#for CRISPR variants script
library(Rsamtools)
library(CrispRVariants)
library(Biostrings)
library(xlsx)
library(tools)
library(tcltk)
library(plyr)
library(reshape2)
library(ggrepel)

## Paths to define before running script

# Path to metadata Excel file
md_fname <- "/Users/jon/zon/clonal_hematopoiesis/204590r_2019_11_20/MGH_sample_key_204590.xlsx"

# Path to file with barcodes
barcode_fw_fasta_file <- "/Users/jon/zon/clonal_hematopoiesis/barcode_fasta/barcode_ppm1d_fw.fasta"

# Path to custom genome file
genome_file <- "/Users/jon/zon/clonal_hematopoiesis/Jon-31/Jon-31.fa"

# Path to cutadapt installation
cutadapt_path = "/Users/jon/Library/Python/3.7/bin/cutadapt"

# Path to Bowtie2 installation
bowtie2_path = "/usr/local/bin/bowtie2"

# Path to samtools installation
samtools_path <- "/usr/local/bin/samtools"


# Function definitions
find_fastq_length <- function(f_name){
  total_lines <- length(count.fields(f_name))
  return(total_lines[1]/4)
}

# Script start
md <- read.xlsx(md_fname,1)
md <- md[!is.na(md$fish),, drop=FALSE]

#Convert metadata columns to characters
md$fish <- as.character(md$fish)
md$tissue <- as.character(md$tissue)
md$amp <- as.character(md$amp)
md$tube <- as.character(md$tube)
md$barcode <- as.character(md$barcode)

tube_list <- unique(md$tube)
tube_list_with_hyphen <- paste(tube_list, "_", sep="")

barcode_fw_fasta <- read.table(barcode_fw_fasta_file, header = FALSE, sep = "\n", stringsAsFactors = FALSE)

barcode_rv_fasta <- barcode_fw_fasta # make a copy of the FW barcodes
for(q in seq(2, nrow(barcode_rv_fasta), by = 2)){ #go by 2 because every other line is barcode ID 
  barcode_rv_fasta[q,1] <- as.character(reverseComplement(DNAString(barcode_fw_fasta[q,1])))
}


genome_file_no_ext <- tools::file_path_sans_ext(genome_file)

#Load sample fastq file names and get directory
fq_dir <- grep('FASTQ$', list.dirs(dirname(md_fname)), value = TRUE) #fastq directory must be labeled 'FASTQ' in same directory as sample key

if(length(fq_dir) < 1){
  message("Error: Could not identify FASTQ directory. Make sure that the FASTQ directory is in same folder as sample key")
  stop()
}

if(length(fq_dir) > 2){
  message("Error: Could not uniquely identify FASTQ directory. Make sure it is the only folder labeled 'FASTQ' in the directory")
  stop()
}

fq_fnames <- list.files(fq_dir, pattern = "*.fastq$", full.names = TRUE)

tube_idx <- sapply(tube_list_with_hyphen, function(x) grep(x, fq_fnames))
tube_idx <- unlist(tube_idx) #in case any tube names are wrong and are not found
fq_fnames <- fq_fnames[tube_idx] #restrict fastq files to only ones that we can find

# fq_dir <- dirname(fq_fnames[1]) #MAC
fq_basename <- basename(fq_dir)
fq_dir_parent <- dirname(fq_dir) #MAC
fq_fnames <- basename(fq_fnames)
sample_names <- gsub(".fastq", "", x=fq_fnames)


#Create trimmed_fastq dir, which will be deleted later
if(dir.exists(file.path(fq_dir_parent,"trimmed_fastq"))) {
  trimmed_fastq_dir <- file.path(fq_dir_parent,"trimmed_fastq",fsep="/") #MAC
} else {
  dir.create(file.path(fq_dir_parent,"trimmed_fastq"))
  trimmed_fastq_dir <- file.path(fq_dir_parent,"trimmed_fastq",fsep="/") #MAC
}

#Create qc folder
if(dir.exists(file.path(fq_dir_parent,"qc"))) {
  qc_dir <- file.path(fq_dir_parent,"qc",fsep="/") #MAC
} else {
  dir.create(file.path(fq_dir_parent,"qc"))
  qc_dir <- file.path(fq_dir_parent,"qc",fsep="/") #MAC
}

#Generate file names for trimmed fastq files
fq_fnames1 <- gsub(".fastq",".trimmed-{name}.1.fastq",fq_fnames)
fq_fnames1 <- file.path(trimmed_fastq_dir,fq_fnames1,fsep="/") #WINDOWS "\\" MAC"/"
fq_fnames2 <- gsub(".fastq",".trimmed-{name}.2.fastq",fq_fnames)
fq_fnames2 <- file.path(trimmed_fastq_dir,fq_fnames2,fsep="/") #WINDOWS "\\" MAC"/"

fq_fnames_full <- file.path(fq_dir_parent,fq_basename,fq_fnames,fsep="/") #WINDOWS "\\" MAC"/"

#Trim fastq by quality using cutadapt method #MAC
before_trim_read_count <- vector()
no_fastq <- list()

for(i in 1:length(fq_fnames)){
  #find corresponding fastq file and barcodes to use in list
  barcode_list <- md$barcode[md$tube == tube_list[i]]
  
  # If all barcodes within a tube are 'NA' (i.e. unbarcoded), then we can still trim the reads without a barcode instead of having to run a different script
  # This 'if' statement tests whether all the barcodes in a given tube are 'NA'.
  # @IMPORTANT! If the tube contains a mix of unbarcoded and barcoded, then this will not work properly. Might output the unbarcoded to 'unknown'.
  
  has_barcode_flag <- FALSE
  if(!all(is.na(barcode_list))){
    has_barcode_flag <- TRUE
    barcode_fasta_ident_idx <- sapply(barcode_list, function(x) grep(x, barcode_fw_fasta[,1]))
    barcode_fasta_seq_idx <- barcode_fasta_ident_idx + 1
    temp_barcode_fasta_idx <- sort(union(barcode_fasta_ident_idx, barcode_fasta_seq_idx))
    temp_barcode_fw_fasta <- barcode_fw_fasta[temp_barcode_fasta_idx,]
    temp_barcode_rv_fasta <- barcode_rv_fasta[temp_barcode_fasta_idx,]
    
    temp_barcode_fw_fasta_file <- file.path(dirname(barcode_fw_fasta_file), "temp_barcode_fw_fasta.fasta", fsep = "/") #make temp fasta file only with barcodes we are looking for
    temp_barcode_rv_fasta_file <- file.path(dirname(barcode_fw_fasta_file), "temp_barcode_rv_fasta.fasta", fsep = "/") #make temp fasta file only with barcodes we are looking for
    
    write.table(temp_barcode_fw_fasta, file = temp_barcode_fw_fasta_file, quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)
    write.table(temp_barcode_rv_fasta, file = temp_barcode_rv_fasta_file, quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)
    
    before_trim_read_count[i] <- find_fastq_length(fq_fnames_full[i])/2
    
    # this is currently written to find the barcode on the 5' end or the reverse complement of the barcode on the 3' end
    cmd <- paste0(cutadapt_path, " --interleaved",
                  " -q 30 -m 50 -g file:", temp_barcode_fw_fasta_file, " -a file:", temp_barcode_rv_fasta_file,
                  " -o ", fq_fnames1[i], " -p ", fq_fnames2[i]," ", fq_fnames_full[i])
  } else if(all(is.na(barcode_list))) {
    # unbarcoded
    before_trim_read_count[i] <- find_fastq_length(fq_fnames_full[i])/2
    
    cmd <- paste0(cutadapt_path, " --interleaved -A XXX -q 30 -m 50",
                  " -o ", fq_fnames1[i]," -p ", fq_fnames2[i]," ", fq_fnames_full[i])
    
  } else {
    print(fq_fnames1[i])
    stop("This tube has a mix of barcoded and non-barcoded samples")
  }
  
  # message(cmd, "\n"); 
  system(cmd, wait = TRUE)
  
  if(has_barcode_flag){
    unlink(temp_barcode_fw_fasta_file)
    unlink(temp_barcode_rv_fasta_file)  
  }
}

# Get trimmed file names and output metadata on trimming
trim_fnames <- list.files(trimmed_fastq_dir, full.names = TRUE)
after_trim_read_count <- sapply(trim_fnames, find_fastq_length)

before_trim_df <- data.frame(fq_fnames_full, before_trim_read_count, row.names = NULL)
colnames(before_trim_df) <- c("file", "count")
after_trim_df <- data.frame(trim_fnames, after_trim_read_count, row.names = NULL)
colnames(after_trim_df) <- c("file", "count")

info_fnames <- file.path(qc_dir, gsub(".fastq", ".qc_info.tsv", fq_fnames), fsep = "/")

for(j in 1:nrow(before_trim_df)){
  trim_meta_output <- before_trim_df[j,]
  trim_idx <- grep(sample_names[j], x = after_trim_df$file)
  
  trim_meta_output <- rbind(trim_meta_output, after_trim_df[trim_idx,])
  
  total_before_read_count <- before_trim_df[j, 2]
  
  trim_percentages <- round(100*trim_meta_output$count/total_before_read_count,2)
  trim_meta_output <- cbind(trim_meta_output, trim_percentages)
  colnames(trim_meta_output) <- c("file", "count", "percentage")
  
  for(k in seq(1, nrow(trim_meta_output), by = 2)){
    if(k != 1){
      trim_meta_output[k, 2:3] <- ""
    }
  }
  
  write.table(trim_meta_output, file = info_fnames[j], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}


if(dir.exists(file.path(fq_dir_parent,"sorted_bam"))) {
  sorted_bam_dir <- file.path(fq_dir_parent,"sorted_bam",fsep="/") #WINDOWS \\ MAC /
} else {
  dir.create(file.path(fq_dir_parent,"sorted_bam"))
  sorted_bam_dir <- file.path(fq_dir_parent,"sorted_bam",fsep="/") #WINDOWS \\ MAC /
}

trim_fnames1 <- trim_fnames[seq(1, length(trim_fnames), by = 2)]
#trim_fnames1 <- trim_fnames1[-grep("unknown", trim_fnames1)] #get rid of unknown barcodes before alignment
trim_fnames2 <- trim_fnames[seq(2, length(trim_fnames), by = 2)]
#trim_fnames2 <- trim_fnames2[-grep("unknown", trim_fnames2)] #get rid of unknown barcodes before alignment

#Map, sort, and index the bam files
sm_fnames <- gsub(".1.fastq",".paired.sam", basename(trim_fnames1))
sm_fnames <- file.path(sorted_bam_dir,sm_fnames,fsep="/") #WINDOWS \\ MAC /
bm_fnames <- gsub(".1.fastq",".paired.bam",basename(trim_fnames1))
bm_fnames <- file.path(sorted_bam_dir,bm_fnames,fsep="/") #WINDOWS \\ MAC /
srt_bm_fnames <- gsub(".paired.bam",".sorted.paired.bam",bm_fnames)

bowtie_info_fnames <- file.path(qc_dir, gsub(".1.fastq", ".bowtie2.info.txt", basename(trim_fnames1)), fsep = "/")

for(i in 1:length(trim_fnames1)) {
  
  cmd1 <- bowtie2_path
  system_args1 <- c("-p", 4, "-x", genome_file_no_ext, "-1", trim_fnames1[i], "-2", trim_fnames2[i],
                    "-S", sm_fnames[i], "--very-sensitive",
                    "2>", bowtie_info_fnames[i])
  
  cmd2 <- samtools_path
  system_args2 <- c("view","-h", "-b", sm_fnames[i], "-o", bm_fnames[i])
  
  cmd3 <- samtools_path
  system_args3 <- c("sort", "-o", srt_bm_fnames[i], bm_fnames[i])
  
  cmd4 <- samtools_path
  system_args4 <- c('index', srt_bm_fnames[i])
  
  system2(cmd1, system_args1, wait = TRUE)
  system2(cmd2, system_args2, wait = TRUE)
  system2(cmd3, system_args3, wait = TRUE)
  system2(cmd4, system_args4, wait = TRUE)
  unlink(bm_fnames[i])
  unlink(sm_fnames[i])
}

# Delete trimmed fastq files. We will only save the aligned files since trimming is fast.
unlink(trimmed_fastq_dir, recursive = TRUE)

print('----------------------------------------')
print('Completed')
print('')

