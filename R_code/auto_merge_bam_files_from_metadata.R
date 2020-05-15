## Paths to define
md_fname <- "/Users/jon/zon/clonal_hematopoiesis/204590r_2019_11_20/MGH_sample_key_204590.xlsx"
samtools_path <- "/usr/local/bin/samtools"

# Script beginning
md <- read.xlsx(md_fname,1)
md <- md[!is.na(md$fish),,drop=FALSE] #if there are weird blank spaces from excel document

bam_dir <- grep('sorted_bam$', list.dirs(dirname(md_fname), recursive = FALSE), value = TRUE) #bam directory must be labeled 'sorted_bam' in same directory as sample key

#make directory for merged output
if(dir.exists(file.path(bam_dir,"merged_bam"))) {
  merge_output_dir <- file.path(bam_dir,"merged_bam")
} else {
  dir.create(file.path(bam_dir,"merged_bam"))
  merge_output_dir <- file.path(bam_dir,"merged_bam")
}


bam_fnames <- list.files(path = bam_dir, pattern = "*.bam$", full.names = TRUE)

#Convert metadata columns to characters
md$fish <- as.character(md$fish)
md$tissue <- as.character(md$tissue)
md$amp <- as.character(md$amp)
md$tube <- as.character(md$tube)
md$barcode <- as.character(md$barcode)

#Get unique fish samples/tissues and amplicons
md$new_sample_names <- paste(md$fish, md$tissue, sep = "_")
unique_sample_names <- unique(md$new_sample_names)
num_of_samples <- length(unique_sample_names)

unique_amps <- unique(md$amp)
unique_amps <- sort(unique_amps)
num_of_unique_amps <- length(unique_amps)

num_of_iters <- length(md$fish)

no_bam <- NULL
no_bam_idx <- NULL
count <- 1
for(i in 1:num_of_iters){
  if(is.na(md$barcode[i]) == TRUE){
    bam_file_idx <- intersect(grep(md$tube[i], bam_fnames), grep("unknown", bam_fnames))  # when the sample came from a non-barcoded sample, the bam file is labeled as "unknown". If you don't include this, you can sometimes also find the "Tube" ID in the barcodes, which is a problem!
  } else {
    bam_file_idx <- intersect(grep(md$tube[i], bam_fnames), grep(md$barcode[i], bam_fnames))
  }
  
  if(length(bam_file_idx)> 0){ #report if bam file isn't found
    md$bam_file[i] <- bam_fnames[bam_file_idx]
  } else {
    no_bam[count] <- paste0(md$tube[i], "_", md$barcode[i])
    no_bam_idx[count] <- i
    md$bam_file[i] <- NA
    count <- count + 1
  }
}

# sort unique samples into different groups based on amplicons
amps_in_samples <- setNames(data.frame(matrix(ncol = num_of_unique_amps, nrow = num_of_samples)), unique_amps) #this will make logical vector with amplicons as columns and unique samples as rows
for(k in 1:num_of_samples){
  for(n in 1:num_of_unique_amps){
    amps_in_samples[k, n] <- any(md$amp[md$new_sample_names == unique_sample_names[k]] == unique_amps[n])
  }
}


#identify unique sets of samples that should be grouped together
unique_amps_in_samples <- unique(amps_in_samples)
num_of_groups <- nrow(unique_amps_in_samples)

amps_in_samples$group <- NA
unique_amps_in_samples$group <- NA
for(m in 1:num_of_samples){
  for(b in 1:num_of_groups){
    ident_test <- all(amps_in_samples[m,] == unique_amps_in_samples[b,], na.rm = TRUE)
    if(ident_test == TRUE){
      amps_in_samples$group[m] <- b
    }
  }
}

#get indices of groups 
unique_output_groups <- list()
unique_output_group_genes <- list()
for(t in 1:num_of_groups){
  unique_output_groups[[t]] <- which(amps_in_samples$group == t)
  unique_output_group_genes[[t]] <- colnames(unique_amps_in_samples[,which(unique_amps_in_samples[t,] == TRUE), drop = FALSE])
}


for(i in 1:length(unique_output_groups)){
  md_output <- setNames(data.frame(matrix(ncol = 4, nrow = num_of_samples)), c("bamfile", "short.name", "group", "genes"))
  
  count <- 1
  for(j in unique_output_groups[[i]]){
    sample_idx <- grep(unique_sample_names[j], md$new_sample_names)
    sample_idx <- sample_idx[!is.na(md$bam_file[sample_idx])] #get rid of samples that don't have bam files to avoid errors
    
    if(length(sample_idx) > 0){
      combined_fname <- file.path(merge_output_dir, paste0(unique_sample_names[j], "_merged.sorted.bam"), fsep = "/")
      # print(combined_fname)
      # print("")
      # print(unique(md$bam_file[sample_idx]))
      # print("")
      
      cmd <- samtools_path
      system_args <- c('merge', '-f', combined_fname, unique(md$bam_file[sample_idx])) #unique is here so that we don't merge 2 bamfiles that are identical, which may happen if there are no barcodes associated with the experiment
      cmd2 <- samtools_path
      system_args2 <- c('index', combined_fname)
      
      system2(cmd, system_args, wait = TRUE)
      system2(cmd2, system_args2, wait = TRUE)
      
      #autogenerate metadata output
      md_output[count, "bamfile"] <- combined_fname
      md_output[count, "short.name"] <- unique_sample_names[j]
      md_output[count, "group"] <- "xx"
      
      count <- count + 1
    }
    
    
  }
  
  excel_output_fname <- file.path(dirname(bam_dir), paste0("metadata_output_", i,".xlsx"), fsep = "/")
  md_output[1:length(unique_output_group_genes[[i]]), "genes"] <- unique_output_group_genes[[i]]
  
  write.xlsx(md_output,excel_output_fname, row.names=FALSE, showNA=FALSE)
}

print(no_bam)
write.table(no_bam, file = file.path(bam_dir, "qc_no_bam_file_found.txt", fsep = "/"), sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)

print('---------------------------------------------')
print('Completed')
print('')
