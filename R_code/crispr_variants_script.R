## Paths to define
parent_dir <- "/Users/jon/zon/clonal_hematopoiesis/sequencing_analysis/213385r_2020_12_30"
genome_file_path <- "/Users/jon/zon/clonal_hematopoiesis/Jon-33/Jon-33.fa"
genome_target_file_path <- "/Users/jon/zon/clonal_hematopoiesis/Jon-33/R_custom_genome_target_locations_Jon33.bed"

# Function definition
run_crispr_variants <- function(batch_path_name, genome_file_path, genome_target_file_path){
  # Load genome
  genome_file_no_ext <- tools::file_path_sans_ext(genome_file_path)
  
  #Import BED file information of CRISPR target sites
  gd_fname <- genome_target_file_path
  gd <- rtracklayer::import(gd_fname)
  gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = "center") #add 5 bp to each side for visualization
  
  #Get reference sequences for each gene
  genome_seq <- readDNAStringSet(genome_file_path)
  reference <- BSgenome::getSeq(genome_seq,gdl)
  
  #Load metadata table
  f_names <- list.files(path = batch_path_name, pattern = "^metadata_output.*.xlsx$")
  if(length(f_names) < 1){
    stop('Could not find metadata output files in the folder')
  }
  
  sample_names <- lapply(f_names, function(x) file_path_sans_ext(x))
  failed_no_reads <- list()
  failed_no_read_count <- 1
  
  for(i in 1:length(f_names)){
    md <- read.xlsx(file.path(batch_path_name, f_names[i]), 1)
    #Convert metadata columns to characters
    md$bamfile <- file.path(md$bamfile)
    md$short.name <- as.character(md$short.name)
    md$genes <- as.character(md$genes)
    md$group <- as.character(md$group)
    
    #Get targeted regions and delete empty spaces
    injected_gRNA <- md$genes
    md <- subset(md, select = -genes)
    md[md == ""] <- NA
    md <- na.omit(md)
    injected_gRNA[injected_gRNA == ""] <- NA
    injected_gRNA <- na.omit(injected_gRNA)
    
    #Get total number of aligned reads per sample
    total_aligned_reads <- lapply(md$bamfile, function(x) countBam(x,param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))))
    md$total_reads <- sapply(total_aligned_reads, function(x) x$records)
    md <- subset(md,total_reads>0)
    
    if(nrow(md) < 1){
      failed_no_reads[failed_no_read_count] <- file.path(batch_path_name, f_names[i])
      failed_no_read_count <- failed_no_read_count + 1
    } else {
      #Find indices of targeted regions in the gdl GRanges
      gene_idx <- pmatch(injected_gRNA, gdl$name)
      
      #Rearrange groups to individual fish
      samples <- split(md, md$short.name)
      num_of_samples <- length(samples)
      
      #Get or make directory for excel outputs
      bam_dir <- dirname(dirname(md$bamfile[1]))
      bam_dir <- dirname(bam_dir) #for subdirectory with merged files!
      
      if(!dir.exists(file.path(bam_dir,"crispr_sets"))) {
        dir.create(file.path(bam_dir,"crispr_sets"))
        dir.create(file.path(bam_dir,"crispr_sets", sample_names[[i]]))
      }
      
      if(!dir.exists(file.path(bam_dir,"crispr_sets", sample_names[[i]]))){
        dir.create(file.path(bam_dir,"crispr_sets", sample_names[[i]]))
      }
      
      crispr_sets_output_dir <- file.path(bam_dir,"crispr_sets", sample_names[[i]])  
      
      # generate excel outputs
      crispr_set <- list(length(gene_idx))
      wb <- list()
      excel_fname <- list()
      
      min_count <- 25
      min_freq <- 0
      
      for(i in 1:num_of_samples) {
        num_of_groups <- length(samples[[i]]$group)
        crispr_set_fname <- file.path(crispr_sets_output_dir, paste(unique(samples[[i]]$short.name),".RData",sep=""),fsep="/")
        crispr_set <- readsToTargets(file.path(samples[[i]]$bamfile), targets = gdl[gene_idx], references = reference[gene_idx],
                                     names = samples[[i]]$group, target.loc = 22, verbose=TRUE)
        
        #save crispr_set to file
        saveRDS(crispr_set,
                file = crispr_set_fname,
                compress = TRUE)
        
      }
    }
  }
  
  if(length(failed_no_reads) > 0){
    print('Failed: ')
    print(failed_no_reads)  
  }
}


run_crispr_variants(parent_dir, genome_file_path, genome_target_file_path)

print('-----------------------------------------')
print('Completed')
print('')

