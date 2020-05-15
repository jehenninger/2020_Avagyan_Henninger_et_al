## Paths to define
parent_dir <- '/Users/jon/zon/clonal_hematopoiesis/204590r_2019_11_20/crispr_sets'

# Function definition
group_crispr_sets_by_common_alleles <- function(batch_path_name){
  dir_list <- list.dirs(batch_path_name, recursive=FALSE)
  
  num_of_dir <- length(dir_list)
  
  if(num_of_dir > 1){ #if you just want the combination of all files in the crispr_set directory, then the dir_list will only be the parent path
    k_start <- 1
  } else{
    if(length(dir_list) == 0){
      dir_list <- list(batch_path_name)
      k_start <- 1
    } else {
      dir_list <- list(dir_list)
      k_start <- 1
    }
  }
  
  no_mutation_count <- 1
  no_mutation_found_output <- list()
  
  for(k in k_start:length(dir_list)){
    path_name <- dir_list[[k]]
    sample_name <- basename(path_name)
    
    f_names <- list.files(path = path_name, pattern = "*.RData")
    sample_names <- lapply(f_names, function(x) file_path_sans_ext(x))
    
    variant_data <- list()
    variant_data_freq <- list()
    variant_data_initial <- list()
    joined_crispr_set_counts <- list()
    joined_crispr_set_freqs <- list()
    gene_name <- NULL
    joined_data_counts <- NULL
    joined_data_freqs <- NULL
    failed_crispr_sets <- list()
    
    failed_count <- 1
    for(i in 1:length(f_names)){
      crispr_data <- readRDS(file.path(path_name, f_names[i]))
      num_of_crispr_sets <- length(crispr_data)
      
      if(num_of_crispr_sets < 1) {
        failed_crispr_sets[failed_count] <- file.path(path_name, f_names[i])
        failed_count <- failed_count + 1
      } else {
        count <- 1
        for(j in 1:num_of_crispr_sets){
          gene_name[j] <- crispr_data[[j]]$target$name
          
          variant_data_initial[[j]] <- variantCounts(crispr_data[[j]], min.count = 10, min.freq = 0.05, include.nonindel = TRUE, result = "counts")
          
          no_var_test <- grep('no var', rownames(variant_data_initial[[j]]))
          if(length(no_var_test) < 1){ #this is if there are no "no variant" reads but there are indel or SNP reads
            variant_data_initial[[j]] <- rbind(0, variant_data_initial[[j]])
            rownames(variant_data_initial[[j]])[1] <- "no variant"
          }
          
          # get sum of no var counts
          no_var_read_count <- sum(variant_data_initial[[j]][grep('no var',rownames(variant_data_initial[[j]]))])
          
          # get sum of SNV counts (to add into 'no var' counts)
          snv_index <- grep('SNV', rownames(variant_data_initial[[j]]))
          snv_read_count <- sum(variant_data_initial[[j]][snv_index])
          
          total_wt_count <- no_var_read_count + snv_read_count
          variant_data_initial[[j]][grep('no var', rownames(variant_data_initial[[j]])), 1] <- total_wt_count #no var will now add SNV's to total read count
          
          # get rid of SNV and "Other" reads
          other_snv_index <- grep('SNV|Other', rownames(variant_data_initial[[j]]))
          if(length(other_snv_index) > 0){
            variant_data[[j]] <- variant_data_initial[[j]][-other_snv_index,, drop = FALSE]
          }else{
            variant_data[[j]] <- variant_data_initial[[j]]
          }
          
          # get total indel count
          total_indel_count <- sum(variant_data[[j]][-grep('no var', rownames(variant_data[[j]]))])
          
          # distinguish frame-shift
          indel_cigar <- rownames(variant_data[[j]][-grep('no var', rownames(variant_data[[j]])), , drop = FALSE])
          
          if(length(indel_cigar) > 0){ #if there are indels
            #Get rid of locations, colons, and commas from modified CIGAR strings
            #Isolate numbers from strings to detect if frameshift or in-frame
            indel_cigar <- gsub("^[^:]*:|,.*:","",indel_cigar)
            indel_cigar <- strsplit(indel_cigar,"[DI]")
            indel_cigar <- lapply(indel_cigar,function(x) sum(as.numeric(x)))
            indel_cigar <- unlist(indel_cigar)
            fs_test <- indel_cigar%%3>0
            fs_test <- c(FALSE, fs_test) #need to add a FALSE because we took out no_var before
            
            # append gene name to each variant
            gene_specific_allele_name <- rownames(variant_data[[j]])
            for(n in 1:length(gene_specific_allele_name)){
              if(n == 1){
                gene_specific_allele_name[n] <- paste(gene_name[j], '0', gene_specific_allele_name[n], sep = "_")
              } else {
                if(fs_test[n] == TRUE){
                  gene_specific_allele_name[n] <- paste(gene_name[j], 'fs', gene_specific_allele_name[n], sep = "_")
                } else {
                  gene_specific_allele_name[n] <- paste(gene_name[j], 'nfs', gene_specific_allele_name[n], sep = "_")
                }
              }
              
            }
          } else { #if there are not indels
            gene_specific_allele_name <- rownames(variant_data[[j]])
            gene_specific_allele_name <- paste(gene_name[j], '0', gene_specific_allele_name, sep = "_")
          }
          
          variant_data[[j]] <- cbind(gene_specific_allele_name, as.data.frame(variant_data[[j]]))
          
          # make column header
          variant_column_names <- colnames(variant_data[[j]])
          for(k in 1:length(variant_column_names)){
            if(k > 1){
              variant_column_names[k] <- paste("fish", sample_names[[i]], variant_column_names[k], sep = "_")
            } else {
              variant_column_names[k] <- "allele"
            }
          }
          colnames(variant_data[[j]]) <- variant_column_names
          
          # calculate frequencies and make frequency table
          variant_data_freq[[j]] <- variant_data[[j]]
          freq <- round(variant_data_freq[[j]][,2] * 100 / (total_wt_count + total_indel_count), 3)
          variant_data_freq[[j]][,2] <- freq
          
          # merge crispr sets across genes
          if(count > 1){
            joined_crispr_set_counts[[i]] <- rbind(joined_crispr_set_counts[[i]], variant_data[[j]], make.row.names = FALSE, stringsAsFactors = FALSE)
            joined_crispr_set_freqs[[i]] <- rbind(joined_crispr_set_freqs[[i]], variant_data_freq[[j]], make.row.names = FALSE, stringsAsFactors = FALSE)
          } else {
            joined_crispr_set_counts[[i]] <- variant_data[[j]]
            joined_crispr_set_freqs[[i]] <- variant_data_freq[[j]]
            count <- count + 1
          }
          #} 
        }
        ############# Joining data ###################
        #Counts
        if(i > 1){
          joined_data_counts <- merge(joined_data_counts, joined_crispr_set_counts[[i]], all = TRUE)
        } else {
          joined_data_counts <- joined_crispr_set_counts[[i]]
        }
        
        #Freqs
        if(i > 1){
          joined_data_freqs <- merge(joined_data_freqs, joined_crispr_set_freqs[[i]], all = TRUE)
        } else {
          joined_data_freqs <- joined_crispr_set_freqs[[i]]
        }
        
      }
    }
  }

  joined_data_counts <- joined_data_counts[order(as.character(joined_data_counts$allele)),]
  joined_data_freqs <- joined_data_freqs[order(as.character(joined_data_freqs$allele)),]
  
  if(length(crispr_data) > 0){
    for(m in 1:length(gene_name)){
      max_gene_idx <- max(grep(gene_name[m], joined_data_counts$allele)) #may be -Inf if the gene can't be found
      if(!is.infinite(max_gene_idx)){
        joined_data_counts <- rbind(joined_data_counts[1:max_gene_idx,], NA, NA, joined_data_counts[-(1:max_gene_idx),]) #adds empty rows between genes
        joined_data_freqs <- rbind(joined_data_freqs[1:max_gene_idx,], NA, NA, joined_data_freqs[-(1:max_gene_idx),])
      }
    }
    
    write.table(joined_data_counts, file = file.path(path_name, paste(sample_name, "combined_data_counts.tsv", sep = "_")), quote = FALSE, row.names = FALSE, sep = "\t", na = "")
    write.table(joined_data_freqs, file = file.path(path_name, paste(sample_name, "combined_data_freqs.tsv", sep = "_")), quote = FALSE, row.names = FALSE, sep = "\t", na = "")
    
    excel_output_file <- file.path(path_name, paste(sample_name, "freq_and_counts.xlsx", sep = "_"))
    write.xlsx(joined_data_freqs, file = excel_output_file, sheetName = "freq", col.names = TRUE, row.names = FALSE, append = FALSE, showNA = FALSE)
    write.xlsx(joined_data_counts, file = excel_output_file, sheetName = "count", col.names = TRUE, row.names = FALSE, append = TRUE, showNA = FALSE)
  }
  
  if(length(failed_crispr_sets) > 0) {
    print('Failed crispr_sets: ')
    print(failed_crispr_sets)
  } else{
    print('No failed crispr sets')
    
  }
}

# Begin script


dir_list <- list.dirs(parent_dir, recursive=FALSE)

if(length(dir_list) == 0){
  dir_list = list(parent_dir)
}

for(f in c(1:length(dir_list))){
  group_crispr_sets_by_common_alleles(dir_list[[f]], include_long_graphs=include_long_graphs, time_points=time_points, my_order=order)
}

print('----------------------------------------------------')
print('Completed')
print('')