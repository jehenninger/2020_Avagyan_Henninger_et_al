# 2020_Avagyan_Henninger_et_al
MATLAB, Python, and R code for "Distinct paths to clonal hematopoiesis with acquired or germline mutations"


## R code notes for CRISPR mutation analysis

### Order of running scripts

1. *barcode_trimming_and_aligning.R*  
This script will read an Excel metadata file, trim reads for barcodes and quality, and then align to a custom genome to produce sorted bam files
        
    Filepaths you must define:  
    *md_fname* = Path to metadata Excel file (example in ./R_code/bin). Fastq files *must* be in a directory called "FASTQ"  
    *barcode_fw_fasta_file* = Path to fasta file with barcodes (example in ./R_code/bin)  
    *genome_file* = Path to fasta file of custom genome that has bowtie2 indexes in same folder (example in ./R_code/bin/Jon-31)  
    *cutadapt_path* = Path to cutadapt installation  
    *bowtie2_path* = Path to bowtie2 installation  
    *samtools_path* = Path to samtools installation  

2. *auto_merge_bam_files_from_metadata.R*  
This script will use the same metadata file to merge reads from different tubes according to the metadata, and it will rename sorted bam files to IDs given in the metadata.
        
    Filepaths you must define:          
    *md_fname* = Path to metadata Excel file (example in ./R_code/bin)  
    *samtools_path* = Path to samtools installation  
            
3. *crispr_variants_script.R*  
This script will use the merged, sorted bam files and run CrispRVariants to detect mutants.
    
    Filepaths you must define:  
    *parent_dir* = Parent directory of data (containing folders of FASTQ, sorted_bam, etc...)  
    *genome_file_path* = Path to fastafile of custom genome (example in ./R_code/bin/Jon-31)  
    *genome_target_file_path* = Path to .bed file that has target regions of gRNAs (example in ./R_code/bin/Jon-31)  
    
4. *group_crisprSets_by_common_alleles.R*  
This script will read the RData objects produced by crispr_variants_script and output Excel files that contain frequencies and counts of all variants found in a group of crisprSets.
    
    Filepaths you must define:  
    *parent_dir* = Parent directory of crispr_sets data  