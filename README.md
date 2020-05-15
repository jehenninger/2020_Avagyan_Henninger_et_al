# 2020_Avagyan_Henninger_et_al
MATLAB, Python, and R code for "Distinct paths to clonal hematopoiesis with acquired or germline mutations"

## MATLAB code notes for Zebrabow color analysis
### Installation  

1. Place the folder containing the code in your Matlab path and add folder to path (Use the 'Set Path' icon in the main Matlab toolbar). Alternatively, extract to any folder and set the Matlab path to include this folder using the same 'Set Path' tool.
	
### Dependencies  

1. Prior to running, you should have the following Matlab toolboxes installed:
    'Image Processing Toolbox'
    'Statistics and Machine Learning Toolbox'

2. Also need to install:  
    logicleTransform - https://www.mathworks.com/matlabcentral/fileexchange/45022-logicle-transformation/content/run_logicle.m  
    
    Follow instructions on that website to build the MEX-functions.

### Before running the program for the first time  

Change FACS parameter names:  

In command console, type 'open zbow_analysisV2'

Under 'loadFSCFile' function, change parameter names corresponding to your FACS machine. Also keep in mind whether you are exporting compensated data. (i.e. our parameters are 'FJComp-PE-A', 'FJComp-FITC-A',etc...). The parameter must match exactly. If you are unsure what the parameter names are, do the following:    

  * In command console, type 'findFACSParameterNames'. Then choose to open a fcs file with the proper parameters. The names of the parameters will then be displayed on the console.  
    
In this function, you can also modify the red, green, and blue width parameters for the biexponential (logicle) transform. (our custom parameters are set for R = 1.5, G = 1.75, and B = 1.75).  
	
### Suggestions prior to running function:

1. Export high quality (i.e. single, live cells of population of interest) FCS files from FlowJo or similar program. This will make clustering much easier. We generally only export the parameters that we will need for clustering (i.e. the color data). Exporting all parameters won't affect anything.

2. If at any time there is an error, it is best to try restarting the GUI as the data structures may be modified.
	
### Running the GUI

1. In MATLAB command console, type 'run_zbow'. This will bring up the GUI.

2. Before loading a file, you can change how much data is sub-sampled from the FACS data. This is also the data that will be used for clustering. It is by default set to 10000 data points. This number provides a good balance between performance and accuracy. On better machines (with more RAM and better CPUs), the algorithms will work on 20-30,000 points, but it will take longer to analyze.

3. To load FCS file, click the "Load" button and select a FCS file. The loaded file name will appear under 'Loaded File'. A 3D scatter and ternary plot should be generated in separate figures. You can use Matlab's figure tools (like rotation) to manipulate viewing the data.

4. In the bottom left of the GUI, there are 3D Scatter and Ternary Graph Controls. These control both the color and display type for each graph.

   Options for color control include:
    * Custom normalized - normalized color from Red, Green, and Blue values using the logicle transformation with custom width parameters (see above for custom values).
    * Default normalized - normalized color from Red, Green, and Blue values using the logical transformation with default width parameters (calculated for each parameter; these will generally match what you would see on the FACS machine).
    * Linear - normalized color from Red, Green, and Blue using the linear intensity values from FACS. Because the dynamic range of FACS is so large, this option isn't very helpful since most cells will look dim in color. Could be more appropriate for microscopy intensity values (although the function currently does not support this. Could add later).
    * Cluster color - pseudocolored clusters that are easy to distinguish. If clustering hasn't been done yet, all cells are considered part of the same cluster.
    
   Options for display control include:
    * Custom normalized - normalized values using the logicle transformation with custom width parameters (see above for custom values).
    * Default normalized - normalized values using the logical transformation with default width parameters (calculated for each parameter; these will generally match what you would see on the FACS machine).
    * Linear - normalized values using the linear intensity from FACS.
        
5. If there are outliers in the data, these can be removed using the 3D Scatter Plot. In the 3D Scatter plot, use the 'Brush/Select data' tool to highlight the outliers. Then, right click the highlighted points and select 'Remove' (or simply hit the Delete key). NOTE: This is irreversible, so if you make a mistake you have to 'Reset' and load the data again. After removing the outliers, select the 'Remove Outliers' button in the GUI. This will update the data, and it will tell you how many cells are remaining after outlier removal.

6. To start clustering the data, select the "Make Decision Graph" button in the GUI. It will calculate the rho and delta values for the 'Decision Graph' (see Clustering by fast search and find of density peaks, Rodriguez and Laio, Science 2014). Clustering can be modified with the "Pre-clustering panel". Here, you can decide on the transformation type to cluster the data (default is the custom logicle transformation), and you can decide whether to cluster on the 3D values (R,G,B) or the 2D values (x,y of the Ternary plot). Although clustering in 2D removes a degree of freedom, we have found that clustering works best with this method (set as default). If these options are changed after clustering, then you must re-select "Make Decision Graph" to re-cluster.

7. Cluster centers must be chosen manually. To add cluster centers, select the "Add Center" button. This will change your cursor to a cross-hair, and it will allow you to draw a rectangle on the Decision Graph to choose centers. It is not necessary to try to get all the centers at once, as you can click "Add Center" as many times as you need. It will only add new, unique centers to the list. The "Cluster List" will populate with the cluster numbers (and # of cells in that cluster). The bar graph next to the ternary plot will also update to show the cluster distribution. To remove a cluster from the list, highlight it in the Cluster List and select "Remove Center". To highlight a specific cluster in the 3D Scatter and Ternary Graphs, select a cluster in the Cluster List and then select the "Highlight cluster" button. This will set the cluster of choice to the color determined in the Color control panels, and it will set every cell not in the cluster to a gray value. This function is useful for determining if a cluster should be split into multiple clusters.

8. To evaluate the clustering solution, select the "Evaluate clustering (silhouette)" button. This will generate a new graph, which plots the silhouette value for each cell in a cluster. The black bar represents the mean silhouette value for that cluster.

9. The "Copy cluster data to clipboard" button will copy the cluster number, # cells in the cluster, and the percentage of total cells to the system clipboard of the cluster highlighted in the cluster list. This is a quick method to copy the data for a specific cluster to Excel.

10. To save the files, select the "Save Session" button. It will ask you for a folder to save to. In this folder, it will generate high quality JPEGs and EPS files (vector format compatible with Adobe Illustrator) of the decision graph, ternary graph, bar graph, and silhouette graph (if it exists). It will also generate 2 Excel (.xlsx) files: a Summary file that contains a summary of the clusters, including number of cells, percentage of total, and mean RGB color; a RawData file that contains all of the normalized, raw data (with headers).
    NOTE: It automatically uses the sample name to generate file names. This will overwrite any JPEG, EPS, or XLSX files with identical names in the folder you choose. 
    NOTE: The graphs are copied exactly as they appear. So you will want to set the Color and Display controls to how you want it to look prior to saving.
    
11. The "Reset" button will close all associated figures and re-launch the GUI to load a different set of data.


## Python code notes for Zebrabow color analysis

### README under development

Running main.py in a Terminal window should bring up a GUI where you can load .fcs files. Clustering is done using HDBSCAN. Clusters can be split or joined, and you can double click a cluster ID in the cluster table to specifically highlight that cluster.

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