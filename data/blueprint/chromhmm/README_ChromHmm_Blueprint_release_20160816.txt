Analysis name: ChromHmm Segmentation of ChIP-Seq data

Analysis description:
  A model with 12 states was generated for 91 samples from healthy donors that had the complete histone modifications set and the input (information about samples used to generate the model summarized in table "index.segmented.healthy.model.csv").
  This model was then used to segment the 91 cell types used in model plus 82 additional cell types with the complete histone modifications set and the input, including 7 cell lines, 9 additional healthy donors and samples from 66 patients with diseases of the haematopietic system.


Data format:

  Files:
    index.segmented.healthy.model.csv			   : Contains information about segmented cell types from healthy donors included in model
    index.segmented.healthy.csv				   : Contains information about segmented cell types from healthy donors not included in model
    index.segmented.disease.csv				   : Contains information about segmented samples from diseases
    index.segmented.cell_lines.csv			   : Contains information about segmented cell lines
    ChromHMM_model_Blueprint_release_201508.md5sum	   : A md5 file to check files integrity 
    model_12_12_Blueprint_release_201608.txt               : The model generated with 12 states
    emissions_12_12_Blueprint_release_201608.txt           : The emission probabilities file tab-delimiter
    emissions_12_12_Blueprint_release_201608.png/svg       : Figures with the emission probabilities
    transitions_12_12_Blueprint_release_201608.txt         : The transition probabilities file tab-delimiter
    transitions_12_12_Blueprint_release_201608.png/svg     : Figures with the transition probabilities
    POSTERIOR.tar.gz     	   			   : Includes 4 folders: POSTERIOR_healthy_model, POSTERIOR_healthy, POSTERIOR_disease and POSTERIOR_cell_lines, and a md5 file for each folder
							     Folders include all the posterior probabilities for each 200bp windows,
                                                             25 files (one each for chr1-22,X,Y and M) per sample with posterior probabilities. 
                                                             Format: TAB-delimited file of the probability of each state in each bin
    STATEBYLINE.tar.gz   				   : Includes 4 folders: STATEBYLINE_healthy_model, STATEBYLINE_healthy, STATEBYLINE_disease and STATEBYLINE_cell_lines, and a md5 file for each folder
							     Folders include the state with the highest posterior probability for each 200bp windows 
    SEGMENTATION.tar.gz  				   : Includes 4 folders: SEGMENTATION_healthy_model, SEGMENTATION_healthy, SEGMENTATION_disease and SEGMENTATION_cell_lines, and a md5 file for each folder
							     Folders include the segmentation files, one bed file per sample with genomic segmentation, and a md5 file for each folder
                                                             Format: 9 coloums BED file
    VISUALIZATION.tar.gz 				   : Includes 4 folders: VISUALIZATION_healthy_model, VISUALIZATION_healthy, VISUALIZATION_disease and VISUALIZATION_cell_lines, and a md5 file for each folder
 					                     Folders include two bed files per sample for UCSC visualization. 
                                                             Format: 11 columns BED file

    More details about formats in: http://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf 



#####
#####

Command line:


	1. Data Binarization: Chip-seq data is binarized in 200bp with BinarizeBam function from ChromHmm package v1.11 (BinarizeBam and BinarizeBed functions produce the same binary output)
	
	 java -mx20G -jar /<path>/ChromHMM.jar BinarizeBam -c <controldir> /<path>/ChromHMM/CHROMSIZES/hg38.txt <inputbeddir> <datafile>  <outputbinarydir>


   Files required:
      A. Control aligment files
      B. Histone modification aligment files
      C. File with chromosome sizes. A two column tab delimited file with the first column being the chromosome and the 
         second being the chromosome length (hg38.txt can be found in CHROMSIZES directory where ChromHmm was installed)
      D. Data file. A tab delimited file each row contains the cell type or other identifier for a groups of marks, 
         then the associated mark, then the name of a bed file, and optionally a corresponding control bed file. Example:

           C000S5H1	H3K4me3	C000S5H1.ERX149093.H3K4me3.dedup.bwa.GRCh38.20150504.bam	C000S5H1.ERX675826.Input.dedup.bwa.GRCh38.20150503.bam
	   C000S5H1	H3K36me3	C000S5H1.ERX149094.H3K36me3.dedup.bwa.GRCh38.20150503.bam	C000S5H1.ERX675826.Input.dedup.bwa.GRCh38.20150503.bam
	   C000S5H1	H3K27ac	C000S5H1.ERX149095.H3K27ac.dedup.bwa.GRCh38.20150504.bam	C000S5H1.ERX675826.Input.dedup.bwa.GRCh38.20150503.bam


	 
	2. Learn Model and genomic segmentation of the samples used in model: ChromHmm package v1.11

        java -mx20G -jar ChromHMM.jar LearnModel -i <outfileID> -l /<path>/ChromHMM/hg38.txt -printposterior -printstatebyline  <BINARIZED_files_dir> <outputdir> numstates <assembly_dir>

    Files required:

      A. File with chromosome sizes. A two column tab delimited file with the first column being the chromosome and the 
         second being the chromosome length (hg38.txt can be found in CHROMSIZES directory where ChromHmm was installed)
      B. The Binarized files generated in the 1st step
      C. The assembly directory path with the annotations
	
	
	3. Segmentation of the remaining samples using the 12-state model previously generated:
	java -mx20G -jar ChromHMM.jar  MakeSegmentation -printposterior -printstatesbyline /<path>/modelfile <inputdir> <outputdir>

	File required: 
	A. File with the model parameters to produce the segmentation

	
	4. Generation of browser_visualization files for the cell types(?) segmented in step3:
	java -mx20G -jar ChromHMM.jar MakeBrowserFiles /<path>/segmentfile segmentationname /<path>/outputfileprefix

	Files required: 
	A. Segmentation files generated in previous step 





ANNOTATION FOR THE 12 STATES

State1.- Repressed Polycomb High signal H3K27me3
State2.- Repressed Polycomb Low signal H3K27me3 
state3.- Low signal
State4.- Heterochromatin High Signal H3K9me3
State5.- Transcription High signal H3K36me3
State6.- Transcription Low signal H3K36me3
State7.- Genic Enhancer High Signal H3K4me1 & H3K36me3
State8.- Enhancer High Signal H3K4me1
State9.- Active Enhancer High Signal H3K4me1 & H3K27Ac
State10.- Distal Active Promoter (2Kb) High Signal H3K4me3 & H3K27Ac & H3K4me1
State11.- Active TSS High Signal H3K4me3 & H3K4me1
state12.- Active TSS High Signal H3K4me3 & H3K27Ac

    More details about ChromHmm software in: http://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf
    
Source code: 

    ChromHMM package v1.11
    http://compbio.mit.edu/ChromHMM/

Publication url: 

    http://www.nature.com/nmeth/journal/v9/n3/full/nmeth.1906.html

Centre:
 
    Spanish National Cancer Research Center- Centro Nacional de Investigaciones Oncológicas (CNIO). Madrid, Spain.

Contact:

	Felipe Were <fnicolau@cnio.es>
	Enrique Carrillo <ecarrillo@cnio.es>    
	 

