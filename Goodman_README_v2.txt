This package was developed to process the raw sequencing data from INSeq experiments generated on an Illumina sequencer, map the transposon insertion locations, count the numbers of reads associated with respective locations, and generate output files of the analyzed data.

Copyright 2011 Jeffrey Gordon laboratory.
This program is distributed under the terms of the GNU General Public License.

************
Requirements:
************
(i) A computational cluster or a desk-top/lap-top running Mac/Unix with minimum 4GB memory. Installation of perl and Bowtie (0.12.7 or higher) is required. 
(ii) The reference genome sequence file in fasta format and gene annotation file (e.g. ptt file from NCBI that provides gene locations along the genome for protein coding genes) for the target species. The genome sequence file can contain sequences for more than one DNA molecule. Each .ptt file contains annotation information for only one chromosome or plasmid. This should be the case if the .ptt file(s) are downloaded from NCBI. The name of the .ptt file should be consistent with the name of the corresponding DNA molecule in the fasta header in the fasta file. For example:
        In the genome fasta file of B.thetaVPI5482, there are sequences for both the main chromosome and the plasmid with the following fasta header lines:
         >Bacteroides_thetaiotaomicron_VPI-5482
         >NC_004703

The corresponding .ptt files present in the indexes directory should be named: Bacteroides_thetaiotaomicron_VPI-5482.ptt,
NC_004703.ptt

NOTE: Please make sure the names are consistent between the ptt files and the header names in the fasta files. You can change either of them to make them match. Please avoid spaces in names.

**************
Getting started (note that a description of the files in INSEQ_analysis and INSEQ_demo is provided at the bottom of this document):
**************   
  (1) Download the analysis package from the website (URL),unzip it (which you probably have already done since you are reading this Readme.txt file) 

##### URL will be provided at time of publication ######

  (2) Download bowtie from http://bowtie-bio.sourceforge.net/index.shtml and install properly to your system, then edit the "config.txt" in the analysis package. Change the bowtie_dir variable to the path where your bowtie installation is.

      e.g. /srv/cgs/local/bowtie/latest (UNIX)
      e.g. /Users/mwu/bowtie-0.12.7  (Mac)

  If you are not sure what the path is, at the command line, go into the bowtie directory, then type 'pwd', the output of that command is the path of your bowtie directory.

  (3) Go to the directory where you unpacked the analysis package (the INSeq_analysis package, not bowtie package), go to the "indexes" directory, create a new directory with your organism's name, e.g. BthetaVPI_5482. Put the multifasta file and the .ptt files in the new directory, go into the new directory and construct a bowtie index by typing the command:

   "bowtie_dir"/bowtie-build <fasta file>  <name of index, same as the new directory>

   one example:
   /Users/mwu/bowtie-0.12.7/bowtie-build Bacteroides_thetaiotaomicron_VPI-5482_wP.fna BthetaVPI_5482
  
  (4) Go to the "Example" directory, type the command: perl ../INSeq_pipeline.pl -i Sample_rawreads.txt -m Sample_mapping_file.txt -s BthetaVPI_5482 -d 0.8. Two mappingjobsXXX.job files will be generated and please run them in series by typing "sh mappingjobsXXX.job" for each.

  If there is no error message, you are set to run the analysis package. 

  NOTE: The data stored in "Example" is a small (artificial) file and a test analysis using this data will produce a large number of genes with no hits. There is a directory in Example called "Results_expected". You can check your results by comparing the output files in that directory to make sure everything function properly. To test the pipeline on a larger dataset, sample data is available from http://gordonlab.wustl.edu/SuppData.html

  (5) Create a directory for the experiment where you would like to store your analysis results.  For example:

      mkdir /Users/mwu/Documents/Experiment_1
      
  (6) Create a mapping file that contains the barcode for each sample in tab-delimited format in the experiment directory:
     <barcode>       <Sample Name>

  A sample mapping file (Sample_mapping_file.txt) is included in the "Example" directory. Tip: If the mapping file is imported from document or spreadsheet software such as Microsoft Word or Excel, be sure to remove any hidden characters.

  (7) From the experimental folder, run the analysis pipeline by running the wrapper script INSeq_pipeline.pl as in the usage:

perl "path for analysis package"/INSeq_pipeline.pl -i <the raw reads file>  -m <Barcodes mapping file>  -s <indexed genome name>  -d <length_disrupt_percent (max=1)>[-operon -c <operon_probability_cutoff (max=1)>] [-arrayed] 

Required arguments: -i gives the input raw reads file, -m gives the mapping file, -s gives the name of the indexed genome that the reads should be mapped to 
Optional arguments: -d gives the region of the gene in which insertions are expected to disrupt gene function. The default is 1, which means when insertion falls anywhere in the gene (100%), the gene function will be disrupted. Setting the -d argument to 0.9, for example, would exclude insertions in the distal 10% of the gene when calculating the total number of reads/insertions for that gene. -operon (no argument) specifies that putative downstream (polar) insertions should be calculated based on a user-provided operon probability file. -c is the cutoff for operon probability, default is 1, which means only when the probability of two genes being in an operon is equal to 1 (100%), the polar effect will be considered for the downstream gene. Setting the -c argument to 0.8, for example, will calculate a polar effect for the downstream gene if the probability of the genes being in an operon is at least 0.8 (80%). -arrayed (no argument) is the option for the arrayed library.

For example:
   perl /Users/mwu/INSEQ_analysis/INSeq_pipeline.pl -i  ../Example/Sample_raw_data.txt -m ../Example/Sample_mapping_file.txt -s BthetaVPI_5482 -d 0.8

Be sure that your raw reads file matches the format of the example raw reads files provided with the software pipeline.

(8) The INSeq_pipeline.pl will generate several mappingjobsXXX.job under the experiment directory. On a desktop/laptop, run them in series by typing "sh mappingjobsXXX.job" for each job, or (on a cluster) run in parallel.


  ************
  Output files
  ************

Several output files will be produced in the results directory:
   INSEQ_experiment.scarf_assigned.txt  :Reads are assigned to different samples based on the sample-specific barcode.
   INSEQ_experiment.scarf.log           :Some statistics for the analysis process, including 1) the total number of reads, 2) the percentage of reads being mapped and trimmed (that have the transposon), 3) how many insertions in the sample (coverage) with how many raw reads (reads after filtering in the normalization step), 4) the scale factor being used for normalization.

   Several output files will be placed in the results folder:
   INSEQ_experiment.scarf_Samplename.bowtiemap :Raw mapping output file from bowtie
   INSEQ_experiment.scarf_Samplename.bowtiemap_processed.txt_chromosomename: Txt file containing the processed mapping output, with format
<chromosome name> <insertion position> < reads mapped to the left side of the insertion> <reads mapped to the right side of the insertion> <the total number of reads mapped to that position>
   INSEQ_experiment.scarf_Samplename.bowtiemap_processed.txt_chromosomename_filter_cpm.txt: Txt file with positions having more than 3 total reads, and normalized to counts per million reads.
   INSEQ_experiment.scarf_Samplename.bowtiemap_processed.txt_chromosomename_filter_cpm.txt_mapped: Txt file that lists genes mapped by insertions, only when the insertions located at the proximate XX percentage (specify by -d option) are considered interrupt the function of genes. The format is
     <gene name>\t<the total number of unique insertions in the gene>\t<the sum of normalized read counts in that gene>\t<gene annotation from the ptt file>

(9) Clean-up
Since there are many intermediate files generated along the pipeline, if storage space is a big concern, deleting them is optional.

You can run clean_up.sh from the experiment folder by using command "sh clean_up.sh", which will remove the sorted barcoded assigned reads file, the folder "bcsortedseqs" and the mappingjob files.

****************
Troubleshooting
****************
(1) If you see the error message "bowtie-build: No such file or directory" when your run build-bowtie, check that you included the entire path for bowtie when you run the command, e.g."/Users/mwu/bowtie-0.12.7/bowtie-build". If you get the error message "bowtie: No such file or directory" when you run the pipeline, please check the "config.txt" in step (2) to make sure you changed the bowtie_dir variable. It should be in the format 
       bowtie_dir="/Users/mwu/bowtie-0.12.7".
(2) If you see the error message "could not find XXXX (name of the reference genome)", please check that you made the index files for your mutant library UNDER the "indexes" directory (should be in the position exactly inside of the "indexes" folder), and please check that the name of the index, which is the prefix of the .X.ebwt file, is consistent with the name of the directory. If not, please re-make the ebwt file using bowtie-build and repeat step (3).                                             
(3) If you find in the results folder, the "INSEQ_experiment.scarf_Samplename.bowtiemap_processed.txt_chromosomename_filter_cpm.txt" lists all the different insertions with normalized read counts, but the "INSEQ_experiment.scarf_Samplename.bowtiemap_processed.txt_chromosomename_filter_cpm.txt_mapped" shows all genes having zero insertions, please check if the header names of the chromosome are consistent with the names of the ptt file as described in Requirements (ii). Check the NOTE there to make sure you did not introduce any spaces in any names and keep the naming consistent.
(4) If the analysis pipeline cannot recognize your raw reads file, check that the format of the file matches that of the sample raw reads files provided with the data analysis pipeline.



Additional features of the analysis pipeline:

(A) Analyze the insertions with operon information. This allows polar effects to be calculated for downstream genes within an operon.

If users have operon information for the genome, this information should be organized as follows:

Operon files must follow the same organization as gene annotation files: Separate files for each chromosome with the same name as the .ptt files and using .operons as the extension. e.g. Bacteroides_thetaiotaomicron_VPI-5482.operons
NC_004703.operons

The format for the operon file is:
   Tab-delimited operon file in format
        Gene A          PA,B   
        Gene B          PB,C
        Gene C          PC,D
   Where PA,B is a value between 0 and 1 (inclusive) reflecting the probability that Gene A and Gene B belong to the same operon.

Put the operon files in the same directory as .ptt files and then invoke the operon function by adding -operon -c <operon_probability_cutoff (max=1)> to the command line. One example is:
 perl /Users/mwu/INSEQ_analysis/INSeq_pipeline.pl -i  ../example/Sample_raw_data.txt -m ../example/Sample_mapping_file.txt -s BthetaVPI_5482 -d 0.8 -operon -c 0.2

The pipeline will generate some additional output files based on the operon information, which are also saved in the results folder:
   INSEQ_experiment.scarf_Samplename.bowtiemap_processed.txt_chromosomename_filter_cpm.txt_genes_direct.txt: The txt file that lists genes mapped by insertions and disrupt_percentage, only the insertions located within the designated 5' portion of the genes are considered to disrupt the sequence. It's same as the output without operon information in **_mapped. 
   INSEQ_experiment.scarf_Samplename.bowtiemap_processed.txt_chromosomename_filter_cpm.txt_genes_sites.txt: The txt file lists all the unique insertions in each gene.
   INSEQ_experiment.scarf_Samplename.bowtiemap_processed.txt_chromosomename_filter_cpm.txt_genes_polar.txt: The txt file lists the genes potentially disrupted by polar effect (that is, genes downstream of an insertion in an operon are counted).
   INSEQ_experiment.scarf_Samplename.bowtiemap_processed.txt_chromosomename_filter_cpm.txt_insertions.txt: The txt file lists all the insertions, and the locations in the genome, including intergenic insertions.

(B) Arrayed library
As described in Box 1 in the main text, when you construct an arrayed library, you will pool these mutant strains into 24 pools, and then sequence these 24 pools with a different barcode assigned to each pool. To match the insertions from the sequencing results to the mutants in the arrayed library, the files you need are:
 (1) The raw reads file
 (2) The sample mapping file, please give the sample name as PoolXX, eg Pool1, Pool2, ..., Pool24. (refer to arrayed_mapping.txt in the directory "Example")
 (3) The arrayed library mapping file Write_dws_X_programs map.txt generated by write_dws.pl when generating the pooling commands for the robot. (refer to Write_dws_27_programs map.txt in the directory "Example")
 
Run the analysis pipeline with the arrayed function turned on by adding "- arrayed" to the command line. For example:
 perl /Users/mwu/INSEQ_analysis/INSeq_pipeline.pl -i AGS_arrayed.txt -m Arrayed_library_mapfile.txt -s BthetaVPI_5482 -arrayed

After all the mappingjob files are finished and all the normalized count reads (.cpm) files are generated, run the perl script arrayed_library.pl in "Arrayed_library" as follows:
   perl arrayed_library.pl <Write_dws_X_programs map.txt > <all processed X_filter_cpm.txt files from the 24 pooled samples>

The first argument is the arrayed library mapping file generated by write_dws.pl when generating the pooling command for the robot, the remaining arguments are the output .cpm files in the "results" directory.

One example command is:
perl arrayed_library.pl Write_dws_20_programs map.txt results/*_filter_cpm.txt 

There are two output files in the "results" directory:
  (1) The arrayed_match.txt with the format
Chromosome      Coordinate      StrainID        Plate   Well    Mismatches

  (2) The mapping_arrayed.log file with the format
<Pool number (sample number)> <number of the wells have mismatches between the mapping file and the actual location>

with the record of how many bad_pools (the number of mismatches between the mapping file and the actual sequencing result in each of 24 pools).

*************************************************************
Description of the scripts and files in the analysis package:
*************************************************************
Typing the command "perl XXXXX.pl" will print out the usage statement for each script.

  (i) "INSeq_pipeline.pl" is the script that will split and sort the raw reads based on the barcodes, find the transposon sequence in the reads, trim, store the remaining chromosome sequences based on their barcodes in folder "bcsortedseqs", and generate several mappingjobsXXXX files. The jobfiles include mapping with bowtie using parameters as "-m 1 -l 17 -n 1 -a --best --strata", which requests only the best unique hit with no more than 1 mismatch will print out into mapping results by bowtie.

  (ii) "process_bowtie_output.pl" is the script that will parse the mapping results, and print out all the insertions with the left and right reads in the format "chromosome name \t Coordinate \t Left Reads Number \t Right Reads number \t Total Reads number \n".

  (iii) "normalize_processed_filter.pl" is the script that will remove insertions with the "Total Reads Number" less than or equal to 3, then normalize the "Total Reads Number" to a million, output the normalized "counts per million" (cpm) for each coordinate in the format "chromosome name \t Coordinate \t Normalized Reads number (cpm) \n".

  (iv) "map_genes.pl" maps the insertion coordinate to genes based on operon data.

  (v) "map_genes_v2.pl" maps the insertion coordinate without operon data.

  (vi) "config.txt" is the txt file to store where the bowtie package is located in the system, and needs to be edited by user when installing the package.

  (vii) Directory "Arrayed_library" includes additional files needed for arrayed library mapping.
  
    (1) "16384_strings.txt" is the list of 24-bit strings with a minimum Hamming distance of 6. This set can be used to map up to 170 96-well trays.

    (2) "13000_strings.txt" is the subset of these that are less likely to be mistaken for each other if multiple archived strains happen to carry transposons at the identical genome coordinate. This set is preferable for arrayed libraries of under 135 trays. 

    (3) "write_dws.pl" is the script that translates the pooling patterns like (vii) (1) and (2) into programs for the EpMotion robot.

    (4) "arrayed_library.pl" is the script that matches the insertions identified by sequencing to the mutants in the arrayed library.

  (viii) Directory "indexes"
  
   The "indexes" folder contains genome reference sequences to map the INSeq reads against. Currently, it contains the Bacteroides_thetaiotaomicron_VPI-5482 genome in the directory "BthetaVPI_5482", which includes the genome sequence formatted by bowtie, ptt file (a file from NCBI that provides gene locations along the genome for protein coding genes), and the operon information for the main chromosome in B.thetaiotaomicron. Users will create their own directory under "index" for their own organism as described in "Getting started" Step (3).

   (ix) Example
   The "Example" folder contains a small set of sample files: one small raw reads file "Sample_rawreads.txt" and mapping file "Sample_mapping_file.txt" and all output files generated by the package as "Results_expected". Users can run the analysis package on this data to check if the package installed properly. There are also two example files for the arrayed library: arrayed_mapping.txt" and "Write_dws_27_programs map.txt".

NOTE: There are also four "demo"-related folders, which can be download as a single .zip file from http://gordonlab.wustl.edu/SuppData.html. "demo" contains the samples of raw reads files and mapping file, "demo_test" is an example of a regular run; "demo_operon" is the example of same data in demo_test but with operon information;  "demo_arrayedlibrary" contains an example of arrayed library mapping.


