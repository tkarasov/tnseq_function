# tnseq_function
Association between sequence homology and functional similarity

## code directory
This directory contains scripts for processing the datatables from JGI and output from panX.

#### compare_gammproteobacteria.py
This script takes in the fitness tables, and the geneCluster.json output by panX. It outputs the gammaproteobacteria_mappings.txt (seen in fitness_tables directory) which maps gene_id in panx to all correpsonding orthologs. The output from this script connects the fitness data with the panX data. 
The script then builds the full_tag_dict_file.cpk dictionary (takes too long). Finally, **the script produces the plots that compare sequence divergence to correlation in fitness values between orthologs**.

#### pan_genome_stats_step2.py
Build fitness object which is then output into tnseq_panx_fitness_data.cpk that has the diversity statistics for an orthology cluster

#### pan_genome_anc_reconstruction_step3.py
Script that uses tidbits of code from panX to calculate the number of times genes were gained and lost along the tree

#### deprecated scripts
map_tnseq_pangenome.py: trying to map gene ids between JGI fitness tables and genomes via blast
tnseq_fitness_class_pan_genome_stats_step3.py: 
master_script_tnseq_pipeline.py: Old pipeline for processing bad mappings
co_occur_step4.py: Looking at correlated co-occurence data. This analysis was essentially scooped in the Nature manuscript.


## fitness_tables directory
This directory contains JGI generated datasets on the tnseq fitness values in different host species and environments. When conducting experiments, JGI gave each selection experiment a set number. In each fitness dataset, a condition is referred to both by its set number and a short description of the condition. This description is not detailed, however. The full description can be found in the file all_conditions.txt.

#### all_conditions.txt:
Full description of conditions in which tnseq experiment was conducted

#### full_tag_dict_file.cpk:
In this pickled list of lists the first element of a row is a tuple (set number (which is the experiment), Gene Cluster ) and the second element is the fitness value of the gene cluster in a given experiment. With this file, one call pull the fitness value for a given gene cluster across every experiment in which that gene cluster was present. 

#### gamma_proteobacteria_mappings.txt
This file identifies the gene ids for the genes assigned to a gene cluster in panX. This was a central issue in analyzing the data. The gene ids in the JGI datasets were different from the gene ids in the full genomes. I must look up which script I used to make these mappings.

#### SNP_whole_matrix.aln
This is the output from the core genome phylogeny of panX run on gammaproteobacteria.

#### distance_matrix_across_strains.txt
Pairwise sequence distance matrix between strains in gammaproteobacteria (calculated from panX output--I must check if renormalized for genome size or instead calculated just on SNP sites)

## Other Relevant data on burrito
path_to_pangenome_dir='/ebio/abt6_projects9/tnseq/wdata/panX-tnseq/data/Gammaproteobacteria/'
path_to_fitness_tables='/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables'
path_to_faa_dir='/ebio/abt6_projects9/tnseq/data/fitness_datasets/tnseq_ref_genomes'
clusters = path_to_pangenome_dir + "/allclusters_final.tsv"
