# tnseq_function
Association between sequence homology and functional similarity

## code directory
This directory contains scripts for processing the datatables from JGI and output from panX.

## fitness_tables directory
This directory contains JGI generated datasets on the tnseq fitness values in different host species and environments. When conducting experiments, JGI gave each selection experiment a set number. In each fitness dataset, a condition is referred to both by its set number and a short description of the condition. This description is not detailed, however. The full description can be found in the file all_conditions.txt.

#### all_conditions.txt:
Full description of conditions in which tnseq experiment was conducted

#### full_tag_dict_file.cpk:
This pickled dictionary has they key of (experiment--set number, Gene Cluster ) and value of the fitness value of the gene cluster in a given experiment. With this file, one call pull the fitness value for a given gene cluster across every experiment in which that gene cluster was present. 

#### SNP_whole_matrix.aln
This is the output from the core genome phylogeny of panX run on gammaproteobacteria.
