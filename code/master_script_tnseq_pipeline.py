#!/usr/bin/env python
#first script maps between tnseq and panX results
python map_tnseq_pangenome.py

#performs anc reconstruction to infer events
python pan_genome_anc_reconstruction_step3.py

#creates tnseq gene object with a ton of attributes
python pan_genome_stats_step2.py

#statistical analysis with fitness data
python tnseq_fitness_class_pan_genome_stats_step3.py

#co-occurence analysis
python co_occur_step4.py