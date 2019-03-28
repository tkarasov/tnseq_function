#!/usr/bin/env python
'''the goal of this script is to take the pan genome for a given bacterial group, and to assign diversity statistics to each ortholog group'''

import sys, os, pickle, ete3, pandas
from Bio import SeqIO

#now build object for every gene in tnseq genome

class tnseq_gene:

	def __init__(self, tnid, pgid, pggc, pg_num, ttree, gene_ev, gene_div):
		self.tnid = tnid
		self.pgid = pgid
		self.pggc = pggc
		self.pg_num = pg_num
		#self.pident = pident
		#self.len_gene = len_gene
		#self.gc_content = gc_content
		#self.gc_bias =
		self.ttree = ttree
		self.gene_ev = gene_ev
		self.gene_div = gene_div

def getGC(seq):
    return(str(int(round((sum([1.0 for nucl in seq if nucl in ['G', 'C']]) / len(seq)) * 100))))

def fit_object(file_name):
    fit = pandas.read_table(file_name, index_col =2)
    #Time spent in phylogenetic tree
    tnseq_stats = {}
    for key,value in tnseq_panx.iteritems():
        try:
            tnid = key.split(":")[1]
        except IndexError:
            print("ERROR")
        tnid=tnseq_panx[key][0].split("_")[2]
        pgid = value[0]
        pggc = value[1]
        pgnum = value[2]
        pident = gene_events[pgnum]
        ttree = ts_tree[pgnum]
        gene_div = genediv[pggc]
        gene_ev = gene_events[pgnum]
        tnseq_stats[tnid] = tnseq_gene(tnid, pgid, pggc, pgnum, ttree, gene_ev, gene_div)
    fit.index=[str(rec) for rec in fit.index]
    fit['ttree'] = [tnseq_stats[str(rec)].ttree for rec in fit.index]
    fit['gene_ev'] = [tnseq_stats[rec].gene_ev for rec in fit.index]
    fit['gene_div'] = [tnseq_stats[rec].gene_div for rec in fit.index]
    with open("/ebio/abt6_projects9/tnseq/tnseq_function/data/"+file_name+"_tnseq_panx_fitness_data.cpk", 'wb') as handle:
        pickle.dump(fit, handle)

#path_to_pangenome_dir = '/ebio/ag-neher/share/users/wding/panX-refseq/data/Pseudomonadales'#sys.argv[1]
path_to_pangenome_dir='/ebio/abt6_projects9/tnseq/wdata/panX-tnseq/data/Gammaproteobacteria/'
clusters = path_to_pangenome_dir + "/allclusters_final.tsv"
#tnseq_genome = list(SeqIO.parse(os.popen("ls | grep _cds_from_genomic.fna").read().strip(), "fasta"))

#tnseq_panx = pickle.load(open("tnseq_panx_mapped_dict.cpk", 'rb')) '''I think this is the deprecated mapping'''
tnseq_panx = [line.strip().split() for line in open("/ebio/abt6_projects9/tnseq/tnseq_function/fitness_tables/gammaproteobacteria_mappings.txt").readlines()]



#Nucleotide diversity of gene
genediv = pickle.load(open(path_to_pangenome_dir+"/geneClusterx/gene_diversity.cpk", 'rb'))

#GC content

#Deletion/Duplication events pickle.load(open(path_to_pangenome_dir+"/geneCluster/dt_geneEvents.cpk"))
gene_events = pickle.load(open(path_to_pangenome_dir + "/geneClusterx/dt_geneEvents.cpk", 'rb'))

#time_spent in tree
ts_tree = pickle.load(open("/ebio/abt6_projects9/tnseq/tnseq_function/data/branch_gene.cpk", 'rb'))

#fitness info





#fit = pandas.read_table("fit_organism_pseudo1_N1B4.txt", index_col =2)
fit = pandas.read_table("NZ_CP009273_fit_organism_Keio.txt", index_col =2)

#Time spent in phylogenetic tree
tnseq_stats = {}
for key,value in tnseq_panx.iteritems():
	try:
		tnid = key.split(":")[1]
	except IndexError:
		print("ERROR")
		tnid=tnseq_panx[key][0].split("_")[2]
	pgid = value[0]
	pggc = value[1]
	pgnum = value[2]
	pident = gene_events[pgnum]
	ttree = ts_tree[pgnum]
	gene_div = genediv[pggc]
	gene_ev = gene_events[pgnum]
	tnseq_stats[tnid] = tnseq_gene(tnid, pgid, pggc, pgnum, ttree, gene_ev, gene_div)

fit.index=[str(rec) for rec in fit.index]
fit['ttree'] = [tnseq_stats[str(rec)].ttree for rec in fit.index]
fit['gene_ev'] = [tnseq_stats[rec].gene_ev for rec in fit.index]
fit['gene_div'] = [tnseq_stats[rec].gene_div for rec in fit.index]


with open("/ebio/abt6_projects9/tnseq/tnseq_function/data/Keio_tnseq_panx_fitness_data.cpk", 'wb') as handle:
	pickle.dump(fit, handle)

