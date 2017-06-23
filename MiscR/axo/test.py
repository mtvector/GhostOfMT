import GENIE3
from numpy import loadtxt
import sys
exprsfile=sys.argv[1]
data = loadtxt(exprsfile,skiprows=1)
infile = open(exprsfile)
gene_names = infile.readline()
infile.close()
gene_names = gene_names.rstrip('\n').split('\t')
#print gene_names

regulators = gene_names
#regulators = ['CD19', 'CDH17','RAD51','OSR2','TBX3']
VIM = genie3(data,gene_names=gene_names,regulators=regulators)
get_link_list(VIM)
#print VIM
