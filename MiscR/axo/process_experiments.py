import subprocess


subprocess.call(["~/lib/Merlin/merlin", "-d", "~/code/data/axo/ExpressionSet/AllAmby002.txt", "-o", "~/code/data/axo/AllAmby002Net", "-l" ,"~/code/data/axo/ExpressionSet/RegulatorProbes.txt", "&" ])    
subprocess.call(["~/lib/Merlin/merlin", "-d", "~/code/data/axo/ExpressionSet/normAllAmby002.txt", "-o", "~/code/data/axo/normAllAmby002Net", "-l" ,"~/code/data/axo/ExpressionSet/RegulatorProbes.txt", "&" ])    

#~/lib/Merlin/merlin -d ~/code/data/axo/ExpressionSet/AllAmby002.txt -o ~/code/data/axo/AllAmby002Net -l ~/code/data/axo/ExpressionSet/RegulatorProbes.txt &


"""~/lib/merlin/learnMERLIN -m example/in/yeast_expr_interpolated -o example/out/ -r4 -p-5 -h 0.6 -l example/in/regulators.txt -v1 -c example/in/clusterassign.txt -k300
 
There are several arguments merlin takes, however, the most important ones are r, p, h, v, c and l described below. We are working on improving the documentation and usability of merlin :
-m model
-o outputdir
-k maxfactorsize
-v cross_validation_cnt
-l restricted_regulator_fname
-p sparsity_prior
-r modularity_prior
-c initial_cluster_assignments"""

"""d [required]
Text file with the expression data using the DREAM Challenge format (see File Formats)
l: Text file with the list of regulators (see File Formats). Default: all genes are potential regulators.
c: Text file that specifies the initial module assignments. Default: performs a random partitioning of genes into max[squareroot(n/2),30] clusters, where n is the number of genes in the data. That is the default option will have no more than 30 initial clusters.
h: hierarchical clustering threshold (default 0.6)
p: parameter for sparsity (default -5)
r: parameter for module prior (default 4)
k: max number of regulators that a gene can have (default 300)
o [required] : name of output directory that must exist before running the program
v [required] : number of folds for cross validation (default 1)"""