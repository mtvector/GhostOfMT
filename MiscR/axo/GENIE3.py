from sklearn.tree.tree import BaseDecisionTree
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import ExtraTreesRegressor
from numpy import *
import time
from operator import itemgetter
from multiprocessing import Process
from multiprocessing import Queue

def compute_feature_importances(estimator):
    if isinstance(estimator, BaseDecisionTree):
        return estimator.tree_.compute_feature_importances(normalize=False)
    else:
        importances = [e.tree_.compute_feature_importances(normalize=False)
                for e in estimator.estimators_]
        importances = asarray(importances)
        return sum(importances,axis=0) / len(estimator)



def get_link_list(VIM,gene_names=None,regulators='all',maxcount='all',file_name=None):

    """Gets the ranked list of (directed) regulatory links.

    Parameters
    ----------

    VIM: numpy array
        Array as returned by the function genie3(), in which the element (i,j) is the score of the edge directed from the i-th gene to the j-th gene. 

    gene_names: list of strings, optional
        List of length p, where p is the number of rows/columns in VIM, containing the names of the genes. The i-th item of gene_names must correspond to the i-th row/column of VIM. When the gene names are not provided, the i-th gene is named Gi.
        default: None

    regulators: list of strings, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names), and the returned list contains only edges directed from the candidate regulators. When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'

    maxcount: 'all' or positive integer, optional
        Writes only the first maxcount regulatory links of the ranked list. When maxcount is set to 'all', all the regulatory links are written.
        default: 'all'

    file_name: string, optional
        Writes the ranked list of regulatory links to the file file_name.
        default: None




    Returns
    -------

    The list of regulatory links, ordered according to the edge score. Auto-regulations do not appear in the list. Regulatory links with a score equal to zero are randomly permuted. In the ranked list of edges, each line has format:

        regulator   target gene     score of edge



    """

    # Check input arguments      
    if not isinstance(VIM,ndarray):
        raise ValueError('VIM must be a square array')
    elif VIM.shape[0] != VIM.shape[1]:
        raise ValueError('VIM must be a square array')

    ngenes = VIM.shape[0]

    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')

    if regulators is not 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')

    if maxcount is not 'all' and not isinstance(maxcount,int):
        raise ValueError('input argument maxcount must be "all" or a positive integer')

    if file_name is not None and not isinstance(file_name,str):
        raise ValueError('input argument file_name must be a string')



    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = range(ngenes)
    else:
        input_idx = []
        for i in range(ngenes):
            if gene_names[i] in regulators:
                input_idx.append(i)

    nTFs = len(input_idx)

    # Get the non-ranked list of regulatory links
    vInter = []
    for i in input_idx:
        for j in range(ngenes):
            if i != j:
                vInter.append((i,j,VIM[i,j]))

    # Rank the list according to the weights of the edges        
    vInter_sort = sorted(vInter,key=itemgetter(2),reverse=True)
    nInter = len(vInter_sort)

    # Random permutation of edges with score equal to 0
    flag = 1
    i = 0
    while flag and i < nInter:
        (TF_idx,target_idx,score) = vInter_sort[i]
        if score == 0:
            flag = 0
        else:
            i += 1

    if not flag:
        items_perm = vInter_sort[i:]
        items_perm = random.permutation(items_perm)
        vInter_sort[i:] = items_perm

    # Write the ranked list of edges
    nToWrite = nInter
    if isinstance(maxcount,int) and maxcount >= 0 and maxcount < nInter:
        nToWrite = maxcount

    if file_name:

        outfile = open(file_name,'w')

        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('%s\t%s\t%.6f\n' % (gene_names[TF_idx],gene_names[target_idx],score))
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('G%d\tG%d\t%.6f\n' % (TF_idx+1,target_idx+1,score))


        outfile.close()

    else:

        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print '%s\t%s\t%.6f' % (gene_names[TF_idx],gene_names[target_idx],score)
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print 'G%d\tG%d\t%.6f' % (TF_idx+1,target_idx+1,score)



def genie3(expr_data,gene_names=None,regulators='all',tree_method='RF',K='sqrt',ntrees=1000):

    '''Computation of tree-based scores for all putative regulatory links.

    Parameters
    ----------

    expr_data: numpy array
        Array containing gene expression values. Each row corresponds to a condition and each column corresponds to a gene.

    gene_names: list of strings, optional
        List of length p, where p is the number of columns in expr_data, containing the names of the genes. The i-th item of gene_names must correspond to the i-th column of expr_data.
        default: None

    regulators: list of string, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names). When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'

    tree-method: 'RF' or 'ET', optional
        Specifies which tree-based procedure is used: either Random Forest ('RF') or Extra-Trees ('ET')
        default: 'RF'

    K: 'sqrt', 'all' or a positive integer, optional
        Specifies the number of selected attributes at each node of one tree: either the square root of the number of candidate regulators ('sqrt'), the number of candidate regulators ('all'), or any positive integer.
        default: 'sqrt'

    ntrees: positive integer, optional
        Specifies the number of trees grown in an ensemble.
        default: 1000


    Returns
    -------

    An array in which the element (i,j) is the score of the edge directed from the i-th gene to the j-th gene. All diagonal elements are set to zero (auto-regulations are not considered). When a list of candidate regulators is provided, all the edges directed from a gene that is not a candidate regulator are set to zero.

    '''

    time_start = time.time()

    # Check input arguments
    if not isinstance(expr_data,ndarray):
        raise ValueError('expr_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')

    ngenes = expr_data.shape[1]

    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')

    if regulators is not 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')        

    if tree_method is not 'RF' and tree_method is not 'ET':
        raise ValueError('input argument tree_method must be "RF" (Random Forests) or "ET" (Extra-Trees)')

    if K is not 'sqrt' and K is not 'all' and not isinstance(K,int): 
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')

    if isinstance(K,int) and K <= 0:
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')

    if not isinstance(ntrees,int):
        raise ValueError('input argument ntrees must be a stricly positive integer')
    elif ntrees <= 0:
        raise ValueError('input argument ntrees must be a stricly positive integer')


    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = range(ngenes)
    else:
        input_idx = []
        for i in range(ngenes):
            if gene_names[i] in regulators:
                input_idx.append(i)


    # Learn an ensemble of trees for each target gene, and compute scores for candidate regulators
    VIM = zeros((ngenes,ngenes))

    proc_array = []
    number_queue = Queue()
    result_queue = Queue()
    for i in range(ngenes):

        print 'Gene %d...' % (i+1)

        p = Process(target=genie3_single, args=(expr_data,i,input_idx,tree_method,K,ntrees,i,number_queue,result_queue))
        p.start()
        proc_array.append(p)

    for i in range(ngenes):
        VIM[number_queue.get(),:] = result_queue.get()

    for p in proc_array:
        p.join()

    VIM = transpose(VIM)

    time_end = time.time()
    print "Elapsed time: %.2f seconds" % (time_end - time_start)

    return VIM

def genie3_single(expr_data,output_idx,input_idx,tree_method,K,ntrees,number,number_queue,result_queue):

    ngenes = expr_data.shape[1]

    # Normalize output data
    output = expr_data[:,output_idx]
    output = output / std(output)

    # Remove target gene from candidate regulators
    input_idx = input_idx[:]
    if output_idx in input_idx:
        input_idx.remove(output_idx)

    expr_data_input = expr_data[:,input_idx]

    # Parameters of the tree-based method
    if tree_method == 'RF':
        print_method = 'Random Forests'
        if K == 'sqrt':
            print_K = round(sqrt(len(input_idx)))
            treeEstimator = RandomForestRegressor(n_estimators=ntrees,max_features="sqrt")
        elif K == 'all':
            print_K = len(input_idx)
            treeEstimator = RandomForestRegressor(n_estimators=ntrees,max_features="auto")
        else:
            if K < len(input_idx):
                print_K = K
                treeEstimator = RandomForestRegressor(n_estimators=ntrees,max_features=K)
            else:
                print_K = len(input_idx)
                treeEstimator = RandomForestRegressor(n_estimators=ntrees,max_features="auto")

    elif tree_method == 'ET':
        print_method = 'Extra Trees'
        if K == 'sqrt':
            print_K = round(sqrt(len(input_idx)))
            treeEstimator = ExtraTreesRegressor(n_estimators=ntrees,max_features="sqrt")
        elif K == 'all':
            print_K = len(input_idx)
            treeEstimator = ExtraTreesRegressor(n_estimators=ntrees,max_features="auto")
        else:
            if K < len(input_idx):
                print_K = K
                treeEstimator = ExtraTreesRegressor(n_estimators=ntrees,max_features=K)
            else:
                print_K = len(input_idx)
                treeEstimator = ExtraTreesRegressor(n_estimators=ntrees,max_features="auto")

    print "Tree method = %s, K = %d, %d trees" % (print_method,print_K,ntrees)

    # Learn ensemble of trees
    treeEstimator.fit(expr_data_input,output)

    # Compute importance scores
    feature_importances = compute_feature_importances(treeEstimator)
    vi = zeros(ngenes)
    vi[input_idx] = feature_importances

    number_queue.put(number)
    result_queue.put(vi)

    return vi
