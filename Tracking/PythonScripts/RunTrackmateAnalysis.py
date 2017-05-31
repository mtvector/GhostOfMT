# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 18:02:05 2016

@author: mschmitz
"""
import math
import os
import subprocess
import random
import networkx as nx
import pandas
import numpy as np
import csv
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
import multiprocessing as multi
from multiprocessing import Manager
from threading import Timer
manager = Manager()
random.seed(64)


def timedRun(cmd, timeout_sec):
    #print cmd
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
    stderr=subprocess.PIPE)
    kill_proc = lambda p: p.kill()
    timer = Timer(timeout_sec, kill_proc, [proc])
    try:
        timer.start()
        stdout,stderr = proc.communicate()
        return stdout,stderr

    finally:
        timer.cancel()



def dist(a1,a2,b1,b2):
    return np.sqrt((a1-b1)**2 + (a2-b2)**2)

def getNumToRoot(G,n,roots):
    count= None
    for r in roots:
        #print("counting")
        try:
            p=nx.shortest_path(G,r,n)
        except: 
            p=[]
        if len(p)>0:
            count=0
            for pathnode in p:
                count += len(G.successors(pathnode))-len(G.predecessors(pathnode))        
            return count
    return count


def evaluateFitness(params):
    #print([os.path.expanduser("~/FIJI/Fiji.app/ImageJ-linux64"), os.path.expanduser("~/code/cb/trackmater.py"), srcpath] + params)
    #subprocess.call(map(str,[os.path.expanduser("~/FIJI/Fiji.app/ImageJ-linux64"), os.path.expanduser("~/code/cb/trackmater.py"), srcpath] + params))
    out,err = timedRun(map(str,[os.path.expanduser("~/FIJI/Fiji.app/ImageJ-linux64"), os.path.expanduser("~/code/cb/trackmater.py"), srcpath] + params), 6000)
    G,ref = loadNets(srcpath,params) 
    if len(G.nodes()) < 1:
        print("Tracking Failed!")
        testpath  = srcpath+'_' +'_'.join(map(str,params))
        f=open(testpath+'/fail.txt','wb')
        writer=csv.writer(f, delimiter=',')
        writer.writerow([err])
        f=open(testpath+'/out.txt','wb')
        writer=csv.writer(f, delimiter=',')
        writer.writerow([out])
        f.close()
        return 0.0,0.0,0.0,
    else:
        return calcErrors(G,ref)

def mergeNodes(G, n1, n2):
    if G.node[n1]['t'] != G.node[n2]['t']:
        return None
    out = {}    
    for attribute in G.node[n1].keys():
        if attribute == 'track':
            out[attribute] = G.node[n1][attribute]
        else:
            out[attribute] = (G.node[n1][attribute] + G.node[n2][attribute])/2
    ended = False
    outnum = 0
    while not ended:
        outnum = outnum +1
        if outnum not in G.nodes():
            ended=True      
    G.add_node(outnum,out)
    nei1=G.neighbors(n1)
    nei2=G.neighbors(n2)
    neibs = nei1+nei2
    for n in neibs:
        if G.node[n]['t'] > G.node[outnum]['t']:
            G.add_edge(outnum, n)
        if G.node[n]['t'] < G.node[outnum]['t']:
            G.add_edge(n,outnum)
    G.remove_nodes_from([n1,n2])

def loadNets(srcpath, params):
    testpath  = srcpath+'_' +'_'.join(map(str,params))
    print "testpath="+testpath
    ref= nx.DiGraph()
    trainset= pandas.read_csv(srcpath+'tracked/trainset.csv')
    remov = trainset.loc[:,'divs'] > -1
    trainset = trainset.loc[remov,:]
    for i in range(trainset.shape[0]):
        ref.add_node(i,dict(trainset.loc[trainset.index[i],:]))
        if(i>0):
            if(trainset.loc[trainset.index[i],'track'] == trainset.loc[trainset.index[i-1],'track'] ):
                ref.add_edge(i-1,i)
    for i in ref.nodes():
        if i not in ref.nodes():
            continue
        for j in ref.nodes():
            if i != j:
                #print i 
                #print j
                A = ref.node[i]
                B = ref.node[j]
                if A['t'] == B['t']:
                    if dist(A['x'],A['y'],B['x'],B['y']) < 3:
                        mergeNodes(ref, i, j)
                        break 
                        
    edgeFrames={}
    trackFrames={}
    if 1==1:        
    #if  os.path.isfile(os.path.expanduser(srcpath+'/graphtable.csv')):
        for root, directories, filenames in os.walk(testpath):
            for filename in filenames:
                if 'edge' in filename:
                    edgeFrames[int(filename.replace('edge','').replace('.txt',''))] = pandas.read_table(os.path.join(root,filename),sep = '\t')
                    #print pandas.read_csv(os.path.join(root,filename))
                if 'track' in filename:
                    trackFrames[int(filename.replace('track','').replace('.txt',''))] = pandas.read_table(os.path.join(root,filename),sep = '\t')
                
        
        G=nx.DiGraph()
        print "edges"
        for edge in edgeFrames.keys():
            source = [int(s.replace('ID','')) for s in edgeFrames[edge].loc[:,'source']] 
            target = [int(s.replace('ID','')) for s in edgeFrames[edge].loc[:,'target']] 
            G.add_edges_from([(source[x],target[x])for x in range(len(source))])
        print "spots"
        for track in trackFrames.keys():
                spotID =trackFrames[track].loc[:,'spotID']        
                t = trackFrames[track].loc[:,'t']
                fluor = trackFrames[track].loc[:,'fluorescence']
                xc = trackFrames[track].loc[:,'x']
                yc = trackFrames[track].loc[:,'y']
                ed = trackFrames[track].loc[:,'estdiam']
                G.add_nodes_from([(spotID[x],{'fluorescence': fluor[x] , 't' : t[x], 'track':track, 'x':xc[x] ,'y':yc[x],'estdiam':ed[x]} )for x in range(len(spotID))])
        
    
        roots = []
        for n in G.nodes():
            if len(G.predecessors(n)) == 0: #and nx.get_node_attributes(G,'t')[n] < 50:
                roots.append(n)
        print len(G.nodes())
        print("screening rootless nodes")
        for n in G.nodes():
            div = getNumToRoot(G,n,roots)
            if div == None:
                G.remove_node(n)
            else:
                G.node[n]['divs'] = div
            
    return G, ref    
  
def getDistOfSpots(x,y,t,net):
    which = lambda lst:list(np.where(lst)[0])
    sameframe = which([q == t for q in nx.get_node_attributes(net,'t').values()])
    xs=nx.get_node_attributes(net,'x').values()
    ys=nx.get_node_attributes(net,'y').values()
    if len(sameframe)>0:
        inds = np.array(net.nodes())[sameframe]
        xs=np.array(xs)[sameframe]
        ys=np.array(ys)[sameframe]
        return {inds[i] : dist(x,y ,xs[i],ys[i]) for i in range(len(xs))}
    else:
        return None
 
 
def calcErrors(G,ref):
    spotscore = 0.0
    totspotscore = .00001
    divscore = 0.0
    totdivscore = .00001
    connexscore = 0.0
    totconnexscore = .00001
    for s in set(nx.get_node_attributes(ref,'track')):
        which = lambda lst:list(np.where(lst)[0])
        trackinds = which([q == s for q in nx.get_node_attributes(ref,'track').values()])
        trackroots = []
        for n in np.array(ref.nodes())[trackinds]:
            if len(ref.predecessors(n)) == 0:
                trackroots.append(n)

        thisroot=[]
        for i in trackroots:
            x=ref.node[i]['x']
            y=ref.node[i]['y']
            t=ref.node[i]['t']-1
            rootframe = which([q == t for q in nx.get_node_attributes(G,'t').values()])
            d = getDistOfSpots(x,y,t,G)
            #get inds of sorted values
            for j in np.array(G.nodes())[rootframe]:
                if d[j] < radiuslim:
                    thisroot.append(j)

        connex = []
        for i in np.array(ref.nodes())[trackinds]:
            totdivscore =totdivscore + 1
            totspotscore = totspotscore +1
            x=ref.node[i]['x']
            y=ref.node[i]['y']
            t=ref.node[i]['t']-1
            d = getDistOfSpots(x,y,t,G)
            #print d
            connex.append([i,None,None])
            if d != None:
                #get inds of sorted values
                inds = sorted(d, key=d.get)
                #print inds
                if(d[inds[0]]<radiuslim):
                    spotscore = spotscore+1
                toroot=[]
                for n in inds[0:2]:
                    toroot.append(getNumToRoot(G,n,thisroot))
                if toroot[0] != None or len(toroot)>1:
                    connex[len(connex)-1][1] = inds[0]
                    if len(toroot)>1:
                        connex[len(connex)-1][2] = inds[1]
                    if toroot[0] != None:
                        divscore =divscore + (1/(1 + abs(toroot[0] - ref.node[i]['divs']))**2)
                    elif toroot[1] != None:
                        divscore= divscore + (1/(1 + abs(toroot[1] - ref.node[i]['divs']))**2)
        GU=G.to_undirected() 
        #print connex
        for oo in range(len(connex)-1):
            for pp in range(len(connex)-1):
                if oo != pp:
                    totconnexscore=totconnexscore+1
                    try: 
                        nx.bidirectional_dijkstra(GU,G.nodes()[connex[oo][1]],G.nodes()[connex[pp][1]])
                        connexscore=connexscore+1   
                        continue
                    except:
                        pass
                    try:
                        nx.bidirectional_dijkstra(GU,G.nodes()[connex[oo][1]],G.nodes()[connex[pp][2]])
                        connexscore=connexscore+1 
                        continue
                    except:
                        pass    
                        
    print(str(divscore/totdivscore)+" "+str(spotscore/totspotscore)+" "+str(connexscore/totconnexscore))
    return divscore/totdivscore, spotscore/totspotscore ,connexscore/totconnexscore,
    


def mutAddInt(individual, lo, hi, indpb):
    for i in range(len(individual)-2):
        individual[i] += random.randint(lo, hi)
        if individual[i]<0:
            individual[i] = 0
        if individual[0]<5:
            individual[0] = 5
        individual[len(individual)-1]=random.choice([True,False])
    return individual,

srcpath = os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x-04572y001856-10x-FLsorted")
#Q,ref=loadNets(srcpath,['tracked'])

#G,E=loadNets(srcpath,[6,5,2,20,20,5,30,True])


#srcpath = os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x-04572y001856-10x-FLsorted_7_7_4_14_9_10_3_True")
'''
Params to trackmater
radius = sys.argv[2]
threshold = sys.argv[3]
maxGap=sys.argv[4]
maxLinkDist=sys.argv[5]
maxGapDist = sys.argv[6]
trackDuration=sys.argv[7]
qualityMin=sys.argv[8]
'''
guess = [8,7,3,20,20,3,20,True]
radiuslim=10

creator.create("FitnessMax", base.Fitness, weights=(1.0,1.0,1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

def initIndividual(icls, content):
    return icls(content)

def initPopulation(pcls, ind_init,num_ind):
    #contents = [[8,7,3,16,10,11,5]for x in range(num_ind)]
    contents = [[x + random.randint(-3,5) for x in guess[0:(len(guess)-1)]] + [random.choice([True,False])] for x in range(num_ind) ]
    return pcls(ind_init(c) for c in contents)

toolbox = base.Toolbox()
toolbox.register("individual", initIndividual, creator.Individual)
toolbox.register("population", initPopulation, list, toolbox.individual)

toolbox.register("evaluate", evaluateFitness)
toolbox.register("mate", tools.cxOnePoint)
toolbox.register("mutate", mutAddInt,lo= -2, hi = 2, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)
NPOP, CXPB, MUTPB, NGEN = 60 ,0.2, 0.3, 20
print("Start of evolution")
pop = toolbox.population(NPOP)
# Evaluate the entire population
'''
hof = tools.HallOfFame(8)
stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", np.mean)
stats.register("std", np.std)
stats.register("min", np.min)
stats.register("max", np.max)
    
pop, log = algorithms.eaSimple(pop, toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=NGEN, 
                   stats=stats, halloffame=hof, verbose=True)


'''
p=multi.Pool(processes=9)
fitnesses = p.map(toolbox.evaluate, pop)

#fitnesses = list(map(toolbox.evaluate, pop))
for ind, fit in zip(pop, fitnesses):
    ind.fitness.values = fit

print("  Evaluated %i individuals" % len(pop))

# Begin the evolution
for g in range(NGEN):
    print("-- Generation %i --" % g)
    
    # Select the next generation individuals
    offspring = toolbox.select(pop, len(pop))
    # Clone the selected individuals
    offspring = list(map(toolbox.clone, offspring))

    # Apply crossover and mutation on the offspring
    for child1, child2 in zip(offspring[::2], offspring[1::2]):

        # cross two individuals with probability CXPB
        if random.random() < CXPB:
            toolbox.mate(child1, child2)

            # fitness values of the children
            # must be recalculated later
            del child1.fitness.values
            del child2.fitness.values

    for mutant in offspring:

        # mutate an individual with probability MUTPB
        if random.random() < MUTPB:
            toolbox.mutate(mutant)
            del mutant.fitness.values

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        
    p=multi.Pool(processes=8)
    #p.map(func, pop)
    fitnesses = p.map(toolbox.evaluate, pop)
    #fitnesses = map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit
    
    print("  Evaluated %i individuals" % len(invalid_ind))
    
    # The population is entirely replaced by the offspring
    pop[:] = offspring
    
    # Gather all the fitnesses in one list and print the stats
    fits = [ind.fitness.values[0] for ind in pop]
    
    length = len(pop)
    mean = sum(fits) / length
    sum2 = sum(x*x for x in fits)
    std = abs(sum2 / length - mean**2)**0.5
    
    print(["  Min %s" % min(fits)])
    print(["  Max %s" % max(fits)])
    print(["  Avg %s" % mean])
    print(["  Std %s" % std])
    best_ind = tools.selBest(pop, 5)
    csvfile = open(os.path.expanduser('~/Desktop/')+'paramtest.csv', 'a')
    spamwriter = csv.writer(csvfile, delimiter=',')
    spamwriter.writerow(["  %i Generation" % g])
    spamwriter.writerow(["  Evaluated %i individuals" % len(invalid_ind)])
    spamwriter.writerow(["  Min %s" % min(fits)])
    spamwriter.writerow(["  Max %s" % max(fits)])
    spamwriter.writerow(["  Avg %s" % mean])
    spamwriter.writerow(["  Std %s" % std])
    for bi in best_ind:
        print bi, bi.fitness.values
        spamwriter.writerow( bi + list(bi.fitness.values))
    csvfile.close()
print("-- End of (successful) evolution --")
best_ind = tools.selBest(pop, 1)[0]
print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))
