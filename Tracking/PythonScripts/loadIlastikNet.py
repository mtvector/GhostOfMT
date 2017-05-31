# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 11:46:56 2016

@author: mschmitz
"""

import math
import os
import random
import networkx as nx
import pandas
import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import colorsys

random.seed(64)

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


def dist(a1,a2,b1,b2):
    return np.sqrt((a1-b1)**2 + (a2-b2)**2)


def loadIlastikNet(srcpath,fn):
    testpath  = srcpath+fn
    print "testpath="+testpath
    G = nx.DiGraph()
    which = lambda lst:list(np.where(lst)[0])
    netframe= pandas.read_csv(testpath+"_table.csv")
    divframe= pandas.read_csv(testpath+"_divisions.csv")
    print divframe
    #netframe = netframe.loc[range(300),:]
    print netframe.columns.values
    print netframe.loc[:,np.array(['object_id','track_id1','timestep','RegionCenter_0','RegionCenter_1'])]
    nameflip = {'object_id':'id','track_id1':'track','timestep':'t','RegionCenter_0':'x','RegionCenter_1':'y'}
    netframe.rename(columns=nameflip, inplace=True)
    for i in range(netframe.shape[0]):
        if netframe.loc[netframe.index[i],'track'] !=0:
            #print netframe.index[i]
            G.add_node(netframe.index[i],dict(netframe.loc[netframe.index[i],:]))
    #print which([q == 12 for q in nx.get_node_attributes(G,'track_id1').values()])
    for s in set(nx.get_node_attributes(G,'track').values()):
        #print 'Track '+str(s)
        trackinds = which([q == s for q in nx.get_node_attributes(G,'track').values()])
        ti = np.array(G.nodes())[trackinds] 
        for n in range(len(ti)-1):   
            G.add_edge(ti[n],ti[n+1])
            #print G.node[ti[n]]['timestep']
            #print G.node[ti[n+1]]['timestep']
    for i in divframe.index:
        G.add_edge(divframe.loc[i,'parent_oid'], divframe.loc[i,'child1_oid'])
        G.add_edge(divframe.loc[i,'parent_oid'], divframe.loc[i,'child2_oid'])
    return G

def loadTrainNet(srcpath):
    ref= nx.DiGraph()
    trainset= pandas.read_csv(srcpath+'tracked/trainsetAdj.csv')
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
    return ref

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
    for s in set(nx.get_node_attributes(ref,'track').values()):
        #print "Track " + str(s)
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
        print thisroot
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
                    print toroot
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

srcpath=os.path.expanduser("~/Desktop/")
fn= "-04572_001856_10_16051_____________-exported_data"

radiuslim=15
srcpath = os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x-04572y001856-10x-FLsorted")
ref = loadTrainNet(srcpath)
srcpath = os.path.expanduser("~/Desktop/")
G=loadIlastikNet(srcpath,fn)
nx.write_edgelist(G, os.path.expanduser('~/Desktop/edgelist.txt'))
nx.write_edgelist(ref, os.path.expanduser('~/Desktop/refedgelist.txt'))

#print calcErrors(G,ref)
outDF = pandas.DataFrame(columns=G.node[G.nodes()[1]].keys())
interact = False
fn = os.path.basename(srcpath)

for n in G.nodes():
    outDF.loc[n] = G.node[n].values()
#f=open( os.path.expanduser(srcpath+'/graphtable.csv'),'wb')
#outDF.to_csv(f, mode='w',sep='\t')
#f.close()
print outDF

rgb = [colorsys.hsv_to_rgb(i ,1.,1.) for i in np.array(nx.get_node_attributes(G,'divs').values())/ (float(max(nx.get_node_attributes(G,'divs').values()))*1.3)]
pp=PdfPages(os.path.expanduser('~/Desktop/'+fn+'analysis.pdf'))
'''
plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'adjFluorescence'],c=rgb)#outDF.loc[:,'divs'])
plt.xlabel('time')
plt.ylabel('fluorescence')
#plt.show()
#plt.draw()
pp.savefig()
plt.close()
positions={x : np.array([G.node[x]['t'], G.node[x]['adjFluorescence']]) for x in G.nodes()}
nx.draw_networkx(G,pos=positions, node_size =15,edge_color='black',linewidths= .1, width = .3,node_color = rgb, with_labels=False)
pp.savefig()
plt.close()
plt.scatter(outDF.loc[:,'divs'] ,outDF.loc[:,'adjFluorescence'],c=rgb)#outDF.loc[:,'divs'])
plt.xlabel('divisions')
plt.ylabel('fluorescence')
#plt.show()
#plt.draw()
pp.savefig()
plt.close()
'''
plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'divs'],c=rgb)
plt.xlabel('time')
plt.ylabel('divisions')
#plt.show()
#plt.draw()
pp.savefig()
plt.close()
pp.close()