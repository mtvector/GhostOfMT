#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:36:22 2016

@author: mschmitz

Usage: python GraphReferenceColorReporter.py continuefile1 continuefile2...
"""
#matplotlib.use('WxAgg')
import os
import numpy as np
#import PIL
import pandas
import networkx as nx
#import shutil
import matplotlib.pyplot as plt
#import wx
import math
import re
from matplotlib.backends.backend_pdf import PdfPages
import colorsys
#import seaborn as sns
import sys
import colour
import csv

#radius = sys.argv[0]
spaths = sys.argv[1:]
#spaths = [os.path.expanduser("~/Desktop/Tracking/TestImg_ContinueFile.csv")]
srcpath=spaths[0]
if len(spaths)>1:
    srcpath = srcpath + "Aggregated"

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

def dist(a1,a2,b1,b2):
    return np.sqrt((a1-b1)**2 + (a2-b2)**2)

#G is directed, Q is undirected version of G
def mergeNodes(G,Q, n1, n2):
    if G.node[n1]['t'] != G.node[n1]['t']:
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
    Q = G.to_undirected()
    G.add_node(outnum,out)
    nei1=Q.neighbors(n1)
    nei2=Q.neighbors(n2)
    neibs = nei1+nei2
    for n in neibs:
        if G.node[n]['t'] > G.node[outnum]['t']:
            G.add_edge(outnum, n)
        if G.node[n]['t'] < G.node[outnum]['t']:
            G.add_edge(n,outnum)
    G.remove_nodes_from([n1,n2])

#other option for above function (DFS)
def getNumToRootRecursive(G,n):
    c = 0
    return getNumToRootHelper(G,n,c)
    
def getNumToRootHelper(G,n,c):
    preds= G.predecessors(n)
    if len(preds) < 1:
        return c
    for x in preds:
        c += len(G.successors(x))-len(G.predecessors(x))
        return getNumToRootHelper(G,x,c)        

def loadNets(trainset,radius):
    ref= nx.DiGraph()
    remov = trainset.loc[:,'divs'] > -2
    trainset=trainset.fillna(0)
    trainset = trainset.loc[remov,:]
    for i in range(trainset.shape[0]):
        ref.add_node(i,dict(trainset.loc[trainset.index[i],:]))
        if(i>0):
            if(trainset.loc[trainset.index[i],'track'] == trainset.loc[trainset.index[i-1],'track'] ):
                ref.add_edge(int(i-1),int(i))   
    unD = ref.to_undirected()
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
                    if dist(A['x'],A['y'],B['x'],B['y']) <= radius:
                        print i , j                        
                        mergeNodes(ref,unD, i, j)
                        break
    return ref

channelnames=[]
nets=[]
for s in spaths:
    s = os.path.expanduser(s)
    csvfile = open(s, 'rU')
    reader = csv.reader(csvfile, delimiter=',')
    diameter=float(reader.next()[1])
    trainset= pandas.read_csv(reader.next()[0],header=0)
    reader.next()
    channelnames = [os.path.basename(os.path.normpath(x[1])) for x in reader]
    csvfile.close()
    nets.append(loadNets(trainset,diameter/2))

ref=nx.DiGraph()
for n in range(len(nets)):
    if len(nx.get_node_attributes(ref,'track').values())>0:
        mx=max(nx.get_node_attributes(ref,'track').values())
        nx.set_node_attributes(nets[n],'track', dict( zip(nx.get_node_attributes(nets[n],'track').keys(), [x+mx+1 for x in nx.get_node_attributes(nets[n],'track').values()])) )
    ref=nx.disjoint_union(ref,nets[n])

print(ref.node[ref.nodes()[0]].keys())

outDF = pandas.DataFrame(columns=ref.node[ref.nodes()[0]].keys())

fn = re.sub(".csv|.txt","",os.path.basename(srcpath))

for n in ref.nodes():
    outDF.loc[n,:] = [ref.node[n][x] for x in ref.node[ref.nodes()[0]].keys()]
#f=open( os.path.expanduser(srcpath+'/graphtable.csv'),'wb')
#outDF.to_csv(f, mode='w',sep='\t')
#f.close()

#Count the number of cells at each division number at each time point
regressDF=[]
for t in set(outDF.loc[:,'track']):
    curDF=outDF.loc[outDF.loc[:,'track']==t]
    #print curDF
    curDF.sort(['t','divs'],ascending=[True,True])
    for i in range(curDF.shape[0]-1):
        r1=curDF.loc[curDF.index[i],:]       
        r2=curDF.loc[curDF.index[i+1],:]
        t1=int(r1['t'])
        t2=int(r2['t'])
        for newT in range(t1,t2+1):
            if r1['divs'] < r2['divs'] and r1['divs']+2 > r2['divs'] and newT> (t1+t2)/2 :
                regressDF.append({'track':r1['track'], 'divs':r2['divs'],'t':newT})
                #print {'track':r1['track'], 'divs':r2['divs'],'t':newT}
            elif r1['divs'] == r2['divs']:
                regressDF.append({'track':r1['track'], 'divs':r1['divs'],'t':newT})

regressDF=pandas.DataFrame(regressDF,columns=('track', 'divs', 't'))            


DFch=outDF.columns[["ch" in x for x in outDF.columns]]
channelOrder=[int(re.sub("ch|Mean","", y)) for y in DFch]


outDir=os.path.join(os.path.dirname(os.path.realpath(srcpath)),'')
rgb = [colorsys.hsv_to_rgb(i ,1.,1.) for i in np.array(nx.get_node_attributes(ref,'divs').values())/ (float(max(nx.get_node_attributes(ref,'divs').values()))*1.3)]
rgbregress = [colorsys.hsv_to_rgb(i ,1.,1.) for i in np.array(regressDF.loc[:,'divs'])/ (float(max(regressDF.loc[:,'divs']))*1.3)]
nx.write_gml(ref, os.path.expanduser(outDir+fn+'.gml'))
nx.drawing.nx_pydot.write_dot(ref,os.path.expanduser(outDir+fn+'.dot'))

def graphviz_fixy(ref,yvals, prog='dot'):
    pos=nx.drawing.nx_pydot.graphviz_layout(ref, prog='dot')
    maxY = float(max([pos[x][1] for x in pos]))
    yvalsDict={pos.keys()[x]:yvals[x] for x in range(len(pos.keys()))}
    posy={x:maxY*(((yvalsDict[x] -min(yvals))/(max(yvals)-min(yvals)))) for x in pos.keys()}
    for x in pos:
        pos[x]=(pos[x][0],posy[x])
    return pos
    

pp=PdfPages(os.path.expanduser(outDir+fn+'analysis.pdf'))
for i in range(len(channelOrder)):
    plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,DFch[i]],c=rgb)#outDF.loc[:,'divs'])
    plt.xlabel('time')
    plt.ylabel(channelnames[channelOrder[i]-1]+'fluorescence')
    pp.savefig()
    plt.close()
    positions={x : np.array([ref.node[x]['t'], ref.node[x][DFch[i]]]) for x in ref.nodes()}
    nx.draw_networkx(ref,pos=positions, node_size =15,edge_color='black',linewidths= .1, width = .1,node_color = rgb, with_labels=False)
    plt.xlabel('time')
    plt.ylabel(channelnames[channelOrder[i]-1]+'fluorescence')    
    pp.savefig()
    plt.close()
    plt.scatter(outDF.loc[:,'divs'] ,outDF.loc[:,DFch[i]],c=rgb)#outDF.loc[:,'divs'])
    plt.xlabel('divisions')
    plt.ylabel(channelnames[channelOrder[i]-1]+'fluorescence')
    pp.savefig()
    plt.close()
    plt.close()
    plt.title(channelnames[channelOrder[i]-1]+'fluorescence')
    black = colour.Color("black")
    fluors = nx.get_node_attributes(ref,DFch[i]).values()
    maxFluors=int(max(fluors))
    minFluors=int(min(fluors))
    colors = list(black.range_to(colour.Color("red"), maxFluors-minFluors ))
    pos= graphviz_fixy(ref,yvals=[ref.node[x]['t'] for x in ref.nodes()])
    ax=plt.yticks(range(int(min([ref.node[x]['t'] for x in ref.nodes()])), int(max([ref.node[x]['t'] for x in ref.nodes()]))))    
    #pos=nx.drawing.nx_pydot.graphviz_layout(ref, prog='dot')
    nx.draw(ref, pos, with_labels=False,node_size =15,edge_color='black',linewidths= .1, width = .1,node_color=[colors[c-minFluors-1].get_rgb() for c in [int(x)for x in fluors]]  )  
    plt.ylabel("Frame")  
    #nx.drawing.draw_spring(ref, prog='dot',node_color = [colors[c].get_rgb() for c in [int(x)for x in nx.get_node_attributes(ref,DFch[i]).values()]])
    #nx.draw_networkx(ref,pos=positions, node_size =15,edge_color='black',linewidths= .1, width = .3,node_color = [colors[c].get_rgb() for c in [int(x)for x in nx.get_node_attributes(ref,DFch[i]).values()]],with_labels=False)
    pp.savefig()
    plt.close()    


plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'divs'],c=rgb)
plt.xlabel('time')
plt.ylabel('divisions')
pp.savefig()
plt.close()
plt.scatter(regressDF.loc[:,'t'] ,regressDF.loc[:,'divs'],c=rgbregress)
plt.xlabel('time')
plt.ylabel('divisions (Interploated)')
pp.savefig()
plt.close()
pp.close()
'''
'''


'''
#PRINTING A FIGURE
pp = PdfPages(os.path.expanduser('~/Desktop/'+sorcpath+'Networks.pdf'))        
fluors = nx.get_node_attributes(G,ch).values()
black = colour.Color("black")
colors = list(black.range_to(Color("red"), int(max(fluors)+1)))
nx.nx_pydot.write_dot(G,os.path.expanduser('~/Desktop/'+sorcpath+'Networks.dot'))        
#print nx.get_node_attributes(T,'track').values()
#print nx.get_node_attributes(T,'t').values()
print G.number_of_edges()
print G.number_of_nodes()
#positions= {nx.get_node_attributes(G,'t').keys()[x] : np.array([nx.get_node_attributes(G,'track').values()[x],nx.get_node_attributes(G,'t').values()[x]]) for x in range(len(nx.get_node_attributes(G,'t').values())) }        
print "building fig"
positions = nx.nx_pydot.graphviz_layout(G, prog='dot')
print "positions calculated"
#nx.get_node_attributes(T,'t')
pp.savefig( nx.draw_networkx(G,pos=positions, node_size =15,edge_color='black',linewidths= .1, width = .3,node_color = [colors[c].get_rgb() for c in [int(x)for x in nx.get_node_attributes(G,'fluorescence').values()]], with_labels=False))
pp.close()        
print 'saved'
sns.violinplot(regressDF.loc[regressDF.loc[:,'t']%50==1,'t'] , regressDF.loc[regressDF.loc[:,'t']%50==1,'divs'])
plt.xlabel('TimePoint')
plt.ylabel('divisions')
pp.savefig()
plt.close()

'''
