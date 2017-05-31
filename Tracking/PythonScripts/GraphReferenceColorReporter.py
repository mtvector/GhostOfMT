# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:36:22 2016

@author: mschmitz
"""

# -*- coding: utf-8 -*-
#matplotlib.use('WxAgg')
import os
import numpy as np
#import PIL
import pandas
import networkx as nx
#import shutil
import matplotlib.pyplot as plt
#import wx
#import math
from matplotlib.backends.backend_pdf import PdfPages
import colorsys
import seaborn as sns


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
    
def mergeNodes(G, n1, n2):
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

def loadNets(srcpath):
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

spaths=[os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x-04572y001856-10x-FLsorted"),
os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x004790y003248-10x-FLsorted"),
os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x-03952y000264-10x-FLsorted")]
   
srcpath = os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x-04572y001856-10x-FLsorted")
#srcpath = os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x004790y003248-10x-FLsorted")
#srcpath = os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x-03952y000264-10x-FLsorted")


#ref = loadNets(srcpath)
nets=[]
for s in spaths:
    nets.append(loadNets(s))
ref=nx.DiGraph()

for n in range(len(nets)):
    if len(nx.get_node_attributes(ref,'track').values())>0:
        mx=max(nx.get_node_attributes(ref,'track').values())
        nx.set_node_attributes(nets[n],'track', dict( zip(nx.get_node_attributes(nets[n],'track').keys(), [x+mx+1 for x in nx.get_node_attributes(nets[n],'track').values()])) )
    ref=nx.disjoint_union(ref,nets[n])

print(ref.node[ref.nodes()[0]].keys())

outDF = pandas.DataFrame(columns=ref.node[ref.nodes()[0]].keys())

fn = os.path.basename(srcpath)

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
print regressDF

rgb = [colorsys.hsv_to_rgb(i ,1.,1.) for i in np.array(nx.get_node_attributes(ref,'divs').values())/ (float(max(nx.get_node_attributes(ref,'divs').values()))*1.3)]
rgbregress = [colorsys.hsv_to_rgb(i ,1.,1.) for i in np.array(regressDF.loc[:,'divs'])/ (float(max(regressDF.loc[:,'divs']))*1.3)]

pp=PdfPages(os.path.expanduser('~/Desktop/'+fn+'analysis.pdf'))
plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'adjFluorescence'],c=rgb)#outDF.loc[:,'divs'])
plt.xlabel('time')
plt.ylabel('fluorescence')
#plt.show()
#plt.draw()
pp.savefig()
plt.close()
positions={x : np.array([ref.node[x]['t'], ref.node[x]['adjFluorescence']]) for x in ref.nodes()}
nx.draw_networkx(ref,pos=positions, node_size =15,edge_color='black',linewidths= .1, width = .3,node_color = rgb, with_labels=False)
pp.savefig()
plt.close()
plt.scatter(outDF.loc[:,'divs'] ,outDF.loc[:,'adjFluorescence'],c=rgb)#outDF.loc[:,'divs'])
plt.xlabel('divisions')
plt.ylabel('fluorescence')
#plt.show()
#plt.draw()
pp.savefig()
plt.close()
plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'divs'],c=rgb)
plt.xlabel('time')
plt.ylabel('divisions')
pp.savefig()
plt.close()
sns.violinplot(regressDF.loc[regressDF.loc[:,'t']%50==1,'t'] , regressDF.loc[regressDF.loc[:,'t']%50==1,'divs'])
plt.xlabel('TimePoint')
plt.ylabel('divisions')
pp.savefig()
plt.close()
plt.scatter(regressDF.loc[:,'t'] ,regressDF.loc[:,'divs'],c=rgbregress)
plt.xlabel('time')
plt.ylabel('divisions (Interploated)')
#plt.show()
#plt.draw()
pp.savefig()
plt.close()
pp.close()


        
'''
#PRINTING A FIGURE
pp = PdfPages(os.path.expanduser('~/Desktop/'+sorcpath+'Networks.pdf'))        
fluors = nx.get_node_attributes(G,'fluorescence').values()
black = Color("black")
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
'''
