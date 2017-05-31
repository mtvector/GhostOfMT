# -*- coding: utf-8 -*-
import matplotlib
#matplotlib.use('WxAgg')
import os
import numpy as np
#import PIL
import pandas
import networkx as nx
import sys
#import shutil
import matplotlib.pyplot as plt
from colour import Color
#import wx
#import math
from matplotlib.backends.backend_pdf import PdfPages
import re
from GraphToImage import *


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



srcpath = os.path.expanduser(sys.argv[1])
if len(sys.argv) < 2:
    interact=False    
elif sys.argv[2] == '0' :
    interact = False
else:
    interact = True

#python GraphColorReporter.py "~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x-04572y001856-10x-FLsorted_5_11_3_17_21_3_21_True" 0
#srcpath=os.path.expanduser("~/code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x-06080y001764-10x-FLsortedtracked")
#overpath = 'code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/'
#for sorcpath in os.listdir(os.path.join(hmpath,overpath)):
#for sorcpath in os.listdir(overpath):
#for sorcpath in ['CB H9pax6td NBnog may10-2_2016051000002_x-06080y001764-10x-FLsortedtracked']:
edgeFrames={}
trackFrames={}
fn = os.path.basename(srcpath)
srcDir = os.path.dirname(srcpath)
print sys.argv
print fn
print srcpath
if 1==1:        
#if  os.path.isfile(os.path.expanduser(srcpath+'/graphtable.csv')):
    for root, directories, filenames in os.walk(srcpath):
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
    

    #print G.nodes()  
    #f=open('~/Desktop'+sorcpath+'.txt','wb')
    #writer=csv.writer(f, delimiter=',')
    #plt.hist(nx.get_node_attributes(G,'t').values(),25)
    #plt.show()
    roots = []
    for n in G.nodes():
        #node = G.nodes()[n]
        if len(G.predecessors(n)) == 0 and nx.get_node_attributes(G,'t')[n] < 900:
            roots.append(n)
    print("screening rootless nodes")
    for n in G.nodes():
        div = getNumToRoot(G,n,roots)
        if div == None:
            G.remove_node(n)
        else:
            G.node[n]['divs'] = div
            
    #print("getting num to root")         
    #for n in G.nodes():
    #    div = getNumToRootRecursive(G,n)
    #    G.node[n]['divs'] = div
    
    #f=open('~/Desktop'+sorcpath+'.txt','wb')
    #writer=csv.writer(f, delimiter=',')
    #print G.nodes()
    #if len(G.nodes())<1:
    #    continue
    print "past continue"
    outDF = pandas.DataFrame(columns=G.node[G.nodes()[0]].keys())

    for n in G.nodes():
        outDF.loc[n] = G.node[n].values()
    f=open( os.path.expanduser(srcpath+'/graphtable.csv'),'wb')
    outDF.to_csv(f, mode='w',sep='\t')
    f.close()
    print outDF

print("doneif")
print os.path.expanduser(srcpath+'/graphtable.csv')
outDF = pandas.read_csv(os.path.expanduser(srcpath+'/graphtable.csv'),sep='\t',index_col=0)
print outDF
#plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'fluorescence'],c=outDF.loc[:,'divs'])
#plt.xlabel('time')
#plt.ylabel('fluorescence')
#plt.show()
#if interact==True:        
if False:
    annotes = list(outDF.index)
    fig, ax = plt.subplots()
    ax.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'fluorescence'],c=outDF.loc[:,'divs'])
    ax.set_xlabel('time')
    ax.set_ylabel('fluorescence')
    af =  ImageFinder(outDF.loc[:,'t'] ,outDF.loc[:,'fluorescence'], annotes, dataDF=outDF,path=re.sub('tracked','', srcpath), ax=ax)
    print("theconnection")        
    fig.canvas.mpl_connect('button_press_event', af)
    print("theshow")
    plt.show()
print "open pdf"
pp=PdfPages(os.path.expanduser('~/Desktop/'+fn+'analysis.pdf'))
print "pdfopen"
plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'fluorescence'],c=outDF.loc[:,'divs'])
plt.xlabel('time')
plt.ylabel('fluorescence')
#plt.show()
#plt.draw()
pp.savefig()
plt.close()
plt.scatter(outDF.loc[:,'divs'] ,outDF.loc[:,'fluorescence'],c=outDF.loc[:,'divs'])
plt.xlabel('divisions')
plt.ylabel('fluorescence')
#plt.show()
#plt.draw()
pp.savefig()
plt.close()
plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'divs'])
plt.xlabel('time')
plt.ylabel('divisions')
#plt.show()
#plt.draw()
pp.savefig()
plt.close()
plt.scatter(outDF.loc[:,'t'] ,outDF.loc[:,'estdiam'])
plt.xlabel('time')
plt.ylabel('diameter')
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
