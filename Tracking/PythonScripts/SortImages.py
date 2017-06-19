# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 11:40:26 2016

@author: mschmitz
Need to install pandas with sudo pip install pandas before using
Command line usage:
python PATH_TO_sortImages.py PATH_TO_DIRECTORY 
or to do multiple
python PATH_TO_sortImages.py PATH_TO_OUTER_DIRECTORY/*
"""
import pandas
import os
import sys
from shutil import copyfile

overpath = sys.argv[1]
#overpath = '~/Desktop/ss HEP1 Jan 27-B05_2016012700205_x000000y000000-04x-FL'
overpath = os.path.expanduser(overpath)
micro = pandas.DataFrame.from_csv(overpath+'/micro.csv',sep='\t', encoding='utf-16')
files= micro.loc[:,'File Name']
for f in files:
    f = f.replace('\\','/')
    print(f)
    if not os.path.exists(overpath+'/sorted'):
        os.makedirs(overpath+'/sorted')
    plist=f.split('/')
    if not os.path.exists(overpath+'/sorted/'+ plist[1]):
        os.makedirs(overpath+'/sorted/'+plist[1])
    copyfile(overpath+'/'+f, overpath+'/sorted/'+ plist[1]+ '/' + plist[len(plist)-1])