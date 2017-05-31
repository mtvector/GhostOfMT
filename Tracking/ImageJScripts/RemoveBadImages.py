#Start by clicking on folder
#!!!Will delete original image files

import math
from ij import IJ, ImagePlus
from ij.process import FloatProcessor
from ij.process import ImageProcessor as IP
import ij
import csv
import os
from time import sleep
import sys
from ij.gui import Toolbar
from ij import WindowManager
import ij.gui.OvalRoi as OvalRoi
import ij.gui.TextRoi as TextRoi
import re
myFile = IJ.getDirectory("Choose a directory")
print str(myFile)
def sort_human(l):
  convert = lambda text: float(text) if text.isdigit() else text
  alphanum = lambda key: [ convert(c) for c in re.split('([-+]?[0-9]*\.?[0-9]*)', key) ]
  l.sort( key=alphanum )
  return l
  
dirname = myFile
print dirname

ch2list=  sort_human(os.listdir(dirname+'Ch2'))
phlist= sort_human(os.listdir(dirname+'Ph'))
ch3list= sort_human(os.listdir(dirname+'Ch3'))


def trimmedMean(lst,num):
    trimmed_lst = sorted(lst)[num:-num]
    return sum(trimmed_lst) / len(trimmed_lst)

def removeOutliers(stack):
    hists = []
    for i in range(1,stack.getNSlices()):
    	stack.setPosition(i)
        hists.append(stack.getProcessor().getHistogram())
    meanHist = [0]*len(hists[1])
    totalHist = [0]*len(hists[1])
    lenHist = [0]*len(hists[1])
    for h in hists:
        for i in range(len(h)):
            totalHist[i] = totalHist[i] + h[i]
            lenHist[i] = lenHist[i] + 1
    for i in range(len(totalHist)):
        meanHist[i] = totalHist[i]/lenHist[i]
    histDiffs = [sum([abs(a) for a in [a-b for a,b in zip(meanHist,h)]]) for h in hists]
    #trimMean = trimmedMean(histDiffs, 12)
    #inds = [q < (trimMean + 4*max(histDiffs))/6 for q in histDiffs]
    #inds_selection= [elem for elem in range(len(inds)) if inds[elem]]
    return(histDiffs)

imp = IJ.run("Image Sequence...", "open=["+ str(myFile)+"Ph] sort")
imp = IJ.getImage()
imp.show()

histDiffs= removeOutliers(imp)
sortedinds=sorted(range(len(histDiffs)),key=lambda x:histDiffs[x],reverse=True)
#print sorted(histDiffs)

i=0
while i < len(sortedinds):
	x = sortedinds[i]
	imp.setPosition(x+1)
	string = IJ.getString("Delete?",'n')
	if string == 'y':
		print 'deleted'
		print x
		delly2 = IJ.openImage(dirname + 'Ch2/'+ ch2list[x] )
		#dellyph = IJ.openImage(dirname + 'Ph/'+ phlist[x])
		delly3 = IJ.openImage(dirname + 'Ch3/' + ch3list[x] )
		delly2.show()
		delly3.show()
		#dellyph.show()
		surestring = IJ.getString("Are you sure?",'n')
		if surestring == 'y':
			print "OUCH"
			os.remove(dirname + 'Ch2/'+ch2list[x])
			os.remove(dirname + 'Ch3/'+ch3list[x])
			os.remove(dirname + 'Ph/'+phlist[x])
		delly2.close()
		delly3.close()
		#dellyph.close()
	elif string == 'b':
		i = i - 2
	elif string == "quit":
		imp.close()
		break
	i = i + 1

	

 
