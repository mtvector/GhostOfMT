#Run from FIJI by pressing entering the script editor with "[" and opening this file

# File(label="Select a file") myFile
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
import shutil
import os
import re
from ij import IJ
from ij.io import Opener
from ij.plugin import Concatenator
from jarray import array
import ij.gui.GenericDialog as GenericDialog 

def setColor(color):
    IJ.run('Colors...', 'foreground=255')

def draw(i, j, array):
    path(i, j, array)
    setColor(255)
    IJ.run("Fill")
    setColor(255)
    IJ.run("Draw")

def drawDivNumber(imp,string, frame, x,y,col=255):
	if float(string) > 0:
  		imp.setRoi(TextRoi( int(x)-diameter/2, int(y)-diameter/2,"X"))
  		IJ.run(imp, "Fill", "slice")
  	else:
		imp.setPosition(1,1,int(frame))
		setColor(col)
		imp.setRoi(TextRoi( int(x), int(y),string))
	  	IJ.run(imp, "Fill", "slice")
	  	imp.setRoi(OvalRoi(int(x)-diameter/2, int(y)-diameter/2,diameter,diameter))
	  	IJ.run(imp, "Fill", "slice")
  	

def importDirAsStack(specdir):
	timeseries = []
	for f in sorted(os.listdir(specdir)):
		imp = Opener.openUsingBioFormats(os.path.join(specdir, f))
		if imp != None:
			imp.setOpenAsHyperStack(False)
			timeseries.append(imp)
	jaimp = array(timeseries, ImagePlus)
	ccc = Concatenator()
	allimp = ccc.concatenate(jaimp, False)
	allimp.setOpenAsHyperStack(True)
	return allimp



tab = []
track = 0
imps=[]
imPaths=[]

gd = GenericDialog("Continue a track?")
gd.addChoice("Yes or No?",["Yes","No","GenerateResults"],"No")
gd.showDialog()
if (gd.wasCanceled()):
	sys.exit()
continuing = gd.getChoices().get(0).getSelectedItem()

if continuing=="Yes":
	continueFilePath=IJ.getFilePath("Choose ContinueFile")
	csvfile = open(continueFilePath, 'r')
	spamreader = csv.reader(csvfile, delimiter=',')
	lll = spamreader.next()
	chosen=lll[0]
	diameter=float(lll[1])
	spamreader.next()
	try:
		for row in spamreader:
			imPaths.append(row[1] )
		csvfile.close()
	except:
		pass

elif continuing=="No":
	gd = GenericDialog("What type of image for channel")
	gd.addMessage("Pick tracking channel first!")
	gd.addChoice("TIFF or Directory of Images?",["TIFF","DIR"],"TIFF")
	gd.showDialog()
	
	if gd.wasCanceled():
		sys.exit()
	chosen = gd.getChoices().get(0).getSelectedItem()

	diameter = IJ.getNumber("# Pixel Diameter Of Nuclei (Suggested=10)",10)
	
	nextImage = True
	while(nextImage):
		if len(imPaths)>0:
			gd = GenericDialog("Add channel "+ str(len(imPaths)))
			gd.addMessage("Add Another Channel?")
			gd.showDialog()
			if gd.wasCanceled():
				break
		if(chosen=="TIFF"):
			myFile=IJ.getFilePath("Choose another Image?")
			if(myFile!=None):
				imPaths.append(myFile)
			else:
				nextImage=False
		else:
			myFile = IJ.getDirectory("Choose another directory?")
			if(myFile!=None):
				imPaths.append(myFile)
			else:
				nextImage=False
			
	
	gd = GenericDialog("Where put file?")
	gd.addMessage("Pick a directory for ContinueFile!")
	gd.showDialog()
	continueDir = IJ.getDirectory("Choose a directory")
	continueFN = IJ.getString("ContinueFile Name?",str(imPaths[0].split("/")[len(imPaths[0].split("/"))-3]) )
	csvfile = open(str(continueDir)+ continueFN+'_ContinueFile.csv', 'w')
	spamwriter = csv.writer(csvfile, delimiter=',')
	spamwriter.writerow([chosen,diameter])
	spamwriter.writerow(["ch","path"])
	for x in range(len(imPaths)):
		spamwriter.writerow([x+1,imPaths[x]])
	csvfile.close()


print("loading files, (This may take a while)")
print imPaths
for p in imPaths:
	if chosen=="TIFF":
		print p
		imps.append(IJ.open(p))
	elif chosen=="DIR":
		print p
		imps.append(importDirAsStack(p))
				

#Process a path for 
myFile= imPaths[0].split("/")[len(imPaths[0].split("/"))-3]
inds=range(len(imPaths[0].split("/"))-2)
myPath=[imPaths[0].split("/")[i] for i in inds]
myPath= os.path.join(*myPath)
trackingFilename='/'+os.path.join(myPath,myFile+'trackingFile.csv')

if continuing =="GenerateResults":
	print "run results script"


	
print trackingFilename
if os.path.exists(trackingFilename):
	csvfile = open(trackingFilename, 'r')
	spamreader = csv.reader(csvfile, delimiter=',')
	keys=spamreader.next()
	try:
		for row in spamreader:
			tab.append( {keys[i]:row[i] for i in range(len(keys)) } )
		csvfile.close()
		track = int(tab[len(tab)-1]['track']) +1
		print "reading"
	except:
		pass
else:
	csvfilenew = open(trackingFilename, 'w')
	spamwriter = csv.writer(csvfilenew, delimiter=',')
	spamwriter.writerow(["t","track","x","y","divs"] + ["ch"+str(x)+"Mean" for x in range(1,len(imPaths)+1)])
	csvfilenew.close()
print(trackingFilename)

print "Starting at track=" + str(track)

imp = imps[0]
imp.show()
frames = imp.getNFrames()


for i in range(len(tab)-1):
	x=tab[i]
	print x
	if  float(tab[i+1]['divs']) == (float(tab[i]['divs']) + 1):
		drawDivNumber(imp, str(x['divs']),float(x['t']),float(x['x']), float(x['y']),col=250)
	else:
		drawDivNumber(imp, str(x['divs']),float(x['t']),float(x['x']), float(x['y']))


#Open the file to write to later
csvfile = open(trackingFilename, 'a')
spamwriter = csv.writer(csvfile, delimiter=',')


#IJ.setTool(Toolbar.HAND)
IJ.setTool(Toolbar.POINT)
canvas = WindowManager.getCurrentImage().getCanvas()
clicked = 0
clickedFrames = []

#HIT CANCEL ON Number prompt to save place
while imp.getNFrames()==frames:
	newClicked = canvas.getModifiers() == 16 
	if clicked and not newClicked and imp.getFrame() not in clickedFrames:
		p = canvas.getCursorLoc()
		#PRESS CANCEL TO EXIT (RETURNS sys.minint)
		number = IJ.getNumber("Number Of Divisions",0)
		#Enter 999 if you made a mistake
		if number == 999:
			clicked = newClicked
			continue
		if number == sys.minint:
			#IJ.saveAs(imp, "TIFF", str(myPath) + 'Working.tif' )
			imp.close()
			csvfile.close()
			break
		else:
			clickedFrames.append(imp.getFrame())
			imp.setRoi(OvalRoi(p.x-diameter/2,p.y-diameter/2, diameter, diameter))
			fluor = []
			for impF in imps:
				impF.setPosition(1,1,imp.getFrame())
				impF.setRoi(OvalRoi(p.x-diameter/2,p.y-diameter/2, diameter, diameter))
				fluor.append(impF.getStatistics().mean)
			spamwriter.writerow([imp.getFrame(),track,p.x,p.y,number] + fluor )
			setColor(255)
			IJ.run(imp, "Fill", "slice")
			IJ.run("Select None")
			#Enter -1 if a cell dies, -2 if you lost a cell
			if imp.getFrame() == imp.getNFrames() or number < 0:
				clickedFrames=[]
				clicked = False
				track = track +1
				imp.setPosition(1,1,1)
				sleep(.75)

	clicked = newClicked
	sleep(.06)

for impF in imps:
	impF.close()



