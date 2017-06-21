#Run from FIJI by pressing "[" to open the script editor. Open this file and run.

# File(label="Select a file") myFile
import math
from ij import IJ, ImagePlus
from ij.process import FloatProcessor
from ij.process import ImageProcessor as IP
import ij
import csv
import os
import ij.plugin.HyperStackConverter
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
import subprocess
import java.awt.Font as Font
import java.awt.Color as Color
import java.lang.Runtime as Runtime
import java.io.InputStreamReader
import java.io.BufferedReader

def setColor(color):
    IJ.run('Colors...', 'foreground=255')

def draw(i, j, array):
    path(i, j, array)
    setColor(255)
    IJ.run("Fill")
    setColor(255)
    IJ.run("Draw")

def drawDivNumber(imp,string, frame, x,y,col=0):
	if float(string) < 0:
		imp.setPosition(1,1,int(frame))
		setColor(col)
  		imp.setRoi(TextRoi( int(x)-diameter/2, int(y)-diameter/2,"X"))
  		IJ.run(imp, "Fill", "slice")
  	else:
		imp.setPosition(1,1,int(frame))
	  	IJ.setForegroundColor(col,col,col)
	  	imp.setRoi(OvalRoi(int(x)-diameter/2, int(y)-diameter/2,diameter,diameter))
	  	IJ.run(imp, "Fill", "slice")
	  	colFill = abs(col-255)
		IJ.setForegroundColor(colFill,colFill,colFill)
		imp.setRoi(TextRoi( int(x)-diameter/3, int(y)-diameter/1.5,re.sub("\\.0","",string),Font("Arial",Font.BOLD,int(diameter))))
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
gd.addChoice("Yes or No?",["Yes","No","GenerateResults"],"Yes")
gd.showDialog()
if (gd.wasCanceled()):
	sys.exit()
continuing = gd.getChoices().get(0).getSelectedItem()

if continuing=="Yes":
	continueFilePath=IJ.getFilePath("Choose ContinueFile")
	csvfile = open(continueFilePath, 'rU')
	spamreader = csv.reader(csvfile, delimiter=',')
	lll = spamreader.next()
	chosen=lll[0]
	diameter=float(lll[1])
	trackingFilename=spamreader.next()[0]
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
	trackFile= imPaths[0].split("/")[len(imPaths[0].split("/"))-3]
	inds=range(len(imPaths[0].split("/"))-2)
	trackPath=[imPaths[0].split("/")[i] for i in inds]
	trackPath= os.path.join(*trackPath)
	trackingFilename='/'+os.path.join(continueDir,continueFN+'_TrackingFile.csv')
	csvfile = open(str(continueDir)+ continueFN+'_ContinueFile.csv', 'w')
	spamwriter = csv.writer(csvfile, delimiter=',')
	spamwriter.writerow([chosen,diameter])
	spamwriter.writerow([trackingFilename])
	spamwriter.writerow(["ch","path"])
	for x in range(len(imPaths)):
		spamwriter.writerow([x+1,imPaths[x]])
	csvfile.close()


print("loading files, (This may take a while)")
print imPaths
for p in imPaths:
	if chosen=="TIFF":
		print p
		#imp=IJ.openImage(p)
		#imp.setOpenAsHyperStack(True)
		#imps.append(imp)
		o=Opener()
		o.setSilentMode(True)
		imp =o.openUsingBioFormats(p)
		imp.setOpenAsHyperStack(True)
		imps.append(imp)
	elif chosen=="DIR":
		print p
		imps.append(importDirAsStack(p))
				
print imps

if continuing =="GenerateResults":
	print "run results script"
	cFiles=[]
	gd = GenericDialog("Choose ContinueFile")
	gd.addMessage("Choose ContinueFile")
	gd.showDialog()
	continueFilePath=IJ.getFilePath("Choose ContinueFile")
	cFiles.append(continueFilePath)
	nextFile=True
	while(nextFile):
		gd = GenericDialog("Add another file to aggregate results?")
		gd.addMessage("Add another file to aggregate results?")
		gd.showDialog()
		if gd.wasCanceled():
				break
		else:
			newCFP=IJ.getFilePath("Choose ContinueFile")
			cFiles.append(newCFP)
	cFilesString= ''.join([" "+str(x) for x in cFiles])
	gd = GenericDialog("Choose Script")
	gd.addMessage("Find GraphReferenceColorReporter.py (or other)")
	gd.showDialog()
	if gd.wasCanceled():
				sys.exit()
	scriptFile=IJ.getFilePath("Choose Script")
	IJ.getString("Paste into Terminal: ", " ".join(["python " , scriptFile , cFilesString]) )
	sys.exit()
		


if os.path.exists(trackingFilename):
	csvfile = open(trackingFilename, 'rU')
	spamreader = csv.reader(csvfile, delimiter=',')
	keys=spamreader.next()
	try:
		for row in spamreader:
			tab.append( {keys[i]:row[i] for i in range(len(row)) } )
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
	#print x
	if  float(tab[i+1]['divs']) == (float(tab[i]['divs']) + 1):
		drawDivNumber(imp, str(x['divs']),float(x['t']),float(x['x']), float(x['y']))
	else:
		drawDivNumber(imp, str(x['divs']),float(x['t']),float(x['x']), float(x['y']))

IJ.run("Select None")
imp.setPosition(1,1,1) 
#Open the file to write to later
csvfile = open(trackingFilename, 'a')
spamwriter = csv.writer(csvfile, delimiter=',')

#IJ.setTool(Toolbar.HAND)
IJ.setTool(Toolbar.POINT)
canvas = WindowManager.getCurrentImage().getCanvas()
clicked = 0
clickedFrames = []

skipFrames = int(IJ.getNumber("Number Of Frames Between Observations",5))

lastFrame=1
nextJump=False
lastDiv=0
#HIT CANCEL ON Number prompt to save place
while imp.getNFrames()==frames:
	newClicked = canvas.getModifiers() == 16
	if not clicked and newClicked and imp.getFrame() not in clickedFrames:
		p = canvas.getCursorLoc()
		#PRESS CANCEL TO EXIT (RETURNS sys.minint)
		number = IJ.getNumber("Number Of Divisions",lastDiv)
		#Enter 999 if you made a mistake
		if number == 999:
			clicked = newClicked
			continue
		if number == sys.minint:
			#IJ.saveAs(imp, "TIFF", str(myPath) + 'Working.tif' )
			imp.close()
			csvfile.close()
			break
		#Get rid of last line of file
		if number == 989:
			#You have to close the file, read the whole thing, remove the last line, then rewrite the file to remove a line :(
			clicked = newClicked
			csvfile.close()
			csvfile = open(trackingFilename, 'rU')
			spamreader = csv.reader(csvfile, delimiter=',')
			lines = [row for row in spamreader]
			print lines
			lines = lines[:-1] 
			csvfile.close()
			csvfile = open(trackingFilename, 'w')
			spamwriter = csv.writer(csvfile, delimiter=',')
			for line in lines:
				spamwriter.writerow(line)
			csvfile.close()
			csvfile = open(trackingFilename, 'a')
			spamwriter = csv.writer(csvfile, delimiter=',')
			continue
		else:
			clickedFrames.append(imp.getFrame())
			fluor = []
			for impF in imps:
				impF.setPosition(1,1,imp.getFrame())
				impF.setRoi(OvalRoi(p.x-diameter/2,p.y-diameter/2, diameter, diameter))
				fluor.append(impF.getStatistics().mean)
			spamwriter.writerow([imp.getFrame(),track,p.x,p.y,number] + fluor )
			tlf=lastFrame
			lastFrame=imp.getFrame()
			lastDiv=number
			drawDivNumber(imp, str(number),float(imp.getFrame()),float(p.x), float(p.y))
			print [imp.getFrame(),track,p.x,p.y,number] + fluor
			IJ.run("Select None")

			#Enter -1 if a cell dies, -2 if you lost a cell, 989 to undo last action
			if imp.getFrame() == imp.getNFrames() or number < 0:
				clickedFrames=[]
				clicked = False
				lastFrame=1
				lastDiv=0
				nextJump=False
				track = track +1
				imp.setPosition(1,1,1)
				sleep(.5)
			else:
				if (tlf != imp.getFrame()-skipFrames) == nextJump:
					nextJump=False
					for i in range(skipFrames):
						sleep(.14)
						imp.setPosition(1,1,int(imp.getFrame()+1))
						#amount of time between each frame
						
				else:
					imp.setPosition(1,1,int(imp.getFrame()+1))
					nextJump=True	
					
	clicked = newClicked
	sleep(.06)

for impF in imps:
	impF.close()
csvfile.close()