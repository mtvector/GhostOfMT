from ij import IJ, ImagePlus
import os
import ij.gui.GenericDialog as GenericDialog 

#This is how to make a dialog to control the if, else based on structure of input
gd = GenericDialog("What type of image")
gd.addChoice("TIFF or Directory of Images?",["TIFF","DIR"],"TIFF")
gd.showDialog()

if (gd.wasCanceled()):
	sys.exit()
chosen = gd.getChoices().get(0).getSelectedItem()

if(chosen=="TIFF"):
	imPath = IJ.getFilePath("Choose and Image")
	print(imPath)
	if(imPath!= None):
		imp=IJ.open(imPath)
		IJ.run(imp, "Subtract Background...", "rolling=100 disable")
		IJ.saveAs(imp, "TIFF", imPath.split(".")[0]+"SubBack")
else:
	myFile = IJ.getDirectory("Choose a directory")
	if not os.path.exists(str(os.path.normpath(myFile))+"SubBack"):
		os.mkdir(str(os.path.normpath(myFile))+"SubBack")
	if os.path.exists(str(myFile)):
		for f in os.listdir(myFile):
			if f !=".DS_Store":
				print f
				imp = IJ.openImage(os.path.join(str(myFile),f))
				IJ.run(imp, "Subtract Background...", "rolling=100 disable")
				IJ.saveAs(imp, "PNG", os.path.join(str(str(os.path.normpath(myFile))+"SubBack"),f))

