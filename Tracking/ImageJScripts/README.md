#The ManTracker
This is a Jython ImageJ Macro. Jython uses Python 2.7 syntax, however most python libraries can't be loaded (although all Java and [Java] ImageJ libraries can, along with all the standard ImageJ Macro language functions). There shouldn't be any setup required for just tracking, as everything it needs should be included in the FIJI install. To generate the figures, you will need to install a few python packages, which will be explained later.

To use ManTrack.py, first open the FIJI script editor by pressing "[" when the main bar is open. Open the file ManTrack.py and click the run button to start the script.

The first dialog box that will pop up will ask if you want to continue a track, start a new track, or run a script to generate results.

####Starting a new tracking dataset for a FOV (field of view)

If you are starting a new tracking dataset, the first thing it will ask for is the structure of your image stacks. You can either have a single multi-image TIF stack, or a directory of individual images for each channel. You'll then be asked to pick a diameter of nuclei, which can change depending on your magnification, cell type etc. You won't be able to change this later so choose wisely! Immediately after that you pick the file/directory that you will be tracking nuclei/cells on. If you want to also record fluorescence for this channel, when prompted to pick another channel, you should select this channel again. When you have selected all the channels you wish to record mean fluorescence values for, hit cancel. You will then pick a directory in which to put your output files. I recommend making a new directory, named in a way that you know exactly which stack and FOV you are tracking on. It will then make a TrackingFile which will record the data about your tracking, and a ContinueFile, which stores metadata and the paths to the images. **If you move your image files, you will have to change the paths in the continuefile or the script will be unable to find them.**

####Continuing an old tracking dataset for a FOV

If you wish to continue tracking, say yes on the first menu. Then all you have to do is find the continuefile and the script should load up everything you've done so far. 

#Tracking

####Tracking a cell

Pick the number of frames you want to have between observations. The program will try to help you to be observing at the same timepoints within a lineage.

The first step is to click in the center of the nucleus of a cell. For the first one, I recommend doing this on the very first frame of the image. This will be the root of the lineage tree. The program should automatically step forward one image once and then the number of images you have selected. If you obeserve a cell dividing, go to the last frame where the cell is one and click it. It will then step forward one frame, which should be the first frame where you see two nuclei. At each division, pick only one daughter cell to follow! You must track only one cell from your starting point to the end of the time series (or the point at which the cell is lost). The reason why is discussed in limitations.

You will track the other daughter cell by going back to the frame before the cell divided, and clicking at the center of the dot you clicked before. Enter the same division number you entered before (# shown in the dot).

####Secret Codes

**999** If you accidentally click on the image DON'T PRESS CANCEL. It will exit out everything. Instead enter **999** and hit okay and it will return

**989** If you think your previous data entry was incorrect, enter **989** and it will delete the last entry to the TrackingFile.

**-1** Signifies that the cell you are tracking was observed to have died at this point. It will end the track and return to the beginning of the stack.

**-2** Signifies that you have lost track of the cell you were tracking. It will end the track and return to the beginning of the stack.

**Cancel** Closes everything and finalizes output file.

####Things to know

Everything should save as you go along. To exit, just click on the image and hit cancel. Or exit the image window and don't save the open image (it can be reconstructed each time with the data from the continuefile and the trackingfile). It might take a little bit for the tracking data to appear in the file. Often it helps to hit run on the ManTrack script, then cancel out (something about ImageJ not finalizing the trackingfile until all ties to the previous run are broken).

####Limitations

This program is really dumb. It creates the lineages by detecting that two clicks on the same frame are within the diameter you entered, and merging them in processing. Also it's very hard to synchronize observations when the cell divisions are not synchronize unless you are clicking in observations on every frame. 

It is best to do something to normalized the fluorescence in the channels. It seems the BioCT can do this for you when you output your files, however the RemoveBackgrounds.py ImageJ script can help you do this as well. There is also an imageJ function called Subtract Background that can help with this. If the fluorescence looks skewed across the image, your data will include this bias unless you do something to correct it!

#Generate Results

Before you run the results you'll want to make sure you have a number of python packages installed. To do this you will use pip, the python package manager. Open your terminal app (on macs it is inside the utilities folder in the applications folder) and enter:

sudo easy_install pip 

This should install pip (you'll need to enter your password). Next you use pip to install a few crucial packages for the plotting by entering the following on your command line:

pip install pandas

pip install matplotlib

pip install networkx

pip install colour

pip install colorsys

pip install numpy

(You should have python 2.7 automatically installed, but you can make sure you have python installed by typing "whereis python" in the command line, and it should give you a path) If you don't have python installed, install python first.

Once you have done this, you can run ManTrack.py through the script editor and select GenerateResults. It will ask you first to select your continue file. If you want to include multiple FOVs from the same experiment, you can then select as many continuefiles as you want, if you have tracked cells in many FOVs and you want to aggregate them all into a single figure. Then you will personally locate the file GraphReferenceColorReporter.py and then it will give you a line to copy and paste into the terminal.

**I apologize for the lack of user friendliness of running the results script.** It is impossible to run real python scripts from a FIJI jython script, so you get this duct-taped code mess. My B.

The output files should go to the same directory as your first continuefile was in. These include a pdf of the graphs of the graphs https://en.wikipedia.org/wiki/Graph_(discrete_mathematics) as well as .dot and .gml files so you can take your analysis to an external platform, like Cytoscape. 

#Automated tracking

I experimented extensively with the trackmate plugin in the FIJI package. According to my analysis and the word of the plugin's author, I have concluded that while it is very good at detecting round objects and tracking them between frames, it is not accurate at detecting divisions, so if you are hoping to draw cell lineages. This plugin is not viable. If you want a machine to do the tracking for you, I recommend trying Ilastik (http://ilastik.org/) and running it on a VNC server. Make sure you remove condensation ruined images before tracking though! They will totally destroy results from any program that can't detect them. 