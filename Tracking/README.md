# ImageJ Scripts
**Used to track cells over time, counting divisions, and monitoring reporter fluorescence.**

See The ImageJScripts README.md file for instructions

Scripts in the ImageJScripts should be run through the script editor of the FIJI software package. (Opened by pressing "[")

The main user-friendly script is ManTrack.py in the ImageJScripts directory. This can launch the manual tracker, as well as a script to generate output analysis images. See README.md in the ImageJScripts Folder for further instructions.

Preprocessing to remove background of TIFF or directories of image files can be done with RemoveBackgrounds.py

RemoveBadImages **(BE CAREFUL WITH THIS ONE)** can sort and delete multiple channels' failed images (condensation) simultaneously.

Images from BioCT output can be sorted into directories for proper channels with SortImages.py

trackmater.py and run TrackmateAnalysis were used to search for parameters with genetic algorithm