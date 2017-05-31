#To be used with RunTrackmateAnalysis
#Trackmate found to be ineffective by genetic algorithm search of parameters
print 'start\nstart\nstart'
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.features.spot import SpotRadiusEstimatorFactory
from ij import IJ
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import sys
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer as TrackDurationAnalyzer
import fiji.plugin.trackmate.features.track as tk 
import os
import re
from ij import ImagePlus, IJ
from net.imagej import ImgPlus as ImgPlus
from ij.io import Opener
from ij.plugin import Concatenator
from jarray import array
import csv
import fiji.plugin.trackmate.visualization as visualization
import fiji.plugin.trackmate.visualization.trackscheme.TrackScheme as TrackScheme
import fiji.plugin.trackmate.visualization.SpotColorGeneratorPerTrackFeature as SpotColorGeneratorPerTrackFeature
import fiji.plugin.trackmate.visualization.PerTrackFeatureColorGenerator as PerTrackFeatureColorGenerator
import fiji.plugin.trackmate.features.ModelFeatureUpdater as ModelFeatureUpdater
import fiji.plugin.trackmate.features.track.TrackIndexAnalyzer as TrackIndexAnalyzer
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzer as SpotIntensityAnalyzer
import shutil
from ij.process import FloatProcessor
import net.imglib2.img.ImagePlusAdapter as ImagePlusAdapter
import copy
import math
from ij.plugin import RGBStackMerge
from ij.process import ImageConverter
import subprocess

def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

def trimmedMean(lst,num):
    trimmed_lst = sorted(lst)[num:-num]
    return sum(trimmed_lst) / len(trimmed_lst)

def removeOutliers(stack):
    hists = []
    for i in stack:
        hists.append(i.getProcessor().getHistogram())
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
    trimMean = trimmedMean(histDiffs, 12)
    print max(histDiffs)
    print trimMean
    #inds = [q < trimMean*2.8 for q in histDiffs]
    #Hacked... Need to replace with real clustering method
    inds = [q < (trimMean + 4*max(histDiffs))/6 for q in histDiffs]
    inds_selection= [elem for elem in range(len(inds)) if inds[elem]]
    #print(inds_selection)
    return(inds_selection)

def meanInRadius(x,y,imp,radius,xlen=1000,ylen=1000):
    x = x
    y = y
    grabCoord = []
    pixVals = []
    for i in range(max([1,x-radius]) , min([xlen, x+radius])+1 ):
        for j in range(max([1,y-radius]) , min([ylen, y+radius])+1 ):
            if math.sqrt((i-x)**2 + (j-y)**2) <=radius:
                grabCoord.append([i,j])
        #print i,j
    for i,j in grabCoord:
        pixVals.append(imp.getPixel(i,j)[0])
        #print dir(imp.getProcessor())
        #imp.getProcessor().putPixel(i,j,0)
    if min(pixVals) < 0:
        print min(pixVals)
    #print pixels[int(x+(y%ylen -1)*xlen)]
    #print pixVals
    return sum(pixVals)/len(pixVals)
    

#os.chdir()
#os.listdir()
#os.getcwd()

#srcpath = IJ.getFilePath('Choose the first file')
#mpath = os.path.realpath('~')
#hmpath = hmpath.replace('~','')
#srcpath=os.path.join(hmpath,'code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/CB H9pax6td NBnog may10-2_2016051000002_x004790y003248-10x-FLsorted/')
#overpath = 'code/data/cb/imageAnalysis/CB H9pax6td NBnog may10-2_2016051000002/'


outpath ='_'.join(map(str,sys.argv[1:]))
print sys.argv
print "outpath= "+ outpath

radius = float(sys.argv[2])
threshold = float(sys.argv[3])
maxGap=int(sys.argv[4])
maxLinkDist=float(sys.argv[5])
maxGapDist = float(sys.argv[6])
trackDuration=int(sys.argv[7])
qualityMin=5/float(sys.argv[8])
subBackground = sys.argv[9]
if subBackground:
    threshold = threshold/30


if not os.path.exists(outpath):
    os.mkdir(outpath)
    if len(os.listdir(outpath))<2:
        #for srcpath in os.listdir(os.path.join(hmpath,overpath)):
        for srcpath in [sys.argv[1]]:
            if 'FLsorted' in srcpath and 'tracked' not in srcpath:
                #srcpath = os.path.join(hmpath,overpath,srcpath)
                filename = os.path.basename(srcpath)
                srcDir = os.path.dirname(srcpath)
                
                GRlist = {}
                
                for root, directories, filenames in os.walk(srcpath):
                    #print root, directories, filenames, "ENDENDENDENDEND\n\n\n\n\n\n\n\n\n"
                    #print os.path.basename(root)
                    #print len(filenames)
                    if ( os.path.basename(root) not in GRlist and len(filenames)>100 and os.path.basename(root) != 'tracks' ):
                        GRlist[os.path.basename(root)] = []
                    for filename in filenames:
                        if '.png' in filename:
                            GRlist[os.path.basename(root)].append(root+'/'+filename)
                 
                print srcpath
                print 'channels: ', len(GRlist)
                print dict.keys(GRlist)
                
                from collections import defaultdict
                #GRlist = sorted(GRlist)
                #timeseries = dict.fromkeys(GRlist.keys(),[])
                timeseries = defaultdict(list)
                for x in sorted(GRlist.keys()):
                    if not os.path.isfile(srcpath+'/'+x+".tif"):
                        print x
                        for timepoint in sorted(GRlist[x]):
                            thisfile = timepoint
                            #print timepoint
                            imp = Opener.openUsingBioFormats(os.path.join(srcDir, thisfile))
                            if imp != None:
                                imp.setOpenAsHyperStack(False)
                                #timeseries[x].append(imp)
                                timeseries.setdefault(x, []).append(imp)
                        print "removing outliers"
                        print len(timeseries[x])
                        if not os.path.isfile(srcpath+"/goodImages.txt"):
                            good_inds = removeOutliers(timeseries[x])
                            f=open(srcpath+'/goodImages.txt','wb')
                            writer=csv.writer(f, delimiter=',')
                            for i in good_inds:
                                writer.writerow([i])
                            f.close()
                        else:
                            f=open(srcpath +'/goodImages.txt','rb')
                            reader = csv.reader(f, delimiter=',')
                            good_inds = [int(row[0]) for row in reader]
                            f.close()
                        timeseries[x] = map(timeseries[x].__getitem__, good_inds)
                        print "done removing outliers"
                        print len(timeseries[x])
                        jaimp = array(timeseries[x], ImagePlus)
                        ccc = Concatenator()
                        allimp = ccc.concatenate(jaimp, False)
                        IJ.saveAs(allimp, "TIFF", srcpath +'/'+ x)
                imp = IJ.openImage(srcpath+"/Ch2.tif")
                impch2= imp
                if not subBackground:
                    if not os.path.isfile(srcpath+"/Ch2Mean.tif"):
                        
                        from fiji.threshold import Auto_Local_Threshold as ALT    
                        #imp = IJ.getImage()
                        #IJ.run(imp, "Auto Local Threshold", "method=Bernsen radius=45 parameter_1=0 parameter_2=0 white");
                        
                        for i in range(imp.getNFrames()+1):
                            imp.setPosition(1,1,i)
                            thimp = ALT().exec(imp, "Mean", 11, 0, 0, True)
                        #IJ.run(imp, "Subtract Background...", "rolling=50 disable stack")
                        #IJ.run(imp, "Subtract...", "value=3 stack")
                        #IJ.run(imp, "Auto Local Threshold", "method=Mean radius=12 parameter_1=0 parameter_2=0 white stack")
                        #IJ.run(imp, "Invert LUT", "")
                        #IJ.run(imp, "Watershed", "stack")
                        #IJ.run(imp, "Invert LUT", "")                    
                        IJ.saveAs(imp, "TIFF", srcpath + '/Ch2Mean' )
                    imp = IJ.openImage(srcpath+'/Ch2Mean.tif')
                else:
                    if not os.path.isfile(srcpath+"/Ch2SubBack.tif"):
                        IJ.run(imp, "Subtract Background...", "rolling=50 disable stack")
                        IJ.saveAs(imp, "TIFF", srcpath + '/Ch2SubBack' )
                    else:
                        imp = IJ.openImage(srcpath+'/Ch2SubBack.tif')
                

                print('frames:'+str(imp.getNFrames()))
                impF = IJ.openImage(srcpath+"/Ch3.tif")
                impch2= IJ.openImage(srcpath+'/Ch2.tif')
                
                #colorimp = RGBStackMerge.mergeChannels(array([impF, impch2], ImagePlus), True)
                
                #IJ.saveAs(ImageConverter(colorimp).convertToRGB(), "TIFF", srcpath + '/Ch2Ch3.tif' )
                
                #----------------------------
                # Create the model object now
                #----------------------------
                    
                # Some of the parameters we configure below need to have
                # a reference to the model at creation. So we create an
                # empty model now.
                    
                model = Model()
                    
                # Send all messages to ImageJ log window.
                model.setLogger(Logger.IJ_LOGGER)
                    
                    
                       
                #------------------------
                # Prepare settings object
                #------------------------
                       
                settings = Settings()
                settings.setFrom(imp)
                
                # Configure detector - We use the Strings for the keys
                settings.detectorFactory = LogDetectorFactory()
                settings.detectorSettings = { 
                    'DO_SUBPIXEL_LOCALIZATION' : True,
                    'RADIUS' : radius,
                    'TARGET_CHANNEL' : 1,
                    'THRESHOLD' : threshold,
                    'DO_MEDIAN_FILTERING' : False,
                }  
                    
                # Configure spot filters - Classical filter on quality
                filter1 = FeatureFilter('QUALITY', qualityMin, True)
                settings.addSpotFilter(filter1)
        
        
                # Configure tracker - We want to allow merges and fusions
                settings.trackerFactory = SparseLAPTrackerFactory()
                settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap() # almost good enough
                settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = True
                settings.trackerSettings['ALLOW_TRACK_MERGING'] = True
                settings.trackerSettings['ALLOW_GAP_CLOSING'] = True
                settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = maxGapDist
                settings.trackerSettings['LINKING_MAX_DISTANCE'] = maxLinkDist
                settings.trackerSettings['MAX_FRAME_GAP'] = maxGap
                settings.addSpotAnalyzerFactory(SpotRadiusEstimatorFactory())
                
                
                # Configure track analyzers - Later on we want to filter out tracks 
                # based on their displacement, so we need to state that we want 
                # track displacement to be calculated. By default, out of the GUI, 
                # not features are calculated. 
                    
                # The displacement feature is provided by the TrackDurationAnalyzer.
                    
                settings.addTrackAnalyzer(TrackDurationAnalyzer())
                 
                # Configure track filters - We want to get rid of the two immobile spots at 
                # the bottom right of the image. Track displacement must be above 10 pixels.
                    
                filter2 = FeatureFilter('TRACK_DISPLACEMENT', 10 , True)
                settings.addTrackFilter(filter2)
                
                filter3 = FeatureFilter('TRACK_DURATION', trackDuration, True)
                settings.addTrackFilter(filter3)
                
                #print dir(tk)
                #sys.exit()
                
                #-------------------
                # Instantiate plugin
                #-------------------
                    
                trackmate = TrackMate(model, settings)
                
                #--------
                # Process
                #--------
                    
                ok = trackmate.checkInput()
                if not ok:
                    sys.exit(str(trackmate.getErrorMessage()))
                    
                ok = trackmate.process()
                if not ok:
                    sys.exit(str(trackmate.getErrorMessage()))
                print "Okay"
                 
                #----------------
                # Display results
                #----------------
                     
                #selectionModel = SelectionModel(model)
                #displayer =  HyperStackDisplayer(model, selectionModel, imp)
                #color2 = visualization.SpotColorGeneratorPerTrackFeature(model, 'SPOT_FLUORESCENCE')
                #displayer.setDisplaySettings('SpotColoring', color2)
                #displayer.render()
                #displayer.refresh()
                
                if os.path.exists(os.path.join(srcpath ,'Ch2/tracks')):
                        shutil.rmtree(os.path.join(srcpath ,'Ch2/tracks'))
                #if os.path.exists(srcpath+'tracked'):
                #        shutil.rmtree(srcpath+'tracked')
                #os.makedirs(srcpath+'tracked')
                
                fm = model.getFeatureModel()
                #print fm.getEdge
                #print fm.getEdgeFeature(1,'SPOT_SOURCE_ID')
                #print fm.getEdgeFeatureNames()
                for id in model.getTrackModel().trackIDs(True):
                    # Fetch the track feature from the feature model.
                    v = fm.getTrackFeature(id, 'TRACK_MEAN_SPEED')
                    spl = fm.getTrackFeature(id, 'TRACK_NUMBER_SPLITS') 
                    mrg = fm.getTrackFeature(id, 'NUMBER_MERGES')
                    #spl = 0
                    #mrg = 0
                    
                    model.getLogger().log('')
                    #print dir(fm)
                    #print fm.getSpotFeatureNames()
                    
                    #model.getLogger().log('Track ' + str(id) + ': mean velocity = ' + str(v) + ' ' + model.getSpaceUnits() + '/' + model.getTimeUnits()+' splits:'+str(spl) + ' merges:'+str(mrg))
                    track = model.getTrackModel().trackSpots(id)
                    edges = model.getTrackModel().trackEdges(id)
                    
                    f=open(outpath+ '/track'+ str(id) +'.txt','wb')
                    writer=csv.writer(f, delimiter='\t')
                    writer.writerow(['spotID', 'x', 'y', 't', 'q', 'snr', 'mean','fluorescence','estdiam'])
                    for spot in track:
                        sid = spot.ID()
                        # Fetch spot features directly from spot. 
                        x=spot.getFeature('POSITION_X')
                        y=spot.getFeature('POSITION_Y')
                        t=spot.getFeature('FRAME')
                        impF.setPosition(int(x),int(y),int(t))
                        spot.putFeature('FLUORESCENCE', meanInRadius(x,y,impF,radius))
                        #if t == 500:
                        #    impF.show()
                        #    asdlkjsfda
                        #print spot.getFeature('SNR')
                        q=spot.getFeature('QUALITY')
                        snr=spot.getFeature('SNR')
                        mean=spot.getFeature('MEAN_INTENSITY')
                        writer.writerow([str(sid), str(x), str(y), str(t), str(q), str(snr), str(mean),str(spot.getFeature('FLUORESCENCE')),str(spot.getFeature('ESTIMATED_DIAMETER'))])
                        
                    f.close()
                    tm = model.getTrackModel()
                    f=open(outpath+ '/edge'+ str(id) +'.txt','wb')
                    writer=csv.writer(f, delimiter='\t')
                    writer.writerow(['weight','source', 'target'])
                    for edge in edges:
                        weight= tm.getEdgeWeight(edge)
                        source= tm.getEdgeSource(edge)
                        target= tm.getEdgeTarget(edge)
                        writer.writerow([str(weight), str(source), str(target)])
                    f.close()
        
                #print(os.path.exists(os.path.expanduser('~/Desktop/GraphColorReporter.py')))
                
                #print(['python',os.path.expanduser('~/Desktop/GraphColorReporter.py'),'"'+outpath+'"','0'] )
                #f=open(os.path.expanduser('~/Desktop/'+os.path.basename(outpath)+'errorlog.txt'),'w')
                #q=open(os.path.expanduser('~/Desktop/'+os.path.basename(outpath)+'log.txt'),'w')
                #subprocess.Popen(['/w4home/mschmitz/python_libs/bin/python',os.path.expanduser('~/Desktop/GraphColorReporter.py'),'"'+outpath+'"','0'],stderr=f,stdout=q ).wait()
                #f.close()
                #q.close()
                
                #sm = SelectionModel(model)
                #trackscheme = TrackScheme(model, sm)
                #color1 = PerTrackFeatureColorGenerator(model, 'TRACK_ID')
                
                #print dir(visualization)
            
                #for attr in dir(color):
                #    if attr not in ['currentTrackID','from']:
                #        print "obj.%s = %s" % (attr, getattr(color,attr))
                
                #color2 = visualization.SpotColorGeneratorPerTrackFeature(model, 'SPOT_FLUORESCENCE')
                #for attr in dir(color):
                #    if attr not in ['currentTrackID','from']:
                #        print "obj.%s = %s" % (attr, getattr(color,attr))
                    
                
                #trackscheme.setDisplaySettings('TrackColoring', color1) 
                #trackscheme.setDisplaySettings('SpotColoring', color2)
                #trackscheme.setDisplaySettings('SpotsVisible', True)
                #trackscheme.setDisplaySettings('TracksVisible', False)
                #trackscheme.render()
                #trackscheme.refresh()
        
        # Echo results with the logger we set at start:
        model.getLogger().log(str(model))