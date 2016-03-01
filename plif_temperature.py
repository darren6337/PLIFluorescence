# -*- coding: utf-8 -*-
"""
PLIF Temperature Calculator
    Created on Wed Feb  3 16:51:48 2016

    @author: Darren Banks

plif_temperature calculates temperature in a plane of rhodamine-B solution
based on the intensity at which the rhodamine fluoresces under planar laser
irradiation. plif_temperature requires the module plif_tools to run.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import PIL as pil
import plif_tools as pt
from time import time
from winsound import Beep


timeZero = time()
""" Simple profiling. """

plt.ioff
""" Suppressing graph output to the iPython console. There likely will be way
    too many graphs produced for that to be practical or time-efficient. """


""" LITERALS AND TEST VALUES """

frameLimit = 5000
""" Maximum number of frames to iterate over. """

numberReferenceImages = 100   
""" Number of frames used to establish baseline fluorescence at each location
    within the images. """
    
frameNumber = numberReferenceImages
""" A counting variable used while iterating and plotting. No plots are
    produced for calibration frames, so the iteration starts at the first 
    non-calibration frame. """
    
soundOn = True
""" A beep is emitted upon request for input and upon script completion. """

gridNumber = 40
""" The number of grid cells applied to the images for analysis in the 
    x and y directions. The number of cells overall is gridNumber^2. """
    
plotLimitRequest = False
""" If plotLimitRequest is false, the code uses preset values of 25 and 250 as
    the minimum and maximum values of the z-axis when plotting results. """
    
plotPath = 'figures'
""" The folder name that will contain plotted temperatures after code runs. """

plotType = '.png'
""" Image file extension for saving results. """

resultsName = 'temperatures.xlsx'
""" Name of MS Excel file to save results. """

plotWidth = 4
""" The base width for output plots to be saved. Units of inches. """


""" IMPORT IMAGES """

rootDirectory = ('I:\\PLIF\\test 12b')
""" Directory containing experiment images and calibration. """

figurePath = rootDirectory + '\\' + plotPath
""" Directory for result figures to be saved. """

if not os.path.exists(figurePath):
    os.makedirs(figurePath)

imageDirectories = pt.experimentDirectory(rootDirectory, '', 'cal')

imagePath = imageDirectories[0]
calibPaths = imageDirectories[1]

def imageImport(imagePath, gridNumber):
    """ Returns a list of all image objects with the path and the RGB averages
        of the gridNumber by gridNumber grid applied to those images. """
        
    allImages = pt.listImages(imagePath)

    allImageAverages = pt.gridAverage(allImages, gridNumber)
    
    return allImages, allImageAverages

    
allImages, allAverages = imageImport(imagePath, gridNumber)

referenceAverages = allAverages[:numberReferenceImages]
imageAverages = allAverages[numberReferenceImages:]
""" Take the RGB mean value for the collections of images with regard to
    each grid square. """


def getAspectRatio(imagePath, decimalPoint = 1):
    """ Returns the input image's aspect ratio (width/height) rounded to a
        default of 1 decimal point. """
        
    referenceImage = pil.Image.open(imagePath)

    imageWidth, imageHeight = referenceImage.size

    aspectRatio = round(imageWidth/imageHeight, 1)
    
    return aspectRatio
 
   
aspectRatio = getAspectRatio(allImages[0])

timeImport = round(time() - timeZero, 0)
print('Import time: {} s.'.format(timeImport))
""" Profiling. """


""" REFERENCE VALUES AND CALIBRATION """

meanReferenceAverages = np.mean(referenceAverages)
""" Take the average of each grid square over the collection of calibration
    images. """
    
calibTemperatures = [path[-2:] for path in calibPaths]   
intTemperatures = [int(temp) for temp in calibTemperatures]
""" Read the calibration temperatures from the calibration folder names, and
    for convenience create a list of integer values of those temperatures. """
  
calibRange = [[calibTemperatures[i], calibTemperatures[i+1],
               intTemperatures[i] - intTemperatures[i+1]] 
               for i in range(len(calibTemperatures)-1)]
""" Pairs of temperatures and the difference between each pair. """

calibImageSets = [pt.listImages(path) for path in calibPaths]
""" Gather the images located in the calibration directories. """

calibAverages = pt.getCalibrationAverages(calibImageSets, calibTemperatures, 
                                          gridNumber)
""" Apply the grid and get RGB averages for each calibration temperature. """

gridSlopesSet = [np.mean(calibAverages.ix[temp[0]] - calibAverages.ix[temp[1]]) 
                 / temp[2] for temp in calibRange]                         
gridSlopes = np.mean(pd.DataFrame(gridSlopesSet))
""" The slope is defined as the change in intensity divided by the change in
    temperature. The average delta-I / delta-T for each grid square over all
    temperatures is used to provide a slope value of that square. """

timeCalibration = round(time() - timeImport - timeZero, 0)
print('Calibration time: {} s.'.format(timeCalibration))
""" Profiling. """


""" CALCULATING TEMPERATURE """

deltaIntensity = imageAverages - meanReferenceAverages

plotTemperatures = deltaIntensity / gridSlopes + intTemperatures[0]
""" Calculate the temperature based on the difference between the calibration
    and the image's grid RGB averages. """
    
plotTemperatures.to_excel(imagePath + '\\' + resultsName)
""" Save the calculated temperatures for analysis. """

if soundOn:
    Beep(600, 500)
    

""" ANALYSIS """

figAnalysis = plt.figure()

""" Plotting and saving the maximum temperature in each frame. """
maxTemperatures = pd.Series([max(T) for i, T in plotTemperatures.iterrows()])
plt.plot(maxTemperatures)
plt.title('Maximum temperature per frame.')
plt.savefig(imagePath + '\\max_temperatures' + plotType, dpi = 100)
plt.clf()

""" Plotting and saving the average temperature in each frame. """
meanTemperatures = pd.Series([np.mean(T) 
                              for i, T in plotTemperatures.iterrows()])
plt.plot(meanTemperatures)
plt.title('Average temperature per frame.')
plt.savefig(imagePath + '\\mean_temperatures' + plotType, dpi = 100)
plt.clf()

""" Plotting and saving standard deviation of temperature in each frame. """
deviationTemperatures = pd.Series([np.std(T) 
                                   for i, T in plotTemperatures.iterrows()])
plt.plot(deviationTemperatures)
plt.title('Standard deviation in temperature per frame.')
plt.savefig(imagePath + '\\std_dev_temperatures' + plotType, dpi = 100)
plt.clf()

print('Maximum value: {}'.format(max(plotTemperatures.max())))
print('Minimum value: {}'.format(min(plotTemperatures.min())))
print('Median value: {}'.format(np.median(plotTemperatures.median())))
print('{} frames'.format(len(imageAverages)))
""" Report the statistics for the user to be able to set the plot limits. """

    
""" PLOTTING TEMPERATURE """

if plotLimitRequest:
    zMinimum = float(input('Min graph? '))
    zMaximum = float(input('Max graph? '))
else:
    zMinimum = 25
    zMaximum = 100
""" User sets the graph maximum and minimum temperature values. """

if zMinimum > zMaximum:
    zMinimum, zMaximum = zMaximum, zMinimum
""" Corrects if the upper limit of the graph is lower than the lower limit. """

plotRange = np.arange(gridNumber)
xGrid, yGrid = np.meshgrid(plotRange, plotRange)
""" Setting up the X and Y array for plotting purposes. """

temperatureIntervals = np.arange(zMinimum, zMaximum, 1)
""" The temperature range for the graph to use in scaling its color map. """

fig = plt.figure(figsize = (plotWidth, plotWidth/aspectRatio))

for index, row in (plotTemperatures).iterrows():

    frameTitle = 'Frame {}'.format(frameNumber - 1)
    """ Title of each plot corresponds to its frame number in video. """
    
    plotTemperatureArray = np.reshape(row, (gridNumber, gridNumber))
    """ plotTemperatureArray is the calculated temperature for a 3-D surface
        plot. It takes the row of the temperature dataFrame and fits it to the
        x- and y-grid set on the image during analysis. """

    plt.contourf(xGrid, yGrid, plotTemperatureArray, temperatureIntervals, 
                 cmap = 'jet', vmin = zMinimum, vmax = zMaximum, extend='both')            
    plt.title(frameTitle)
    plt.xticks(np.arange(0, gridNumber, 1))
    plt.yticks(np.arange(0, gridNumber, 1))
    plt.colorbar()
    plt.grid(color = 'k', linestyle = 'solid', which='both')
    """ Creating and formatting the plot with a colormap, the previously set
        Z limits, ticks with intervals of 1, and a black grid. """

    plt.savefig(figurePath + '\\' + frameTitle + plotType, dpi = 100)
    plt.clf()
    
    frameNumber += 1
    if np.mod(frameNumber, 100) == 0:
        print('Frame {}'.format(frameNumber))
    if frameNumber > frameLimit: 
        print('Frame limit reached.')        
        break
    """ Iterating over the frames. Ends the loop if past the frame limit. """
    """ Save the figure within a subfolder of the initial directory, and then
        clear the figure for the next plot. """


""" CLEANUP """

timeComplete = round(time() - timeZero, 0)
print('Completed in {} s.'.format(timeComplete))
""" Profiling. """

if soundOn:
    Beep(600, 500)


""" END OF FILE """