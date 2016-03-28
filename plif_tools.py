# -*- coding: utf-8 -*-
"""
PLIF Toolset
    Created on Thu Jan 21 13:56:54 2016

    @author: Darren Banks

PLIF Temperature Calculator uses grayscale-average values from video recording
of rhodamine fluorescence to estimate the temperature field from the images.

It assumes an exponential relationship between temperature and fluorescent
intensity, as spelled out by Lemoine, et al, 'Simultaneous temperature and 2D
velocity measurements in a turblent heated jet using combined laser-induced
fluorescence and LDA', Experiments in Fluids, 26, p.315-323, 1999.
"""


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import pandas as pd
import PIL as pil
from warnings import warn


gridNum = 32
""" *Grid specifies the shape of the grid used for image analysis. Setting
    xGrid and yGrid to 10 means the images were divided into a 10x10 set of
    subdivisions. 
    
    Addition: need to verify that the gridNum is evenly divisible into the
    images being used. If not, warn the user that some pixels on the right and/
    or the bottom of the image are excluded from the analysis. """
    
experimentPath = 'H:\\test 06'
""" User-set path containing relevant experimental files. """

frameLimit = 1000


def listImages(imageDirectory, extension = '.tif'):
    """ Returns a list of images within imageDirectory that have the specified
        extension. Current camera software saves video to .tif files, so that
        is assigned as the default argument. As long as the filetype is later
        readable by PIL, any extension should work. """
    
    filePath = imageDirectory + '\\'

    imageNames = [filePath + file for file in os.listdir(imageDirectory) 
                  if extension in file]
                      
    return imageNames


def experimentDirectory(rootDirectory, primeDirectory = 'images', subDirectory = 'calibration'):
    """ Returns imageDir, a list with a primary image directory in its first
        entry, and a sub-list containing calibration directories second. It
        assumes the images are contained in a subfolder appropriately named
        'images', and that any subfolders with the word 'calibration' in their
        names contain calibration images. """
        
    if rootDirectory != primeDirectory and primeDirectory != '':
        """ If rootDirectory and primeDirectory are not the same, the prime
            should be a subfolder of root and the subDirectories contained
            within the prime. """
        
        rootPath = rootDirectory + '\\'
    
        imageDir = [rootPath + entry for entry in os.listdir(rootDirectory)
                    if os.path.isdir(rootPath + entry) 
                    if primeDirectory in entry]
        """ Creates imageDir, the directory within the rootDirectory that
            contains the name listed in primeDirectory. """
            
    elif rootDirectory == primeDirectory or primeDirectory == '':
        """ If rootDirectory and primeDirectory are the same, or primeDirectory
            is specified as an empty string, the images are located within the
            rootDirectory itself. """
        
        imageDir = [rootDirectory]

    imageDir.append([imageDir[0] + '\\'  + entry 
                     for entry in os.listdir(imageDir[0])
                     if subDirectory in entry])
    """ Appends a sublist of folders to the second entry in imageDir, which
        are the subfolders containing subDirectory name. """
                         
    return imageDir


def gridAverage(images, gridNum = 32):
    """ Using PIL/pillow, divides the input images into a gridNum x gridNum
        square grid. Then calculates the average RGB values for each of the
        cells of that grid. Returns a pandas DataFrame containing a row for
        each image in images and a column for each grid cell. """
        
    imageAverageList = []

    for image in images:
        
        currentImage = pil.Image.open(image)
    
        width = currentImage.size[0]
        xStep = width/gridNum
        xGridCoords = np.arange(0, width, xStep)
        height = currentImage.size[1]
        yStep = height/gridNum
        yGridCoords = np.arange(0, height, yStep)
        """ Based on the image's size, determine the width and coordinates of
            each grid square in pixels. """
       
        gridBoxSet = [(xCoord, yCoord, xCoord+xStep-1, yCoord+yStep-1)
                      for yCoord in yGridCoords for xCoord in xGridCoords]
        """ gridBoxSet is the collection of left, upper, right, and lower 
            bounds of each grid cell based on the image size and the desired
            number of cells. """
        
        gridAvgs = []
        """ Pre-defining gridAvgs as a list, also, clearing it out so fresh
            values are appended for each image. """
    
        for gridBox in gridBoxSet:
        
            gridBox = [int(num) for num in gridBox]
            """ gridBox is the iterating collection of each coordinate of a
                cell on the grid. The values are forced to integers because
                the image cropping function doesn't accept float values. """

            currentBox = currentImage.crop((gridBox[0:4]))
            gridAvgs.append(np.mean(currentBox))
            """ The orignal image is cropped to a single grid cell, and the
                average RGB value is added to the list of all cell values in
                the working image. """
        
        imageAverageList.append(gridAvgs)
        """ imageAverageList collects the averages of each analyzed image. """
    
#    if np.mod(width, gridNum) != 0: 
#        warn('Grid does not contain all pixels. Some will be omitted.')
#   USE THIS SECTION TO ADD GRID SQUARES TO THE EXCESS AREA - THEY WILL BE A
#   DIFFERENT SIZE THAN THE REST BUT THAT IS OKAY.
    
    imageAverages = pd.DataFrame(imageAverageList)
    
    return(imageAverages)


def getCalibrationAverages(calibImgSets, calTemperatures, gridNum=32):
    """ Using the gridAverage function, returns a DataFrame containing the RGB
        averages for each set of images denoted as calibrations. The DataFrame
        is indexed by the calibration temperatures. """
        
    calAvgList = [gridAverage(calSet, gridNum) for calSet in calibImgSets]
                       
    calAverages = pd.concat(calAvgList, keys = calTemperatures)
    
    return calAverages


def referenceIntensity(calAverages, calTemperatures):
    """ Assuming the lowest of the calibration temperatures is the reference
        value, sets that value and the grid-based RGB averages corresponding to
        images at that temperature. """
        
    refTemperature = min(calTemperatures)
    
    refAverages = calAverages.ix[refTemperature]
    
    return refAverages, refTemperature

    
def logCalibAverages(calibAvg, calTemps, refAvg, refTemp):
    """ Takes the natural logarithm of the calibration grid intensities using
        the reference temperature and intensity. """
    
    calTemps.remove(refTemp)
    """ For convenience, removes the reference temperature from the list of
        calibration temperatures. The ref temperature and corresponding image
        values are used to normalize the other calibration values; including
        them in the calculation causes problems with the logarithm. """
    
    calRatiosList = [calibAvg.ix[temp] / np.mean(refAvg) for temp in calTemps]
                                          
    calIntensityRatios = pd.concat(calRatiosList, keys = calTemps)
    
    calLogRatiosList = [np.log(calIntensityRatios.ix[key]) 
                        for key in calTemps]
    
    calLogRatios = pd.concat(calLogRatiosList, keys = calTemps)
    
    return calTemps, calLogRatios
    
def coefficientSolver(calTemperatures, refTemp, calLogRatios):
    """ for each value in calLogRatios, solve for for the coefficient of the 
        log-linear temperature-intensity curve by the following formula:
        
                C = ln(I/I_0) / ((1/T) - (1/T_ref))        
        
        calTemperatures contains the T values
        refTemperature is T_ref
        ln(I/I_0) is contained within the calLogRatios DataFrame
        """
        
    coeffValuesList = [calLogRatios.ix[temp] / ((1/temp) - (1/refTemp))
                       for temp in calTemperatures]
    
    coeffValues = pd.concat(coeffValuesList)
    
    coeffMean = pd.DataFrame(np.mean(coeffValues)).T
    
    return coeffMean


def logTestAverages(testImageAvg, refAvgs):
    """ Takes the ratio of the test image grid averages to the mean of the
        reference image grid averages. Then takes the logarithm of that ratio.
        Used to calculate the temperature. Simple but for convenience. """
        
    logTestRatios = np.log(testImageAvg / np.mean(refAvgs))
    
    return logTestRatios


def solveTemperature(betaCoeff, refTemp, logRatios):
    """ From the calculated coefficient and the logarithmic ratio averages from
        each experimental frame, returns the temperature at that location by
        this expression:
        
            T = 1 / ((1/C)ln(I/I_ref) + (1/T_ref))
        
        T is the temperature, C is the coefficient for each grid square found
        by coefficientSolver, I and I_ref are the instantaneous and reference
        RGB averages for each square, and T_ref the reference temperature.
        """
    
    invCoeff = 1 / betaCoeff
    
    invTemp = 1 / refTemp
    
    estTemperature = pd.concat([1 / ((invCoeff * row) + invTemp) 
                                for i, row in logRatios.iterrows()])
    
    return estTemperature
    

def plotTemperatures(tempValues, gridNum):
    
    plotRange = np.arange(gridNum)
    
    X, Y = np.meshgrid(plotRange, plotRange)
    
    frame = 0
    
    for index, row in tempValues.iterrows():
        
        Z = np.reshape(row, (gridNum, gridNum))
        """ Z is a height array based on the calculated temperature for a 3-D
            surface plot. It is takes the row of the temperature dataFrame and
            fits it to the x- and y-grid set on the image during analysis. """
            
        frame += 1
        
        if frame > frameLimit: break        
        
        frameTitle = 'Frame ' + str(frame)
        """ Title of each plot corresponds to its frame number in video. """
    
        fig = plt.figure()
        fig.suptitle(frameTitle)
        ax = fig.gca(projection='3d')

        ax.plot_surface(X, Y, Z, cmap = 'jet', rstride = 1, cstride = 1, 
                        vmin = 20, vmax = 100)
        """ cmap is the color map being used. 'jet' ranges from blue at minimum to
            red at maximum. rstride and cstride define the plotted grid as having
            increments of 1 in the x and y directions. vmin and vmax set the limits
            for the color map - so blue is defined as approximately 20 C and below,
            while red is 100 C and above. """                 
                    
        ax.view_init(elev = 90, azim = 0)
        ax.set_zlim3d(20, 100)
        """ view_init sets the viewing elevation as 60 degrees above the xy plane,
            and 45 degrees from orthogonal to the xy plane. zlim3d sets the limits
            on the z axis to the expected temperatures. """
        ax.set_zlabel('Deg. C')
        
        figPath = (experimentPath + '\\Figures')
        if not os.path.exists(figPath):
            os.makedirs(figPath)
        plt.savefig(figPath + '\\' + frameTitle + '.png')
        
        plt.show()
    pass