# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:51:48 2016

@author: Darren
"""

import matplotlib.pyplot as plt

import numpy as np
import os
import pandas as pd
import PIL as pil
import plif_tools
import winsound

plt.ioff
""" Suppressing graph output to the iPython console. There liekly will be way
    too many graphs produced for that to be practical. """

frameLimit = 500
""" Maximum number of frames to iterate over. """
referenceImages = 100   
""" Number of frames used to establish baseline fluorescence
    at each location within the images. """
frame = referenceImages
""" A counting variable used while iterating and plotting. No plots are
    produced for calibration frames, so the iteration starts at the first 
    non-calibration frame. """

sound = True
""" A beep is emitted upon request for input and upon script completion. """

gridNum = 40
""" The number of grid cells applied to the images for analysis in the 
    x and y directions. The number of cells overall is gridNum^2. """

imgpath = ('I:\\PLIF\\test 11\\images 2')
""" Directory containing images for analysis. """
calpath = ('I:\\PLIF\\test 11\\images 2')

refimgs = plif_tools.listImages(imgpath)[:referenceImages]
""" The first few images are used as unheated references. """

imgs = plif_tools.listImages(imgpath)[referenceImages:frameLimit]
""" The images taken during heating to be analyzed. """

imgavgs = plif_tools.gridAverage(imgs, gridNum)
refavgs = plif_tools.gridAverage(refimgs, gridNum)
""" Taking the RGB mean value for the collections of images with regard to
    each grid square. """

refavgsavg = np.mean(refavgs)
""" Take the average of each grid square over the collection of calibration
    images. """
    
caltemps = ['25', '30', '35', '40', '45', '50']
calpaths = [calpath + '\\' + temp for temp in caltemps]

calrange = []
i = 0
while i < len(caltemps)-1:
    calrange.append([int(caltemps[i]), int(caltemps[i+1]),
                     int(caltemps[i]) - int(caltemps[i+1])])
    i+=1

calimgsets = [plif_tools.listImages(path) for path in calpaths]

calavgs = plif_tools.getCalibrationAverages(calimgsets, caltemps, gridNum)

slopes = [np.mean((calavgs.ix[str(p[0])] - calavgs.ix[str(p[1])])) / p[2] 
                      for p in calrange]
slopes = pd.DataFrame(slopes)
    
#slopes1 = np.mean((calavgs.ix[str(calrange1[0])] - calavgs.ix[str(calrange1[1])]/(calrange1[1] - calrange1[0])))
#slopes2 = np.mean((calavgs.ix[str(calrange2[0])] - calavgs.ix[str(calrange2[1])]/(calrange2[1] - calrange2[0])))
#slopes3 = np.mean((calavgs.ix[str(calrange3[0])] - calavgs.ix[str(calrange3[1])]/(calrange3[1] - calrange3[0])))
#slopes4 = np.mean((calavgs.ix[str(calrange4[0])] - calavgs.ix[str(calrange4[1])]/(calrange4[1] - calrange4[0])))
#slopes5 = np.mean((calavgs.ix[str(calrange5[0])] - calavgs.ix[str(calrange5[1])]/(calrange5[1] - calrange5[0])))
#
#slopesset = [slopes2, slopes3, slopes4, slopes1]
#
#slopescon = pd.concat(slopesset)
#
#slopes = np.mean(slopesset)

deltaintensity = refavgsavg - imgavgs

deltaintensity[deltaintensity<0] = 0

scale = 1

#plottemps = imgavgs / refavgsavg
plottemps = slopes*scale*deltaintensity + int(caltemps[0])
""" Calculates the temperature based on the difference between the calibration
    and the image's grid RGB averages. Need to use temperature calibrations to
    estimate how this differenc relates to temperature. """

if sound:
    winsound.Beep(600, 500)
    
#looper = 1
#
#while looper < 100:
#    print('Maximum value: {}'.format(max(plottemps.max())))
#    print('Minimum value: {}'.format(min(plottemps.min())))
#    print('Median value: {}'.format(np.median(plottemps.median())))
#    """ Report the statistics of the calculated temperatures so that graph limits
#        can be set by the user. """
#
#    scale = float(input('scale? '))
#    
#    plottemps = plottemps * scale
#
#    looper = int(input('Loop >100 will escape.'))


print('Maximum value: {}'.format(max(plottemps.max())))
print('Minimum value: {}'.format(min(plottemps.min())))
print('Median value: {}'.format(np.median(plottemps.median())))

print('{} frames'.format(len(imgs)))
#zmin = float(input('Min graph? '))
#zmax = float(input('Max graph? '))
zmin = 25
zmax = 250
""" User sets the graph maximum and minimum temperature values. """

if zmin > zmax:
    re_zmin = zmax
    re_zmax = zmin
    zmin = re_zmin
    zmax = re_zmax
""" Might be superfluous but this corrects if the user sets the upper limit of
    the graph to a lower value than the lower limit. During debugging this
    happened to me and wasted a lot of computing time. """

plotRange = np.arange(gridNum)
X, Y = np.meshgrid(plotRange, plotRange)
""" Setting up the X and Y array for plotting purposes. """

levels = np.arange(zmin, zmax, 1)
""" The temperature range for the graph to use in scaling its color map. """

fig = plt.figure(figsize = (4, 6))
""" Opening a figure. Currently the aspect ratio is set to a portrait value,
    reflecting the general aspect of the images on video. """

for index, row in (plottemps).iterrows():

    Z = np.reshape(row, (gridNum, gridNum))
    """ Z is a height array based on the calculated temperature for a 3-D
        surface plot. It is takes the row of the temperature dataFrame and fits
        it to the x- and y-grid set on the image during analysis. """

    frame += 1
    if np.mod(frame, 100) == 0:
        print('Frame {}'.format(frame))
    if frame > frameLimit: 
        print('Frame limit reached.')        
        break
    """ Iterating over the frames, and escaping the loop if the frame limit
        is exceeded. """

    frameTitle = 'Frame ' + str(frame-101)
    """ Title of each plot corresponds to its frame number in video. """

    plt.contourf(X, Y, Z, levels, cmap = 'jet', vmin = zmin, vmax = zmax, 
                 antialiased = False, extend='both')            
    plt.title(frameTitle)
    plt.xticks(np.arange(0,gridNum,1))
    plt.yticks(np.arange(0,gridNum,1))
    plt.colorbar()
    plt.grid(color = 'k', linestyle = 'solid', which='both')
    """ Creating and formatting the plot with a colormap, the previously set
        Z limits, ticks with intervals of 1, and a black grid. """

    plt.savefig(imgpath + '\\figures\\' + frameTitle + '.png', dpi=100)
    plt.clf()
    """ Save the figure within a subfolder of the initial directory, and then
        clear the figure for the next plot. """

if sound:
    winsound.Beep(600, 500)