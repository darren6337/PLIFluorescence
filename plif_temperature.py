# -*- coding: utf-8 -*-
"""
PLIF Temperature Calculator
    Created on Wed Feb  3 16:51:48 2016

    @author: Darren Banks

plif_temperature calculates temperature in a plane of rhodamine-B solution
based on the intensity at which the rhodamine fluoresces under planar laser
irradiation. plif_temperature requires the module plif_tools to run.
"""

import logging
import matplotlib.pyplot as plt
import numpy as np
from os import makedirs
from os.path import exists
import pandas as pd
import plif_tools as pt


""" Logging setup """

logger = logging.getLogger('plif')
logger.setLevel(logging.DEBUG)

console_format = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s',
                                   datefmt = '%H:%M:%S')
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(console_format)
logger.addHandler(console_handler)
""" Create console handler. """

logfile = 'C:\\Users\\Darren\\Documents\\GitHub\\PLIFluorescence\\info.log'
logfile_format = logging.Formatter('%(asctime)s %(name)-24s %(levelname)-8s %(message)s', 
                                   datefmt = '%Y-%m-%d %H:%M:%S')
file_handler = logging.FileHandler(logfile)
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logfile_format)
logger.addHandler(file_handler)
""" Create logging file handler. """

logger.debug('Starting.')

plt.ioff
""" Suppressing graph output to the iPython console. """


""" Literals """

num_reference_images = 100   
""" Number of frames to establish base fluorescence within images. """

grid_number = 40
""" Number of grid cells applied to the images for analysis.
    The number of cells overall is grid_number^2. 
    """
    
want_plots = True
""" If want_plots is False, the temperature surface plots will 
    not be produced. Generally a time-saving value if False.
    """
    
plot_path = 'figures 2'
""" The folder name that will contain temperatures plots. """

plot_type = '.png'
""" Image file extension for saving results. """

results = 'temperatures 2.xlsx'
""" Name of MS Excel file to save results. """

statistics = 'statistics 2.xlsx'
""" Name of MS Excel file to save summarizing statistics. """

plot_width = 4
""" The base width in inches for output plots. """

plt.rc('font', family = 'serif', size = 24.0, serif = 'Times New Roman') 
""" Set the default font for plotting to Times New Roman, so it
    matches that used in the paper. 
    """


""" Image import """

root_directory = ('I:\\PLIF\\test 11\\images 2 - Copy')

logger.info('Directory: ' + root_directory)
""" Directory containing experiment images and calibration. """

figure_path = root_directory + '\\' + plot_path
""" Directory for result figures to be saved. """

if not exists(figure_path):
    makedirs(figure_path)

image_directories = pt.exptDirectory(root_directory, '', 'cal')

image_path = image_directories[0]
calib_paths = image_directories[1]

all_images = pt.listImages(image_path)

all_averages = pt.gridAverage(all_images, grid_number)

reference_averages = all_averages[:num_reference_images]
image_averages = all_averages[num_reference_images:]
logger.debug('First {} images used as reference'.format(num_reference_images))
""" Take the RGB mean value for the images in each grid square. """

aspect_ratio = pt.getAspectRatio(all_images[0])

logger.debug('File import complete')


""" Calibration of intensity to temperature """

mean_reference_averages = np.mean(reference_averages)
""" Take the average of each grid square over the collection of
    calibration images.
    """
    
calib_temperatures = [path[-2:] for path in calib_paths]   
int_temperatures = [int(temp) for temp in calib_temperatures]
""" Read the calibration temperatures from the calibration folder 
    names, and for convenience create a list of integer values of 
    those temperatures. 
    """
  
calib_range = [[calib_temperatures[i], calib_temperatures[i+1],
                int_temperatures[i] - int_temperatures[i+1]] 
                for i in range(len(calib_temperatures)-1)]
""" Pairs of temperatures and the difference between each pair. """

calib_image_sets = [pt.listImages(path) for path in calib_paths]
""" Gather the images located in the calibration directories. """

calib_averages = pt.getCalibrationAverages(calib_image_sets, 
                                           calib_temperatures, 
                                           grid_number)
""" Apply the grid and get RGB averages for each calibration 
    temperature. 
    """

grid_slopes_set = [np.mean(calib_averages.ix[temp[0]] - 
                           calib_averages.ix[temp[1]]) 
                   / temp[2] for temp in calib_range]
                   
grid_slopes = np.mean(pd.DataFrame(grid_slopes_set))
""" The slope is defined as the change in intensity divided by
    the change in temperature. The average delta-I / delta-T 
    for each grid square over all temperatures is used to provide
    a slope value of that square.
    """

logger.debug('Temperature calibration complete.')


""" Calculating temperature """

delta_intensity = image_averages - mean_reference_averages

delta_temperature = delta_intensity / grid_slopes

delta_temperature.to_excel(image_path+'\\temperature_deltas.xlsx')

plot_temperatures = delta_intensity / grid_slopes + int_temperatures[0]
""" Calculate the temperature based on the difference between the 
    calibration and the image's grid RGB averages. 
    """
    
plot_temperatures.to_excel(image_path + '\\' + results)
""" Save the calculated temperatures for analysis. """
    

""" Analysis, showing statistics on results. """

fig_analysis = plt.figure()

""" Maximum temperature in each frame. """
max_temperatures = pd.Series([max(T) for i, T in plot_temperatures.iterrows()])

plt.plot(max_temperatures)
plt.title('Maximum temperature per frame.', fontname = 'Times New Roman')
plt.ylabel('Deg. C')
plt.xlabel('Frame')
plt.savefig(image_path + '\\max_temperatures' + plot_type, dpi = 100)
plt.clf()

""" Average temperature in each frame. """
mean_temperatures = pd.Series([np.mean(T) 
                               for i, T in plot_temperatures.iterrows()])
                                   
plt.plot(mean_temperatures)
plt.title('Average temperature per frame.')
plt.ylabel('Deg. C')
plt.xlabel('Frame')
plt.savefig(image_path + '\\mean_temperatures' + plot_type, dpi = 100)
plt.clf()

""" Standard deviation of temperature in each frame. """
deviation_temperatures = pd.Series([np.std(T) 
                                    for i, T in plot_temperatures.iterrows()])
                                        
plt.plot(deviation_temperatures)
plt.title('Standard deviation in temperature per frame.')
plt.ylabel('Deg. C')
plt.xlabel('Frame')
plt.savefig(image_path + '\\std_dev_temperatures' + plot_type, dpi = 100)
plt.clf()

logger.info('Maximum: {}'.format(round(max(plot_temperatures.max()))))
logger.info('Minimum: {}'.format(round(min(plot_temperatures.min()))))
logger.info('Median: {}'.format(round(np.median(plot_temperatures.median()))))
logger.debug('{} frames to be analyzed'.format(len(plot_temperatures)))

if min(plot_temperatures.min()) < 25:
    logger.warn('Subcooled, possibly erroneous temperatures found')
""" Report the temperature statistics. """

thermal_stats_list = [max_temperatures, mean_temperatures, 
                      deviation_temperatures]
thermal_statistics = pd.concat(thermal_stats_list,
                               keys= ['max', 'mean', 's.dev'])
thermal_statistics.to_frame().to_excel(image_path + '\\' + statistics)

    
""" Plot temperature contour. """

if want_plots:
    
    z_minimum = 25
    z_maximum = 100
    """ User sets the graph maximum and minimum temperature values. """

    plot_range = np.arange(grid_number)
    x_grid, y_grid = np.meshgrid(plot_range, plot_range)
    """ Setting up the X and Y array for plotting purposes. """

    temperature_intervals = np.arange(z_minimum, z_maximum, 1)
    """ The temperature range for the graph to use in scaling its
        color map.
        """

    fig = plt.figure(figsize = (2.5*plot_width, 2.0*plot_width/aspect_ratio))

    for index, row in plot_temperatures.iterrows():

        frame_title = 'Frame {}'.format(index-99)
        """ Title of each plot corresponds to its frame number in video. """
    
        plot_temperature_array = np.reshape(row, (grid_number, grid_number))
        """ plotTemperatureArray is the calculated temperature for a 
            3-D surface plot. It takes the row of the temperature 
            dataFrame and fits it to the x- and y-grid set on the 
            image during analysis. 
            """

        plt.contourf(x_grid, y_grid, plot_temperature_array, 
                     temperature_intervals, cmap = 'jet', extend='both',
                     vmin = z_minimum, vmax = z_maximum)            
        plt.title(frame_title)
        plt.xticks(np.arange(0, grid_number, 1))
        plt.yticks(np.arange(0, grid_number, 1))
        plt.colorbar()
        plt.grid(color = 'k', linestyle = 'solid', which='both')
        """ Creating and formatting the plot with a colormap, the 
            previously set Z limits, ticks with intervals of 1, 
            and a black grid. 
            """

        """ Save the figure within a subfolder of the initial 
            directory, and then clear the figure. 
            """
        plt.savefig(figure_path + '\\' + frame_title + plot_type, dpi = 50)
        plt.clf()
    
        if np.mod(index-99, 100) == 0:
            logger.debug('Frame {} graphed'.format(index-99))
        """ Iterating over the frames. """
        
plt.close()
plt.close()

logger.info('Complete\n')