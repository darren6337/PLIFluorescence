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
import plif_tools as pt
import sys


""" Logging setup """

logger = logging.getLogger('plif')
logger.setLevel(logging.DEBUG)

con_format = '%(asctime)s - %(name)s -  %(levelname)-8s: %(message)s'
console_format = logging.Formatter(con_format, datefmt='%H:%M:%S')

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(console_format)
logger.addHandler(console_handler)
""" Create console handler. """

log_file = 'C:\\Users\\Darren\\Documents\\GitHub\\PLIFluorescence\\debug.log'
if not exists(log_file):
    file = open(log_file, 'a')
    file.close()

log_format = '%(asctime)s %(name)-24s %(levelname)-8s %(message)s'
logfile_format = logging.Formatter(log_format, datefmt='%Y-%m-%d %H:%M:%S')

file_handler = logging.FileHandler(log_file)
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logfile_format)
logger.addHandler(file_handler)
""" Create debug logging file handler. """

info_file = 'C:\\Users\\Darren\\Documents\\GitHub\\PLIFluorescence\\info.log'
if not exists(info_file):
    file = open(info_file, 'a')
    file.close()

info_handler = logging.FileHandler(info_file, mode='w')
info_handler.setLevel(logging.INFO)
info_handler.setFormatter(logfile_format)
logger.addHandler(info_handler)
""" Creating info logging file handler. """

logger.debug('Starting.')

plt.ioff
""" Suppressing graph output to the iPython console. """


""" Literals """

num_reference_images = 100
""" Number of frames to establish base fluorescence within images. """

grid_number = 40
""" Number of grid cells applied to the images for analysis. """

want_plots = False
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

plt.rc('font', family='serif', size=24.0, serif='Times New Roman')
""" Set the default font for plotting to Times New Roman, so it
    matches that used in the paper.
    """


""" Image import """

root_directory = ('I:\\PLIF\\test 11\\images 2 - Copy')

if not exists(root_directory):
    logger.error('Experiment directory does not exist!')
    sys.exit()

logger.info('Directory: ' + root_directory)
""" Directory containing experiment images and calibration. """

figure_path = root_directory + '\\' + plot_path
""" Directory for result figures to be saved. """

if not exists(figure_path):
    makedirs(figure_path)

[image_path, calib_paths] = pt.exptDirectory(root_directory, '', 'cal')

all_images = pt.listImages(image_path)

all_averages = pt.gridAverage(all_images, grid_number)

reference_averages = all_averages[:num_reference_images]
image_averages = all_averages[num_reference_images:]
logger.debug('First {} images used as reference'.format(num_reference_images))
""" Take the RGB mean value for the images in each grid square. """

aspect_ratio = pt.getAspectRatio(all_images[0])

logger.info('File import complete')


""" Calibration of intensity to temperature """

mean_reference_averages = np.mean(reference_averages)
""" Take the average of each grid square over the collection of
    calibration images.
    """

calib_temperatures = [path[-2:] for path in calib_paths]

calib_image_sets = [pt.listImages(path) for path in calib_paths]
""" Gather the images located in the calibration directories. """

calib_averages = pt.getCalibrationAverages(calib_image_sets,
                                           calib_temperatures, grid_number)
""" Apply grid and get RGB averages for each calibration temperature. """

grid_slopes = pt.getGridSlopes(calib_averages, calib_temperatures)

logger.info('Temperature calibration complete.')


""" Calculating temperature """

delta_intensity = image_averages - mean_reference_averages

delta_temperature = delta_intensity / grid_slopes

delta_temperature.to_excel(image_path+'\\temperature_deltas.xlsx')

plot_temperatures = delta_intensity / grid_slopes + int(calib_temperatures[0])
""" Calculate the temperature based on the difference between the
    calibration and the image's grid RGB averages.
    """

if min(plot_temperatures.min()) < 25:
    logger.warn('Subcooled, possibly erroneous temperatures')


plot_temperatures.to_excel(image_path + '\\' + results)
""" Save the calculated temperatures for analysis. """


""" Reporting the temperature statistics. """

stats_list = pt.getTemperatureStats(plot_temperatures, image_path, statistics)

pt.plotTemperatureStats(stats_list, image_path, plot_type)


""" Plotting temperature contour in each video frame. """

if want_plots:

    z_minimum = 25
    z_maximum = 100
    """ User sets the graph maximum and minimum temperature values. """

    plot_range = np.arange(grid_number)
    x_grid, y_grid = np.meshgrid(plot_range, plot_range)
    """ Setting up the X and Y array for plotting purposes. """

    temperature_intervals = np.arange(z_minimum, z_maximum, 1)
    """ The temperature range to scale the color map. """

    fig = plt.figure(figsize=(2.5*plot_width, 2.0*plot_width/aspect_ratio))

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
                     temperature_intervals, cmap='jet', extend='both',
                     vmin=z_minimum, vmax=z_maximum)
        plt.title(frame_title)
        plt.xticks(np.arange(0, grid_number, 1))
        plt.yticks(np.arange(0, grid_number, 1))
        plt.colorbar()
        plt.grid(color='k', linestyle='solid', which='both')
        """ Creating and formatting the plot with a colormap, the
            previously set Z limits, ticks with intervals of 1,
            and a black grid.
            """

        """ Save the figure within a subfolder of the initial
            directory, and then clear the figure.
            """
        plt.savefig(figure_path + '\\' + frame_title + plot_type, dpi=50)
        plt.clf()

        if np.mod(index-99, 100) == 0:
            logger.debug('Frame {} graphed'.format(index-99))
        """ Iterating over the frames. """

plt.close('all')

if not want_plots:
    logger.info('Temperatures not plotted.')

logger.info('Complete\n')
