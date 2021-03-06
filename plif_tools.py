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


import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import PIL as pil


logger = logging.getLogger('plif.tools')


def listImages(image_directory, extension='.tif'):
    """ Returns a list of images within image_directory that have the
        specified extension.

        Current camera software saves video to '.tif' files, hence the
        default argument. As long as files are readable by PIL, any
        extension should work.
        """

    file_path = image_directory + '\\'

    image_names = [file_path + file for file in os.listdir(image_directory)
                   if extension in file]

    num_images = len(image_names)

    logger = logging.getLogger('plif.tools.listImages')
    if len(image_directory) > 16:
        short_dir = '...' + image_directory[-16:]
    else:
        short_dir = image_directory
    logger.debug('{} images loaded from {}'.format(num_images, short_dir))

    return image_names


def exptDirectory(root_directory, prime_directory='images',
                  sub_directory='calibration'):
    """ Returns image_dir, a list with a primary image directory in
        its first entry, and a sub-list containing calibration
        directories second. It assumes the images are contained in
        a subfolder appropriately named 'images', and that any
        subfolders with the word 'calibration' in their names
        contain calibration images.
        """

    if root_directory != prime_directory and prime_directory != '':
        """ If root_directory and prime_directory are not the same,
            the prime should be a subfolder of root and the
            subDirectories contained within the prime.
            """

        root_path = root_directory + '\\'

        image_dir = [root_path + entry for entry in os.listdir(root_directory)
                     if os.path.isdir(root_path + entry)
                     if prime_directory in entry]
        """ Creates image_dir, the directory within the root_directory
            that contains the name listed in prime_directory.
            """

    elif root_directory == prime_directory or prime_directory == '':
        """ If root_directory and prime_directory are the same, or
            prime_directory is specified as an empty string, the
            images are located within the root_directory itself.
            """

        image_dir = [root_directory]

    image_dir.append([image_dir[0] + '\\' + entry
                     for entry in os.listdir(image_dir[0])
                     if sub_directory in entry])
    """ Appends a sublist of folders to the second entry in image_dir,
        which are the subfolders containing subDirectory name.
        """

    return image_dir


def gridAverage(images, grid_num=32):
    """ Using PIL/pillow, divides the input images into a gridNum x
        gridNum square grid. Then calculates the average RGB values
        for each of the cells of that grid. Returns a pandas
        DataFrame containing a row for each image in images and a
        column for each grid cell.
        """

    image_average_list = []

    for image in images:

        current_image = pil.Image.open(image)

        width, height = current_image.size[:]
        x_step = width / grid_num
        x_coords = np.arange(0, width, x_step)
        y_step = height / grid_num
        y_coords = np.arange(0, height, y_step)
        """ Based on the image's size, determine the width and
            coordinates of each grid square in pixels.
            """

        grid_set = [(x_coord, y_coord, x_coord + x_step - 1,
                     y_coord + y_step - 1)
                    for y_coord in y_coords for x_coord in x_coords]
        """ gridBoxSet is the collection of left, upper, right,
            and lower bounds of each grid cell based on the image
            size and the desired number of cells.
            """

        grid_averages = []
        """ Pre-defining gridAvgs as a list, also, clearing it out
            so fresh values are appended for each image.
            """

        for grid_box in grid_set:

            grid_box = [int(num) for num in grid_box]
            """ gridBox is the iterating collection of each
                coordinate of a cell on the grid. The values are
                forced to integers because the image cropping
                function doesn't accept float values.
                """

            current_box = current_image.crop((grid_box[0:4]))
            grid_averages.append(np.mean(current_box))
            """ The orignal image is cropped to a single grid cell,
                and the average RGB value is added to the list of
                all cell values in the working image.
                """

        image_average_list.append(grid_averages)
        """ imageAverageList collects averages of each image. """

    image_averages = pd.DataFrame(image_average_list)

    return(image_averages)


def getCalibrationAverages(calib_images, cal_temperatures, grid_num=32):
    """ Using the gridAverage function, returns a DataFrame containing
        the RGB averages for each set of images denoted as
        calibrations. The DataFrame is indexed by the calibration
        temperatures.
        """

    cal_avg_list = [gridAverage(cal_set, grid_num) for cal_set in calib_images]

    cal_averages = pd.concat(cal_avg_list, keys=cal_temperatures)

    return cal_averages


def getGridSlopes(cal_averages, cal_temperatures):

    calib_range = [[cal_temperatures[i], cal_temperatures[i+1],
                    int(cal_temperatures[i]) - int(cal_temperatures[i+1])]
                   for i in range(len(cal_temperatures)-1)]

    grid_slopes_set = [np.mean(cal_averages.ix[temp[0]] -
                       cal_averages.ix[temp[1]]) / temp[2]
                       for temp in calib_range]

    grid_slopes = np.mean(pd.DataFrame(grid_slopes_set))

    return(grid_slopes)


def getAspectRatio(image_path, decimal_point=1):
    """ Returns the input image's aspect ratio.

        (width/height) rounded to a default of 1 decimal point.
        """

    reference_image = pil.Image.open(image_path)

    image_width, image_height = reference_image.size

    aspect_ratio = round(image_width/image_height, decimal_point)

    logger = logging.getLogger('plif.tools.aspectRatio')
    logger.debug('Aspect ratio {}'.format(aspect_ratio))

    return aspect_ratio


def getTemperatureStats(temperatures, image_path, statistic='statistics.xlsx'):
    """ Returns stats_list, a list of the maximum, average, and standard
        deviation of temperatures within each row of the temperatures
        DataFrame.
    """

    logger = logging.getLogger('plif.tools.getTempStats')

    temperature_rows = temperatures.iterrows()

    max_temperatures = pd.Series([max(T) for i, T in temperature_rows])
    mean_temperatures = pd.Series([np.mean(T) for i, T in temperature_rows])
    stddev_temperatures = pd.Series([np.std(T) for i, T in temperature_rows])

    logger.info('Maximum: {}'.format(round(max(temperatures.max()))))
    logger.info('Minimum: {}'.format(round(min(temperatures.min()))))
    logger.info('Median: {}'.format(round(np.median(temperatures.median()))))
    logger.debug('{} frames to be analyzed'.format(len(temperatures)))

    stats_list = [max_temperatures, mean_temperatures, stddev_temperatures]
    thermal_statistics = pd.concat(stats_list, keys=['max', 'mean', 's.dev'])
    thermal_statistics.to_frame().to_excel(image_path + '\\' + statistic)

    return(stats_list)


def plotTemperatureStats(stats_list, image_path, plot_type='.png'):
    """ Using the stats_list produced by getTemperatureStats, plots the
        maximum, mean, and standard deviation of temperatures in each row
        of the temperatures in the DataFrame fed to getTemperatureStats.
    """

    plt.plot(stats_list[0])
    plt.title('Maximum temperature per frame.', fontname='Times New Roman')
    plt.ylabel('Deg. C')
    plt.xlabel('Frame')
    plt.savefig(image_path + '\\max_temperatures' + plot_type, dpi=100)
    plt.clf()

    plt.plot(stats_list[1])
    plt.title('Average temperature per frame.')
    plt.ylabel('Deg. C')
    plt.xlabel('Frame')
    plt.savefig(image_path + '\\mean_temperatures' + plot_type, dpi=100)
    plt.clf()

    plt.plot(stats_list[2])
    plt.title('Standard deviation in temperature per frame.')
    plt.ylabel('Deg. C')
    plt.xlabel('Frame')
    plt.savefig(image_path + '\\std_dev_temperatures' + plot_type, dpi=100)
    plt.clf()

    pass
