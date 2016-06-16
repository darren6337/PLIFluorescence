# -*- coding: utf-8 -*-
"""
Storage script for old/unused plif_tools code

Created on Tue Jun 14 10:29:17 2016

@author: Darren
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import pandas as pd


grid_num = 32
""" *Grid specifies the shape of the grid used for image analysis. 
    Setting xGrid and yGrid to 10 means the images were divided
    into a 10x10 set of subdivisions. 
    
    Addition: need to verify that the gridNum is evenly divisible
    into the images being used. If not, warn the user that some
    pixels on the right and/or the bottom of the image are excluded 
    from the analysis. 
    """
    
experiment_path = 'H:\\test 06'
""" User-set path containing relevant experimental files. """

frame_limit = 1000

def referenceIntensity(cal_averages, cal_temperatures):
    """ Assuming the lowest of the calibration temperatures is the 
        reference value, sets that value and the grid-based RGB 
        averages corresponding to images at that temperature.
        """        

    reference_temperature = min(cal_temperatures)
    
    reference_averages = cal_averages.ix[reference_temperature]
    
    return reference_averages, reference_temperature
    """ UNUSED """

    
def logCalibAverages(calib_average, cal_temps, ref_average, ref_temperature):
    """ Takes the natural logarithm of the calibration grid 
        intensities using the reference temperature and intensity. 
        """
    
    cal_temps.remove(ref_temperature)
    """ For convenience, removes the reference temperature from the 
        list of calibration temperatures. The ref temperature and 
        corresponding image values are used to normalize the other 
        calibration values; including them in the calculation causes
        problems with the logarithm. 
        """
    
    cal_ratios_list = [calib_average.ix[temp] / np.mean(ref_average) 
                       for temp in cal_temps]
                                          
    cal_intensity_ratios = pd.concat(cal_ratios_list, keys = cal_temps)
    
    cal_log_ratios_list = [np.log(cal_intensity_ratios.ix[key]) 
                           for key in cal_temps]
    
    cal_log_ratios = pd.concat(cal_log_ratios_list, keys = cal_temps)
    
    return cal_temps, cal_log_ratios
    """ DEPRECATED """
    
def coefficientSolver(cal_temperatures, ref_temperature, cal_log_ratios):
    """ for each value in calLogRatios, solve for for the 
        coefficient of the log-linear temperature-intensity curve
        by the following formula:
        
                C = ln(I/I_0) / ((1/T) - (1/T_ref))        
        
        calTemperatures contains the T values
        refTemperature is T_ref
        ln(I/I_0) is contained within the calLogRatios DataFrame
        """
        
    coeff_values_list = [cal_log_ratios.ix[temp] / 
                         ((1 / temp) - (1 / ref_temperature))
                         for temp in cal_temperatures]
    
    coeff_values = pd.concat(coeff_values_list)
    
    coeff_mean = pd.DataFrame(np.mean(coeff_values)).T
    
    return coeff_mean
    """DEPRECATED"""


def logTestAverages(test_image_averages, ref_averages):
    """ Takes the ratio of the test image grid averages to the 
        mean of the reference image grid averages. Then takes 
        the logarithm of that ratio. Used to calculate the 
        temperature. Simple but for convenience.
        """
        
    log_test_ratios = np.log(test_image_averages / np.mean(ref_averages))
    
    return log_test_ratios
    """ DEPRECATED """


def solveTemperature(beta_coeff, ref_temperature, log_ratios):
    """ From the calculated coefficient and the logarithmic ratio
        averages from each experimental frame, returns the 
        temperature at that location by this expression:
        
            T = 1 / ((1/C)ln(I/I_ref) + (1/T_ref))
        
        T is the temperature, C is the coefficient for each grid 
        square found by coefficientSolver, I and I_ref are the 
        instantaneous and reference RGB averages for each square, 
        and T_ref the reference temperature.
        """
    
    inv_coeff = 1 / beta_coeff
    
    inv_temperature = 1 / ref_temperature
    
    est_temperature = pd.concat([1 / ((inv_coeff * row) + inv_temperature) 
                                for index, row in log_ratios.iterrows()])
    
    return est_temperature
    """" DEPRECATED """
    

def plotTemperatures(temperatures, grid_num):
    
    plot_range = np.arange(grid_num)
    
    x, y = np.meshgrid(plot_range, plot_range)
    
    frame = 0
    
    for index, row in temperatures.iterrows():
        
        z = np.reshape(row, (grid_num, grid_num))
        """ Z is a height array based on the calculated temperature
            for a 3-D surface plot. It is takes the row of the 
            temperature dataFrame and fits it to the x- and y-grid 
            set on the image during analysis.
            """
            
        frame += 1
        
        if frame > frame_limit: break    # uhhhhhh     
        
        frame_title = 'Frame ' + str(frame)
        """ Title of each plot corresponds to its frame number in video. """
    
        fig = plt.figure()
        fig.suptitle(frame_title)
        ax = fig.gca(projection='3d')

        ax.plot_surface(x, y, z, cmap = 'jet', rstride = 1, cstride = 1, 
                        vmin = 20, vmax = 100)
        """ cmap is the color map being used. 'jet' ranges from blue at 
            minimum to red at maximum. rstride and cstride define the 
            plotted grid as having increments of 1 in the x and y 
            directions. vmin and vmax set the limits for the color map
            - so blue is defined as approximately 20 C and below, while
            red is 100 C and above. 
            """                 
                    
        ax.view_init(elev = 90, azim = 0)
        ax.set_zlim3d(20, 100)
        """ view_init sets the viewing elevation as 60 degrees above the
            xy plane, and 45 degrees from orthogonal to the xy plane. 
            zlim3d sets the limits on the z axis to the expected 
            temperatures. 
            """
        ax.set_zlabel('Deg. C')
        
        fig_path = (experiment_path + '\\Figures')
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        plt.savefig(fig_path + '\\' + frame_title + '.png')
        
        plt.show()
    pass
""" So, this one is kinda deprecated, because a replacement has been developed
    within the plif_temperature script... """