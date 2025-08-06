import matplotlib.pyplot as plt
import numpy as np
import struct
import skimage as io
import cv2 as cv
import cv2
import logging
import sys
import glob
from os import path
import colorama 
import scipy
from scipy.interpolate import LinearNDInterpolator
from IPython.display import clear_output
from skimage import measure
import skimage as io
from scipy import ndimage
from skimage.morphology import (square, rectangle, diamond, disk)
from scipy.stats import binned_statistic
import pickle

def wrap_ang2pio2(X):
    import numpy as np
    pi =np.pi
    
    """
    wrap_ang2pio2 unwraps the input angle to an angle between -pi/2 to pi/2.
    input: 
    the angle(s) to be unwraped. eitheer in pi domain or 2pi
    output: 
    it returns the result in radians.
    As adapted from: LR 2007.
    """
    theta_wraped = abs((abs(X) -((pi/2))) % ((pi)) - (pi/2))
    return theta_wraped
    
def interpolate(X,Y,Z):
    '''
    creates an interpolator
    '''
    from scipy import interpolate
    Xl = list(X.flatten())
    Yl = list(Y.flatten())
    Zl = list(Z.flatten())
    interpolator=interpolate.LinearNDInterpolator(np.array([Xl,Yl]).T,Zl)
    return interpolator
    
def extract_contours_from_plot(contour_plot):
    # Extract the contour paths from matplotlib's contour plot
    contours = []
    for collection in contour_plot.collections:
        for path in collection.get_paths():
            contour = path.vertices  # Get contour points
            contours.append(contour)
    return contours
    
def smooth_contour(contour, window_size):
    # Smoothing function using Hann window
    hann_window = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(window_size) / (window_size - 1))
    hann_window /= hann_window.sum()  # Normalize the Hann window

    # Separate x and y coordinates from the contour
    x_vals = contour[:, 0]
    y_vals = contour[:, 1]

    # Pad the contour data to apply the moving window (to handle the edges)
    x_padded = np.pad(x_vals, (window_size // 2,), mode='wrap')
    y_padded = np.pad(y_vals, (window_size // 2,), mode='wrap')

    # Apply the Hann window using convolution
    x_smoothed = np.convolve(x_padded, hann_window, mode='valid')
    y_smoothed = np.convolve(y_padded, hann_window, mode='valid')

    # Combine the smoothed x and y values
    smoothed_contour = np.vstack((x_smoothed, y_smoothed)).T

    return smoothed_contour

def wrap_angpi(X):
    import numpy as np
    pi =np.pi
    
    """
    wrap_angpi unwraps the input angle to an angle between -pi to pi.
    input: 
    the angle(s) to be unwraped. eitheer in pi domain or 2pi
    output: 
    it returns the result in radians.
    As adapted from: LR 2007.
    """
    theta_wraped = (((X-pi) ) % ((2*pi)) - pi)
    return theta_wraped

def compute_dl(theta, dx,dy):

    """
    Takes in the angle (theta) from the breaking front velocities [u and v], 
    dx & dy are the grid resolution.
    computes the dl [which is the length of the breaking front] per pixel
    Returns the individual breaker crest length 
     NB: This is  wraped to 90deg
    """
    import numpy as np
    pi = np.pi
    dx = abs(np.asarray(dx))
    dy = abs(np.asarray(dy))
    
    my_dls = list()
    for i in range(len(dx)):
        
        theta_critical = np.arctan(dy[i]/dx[i])
        #theta = np.arctan(abs(v)/abs(u))
        theta_o = abs(theta[i] + (pi/2))     

        theta_o = wrap_ang2pio2(theta_o) #wrap to 90
        gamma = (pi/2) - theta_o 

        if theta_o <= theta_critical:
            hypothenus = dx[i]/(np.cos(theta_o))
            my_dls.append(hypothenus)

        else:
            hypothenus = dy[i]/(np.cos((gamma)))
            my_dls.append(hypothenus) 
    return my_dls

def remove_background_brigthness(img, kernel_size = 400,left = 750, top = 350 , height = 1650,width = 1600,  plot_results=False):
    ''' based on KM11'''

    
    
    #kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (kernel_size, kernel_size))
    #kernel = io.morphology.diamond(300);
    
    kernel = io.morphology.diamond(kernel_size);
    #kernel = kernel[100:500,100:500]
    kernel = kernel[int(kernel_size/4):int(2*kernel_size  - kernel_size/4), int(kernel_size/4) : int(2*kernel_size - kernel_size/4)]
    
    background = cv.morphologyEx(img, cv.MORPH_OPEN, kernel)
    background_cropped= background[top:top+height, left:left+width]
    mean_background = np.nanmean(background_cropped)
    
    pixel_division = img.astype(float) / background.astype(float)
    corrected_float0 =  pixel_division * mean_background
    corrected_float = np.nan_to_num(corrected_float0, nan=0)
    
    corrected_clipped = np.clip(corrected_float, 0, 255)
    corrected = np.round(corrected_clipped).astype(int)
    
    corrected[img==255] = 255
    
    if plot_results:
        fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
        img_plot = axes[0, 0].imshow(img, cmap='gray', vmin=0, vmax=255)
        axes[0, 0].set_title('img')
        cropped_img_plot = axes[0, 1].imshow(pixel_division, cmap='gray', vmin=0, vmax=255)
        axes[0, 1].set_title('pixel_division')
        background_plot = axes[1, 0].imshow(background, cmap='gray', vmin=0, vmax=255)
        axes[1, 0].set_title('background')
        #corrected_plot = axes[1, 1].imshow(corrected, cmap='gray', vmin=255, vmax=256)
        corrected_plot = axes[1, 1].imshow(corrected, cmap='gray')

        axes[1, 1].set_title('corrected')
    
        # Add colorbars to each subplot
        fig.colorbar(img_plot, ax=axes[0, 0])
        fig.colorbar(cropped_img_plot, ax=axes[0, 1])
        fig.colorbar(background_plot, ax=axes[1, 0])
        fig.colorbar(corrected_plot, ax=axes[1, 1])
    
        plt.tight_layout()
        plt.show()
   
    return  corrected
    
def filterout_lowspeeds(CXo, CYo, inpo, jnpo, my_dls,slowspeed=1):
    """
    Filters out the very slow detected breaking that may not be trusted from visible cameras. Then constrains the breaking direction 
    Following sutherland and melville 2013. the breaking direction is restricted to plus or minus pi/2 in the direction of the dominant wind/waves.
    """
    
    speeds = np.sqrt(CXo**2 +CYo**2)
    index_speeds = np.where(speeds > slowspeed) ###- filter very slow velocities --
    CXo, CYo, inpo, jnpo, my_dls = CXo[index_speeds], CYo[index_speeds], inpo[index_speeds], jnpo[index_speeds], np.array(my_dls)[index_speeds]
    
    return inpo, jnpo, CXo, CYo, my_dls
    
def save_loc_data(datapath, fname_savestats, C, n,aTOT,  **data_dict):
    import pickle
    with open(datapath, 'wb') as file:
        pickle.dump(data_dict, file)
        
    C_min = np.nanmin(C) if len(C) > 0 else np.nan
    C_max = np.nanmax(C) if len(C) > 0 else np.nan
    C_mean = np.nanmean(C) if len(C) > 0 else np.nan
    C_median = np.nanmedian(C) if len(C) > 0 else np.nan

    np.savetxt(fname_savestats, (aTOT, C_min, C_mean, C_median, C_max, n, n+1),
               delimiter=", ", fmt='%s')

    print("--------------------+--------------------------")
    print("LOC stats in m/s:   |")
    print("--------------------+--------------------------")
    print("    Cmin:           | {:.2f}".format(C_min))
    print("    Cmax:           | {:.2f}".format(C_max))
    print("    Cmean:          | {:.2f}".format(C_mean))
    print("    Cmedian:        | {:.2f}".format(C_median))
    print("--------------------+--------------------------")

def triangle_smoothing_1d(data, kernel_size):
    
    # Ensure that kernel size is odd for symmetry
    if kernel_size % 2 == 0:
        raise ValueError("Kernel size must be odd for symmetry")

    # Create a uniform triangular kernel
    triangle_kernel = np.arange(1, kernel_size + 1)
    triangle_kernel = triangle_kernel / np.sum(triangle_kernel)
    
    smoothed_data = np.convolve(data, triangle_kernel, mode='same')
    return smoothed_data
