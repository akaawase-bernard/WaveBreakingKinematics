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
from scipy.interpolate import interp1d


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
    
def remove_nesting(my_nested_list):
    """
    Flattens list of list.
    """
    try:
        from collections.abc import Iterable
    except ImportError:
        from collections import Iterable
    for item in my_nested_list:
        if isinstance(item, Iterable) and not isinstance(item, str):
            
            for sub_element in remove_nesting(item):
                 yield sub_element
        else:        
             yield item

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


def zero2nans_2d(data):
    '''
    sets all zeros as nans
    '''
    data = data.astype('float')
    data [data  == 0] = 'nan' 
    return data



def get_KfromC_intermediate(c, depth, g = 9.81, check=0):
    
    '''
    This function computes wavenumber(K) from wave phase speed (c) and water depth (h)
    based on intermediate water dispersion relationship.
    
    Adapted from the explanation of LR 2023
    coded by Bernard Akaawase
    '''
    
    pi = np.pi

    #initial guess based on deepwater 
    Kdeep = g/(c**2)
    
    for item in range(100):
        if item == 0:
            Kinitial = Kdeep
        else:
            Kinitial = Kintermediate
        Kfunction = c**2 - g*np.tanh(Kinitial*depth)/Kinitial
        Kderivative = g * (np.tanh(Kinitial * depth) - Kinitial * depth * 1/np.cosh(Kinitial * depth)**2) / Kinitial**2
        Kintermediate = Kinitial - (Kfunction/Kderivative)

    if check:
        print('The calculated value of C:', np.sqrt(g/Kintermediate),'m/s')
        print('The inputed value of C:', c,'m/s')

    return Kintermediate


def generate_img_pairs(base_dir, n):

    '''
    This function fetches the stabilized then processed wass data corresponding to an instance (n). 
    '''
    file_pairs = []

    # Generate the file pairs directly
    for i in range(n, n + 2):
        file_string = path.join(base_dir, f"{i:06d}_{n:06d}.png")
        file_pairs.append(file_string)

    return sorted(file_pairs)
    
def generate_file_pairs(base_dir, n):

    '''
    This function fetches the stabilized then processed wass data corresponding to an instance (n). 
    '''
    file_pairs = []

    # Generate the file pairs directly
    for i in range(n, n + 2):
        file_string = path.join(base_dir, f"{i:06d}_{n:06d}_wd/")
        file_pairs.append(file_string)

    return sorted(file_pairs)
    
def mean_dir_spread(k, sigma1_brk_r19):

    # Filter k and spreading_deg
    kvec = np.array(k)
    mask = (k > 0.1) & (k < 4.1)

    filtered_kvec = kvec[mask]
    filtered_spreading_deg = sigma1_brk_r19[mask]

    # Define a constant kr grid
    kr = np.linspace(filtered_kvec.min(), filtered_kvec.max(), len(kvec))

    # Interpolate spreading values onto the constant dk grid
    fc = interp1d(filtered_kvec, filtered_spreading_deg, kind='linear')
    spreading_interp = fc(kr)
    #spreading_interp = fc(filtered_kvec )

    mean_spreading = np.nansum(spreading_interp * np.gradient(kr)) / np.nansum(np.gradient(kr))
    #mean_spreading = np.nansum(spreading_interp * np.gradient(filtered_kvec)) / np.nansum(np.gradient(filtered_kvec))
    print('The mean dir. spread is:', mean_spreading) 


def draw_contours(frame, threshold ): 
    """
    input = first image, second image
    output = contours of the detected breakers
    """
    
    import cv2

    _, binary_mask = cv2.threshold(frame, threshold * 255, 255, cv2.THRESH_BINARY)
    contours, _ = cv2.findContours(binary_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

    return contours   
    
def buffer_image_edges(image_a, image_b, ip1, jp1, ip2, jp2):
    '''
    Create a pad of 255 around images. Receives single-channel images (e.g grayscale image)
    Outputs the padded images.
    '''
    Ih, Iw = np.shape(image_a)


    # Blank pixels outside the crop window
    image_a[:, :ip1] = 0
    image_a[:jp1, :] = 0
    image_a[:, ip2:] = 0
    image_a[jp2:, :] = 0

    image_b[:, :ip1] = 0
    image_b[:jp1, :] = 0
    image_b[:, ip2:] = 0
    image_b[jp2:, :] = 0

    return image_a, image_b

def mask_lower_corner(image_a, image_b, y2, x2):

    '''
    y2 and x2 are the pixels values of the upper corner that is to be masked.
    '''
    
    #mask region
    corner_width =  x2
    corner_height = np.shape(image_a)[0] - y2 # just suppy the ypixel value here
    
    image_a[-corner_height:, :corner_width] = 0 #255 #0
    image_b[-corner_height:, :corner_width] = 0 #255 #0

    return image_a, image_b
    
def my_guassian_blur(frame, sigma_x =2, sigma_y=2):
    import cv2
    
    return cv2.GaussianBlur(frame, (7, 7), sigma_x, sigma_y)


def load_camera_mesh( meshfile ):
    
    with open(meshfile, "rb") as mf:
        npts = struct.unpack( "I", mf.read( 4 ) )[0]
        limits = np.array( struct.unpack( "dddddd", mf.read( 6*8 ) ) )
        Rinv = np.reshape( np.array(struct.unpack("ddddddddd", mf.read(9*8) )), (3,3) )
        Tinv = np.reshape( np.array(struct.unpack("ddd", mf.read(3*8) )), (3,1) ) 
                
        data = np.reshape( np.array( bytearray(mf.read( npts*3*2 )), dtype=np.uint8 ).view(dtype=np.uint16), (3,npts), order="F" )
        mesh_cam = data.astype( np.float32 )
        mesh_cam = mesh_cam / np.expand_dims( limits[0:3], axis=1) + np.expand_dims( limits[3:6], axis=1 );
        mesh_cam = Rinv@mesh_cam + Tinv;
    
        return mesh_cam


def load_pivlab_data(file_path):
    """
    Load PIVLab data from a .txt file and extract 'x' (as integers), 'y' (as integers), 'u' (as floats), and 'v' (as floats).

    Parameters:
    file_path (str): The path to the PIVLab data file in .txt format.

    Returns:
    x (list of int): List of x-coordinate values as integers.
    y (list of int): List of y-coordinate values as integers.
    u (list of float): List of u-velocity values as floats (px/frame).
    v (list of float): List of v-velocity values as floats (px/frame).
    """
    x = []
    y = []
    u = []
    v = []

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if not line.strip():
                continue  # Skip empty lines
            parts = line.strip().split('\t')
            if len(parts) == 5 and parts[0].replace(".", "").isdigit():
                x.append(int(parts[0]))
                y.append(int(parts[1]))
                u.append(float(parts[2]))
                v.append(float(parts[3]))
    print('Loading PIVLab data:', file_path)
    return x, y, u, v
    
def normalize_1d(arr):
    max_value = np.nanmax(arr)

    # Check if the maximum value is not zero to avoid division by zero
    if max_value != 0:
        normalized_arr = arr / max_value
    else:
        # Handle the case where all values are zero
        normalized_arr = arr.copy()  # Create a copy of the input array
        
    return normalized_arr

def compute_sea_plane_RT( plane  ):
    assert len(plane)==4, "Plane must be a 4-element vector"
    a=plane[0]
    b=plane[1]
    c=plane[2]
    d=plane[3];
    q = (1-c)/(a*a + b*b)
    R=np.array([[1-a*a*q, -a*b*q, -a], [-a*b*q, 1-b*b*q, -b], [a, b, c] ] )
    T=np.expand_dims( np.array([0,0,d]), axis=1)
    
    return R, T
    
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


def directional_filter(CXo, CYo, inpo, jnpo, wind_wave_dir_deg=22.5):
    """ Filter data based on directional constraints relative to wind/wave direction.

    Args:
        CXo (np.array): X-coordinates or corresponding component array.
        CYo (np.array): Y-coordinates or corresponding component array.
        inpo (np.array): Array of indices or any other related data.
        jnpo (np.array): Array of indices or any other related data.
        wind_wave_dir_deg (float): Wind/wave direction in degrees.

    Returns:
        tuple: Filtered CXo, CYo, inpo, jnpo arrays.
    """
    lower_limit = -np.pi / 2
    upper_limit = np.pi / 2
    wind_wave_dir = np.deg2rad(wind_wave_dir_deg)

    # Calculate the median direction
    dir = wrap_angpi(np.arctan2(CYo, CXo) - np.pi / 2)
    med_dir = np.nanmedian(dir)

    # Calculate the absolute angular difference and adjust to [-pi, pi]
    angular_diff = wrap_angpi(med_dir - wind_wave_dir)

    # Check if the median direction is within ±120 degrees of the wind/wave direction
    if np.abs(angular_diff) <= np.deg2rad(120):
        print("med_dir is within ±120 degrees of wind_wave_dir.")
        angular_difference = np.abs(wrap_angpi(dir - med_dir))
        indices_within_range = np.where(angular_difference <= np.pi / 2)
    else:
        print("med_dir is NOT within ±120 degrees of wind_wave_dir.")
        lower_limit = -np.deg2rad(120)
        upper_limit = abs(lower_limit)
        theta_ind = np.where((dir >= lower_limit) & (dir <= upper_limit))
        CXo = CXo[theta_ind]
        CYo = CYo[theta_ind]
        jnpo = jnpo[theta_ind]
        inpo = inpo[theta_ind]
        return CXo, CYo, inpo, jnpo

    # Return filtered data for the ±90 degrees case
    CXo = CXo[indices_within_range]
    CYo = CYo[indices_within_range]
    jnpo = jnpo[indices_within_range]
    inpo = inpo[indices_within_range]
    return CXo, CYo, inpo, jnpo
    
def align_on_sea_plane_RT( mesh, R, T ):
    # Rotate, translate
    mesh_aligned = R@mesh + T;
    # Invert z axis
    mesh_aligned[2,:]*=-1.0;
    return mesh_aligned


def align_on_sea_plane( mesh, plane ):
    assert mesh.shape[0]==3, "Mesh must be a 3xN numpy array"
    R,T = compute_sea_plane_RT( plane )
    return align_on_sea_plane_RT(mesh, R,T)

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
    
def filter_contours_by_size(contours, min_size):
    filtered_contours = [contour for contour in contours if cv2.contourArea(contour) > min_size]
    return filtered_contours


def fetch_contour_coords(data):
    data_squeezed = np.squeeze(data)
    
    if data_squeezed.size < 2:
        # Handle the case where the contour has insufficient points
        return None, None
    
    x_coordinates = data_squeezed[:, 0]
    y_coordinates = data_squeezed[:, 1]
    
    return x_coordinates, y_coordinates



def get_outward_vectors(ipc, jpc, CXi, CYi):
    """
    Calculate the outward vectors for the contour points based on their coordinates and velocity components.

    Args:
    ipc (np.array): Array of x-coordinates.
    jpc (np.array): Array of y-coordinates.
    CXi (np.array): Array of x-components of velocity.
    CYi (np.array): Array of y-components of velocity.

    Returns:
    tuple: Arrays of outward contour point coordinates, their normal vectors, and velocity components.
    """
    dXipc = []
    dYjpc = []

    # Check if jpc has elements and compute the differences with wrapping
    if len(jpc) > 0:
        dYjpc = jpc[1:] - jpc[:-1]
        tmpy = jpc[0] - jpc[-1]
        dYjpc = np.append(dYjpc, tmpy)  # appending looped values during differentiation

    # Check if ipc has elements and compute the differences with wrapping
    if len(ipc) > 0:
        dXipc = ipc[1:] - ipc[:-1]
        tmpx = ipc[0] - ipc[-1]
        dXipc = np.append(dXipc, tmpx)
    print('Finding the outward vectors')

    # Normal outward vector to the contours
    nCX = -1 * dYjpc
    nCY = dXipc

    # Project our velocities (as evaluated on the contours) onto the normal vectors
    nxp = nCX * CXi + nCY * CYi  
    outward_vec = np.argwhere(nxp > 0)  # find the projection that is greater than zero == outward vector

    # Extracting information for outward vectors
    inpo, jnpo = ipc[outward_vec], jpc[outward_vec]
    nCXo, nCYo = nCX[outward_vec], nCY[outward_vec]
    CXo, CYo = CXi[outward_vec], CYi[outward_vec]

    return inpo, jnpo, CXo, CYo
    
def load_camXYZp_map2IJ_without_putin2grid(wass_frame, fidstr, all_planes, baseline =2.5, checks=0):
    '''
    input: the working directory, as generated by during 3d reconstruction of sea surface eg with 
            WASS. This typically lives in the out directory. The should be of the 
            form 0000023_wd
            This functions loads xyz, removes plane, and map to IJ without interpolation. Uses individual plane files
    output: brings out the x,y,z and the projection in pixel. if checks is set as 1, plots 
            will be made.    
    '''
    import sys
    import glob
    import matplotlib.pyplot as plt
    import cv2 as cv
    from os import path
    import numpy as np
    import colorama #pip install colorama
    colorama.init()
    
    meshname = path.join(wass_frame,"mesh_cam.xyzC")
    print("Loading ", meshname )

    R = np.loadtxt(path.join(wass_frame,'Cam0_poseR.txt'))
    T = np.loadtxt(path.join(wass_frame,'Cam0_poseT.txt'))
    P0Cam =  np.vstack( (np.loadtxt( path.join(wass_frame,"P0cam.txt"))  ,[0, 0, 0, 1] ) )
    P1Cam =  np.vstack( (np.loadtxt( path.join(wass_frame,"P1cam.txt"))  ,[0, 0, 0, 1] ) )

    I = cv.imread(  path.join(wass_frame,"undistorted","00000001.png"), cv.IMREAD_ANYCOLOR )
    Iw,Ih = I.shape[1],I.shape[0]
    if Iw is None or Ih is None:
        #print("Unable to determine the camera image size. Please set it manually with -Iw,-Ih program arguments")
        sys.exit(-3)

    Rpl, Tpl = compute_sea_plane_RT( all_planes )
    mesh = load_camera_mesh(meshname)
    mesh_aligned = align_on_sea_plane( mesh, all_planes) * baseline
    
    Ri = Rpl.T
    Ti = -Rpl.T@Tpl
    RTplane = np.vstack( (np.hstack( (Ri,Ti) ),[0,0,0,1]) )
    toNorm = np.array( [[ 2.0/Iw, 0     , -1, 0],
                        [ 0     , 2.0/Ih, -1, 0],
                        [ 0,      0,       1, 0],
                        [ 0,      0,       0, 1]], dtype=np.float64 )

    SCALEi = 1.0/baseline
    P0plane = toNorm @ P0Cam @ RTplane @ np.diag((SCALEi,SCALEi,-SCALEi, 1))
    Npts = len(mesh[2,:]);
    my_ones = np.ones(len(mesh[2,:]))
    #my_ones = np.ones(len(mesh_aligned[2,:]))
    mesh_reshaped = np.vstack([mesh,my_ones])
    #mesh_reshaped = np.vstack([mesh_aligned,my_ones])
    
    P1 = np.resize(P1Cam,(3,4))
    pt2d = np.matmul(P1,mesh_reshaped);
    pt2d = pt2d / np.tile( pt2d[2,:],(3,1));
    
    plane_reshaped = np.tile(all_planes,(1,Npts) )
    my_plane = np.resize(plane_reshaped, (4,Npts))
    zp = mesh_aligned[2,:]
    
    if checks:
        plt.figure(figsize=(11,8))
        plt.imshow(I, cmap = 'gray', origin = 'upper')
        c= plt.scatter(pt2d[0,:],pt2d[1,:],
                       c=mesh_aligned[2,:], s =0.02, 
                       cmap= 'jet',
                        vmin=-1.5, vmax = 1.5, 
                       marker = '.')
        
        plt.title('Frame ' + fidstr)
        c=plt.colorbar(c)
        c.set_label('elevation (m)')
        plt.show()
        
    i = pt2d[0,:]
    j = pt2d[1,:]
    x = mesh_aligned[0,:]
    y = mesh_aligned[1,:]
    
    return i, j, x, y, zp


def remove_background_brigthness(img, kernel_size = 400,left = 750, top = 350 , height = 1650,width = 1600,  plot_results=False):
    ''' based on KM11'''

    #kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (kernel_size, kernel_size))

    kernel = io.morphology.diamond(300);
    kernel = kernel[100:500,100:500]
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
    
def filterout_lowspeeds(CXo, CYo, inpo, jnpo, my_dls,slowspeed=1):
    """
    Filters out the very slow detected breaking that may not be trusted from visible cameras. Then constrains the breaking direction 
    Following sutherland and melville 2013. the breaking direction is restricted to plus or minus pi/2 in the direction of the dominant wind/waves.
    """
    #print('Running a simple filter')#simple filter based on winddir
    
    speeds = np.sqrt(CXo**2 +CYo**2)
    index_speeds = np.where(speeds > slowspeed) ###- filter very slow velocities --
    CXo, CYo, inpo, jnpo, my_dls = CXo[index_speeds], CYo[index_speeds], inpo[index_speeds], jnpo[index_speeds], np.array(my_dls)[index_speeds]


    return inpo, jnpo, CXo, CYo, my_dls
    
def ffgrid(x,y,z, xvec, yvec, dx,dy):
    
    """
    The input:
    The function recieves ungridded x, y, z as 1-D arrays, 
    generates a grid of resolution specified by the user as 
    xvec & yvec. Note, xvec & yvec should be center bins.

    The ouput: 
    The gridded Z on the X,Y space 
    """
        
    import numpy as np    
    from scipy.stats import binned_statistic_2d

    X,Y = np.meshgrid(xvec, yvec)

    #establish edge bins
    x_edge = xvec - (dx/2.)
    x_edge = np.append(x_edge, max(x_edge)+ dx)

    y_edge = yvec - (dy/2.)
    y_edge = np.append(y_edge, max(y_edge)+ dy)

    #call binning function
    ret = binned_statistic_2d(x, y, z, 'mean', bins=[x_edge, y_edge], 
        expand_binnumbers=False)

    Z = ret.statistic
    Z = Z.T #need to transpose to match the meshgrid orientation with is (J,I)
    return X, Y, Z 

def clear_borders(frame_b, inpo, jnpo, CXo, CYo, my_dls, radius=4, plot=0):
    """
    The radius sets how far from the masked region you want your data to start
    """
    # Create a binary mask based on zero intensity points in the image
    mask = (frame_b == 0)

    # Initialize lists to store filtered vector data
    filtered_x = []
    filtered_y = []
    filtered_v = []
    filtered_u = []
    filtered_dls = []

    for xi, yi, ui, vi, my_dl in zip(inpo, jnpo, CXo, CYo, my_dls):
        xi_int = int(xi)  # Convert x coordinate to integer
        yi_int = int(yi)  # Convert y coordinate to integer

        if not np.any(mask[yi_int - radius:yi_int + radius + 1, xi_int - radius:xi_int + radius + 1]):
            filtered_x.append(xi)
            filtered_y.append(yi)
            filtered_v.append(vi)
            filtered_u.append(ui)
            filtered_dls.append(my_dl)

    # Convert the filtered vectors to NumPy arrays
    filtered_x = np.array(filtered_x)
    filtered_y = np.array(filtered_y)
    filtered_v = np.array(filtered_v)
    filtered_u = np.array(filtered_u)
    filtered_dls = np.array(filtered_dls)

    if plot:
        plt.figure()
        # Plot the vectors
        plt.imshow(frame_b, cmap='gray')
        plt.quiver(filtered_x, filtered_y, filtered_u, filtered_v, color='r', angles='xy', scale_units='xy', scale=1)
        plt.show()

    return filtered_x, filtered_y, filtered_u, filtered_v, filtered_dls

def save_loc_data(datapath, fname_savestats, C, n,aTOT,  **data_dict):
    import pickle
    # Save the data using pickle
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
    