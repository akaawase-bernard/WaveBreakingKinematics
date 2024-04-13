import sys
sys.path.append('../code/auxiliary/')
import detector
from os import path
import os
import cv2
from concurrent.futures import ProcessPoolExecutor
import warnings
warnings.filterwarnings("ignore")


def remove_background(input_path, output_path, kernel_size=400):
    frame = cv2.imread(input_path)
    if frame is not None:
        img = detector.remove_background_brigthness(frame, kernel_size=kernel_size, plot_results=False)

        # Save the result with the same name as the input file in the output folder
        cv2.imwrite(output_path, img)
        print(f"Processed {input_path} and saved result to {output_path}")

def process_images(input_folder, output_folder, num_processors=40):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Get a list of input image paths
    input_paths = [os.path.join(input_folder, file) for file in os.listdir(input_folder) if file.endswith('.png')]

    with ProcessPoolExecutor(max_workers=num_processors) as executor:
        # Create the output paths based on the input paths and output folder
        output_paths = [os.path.join(output_folder, os.path.basename(path)) for path in input_paths]

        # Submit the tasks to the executor
        executor.map(remove_background, input_paths, output_paths)


if __name__ == "__main__":
    
    data_pth = '/d2/lab/data/bernard/wass/AA_2015/2015-03-05_10-35-00_12Hz/output_stbv2/'
    write_pth = '/scratch/lab/data/bernard/wass/AA_2015/2015-03-05_10-35-00_12Hz/output_stbv2/' 
    
    input_folder = data_pth + "stabilized_undistorted_imgs"
    output_folder = write_pth + "background_removed_matlab_stbv2"
    num_processors = 1 
    
    process_images(input_folder, output_folder, num_processors)


#visualize
if 0:
    import glob
    import matplotlib.pyplot as plt    
    
    b4backrmd = sorted(glob.glob(input_folder + '/*.png')) 
    backrmd = sorted(glob.glob(output_folder + '/*.png')) 
    
    frameb4 = cv2.imread(b4backrmd[2],0)
    frame_after = cv2.imread(backrmd[2],0)

    plt.figure()
    plt.imshow(frameb4, cmap = 'gray')
    
    plt.figure()
    plt.imshow(frame_after , cmap = 'gray')
