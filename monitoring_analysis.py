#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:28:07 2023

@author: matthew

# initialise something to store:


  dim 1: day
  dim 2; volcano
  dim 3: which value

  plus: 
    list of days
    list of volcano names.  

  # need to get number of days 
  # need to get number of volcanoes.  


# get the region directories.  
# loop through these.  
  # get the volcano directories.  
  # loop through these.  


2024 update:
    
    jasmin new LiCSAlert version.  
    list of volcanoes that are public.  
    copy all from Jasmin
    


"""

import numpy as np
import matplotlib.pyplot as plt 
from pathlib import Path
import pdb
import os
from glob import glob
import sys






#%% from other modules


def matrix_show(matrix, title=None, ax=None, fig=None, save_path = None, vmin0 = False):
    """Visualise a matrix
    Inputs:
        matrix | r2 array or masked array
        title | string
        ax | matplotlib axes
        save_path | string or None | if a string, save as a .png in this location.  
        vmin0 | boolean | 
        

    2017/10/18 | update so can be passed an axes and plotted in an existing figure
    2017/11/13 | fix bug in how colorbars are plotted.
    2017/12/01 | fix bug if fig is not None
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if ax is None:
        fig, ax = plt.subplots()
    matrix = np.atleast_2d(matrix)                   # make at least 2d so can plot column/row vectors

    if isinstance(matrix[0,0], np.bool_):           # boolean arrays will plot, but mess up the colourbar
        matrix = matrix.astype(int)                 # so convert

    if vmin0:
        matrixPlt = ax.imshow(matrix,interpolation='none', aspect='auto', vmin = 0)
    else:
        matrixPlt = ax.imshow(matrix,interpolation='none', aspect='auto')
    fig.colorbar(matrixPlt,ax=ax)
    if title is not None:
        ax.set_title(title)
        fig.canvas.manager.set_window_title(f"{title}")

    if save_path is not None:                                                   # possibly save the figure
        if title is None:                                                       # if not title is supplied, save with a default name
            fig.savefig(f"{save_path}/matrix_show_output.png")
        else:
            fig.savefig(f"{save_path}/{title}.png")                             # or with the title, if it's supplied 
            
    
    plt.pause(1)                                                                    # to force it to be shown when usig ipdb



def pngs_to_gif(input_folder, output_file, image_duration = 1000):
    """A function to conbime a folder of .pngs into one gif.
    Inputs:
        input_folder | pathlib Path | the location of the pngs.  E.g. "output_data".  Must be a pathlib Path
        output_file  | string | e.g. "output.gif"
        image_duration | float | the time each png is displayed for in milliseconds.  e.g. 0.5
    Returns:
        gif
    History:
        2018/07/23 | MEG: adapted from a script.
        2018/12/05 | paths must now be absolute.
        2020/08/18 | MEG | Update to handle paths formatted using the pathlib module.  
        2021_05_04 | MEG | Update to handle paths better.  
        2023_09_04 | MEG | Add funtion to make sure all images are the same size.  
        2023_09_09 | MEG | GIF time units have changed in milliseconds in new version of imageio - update here.  
    """

    import imageio
    import os
    from PIL import Image
    import re

    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [ atoi(c) for c in re.split('(\d+)', text) ]
    
    def determine_minimim_image_size(images):
        """ Iterate along the list of images and find the size of the smallest image
        """
        for image_n, image in enumerate(images):
            ny, nx, _ = image.shape
            if image_n == 0:
                min_y = ny
                min_x = nx
            else:
                if ny < min_y:
                    min_y = ny
                if nx < min_x:
                    min_x = nx
        return min_y, min_x

           
    def crop_if_larger(images, min_y, min_x):
        """ If an image is larger than the minimum, crop it.  
        """
        for image_n, image in enumerate(images):
            ny, nx, _ = image.shape
            if (ny > min_y) or (nx > min_x):
                print(f"Resizing image from {ny} to {min_y} in y, and {nx} to {min_x} in x ")
                images[image_n] = image[:min_y, :min_x, :]
        return images

    # sort the files into a human order
    image_names = os.listdir(input_folder)              # get 
    image_names.sort(key=natural_keys)                  # sort into the correct (human) order

    # write the images to the .gif
    images = []
    for counter, filename in enumerate(image_names):
        images.append(imageio.imread(input_folder / filename))                                  # input_folder must be a pathlib Path to allow simple extension with /
        print(f"Processed {counter} of {len(image_names)} images.  ")
    
    min_y, min_x = determine_minimim_image_size(images)                             # fnid the size of the minimum image
    images = crop_if_larger(images, min_y, min_x)                                   # make sure there are none bigger than that.  

    # for image in images:
    #     print(image.shape)
    
    print('Writing the .gif file...', end = "")
    imageio.mimsave(output_file, images, format = 'GIF', duration = image_duration,  loop = 1)
    print('Done!')


#%%


from licsalert.status import extract_licsalert_status, get_all_volcano_dirs
from licsalert.status import find_volcano

from licsalert.temporal import create_day_list

from licsalert.plotting import status_fig_one_volc, status_fig_all_volcs
        
#%% Things to set


################## jasmin outputs 2023 (for fringe)
# #licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/01_test/")
# licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/02_all/")
# out_dir = Path("./monitoring_jasmin/one_volc/ale_bagu")
# d_start = "20200101"
# d_stop = "20230901"
# test_volcano = 'ale_bagu*'
# regions = True
################## end jasmin outputs 2023 (for fringe)


################## jasmin outputs 2024 (for Cities on Volcanoes)

licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/05_jasmin_clones/02_2024_01_19_cov_only")
out_dir = Path("./status_outputs/cov")
d_start = "20150101"
d_stop = "20230901"
test_volcano = 'sabancaya*'
regions = True
################## end jasmin outputs 2024 (for Cities on Volcanoes)



#################  local test volcs (processed with LiCSBAS by me?  )
# all_volcs_figs_dir = Path('monitoring_steps_all_volcs_test')
# licsalert_dir = Path("/home/matthew/university_work/31_from_bright/2023_09_08/01a_LiCSAlert_batch_mode_2/")
# regions = False

# # Option for the 1 volc we plot in detail.  
# # option 1
# out_dir = Path("./monitoring_test/one_volc/sierra_negra")
# d_start = "20170101"
# d_stop = "20230901"
# test_volcano = '*sierra*'

# # option 2
# # out_dir = Path("./monitoring_test/one_volc/erta_ale")
# # d_start = "20160601"
# # d_stop = "20230901"
# # test_volcano = '*erta*'

################# end local test volcs            


#%% ?

# get path to outputs for each volcanoes, and their names (snake case)
volc_frame_dirs, volc_frame_names = get_all_volcano_dirs(licsalert_dir, regions)                                               

# simple list of datetime for each status day
day_list = create_day_list(d_start, d_stop)

# get the licsalert status for each volcano frame, 
# also remove dates where no statuses change.  
output = extract_licsalert_status(volc_frame_dirs, volc_frame_names, day_list)           
(licsalert_status, volc_frame_names, volc_frame_dirs) = output; del output

#%% figure for one volcano at all times.  

name_and_indexes = find_volcano(volc_frame_names, test_volcano)



out_dir_1_volc = out_dir / name_and_indexes[0][0]

status_fig_one_volc(name_and_indexes[-1], licsalert_status, day_list, 
                fig_type = 'png', out_dir = out_dir_1_volc)

status_fig_one_volc(name_and_indexes[-1], licsalert_status, day_list, 
                    fig_type = 'gif', volc_dirs = volc_frame_dirs, 
                    out_dir = out_dir_1_volc/ "gif") 
#pngs_to_gif(out_dir / "gif", out_dir  / "gif" / "animation.gif", image_duration = 0.5)






#%% figure for all volcs

status_fig_all_volcs(licsalert_status, volc_frame_names, day_list,  
                     out_dir = out_dir / "all_volcs", 
                     figsize = (20, 12), plot_frequency = 6, label_fs = 16)

pngs_to_gif(out_dir / "all_volcs", out_dir / "all_volcs" / "animation.gif", image_duration = 250)                                                 # convert all times into a gif


#%%

    
    
    
    
    
    
    
    
    
    
    
    





