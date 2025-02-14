#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:28:07 2023

@author: matthew

"""

print("Succesfully started the script.  ")

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




#%%
# def find_directories_with_more_than_n_subdirectories(root_dir, n=2):
#     # List to store directories with more than n subdirectories
#     result_directories = []

#     # Walk through the directory tree
#     for root, dirs, files in os.walk(root_dir):
#         # Count the number of subdirectories in the current directory
#         num_subdirectories = len(dirs)
        
#         # If the number of subdirectories exceeds n, add the directory to the result list
#         if num_subdirectories > n:
#             result_directories.append(root)

#     return result_directories

# # Example usage:
# directory_to_search = '/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/05_jasmin_clones/02_2024_01_19_cov_only/'

# # Call the function
# directories_with_more_than_2_subdirectories = find_directories_with_more_than_n_subdirectories(directory_to_search, n=2)

# directories_with_more_than_2_subdirectories = sorted(directories_with_more_than_2_subdirectories)

# # Write the result to a text file
# with open('directories_with_more_than_2_subdirectories.txt', 'w') as f:
#     for directory in directories_with_more_than_2_subdirectories:
#         f.write(directory + '\n')
    


#%%

from licsalert.jasmin_tools import open_comet_frame_files, volcano_name_to_comet_frames
from licsalert.jasmin_tools import write_jasmin_download_shell_script
from licsalert.jasmin_tools import update_volcs_with_data, get_lon_lat_of_volcs


from licsalert.status import extract_licsalert_status, get_all_volcano_dirs
from licsalert.status import find_volcano, get_volc_names_fron_dir_of_frames

from licsalert.temporal import create_day_list
from licsalert.volcano_portal_tools import get_portal_public_volcanoes


from licsalert.plotting import status_fig_one_volc, status_fig_all_volcs
        
#%% Things to set

#%%  ################# jasmin outputs 2023 (for fringe)
# #licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/01_test/")
# licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/02_all/")
# out_dir = Path("./monitoring_jasmin/one_volc/ale_bagu")
# d_start = "20200101"
# d_stop = "20230901"
# test_volcano = 'ale_bagu*'
# regions = True
################## end jasmin outputs 2023 (for fringe)

#%% ################## jasmin outputs 2024 (for Cities on Volcanoes)

# option 1: volcanoes on public portal  
# cov_volc_names = ['Antisana',
#                  'Cerro Overo',
#                  'Cordon del Azufre',
#                  'Cotopaxi',
#                  'Domuyo',
#                  'Fernandina',
#                  #'Guagua Pichinca',                    # volcano portal spelling
#                  'Guagua Pichincha',                     # licsbas spelling
#                  'Laguna del Maule',
#                  'Masaya',
#                  'Nevados de Chillan',
#                  'Planchón-Peteroa',                    # 
#                  'Sabancaya',
#                  'San Miguel',
#                  'San Salvador',
#                  'Sangay',
#                  'Santa Ana',
#                  'Sierra Negra',
#                  'Tungurahua',
#                  'Turrialba',
#                  'Villarrica']

# # option 2: Camila deforming volcanoes list
# cov_volc_names = [ # camilla deforming list
#                   "Cordon Caulle",
#                   "Wolf",
#                   "Sabancaya",
#                   "Laguna del Maule",
#                   "Cordon del Azufre",
#                   "Sangay",
#                   "Domuyo",
#                   "Alcedo",
#                   "Mevados de Chillan",
#                   "Cerro Azul",
#                   "Fernandina",
#                   "Sierra Negra",
#                   # below are Visible on portal, def observed, central america
#                   "Masaya",                 
#                   "San Miguel",
#                   "Santa Ana",
#                   "Turrialba"]



# # ?
# # licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/"
# #                      "06_LiCSAlert/05_jasmin_clones/02_2024_01_19_cov_only")

# # current clone of Jasmin data
# licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/"
#                       "06_LiCSAlert/05_jasmin_clones/2024_06_11_comet_talk")

# # print(f"USING THE TEMPORARY TEST DIRECTORY")
# # licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/"
# #                       "06_LiCSAlert/05_jasmin_clones/test")


# out_dir = Path("./status_outputs/03_comet")
# # d_start = "20200101"                # LiCSAlert baseline ends at end of 2019?  
# d_start = "20180101"                # LiCSAlert baseline ends at end of 2019?  
# d_stop = "20240621"
# test_volcano = 'sabancaya*'
# regions = True

# # volcs_omit = ["guagua_pichincha_121A_08908_131313",  # no coherence
# #               "planchon-peteroa_018A_12472_131313",  # single pixels of large value
# #               ]
# volcs_omit = []

################## end jasmin outputs 2024 (for Cities on Volcanoes)
#%% ################## jasmin outputs 2024 Chile 

# # current clone of Jasmin data
# licsalert_dir = Path(
#     "/home/matthew/university_work/03_automatic_detection_algorithm/"
#     "06_LiCSAlert/05_jasmin_clones/2024_11_14_chile_v2"
#     )

# # names of all comet frames.  
# comet_volcano_frame_index_dir = Path('./comet_volcano_frames/')

# out_dir = Path("./status_outputs/04_chile_november_2024")
# d_start = "20210101"                
# d_stop = "20241101"
# #test_volcano = 'sabancaya*'
# regions = False
# volcs_omit = []
# generate_bash_dl_file = True
# ################## end jasmin outputs 2024 for Chile

#%% ################## 10 years of Sentinel paper (2025_02_12)

# current clone of Jasmin data
# licsalert_dir = Path(
#     "/home/matthew/university_work/03_automatic_detection_algorithm/"
#     "06_LiCSAlert/05_jasmin_clones/2024_11_14_chile_v2"
#     )

# names of all comet frames.  
comet_volcano_frame_index_dir = Path('./comet_volcano_frames/')

#out_dir = Path("./status_outputs/04_chile_november_2024")
d_start = "20210101"                
d_stop = "20241101"
#test_volcano = 'sabancaya*'
regions = False
volcs_omit = []
generate_bash_dl_file = True
################## end jasmin outputs 2024 for Chile


#%%  #################  local test volcs (processed with LiCSBAS by me?  )
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

#%% Step 00: volcano names:
   
try:
    
    volc_names = get_volc_names_fron_dir_of_frames(
        licsalert_dir, regions = regions
        )
    
    print(
        "Succesfully created a list of volcanoes (volc_names) from a "
        "licsalert directory.  Continuing.  "
        )
except:
    print(
        "Failed to get the volcano names from a licsalert directory.  Perhaps "
        "this is the first time this has been run?  You could manually "
        "provide a list of volc_names, but for now we'll try to generate "
        "one from the volcanoes that are visible on the COMET volcano portal"
        )

    try:
        # get the names of the COMET portal public volcanoes
        volc_names = get_portal_public_volcanoes()
        print(
            "Succesfully generated a list of volcanoes that are visible "
            "publicly on the COMET volcano portal.  Continuing.  "
            )
    except:
        print(
            "Failed to generate a list of volcanoes from the COMET volcano "
            "portal.  Exiting.  ")


# open the info on COMET volcano frames and what region they're in.  
# region is a key, and each value isa list of comet frames in that region
comet_volcano_frame_index = open_comet_frame_files(comet_volcano_frame_index_dir)

# convert volc_names to volcs, which is a list of comet_volcano objects
# output is jasmin_sync_script.sh
volcs = volcano_name_to_comet_frames(volc_names, comet_volcano_frame_index)    

# tidy up.  
del volc_names



#%% Step 01: Generate bash file to copy selected data from Jasmin



if generate_bash_dl_file:
    jasmin_dir = Path(
        "mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/"
        "projects/LiCS/volc-portal/processing_output/licsalert"
        )
    local_dir = Path('./licsalert_sync/')
    bash_sync_fle = "jasmin_sync_script.sh"
    

    # write a shell script to download the LiCSAlert data from jasmin 
    write_jasmin_download_shell_script(
        jasmin_dir, local_dir, bash_sync_fle, volcs, exclude_json_gz = True,
        exclude_original_ts = False, exclude_ICASAR_data = True,
        exclude_epoch_images = True, exclude_epoch_data = True)



pdb.set_trace()



#%% Step 02:
    
# run jasmin sync file on FOE-linux (to copy data from jasmin to school server)
# copy from FOE-linux to local.  

#%% Step 02a: possibly move directories that don't have LiCSAlert results
# to a different directory

from licsalert.jasmin_tools import move_small_directories

move_small_directories(
    licsalert_dir, licsalert_dir.parent / "to_delete", n = 6)

pdb.set_trace()

#%% Step 03: Compile licalert status for all frames at all times.  

# get path to licsalert outputs for each volcanoes, and their names (snake case)
volc_frames, volc_frame_names = get_all_volcano_dirs(licsalert_dir, volcs_omit,
                                                         regions)                                               

# update volcs so it only contains volcanoes we have data for.  
volcs = update_volcs_with_data(volcs, volc_frames, volc_frame_names,
                               verbose = False)

# add the lon lat info for each volcano
volcs = get_lon_lat_of_volcs(volcs)

# simple list of datetime for each status day
day_list = create_day_list(d_start, d_stop)

# get the licsalert status for each volcano frame (update volcs)
extract_licsalert_status(volcs, day_list)           





#%% Step 04: figure for one volcano at all times.  

# name_and_indexes = find_volcano(volc_frame_names, test_volcano)



# out_dir_1_volc = out_dir / name_and_indexes[0][0]

# status_fig_one_volc(name_and_indexes[-1], licsalert_status, day_list, 
#                 fig_type = 'png', out_dir = out_dir_1_volc)

# status_fig_one_volc(name_and_indexes[-1], licsalert_status, day_list, 
#                     fig_type = 'gif', volc_dirs = volc_frame_dirs, 
#                     out_dir = out_dir_1_volc/ "gif") 
# #pngs_to_gif(out_dir / "gif", out_dir  / "gif" / "animation.gif", image_duration = 0.5)






#%% Step 05: figure for all volcs at all times.  

# status_fig_all_volcs(licsalert_status, volc_frame_names, day_list,  
#                      out_dir = out_dir / "all_volcs", 
#                      figsize = (20, 12), plot_frequency = 6, label_fs = 16)

# # convert the directory of pngs to a single gif
# pngs_to_gif(out_dir / "all_volcs", out_dir / "all_volcs" / "animation.gif",
#             image_duration = 250)                                                 

    
#%% Step 06: New worldmap figure

from licsalert.plotting import offset_volc_lls

# adjust the lon lats of the volcanoes so when plotted they don't overlap
offset_volc_lls(volcs,  threshold = 2., offset = 2.1, attempts = 200)
    
    
from licsalert.plotting import licsalert_status_map
    
licsalert_status_map(volcs, day_list = day_list[-2:-1], 
                     outdir = Path("./status_outputs/03_comet/"),
                     figure_type = 'window', plot_frequency = 'daily')

pngs_to_gif(out_dir / "short_status_map", out_dir / "status_map" / "status_map_animation.gif",
            image_duration = 750)                                                 
    
    
    
    
    
    
    
    
    
    





