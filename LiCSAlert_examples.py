#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 10:55:48 2020

@author: matthew
"""

import sys
import pickle
from pathlib import Path
import copy
import numpy as np

import licsalert
from licsalert.monitoring_functions import LiCSAlert_monitoring_mode




# thinkpad
LiCSAlert_pkg_dir = Path("/home/matthew/university_work/crucial_1000gb/03_automatic_detection_algorithm/06_LiCSAlert/00_LiCSAlert_GitHub")          # path to LiCSAlert.  Could also be this directory.  
debug_scripts = "/home/matthew/university_work/python_stuff/python_scripts"

# bright
#LiCSAlert_pkg_dir = Path("/home/matthew/university_work/crucial_1000gb/03_automatic_detection_algorithm/06_LiCSAlert/00_LiCSAlert_GitHub")          # path to LiCSAlert.  Could also be this directory.  
#debug_scripts = "/home/matthew/university_work/crucial_1000gb/python_stuff/python_scripts"




if debug_scripts not in sys.path:                                                                             # check if already on path
    sys.path.append(debug_scripts)
from small_plot_functions import matrix_show, quick_linegraph

#%%

# from licsalert.aux import col_to_ma



# path = "/home/matthew/university_work/crucial_1000gb/03_automatic_detection_algorithm/06_LiCSAlert/00_LiCSAlert_GitHub/campi_flegrei_example/20210912/epoch_images_data.pkl"
# with open(path, 'rb') as f:
#     test = pickle.load(f)
# matrix_show([col_to_ma(test['cumulative'], test['mask'])])        
# matrix_show([col_to_ma(test['incremental'], test['mask'])])        
# matrix_show([col_to_ma(test['reconstruction'], test['mask'])])        
# matrix_show([col_to_ma(test['residual'], test['mask'])])        
    

# path = "/home/matthew/university_work/crucial_1000gb/03_automatic_detection_algorithm/06_LiCSAlert/00_LiCSAlert_GitHub/campi_flegrei_example/aux_figures/aux_images_data.pkl"

# with open(path, 'rb') as f:
#     test = pickle.load(f)
# matrix_show([col_to_ma(test['icasar_sources'][0,:], test['mask'])])        
# matrix_show([col_to_ma(test['dem'], test['mask'])])        



#%% Settings 


licsalert_settings = {"baseline_end" : "20170101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
                      "figure_intermediate" : False,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "figure_type"         : 'both',                             # either 'window' or 'png' (to save as pngs), or 'both'
                      "downsample_run" : 0.5,                                     # data can be downsampled to speed things up
                      "downsample_plot"      : 0.5,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "residual_type"        : 'cumulative'}                      # controls the type of residual used in the lower plot.  Either cumulative or window   


icasar_settings = {"n_comp" : 5,                                                  # number of components to recover with ICA (ie the number of PCA sources to keep)
                   "bootstrapping_param" : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                   "tsne_param" : (30, 12),                                       # (perplexity, early_exaggeration)
                   "ica_param" : (1e-2, 150),                                     # (tolerance, max iterations)
                   "hdbscan_param" : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
                   "ifgs_format"        : 'cum'}                                  # can be 'all', 'inc' (incremental - short temporal baselines), or 'cum' (cumulative - relative to first acquisition)



licsalert_dir = Path("./")                                                        # outdir parent
volcano = 'campi_flegrei_example'                                                 # outdir final                                                                                                                      
licsbas_dir = Path("./022D_04826_121209")                                         # input data


LiCSAlert_monitoring_mode(region = None, volcano = volcano, LiCSAlert_pkg_dir = LiCSAlert_pkg_dir,
                          licsbas_dir = licsbas_dir, licsalert_dir = licsalert_dir,
                          licsalert_settings = licsalert_settings, icasar_settings = icasar_settings)

    

#%% We can also make the LiCSAlert figure for all times, but this is slow

# licsalert_settings['figure_intermediate'] = True                                                    # Turn on figures for all times 
# licsalert_settings['figure_type'] = 'png'                                                           # There will be too many figure windows if this is 'window' or 'both' ! 

# volcano = 'campi_flegrei_example_all_times'                                                         # new outdir

# LiCSAlert_monitoring_mode(region = None, volcano = volcano, LiCSAlert_pkg_dir = LiCSAlert_pkg_dir,
#                           licsbas_dir = licsbas_dir, licsalert_dir = licsalert_dir,
#                           licsalert_settings = licsalert_settings, icasar_settings = icasar_settings)
