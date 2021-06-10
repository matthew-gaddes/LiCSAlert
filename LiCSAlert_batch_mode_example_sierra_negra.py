#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 10:55:48 2020


Things to address:
    
    downsampling is compound?


@author: matthew
"""

import sys
import pickle
from pathlib import Path
import copy
import numpy as np

sys.path.append("./lib")
from LiCSAlert_functions import LiCSAlert_batch_mode

#%% Load Sentinel-1 data for Sierra Negra, note that this wasn't processed with LiCSBAS, and doesn't include the DEM.  

displacement_r2 = {}                                                                    # initiate a dictionary to store some info in.  
with open(f"sierra_negra_example_data.pkl", 'rb') as f:
    phUnw_files = pickle.load(f)                                                        # these are the names of the interferograms used (e.g. YYYYMMDD_YYYYMMDD.unw)
    displacement_r2["incremental"] = pickle.load(f)                                     # the incremental interferograms as row vectors.
    displacement_r2["mask"] = pickle.load(f)                                            # a mask to convert a row vector back to a rank 2 array (col_to_ma in LiCSAlert_functions.py does this easily.  )
    cumulative_baselines = pickle.load(f)                                               # the cumulative baselines. ie if there are acquisitions every 12 days, these would be 12,24,36 etc.  
    acq_dates = pickle.load(f)                                                          # acquisition dates.  This should be one longer than the names of the interferograms
    lons = pickle.load(f)                                                               # matrix of longitudes for bottom row (ie rank 1)
    lats = pickle.load(f)                                                               # matrix of latitdues for left column (ie rank 1)
f.close()

displacement_r2['lons'] = np.repeat(lons[np.newaxis, :],lats.shape[0], 0 )              # convert from rank 1 to rank 2
displacement_r2['lats'] = np.repeat(lats[::-1, np.newaxis],lons.shape[0], 1 )           # convert from rank 1 to rank 2
del lons, lats


displacement_r2['ifg_dates'] = [phUnw_file.split('.')[0] for phUnw_file in phUnw_files] # use the file names to get the dates that each interferogram spans.  
displacement_r2_copy = copy.deepcopy(displacement_r2)                                   # make a copy for use with the 2nd example  



#%% Example 1, creating one LiCSAlert figure for the complete time series
ICASAR_path = Path("/home/matthew/university_work/15_my_software_releases/ICASAR-2.2.0/lib/")                               # location of ICASAR functions
#ICASAR_path = Path("/home/matthew/university_work/01_blind_signal_separation_python/13_ICASAR/ICASAR_GitHub/lib")           # development version
vudl_net_21_path = Path("/home/matthew/university_work/02_neural_networks_python/01_phd_script_approach/04_vgg16/03_best_2head_synth_real_v2/01_vgg16_2head_block5_trained.h5")

LiCSAlert_settings = {"n_baseline_end" : 35,                                    # n_ifgs that are used in the baseline stage (i.e. by ICASAR)
                      "out_folder" : "LiCSAlert_01_Sierra_Negra_no_intermediate",    # no spaces, snake or camel case
                      "run_ICASAR" : True,                                # If False, attempt to load results from previous run.  If True, run (which can be slow)
                      "intermediate_figures" : False,                      # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "downsample_run" : 0.5,                              # data can be downsampled to speed things up
                      "downsample_plot" : 0.5}                             # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      
ICASAR_settings = {"n_comp" : 6,                                    # number of components to recover with ICA (ie the number of PCA sources to keep)
                    "bootstrapping_param" : (200, 0),               # number of runs with bootstrapping, number of runs without bootstrapping
                    "hdbscan_param" : (35, 10),                     # (min_cluster_size, min_samples)
                    "tsne_param" : (30, 12),                        # (perplexity, early_exaggeration)
                    "ica_param" : (1e-2, 150),                      # (tolerance, max iterations)
                    "figures" : "png",                       # if png, saved in a folder as .png.  If window, open as interactive matplotlib figures, if window+png then both.  
                    "create_all_ifgs_flag" : False}                 # Creates all possible pairs of ifgs between all acquisitions.  Results can be more complex with this set to True, but also better at recovering small magnitude signals
                    


LiCSAlert_batch_mode(displacement_r2, ICASAR_settings = ICASAR_settings, **LiCSAlert_settings, ICASAR_path = ICASAR_path, vudl_net_21_path = vudl_net_21_path)


#%% Example 2, creating the LiCSAlert figure at all time steps

LiCSAlert_settings["intermediate_figures"] = True                # change one key/variable in the dict so that all the intermediate ifgs are made by LiCSAlert (which can be slow, but is good for animations)
LiCSAlert_settings["out_folder"] = "LiCSAlert_02_Sierra_Negra_intermediate"
LiCSAlert_batch_mode(displacement_r2_copy, ICASAR_settings = ICASAR_settings, **LiCSAlert_settings,ICASAR_path = ICASAR_path, vudl_net_21_path = vudl_net_21_path)


#%%


    
    