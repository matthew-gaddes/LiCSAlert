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

sys.path.append("./lib")
from LiCSAlert_functions import LiCSAlert_batch_mode

#%% Load Sentinel-1 data for Sierra Negra

displacement_r2 = {}                                                                    # initiate a dictionary to store some info in.  
with open(f"sierra_negra_example_data.pkl", 'rb') as f:
    phUnw_files = pickle.load(f)                                                        # these are the names of the interferograms used (e.g. YYYYMMDD_YYYYMMDD.unw)
    displacement_r2["incremental"] = pickle.load(f)                                     # the incremental interferograms as row vectors.
    displacement_r2["mask"] = pickle.load(f)                                            # a mask to convert a row vector back to a rank 2 array (col_to_ma in LiCSAlert_functions.py does this easily.  )
    cumulative_baselines = pickle.load(f)                                               # the cumulative baselines. ie if there are acquisitions every 12 days, these would be 12,24,36 etc.  
    acq_dates = pickle.load(f)                                                          # acquisition dates.  This should be one longer than the names of the interferograms
    displacement_r2['lons'] = pickle.load(f)                                            # matrix of longitudes for each pixel.  Should be the same size as mask
    displacement_r2['lats'] = pickle.load(f)                                            # matrix of latitdues for each pixel.  Should be the same size as mask
f.close()

displacement_r2_copy = copy.deepcopy(displacement_r2)                                   # make a copy for use with the 2nd example  


#%% Example 1, creating one LiCSAlert figure for the complete time series
ICASAR_path = Path("/home/matthew/university_work/01_blind_signal_separation_python/13_ICASAR/ICASAR-1.0")                  # location of ICASAR functions

LiCSAlert_settings = {"n_baseline_end" : 35,                              # n_ifgs that are used in the baseline stage (i.e. by ICASAR)
                      "out_folder" : "01_Sierra_Negra_no_intermediate",    # no spaces, snake or camel case
                      "run_ICASAR" : False,                                # If False, attempt to load results from previous run.  If True, run (which can be slow)
                      "intermediate_figures" : False,                      # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "downsample_run" : 0.5,                              # data can be downsampled to speed things up
                      "downsample_plot" : 0.5}                             # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      
ICASAR_settings = {"n_comp" : 6,                                    # number of components to recover with ICA (ie the number of PCA sources to keep)
                    "bootstrapping_param" : (200, 0),               # number of runs with bootstrapping, number of runs without bootstrapping
                    "hdbscan_param" : (35, 10),                     # (min_cluster_size, min_samples)
                    "tsne_param" : (30, 12),                        # (perplexity, early_exaggeration)
                    "ica_param" : (1e-2, 150),                      # (tolerance, max iterations)
                    "ge_kmz"    :  True,                            # make a google earth .kmz of the ICs
                    "figures" : "png+window"}                       # if png, saved in a folder as .png.  If window, open as interactive matplotlib figures, if window+png then both.  
                    


LiCSAlert_batch_mode(displacement_r2, cumulative_baselines, acq_dates,
                      ICASAR_settings = ICASAR_settings, **LiCSAlert_settings, ICASAR_path = ICASAR_path)


#%% Example 2, creating the LiCSAlert figure at all time steps

LiCSAlert_settings = {"n_baseline_end" : 35,                              # n_ifgs that are used in the baseline stage (i.e. by ICASAR)
                      "out_folder" : "02_Sierra_Negra_intermediate",    # no spaces, snake or camel case
                      "run_ICASAR" : False,                                # If False, attempt to load results from previous run.  If True, run (which can be slow)
                      "intermediate_figures" : True,                     # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "downsample_run" : 0.5,                             # data can be downsampled to speed things up
                      "downsample_plot" : 0.5}                            # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      
ICASAR_settings = {"n_comp" : 6,                                    # number of components to recover with ICA (ie the number of PCA sources to keep)
                    "bootstrapping_param" : (200, 0),               # number of runs with bootstrapping, number of runs without bootstrapping
                    "hdbscan_param" : (35, 10),                     # (min_cluster_size, min_samples)
                    "tsne_param" : (30, 12),                        # (perplexity, early_exaggeration)
                    "ica_param"  : (1e-2, 150),                      # (tolerance, max iterations)
                    "ge_kmz"     :  True,                            # make a google earth .kmz of the ICs
                    "figures"    : "png"}                       # if png, saved in a folder as .png.  If window, open as interactive matplotlib figures, if window+png then both.  


LiCSAlert_batch_mode(displacement_r2_copy, cumulative_baselines, acq_dates,
                     ICASAR_settings = ICASAR_settings, **LiCSAlert_settings,ICASAR_path = ICASAR_path)

#%%