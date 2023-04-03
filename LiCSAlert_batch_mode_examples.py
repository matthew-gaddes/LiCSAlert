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
from licsalert.licsalert import LiCSAlert_batch_mode

#ICASAR_path = Path("/home/matthew/university_work/15_my_software_releases/ICASAR-2.7.2/")                               # location of ICASAR functions
ICASAR_path = Path("/home/matthew/university_work/crucial_1000gb/01_blind_signal_separation_python/13_ICASAR/ICASAR_GitHub/")           # development version




#%% Example 3, Running LiCSAlert with a smaller signal (Campi Flegrei, processed with LiCSBAS)

LiCSBAS_out_folder_campi_flegrei = Path('./022D_04826_121209')
sys.path.append(str(ICASAR_path))                                                                       # Add ICASAR to the path so we can use one of its functions.  
import icasar
from icasar.icasar_funcs import LiCSBAS_to_ICASAR

LiCSAlert_settings = {"n_baseline_end" : 55,                                         # n_ifgs that are used in the baseline stage (i.e. by ICASAR)
                      "out_folder" : Path("LiCSAlert_campi_flegrei"),                # pathlib Path
                      "run_ICASAR" : False,                                           # If False, attempt to load results from previous run.  If True, run (which can be slow)
                      "figure_intermediate" : False,                                # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "figure_type"         : 'png',                                 # either 'window' or 'png' (to save as pngs)
                      "downsample_run" : 0.5,                                        # data can be downsampled to speed things up
                      "downsample_plot" : 0.5,                                       # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "residual_type"        : 'cumulative'}                              # controls the type of residual used in the lower plot.  Either cumulative or window   


ICASAR_settings = {"n_comp" : 5,                                         # number of components to recover with ICA (ie the number of PCA sources to keep)
                   "bootstrapping_param" : (200, 0),                    # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                   "tsne_param" : (30, 12),                             # (perplexity, early_exaggeration)
                   "ica_param" : (1e-2, 150),                           # (tolerance, max iterations)
                   "hdbscan_param" : (100,10),                           # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
                   "ifgs_format"        : 'cum',
                   "load_fastICA_results" : True,                      # If all the FastICA runs already exisit, setting this to True speeds up ICASAR as they don't need to be recomputed.  
                   "figures" : "png+window"}                            # if png, saved in a folder as .png.  If window, open as interactive matplotlib figures,


displacement_r2, tbaseline_info, _ = LiCSBAS_to_ICASAR(LiCSBAS_out_folder_campi_flegrei, figures=True)        # open various LiCSBAS products, spatial ones in displacement_r2, temporal ones in tbaseline_info
displacement_r2['ifg_dates'] = tbaseline_info['ifg_dates']                                                  # Unlike ICASAR, LiCSAlert always needs the ifg_dates too.  

LiCSAlert_batch_mode(displacement_r2, ICASAR_settings = ICASAR_settings, **LiCSAlert_settings, ICASAR_path = ICASAR_path)
    
    