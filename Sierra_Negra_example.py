#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 10:55:48 2020

Example for Github with Sierra Negra Data

@author: matthew
"""


import numpy as np
import numpy.ma as ma
import sys
import os

import pickle


# LiCSAlert imports
sys.path.append("./lib")
from LiCSAlert_functions import LiCSAlert, LiCSAlert_figure, save_pickle, shorten_LiCSAlert_data, LiCSAlert_preprocessing
from downsample_ifgs import downsample_ifgs
from LiCSAlert_aux_functions import col_to_ma

# Which includes ICASAR
sys.path.append("/home/matthew/university_work/01_blind_signal_separation_python/13_ICASAR/ICASAR_GitHub")                  # location of ICASAR functions
from ICASAR_functions import ICASAR


# sys.path.append("/home/matthew/university_work/python_stuff/python_scripts")
# from small_plot_functions import matrix_show, low_resolution_ifgs



#%% Things to set

LiCSAlert_settings = {"volcano_name" : "Sierra Negra",
                      "run_ICASAR" : True,                          # If False, attempt to load results from previous run
                      "n_baseline_end" : 35,                        # n_ifgs that are used in the baseline stage (i.e. by ICASAR)
                      "verbose" : True,
                      "downsample_run" : 1.,                       # data can be downsampled to speed things up
                      "downsample_plot" : 0.3}                      # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound

ICASAR_settings = {"n_comp" : 6,                                    # number of components to recover with ICA (ie the number of PCA sources to keep)
                    "bootstrapping_param" : (200, 0),
                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                    "tsne_param" : (30, 12),                       # (perplexity, early_exaggeration)
                    "ica_param" : (1e-2, 150)}                     # (tolerance, max iterations)


displacement_r2 = {}                                                            # initiate
with open(f"sierra_negra_example_data.pkl", 'rb') as f:
    phUnw_files = pickle.load(f)    
    displacement_r2["incremental"] = pickle.load(f)    
    displacement_r2["mask"] = pickle.load(f)    
    time_values = pickle.load(f)    
f.close()


#%% Preliminary steps

displacement_r2 = LiCSAlert_preprocessing(displacement_r2, LiCSAlert_settings["downsample_run"], 
                                                           LiCSAlert_settings["downsample_plot"], LiCSAlert_settings["verbose"])
n_baseline_end = LiCSAlert_settings["n_baseline_end"]
LiCSAlert_settings["volcano_name_filename"]  = "_".join(LiCSAlert_settings["volcano_name"].split(" "))                     # make a version of the name with no spaces that can be used to save files

try:
    os.mkdir(f"png_{LiCSAlert_settings['volcano_name_filename']}")                                                              # create a directory for pngs output by LiCSAlert
except:
    print(f"The output folder probably already exisits (png_{LiCSAlert_settings['volcano_name_filename']}). "
          f"The contents will be overwritten.  ")    

#%% run ICASAR

if LiCSAlert_settings["run_ICASAR"]:
    sources, tcs, residual, Iq, n_clusters, S_all_info = ICASAR(phUnw = displacement_r2['incremental'][:n_baseline_end],
                                                                 mask = displacement_r2['mask'], **ICASAR_settings)

    sources_downsampled, _ = downsample_ifgs(sources, displacement_r2["mask"], LiCSAlert_settings["downsample_plot"])                     # downsample for plots
    
    save_pickle(f"{LiCSAlert_settings['volcano_name_filename']}_ICASAR", sources, sources_downsampled, tcs, residual, Iq, n_clusters) 
    del tcs, residual, S_all_info, n_clusters                                                   # these aren't needed if then running LiCSAlert
    if LiCSAlert_settings["verbose"]:
        print(f"The resuls of ICASAR have been saved for future use in '{LiCSAlert_settings['volcano_name']}_ICASAR.pkl'.  ")
    
else:
    try:
        with open(f"{LiCSAlert_settings['volcano_name_filename']}_ICASAR.pkl", 'rb') as f:
            sources = pickle.load(f)    
            sources_downsampled = pickle.load(f)    
    except:
        raise Exception(f"Unable to find the results of ICASAR for {LiCSAlert_settings['volcano_name_filename']} (which are stored in {LiCSAlert_settings['volcano_name']}_ICASAR.pkl) "
                        f"Try re-running and enabling ICASAR with 'run_ICASAR' set to 'True'.  ")


#%% Loop through LiCSAlert as each new ifg is ingested, plotting each time

for ifg_n in np.arange(n_baseline_end+1, displacement_r2["incremental"].shape[0]+1):
    
    displacement_r2_current = shorten_LiCSAlert_data(displacement_r2, n_end=ifg_n)                        # get the current ifgs to simulate monitoring
    time_values_current = time_values[:ifg_n]                                                             # also get current time values


    sources_tcs_monitor, residual_monitor = LiCSAlert(sources, time_values_current, displacement_r2_current["incremental"][:n_baseline_end], 
                                                                                    displacement_r2_current["incremental"][n_baseline_end:], t_recalculate=10)    

    LiCSAlert_figure(sources_tcs_monitor, residual_monitor, sources_downsampled, displacement_r2_current, n_baseline_end, 
                      time_values_current, time_value_end=time_values[-1], out_folder = f"png_{LiCSAlert_settings['volcano_name_filename']}")


#%% Or just run on the entire set of data
    
n_end = 97

displacement_r2_current = shorten_LiCSAlert_data(displacement_r2, n_end)                              # get the current ifgs to simulate monitoring
time_values_current = time_values[:n_end]                                                             # also get current time values


sources_tcs_monitor, residual_monitor = LiCSAlert(sources, time_values_current, displacement_r2_current["incremental"][:n_baseline_end], 
                                                                                displacement_r2_current["incremental"][n_baseline_end:], t_recalculate=10)    

LiCSAlert_figure(sources_tcs_monitor, residual_monitor, sources_downsampled, displacement_r2_current, n_baseline_end, 
                  time_values_current, time_value_end=time_values[-1])    


