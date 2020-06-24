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
import shutil
import pickle


# LiCSAlert imports
sys.path.append("./lib")
from LiCSAlert_functions import LiCSAlert, LiCSAlert_figure, save_pickle, shorten_LiCSAlert_data, LiCSAlert_preprocessing
from downsample_ifgs import downsample_ifgs
from LiCSAlert_aux_functions import col_to_ma

# Which includes ICASAR
sys.path.append("/home/matthew/university_work/01_blind_signal_separation_python/13_ICASAR/ICASAR_GitHub")                  # location of ICASAR functions
from ICASAR_functions import ICASAR


sys.path.append("/home/matthew/university_work/python_stuff/python_scripts")
from small_plot_functions import matrix_show



#%% Things to set (and load the data)

LiCSAlert_settings = {"volcano_name" : "Sierra Negra",
                      "run_ICASAR" : True,                          # If False, attempt to load results from previous run
                      "n_baseline_end" : 35,                        # n_ifgs that are used in the baseline stage (i.e. by ICASAR)
                      "verbose" : True,
                      "downsample_run" : 1.,                       # data can be downsampled to speed things up
                      "downsample_plot" : 0.3,                      # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "out_folder" : 'LiCSAlert_results'}           # figures produced by LiCSAlert will be saved here

ICASAR_settings = {"n_comp" : 6,                                    # number of components to recover with ICA (ie the number of PCA sources to keep)
                    "bootstrapping_param" : (200, 0),
                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                    "tsne_param" : (30, 12),                       # (perplexity, early_exaggeration)
                    "ica_param" : (1e-2, 150),                     # (tolerance, max iterations)
                    "ge_kmz"    :  True,                            # make a google earth .kmz of the ICs
                    "figures" : "png+window"}                       # if png, saved in a folder as .png.  If window, open as interactive matplotlib figures,


displacement_r2 = {}                                                            # initiate
with open(f"sierra_negra_example_data.pkl", 'rb') as f:
    phUnw_files = pickle.load(f)    
    displacement_r2["incremental"] = pickle.load(f)    
    displacement_r2["mask"] = pickle.load(f)    
    cumulative_baselines = pickle.load(f)    
    acq_dates = pickle.load(f)
    lons = pickle.load(f)    
    lats = pickle.load(f)    
f.close()



#%% Preliminary steps

# create a folder that will be used for outputs
try:
    print(f"Trying to remove the existing outputs folder ({LiCSAlert_settings['out_folder']})... ", end = '')
    shutil.rmtree(LiCSAlert_settings['out_folder'])                                                                       # try to remove folder
    print('Done!')
except:
    print("Failed!")                                                                                # 
try:
    print(f"Trying to create a new outputs folder ({LiCSAlert_settings['out_folder']})... ", end = '')                                    # try to make a new folder
    os.mkdir(LiCSAlert_settings['out_folder'])                                                                       
    print('Done!')
except:
    print("Failed!") 


n_baseline_end = LiCSAlert_settings["n_baseline_end"]
LiCSAlert_settings["volcano_name_filename"]  = "_".join(LiCSAlert_settings["volcano_name"].split(" "))                     # make a version of the name with no spaces that can be used to save files





#%% run ICASAR, or load the results from a previous run
displacement_r2 = LiCSAlert_preprocessing(displacement_r2, LiCSAlert_settings["downsample_run"], 
                                                           LiCSAlert_settings["downsample_plot"], LiCSAlert_settings["verbose"])                    # mean centre and downsize the data

if LiCSAlert_settings["run_ICASAR"]:
    sources, tcs, residual, Iq, n_clusters, S_all_info, means = ICASAR(phUnw = displacement_r2['incremental'][:n_baseline_end],
                                                                       mask = displacement_r2['mask'], lons = lons, lats = lats,
                                                                       **ICASAR_settings)

    sources_downsampled, _ = downsample_ifgs(sources, displacement_r2["mask"], LiCSAlert_settings["downsample_plot"])                     # downsample for plots
    

    
else:
    try:
        with open(f"ICASAR_results/ICASAR_results.pkl", 'rb') as f:
            sources = pickle.load(f)    
            tcs  = pickle.load(f)    
            source_residuals = pickle.load(f)    
            Iq_sorted = pickle.load(f)    
            n_clusters = pickle.load(f)    
            
        sources_downsampled, _ = downsample_ifgs(sources, displacement_r2["mask"], LiCSAlert_settings["downsample_plot"])                     # downsample the sources as this can speed up plotting

        
    except:
        raise Exception(f"Unable to open the results of ICASAR (which are usually stored in 'ICASAR_results') "
                        f"Try re-running and enabling ICASAR with 'run_ICASAR' set to 'True'.  ")


#%% Either loop through LiCSAlert as each new ifg is ingested (as if it were becoming available), plotting each time

for ifg_n in np.arange(n_baseline_end+1, displacement_r2["incremental"].shape[0]+1):
    
    displacement_r2_current = shorten_LiCSAlert_data(displacement_r2, n_end=ifg_n)                        # get the current ifgs to simulate monitoring
    cumulative_baselines_current = cumulative_baselines[:ifg_n]                                                             # also get current time values


    sources_tcs_monitor, residual_monitor = LiCSAlert(sources, cumulative_baselines_current, displacement_r2_current["incremental"][:n_baseline_end], 
                                                                                    displacement_r2_current["incremental"][n_baseline_end:], t_recalculate=10)    

    LiCSAlert_figure(sources_tcs_monitor, residual_monitor, sources_downsampled, displacement_r2_current, n_baseline_end, 
                      cumulative_baselines_current, time_value_end=cumulative_baselines[-1], out_folder = f"{LiCSAlert_settings['out_folder']}",
                      day0_date = acq_dates[0], sources_downsampled = True)                                                                                 # main LiCSAlert figure, note that we use downsampled sources to speed things up



#%% Or just run on the entire set of data
    
n_end = 97

displacement_r2_current = shorten_LiCSAlert_data(displacement_r2, n_end)                              # get the current ifgs to simulate monitoring
cumulative_baselines_current = cumulative_baselines[:n_end]                                                             # also get current time values


sources_tcs_monitor, residual_monitor = LiCSAlert(sources, cumulative_baselines_current, displacement_r2_current["incremental"][:n_baseline_end], 
                                                  displacement_r2_current["incremental"][n_baseline_end:], t_recalculate=10)    

LiCSAlert_figure(sources_tcs_monitor, residual_monitor, sources, displacement_r2_current, n_baseline_end, 
                  cumulative_baselines_current, time_value_end=cumulative_baselines[-1], day0_date = acq_dates[0], sources_downsampled = False)                 # note that as we're only plotting once, 
                                                                                                                                                                # this time we just use the full resolution sources, and set the sources_downsampled flag appropraitely (ie to False)




#%%

