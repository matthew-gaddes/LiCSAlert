#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 10:55:48 2020

@author: matthew


To do:
    - update colorbar ticks in results explorer
    - check results explorere works with even number of sources.  
    - colour of labels (blakc, orange, purple).  
"""

import sys
import pickle
from pathlib import Path
import copy
import numpy as np
import matplotlib.pyplot as plt
import pdb

import licsalert
from licsalert.monitoring_functions import LiCSAlert_monitoring_mode
from licsalert.plotting import licsalert_results_explorer
from licsalert.licsalert import reconstruct_ts_from_dir






#%% Example 01: Campi Flegrei, only the final figure.  


licsalert_settings = {"baseline_end"        : "20170101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
                      "figure_intermediate" : False,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "figure_type"         : 'png',                             # either 'window' or 'png' (to save as pngs), or 'both'
                      "downsample_run"      : 0.5,                                     # data can be downsampled to speed things up
                      "downsample_plot"     : 0.5,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "residual_type"       : 'cumulative',                      # controls the type of residual used in the lower plot.  Either cumulative or window   
                      "t_recalculate"       : 40,                               # Number of acquisitions that the lines of best fit are calcualted over.  Larger value makes the algorithm more sensitive
                      'inset_ifgs_scaling'  : 15}                               # scales the size of the incremental and cumulative ifgs in the top row of the figure.  Smaller values gives a bigger figures.  


icasar_settings = {"n_comp"                  : 5,                                                  # number of components to recover with ICA (ie the number of PCA sources to keep)
                   "bootstrapping_param"    : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                    "tsne_param"             : (30, 12),                                       # (perplexity, early_exaggeration)
                    "ica_param"              : (1e-2, 150),                                     # (tolerance, max iterations)
                    "hdbscan_param"          : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
                    "ifgs_format"            : 'cum',                                  # can be 'all', 'inc' (incremental - short temporal baselines), or 'cum' (cumulative - relative to first acquisition)
                    "sica_tica"              : 'tica' }

licsbas_settings = {"filtered"               : False,
                    "date_start"            : None,
                    "date_end"              : None,
                    'mask_type'             : 'licsbas',
                    'crop_pixels'           : None}


outdir = Path("./")
volcano = '001_campi_flegrei_example'                                                 # outdir final                                                                                                                      
licsbas_dir = Path("./022D_04826_121209_campi_flegrei")                                         # input data


LiCSAlert_monitoring_mode(outdir = outdir, region = None, volcano = volcano,
                          licsbas_dir = licsbas_dir,
                          licsalert_settings = licsalert_settings, 
                          icasar_settings = icasar_settings,
                          licsbas_settings = licsbas_settings)

licsalert_out_dir = outdir / volcano

licsalert_results_explorer(outdir / volcano, fig_width = 18)                                                 # use this function to explore the results


#ics_one_hot = [1, 1, 1, 1]                                                                   # One hot encoding of which sources to use in the reconstruction.  1 means used, 0 means not.  list must be the same length as the number of ICs.  
#X_inc_r3, X_cum_r3 = reconstruct_ts_from_dir(ics_one_hot, outdir / volcano)                               # return the cumualtive interferograms reconstrutced using the ICs selected above.  All mean centering has been removed.  




#%% Example 2: make all intermediate figures



licsalert_settings = {"baseline_end"        : "20170101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
                      ### changed here
                      "figure_intermediate" : True,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      ### end change
                      "figure_type"         : 'png',                             # either 'window' or 'png' (to save as pngs), or 'both'
                      "downsample_run"      : 0.5,                                     # data can be downsampled to speed things up
                      "downsample_plot"     : 0.5,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "residual_type"       : 'cumulative',                      # controls the type of residual used in the lower plot.  Either cumulative or window   
                      "t_recalculate"       : 40,                               # Number of acquisitions that the lines of best fit are calcualted over.  Larger value makes the algorithm more sensitive
                      'inset_ifgs_scaling'  : 15}                               # scales the size of the incremental and cumulative ifgs in the top row of the figure.  Smaller values gives a bigger figures.  


icasar_settings = {"n_comp"                  : 5,                                                  # number of components to recover with ICA (ie the number of PCA sources to keep)
                   "bootstrapping_param"    : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                    "tsne_param"             : (30, 12),                                       # (perplexity, early_exaggeration)
                    "ica_param"              : (1e-2, 150),                                     # (tolerance, max iterations)
                    "hdbscan_param"          : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
                    "ifgs_format"            : 'cum',                                  # can be 'all', 'inc' (incremental - short temporal baselines), or 'cum' (cumulative - relative to first acquisition)
                    "sica_tica"              : 'sica' }

licsbas_settings = {"filtered"               : False,
                    "date_start"            : None,
                    "date_end"              : None,
                    'mask_type'             : 'licsbas',
                    'crop_pixels'           : None}


outdir = Path("./")
volcano = '002_campi_flegrei_example_all_figures'                                                 # outdir final                                                                                                                      
licsbas_dir = Path("./022D_04826_121209_campi_flegrei")                                         # input data


LiCSAlert_monitoring_mode(outdir = outdir, region = None, volcano = volcano,
                          licsbas_dir = licsbas_dir,
                          licsalert_settings = licsalert_settings, 
                          icasar_settings = icasar_settings,
                          licsbas_settings = licsbas_settings)


    
#%% Example 3: temporal ICA at Vesuvius

outdir = Path("./")
volcano = '003_vesuvius_example_tica'                                                 # outdir final                                                                                                                      
licsbas_dir = Path("./022D_04826_121209_vesuvius_crop_rationalized")                                         # input data


licsalert_settings = {"baseline_end"        : "20170101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
                      "figure_intermediate" : False,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "figure_type"         : 'png',                             # either 'window' or 'png' (to save as pngs), or 'both'
                      "downsample_run"      : 0.5,                                     # data can be downsampled to speed things up
                      "downsample_plot"     : 0.5,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "t_recalculate"       : 40,                                   # Number of acquisitions that the lines of best fit are calcualted over.  Larger value makes the algorithm more sensitive
                      "residual_type"       : 'cumulative',                      # controls the type of residual used in the lower plot.  Either cumulative or window   
                      'inset_ifgs_scaling'  : 15}                               # scales the size of the incremental and cumulative ifgs in the top row of the figure.  Smaller values gives a bigger figures.  

icasar_settings = {"n_comp"                 : 5,                                                  # number of components to recover with ICA (ie the number of PCA sources to keep)
                    "bootstrapping_param"   : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                    "tsne_param"            : (30, 12),                                       # (perplexity, early_exaggeration)
                    "ica_param"             : (1e-2, 150),                                     # (tolerance, max iterations)
                    "hdbscan_param"         : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
                    #"ifgs_format"        : 'cum',                                  # ifg format is redundant when using tICA as cumulative signals must be considered when working temporally.  
                    ### Note set to tica here for temporal:
                    "sica_tica"             : 'tica' }

licsbas_settings = {"filtered"               : False,
                    "date_start"            : '20141218',
                    "date_end"              : '20211217',
                    'mask_type'             : 'licsbas',
                    'crop_pixels'           : None}


LiCSAlert_monitoring_mode(outdir = outdir, region = None, volcano = volcano,
                          licsbas_dir = licsbas_dir,
                          licsalert_settings = licsalert_settings, 
                          icasar_settings = icasar_settings,
                          licsbas_settings = licsbas_settings)


licsalert_results_explorer(outdir / volcano, fig_width = 18)                                                 # use this function to explore the results


#%% Example 4: spatial ICA at Vesuvius

outdir = Path("./")
volcano = '004_vesuvius_example_sica'                                                 
licsbas_dir = Path("./022D_04826_121209_vesuvius_crop_rationalized")                                         # input data



licsalert_settings = {"baseline_end"        : "20170101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
                      "figure_intermediate" : False,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "figure_type"         : 'png',                             # either 'window' or 'png' (to save as pngs), or 'both'
                      "downsample_run"      : 0.5,                                     # data can be downsampled to speed things up
                      "downsample_plot"     : 0.5,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "t_recalculate"       : 40,                                   # Number of acquisitions that the lines of best fit are calcualted over.  Larger value makes the algorithm more sensitive
                      "residual_type"       : 'cumulative',                      # controls the type of residual used in the lower plot.  Either cumulative or window   
                      'inset_ifgs_scaling'  : 15}                               # scales the size of the incremental and cumulative ifgs in the top row of the figure.  Smaller values gives a bigger figures.  

icasar_settings = {"n_comp"                 : 5,                                                  # number of components to recover with ICA (ie the number of PCA sources to keep)
                    "bootstrapping_param"   : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                    "tsne_param"            : (30, 12),                                       # (perplexity, early_exaggeration)
                    "ica_param"             : (1e-2, 150),                                     # (tolerance, max iterations)
                    "hdbscan_param"         : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
                    ### begin change from example 3
                    "ifgs_format"        : 'cum',                                  # ifg format is redundant when using tICA as cumulative signals must be considered when working temporally.  
                    "sica_tica"             : 'sica' }
                    ### end change from example 3

licsbas_settings = {"filtered"               : False,
                    "date_start"            : '20141218',
                    "date_end"              : '20211217',
                    'mask_type'             : 'licsbas',
                    'crop_pixels'           : None}




LiCSAlert_monitoring_mode(outdir = outdir, region = None, volcano = volcano,
                          licsbas_dir = licsbas_dir,
                          licsalert_settings = licsalert_settings, 
                          icasar_settings = icasar_settings,
                          licsbas_settings = licsbas_settings)


licsalert_results_explorer(outdir / volcano, fig_width = 18)                                                 # use this function to explore the results