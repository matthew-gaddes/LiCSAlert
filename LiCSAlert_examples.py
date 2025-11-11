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
from licsalert.licsalert import reconstruct_ts_from_dir



#%% Example 001:  Cordon Caulle with time-varying coherence

# licsalert_settings = {"baseline_end"        : "20210101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
#                       "figure_intermediate" : False,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
#                       "figure_type"         : 'png',                             # either 'window' or 'png' (to save as pngs), or 'both'
#                       "downsample_run"      : 1.0,                                     # data can be downsampled to speed things up
#                       "downsample_plot"     : 0.5,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
#                       "residual_type"       : 'cumulative',                      # controls the type of residual used in the lower plot.  Either cumulative or window   
#                       "t_recalculate"       : 40,                               # Number of acquisitions that the lines of best fit are calcualted over.  Larger value makes the algorithm more sensitive
#                       'inset_ifgs_scaling'  : 15}                               # scales the size of the incremental and cumulative ifgs in the top row of the figure.  Smaller values gives a bigger figures.  


# icasar_settings = {"n_pca_comp_start"       : 6,                                                  
#                    "n_pca_comp_stop"        : 8,                                                  
#                    "bootstrapping_param"    : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
#                     "tsne_param"             : (30, 12),                                       # (perplexity, early_exaggeration)
#                     "ica_param"              : (1e-2, 150),                                     # (tolerance, max iterations)
#                     "hdbscan_param"          : (35,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
#                     "ifgs_format"            : 'cum',                                  # can be 'all', 'inc' (incremental - short temporal baselines), or 'cum' (cumulative - relative to first acquisition)
#                     "load_fastICA_results"   : True}

# licsbas_settings = {"filtered"               : False,
#                     "date_start"            : None,
#                     "date_end"              : None,
#                     'mask_type'             : 'licsbas',                        # "dem" or "licsbas"
#                     'crop_pixels'           : None}


# outdir = Path("./")
# volcano = '001_cordon_caulle_example' 

# # Load the data
# data_as_arg={}

# with open(Path('./example_data/cordon_caulle_ts_164A.pkl'), 'rb') as f:
#     data_as_arg['displacement_r3'] = pickle.load(f)   
#     data_as_arg['tbaseline_info'] = pickle.load(f)
    

# LiCSAlert_monitoring_mode(
#     outdir = outdir, region = None, volcano = volcano,
#     data_as_arg=data_as_arg,
#     licsalert_settings = licsalert_settings, 
#     icasar_settings = icasar_settings,
#     licsbas_settings = licsbas_settings,
#     )


# #%% Example 002: Campi Flegrei, only the final figure.  


# licsalert_settings = {"baseline_end"        : "20170101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
#                       "figure_intermediate" : False,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
#                       "figure_type"         : 'png',                             # either 'window' or 'png' (to save as pngs), or 'png'
#                       "downsample_run"      : 0.5,                                     # data can be downsampled to speed things up
#                       "downsample_plot"     : 0.2,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
#                       "residual_type"       : 'cumulative',                      # controls the type of residual used in the lower plot.  Either cumulative or window   
#                       "t_recalculate"       : 40,                               # Number of acquisitions that the lines of best fit are calcualted over.  Larger value makes the algorithm more sensitive
#                       'inset_ifgs_scaling'  : 15}                               # scales the size of the incremental and cumulative ifgs in the top row of the figure.  Smaller values gives a bigger figures.  


# icasar_settings = {"n_pca_comp_start"       : 6,                                                  
#                    "n_pca_comp_stop"        : 7,                                                  
#                    "bootstrapping_param"    : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
#                     "tsne_param"             : (30, 12),                                       # (perplexity, early_exaggeration)
#                     "ica_param"              : (1e-2, 150),                                     # (tolerance, max iterations)
#                     "hdbscan_param"          : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
#                     "ifgs_format"            : 'cum',                                  # can be 'all', 'inc' (incremental - short temporal baselines), or 'cum' (cumulative - relative to first acquisition)
#                     "load_fastICA_results"   : True}

# licsbas_settings = {"filtered"               : False,
#                     "date_start"            : None,
#                     "date_end"              : None,
#                     'mask_type'             : 'licsbas',                        # "dem" or "licsbas"
#                     'crop_pixels'           : None}


# outdir = Path("./")
# volcano = '002_campi_flegrei_example'                                                 # outdir final                                                                                                                      
# licsbas_dir = Path("./example_data/022D_04826_121209_campi_flegrei")                                         # input data


# LiCSAlert_monitoring_mode(
#     outdir = outdir, region = None, volcano = volcano,
#     licsbas_dir = licsbas_dir,
#     licsalert_settings = licsalert_settings, 
#     icasar_settings = icasar_settings,
#     licsbas_settings = licsbas_settings,
#     )


#%% Example 2: make all intermediate figures



licsalert_settings = {"baseline_end"        : "20170101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
                      ### changed here
                      "figure_intermediate" : True,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      ### end change
                      "figure_type"         : 'png',                             # either 'window' or 'png' (to save as pngs), or 'both'
                      "downsample_run"      : 0.5,                                     # data can be downsampled to speed things up
                      "downsample_plot"     : 0.25,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "residual_type"       : 'cumulative',                      # controls the type of residual used in the lower plot.  Either cumulative or window   
                      "t_recalculate"       : 40,                               # Number of acquisitions that the lines of best fit are calcualted over.  Larger value makes the algorithm more sensitive
                      'inset_ifgs_scaling'  : 15}                               # scales the size of the incremental and cumulative ifgs in the top row of the figure.  Smaller values gives a bigger figures.  


icasar_settings = {"n_pca_comp_start"       : 6,                                                  
                   "n_pca_comp_stop"        : 7,                                                  
                   "bootstrapping_param"    : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                    "tsne_param"             : (30, 12),                                       # (perplexity, early_exaggeration)
                    "ica_param"              : (1e-2, 150),                                     # (tolerance, max iterations)
                    "hdbscan_param"          : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
                    "ifgs_format"            : 'cum',                                  # can be 'all', 'inc' (incremental - short temporal baselines), or 'cum' (cumulative - relative to first acquisition)
                    }

licsbas_settings = {"filtered"               : False,
                    "date_start"            : None,
                    "date_end"              : None,
                    'mask_type'             : 'licsbas',
                    'crop_pixels'           : None}


outdir = Path("./")
volcano = '003_campi_flegrei_example_all_figures'                                                 # outdir final                                                                                                                      
licsbas_dir = Path("./example_data/022D_04826_121209_campi_flegrei")                                         # input data


LiCSAlert_monitoring_mode(
    outdir = outdir, region = None, volcano = volcano,
    licsbas_dir = licsbas_dir,
    licsalert_settings = licsalert_settings, 
    icasar_settings = icasar_settings,
    licsbas_settings = licsbas_settings
    )



#%% 004: AlginSAR datacube


licsalert_settings = {"baseline_end"        : "20220101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
                      "figure_intermediate" : False,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "figure_type"         : 'png',                             # either 'window' or 'png' (to save as pngs), or 'both'
                      "downsample_run"      : 0.5,                                     # data can be downsampled to speed things up
                      "downsample_plot"     : 0.5,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "residual_type"       : 'cumulative',                      # controls the type of residual used in the lower plot.  Either cumulative or window   
                      "t_recalculate"       : 40,                               # Number of acquisitions that the lines of best fit are calcualted over.  Larger value makes the algorithm more sensitive
                      'inset_ifgs_scaling'  : 15}                               # scales the size of the incremental and cumulative ifgs in the top row of the figure.  Smaller values gives a bigger figures.  


icasar_settings = {"n_pca_comp_start"       : 6,                                                  
                   "n_pca_comp_stop"        : 7,                                                  
                   "bootstrapping_param"    : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                    "tsne_param"             : (30, 12),                                       # (perplexity, early_exaggeration)
                    "ica_param"              : (1e-2, 150),                                     # (tolerance, max iterations)
                    "hdbscan_param"          : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
                    "ifgs_format"            : 'cum',                                  # can be 'all', 'inc' (incremental - short temporal baselines), or 'cum' (cumulative - relative to first acquisition)
                    "sica_tica"              : 'sica',
                    "load_fastICA_results"   : True}

licsbas_settings = {"filtered"               : False,
                    "date_start"            : None,
                    "date_end"              : None,
                    'mask_type'             : 'licsbas',                        # "dem" or "licsbas"
                    'crop_pixels'           : None}


outdir = Path("./004_AlignSAR_data_cube_example")
volcano = 'alignsar_campi_flegrei'                                                 


LiCSAlert_monitoring_mode(
    outdir = outdir, region = None, volcano = volcano,
    alignsar_dc = Path("./example_data/022D_ALIGNSAR_v4.nc"),
    licsalert_settings = licsalert_settings, 
    icasar_settings = icasar_settings,
    licsbas_settings = licsbas_settings,
    )



#%% LiCSBAS COMET Volcano Portal Jasmin file

# LiCSAlert_monitoring_mode(
#     outdir = Path('005_COMET_portal_example'), 
#     region = 'africa', 
#     volcano = "erta_ale_014A_07688_131313",  
#     licsbas_jasmin_dir = Path("example_data/jasmin_simulation"),
#     licsalert_pkg_dir = Path('./'),
#     )
