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
import pdb

import licsalert
from licsalert.monitoring_functions import LiCSAlert_monitoring_mode
#from licsalert.plotting import licsalert_results_explorer




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


# licsalert_settings = {"baseline_end" : "20170101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
#                       "figure_intermediate" : False,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
#                       "figure_type"         : 'both',                             # either 'window' or 'png' (to save as pngs), or 'both'
#                       "downsample_run"      : 0.5,                                     # data can be downsampled to speed things up
#                       "downsample_plot"      : 0.5,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
#                       "residual_type"        : 'cumulative'}                      # controls the type of residual used in the lower plot.  Either cumulative or window   


# icasar_settings = {"n_comp" : 5,                                                  # number of components to recover with ICA (ie the number of PCA sources to keep)
#                    "bootstrapping_param" : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
#                    "tsne_param" : (30, 12),                                       # (perplexity, early_exaggeration)
#                    "ica_param" : (1e-2, 150),                                     # (tolerance, max iterations)
#                    "hdbscan_param" : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
#                    "ifgs_format"        : 'cum',                                  # can be 'all', 'inc' (incremental - short temporal baselines), or 'cum' (cumulative - relative to first acquisition)
#                    "sica_tica"          : 'sica' }



# licsalert_dir = Path("./")                                                        # path to licsalert package
# volcano = '001_campi_flegrei_example'                                                 # outdir final                                                                                                                      
# licsbas_dir = Path("./022D_04826_121209")                                         # input data


# LiCSAlert_monitoring_mode(region = None, volcano = volcano, LiCSAlert_pkg_dir = LiCSAlert_pkg_dir,
#                           licsbas_dir = licsbas_dir, licsalert_dir = licsalert_dir,
#                           licsalert_settings = licsalert_settings, icasar_settings = icasar_settings)


# licsalert_results_explorer(licsalert_dir / volcano, fig_width = 18)                                                 # use this function to explore the results



    
#%% Example 2: temporal ICA at Vesuvius

licsalert_dir = Path("./")                                                        # path to licsalert package
volcano = '002_vesuvius_example'                                                 # outdir final                                                                                                                      
licsbas_dir = Path("./022D_04826_121209_vesuvius_crop_rationalized")                                         # input data


licsalert_settings = {"baseline_end" : "20170101",                               # end baseline stage at YYYYMMDD, need to be before the last acquisition of LiCSAlert will never monitoring anyhting.  
                      "figure_intermediate" : False,                             # if set to True, a figure is produced for all time steps in the monitoring data, which can be time consuming.  
                      "figure_type"         : 'png',                             # either 'window' or 'png' (to save as pngs), or 'both'
                      "downsample_run"      : 0.5,                                     # data can be downsampled to speed things up
                      "downsample_plot"      : 0.5,                               # and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
                      "residual_type"        : 'cumulative'}                      # controls the type of residual used in the lower plot.  Either cumulative or window   

icasar_settings = {"n_comp" : 5,                                                  # number of components to recover with ICA (ie the number of PCA sources to keep)
                   "bootstrapping_param" : (200, 0),                              # (number of runs with bootstrapping, number of runs without bootstrapping)                    "hdbscan_param" : (35, 10),                        # (min_cluster_size, min_samples)
                   "tsne_param" : (30, 12),                                       # (perplexity, early_exaggeration)
                   "ica_param" : (1e-2, 150),                                     # (tolerance, max iterations)
                   "hdbscan_param" : (100,10),                                    # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise. 
                   #"ifgs_format"        : 'cum',                                  # ifg format is redundant when using tICA as cumulative signals must be considered when working temporally.  
                   "sica_tica"          : 'tica' }

volcano = '002_vesuvius_example_tica'                                                         # new outdir


LiCSAlert_monitoring_mode(region = None, volcano = volcano, LiCSAlert_pkg_dir = LiCSAlert_pkg_dir,
                          licsbas_dir = licsbas_dir, licsalert_dir = licsalert_dir,
                          licsalert_settings = licsalert_settings, icasar_settings = icasar_settings)



def licsalert_results_explorer(licsalert_out_dir, fig_width = 18):
    """ The interactive figure for exploring LiCSAlert results.
    
    Inputs:
        licsalert_out_dir | pathlib Path | path to directory of LiCSAlert results.   e.g. Path("./../001_campi_flegrei_example")
        fig_width | int | figure width in inches.  Height is set relative to this and optimised for display on two FHD screens (i.e. 2 x 1920 x 1080)
    
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from matplotlib.widgets import CheckButtons
    import matplotlib.gridspec as gridspec
    
    from licsalert.licsalert import reconstruct_ts
    from licsalert.data_importing import open_aux_data, open_tcs
    from licsalert.plotting import truncate_colormap, xticks_every_nmonths
    from licsalert.aux import r2_to_r3, moving_average, col_to_ma
    
    
    def plot_original_reconstruction_dem(fig, ax_orig, ax_reco, ax_dem, cax_def, cax_dem, cumulative_r2, cumulative_reco_r2, displacement_r2, pixel):
        """ Make the three plots that show the raw signal (from the input time series) and the reconstruction using various 
        IC components, and the DEM.  
        
        Inputs:
            ax_orig | matplotlib axes | axes to plot the orignal data on
            ax_reco | matplotlib axes | axes to plot the reconstruction data on
            ax_dem  | matplotlib axes | axes to plot the DEM on
            cax_dem | matplotlib axes | axes  for the horizontal colorbar for DEM
            cax_def | matplotlib axes | axes  for the horizontal colorbar for deformaiton (original and reconstruction)
            
            cumulative_r2 | r2 array | cumulative ifgs as row vectors.  
            cumulative_reco_r2 | r2 array | cumulative ifgs as row vectors, reconstructed using the ICs
            displacement_r2 | dict | must contain the mask, DEM, lons and lats (all r2 arrays)
            pixel | dict | contains 'x' and 'y', integer values of pixel being plotted in the time series.  
        Returns:
            Figure
        History:
            2023_11_01 | MEG | Written.  
        """
        #from matplotlib import ticker
        from licsalert.plotting import truncate_colormap
        
        for ax in [ax_orig, ax_reco, ax_dem]:
            ax.clear()
            
        vmin = np.min(np.concatenate([cumulative_r2[-1,], cumulative_reco_r2[-1]]))
        vmax = np.max(np.concatenate([cumulative_r2[-1,], cumulative_reco_r2[-1]]))
        
        reconstruction = ax_reco.matshow(col_to_ma(cumulative_reco_r2[-1,], displacement_r2['mask']), vmin = vmin, vmax = vmax)             # Plot the reconstructed last cumulative ifg.  
        ax_reco.scatter(pixel['x'], pixel['y'], s = 20, c = 'r', marker = 'x')
        ax_reco.set_ylabel('Reconstruction')        
            
        original = ax_orig.matshow(col_to_ma(cumulative_r2[-1,], displacement_r2['mask']), vmin = vmin, vmax = vmax)                         # Plot the raw data last cumulative ifg.  
        ax_orig.scatter(pixel['x'], pixel['y'], s = 20, c = 'r', marker = 'x')
        ax_orig.set_ylabel('Cumulative\nInterferogram')
        
        for ax in [ax_reco, ax_orig]:
            ax.set_xticks([])
            ax.set_yticks([])
            fig.add_subplot(ax)
            
        # : Plot the colorbar
        ics_cbar = fig.colorbar(original, cax=cax_def, orientation='horizontal')
        tick_locator = ticker.MaxNLocator(nbins=4)
        ics_cbar.locator = tick_locator
        ics_cbar.update_ticks()
        ics_cbar.set_label('LOS displacement (m)')
        
        # plot the DEM
        terrain_cmap = plt.get_cmap('terrain')                                                                                  # appropriate colours for a dem
        terrain_cmap = truncate_colormap(terrain_cmap, 0.2, 1)                                                                  # but crop (truncate) the blue parts as we are only interested in land
        dem_plot = ax_dem.matshow(displacement_r2["dem"], cmap = terrain_cmap)                                                   # plot the DEM
        ax_dem.tick_params(axis='both', which='major', labelsize=7)                                                             # adjust fontsize
        ax_dem.tick_params(axis='both', which='minor', labelsize=7)
        ax_dem.set_xticks([0, displacement_r2['dem'].shape[1]])                                                                 # tick only the min and max in each direction
        ax_dem.set_yticks([0, displacement_r2['dem'].shape[0]])
        ax_dem.set_xticklabels([str(round(displacement_r2['lons'][-1,0], 2)) + "$^\circ$", str(round(displacement_r2['lons'][-1,-1], 2)) + "$^\circ$" ])
        ax_dem.set_yticklabels([str(round(displacement_r2['lats'][0,0], 2))  + "$^\circ$", str(round(displacement_r2['lats'][-1,0], 2)) + "$^\circ$"])
        ax_dem.yaxis.tick_right()
        ax_dem.scatter(pixel['x'], pixel['y'], s = 20, c = 'r', marker = 'x')
        fig.add_subplot(ax_dem)
        
        fig.colorbar(dem_plot, cax = cax_dem, orientation = 'horizontal')                                    # colorbar, tick only 0 and the max (and check max is not a nan)
        cax_dem.tick_params(axis='both', which='major', labelsize=8, rotation = 315)                                               #
        cax_dem.set_xlabel('DEM (m)')            
    
        
    
    
    def plot_ts(fig, ax_ts, cumulative_r2, cumulative_reco_r2, mask, tbaseline_info, pixel):
        """
        Plot the time series for a pixel in two different datasets.  
        
        Inputs:
            fig | matplotlib figure | figure axes is in.  
            ax_ts | matplotlib axes | axes to draw in.  
            cumulative_r2 | r2 array | cumulative ifgs as row vectors.  
            cumulative_reco_r2 | r2 array | reconstrcution of cumultaive data.  
            mask | r2 boolean | true where maskd.  
            tbaseline_info | dict | must contain baselines_cumulative (i.e. 0,6,12,18 for 6 day acquisitions)
            pixel | dict | if None, pixel with maximum deformation is plotted.  Or can be dict contiing 'x' and 'y' to choose pixel.  
        Returns:
            Plot in axes
        History:
            2023 | MEG | Written
    
        """
        
        from licsalert.aux import r2_to_r3, moving_average
        from licsalert.plotting import xticks_every_nmonths
        
        
        data_colours = ['tab:purple', 'tab:orange']
        data_names = ['Original', 'Reconstruction']
        
        cumulative_r3 = r2_to_r3(cumulative_r2, mask)
        cumulative_reco_r3 = r2_to_r3(cumulative_reco_r2, mask)
        
        
        for i, data in enumerate([cumulative_r3, cumulative_reco_r3]):
    
            ts = data[:, pixel['y'], pixel['y']]                                                                                   # get the time series for the pixel of interest
            ts_smooth, valid = moving_average(ts)                                                           # smooth it 
            ax_ts.scatter(np.concatenate((np.array([0]), tbaseline_info['baselines_cumulative'])),
                          ts, alpha = 0.4, marker = '.', s = 4, 
                          label = data_names[i], c = data_colours[i])                                       # plot each point for the 
            ax_ts.plot(np.concatenate((np.array([0]), tbaseline_info['baselines_cumulative'])),
                       ts_smooth, c = data_colours[i])                                   # and the smoothed one as a line
        
        
        ax_ts.axhline(0, c = 'k')
        ax_ts.grid(True)
        ax_ts.yaxis.tick_right()
        fig.add_subplot(ax_ts)                                                                   # add to figure
        ax_ts.yaxis.set_label_position("right")
        ax_ts.set_ylabel("LOS displacemnt (m)")
        
        xticks_every_nmonths(ax_ts, tbaseline_info['acq_dates'][0], tbaseline_info['baselines_cumulative'], include_tick_labels = True, 
                             major_ticks_n_months = 12, minor_ticks_n_months = 1)
    
        ax_ts.legend()

    def replot_on_pixel_select(event):
        """ When the selected pixel changes, replot the time series for that point and the three images with a point showing where was clicked.  
        """
        
        button=event.button
        if (event.inaxes is ax_reco) or (event.inaxes is ax_cum) or (event.inaxes is ax_dem):                          # check we are in axes as otherwise there's no meaning to click position in data coords.  
            pixel['x'] = int(event.xdata)                                                                               # get the data coords of where we clicked and update hte dict.  
            pixel['y'] = int(event.ydata)
            ax_ts.clear()                                                                                               # remove previous plot
            plot_ts(f, ax_ts, cumulative_r2, cumulative_reco_r2,        
                    displacement_r2['mask'], tbaseline_info,  pixel)                                                    # plot the time series for that pixel in original and reconstruction.  
        
            plot_original_reconstruction_dem(f, ax_cum, ax_reco, ax_dem, cax_def, cax_dem, 
                                             cumulative_r2, cumulative_reco_r2, displacement_r2, pixel)         # also replot the images with a cross on the selected point.  
            plt.draw()
    
    def replot_on_ic_select(label):
        """ When the ICs selected changes, replot the reconstruction and the time series using those ICs.  
        """
        ics_one_hot = [1 if status else 0 for status in check.get_status()]                                                    # convert boolean to ints and work out which buttons are currently selected
        cumulative_reco_r2 = reconstruct_ts(ics_one_hot, sources_tcs, aux_data, displacement_r2)                               # make the time series with those sources.
        cumulative_reco_r2 = np.concatenate((np.zeros((1, n_pixels)), cumulative_reco_r2))                                     # add zero to first acquisition
        plot_original_reconstruction_dem(f, ax_cum, ax_reco, ax_dem, cax_def, cax_dem, 
                                         cumulative_r2, cumulative_reco_r2, displacement_r2, pixel)                    # plot the images with the new reconsrution.   
        ax_ts.clear()                                                                                                          # clear the time series plot ready for new plot
        plot_ts(f, ax_ts, cumulative_r2, cumulative_reco_r2,    
                displacement_r2['mask'], tbaseline_info, pixel)                                                                 # replot the time series using the new reconstruction.  
        plt.draw()
    

      
    # 1: Open data and some simple processing 
    displacement_r2, tbaseline_info, aux_data = open_aux_data(licsalert_out_dir)
    sources_tcs = open_tcs(licsalert_out_dir)    
    
    n_sources = len(sources_tcs)
    n_pixels = np.size(displacement_r2['incremental'], axis = 1)
    
    cumulative_r2 = np.concatenate((np.zeros((1, n_pixels)), np.cumsum(displacement_r2['incremental'], axis = 0)), axis = 0)                                                 # calculate the cumulative displacments, 0 on first acquisition
    cumulative_reco_r2 = np.concatenate((np.zeros((1, n_pixels)), reconstruct_ts([1 for i in range(n_sources)], sources_tcs, aux_data, displacement_r2)), axis = 0)          # reconstruct the data using all the sources, 0 on first acquisition
    cumulative_r3 = r2_to_r3(cumulative_r2, displacement_r2['mask'])                                                                                                         # conver to rank 3
    
    
    #%% debugging
    # pdb.set_trace()
    # plt.switch_backend('qt5agg')
    
    # matrix_show(col_to_ma(np.sum(displacement_r2['incremental'][:,], axis = 0), displacement_r2['mask']))          # noise and zeros
    # matrix_show(col_to_ma(displacement_r2['incremental'][-1,], displacement_r2['mask']))                    # looks like a normal ifg
    
    # matrix_show(col_to_ma(cumulative_r2[-1], displacement_r2['mask']))          # noise and zeros
    # matrix_show(cumulative_r3[-1,])                                             # noise and zeros
    
    # cumulative_reco_r3 = r2_to_r3(cumulative_reco_r2, displacement_r2['mask'])                                                                                            # conver to rank 3

    # matrix_show(cumulative_reco_r3[-1,])                    # signal, but very small?
    
    
    #%%
    
    t, y, x = np.unravel_index(np.argmax(cumulative_r3), cumulative_r3.shape)                                                                                   # which makes finding the 
    pixel = {'x' : x, 'y' : y}
    del t, y, x, cumulative_r3
    
    
    

    # 2: start the figure.      
    f = plt.figure(figsize = (fig_width, fig_width /  (2 * (1920 / 1080))))
    f.canvas.manager.set_window_title('LiCSAlert results visualiser')
    
    if n_sources % 2 == 1:                                                                                  # number of rows is the smallest even number greater than or equal to the number of sources.  
        n_rows = n_sources + 1
    else:
        n_rows = n_sources
    grid = gridspec.GridSpec(n_rows, 20, wspace=0.2, hspace=0.2)                                            # 
    
    # create all the axes from the grid
    ax_cum = plt.Subplot(f, grid[:int(n_sources/2), 10:15])                                                # for the final cumulative ifg.  
    ax_reco = plt.Subplot(f, grid[int(n_sources/2) : , 10:15])                                             # for the reconstruction of the cumulative ifg.  
    ax_ts = plt.subplot(grid[:5,15:])                                                                      # for the time series of a point in both original and reconstruted.  
    ax_dem = plt.Subplot(f, grid[int(n_sources/2) : , 5:10])                                               # for the dem 
    cax_def = f.add_axes([0.535, 0.1, 0.15, 0.02])                                                         # cumulative and reconstruction shared colorbar
    cax_dem = f.add_axes([0.34, 0.1, 0.15, 0.02])                                                          # dem colorbar 
    cax_ics = f.add_axes([0.125, 0.11, 0.03, 0.02])                                                        # ICs colorbar
    rax = plt.axes([0.35, 0.6, 0.1, 0.2])                                                                  # radio buttons
    
    
    # 3: start plotting.  
    plot_original_reconstruction_dem(f, ax_cum, ax_reco, ax_dem, cax_def, cax_dem, 
                                     cumulative_r2, cumulative_reco_r2, displacement_r2, pixel)     # plot the cumulative, reconscructcion, and DEM.   
    
    plot_ts(f, ax_ts, cumulative_r2, cumulative_reco_r2, displacement_r2['mask'], tbaseline_info, pixel)    # plot the time series for the point of interest in the two ways.  
    
    ## Plot the ICs 
    vmin = np.min(aux_data['icasar_sources'])
    vmax = np.max(aux_data['icasar_sources'])
    for source_n, source in enumerate(aux_data['icasar_sources']):
        # plot spatial pattern
        ax_ic = plt.Subplot(f, grid[source_n, 0])                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
        ic = ax_ic.matshow(col_to_ma(source, displacement_r2['mask']), vmin = vmin, vmax = vmax)                         # Plot the raw data last cumulative ifg.  
        ax_ic.set_title(f'IC {source_n}')
        ax_ic.set_xticks([])
        ax_ic.set_yticks([])
        f.add_subplot(ax_ic)
        
        # plot the cumulative tc
        ax_ctc = plt.Subplot(f, grid[source_n, 1:5])                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
        ctc = np.concatenate((np.array([[0]]), sources_tcs[source_n]['cumulative_tc']), axis = 0)
        ctc_smooth, valid = moving_average(ctc)                                                           # smooth it 
        ax_ctc.scatter(np.concatenate((np.array([0]), tbaseline_info['baselines_cumulative'])), ctc, alpha = 0.4, marker = '.', s = 2) 
        ax_ctc.plot(np.concatenate((np.array([0]), tbaseline_info['baselines_cumulative'])), ctc_smooth)                                   # and the smoothed one as a line
        ax_ctc.axhline(0, c = 'k')
        ax_ctc.grid(True)
        ax_ctc.yaxis.tick_right()
        ax_ctc.tick_params(axis='both', which='both', labelsize=8)
        f.add_subplot(ax_ctc)
        
        if source_n == (n_sources - 1):
            include_tick_labels = True
        else:
            include_tick_labels = False
        xticks_every_nmonths(ax_ctc, tbaseline_info['acq_dates'][0], tbaseline_info['baselines_cumulative'], include_tick_labels = include_tick_labels, 
                             major_ticks_n_months = 12, minor_ticks_n_months = 1)
        
    f.colorbar(ic, cax = cax_ics, orientation = 'horizontal')                                    # ICs colorbar, tick only 0 and the max (and check max is not a nan)
    cax_ics.tick_params(axis='both', which='major', labelsize=8, rotation = 315)                                               #
    cax_ics.set_xlabel('IC')            
    
    # 3:  start the interative bit
    button_names = [f"IC{i}" for i in range(n_sources)]
    button_status = [True for i in range(n_sources)]
    check = CheckButtons(rax, button_names, button_status)

    check.on_clicked(replot_on_ic_select)                                                       # check buttons to select which IC
    cid = f.canvas.mpl_connect('button_press_event', replot_on_pixel_select)                    # click on point to plot it.  
    





licsalert_results_explorer(licsalert_dir / volcano, fig_width = 18)                                                 # use this function to explore the results


sys.exit()

#%% We can also make the LiCSAlert figure for all times, but this is slow

# licsalert_settings['figure_intermediate'] = True                                                    # Turn on figures for all times 
# licsalert_settings['figure_type'] = 'png'                                                           # There will be too many figure windows if this is 'window' or 'both' ! 

# volcano = 'campi_flegrei_example_all_times'                                                         # new outdir

# LiCSAlert_monitoring_mode(region = None, volcano = volcano, LiCSAlert_pkg_dir = LiCSAlert_pkg_dir,
#                           licsbas_dir = licsbas_dir, licsalert_dir = licsalert_dir,
#                           licsalert_settings = licsalert_settings, icasar_settings = icasar_settings)
