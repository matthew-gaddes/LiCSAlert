#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 14:49:10 2023

@author: matthew

Thoughts:

    
    plot ICS and cumulative tcs

    add the dem
    
"""

from numpy import pi, sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

from pathlib import Path
import pickle
import sys
from glob import glob
import pdb


debug_scripts = "/home/matthew/university_work/python_stuff/python_scripts"

if debug_scripts not in sys.path:                                                                             # check if already on path
    sys.path.append(debug_scripts)
from small_plot_functions import matrix_show, quick_linegraph

#%%

print("Running")

#%%

#def licsalert_ic_combiner(licsalert_dir):
    
    
    
licsalert_dir = Path("./../001_campi_flegrei_example")
    
    # open the licsbas data.  
    
sys.path.append("./..")
    
import licsalert
from licsalert.aux import col_to_ma
    


    
def reconstruct_ts(ics_one_hot, sources_tcs, aux_data, displacement_r2):
    """ Reconstruct a LiCSAlert cumulative time series form the ICs and cumulative time course using a choice of components.
    
    Inputs:
        ics_one_hot | list | 1 if IC to be used, 0 if not.  Must be same length as number of ICS  e.g. [1,0,0,0]
        sources_tcs | list of dicts | one item in list for each IC.  Contains cumulative time course and associated data.  
        aux_data | dict | dict_keys(['icasar_sources', 'dem', 'mask'])
        displacement_r2 | dict | dict_keys(['dem', 'mask', 'lons', 'lats', 'E', 'N', 'U', 'incremental', 'means', 'incremental_downsampled', 'mask_downsampled'])
    Returns:
        X | r2 array | cumulative ifgs as rows (must be combined with a mask to turn back to masked arrays).  Mean centering has been removed.  
    History:
        2023_10_26 | MEG | Written
        
    """
    if len(ics_one_hot) != len(sources_tcs):
        raise Exception(f"'sources_tcs' is of length {len(sources_tcs)} so contains {len(sources_tcs)} sources.  "
                        f"However, 'ics_one_hot' is {len(ics_one_hot)}, which doesn't agree.  Exiting.  ")
        
    n_sources = len(sources_tcs)
    n_times = sources_tcs[0]['cumulative_tc'].shape[0]
    n_pixels = displacement_r2['incremental'].shape[1]
    
    A = np.zeros((n_times, n_sources))
    for n_source in range(n_sources):
        A[:, n_source] = np.ravel(ics_one_hot[n_source] *  sources_tcs[n_source]['cumulative_tc'])
    S = aux_data['icasar_sources']
    
    if displacement_r2['means'].shape[0] == n_times:                                   # mean for each pixel means sICA was run
        means_r2 = np.repeat(displacement_r2['means'][:,np.newaxis], n_pixels, axis = 1 )
    
    elif displacement_r2['means'].shape[0] == n_pixels:                                   # mean for each pixel means tICA was run
        means_r2 = np.repeat(displacement_r2['means'][np.newaxis, :], n_times, axis = 0)
    
    X = A@S + means_r2
    
    return X

def plot_original_reconstruction_dem(fig, ax_orig, ax_reco, ax_dem, cax_def, cax_dem, cumulative_r2, cumulative_reco_r2, mask, pixel):
    """ Make the two plots that show the raw signal (from the input time series) and the reconstruction using various 
    IC components.  
    
    Inputs:
        ax_orig | matplotlib axes | axes to plot the orignal data on
        ax_reco | matplotlib axes | axes to plot the reconstruction data on
        ax_dem  | matplotlib axes | axes to plot the DEM on
        cax_dem | matplotlib axes | axes  for the horizontal colorbar for DEM
        cax_def | matplotlib axes | axes  for the horizontal colorbar for deformaiton (original and reconstruction)
        
        cumulative_r2 | r2 array | cumulative ifgs as row vectors.  
        cumulative_reco_r2 | r2 array | cumulative ifgs as row vectors, reconstructed using the ICs
        mask | r2 boolean | True where data is masked out.  
        pixel | dict | contains 'x' and 'y', integer values of pixel being plotted in the time series.  
    Returns:
        Figure
    History:
        2023_11_01 | MEG | Written.  
    
    """
    for ax in [ax_orig, ax_reco, ax_dem]:
        ax.clear()
        
    
    vmin = np.min(np.concatenate([cumulative_r2[-1,], cumulative_reco_r2[-1]]))
    vmax = np.max(np.concatenate([cumulative_r2[-1,], cumulative_reco_r2[-1]]))
    
    reconstruction = ax_reco.matshow(col_to_ma(cumulative_reco_r2[-1,], mask), vmin = vmin, vmax = vmax)             # Plot the reconstructed last cumulative ifg.  
    ax_reco.scatter(pixel['x'], pixel['y'], s = 20, c = 'r', marker = 'x')
    ax_reco.set_ylabel('Reconstruction')        
        
    original = ax_orig.matshow(col_to_ma(cumulative_r2[-1,], mask), vmin = vmin, vmax = vmax)                         # Plot the raw data last cumulative ifg.  
    ax_orig.scatter(pixel['x'], pixel['y'], s = 20, c = 'r', marker = 'x')
    ax_orig.set_ylabel('Cumulative\nInterferogram')
    
    for ax in [ax_reco, ax_orig]:
        ax.set_xticks([])
        ax.set_yticks([])
        fig.add_subplot(ax)
        
    # : Plot the colorbar
    ics_cbar = f.colorbar(original, cax=cax_def, orientation='horizontal')
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
    f.add_subplot(ax_dem)
    
    f.colorbar(dem_plot, cax = cax_dem, orientation = 'horizontal')                                    # colorbar, tick only 0 and the max (and check max is not a nan)
    cax_dem.tick_params(axis='both', which='major', labelsize=8, rotation = 315)                                               #
    cax_dem.set_xlabel('DEM (m)')            

    
    #ics_cbar.ax.yaxis.set_label_position('left')
    
    



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
    
    from licsalert.aux import r2_to_r3
    
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
    
    
    #%%
        
from matplotlib import ticker
from matplotlib.widgets import CheckButtons
import matplotlib.gridspec as gridspec

from licsalert.data_importing import open_aux_data, open_tcs
from licsalert.plotting import truncate_colormap, xticks_every_nmonths
from licsalert.aux import r2_to_r3, moving_average
  
displacement_r2, tbaseline_info, aux_data = open_aux_data(licsalert_dir)
sources_tcs = open_tcs(licsalert_dir)    

n_sources = len(sources_tcs)
n_pixels = np.size(displacement_r2['incremental'], axis = 1)


cumulative_r2 = np.concatenate((np.zeros((1, n_pixels)), np.cumsum(displacement_r2['incremental'], axis = 0)), axis = 0)                                                 # calculate the cumulative displacments, 0 on first acquisition
cumulative_reco_r2 = np.concatenate((np.zeros((1, n_pixels)), reconstruct_ts([1 for i in range(n_sources)], sources_tcs, aux_data, displacement_r2)), axis = 0)          # reconstruct the data using all the sources, 0 on first acquisition

cumulative_r3 = r2_to_r3(cumulative_r2, displacement_r2['mask'])                                                                                            # conver to rank 3
t, y, x = np.unravel_index(np.argmax(cumulative_r3), cumulative_r3.shape)                                                                                   # which makes finding the 
pixel = {'x' : x, 'y' : y}
del t, y, x, cumulative_r3


fig_width = 18
f = plt.figure(figsize = (fig_width, fig_width /  (2 * (1920 / 1080))))
f.canvas.manager.set_window_title('LiCSAlert results visualiser')

if n_sources % 2 == 1:
    n_rows = n_sources + 1
else:
    n_rows = n_sources
grid = gridspec.GridSpec(n_rows, 20, wspace=0.2, hspace=0.2)                        # divide into 2 sections, 1/5 for ifgs and 4/5 for components

# create all the axes from the grid
ax_cum = plt.Subplot(f, grid[:int(n_sources/2), 10:15])                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
ax_reco = plt.Subplot(f, grid[int(n_sources/2) : , 10:15])                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
ax_ts = plt.subplot(grid[:5,15:])                                                                                      # create an axes for the IC (spatial source)
ax_dem = plt.Subplot(f, grid[int(n_sources/2) : , 5:10])                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
cax_def = f.add_axes([0.535, 0.1, 0.15, 0.02])                                                              # 
cax_dem = f.add_axes([0.34, 0.1, 0.15, 0.02])                                                              # 
cax_ics = f.add_axes([0.125, 0.11, 0.03, 0.02])                                                              # 
rax = plt.axes([0.35, 0.6, 0.1, 0.2])                                                                       # radio buttons


plot_original_reconstruction_dem(f, ax_cum, ax_reco, ax_dem, cax_def, cax_dem, cumulative_r2, cumulative_reco_r2, displacement_r2['mask'], pixel)


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
    #ax_ts.set_ylabel("LOS displacemnt (m)")
    
    if source_n == (n_sources - 1):
        include_tick_labels = True
    else:
        include_tick_labels = False
    xticks_every_nmonths(ax_ctc, tbaseline_info['acq_dates'][0], tbaseline_info['baselines_cumulative'], include_tick_labels = include_tick_labels, 
                         major_ticks_n_months = 12, minor_ticks_n_months = 1)
    
## ICs colorbar
f.colorbar(ic, cax = cax_ics, orientation = 'horizontal')                                    # colorbar, tick only 0 and the max (and check max is not a nan)
cax_ics.tick_params(axis='both', which='major', labelsize=8, rotation = 315)                                               #
cax_ics.set_xlabel('IC')            


# start the interative bit
button_names = [f"IC{i}" for i in range(n_sources)]
button_status = [True for i in range(n_sources)]
check = CheckButtons(rax, button_names, button_status)




plot_ts(f, ax_ts, cumulative_r2, cumulative_reco_r2, displacement_r2['mask'], tbaseline_info, pixel)









def replot_on_pixel_select(event):
    """ When the selected pixel changes, replot the time series for that point.  
    """
    
    button=event.button

    if (event.inaxes is ax_reco) or (event.inaxes is ax_cum) or (event.inaxes is ax_dem):               # check we are in axes as otherwise there's no meaning to click position in data coords.  
        pixel['x'] = int(event.xdata)
        pixel['y'] = int(event.ydata)
        
        ax_ts.clear()
        plot_ts(f, ax_ts, cumulative_r2, cumulative_reco_r2, 
                displacement_r2['mask'], tbaseline_info,  pixel)
    
        plot_original_reconstruction_dem(f, ax_cum, ax_reco, ax_dem, cax_def, cax_dem, cumulative_r2, cumulative_reco_r2, displacement_r2['mask'], pixel)
        
        plt.draw()
        
        
    
    



def replot_on_ic_select(label):
    """ When the ICs selected changes, replot the reconstruction and the time series using those ICs.  
    """
    ics_one_hot = [1 if status else 0 for status in check.get_status()]                                                 # convert boolean to ints and work out which buttons are currently selected
    cumulative_reco_r2 = reconstruct_ts(ics_one_hot, sources_tcs, aux_data, displacement_r2)                                # make the time series with those sources.
    cumulative_reco_r2 = np.concatenate((np.zeros((1, n_pixels)), cumulative_reco_r2))                                     # add zero to first acquisition
    plot_original_reconstruction_dem(f, ax_cum, ax_reco, ax_dem, cax_def, cax_dem, cumulative_r2, cumulative_reco_r2,
                                     displacement_r2['mask'], pixel)
    ax_ts.clear()                                                                                                           # clear the time series plot ready for new plot
    plot_ts(f, ax_ts, cumulative_r2, cumulative_reco_r2,    
            displacement_r2['mask'], tbaseline_info, pixel)                                                                 # replot the time series using the new reconstruction.  
    plt.draw()

    
        
check.on_clicked(replot_on_ic_select)   
cid = f.canvas.mpl_connect('button_press_event', replot_on_pixel_select)
    


    


    
#%%
    