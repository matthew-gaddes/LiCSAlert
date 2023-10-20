#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 11:41:16 2023

@author: matthew
"""


#%%

def two_spatial_signals_plot(images, mask, dem, tcs_dc, tcs_all, t_baselines_dc, t_baselines_all, 
                             title, ifg_dates_dc, fig_kwargs):
    """
    Product the two plots that show spatial sources (comparison to DEM and ifg baseline, and then sources and cumulative time courses).  
    Inputs:
        images | n_images x n_pixels | spatial signals as row vectors.  
        mask | boolean rank 2 |         mask to convert row vector to a masked array.  
        dem | rank 2 masked array |         The DEM.  Used to compare spatial signals to it.  
        tcs_dcs | n_times x n_images      | daisy chain time courses.  
        tcs_all | n_all_ifgs x n_images or None  | time courses if we have made all possible ifg pairs.  None if not available.  
        t_baselines_dc | n_times         | temporal baselines for the diasy chain interferograms.  
        t_baselines_all | n_all_ifgs     | temporal baselines for all possible ifg pairs, None if not available.  
        title | string | figure title
        fig_kwargs | dict | see other functions.  
    Returns:
        2 figures
        dem_to_sources_comparisons | dict | results of comparison between spatial sources and DEM.  xyzs: x and y and colours for each point from kernel density estimate, line_xys | list | xy points for line of best fit. corr_coefs | pearson correlation coefficeint for points.  
        tcs_to_tempbaselines_comparisons | dict | as above, but correlations between temporal baselines and time courses (ie if long temrpoal baseline, is component used strongly.  )
       
    History:
        2021_11_23 | MEG | Written
        2022_01_14 | MEG | Add return of comparison dicts.  
        
    """
    
    # 1: First figure is simpler as only use daisy chain time courses
    plot_spatial_signals(images.T, mask, tcs_dc.T, mask.shape, title = f"{title}_time",
                         temporal_baselines = t_baselines_dc, ifg_dates_dc = ifg_dates_dc, **fig_kwargs)                      # the usual plot of the sources and their time courses (ie contributions to each ifg)                              

    # 2: Second figure may have access to all interfergram time courses and temporal baselines, but may also not.          
    if t_baselines_all is not None:
        temporal_data = {'tcs'                : tcs_all,
                         'temporal_baselines' : t_baselines_all}
    else:
        temporal_data = {'tcs'                : tcs_dc,
                         'temporal_baselines' : t_baselines_dc}
        
    dem_to_sources_comparisons, tcs_to_tempbaselines_comparisons = dem_and_temporal_source_figure(images, mask, fig_kwargs, dem, temporal_data, fig_title = f"{title}_correlations")        # also compare the sources to the DEM, and the correlation between their time courses and the temporal baseline of each interferogram.                                                                                                              # 
                                                                                                                                                                                            # also note that it now returns information abou the sources and correlatiosn (comparison to the DEM, and how they're used in time.  )
    return dem_to_sources_comparisons, tcs_to_tempbaselines_comparisons


  
  
#%%
def plot_spatial_signals(spatial_map, pixel_mask, tcs, shape, title, ifg_dates_dc, 
                         temporal_baselines, figures = "window",  png_path = './'):
    """
    Input:
        spatial map | pxc matrix of c component maps (p pixels) (i.e. images are column vectors)
        pixel_mask | mask to turn spaital maps back to regular grided masked arrays
        tcs | cxt matrix of c time courses (t long)   
        shape | tuple | the shape of the grid that the spatial maps are reshaped to
        title | string | figure tite and png filename (nb .png will be added, don't include here)
        temporal_baselines | x axis values for time courses.  Useful if some data are missing (ie the odd 24 day ifgs in a time series of mainly 12 day)
        figures | string,  "window" / "png" / "png+window" | controls if figures are produced (either as a window, saved as a png, or both)
        png_path | string | if a png is to be saved, a path to a folder can be supplied, or left as default to write to current directory.  
        
    Returns:
        Figure, either as a window or saved as a png
        
    2017/02/17 | modified to use masked arrays that are given as vectors by spatial map, but can be converted back to 
                 masked arrays using the pixel mask    
    2017/05/12 | shared scales as decribed in 'shared'
    2017/05/15 | remove shared colorbar for spatial maps
    2017/10/16 | remove limit on the number of componets to plot (was 5)
    2017/12/06 | Add a colorbar if the plots are shared, add an option for the time courses to be done in days
    2017/12/?? | add the option to pass temporal baselines to the function
    2020/03/03 | MEG | Add option to save figure as png and close window
    """
    
    import numpy as np
    import numpy.ma as ma  
    import matplotlib.pyplot as plt
    import matplotlib
    import matplotlib.gridspec as gridspec
    import pdb
    

    
    def linegraph(sig, ax, temporal_baselines = None):
        """ signal is a 1xt row vector """
        
        if temporal_baselines is None:
            times = sig.size
            a = np.arange(times)
        else:
            a = np.cumsum(temporal_baselines)                                        # cumulative baselines from temporal baselines.  
        ax.scatter(a,sig,marker='.', color='k', alpha = 0.5)
        ax.plot(a,sig, linestyle = '-', color='k',)
        ax.axhline(y=0, color='k', alpha=0.4) 
        
    def xticks_every_3months(ax_to_update, day0_date, time_values, include_tick_labels):
        """Given an axes, update the xticks so the major ones are the 1st of jan/april/july/october, and the minor ones are the first of the 
        other months.  
        Inputs:
            ax_to_update | matplotlib axes | the axes to update.  
            day0_date | string | in form yyyymmdd
            time_values | rank 1 array | cumulative temporal baselines, e.g. np.array([6,18, 30, 36, 48])
        Returns:
            updates axes
        History:
            2021_09_27 | MEG | Written
        """
        import datetime as dt
        from dateutil.relativedelta import relativedelta                                                    # add 3 months and check not after end
        from matplotlib.ticker import AutoMinorLocator      
        
       
        xtick_label_angle = 315
        
        tick_labels_days = ax_to_update.get_xticks().tolist()                                                # get the current tick labels
        day0_date_dt = dt.datetime.strptime(day0_date, "%Y%m%d")                                            
        dayend_date_dt = day0_date_dt +  dt.timedelta(int(time_values[-1]))                                 # the last time value is the number of days we have, so add this to day0 to get the end.  
    
        # 1: find first tick date (the first of the jan/ april/jul /oct)                        
        date_tick0 = day0_date_dt                                                                           
        while not ( (date_tick0.day) == 1 and (date_tick0.month == 1  or date_tick0.month == 4 or date_tick0.month == 7 or date_tick0.month == 10 )):
            date_tick0 +=  dt.timedelta(1)
            
        # 2: get all the other first of the quarters
        ticks = {'datetimes' : [date_tick0],
                 'yyyymmdd'   : [],
                 'n_day'     : []}
       
        while ticks['datetimes'][-1] < (dayend_date_dt - relativedelta(months=+3)):                         # subtract 3 months to make sure we don't go one 3 month jump too far. 
            ticks['datetimes'].append(ticks['datetimes'][-1] + relativedelta(months=+3))
        
        # 3: work out what day number each first of the quarter is.  
        for tick_dt in ticks['datetimes']:                                                                   # find the day nubmers from this.             
            ticks['yyyymmdd'].append(dt.datetime.strftime(tick_dt, "%Y/%m/%d"))
            ticks['n_day'].append((tick_dt - day0_date_dt).days)
            
        # 4: Update the figure.  
        ax_to_update.set_xticks(ticks['n_day'])                                                                   # apply major tick labels to the figure
        minor_locator = AutoMinorLocator(3)                                                                       # there are three months in each quarter, so a minor tick every month
        ax_to_update.xaxis.set_minor_locator(minor_locator)                                                       # add to figure.  
        if include_tick_labels:
            ax_to_update.set_xticklabels(ticks['yyyymmdd'], rotation = xtick_label_angle, ha = 'left')            # update tick labels, and rotate
            plt.subplots_adjust(bottom=0.15)
            ax_to_update.set_xlabel('Date')
        else:
            ax_to_update.set_xticklabels([])                                                                    # remove any tick lables if they aren't to be used.  
        
        # add vertical lines every year.  
        for major_tick_n, datetime_majortick in enumerate(ticks['datetimes']):
            if datetime_majortick.month == 1:
                ax_to_update.axvline(x = ticks['n_day'][major_tick_n], color='k', alpha=0.1, linestyle='--')           
        
 
    
    # colour map stuff
    ifg_colours = plt.get_cmap('coolwarm')
    cmap_mid = 1 - np.max(spatial_map)/(np.max(spatial_map) + abs(np.min(spatial_map)))          # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
    if cmap_mid < (1/257):                                                                       # this is a fudge so that if plot starts at 0 doesn't include the negative colorus for the smallest values
        ifg_colours_cent = remappedColorMap(ifg_colours, start=0.5, midpoint=0.5, stop=1.0, name='shiftedcmap')
    else:
        ifg_colours_cent = remappedColorMap(ifg_colours, start=0.0, midpoint=cmap_mid, stop=1.0, name='shiftedcmap')
    
    #make a list of ifgs as masked arrays (and not column vectors)
    spatial_maps_ma = []
    for i in range(np.size(spatial_map,1)):
        spatial_maps_ma.append(ma.array(np.zeros(pixel_mask.shape), mask = pixel_mask ))
        spatial_maps_ma[i].unshare_mask()
        spatial_maps_ma[i][~spatial_maps_ma[i].mask] = spatial_map[:,i].ravel()
    tmp, n_sources = spatial_map.shape
#    if n_sources > 5:
#        n_sources = 5
    del tmp
    
    day0_date = ifg_dates_dc[0][:8]
    
    #import pdb; pdb.set_trace()
    
    fig1 = plt.figure(figsize=(14,8))
    #f.suptitle(title, fontsize=14)
    grid = gridspec.GridSpec(n_sources, 6, wspace=0.3, hspace=0.3)                        # divide into 2 sections, 1/5 for ifgs and 4/5 for components
    fig1.canvas.manager.set_window_title(title)
    for i in range(n_sources):    
        # 0: define axes
        ax_source = plt.Subplot(fig1, grid[i,0])                                                                                        # spatial pattern
        ax_ctc = plt.Subplot(fig1, grid[i,1:])                                                                                          # ctc = cumulative time course
        # 1: plot the images (sources)
        im = ax_source.matshow(spatial_maps_ma[i], cmap = ifg_colours_cent, vmin = np.min(spatial_map), vmax = np.max(spatial_map))
        ax_source.set_xticks([])
        ax_source.set_yticks([])

        linegraph(np.cumsum(tcs[i,:]), ax_ctc, temporal_baselines)
        if i != (n_sources-1):
            xticks_every_3months(ax_ctc, day0_date, np.cumsum(temporal_baselines), include_tick_labels = False)                     # no tick labels
        else:
            xticks_every_3months(ax_ctc, day0_date, np.cumsum(temporal_baselines), include_tick_labels = True)                      # unles the last one.  
            
        # if shared ==1:
        #     ax_ctc.set_ylim([np.min(np.cumsum(tcs)), np.max(np.cumsum(tcs))])
        ax_ctc.yaxis.tick_right()
        fig1.add_subplot(ax_source)
        fig1.add_subplot(ax_ctc)
            
    ax_source.set_xlabel('Spatial sources')
    ax_ctc.set_xlabel('Cumulative time courses')
            
    # if shared == 1:                                                             # if the colourbar is shared between each subplot, the axes need extending to make space for it.
    #     #fig1.tight_layout(rect=[0.1, 0, 1., 1])
    cax = fig1.add_axes([0.03, 0.1, 0.01, 0.3])
    fig1.colorbar(im, cax=cax, orientation='vertical')
    #pdb.set_trace()
    
    if figures == 'window':                                                                 # possibly save the output
        pass
    elif figures == "png":
        try:
            fig1.savefig(f"{png_path}/{title}.png")
            plt.close()
        except:
            print(f"Failed to save the figure.  Trying to continue.  ")
    elif figures == 'png+window':
        try:
            fig1.savefig(f"{png_path}/{title}.png")
        except:
            print(f"Failed to save the figure.  Trying to continue.  ")
    else:
        pass
    
    

#%%

def dem_and_temporal_source_figure(sources, sources_mask, fig_kwargs, dem = None, temporal_data = None, fig_title = None,
                                   max_pixels = 1000):
    """ Given sources recovered by a blind signal separation method (e.g. PCA or ICA) compare them in space to hte DEM,
    and in time to the temporal baselines.  
    Inputs:
        sources | rank 2 array | as row vectors.  eg. 5 x999 for 5 sources.  
        sources_mask | rank 2 | boolean array with a True value for any masked pixel.  Number of False pixels should be the same as the number of columns in row_vectors
        fig_kwargs | dict | pass straing to plot_source_tc_correlations, see that for details.  
        dem | masked array | a DEM as a masked array.  It should work if not availble.  
        temporal data | dict | contains the temporal_baselines and the tcs (time courses).  It should work if not available.  
        fig_title | string | sets the window title and the name of the .png produced.  
        
    Returns:
        Figure
        dem_to_sources_comparisons | dict | results of comparison between spatial sources and DEM.  xyzs: x and y and colours for each point from kernel density estimate, line_xys | list | xy points for line of best fit. corr_coefs | pearson correlation coefficeint for points.  
        tcs_to_tempbaselines_comparisons | dict | as above, but correlations between temporal baselines and time courses (ie if long temrpoal baseline, is component used strongly.  )
        
    History:
        2021_09_12 | MEG | Written.  
        2022_01_14 | MEG | Add return of comparison dicts.  
    """
    import numpy as np
    import numpy.ma as ma
    
    from licsalert.monitoring_functions import update_mask_sources_ifgs
    from licsalert.icasar.aux import signals_to_master_signal_comparison

    def reduce_n_pixs(r2_arrays, n_pixels_new):
        """
        """
        n_pixels_old = r2_arrays[0].shape[1]                                # each pixel is a new column
        index = np.arange(n_pixels_old)
        np.random.shuffle(index)
        r2_arrays_new = []
        index = index[:n_pixels_new]
        for r2_array in r2_arrays:
            r2_arrays_new.append(r2_array[:, index])
        return r2_arrays_new
    
    if fig_title is not None:
        print(f"Starting to create the {fig_title} figure:")
    
    if dem is not None:
        dem_ma = ma.masked_invalid(dem)                                                                                                             # LiCSBAS dem uses nans, but lets switch to a masked array (with nans masked)
        dem_new_mask, sources_new_mask, mask_both = update_mask_sources_ifgs(sources_mask, sources,                             # this takes mask and data as row vectors for one set of masked pixels (the sources from pca) 
                                                                             ma.getmask(dem_ma), ma.compressed(dem_ma)[np.newaxis,:])            # and the mask and data as row vectors from the other set of masked pixels (the DEM, hence why it's being turned into a row vector)
        
        [sources_new_mask, dem_new_mask] = reduce_n_pixs([sources_new_mask, dem_new_mask], max_pixels)                                                          # possibly reduce the number of pixels to speed things up (kernel density estimate is slow)
        dem_to_sources_comparisons = signals_to_master_signal_comparison(sources_new_mask, dem_new_mask, density = True)                                        # And then we can do kernel density plots for each IC and the DEM
        
    else:
        dem_to_sources_comparisons = None
        dem_ma = None
    
    if temporal_data is not None:
        tcs_to_tempbaselines_comparisons = signals_to_master_signal_comparison(temporal_data['tcs'].T, 
                                                                               np.asarray(temporal_data['temporal_baselines'])[np.newaxis,:], density = True)               # And then we can do kernel density plots for each IC and the DEM
    else:
        tcs_to_tempbaselines_comparisons = None
                                                 
    plot_source_tc_correlations(sources, sources_mask, dem_ma, dem_to_sources_comparisons, tcs_to_tempbaselines_comparisons, fig_title = fig_title, **fig_kwargs)       # do the atual plotting
    print("Done.  ")
    return dem_to_sources_comparisons, tcs_to_tempbaselines_comparisons

#%%



def plot_source_tc_correlations(sources, mask, dem = None, dem_to_ic_comparisons = None, tcs_to_tempbaselines_comparisons = None,
                                png_path = './', figures = "window", fig_title = None):
    """Given information about the ICs, their correlations with the DEM, and their time courses correlations with an intererograms temporal basleine, 
    create a plot of this information.  
    Inputs:
        sources | rank 2 array | sources as row vectors.  
        mask | rank 2 boolean | to convert a source from a row vector to a rank 2 masked array.  
        dem | rank 2 array | The DEM.  Can also be a masked array.  
        dem_to_ic_comparisons | dict | keys:
                                        xyzs | list of rank 2 arrays | entry in the list for each signal, xyz are rows.  
                                        line_xys | list of rank 2 arrays | entry in the list for each signal, xy are points to plot for the lines of best fit
                                        cor_coefs | list | correlation coefficients between each signal and the master signal.  
        tcs_to_tempbaselines_comparisons| dict | keys as above.  
        png_path | string | if a png is to be saved, a path to a folder can be supplied, or left as default to write to current directory.  
        figures | string,  "window" / "png" / "png+window" | controls if figures are produced (either as a window, saved as a png, or both)
    Returns:
        figure
    History:
        2021_04_22 | MEG | Written.  
        2021_04_23 | MEG | Update so that axes are removed if they are not being used.  
        
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from licsalert.aux import col_to_ma
    from licsalert.icasar.plotting import remappedColorMap, truncate_colormap

    n_sources = sources.shape[0]

    # colour map stuff
    ifg_colours = plt.get_cmap('coolwarm')
    cmap_mid = 1 - np.max(sources)/(np.max(sources) + abs(np.min(sources)))                                    # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
    if cmap_mid < (1/257):                                                                                                 # this is a fudge so that if plot starts at 0 doesn't include the negative colorus for the smallest values
        ifg_colours_cent = remappedColorMap(ifg_colours, start=0.5, midpoint=0.5, stop=1.0, name='shiftedcmap')
    else:
        ifg_colours_cent = remappedColorMap(ifg_colours, start=0.0, midpoint=cmap_mid, stop=1.0, name='shiftedcmap')


    f, axes = plt.subplots(3, (n_sources+1), figsize = (15,7))
    plt.subplots_adjust(wspace = 0.1)
    f.canvas.manager.set_window_title(f"{fig_title}")
    
    # 1: Plot the DEM:
    if dem is not None:
        terrain_cmap = plt.get_cmap('terrain')
        terrain_cmap = truncate_colormap(terrain_cmap, 0.2, 1)    
        dem_plot = axes[1,0].matshow(dem, cmap = terrain_cmap)
        axin = axes[1,0].inset_axes([0, -0.06, 1, 0.05])
        cbar_1 = f.colorbar(dem_plot, cax=axin, orientation='horizontal')
        cbar_1.set_label('Height (m)', fontsize = 8)
        axes[1,0].set_title('DEM')
        axes[1,0].set_xticks([])
        axes[1,0].set_yticks([])
    else:
        axes[1,0].set_axis_off()
        
    # 2: Find the x and y limits for the 2d scatter plots
    if dem_to_ic_comparisons is not None:                                                               # first check that it actually exists.  
        row1_all_xyzs = np.stack(dem_to_ic_comparisons['xyzs'], axis = 2)                               # merge together into a rank 3 numpy array. (3 x n_pixels x n_ics?)
        row1_xlim = (np.min(row1_all_xyzs[0,]), np.max(row1_all_xyzs[0,]))                              # x limits are min and max of the first row
        row1_ylim = (np.min(row1_all_xyzs[1,]), np.max(row1_all_xyzs[1,]))                              # y limits are min and max of the second row
    
    if tcs_to_tempbaselines_comparisons is not None:                                                    # as above.          
        row2_all_xyzs = np.stack(tcs_to_tempbaselines_comparisons['xyzs'], axis = 2)
        row2_xlim = (np.min(row2_all_xyzs[0,]), np.max(row2_all_xyzs[0,]))        
        row2_ylim = (np.min(row2_all_xyzs[1,]), np.max(row2_all_xyzs[1,]))        
    
    # 3: Loop through each IC
    for ic_n in range(n_sources):
        # 2a: Plotting the IC
        im = axes[0,ic_n+1].matshow(col_to_ma(sources[ic_n,:], mask), cmap = ifg_colours_cent, vmin = np.min(sources), vmax = np.max(sources))
        axes[0,ic_n+1].set_xticks([])
        axes[0,ic_n+1].set_yticks([])
        axes[0,ic_n+1].set_title(f"Source {ic_n}")
        
        # 2B: Plotting the IC to DEM scatter, if the data are available
        if dem_to_ic_comparisons is not None:
            axes[1,ic_n+1].scatter(dem_to_ic_comparisons['xyzs'][ic_n][0,:],
                                   dem_to_ic_comparisons['xyzs'][ic_n][1,:], c= dem_to_ic_comparisons['xyzs'][ic_n][2,:])
            axes[1,ic_n+1].plot(dem_to_ic_comparisons['line_xys'][ic_n][0,:], dem_to_ic_comparisons['line_xys'][ic_n][1,:], c = 'r')
            axes[1,ic_n+1].set_xlim(row1_xlim[0], row1_xlim[1])
            axes[1,ic_n+1].set_ylim(row1_ylim[0], row1_ylim[1])
            axes[1,ic_n+1].axhline(0, c='k')
            axes[1,ic_n+1].yaxis.tick_right()                                       # set ticks to be on the right. 
            if ic_n != (n_sources-1):
                axes[1,ic_n+1].yaxis.set_ticklabels([])                             # if it's not the last one, turn the  tick labels off
            else:
                axes[1,ic_n+1].yaxis.set_ticks_position('right')                    # but if it is, make sure they're on the right.  
                axes[1,ic_n+1].set_ylabel(f"IC")
                axes[1,ic_n+1].yaxis.set_label_position('right')
            if ic_n == int(n_sources/2):
                axes[1,ic_n+1].set_xlabel('Height (m)')
            axes[1,ic_n+1].set_title(f"CorCoef: {np.round(dem_to_ic_comparisons['cor_coefs'][ic_n],2)}", fontsize = 7, color = 'r')    
        else:
            axes[1,ic_n+1].set_axis_off()
                
                

        if tcs_to_tempbaselines_comparisons is not None:
            axes[2,ic_n+1].scatter(tcs_to_tempbaselines_comparisons['xyzs'][ic_n][0,:],
                                   tcs_to_tempbaselines_comparisons['xyzs'][ic_n][1,:], c= tcs_to_tempbaselines_comparisons['xyzs'][ic_n][2,:])
            axes[2,ic_n+1].plot(tcs_to_tempbaselines_comparisons['line_xys'][ic_n][0,:], tcs_to_tempbaselines_comparisons['line_xys'][ic_n][1,:], c = 'r')
            axes[2,ic_n+1].set_xlim(row2_xlim[0], row2_xlim[1])
            axes[2,ic_n+1].set_ylim(row2_ylim[0], row2_ylim[1])                                                    # force them to all share a y axis.  Gnerally not good as such varying scales.  
            axes[2,ic_n+1].axhline(0, c='k')
            axes[2,ic_n+1].yaxis.tick_right()

            if ic_n != (n_sources-1):
                axes[2,ic_n+1].yaxis.set_ticklabels([])                                                            # if it's not the last one, turn the  tick labels off
            else:
                axes[2,ic_n+1].yaxis.set_ticks_position('right')                                                   # but if it is, make sure they're on the right.  
                axes[2,ic_n+1].set_ylabel(f"IC usage strength")
                axes[2,ic_n+1].yaxis.set_label_position('right')
            
            if ic_n == int(n_sources/2):                                                                            # on roughly the middle plot....
                axes[2,ic_n+1].set_xlabel('Temporal Baseline (days)')                                               # add an x label.  
                                
            axes[2,ic_n+1].set_title(f"CorCoef: {np.round(tcs_to_tempbaselines_comparisons['cor_coefs'][ic_n],2)}", fontsize = 7, color = 'r')    
        else:
            axes[2,ic_n+1].set_axis_off()


    # 3: The ICs colorbar
    axin = axes[0,0].inset_axes([0.5, 0, 0.1, 1])            
    cbar_2 = f.colorbar(im, cax=axin, orientation='vertical')
    cbar_2.set_label('IC')
    axin.yaxis.set_ticks_position('left')

    # last tidying up
    for ax in [axes[0,0], axes[2,0]]:
        ax.set_axis_off()

    f.tight_layout()    
    
    if figures == 'window':                                                                 # possibly save the output
        pass
    elif figures == "png":
        f.savefig(f"{png_path}/{fig_title}.png")
        plt.close()
    elif figures == 'png+window':
        f.savefig(f"{png_path}/{fig_title}.png")
    else:
        pass  


#%%



def visualise_ICASAR_inversion(interferograms, sources, time_courses, mask, n_data = 10):
    """Given interferograms (that don't need to be mean centered), visualise how the 
    BSS recontruction (A@S) can be used to reconstruct them.  
    Inputs:
        Interferograms | rank 2 array | interferograms as row vectors
        sources | rank 2 array | souces as row vectors(e.g. 5 x 54000)
        time_courses | rank 2 array | time courses as column vectors (e.g. 85 x 5)
        mask | rank 2 array | to convert a row vector to a rank 2 array (image) for plotting.  
        n_data | int | number of data to plot. 
    Returns:
        figure
    History:
        2021_03_03 | MEG | Written.  
        2021_11_25 | MEG | Write the docs
    """
    import numpy as np
    
    def plot_ifg(ifg, ax, mask, vmin, vmax):
        """
        """
        w = ax.matshow(col_to_ma(ifg, mask), interpolation ='none', aspect = 'equal', vmin = vmin, vmax = vmax)                                                   # 
        axin = ax.inset_axes([0, -0.06, 1, 0.05])
        fig.colorbar(w, cax=axin, orientation='horizontal')
        ax.set_yticks([])
        ax.set_xticks([])
    
    import matplotlib.pyplot as plt
    
    interferograms_mc = interferograms - np.mean(interferograms, axis = 1)[:, np.newaxis]
    interferograms_ICASAR = time_courses @ sources
    residual = interferograms_mc - interferograms_ICASAR
    
    if n_data > interferograms.shape[0]:
        n_data = interferograms.shape[0]

    
    fig, axes = plt.subplots(3, n_data, figsize = (15,7))  
    if n_data == 1:    
        axes = np.atleast_2d(axes).T                                                # make 2d, and a column (not a row)
    
    row_labels = ['Data', 'Model', 'Resid.' ]
    for ax, label in zip(axes[:,0], row_labels):
        ax.set_ylabel(label)

    for data_n in range(n_data):
        vmin = np.min(np.stack((interferograms_mc[data_n,], interferograms_ICASAR[data_n,], residual[data_n])))
        vmax = np.max(np.stack((interferograms_mc[data_n,], interferograms_ICASAR[data_n,], residual[data_n])))
        plot_ifg(interferograms_mc[data_n,], axes[0,data_n], mask, vmin, vmax)
        plot_ifg(interferograms_ICASAR[data_n,], axes[1,data_n], mask, vmin, vmax)
        plot_ifg(residual[data_n,], axes[2,data_n], mask, vmin, vmax)

#%%


def plot_2d_interactive_fig(xy, colours, spatial_data = None, temporal_data = None,
                        inset_axes_side = {'x':0.1, 'y':0.1}, arrow_length = 0.1, figsize = (10,6), 
                        labels = None, legend = None, markers = None, 
                        figures = 'window', png_path = './', fig_filename = '2d_interactive_plot'):
    """ Data are plotted in a 2D space, and when hovering over a point, further information about it (e.g. what image it is)  appears in an inset axes.  
    Inputs:
        xy | rank 2 array | e.g. 2x100, the x and y positions of each data
        colours | rank 1 array | e.g. 100, value used to set the colour of each data point
        spatial_data | dict or None | contains 'images_r3' in which the images are stored as in a rank 3 array (e.g. n_images x heigh x width).  Masked arrays are supported.  
        temporal_data | dict or None | contains 'tcs_r2' as time signals as row vectors and 'xvals' which are the times for each item in the timecourse.   
        inset_axes_side | dict | inset axes side length as a fraction of the full figure, in x and y direction
        arrow_length | float | lenth of arrow from data point to inset axes, as a fraction of the full figure.  
        figsize | tuple |  standard Matplotlib figsize tuple, in inches.  
        labels | dict or None | title for title, xlabel for x axis label, and ylabel for y axis label
        legend | dict or None | elements contains the matplotilb symbols.  E.g. for a blue circle: Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4')
                                labels contains the strings for each of these.       
        markers | dict or None | dictionary containing labels (a numpy array where each number relates to a different marker style e.g. (1,0,1,0,0,0,1 etc))) 
                                 and markers (a list of the different Matplotlib marker styles e.g. ['o', 'x'])
        figures | string,  "window" / "png" / "png+window" | controls if figures are produced (either as a window, saved as a png, or both)
        png_path | string | if a png is to be saved, a path to a folder can be supplied, or left as default to write to current directory.  
        fig_filename | string | name of file, if you wish to set one.  Doesn't include the extension (as it's always a png).  
     Returns:
        Interactive figure
    History:
        2020/09/09 | MEG | Modified from a sript in the ICASAR package.  
        2020/09/10 | MEG | Add labels, and change so that images are stored as rank3 arrays.  
        2020/09/10 | MEG | Add legend option.  
        2020/09/11 | MEG | Add option to have different markers.  
        2020/09/15 | MEG | Add option to set size of inset axes.  
        2021_04_16 | MEG | Add figures option (png, png and window, or just window), option to save to a directory, and option to set filename.  
    
    """
    def remove_axes2_and_arrow(fig):
        """ Given a figure that has a second axes and an annotation arrow due to a 
        point having been hovered on, remove this axes and annotation arrow.  
        Inputs:
            fig | matplotlib figure 
        Returns:
        History:
            2020/09/08 | MEG | Written
        """
        # 1: try and remove any axes except the primary one
        try:
            fig.axes[1].remove()                
        except:
            pass
        
        # 2: try and remove any annotation arrows
        for art in axes1.get_children():
            if isinstance(art, matplotlib.patches.FancyArrow):
                try:
                    art.remove()        
                except:
                    continue
            else:
                continue
        fig.canvas.draw_idle()                                          # update the figure
    
    
    def axes_data_to_fig_percent(axes_lims, fig_lims, point):
        """ Given a data point, find where on the figure it plots (ie convert from axes coordinates to figure coordinates) 
        Inputs:
            axes_xlims | tuple | usually just the return of something like: axes1.get_ylim()
            fig_lims | tuple | the limits of the axes in the figure.  usuall (0.1, 0.9)  for an axes made with something like this: axes1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])                  # main axes
            point |float | point in data coordinates
        Returns:
            fig_position | float | where the data point is in the figure.  (0,0) would be the lower left corner.  
        History:
            2020/09/08 | MEG | Written
            
        """
        gradient = (fig_lims[1] - fig_lims[0])/(axes_lims[1] - axes_lims[0])
        y_intercept = fig_lims[0] - (gradient * axes_lims[0])
        fig_position = (gradient * point) + y_intercept
        return fig_position
    
    def calculate_insetaxes_offset(lims, points, offset_length):
        """
        The offsets between the inset axes and the point are different depending on which quadrant of the graph the point is in.  
        Inputs:
            lims | list | length is equal to the number of dimensions.  Filled with tuples of the axes limits.  
            point | list | length is equal to the number of diemsions. Filled with points.  
            offset_length | float | length of the arrow.  
        Returns:
            offsets | list | length is equal to the number of dimensions.  Length of offset for inset axes in each dimension.  
        History:
            2020/09/08 | MEG | Written
        """
        import numpy as np
        offsets = []
        for dim_n in range(len(lims)):                                        # loop through each dimension.  
            dim_centre = np.mean(lims[dim_n])
            if points[dim_n] < dim_centre:
                offsets.append(-offset_length)
            else:
                offsets.append(offset_length)
        return offsets
    
    def hover(event):
        if event.inaxes == axes1:                                                       # determine if the mouse is in the axes
            cont, ind = sc.contains(event)                                              # cont is a boolean of if hoving on point, ind is a dictionary about the point being hovered over.  Note that two or more points can be in this.  
            if cont:                                                                    # if on point
                remove_axes2_and_arrow(fig)                                             # remove the axes and arrow created when hovering on the point (incase cursor moves from one point to next without going off a point)
                point_n = ind['ind'][0]                                                 # get the index of which data point we're hovering on in a simpler form.      
                
                # 1: Add the annotation arrow (from inset axes to data point)
                arrow_lengths = calculate_insetaxes_offset([axes1.get_xlim(), axes1.get_ylim()], 
                                                          [xy[0,point_n], xy[1,point_n]], arrow_length)                               # calculate the length of the arrow, which depends which quadrant we're in (as the arrow always go away from the plot)
                axes1.arrow(xy[0,point_n] + arrow_lengths[0], xy[1,point_n] + arrow_lengths[1],                                       # add the arrow.  Notation is all a bit backward as head is fixed at end, so it has to be drawn backwards.  
                            -arrow_lengths[0], -arrow_lengths[1], clip_on = False, zorder = 999)                                # clip_on makes sure it's visible, even if it goes off the edge of the axes.  

                # 2: Add the inset axes                
                fig_x = axes_data_to_fig_percent(axes1.get_xlim(), (0.1, 0.9), xy[0,point_n] + arrow_lengths[0])                   # convert position on axes to position in figure, ready to add the inset axes
                fig_y = axes_data_to_fig_percent(axes1.get_ylim(), (0.1, 0.9), xy[1,point_n] + arrow_lengths[1])                   # ditto for y dimension
                if arrow_lengths[0] > 0 and arrow_lengths[1] > 0:                                                          # top right quadrant
                    inset_axes = fig.add_axes([fig_x, fig_y,                                                               # create the inset axes, simple case, anochored to lower left forner
                                               inset_axes_side['x'], inset_axes_side['y']], anchor = 'SW')               
                elif arrow_lengths[0] < 0 and arrow_lengths[1] > 0:                                                        # top left quadrant
                    inset_axes = fig.add_axes([fig_x - inset_axes_side['x'], fig_y,                                        # create the inset axes, nudged in x direction, anchored to lower right corner
                                               inset_axes_side['x'], inset_axes_side['y']], anchor = 'SE')     
                elif arrow_lengths[0] > 0 and arrow_lengths[1] < 0:                                                        # lower right quadrant
                    inset_axes = fig.add_axes([fig_x, fig_y - inset_axes_side['y'],                                        # create the inset axes, nudged in y direction
                                               inset_axes_side['x'], inset_axes_side['y']], anchor = 'NW')                 
                else:                                                                                                      # lower left quadrant
                    inset_axes = fig.add_axes([fig_x - inset_axes_side['x'], fig_y - inset_axes_side['y'],                 # create the inset axes, nudged in both x and y
                                               inset_axes_side['x'], inset_axes_side['y']], anchor = 'NE')                
                
                # 3: Plot on the inset axes
                if temporal_data is not None:
                    inset_axes.plot(temporal_data['xvals'], temporal_data['tcs_r2'][point_n,])                            # draw the inset axes time course graph
                    inset_axes.axhline(0)
                if spatial_data is not None:
                    inset_axes.matshow(spatial_data['images_r3'][point_n,])                                                      # or draw the inset axes image
                inset_axes.set_xticks([])                                                                                       # and remove ticks (and so labels too) from x
                inset_axes.set_yticks([])                                                                                       # and from y
                fig.canvas.draw_idle()                                                                                          # update the figure.  
            else:                                                                       # else not on a point
                remove_axes2_and_arrow(fig)                                             # remove the axes and arrow created when hovering on the point                       
        else:                                                                           # else not in the axes
            remove_axes2_and_arrow(fig)                                                 # remove the axes and arrow created when hovering on the point (incase cursor moves from one point to next without going off a point)
    
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np

    # 1: Check some inputs:
    if temporal_data is None and spatial_data is None:                                                                  # check inputs
        raise Exception("One of either spatial or temporal data must be supplied.  Exiting.  ")
    if temporal_data is not None and spatial_data is not None:
        raise Exception("Only either spatial or temporal data can be supplied, but not both.  Exiting.  ")

    # 2: Draw the figure
    fig = plt.figure(figsize = figsize)                                                                # create the figure, size set in function args.  
    axes1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])                                                         # main axes
    if markers is None:                                                                                # if a dictionary about different markers is not supplied... 
        sc = axes1.scatter(xy[0,],xy[1,],c=colours, s=100)                                             # draw the scatter plot, just draw them all with the default marker
    else:                                                                                                                                     # but if we do have a dictionary of markers.  
        n_markers = len(markers['styles'])                                                                                                    # get the number of unique markers
        for n_marker in range(n_markers):                                                                                                     # loop through each marker style
            point_args = np.ravel(np.argwhere(markers['labels'] == n_marker))                                                                 # get which points have that marker style
            try:
                sc = axes1.scatter(xy[0,point_args], xy[1,point_args], c=colours[point_args], s=100, marker = markers['styles'][n_marker])        # draw the scatter plot with different marker styles
            except:
                pass
        sc = axes1.scatter(xy[0,],xy[1,],c=colours, s=100, alpha = 0.0)                                                                       # draw the scatter plot again with all the points (regardless of marker style), but with invisble markers.  As the last to be drawn, these are the ones that are hovered over, and indexing works as all the points are draw this time.  

    # 3: Try and add various labels from the labels dict
    try:
        fig.canvas.manager.set_window_title(labels['title'])
        fig.suptitle(labels['title'])
    except:
        pass
    try:
        axes1.set_xlabel(labels['xlabel'])
    except:
        pass
    try:
        axes1.set_ylabel(labels['ylabel'])
    except:
        pass
    
    # 4: Possibly add a legend, using the legend dict.  
    if legend is not None:
        axes1.legend(handles = legend['elements'], labels = legend['labels'], 
                     bbox_to_anchor=(1., 0.5), loc = 'center right', bbox_transform=plt.gcf().transFigure)                           # Put a legend to the right of the current axis.  bbox is specified in figure coordinates.  
              
    fig.canvas.mpl_connect("motion_notify_event", hover)                                # connect the figure and the function.  
    
    if figures == 'window':
        pass
    elif figures == "png":
        fig.savefig(f"{png_path}/{fig_filename}.png")
        plt.close()
    elif figures == 'png+window':
        fig.savefig(f"{png_path}/{fig_filename}.png")
    else:
        pass


#%%


def plot_temporal_signals(signals, title = None, signal_names = None,
                          figures = "window",  png_path = './'):
    """Plot a set of time signals stored in a matrix as rows.  
    Inputs:
        signals | rank 2 array | signals as row vectors.  e.g. 1x100
        title | string | figure title.  
        signals_names | list of strings | names of each signal
        figures | string,  "window" / "png" / "png+window" | controls if figures are produced (either as a window, saved as a png, or both)
        png_path | string | if a png is to be saved, a path to a folder can be supplied, or left as default to write to current directory.  
    Returns:
        Figure, either as a window or saved as a png
    History:
        2020/09/09 | MEG | Written
    
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    n_signals = signals.shape[0]
    fig1, axes = plt.subplots(n_signals,1, figsize = (10,6))
    if title is not None:
        fig1.canvas.manager.set_window_title(title)
        fig1.suptitle(title)
    for signal_n, signal in enumerate(signals):
        axes[signal_n].plot(np.arange(0, signals.shape[1]), signal)
        if signal_names is not None:
            axes[signal_n].set_ylabel(signal_names[signal_n])
        axes[signal_n].grid(alpha = 0.5)
        if signal_n != (n_signals-1):
            axes[signal_n].set_xticklabels([])
            
    if figures == 'window':                                                                 # possibly save the output
        pass
    elif figures == "png":
        fig1.savefig(f"{png_path}/{title}.png")
        plt.close()
    elif figures == 'png+window':
        fig1.savefig(f"{png_path}/{title}.png")
    else:
        pass
                          # connect the figure and the function.  


#%%


def plot_pca_variance_line(pc_vals, title = '', figures = 'window', png_path = './'):
    """
    A function to display the cumulative variance in each dimension of some high D data
    Inputs:
        pc_vals | rank 1 array | variance in each dimension.  Most important dimension first.  
        title | string | figure title
        figures | string,  "window" / "png" / "png+window" | controls if figures are produced (either as a window, saved as a png, or both)
        png_path | string or None | if a png is to be saved, a path to a folder can be supplied, or left as default to write to current directory.  
    Returns:
        figure, either as window or saved as a png
    History:
        2019/XX/XX | MEG | Written
        2020/03/03 | MEG | Add option to save as png
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    
    
    f, ax = plt.subplots()
    pc_vals_cs = np.concatenate((np.array([0]), np.cumsum(pc_vals)))
    x_vals = np.arange(len(pc_vals_cs)) 
    ax.plot(x_vals, pc_vals_cs/pc_vals_cs[-1])
    ax.scatter(x_vals, pc_vals_cs/pc_vals_cs[-1])
    ax.set_xlabel('Component number')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel('Cumulative Variance')
    ax.set_ylim([0, 1])
    ax.set_title(title)
    f.canvas.manager.set_window_title(title)
    

    if figures == 'window':
        pass
    elif figures == "png":
        f.savefig(f"{png_path}/01_pca_variance_line.png")
        plt.close()
    elif figures == 'png+window':
        f.savefig(f"{png_path}/01_pca_variance_line.png")
    else:
        pass
        
        
        
        
#%%
        
        



#%% taken from insar_tools.py


#%% Copied from small_plot_functions.py


def r2_arrays_to_googleEarth(images_r3_ma, lons, lats, layer_name_prefix = 'layer', kmz_filename = 'ICs',
                             out_folder = './'):
    """ Given one or several arrays in a rank3 array, create a multilayer Google Earth file (.kmz) of them.  
    Inputs:
        images_r3_ma | rank3 masked array |x n_images x ny x nx
        lons | rank 2 array | lons of each pixel in the image.  
        lats | rank 2 array | lats of each pixel in theimage. 
        layer_name_prefix | string | Can be used to set the name of the layes in the kmz (nb of the form layer_name_prefix_001 etc. )
        kmz_filename | string | Sets the name of the kmz produced
        out_folder | pathlib Path | path to location to save .kmz.  
    Returns:
        kmz file
    History:
        2020/06/10 | MEG | Written
        2021/03/11 | MEG | Update to handle incorrectly sized lons and lats arrays (e.g. rank2 arrays instead of rank 1)
    """
    import numpy as np
    import os
    import shutil
    import simplekml
    from pathlib import Path


    n_images = images_r3_ma.shape[0]    
    if type(out_folder) == str:                                                                                     # this should really be a path, but it could easily be a string.  
        out_folder = Path(out_folder)                                                                               # if it is a string, conver it.  
    # 0 temporary folder for intermediate pngs
    try:
        os.mkdir('./temp_kml')                                                                       # make a temporay folder to save pngs
    except:
        print("Can't create a folder for temporary kmls.  Trying to delete 'temp_kml' incase it exisits already... ", end = "")
        try:
            shutil.rmtree('./temp_kml')                                                              # try to remove folder
            os.mkdir('./temp_kml')                                                                       # make a temporay folder to save pngs
            print("Done. ")
        except:
          raise Exception("Problem making a temporary directory to store intermediate pngs" )

    # 1: Initiate the kml
    kml = simplekml.Kml()
        
    # 2 Begin to loop through each iamge
    for n_image in np.arange(n_images)[::-1]:                                           # Reverse so that first IC is processed last and appears as visible
        layer_name = f"{layer_name_prefix}_{str(n_image).zfill(3)}"                     # get the name of a layer a sttring
        r2_array_to_png(images_r3_ma[n_image,], layer_name, './temp_kml/')              # save as an intermediate .png
        
        ground = kml.newgroundoverlay(name= layer_name)                                 # add the overlay to the kml file
        ground.icon.href = f"./temp_kml/{layer_name}.png"                               # and the actual image part
    
        ground.gxlatlonquad.coords = [(lons[-1,0], lats[-1,0]), (lons[-1,-1],lats[-1,-1]),           # lon, lat of image south west, south east
                                      (lons[0,-1], lats[0,-1]), (lons[0,0],lats[0,0])]         # north east, north west  - order is anticlockwise around the square, startign in the lower left
       
    #3: Tidy up at the end
    kml.savekmz(out_folder / f"{kmz_filename}.kmz", format=False)                                    # Saving as KMZ
    shutil.rmtree('./temp_kml')    


#%% Copied from small_plot_functions.py

def r2_array_to_png(r2, filename, png_folder = './'):    
    """ Given a rank 2 array/image, save it as a png with no borders.  
    If a masked array is used, transparency for masked areas is conserved.  
    Designed for use with Google Earth overlays.  
    
    Inputs:
        r2 | rank 2 array | image / array to be saved
        filename | string | name of .png
        png_folder | string | folder to save in, end with /  e.g. ./kml_outputs/
    Returns:
        png of figure
    History:
        2020/06/10 | MEG | Written
        2021_05_05 | MEG | Change colours to coolwarm.  
        
    """
    import matplotlib.pyplot as plt
    
    f, ax = plt.subplots(1,1)
    ax.matshow(r2, cmap = plt.get_cmap('coolwarm'))
    plt.gca().set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.margins(0,0)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.savefig(f"{png_folder}{filename}.png", bbox_inches = 'tight',pad_inches = 0, transparent = True)
    plt.close()


#%%

def prepare_point_colours_for_2d(labels, cluster_order):
    """Given the label for each point (ie 1, 2 or 3 say, or -1 if noise) and the order of importance to the clusters 
    (ie cluster 3 is the most compact and isolated so has the highest Iq value, then cluster 1, then cluster 2), return 
    a list of colours for each point so they can be plotted using a standard  .scatter funtcion.  Ie all the points labelled
    3 have the same colour.  
    
    Inputs:
        label | rank 1 array | the label showing which cluster each point is in.  e.g. (1000)
        cluster_order | rank 1 array | to determine which cluster should be blue (the best one is always in blue, the 2nd best in orange etc.  )
    Returns:
        labels_chosen_colours | np array | colour for each point.  Same length as label.  
    History:
        2020/09/10 | MEG | Written
        2020/09/11 | MEG | Update so returns a numpy array and not a list (easier to index later on.  )
    """
    import numpy as np
    n_clusters = len(cluster_order)
    
    colours = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']           # the standard nice Matplotlib colours
    if n_clusters > 10:                                                         # if we have more than 10 clsuters, generate some random colours
        for i in range(n_clusters - 10):                                        # how many more colours we need (as we already have 10)
            colours.append(('#%02X%02X%02X' % (np.random.randint(0,25), 
                                               np.random.randint(0,25),
                                               np.random.randint(0,25))))       # generate colours randomly (ie a point between 0 and 255 in 3 dimensions.  )
    else:
        colours = colours[:n_clusters]                                          # crop to length if we have 10 or less colours
        
    colours2 = []                                                               # new list of colours, 1st item is the colour that label 0 should be (which is not necesarily blue)
    for i in range(n_clusters):                                                 # loop through each cluster
        colours2.append(colours[int(np.argwhere(cluster_order == i))])          # populate the list
    
    labels_chosen_colours = []                                           # initiate a list where instead of label for each source, we have its colour
    for label in(labels):                                           # Loop through each point's label
        if label == (-1):                                           # if noise, 
            labels_chosen_colours.append('#c9c9c9')                     # colour is grey
        else:               
            labels_chosen_colours.append(colours2[label])           # otherwise, the correct colour (nb colours 2 are reordered so the most imporant clusters have the usual blue etc. colours)
    labels_chosen_colours = np.asarray(labels_chosen_colours)       # convert from list to numpy array
    return labels_chosen_colours


#%%

def prepare_legends_for_2d(clusters_by_max_Iq_no_noise, Iq):
        """Given the cluster order and the cluster quality index (Iq), create a lenend ready for plot_2d_interactive_fig.  
        Inputs:
            clusters_by_max_Iq_no_noise | rank1 array | e.g. (3,2,4,1) if cluster 3 has the highest Iq.  
            Iq | list | Iq for each clusters.  1st item in list is Iq for 1st cluster.  
        Returns:
            legend_dict | dict | contains the legend symbols (a list of complicated Matplotlib 2D line things), and the labels as a list of strings.  
        History:
            2020/09/10 | MEG | Written
            
        """
        import numpy as np
        from matplotlib.lines import Line2D                                  # for the manual legend
        n_clusters = len(Iq)
        
        legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4'), 
                            Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff7f0e'), 
                            Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ca02c'), 
                            Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728'), 
                            Line2D([0], [0], marker='o', color='w', markerfacecolor='#9467bd'), 
                            Line2D([0], [0], marker='o', color='w', markerfacecolor='#8c564b'), 
                            Line2D([0], [0], marker='o', color='w', markerfacecolor='#e377c2'), 
                            Line2D([0], [0], marker='o', color='w', markerfacecolor='#7f7f7f'), 
                            Line2D([0], [0], marker='o', color='w', markerfacecolor='#bcbd22'), 
                            Line2D([0], [0], marker='o', color='w', markerfacecolor='#17becf')]
        if n_clusters > 10:                                                                  # if we have more than 10 clsuters, repeat the same colours the required number of times
            for i in range(n_clusters-10):
                legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='#%02X%02X%02X' % (np.random.randint(0,255),
                                                                                                                  np.random.randint(0,255),
                                                                                                                  np.random.randint(0,255))))
        legend_elements = legend_elements[:n_clusters]                                      # crop to length
    
        legend_labels = []
        for i in clusters_by_max_Iq_no_noise:
                legend_labels.append(f'#: {i}\nIq: {np.round(Iq[i], 2)} ')                                   # make a list of strings to name each cluster
        legend_labels.append('Noise')
        legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='#c9c9c9'))              # but if we have 10 clusters (which is the max we plot), Noise must be added as the 11th
        legend_dict = {'elements' : legend_elements,
                       'labels'   : legend_labels}
        return legend_dict

#%%


def remappedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the median value of a colormap, and scale the
    remaining color range (i.e. truncate the colormap so that it isn't 
    compressed on the shorter side) . Useful for data with a negative minimum and
    positive maximum where you want the middle of the colormap's dynamic
    range to be at zero.
    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and 0.5; if your dataset mean is negative you should leave 
          this at 0.0, otherwise to (vmax-abs(vmin))/(2*vmax) 
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0; usually the
          optimal value is abs(vmin)/(vmax+abs(vmin)) 
          Only got this to work with:
              1 - vmin/(vmax + abs(vmin))
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          0.5 and 1.0; if your dataset mean is positive you should leave 
          this at 1.0, otherwise to (abs(vmin)-vmax)/(2*abs(vmin)) 
          
      2017/??/?? | taken from stack exchange
      2017/10/11 | update so that crops shorter side of colorbar (so if data are in range [-1 100], 
                   100 will be dark red, and -1 slightly blue (and not dark blue))
      '''
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    
    if midpoint > 0.5:                                      # crop the top or bottom of the colourscale so it's not asymetric.  
        stop=(0.5 + (1-midpoint))
    else:
        start=(0.5 - midpoint)
    
    
    cdict = { 'red': [], 'green': [], 'blue': [], 'alpha': []  }
    # regular index to compute the colors
    reg_index = np.hstack([np.linspace(start, 0.5, 128, endpoint=False),  np.linspace(0.5, stop, 129)])

    # shifted index to match the data
    shift_index = np.hstack([ np.linspace(0.0, midpoint, 128, endpoint=False), np.linspace(midpoint, 1.0, 129)])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    #plt.register_cmap(cmap=newcmap)
    return newcmap
    
    
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    import matplotlib.colors as colors
    import numpy as np
    new_cmap = colors.LinearSegmentedColormap.from_list(
    'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
    cmap(np.linspace(minval, maxval, n)))
    return new_cmap 
