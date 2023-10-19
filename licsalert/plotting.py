#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 15:57:28 2023

@author: matthew
"""

import pdb
import matplotlib.pyplot as plt
#%%


def LiCSAlert_figure(sources_tcs, residual, sources, displacement_r2, n_baseline_end, time_values, day0_date=None,
                     time_value_end=None, figure_type = 'both', figure_out_dir=None, ifg_xpos_scaler = 15, sources_labels = None,
                     cmap = plt.get_cmap('coolwarm')):
    """
    The main fucntion to draw the LiCSAlert figure.  
    
    Inputs:
        sorces_tcs | list of dicts | Each source is an item in the list and has its own dictionary, containing information such as the lines of best, the graident learned
                                     in the baseline stage, the lines-of-best-fit to points distances etc.  
         residual | list of dicts | Same structure as above, but as there is only one residual, list is of length 1.  Shuld contain: cumulative timecourse 
                                    (cumualtive_tc), gradient, lines, sigma, distances, and t_recalculate.  
        sources | r2 array or None | sorces (recoverd by ICASAR) as row vectors.  N.b. must be the same size as the downsampled mask in displacement_r2
                                     If set to None, the full resolution source will be plotted        
        displacement_r2 | dict | contains ifgs as row vectors in "incremental" and their mask ("mask"), and also downsampled versions for faster figures,
                                  and their mask.  Downsampled ones used in plotting!  
        mask | r2 array | to convert an ifg (or source) as a row vector into a rank 2 masked array
        n_baseline_end | int | number of ifgs at which we switch from baseline to monitoring
        time_values | r1 array | time values to end date of each ifg.  e.g. [12,24,36] etc.  Also could be called cumulative baselines
        day0_date | string or None |  date of start of time series / first acquisition.  Used along with time values to make x tick labels as dates

        time_value_end | int or None | if an int, axes are extended to this time value, even if it's larrger than the last value of time_values
                                        N.b. time is in days, so e.g. set it to 96, or 192
        
        figure_outtype |
        figure_outdir  | pathlib Path | 
        
        ifg_xpos_scaler | int | To be positioned correctly in the x direction, the ifgs that are plotted on the upper row must not be taller
                                than the axis they lie within.  Increasing this value makes the ifgs smaller, and therefore fit.  
        WIP sources_labels | dict |  keys and variables:
                                    'defo_sources' : ['dyke', 'sill', 'no_def' ]
                                    'Y_class' : n_sources x 3, labels as one hot encoding.  
                                    'Y_loc' : n_sources x 4, location of deformation.  [0,0,0,0] if none present.  
        
    Returns:
        figure
        
    History:
        2020/01/XX | MEG | Written
        2020/01/10 | MEG | update to add "upper_time_values"
        2020/02/16 | MEG | add ifg_xpos_scaler to make sure ifgs are plotted with the correct x value.  
        2020/03/08 | MEG | Change plotting of ifgs and sources to awlays be the downsampled ones.  
        2020/04/20 | MEG | Update so that x tick labels are dates and not numbers since time series started.  
        2020/06/23 | MEG | Write documentation for dates argument.  
        2020/12/15 | MEG | Determine whether sources are downsampled automatically, and also raise exception if the number of pixels doesn't agree.  
        2021_09_28 | MEG | Update various plotting featuers (add cumulative ifgs, add DEM with lons and lats.  )
        2021_09_29 | MEG | Add ifg date
    
    """
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib import ticker
    import matplotlib as mpl
    from matplotlib.ticker import MultipleLocator
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes            
    
    import datetime as dt 
    
    from licsalert.aux import col_to_ma
        
    def calcualte_line_args(n_ifgs, t_recalculate):
        """Lines of best fit are calculated for eahc time step, but we don't want
        to plot them all (as they lie on top of each other, mostly).  
        Therefore, calcaulte the numbers of which to plot so that they don't overlap.  '
        """
        line_args = []    
        for k in range(n_ifgs):
            if k % t_recalculate == 0 and k != 0:                                   # pick which ones, but have to exclue the 0th one
                line_args.append(k-1)
            if k == (n_ifgs-1) and k not in line_args:                              # and always pick the last one
                line_args.append(k)
        return line_args
    

    def plot_ifgs(ifgs, pixel_mask, figure, gridspec_area, time_values, xlim, ylabel = '', day0_date = None, cumulative = False):
        """ Plot all the ifgs (baseline and monitoring) within the grispec_area.  
        """
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes            
        # 1: Create a single wide axes for the inset axes to be plotted on
        ax_ifgs = plt.Subplot(figure, gridspec_area)                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
        fig1.add_subplot(ax_ifgs)                                                                   # add to figure
        ax_ifgs.set_yticks([])                                                                      # no y ticks
        ax_ifgs.set_ylim(bottom = 0, top = 1)
        ax_ifgs.set_xlim(left = 0, right = xlim)                                                    # set x axis upper limit to be the number of acquisitions in the time series
        if day0_date is not None:       
            xticks_every_3months(ax_ifgs, day0_date, time_values, include_tick_labels = False)      # update the ticks (but not labels) to be the same as the time course and residual axis

        for ifg_n, source in enumerate(ifgs):                                                                                # ifgs are rows, loop through
            iax = ax_ifgs.inset_axes([time_values[ifg_n], 0., (xlim/ifg_xpos_scaler), 1.], transform=ax_ifgs.transData)      # [xpos, ypos, xwidth, ywidth], note the x_pos_scaler that reduces the size of the inset axes to make sure it remains in tehe right place
            #ifg_plot = iax.imshow(col_to_ma(source, pixel_mask), cmap = plt.get_cmap('coolwarm'))                                       # plot on the axes
            ifg_plot = iax.matshow(col_to_ma(source, pixel_mask), cmap = cmap)                                       # plot on the axes
            iax.set_xticks([])                                                                                               # images so get rid of x and y ticks (nb this is just for the inset axes)
            iax.set_yticks([])
            # print(time_values[ifg_n])                                                                       # for debugging
            # plt.pause(1)                                                                                    # "
            if ifg_n == (ifgs.shape[0] -1):                                                                                             # if it's the last interferogram
                cbar_ax = inset_axes(iax, width="7%", height="40%",   loc='lower left',  bbox_to_anchor=(1.05, 0.1, 1, 1),              # Colorbar: isnet axes just to left of the main axis
                                     bbox_transform=iax.transAxes,borderpad=0)
                
                #cbar = fig1.colorbar(ifg_plot, cax = cbar_ax, ticks = [np.nanmin(source), np.nanmax(source)])                            # colorbar, tick only 0 and the max (and check max is not a nan)
                cbar = fig1.colorbar(ifg_plot, cax = cbar_ax)                           
                cbar.set_ticks([np.nanmin(source), np.nanmax(source)])
                cbar.set_ticklabels([f"{np.nanmin(source):.3} m", f"{np.nanmax(source):.3} m"])
                cbar_ax.tick_params(labelsize=6)                                                                                         #                                    
                if day0_date is None:                                                                                                   # if we don't have date information,
                    pass                                                                                                                # we can't add it via the colorbar title
                    #cbar_ax.set_title('LOS disp. (m)', fontsize = 6, loc = 'left')
                else:                                                                                                                    # if we do ahve date information
                    day0_date_dt = dt.datetime.strptime(day0_date, "%Y%m%d")                                                             # label the ifg with date interval, start by converting day0 date
                    if cumulative:                                                                                                       # cumulative ifg spans from start to end
                        ifg_start_date = day0_date_dt
                    else:
                        ifg_start_date = day0_date_dt + dt.timedelta(int(time_values[-2]))                                               # but incremental is from previous date
                    ifg_end_date = day0_date_dt + dt.timedelta(int(time_values[-1]))                                                     # to final date
                    cbar_ax.set_title(f"{dt.datetime.strftime(ifg_start_date, '%Y%m%d')}\n"
                                      f"{dt.datetime.strftime(ifg_end_date, '%Y%m%d')}", fontsize = 6, loc = 'left')
        ax_ifgs.set_ylabel(ylabel, fontsize = 7, labelpad = -1)
        

    def colourbar_for_sources(icasar_sources):
        """ Creat a colourbar for the ICA sources that is centered on 0, and cropped so that each side is equal
        (i.e. if data lies in range [-1 10], will only go slightly blue, but up to max red, with grey at 0)
        """
        
        ics_min = np.min(icasar_sources)                                                       # 
        ics_max = np.max(icasar_sources)
        #ic_colours = plt.get_cmap('coolwarm')
        ic_colours = cmap
        cmap_mid = 1 - ics_max/(ics_max + abs(ics_min))                                     # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
        ic_colours_cent = remappedColorMap(ic_colours, start=0.0, midpoint=cmap_mid, stop=1, name='ic_colours_cent')                    # make the colours for plotting the ICs
        return ic_colours_cent
    
    def sigma_bar_plotter(ax_tc, xvals, yvals, sigma_cmap):
        """
        Add the bar plots onto an axes to show how many sigmas from the line each point is.  
        Inputs:
            ax_tcs | axes object | axes on which to plot
            xvals | x values of bars - usually time
            yvals | height of bars - the number of sigmas from the mean that point is
        """
        ax_tc2 = ax_tc.twinx()                                                       # instantiate a second axes that shares the same x-axis
        for point_n, xval in enumerate(xvals):                                       # loop along each bar
            ax_tc2.bar(xval, yvals[point_n], width=10, alpha = 0.3, color = sigma_cmap(yvals[point_n]/5))   # and plot in required colour
        ax_tc2.set_yticklabels([])                                                      # turn off y labels
        ax_tc2.set_ylim(top = 10)                                                       # set so in range 0 to 10

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        import matplotlib.colors as colors
        new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
        return new_cmap 
    
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
                   

    # -2. Check matplotlib backend is set correctly 
    if figure_type == 'png':
        plt.switch_backend('Agg')                                                           #  works when there is no X11 forwarding, and when displaying plots during creation would be annoying.  
    else: 
        if plt.get_backend() != 'Qt5Agg':                                                               # check what the backend is 
            plt.switch_backend('Qt5Agg')                                                           #  and switch to interactive if it wasn't already.  



    # -1: Check that the sizes of the sources and the interferograms agree.  Raise error if not.  
    if sources.shape[1] == displacement_r2['incremental'].shape[1]:
        sources_downsampled = False
    elif sources.shape[1] == displacement_r2['incremental_downsampled'].shape[1]:
        sources_downsampled = True
    else:
        raise Exception(f"There appears to be a mismatch in the number of pixels contained within the sources ({sources.shape[1]} pixels) "
                        f"and the interferograms ({displacement_r2['incremental'].shape[1]} pixels) or the downsampled interferograms "
                        f"({displacement_r2['incremental_downsampled'].shape[1]} pixels).  The sources must have the same number of pixels "
                        f"as one of these so that their mask can be used to turn the sources from row vectors into images.  ")
        
    
    # 0: Start, some definitions that shouldn't need changing (ie hard coded variables)
    #line_best_fit_alpha = 0.7
    dot_marker_size = 12
    if time_value_end is None:
        t_end = time_values[-1]                                                         # last time value - i.e. the right hand x value of all the axes
    else:
        t_end = time_value_end                                                          # or if it was provided, just use that value
    
    
    # 1 set some preliminary stuff
    t_recalculate = sources_tcs[0]["t_recalculate"]
    n_ics = len(sources_tcs)
    n_times = time_values.shape[0] 
    n_ifgs = displacement_r2["incremental"].shape[0]
    line_args= calcualte_line_args(n_times, t_recalculate)                      # which lines of best fit to plot
    c = mpl.colors.ColorConverter().to_rgb                                      
    cmap_discrete = make_colormap(  [c('black'), c('orange'), 0.33, c('orange'), c('yellow'), 0.66, c('yellow'), c('red')])     # custom colorbar for number of sigmas from line
    cmap_sources = colourbar_for_sources(sources)
    figtitle = f'LiCSAlert figure with {(n_ifgs-n_baseline_end):03d} monitoring interferograms'

    # 2 Initiate the figure    
    fig1 = plt.figure(figsize=(18,10))
    fig1.canvas.manager.set_window_title(figtitle)
    grid = gridspec.GridSpec((n_ics + 3), 11, wspace=0.3, hspace=0.1)                        # divide into 2 sections, 1/5 for ifgs and 4/5 for components

    
    # 3: Plot the ifgs along the top
    plot_ifgs(np.cumsum(displacement_r2["incremental_downsampled"], axis = 0), displacement_r2["mask_downsampled"], fig1, grid[0,1:], 
              time_values, t_end, ylabel = 'Cumulative', day0_date = day0_date, cumulative = True)                                                                # cumulative ifgs
    plot_ifgs(displacement_r2["incremental_downsampled"], displacement_r2["mask_downsampled"], fig1, grid[1,1:], 
              time_values, t_end, ylabel = ' Incremental', day0_date = day0_date, cumulative = False)                                                               # incremental ifgs

    # 4: Plot each source and its time course 
    try:
        baseline_monitor_change = np.mean([time_values[n_baseline_end-1], time_values[n_baseline_end]])                                     # Vertical line will be drawn at this time value to show that we switch from baseline to monitoring
    except:
        baseline_monitor_change = np.mean([time_values[n_baseline_end-1], time_values[n_baseline_end-1] + 12])                              # But the above won't work if there are no monitoring ifgs, so just guess next ifg will be after 12 days and draw line as if that were true (ie 6 days after last point)
    for row_n, source_tc in enumerate(sources_tcs):
        # 4a: Plot the source
        ax_source = plt.Subplot(fig1, grid[row_n+2,0])                                                                                      # create an axes for the IC (spatial source)
       
        if sources_downsampled:
            #im = ax_source.imshow(col_to_ma(sources[row_n], displacement_r2["mask_downsampled"]), cmap = cmap_sources, vmin = np.min(sources), vmax = np.max(sources))   # plot the downsampled source
            im = ax_source.matshow(col_to_ma(sources[row_n], displacement_r2["mask_downsampled"]), cmap = cmap_sources) #, vmin = np.min(sources), vmax = np.max(sources))   # plot the downsampled source
        else:
            #im = ax_source.imshow(col_to_ma(sources[row_n], displacement_r2["mask"]), cmap = cmap_sources, vmin = np.min(sources), vmax = np.max(sources))                # or plot the full resolution source
            #im = ax_source.matshow(col_to_ma(sources[row_n], displacement_r2["mask"]), cmap = cmap_sources, vmin = np.min(sources), vmax = np.max(sources))                # or plot the full resolution source
            im = ax_source.matshow(col_to_ma(sources[row_n], displacement_r2["mask"]), cmap = cmap_sources) #, vmin = np.min(sources), vmax = np.max(sources))                # or plot the full resolution source
            
        ax_source.set_xticks([])
        ax_source.set_yticks([])
        ax_source.set_ylabel(f"IC {row_n}")
        
        # 4b: Possible annotate it with VUDL-net-21 output
        ######################## start WIP
        if sources_labels != None:
            if row_n == 0:
                ax_source.set_title('VUDL-Net-21\nprediction', fontsize = 8, color = 'tab:orange')
            from LiCSAlert_neural_network_functions import centre_to_box, add_square_plot
            #start_stop_locs_pred = centre_to_box(locs_predicted[plot_args[n_plot]])                                         # covert from centre width notation to start stop notation, # [x_start, x_stop, Y_start, Y_stop]
            start_stop_locs = centre_to_box(sources_labels['Y_loc'][row_n])                                         # covert from centre width notation to start stop notation, # [x_start, x_stop, Y_start, Y_stop]
            add_square_plot(start_stop_locs[0], start_stop_locs[1], 
                            start_stop_locs[2], start_stop_locs[3], ax_source, colour='tab:orange')                           # box around deformation
        
            ax_orange = ax_source.twinx()
            ax_orange.set_ylabel(f"{sources_labels['defo_sources'][np.argmax(sources_labels['Y_class'][row_n,])]}", fontsize = 8, color = 'tab:orange')
            #ax_orange.yaxis.label.set_color("tab:orange")
            ax_orange.set_yticks([])
            
        fig1.add_subplot(ax_source)
        ######################## end WIP
        
        
        # 4c: plot the time courses for that IC, and the rolling lines of best fit      
        ax_tc = plt.Subplot(fig1, grid[row_n+2,1:])
        ax_tc.scatter(time_values, source_tc["cumulative_tc"], c = source_tc["distances"], marker='o', 
                      s = dot_marker_size, cmap = cmap_discrete, vmin = 0, vmax = 5, )                                      # note that time_values don't start at 0 as the cumualtive tcs also don't start at 0 
        for line_arg in line_args:                                                                                          # line args sets which lines of best fit to plot (there is a line of best fit for each point, but it's too busy if we plot them all)
            ax_tc.plot(time_values, source_tc["lines"][:,line_arg], c = 'k')                                                # ie each column is a line of best fit
    
        # 4d: tidy up some stuff on the axes
        ax_tc.axhline(y=0, color='k', alpha=0.3)  
        ax_tc.axvline(x = baseline_monitor_change, color='k', alpha=0.3)                          #line the splits between baseline and monitoring ifgs
        ax_tc.set_xlim(left = 0, right = t_end)
        if day0_date is not None:
            xticks_every_3months(ax_tc, day0_date, time_values, include_tick_labels = False)
        fig1.add_subplot(ax_tc)
        sigma_bar_plotter(ax_tc, time_values, source_tc["distances"], cmap_discrete)                # draw the bar graph showing sigma values
        ax_tc.yaxis.tick_right()                                                                    # has to be called after sigma_bar_plotter
        
                                                                
    # 5: Plot the residual
    ax_residual = plt.Subplot(fig1, grid[-1,1:])                                                                    # plot on the last row
    ax_residual.scatter(time_values, residual[0]["cumulative_tc"], marker='o', s = dot_marker_size, cmap = cmap_discrete, vmin = 0, vmax = 5, c = residual[0]["distances"])         # 
    for line_arg in line_args:                                                                                      # plot the rolling line of best fit, but not all of them (only those in line_args)
        ax_residual.plot(time_values, residual[0]["lines"][:,line_arg], c = 'k')                                    # each column is a line of best fit
    ax_residual.axhline(y=0, color='k', alpha=0.3)
    ax_residual.axvline(x = baseline_monitor_change, color='k', alpha=0.3)                          #line the splits between baseline and monitoring ifgs
    ax_residual.set_xlim(left = 0, right = t_end)                    # and finaly tidy up axis and labels etc.  
    ax_residual.yaxis.tick_right()
    ax_residual.yaxis.set_label_position("right")
    ax_residual.set_ylabel('RMS\nresidual')
    fig1.add_subplot(ax_residual)
    sigma_bar_plotter(ax_residual, time_values, residual[0]["distances"], cmap_discrete)                    # draw the bar graph showing sigma values
    ax_residual.yaxis.tick_right()                                                                        # has to be called after sigma_bar_plotter
    if day0_date is not None:
        xticks_every_3months(ax_residual, day0_date, time_values, include_tick_labels = True)                   # create the ticks and labels on the 1st of the quater.  
    
    ## 6: add the two colorbars
    cax = fig1.add_axes([0.12, 0.08, 0.005, 0.1])                                      # source strength
    ics_cbar = fig1.colorbar(im, cax=cax, orientation='vertical')
    tick_locator = ticker.MaxNLocator(nbins=4)
    ics_cbar.locator = tick_locator
    ics_cbar.update_ticks()
    ics_cbar.set_label('IC (m)')
    ics_cbar.ax.yaxis.set_label_position('left')
    
    cax2 = fig1.add_axes([0.17, 0.08, 0.005, 0.1])                                            # number of sigmas from the mean
    norm = mpl.colors.Normalize(vmin=0, vmax=5)
    std_cbar = mpl.colorbar.ColorbarBase(cax2, cmap=cmap_discrete, norm=norm,orientation='vertical')
    tick_locator2 = ticker.MaxNLocator(nbins=5)
    std_cbar.locator = tick_locator2
    std_cbar.update_ticks()
    std_cbar.set_label(r'$\sigma$ from trend line')
    std_cbar.ax.yaxis.set_label_position('left')
    
    # 7: Possibly add the DEM
    if 'dem' in displacement_r2.keys():                                                                                         # DEM is not alway included.  
        ax_dem = plt.Subplot(fig1, grid[1,0])                                                                                   # create an axes for the IC (spatial source)
        terrain_cmap = plt.get_cmap('terrain')                                                                                  # appropriate colours for a dem
        terrain_cmap = truncate_colormap(terrain_cmap, 0.2, 1)                                                                  # but crop (truncate) the blue parts as we are only interested in land
        dem_plot = ax_dem.imshow(displacement_r2["dem"], cmap = terrain_cmap)                                                   # plot the DEM
        ax_dem.xaxis.tick_top()                                                                                                 #
        ax_dem.tick_params(axis='both', which='major', labelsize=7)                                                             # adjust fontsize
        ax_dem.tick_params(axis='both', which='minor', labelsize=7)
        ax_dem.set_xticks([0, displacement_r2['dem'].shape[1]])                                                                 # tick only the min and max in each direction
        ax_dem.set_yticks([0, displacement_r2['dem'].shape[0]])

        ax_dem.xaxis.set_label_position('top')
        if ('lons' in displacement_r2.keys()) and ('lats' in displacement_r2.keys()):                                           # if we have lons and lats, we can update the tick lables to be lons and lats.  
            
            ax_dem.set_xticklabels([str(round(displacement_r2['lons'][-1,0], 2)) + "$^\circ$", str(round(displacement_r2['lons'][-1,-1], 2)) + "$^\circ$" ])
            ax_dem.set_yticklabels([str(round(displacement_r2['lats'][0,0], 2))  + "$^\circ$", str(round(displacement_r2['lats'][-1,0], 2)) + "$^\circ$"])
        # ax_dem.set_ylabel('Latitude ($^\circ$)', fontsize = 6)
        # ax_dem.set_xlabel("Longitude ($^\circ$)", fontsize = 6)
            
        #colorbar for the DEM, just gets in the way.              
        # axins = inset_axes(ax_dem, width="7%", height="50%",   loc='lower left',  bbox_to_anchor=(1.05, 0., 1, 1),              # isnet axes just to left of the main axix for a colorbar
        #                     bbox_transform=ax_dem.transAxes,borderpad=0)
        # #fig1.colorbar(dem_plot, cax = axins, ticks = [0, np.nanmax(displacement_r2['dem'])])                                    # colorbar, tick only 0 and the max (and check max is not a nan)
        # fig1.colorbar(dem_plot, cax = axins, ticks = [])                                    # colorbar, tick only 0 and the max (and check max is not a nan)
        # #axins.tick_params(axis='both', which='major', labelsize=6, rotation = 90)                                               #
        fig1.add_subplot(ax_dem)
        
        # work out the size of the ICs/ DEM and add to the DEM bit of the figure.  
        from geopy import distance
        image_size = {}
        image_size['x'] = int(distance.distance((displacement_r2['lats'][-1,0], displacement_r2['lons'][-1,0]),                       # bottom left corner  
                                                (displacement_r2['lats'][-1,-1], displacement_r2['lons'][-1,-1])).meters / 1000)      #  to bottom right, and convert to integere kms
        image_size['y'] = int(distance.distance((displacement_r2['lats'][-1,0], displacement_r2['lons'][-1,0]),                       # bottom left 
                                            (displacement_r2['lats'][0,0], displacement_r2['lons'][0,0])).meters / 1000)              # to to top left, and conver to integer kms
        
        ax_dem.text(-0.5 * displacement_r2['dem'].shape[1], -0.75 * displacement_r2['dem'].shape[0], f"WxH (km): {image_size['x']} x {image_size['y']}\n"              # add these in these labels in the space above the DEM.  
                        f"DEM (m): {int(np.nanmin(displacement_r2['dem'])), int(np.nanmax(displacement_r2['dem']))}", fontsize = 6 )
    
    # 8: Possible save output
    if (figure_type == 'png') or (figure_type == 'both'):
        filename = "_".join(figtitle.split(" "))                                            # figtitle has spaces, but filename must use underscores instead.  
        fig1.savefig(figure_out_dir / f"{filename}.png", bbox_inches='tight')
        
        
    if plt.get_backend() != 'Qt5Agg':                                                           # check if we need to reset the backend (to interactive window)
        plt.switch_backend('Qt5Agg')                                                           #  





#%%

def LiCSAlert_epoch_figures(displacement_r2_current, reconstructions, residuals, tbaseline_info,
                            figure_type, figure_out_dir):
    """
    Plot last cumulative ifg, last incremetnal ifg, reconstrution for last incremental ifg, and residual for last incremental ifg. 
    Also save the data from these figures as a pickle.  
    
    Inputs:
        displacement_r2_current | dict | licsalert dict of ifgs
        reconstrutions          | r2 array | reconstructions as row vectors
        residuals               | r2 array | reconstructions as row vectors
        tbasline_info           | dict | licsalert dict of tbaseline info
        figure_type             | string | png / window / both
        figure_out_dir          | Path   | out dir path.  
        
    Returns:
        figure and pickle file
        
    History: 
        2023_04_03 | MEG | Written
        2023_10_19 | MEG | Written
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from licsalert.aux import col_to_ma
    import pickle
    
    # 0 Check matplotlib backend is set correctly 
    if figure_type == 'png':
        plt.switch_backend('Agg')                                                           #  works when there is no X11 forwarding, and when displaying plots during creation would be annoying.  
    else: 
        if plt.get_backend() != 'Qt5Agg':                                                               # check what the backend is 
            plt.switch_backend('Qt5Agg')                                                           #  and switch to interactive if it wasn't already.  
        
    ifg_n = displacement_r2_current['incremental'].shape[0]
    inc_ifg_date = tbaseline_info['ifg_dates'][ifg_n-1]
    cum_ifg_date = f"{tbaseline_info['ifg_dates'][0][:8]}_{tbaseline_info['ifg_dates'][ifg_n-1][-8:]}"
    
    # 1: do the plots
    plot_1_image(np.sum(displacement_r2_current['incremental'], axis = 0), displacement_r2_current['mask'], f"01_cumulative_{cum_ifg_date}", 
                 figure_type, figure_out_dir)
    plot_1_image(displacement_r2_current['incremental'][-1, :], displacement_r2_current['mask'], f"02_incremental_{inc_ifg_date}",
                 figure_type, figure_out_dir)
    plot_1_image(reconstructions[:,-1], displacement_r2_current['mask'], f"03_reconstruction_{inc_ifg_date}",
                 figure_type, figure_out_dir)
    plot_1_image(residuals[:, -1], displacement_r2_current['mask'], f"04_residual_{inc_ifg_date}",
                 figure_type, figure_out_dir)

    #  save the data that is plotted as png
    epoch_images = {'cumulative'     : np.sum(displacement_r2_current['incremental'], axis = 0),
                    'incremental'    : displacement_r2_current['incremental'][-1, :],
                    'reconstruction' : reconstructions[:,-1],
                    'residual'       : residuals[:,- 1],
                    'mask'           : displacement_r2_current['mask']}
    with open(figure_out_dir / 'epoch_images_data.pkl', 'wb') as f:
        pickle.dump(epoch_images, f)
    f.close()
    
    # if plt.get_backend() != 'Qt5Agg':                                                           # check if we need to reset the backend (to interactive window)
    #     plt.switch_backend('Qt5Agg')        
        
        
#%%

def LiCSAlert_aux_figures(parent_dir, icasar_sources, dem, mask):
    """Plot the DEM and the sources as png images.  
        Also save the data in these images as a pickle.  
    
    Inputs:
        parent_dir | Path | directory of parent dir, aux_figures is then made in this, and the figures in that.  
        icasar_sources | r2 array | sources (images) as row vectors (e.g. 4 x 8533)
        dem | r2 array | dem, nans for water. 
        mask | r2 array | mask, 1 for masked.  
        
    Returns:
        figures and pickle file.  
        
    History:
        2023_10_19 | MEG | Written
    """
    import numpy.ma as ma
    import pickle
    
    aux_fig_dir = parent_dir / "aux_figures"                                                                                                                    # make a directory to save some aux figures.  
    aux_fig_dir.mkdir(parents=True, exist_ok=True)                                   
    dem_ma_r1 = ma.compressed(ma.array(dem, mask = mask))
    plot_1_image(dem_ma_r1, mask, f"DEM", figure_type = 'png', figure_out_dir = aux_fig_dir, figsize = (18,9))
    for source_n, ic in enumerate(icasar_sources):
        plot_1_image(ic, mask, f"IC_{source_n}", figure_type = 'png', figure_out_dir = aux_fig_dir, figsize = (18,9))
        
    aux_images = {'icasar_sources'  : icasar_sources,
                  'dem'             : dem_ma_r1,
                  'mask'            : mask}
    with open(aux_fig_dir / 'aux_images_data.pkl', 'wb') as f:
        pickle.dump(aux_images, f)
    f.close()    
    
        
#%%

def plot_1_image(im_r1, mask, title, figure_type, figure_out_dir, figsize = (18,9)):
    """Plot a single image (column or row vector) when also given its mask.  
    
    """
    import matplotlib.pyplot as plt
    from licsalert.aux import col_to_ma
    
    
    f, ax = plt.subplots(1,1, figsize = figsize)
    im = ax.matshow(col_to_ma(im_r1, mask))
    f.colorbar(im, label = 'Displacement (m)')
    f.suptitle(title)
    
    if (figure_type == 'png') or (figure_type == 'both'):
        f.savefig(figure_out_dir / f"{title}.png", bbox_inches='tight')
        
    # if figure_type == 'png':
    #     plt.close(f)
    
    
#%%


def LiCSAlert_mask_figure(icasar_mask, licsbas_mask, mask_combined, licsbas_date, current_output_dir, figure_type):
    """ Create a .png showing the licsbas mask, the ICASAR mask, and the current combined mask (ie the pixels in both).  
    
    Inputs:
        icasar_mask | r2 array | the mask used by ICASAR
        licsbas_mask | r2 array | the mask produced by the last run of LiCSBAS
        mask_combined | r2 array | the mask that removes any pixels that aren't in bothh the sources and the ifgs
        licsbas_date | string | the date that LiCSAlert is being run to.  
        current_output_dir | Path | the folder that LiCSALert is currently outputting to
    Returns:
        .png figure
    History:
        2020/06/25 | MEG | Written
        2020/07/01 | MEG | Major rewrite to suit directory based structure.  
        2020/07/03 | MEG | continue major rewrite, and write docs.  
        2021_10_20 | MEG Simplify for LiCSAlert 2.0
    """
    import matplotlib.pyplot as plt
    
    # 0 Check matplotlib backend is set correctly 
    if figure_type == 'png':
        plt.switch_backend('Agg')                                                           #  works when there is no X11 forwarding, and when displaying plots during creation would be annoying.  
    else: 
        if plt.get_backend() != 'Qt5Agg':                                                               # check what the backend is 
            plt.switch_backend('Qt5Agg')                                                           #  and switch to interactive if it wasn't already.  
    
    title = f"LiCSBAS_last_date_{licsbas_date}"

    # 1 Figure showing the masks
    f1,axes = plt.subplots(1,3, figsize = (12,6))
    axes[0].imshow(icasar_mask)
    axes[0].set_title('(ICASAR) sources mask')
    axes[1].imshow(licsbas_mask)
    axes[1].set_title('LiCSBAS mask')
    axes[2].imshow(mask_combined)
    axes[2].set_title('Current combined mask')
    f1.suptitle(title)
    
    f1.canvas.manager.set_window_title(title)
    if (figure_type == 'png') or (figure_type == 'both'):
        f1.savefig(current_output_dir / "mask_status.png", bbox_inches='tight')


#%%

                     

def make_colormap(seq):
    """
    Taken from Stackechange - https://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
    Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).

    e.g. useage:  rvb = make_colormap(  [c('black'), c('orange'), 0.33, c('orange'), c('yellow'), 0.66, c('yellow'), c('red')])
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors


    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)



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
