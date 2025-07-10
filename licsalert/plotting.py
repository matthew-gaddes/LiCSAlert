#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 15:57:28 2023

@author: matthew
"""

import pdb
import matplotlib.pyplot as plt



#%%


def create_manual_mask(image):
    """ 
    Use a matploblib interactive figure to draw a mask manually.  Useful if a 
    user would like to omit part of the scene (e.g. a disconnected island)
    
    Inputs:
        image | r2 numpy array | one channel image, usually cumulative ifg.  
        
    Returns:
        mask | r2 numpy array | True for masked areas.  
        
    History:
        2025_02_21 | MEG | Written

    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
            
    class FreehandMaskDrawer:
        def __init__(self, image, mask_mutable):
            # check if we need to reset the backend (to interactive window)
            if plt.get_backend() != 'Qt5Agg':                                                           
                plt.switch_backend('Qt5Agg')
            
            self.image = image
            self.mask_mutable = mask_mutable
            self.mask = np.zeros_like(image, dtype=bool)
            
            # build the colourmap, coolwarm, with 0 on grey
            cmap = plt.get_cmap('coolwarm')
            im_min = np.min(image)
            im_max = np.max(image)
            # get the ratio of the data that 0 lies at 
            cmap_mid = 1 - im_max/(im_max + abs(im_min))                                     
            
            # make the colours for plotting the ICs, 0 data is middle (grey) in colour
            cmap_remapped = remappedColorMap(
                cmap, start=0.0, midpoint=cmap_mid,
                stop=1, name='ic_colours_cent'
                ) 
            
            self.fig, self.ax = plt.subplots()
            self.ax.imshow(image, cmap=cmap_remapped, origin='upper')
            
            # initiliase to store whether to advance the figure or not
            self.continue_status = False
            
            # Variables to store drawing coordinates
            self.xs, self.ys = [], []
            self.is_drawing = False
            
            # Line object to display the drawing path
            self.line, = self.ax.plot([], [], 'r-', linewidth=2)
            
            # Connect mouse events
            self.cid_press = self.fig.canvas.mpl_connect(
                'button_press_event', self.on_press
                )
            self.cid_motion = self.fig.canvas.mpl_connect(
                'motion_notify_event', self.on_motion
                )
            self.cid_release = self.fig.canvas.mpl_connect(
                'button_release_event', self.on_release
                )
    
        def on_press(self, event):
            # Begin drawing if the click is within the axes.
            if event.inaxes != self.ax:
                return
            self.is_drawing = True
            self.xs = [event.xdata]
            self.ys = [event.ydata]
            self.line.set_data(self.xs, self.ys)
            self.fig.canvas.draw()
    
        def on_motion(self, event):
            # Append coordinates as the mouse moves with the button held down.
            if not self.is_drawing or event.inaxes != self.ax:
                return
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            self.line.set_data(self.xs, self.ys)
            self.fig.canvas.draw()
    
        def on_release(self, event):
            # Finish drawing when the mouse button is released.
            if not self.is_drawing:
                return
            self.is_drawing = False
            
            # Ensure the polygon is closed by appending the first point if needed.
            if (self.xs[0], self.ys[0]) != (self.xs[-1], self.ys[-1]):
                self.xs.append(self.xs[0])
                self.ys.append(self.ys[0])
            
            # Build the polygon from the drawn points.
            polygon = np.column_stack((self.xs, self.ys))
            
            # Create a grid corresponding to the pixel centers.
            ny, nx = self.image.shape
            x, y = np.meshgrid(np.arange(nx), np.arange(ny))
            points = np.vstack((x.ravel(), y.ravel())).T
            
            # Use Path to determine which pixel centers are inside the drawn polygon.
            path = Path(polygon)
            # Using a small negative radius helps include points on the boundary.
            grid = path.contains_points(points, radius=-1e-10)
            self.mask = grid.reshape((ny, nx))
            
            # Optionally, overlay the mask on the original image in the current window.
            self.ax.imshow(
                np.ma.masked_where(~self.mask, self.mask), cmap='jet', 
                alpha=0.5
                )
            self.fig.canvas.draw()
    
            # Disconnect event handlers since drawing is finished.
            self.fig.canvas.mpl_disconnect(self.cid_press)
            self.fig.canvas.mpl_disconnect(self.cid_motion)
            self.fig.canvas.mpl_disconnect(self.cid_release)

            # update status that will exit while that pauses execution 
            # whilst the user draws the mask
            self.continue_status = True
    
            # add to mask to mutable list
            self.mask_mutable.append(self.mask)
            print("Mask created.")
    
        def show(self):

            # enter a loop to pause plotting.  self.continue_status is updated
            # in self.on_release
            while not self.continue_status:
                plt.pause(1.)
                print("Waiting for the user to draw the area to retain...")
    
    # initiliase mask
    mask_drawn = []
    
    drawer = FreehandMaskDrawer(image, mask_mutable = mask_drawn)
    drawer.show()
    
    # invert, so True for masked, False for not masked (LiCSAlert convention)
    mask = ~mask_drawn[0]
    
    return mask
    

#%%


def pngs_to_gif(input_folder, output_file, image_duration = 1000):
    """A function to conbime a folder of .pngs into one gif.
    Inputs:
        input_folder | pathlib Path | the location of the pngs.  
                        E.g. "output_data".  Must be a pathlib Path
        output_file  | string | e.g. "output.gif"
        image_duration | float | the time each png is displayed for in 
                        milliseconds.  e.g. 0.5
    Returns:
        gif
    History:
        2018/07/23 | MEG: adapted from a script.
        2018/12/05 | paths must now be absolute.
        2020/08/18 | MEG | Update to handle paths formatted using the pathlib 
                            module.  
        2021_05_04 | MEG | Update to handle paths better.  
        2023_09_04 | MEG | Add funtion to make sure all images are the same 
                            size.  
        2023_09_09 | MEG | GIF time units have changed in milliseconds in new 
                            version of imageio - update here.  
    """

    import imageio
    import os
    from PIL import Image
    import re

    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [ atoi(c) for c in re.split('(\d+)', text) ]
    
    def determine_minimim_image_size(images):
        """ Iterate along the list of images and find the size of the smallest 
        image
        """
        for image_n, image in enumerate(images):
            ny, nx, _ = image.shape
            if image_n == 0:
                min_y = ny
                min_x = nx
            else:
                if ny < min_y:
                    min_y = ny
                if nx < min_x:
                    min_x = nx
        return min_y, min_x

           
    def crop_if_larger(images, min_y, min_x):
        """ If an image is larger than the minimum, crop it.  
        """
        for image_n, image in enumerate(images):
            ny, nx, _ = image.shape
            if (ny > min_y) or (nx > min_x):
                print(f"Resizing image from {ny} to {min_y} in y, and {nx} to {min_x} in x ")
                images[image_n] = image[:min_y, :min_x, :]
        return images

    # sort the files into a human order
    image_names = os.listdir(input_folder)              
    image_names.sort(key=natural_keys)

    # write the images to the .gif
    images = []
    for counter, filename in enumerate(image_names):
        images.append(imageio.imread(input_folder / filename))
        print(f"Processed {counter} of {len(image_names)} images.  ")
    
    # find the minimum image size
    min_y, min_x = determine_minimim_image_size(images)
    # and possibly crop all images to that size 
    images = crop_if_larger(images, min_y, min_x)

    # for image in images:
    #     print(image.shape)
    
    print('Writing the .gif file...', end = "")
    imageio.mimsave(
        output_file, images, format = 'GIF', duration = image_duration,  
        loop = 1
        )
    print('Done!')


#%%


def offset_volc_lls(volcs, threshold = 0.1, offset = 0.2, attempts = 10,
                    verbose = False):
    """ Given the lis of all volcs, update the lons and lats so that
    when plotted, they do not lie on top of each other.  
    """
    from copy import deepcopy
    from licsalert.plotting import update_overlapping_points
    
    def volc_ll_to_array(volcs):
        """ Extract the lon lats of each volcano from the volc list of objects.  
        """
        import math
        
        lls = []
        indexes = []
        for volc_n, volc in enumerate(volcs):
            if (not math.isnan(volc.lon_lat[0])) and (not math.isnan(volc.lon_lat[1])):
                lls.append(volc.lon_lat)
                indexes.append(volc_n)
        # convert from list of tuples to list of lists
        lls = [list(ll) for ll in lls]
        return lls, indexes

    # extract the lon lats of all the volcs    
    lls_original, indexes = volc_ll_to_array(volcs)

    # repeatdely try to move overlapping points
    attempt = 0
    finished = False   
    lls = deepcopy(lls_original)
    while (finished == False) and (attempt < attempts):
        if verbose:
            print(f"Attempt: {attempt}")
        lls, finished = update_overlapping_points(lls, threshold = threshold,
                                                  offset = offset, verbose = verbose)
        attempt +=1 
    
    # updat the volcs list
    for index_n, index in enumerate(indexes):
        setattr(volcs[index], 'lon_lat_offset', tuple(lls[index_n]))
    
#%%


def update_overlapping_points(shifted_points, threshold = 0.01,
                              offset = 0.02, verbose = False):
    """ When plotting scatter points, multiple close volcanoes can plot on
    top of each other.  Calculate new lons and lats for them that intelligently 
    displace them from their true location.  
    Inputs:
        shifted_points | list of lists | [[lon, lat], [lon, lat]] Points will be
                                        updated so they don't overlap
        threshold | float | points within this distance of each other will be
                            moved.  
        offset | float | how much they will be moved.  
        
    Returns:
        shifted_points | list of lists | as above, but updated.  
        
    History:
        2024_06_21 | MEG | Written
    """
    
    import numpy as np
        
    def move_point(x, y, offset, angle_degrees):
        """
        Move a point (x, y) by a given distance (offset) at a given angle (angle_degrees).
    
        Parameters:
        - x, y: Original coordinates of the point
        - offset: Distance to move the point
        - angle_degrees: Angle at which to move the point, where 0 is up the page
    
        Returns:
        - new_x, new_y: New coordinates of the point
        """
        import math
        angle_radians = math.radians(angle_degrees)  # Convert angle to radians
        new_x = x + offset * math.sin(angle_radians)  # Calculate new x position
        new_y = y + offset * math.cos(angle_radians)  # Calculate new y position
        return new_x, new_y
    
    
    
    def calculate_bearing(lat1, lon1, lat2, lon2):
        """
        Calculate the bearing between two points.
    
        Parameters:
        lat1, lon1 : float
            Latitude and longitude of the first point in degrees.
        lat2, lon2 : float
            Latitude and longitude of the second point in degrees.
    
        Returns:
        float
            Bearing from the first point to the second point in degrees.
        """
        import math
        
        # Convert latitude and longitude from degrees to radians
        lat1 = math.radians(lat1)
        lon1 = math.radians(lon1)
        lat2 = math.radians(lat2)
        lon2 = math.radians(lon2)
    
        # Difference in longitudes
        delta_lon = lon2 - lon1
    
        # Calculate the bearing
        x = math.sin(delta_lon) * math.cos(lat2)
        y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1) * math.cos(lat2) * math.cos(delta_lon))
    
        initial_bearing = math.atan2(x, y)
    
        # Convert bearing from radians to degrees
        initial_bearing = math.degrees(initial_bearing)
    
        # Normalize the bearing to 0°-360°
        compass_bearing = (initial_bearing + 360) % 360
    
        return compass_bearing
    
    
    # start main function    
    for point_n1, point1 in enumerate(shifted_points):
        
        # calculate distance from that point to all others.  
        distances = np.zeros(len(shifted_points))
        for point_n2, point2 in enumerate(shifted_points):
            distances[point_n2] = np.linalg.norm(np.array(point1) - np.array(point2))
        
        # the numbers of the points that overlap
        overlap_args = np.argwhere(distances < threshold)[:,0]
        
        # get the number of points that overlap
        n_overlap = len(overlap_args)
        

        # there is always an "overlap" of 0 as its the point to itself
        if n_overlap > 1:
            
            # get the mean of the overlapping points
            overlapping_points = []
            overlapping_xs = []
            overlapping_ys = []
            for overlap_arg in overlap_args:
                overlapping_points.append(shifted_points[overlap_arg])
                overlapping_xs.append(shifted_points[overlap_arg][0])
                overlapping_ys.append(shifted_points[overlap_arg][1])
            mean_x = np.mean(overlapping_xs)
            mean_y = np.mean(overlapping_ys)
            
            # find the bearing from the mean to each point.  
            bearings = []
            for overlap_arg in overlap_args:
                # lat then lon
                bearings.append(calculate_bearing(mean_y, mean_x,
                                                  shifted_points[overlap_arg][1], 
                                                  shifted_points[overlap_arg][0]))
                
            # get what order the bearings are in
            indices = list(range(len(bearings)))
            indices = sorted(indices, key=lambda i: bearings[i])
                
            # debug check
            # f, ax = plt.subplots()
            # for point_n, point in enumerate(overlapping_points):
            #     ax.scatter(point[0], point[1])
            #     ax.text(point[0], point[1], bearings[point_n])
            # ax.scatter(mean_x, mean_y, marker = '*')
               
            # points move evenly spaced around a circle from original point
            degrees = 360 / n_overlap
            # iterate through the overlapping points and move (note that 
            # iterate in the order given by bearins so that points don't overlap)
            for n, overlap_arg in enumerate(overlap_args[indices]):
                
                if n == 0:
                    original_x = shifted_points[overlap_arg][0]
                    original_y = shifted_points[overlap_arg][1]
                
                new_x, new_y = move_point(shifted_points[overlap_arg][0],
                                          shifted_points[overlap_arg][1],
                                          offset, degrees * n)
                
                shifted_points[overlap_arg][0] = new_x
                shifted_points[overlap_arg][1] = new_y
            if verbose:    
                print(f"{n_overlap} overlapping points at ({original_x}, {original_y}) "
                      f"were dispersed by {offset}   ")
            finished = False
            return shifted_points, finished
                
            
        else:
            # move to the next point to see if it overlaps with anything
            pass
    # never have distances below threshold so finished
    finished = True
    if verbose:
        print(f"No overlapping points (below 'threshold') were found, so exiting.  ")
    return shifted_points, finished


#%% licsalert_status_map()



#%% LiCSAlert_figure()


def LiCSAlert_figure(
        sources_tcs, residual, sources, sources_mask, 
        displacement_r3, figure_date, 
        acq_dates, baselines_cs, baseline_end_date, dayend_date = None,
        figure_type = 'both', figure_out_dir=None, ifg_xpos_scaler = 15,
        sources_labels = None, cmap = plt.get_cmap('coolwarm')
        ):
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

    figure_date
    acq_dates | list | acquisition dates as strings in form yyyymmdd
    baselines_cs | r1 array | temporal baselines of each acquisition from the first acquisition (so starts at 0)
    baseline_end_date 
    dayend_date
    

        
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
        2020/01/10 | MEG | update to add "upper_baselines_cs"
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
    from licsalert.licsalert import licsalert_date_obj
        
    def calcualte_line_args(n_acqs, t_recalculate):
        """Lines of best fit are calculated for each time step, but we don't 
        wantto plot them all (as they lie on top of each other, mostly).  
        Therefore, calcaulte the numbers of which to plot so that they don't 
        overlap. """
        line_args = []    
        for n_acq in range(n_acqs):
            # picn_acq which ones, but have to exclue the 0th one
            if n_acq % t_recalculate == 0 and n_acq != 0:                                   
                line_args.append(n_acq-1)
            # and always pick the last one
            if n_acq == (n_acqs-1) and n_acq not in line_args:                              
                line_args.append(n_acq)
        return line_args
    

    def plot_ifgs(ifgs_r3, figure, gridspec_area, baselines_cs, 
                  ylabel, day0_date, figure_date, dayend_date, 
                  cumulative = False):
        """Plot all the ifgs (baseline and monitoring) within the grispec_area,
        up to the date set by figure_date
        """
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes            
        
        # 1: Create a single wide axes for the inset axes to be plotted on
        ax_ifgs = plt.Subplot(figure, gridspec_area)                                                
        fig1.add_subplot(ax_ifgs)
        ax_ifgs.set_yticks([])
        ax_ifgs.set_ylim(bottom = 0, top = 1)
        # set x axis upper limit to be the number of acquisitions (so figure
        # doesn't change size too much as new data added)
        ax_ifgs.set_xlim(left = 0, right = dayend_date.day_n)                                                    

        # update the ticks (but not labels) to be the same as the time course 
        # and residual axis
        xticks_every_3months(
            ax_ifgs, day0_date.date, baselines_cs, include_tick_labels = False
            )      

        # ifgs are rows, loop through.  Make sure we don't try to iterate 
        # through more than we have dates for  
        for ifg_n, ifg_r2 in enumerate(ifgs_r3[:figure_date.acq_n, :]):                                                         
            #make the ax for the ifg,  x pos, y pox, x width, y width,
            iax = ax_ifgs.inset_axes(
                [baselines_cs[ifg_n+1], 0.,
                 (dayend_date.day_n/ifg_xpos_scaler), 1.],
                transform=ax_ifgs.transData
                )      
            ifg_plot = iax.matshow(ifg_r2, cmap = cmap)                                       
            iax.set_xticks([])                                                                                               
            iax.set_yticks([])
            
            # the last ifg has colorbar and label of range.  
            if ifg_n == (figure_date.acq_n -1): 
                # make an axes for the colourbar
                cbar_ax = inset_axes(
                    iax, width="7%", height="40%",   loc='lower left',  
                    bbox_to_anchor=(1.05, 0.1, 1, 1),
                    bbox_transform=iax.transAxes,borderpad=0
                    )

                cbar = fig1.colorbar(ifg_plot, cax = cbar_ax)                           
                cbar.set_ticks([np.nanmin(ifg_r2), np.nanmax(ifg_r2)])
                cbar.set_ticklabels(
                    [f"{np.nanmin(ifg_r2):.3} m", f"{np.nanmax(ifg_r2):.3} m"]
                    )
                cbar_ax.tick_params(labelsize=6)
                # make the label (that shows the date range of the last ifg)
                if cumulative:                                                                                                       
                    ifg_start_date = day0_date.dt
                else:
                    # need to get date of previous acquisition.  
                    # get the number of days to the acquisition before (which
                    # the incremental ifg is made to), then add this to the 
                    # first date to make a date.  
                    ifg_start_date = (day0_date.dt + 
                                      dt.timedelta(int(
                                          baselines_cs[figure_date.acq_n-1]))) 
                    
                ifg_end_date = figure_date.dt
                cbar_ax.set_title(
                    f"{dt.datetime.strftime(ifg_start_date, '%Y%m%d')}\n"
                    f"{dt.datetime.strftime(ifg_end_date, '%Y%m%d')}", 
                    fontsize = 6, loc = 'left'
                    )
        ax_ifgs.set_ylabel(ylabel, fontsize = 7, labelpad = -1)
        

    def colourbar_for_sources(icasar_sources):
        """ Creat a colourbar for the ICA sources that is centered on 0, and 
        cropped so that each side is equal (i.e. if data lies in range [-1 10],
        will only go slightly blue, but up to max red, with grey at 0)
        """
        
        ics_min = np.min(icasar_sources)                                                       # 
        ics_max = np.max(icasar_sources)
        #ic_colours = plt.get_cmap('coolwarm')
        ic_colours = cmap
        cmap_mid = 1 - ics_max/(ics_max + abs(ics_min))                                     # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
        
        # make the colours for plotting the ICs, 0 data is middle (grey) in colour
        ic_colours_cent = remappedColorMap(ic_colours, start=0.0, midpoint=cmap_mid,
                                           stop=1, name='ic_colours_cent')                    
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

    

    # -3: dealing with time.  
    # if not (sources_tcs[0]['cumulative_tc'].shape)[0] == displacement_r2['incremental_downsampled'].shape[0] == baselines_cs.shape[0] == (len(acq_dates) - 1):          # check all the sizes agree.  
    #     raise Exception(f"The time courses have {(sources_tcs[0]['cumulative_tc'].shape)[0]} entries,  "
    #                     f"there are {displacement_r2['incremental_downsampled'].shape[0]} interferograms, "
    #                     f"the cumulative baselines have {baselines_cs.shape[0]} entries, and "
    #                     f"there are {len(acq_dates)} acquisition dates (which should be 1 more than the others). "
    #                     f"As these don't agree, exiting.  ")
    n_acqs = baselines_cs.shape[0] 
    day0_date = licsalert_date_obj(acq_dates[0], acq_dates)                          # time of first acquisition (x = 0 value)

    if dayend_date is None:
        print(f"No dayend_date was provided which sets the time the plots spans to (although "
              f" the data doesn't necessarily fill the whole plot.  Setting it to the last  "
              f"acquisition date.  ")
        dayend_date = licsalert_date_obj(acq_dates[-1], acq_dates)                         # the highest x value (time) of the figure, although data doesn't necesarily plot to here
    else:
        dayend_date = dayend_date



        

    # -2. Check matplotlib backend is set correctly 
    if figure_type == 'png':
        plt.switch_backend('Agg')                                                           #  works when there is no X11 forwarding, and when displaying plots during creation would be annoying.  
    else: 
        if plt.get_backend() != 'Qt5Agg':                                                               # check what the backend is 
            plt.switch_backend('Qt5Agg')                                                           #  and switch to interactive if it wasn't already.  

    
    # 0: Start, some definitions that shouldn't need changing (ie hard coded variables)
    #line_best_fit_alpha = 0.7
    dot_marker_size = 12
    
    # 1 set some preliminary stuff
    t_recalculate = sources_tcs[0]["t_recalculate"]
    n_ics = len(sources_tcs)
    
    #n_ifgs = displacement_r2["incremental"].shape[0]
    # which lines of best fit to plot (units are acq_n), i.e. don't plot 
    # all the lines or they just overlap, only plot every t_recalculate one
    line_args= calcualte_line_args(n_acqs, t_recalculate)                      
    c = mpl.colors.ColorConverter().to_rgb                                      
     # custom colorbar for number of sigmas from line
    cmap_discrete = make_colormap(
        [c('black'), c('orange'), 0.33, c('orange'), c('yellow'), 0.66,
         c('yellow'), c('red')]
        )
    cmap_sources = colourbar_for_sources(sources)
    figtitle = f'LiCSAlert figure on {figure_date.date}'

    # 2 Initiate the figure    
    fig1 = plt.figure(figsize=(18,10))
    fig1.canvas.manager.set_window_title(figtitle)
    grid = gridspec.GridSpec((n_ics + 3), 11, wspace=0.3, hspace=0.1)                        

    
    # 3: Plot the ifgs along the top
    plot_ifgs(displacement_r3['cum_ma_downsampled'],
              fig1, grid[0,1:], baselines_cs, 'Cumulative', 
              day0_date, figure_date, dayend_date, cumulative = True)
    
    plot_ifgs(displacement_r3['inc_ma_downsampled'],
              fig1, grid[1,1:], baselines_cs, 'Incremental', 
              day0_date, figure_date, dayend_date, cumulative = False)
    
    # 4: Plot each source and its time course 
    for row_n, source_tc in enumerate(sources_tcs):
        # 4a: Plot the source
        ax_source = plt.Subplot(fig1, grid[row_n+2,0])
        im = ax_source.matshow(
            col_to_ma(sources[row_n], sources_mask),
            cmap = cmap_sources) #, vmin = np.min(sources), vmax = np.max(sources))  

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
        ax_tc.scatter(baselines_cs[:figure_date.acq_n+1], source_tc["cumulative_tc"][:figure_date.acq_n+1], 
                      c = source_tc["distances"][:figure_date.acq_n+1], marker='o', 
                      s = dot_marker_size, cmap = cmap_discrete, vmin = 0, vmax = 5, )                                      # note that baselines_cs don't start at 0 as the cumualtive tcs also don't start at 0 
        
        for line_arg in line_args:                                                                                          # line args sets which lines of best fit to plot (there is a line of best fit for each point, but it's too busy if we plot them all)
            if line_arg <= figure_date.acq_n:                                                                                # check that line_arg is not in the future relative to the acquisition date we're currently on
                # ensure that the start acquisition doens't go negative 
                start_acq = np.max((0, (line_arg+1) - source_tc['t_recalculate']))
                ax_tc.plot(baselines_cs[start_acq: line_arg+1], 
                           source_tc["lines"][line_arg], c = 'k')      # 
   
        # 4d: tidy up some stuff on the axes
        ax_tc.axhline(y=0, color='k', alpha=0.3)  
        ax_tc.axvline(x = baseline_end_date.day_n, color='k', alpha=0.3)                          #line the splits between baseline and monitoring ifgs
        ax_tc.set_xlim(left = 0, right = dayend_date.day_n)
        if day0_date is not None:
            xticks_every_3months(ax_tc, day0_date.date, baselines_cs, include_tick_labels = False)
        fig1.add_subplot(ax_tc)
        sigma_bar_plotter(ax_tc, baselines_cs[:figure_date.acq_n], 
                          source_tc["distances"][:figure_date.acq_n], cmap_discrete)                # draw the bar graph showing sigma values
        ax_tc.yaxis.tick_right()                                                                    # has to be called after sigma_bar_plotter
        
                                                                
    # 5: Plot the residual
    ax_residual = plt.Subplot(fig1, grid[-1,1:])                                                                    # plot on the last row
    ax_residual.scatter(baselines_cs[:figure_date.acq_n+1], residual[0]["cumulative_tc"][:figure_date.acq_n+1],
                        marker='o', s = dot_marker_size, cmap = cmap_discrete,
                        vmin = 0, vmax = 5, c = residual[0]["distances"][:figure_date.acq_n+1])         # 
    for line_arg in line_args:                                                                                      # plot the rolling line of best fit, but not all of them (only those in line_args)
        if line_arg < figure_date.acq_n:
            # ax_residual.plot(baselines_cs[:figure_date.acq_n], residual[0]["lines"][:figure_date.acq_n,line_arg], c = 'k')                                    # each column is a line of best fit            
            
            start_acq = np.max((0, (line_arg+1) - source_tc['t_recalculate']))
            ax_residual.plot(baselines_cs[start_acq: line_arg+1],
                             residual[0]["lines"][line_arg], c = 'k')                                    # each column is a line of best fit

            
    ax_residual.axhline(y=0, color='k', alpha=0.3)
    ax_residual.axvline(x = baseline_end_date.day_n, color='k', alpha=0.3)                          #line the splits between baseline and monitoring ifgs
    ax_residual.set_xlim(left = 0, right = dayend_date.day_n)                    # and finaly tidy up axis and labels etc.  
    ax_residual.yaxis.tick_right()
    ax_residual.yaxis.set_label_position("right")
    ax_residual.set_ylabel('RMS\nresidual')
    fig1.add_subplot(ax_residual)
    sigma_bar_plotter(ax_residual, baselines_cs[:figure_date.acq_n], 
                      residual[0]["distances"][:figure_date.acq_n], cmap_discrete)                    # draw the bar graph showing sigma values
    ax_residual.yaxis.tick_right()                                                                        # has to be called after sigma_bar_plotter
    if day0_date is not None:
        xticks_every_3months(ax_residual, day0_date.date, baselines_cs, include_tick_labels = True)                   # create the ticks and labels on the 1st of the quater.  
    
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
    if 'dem' in displacement_r3.keys():                                                                                         # DEM is not alway included.  
        ax_dem = plt.Subplot(fig1, grid[1,0])                                                                                   # create an axes for the IC (spatial source)
        terrain_cmap = plt.get_cmap('terrain')                                                                                  # appropriate colours for a dem
        terrain_cmap = truncate_colormap(terrain_cmap, 0.2, 1)                                                                  # but crop (truncate) the blue parts as we are only interested in land
        dem_plot = ax_dem.imshow(
            displacement_r3["dem"],
            cmap = terrain_cmap
        )
        ax_dem.xaxis.tick_top()                                                                                                 #
        ax_dem.tick_params(axis='both', which='major', labelsize=7)                                                             # adjust fontsize
        ax_dem.tick_params(axis='both', which='minor', labelsize=7)
        ax_dem.set_xticks([0, displacement_r3['dem'].shape[1]])                                                                 # tick only the min and max in each direction
        ax_dem.set_yticks([0, displacement_r3['dem'].shape[0]])

        ax_dem.xaxis.set_label_position('top')
        if ('lons_mg' in displacement_r3.keys()) and ('lats_mg' in displacement_r3.keys()):                                           # if we have lons_mg and lats_mg, we can update the tick lables to be lons_mg and lats_mg.  
            
            ax_dem.set_xticklabels([str(round(displacement_r3['lons_mg'][-1,0], 2)) + "$^\circ$", str(round(displacement_r3['lons_mg'][-1,-1], 2)) + "$^\circ$" ])
            ax_dem.set_yticklabels([str(round(displacement_r3['lats_mg'][0,0], 2))  + "$^\circ$", str(round(displacement_r3['lats_mg'][-1,0], 2)) + "$^\circ$"])
        # ax_dem.set_ylabel('Latitude ($^\circ$)', fontsize = 6)
        # ax_dem.set_xlabel("Longitude ($^\circ$)", fontsize = 6)
            
        #colorbar for the DEM, just gets in the way.              
        # axins = inset_axes(ax_dem, width="7%", height="50%",   loc='lower left',  bbox_to_anchor=(1.05, 0., 1, 1),              # isnet axes just to left of the main axix for a colorbar
        #                     bbox_transform=ax_dem.transAxes,borderpad=0)
        # #fig1.colorbar(dem_plot, cax = axins, ticks = [0, np.nanmax(displacement_r3['dem'])])                                    # colorbar, tick only 0 and the max (and check max is not a nan)
        # fig1.colorbar(dem_plot, cax = axins, ticks = [])                                    # colorbar, tick only 0 and the max (and check max is not a nan)
        # #axins.tick_params(axis='both', which='major', labelsize=6, rotation = 90)                                               #
        fig1.add_subplot(ax_dem)
        
        # work out the size of the ICs/ DEM and add to the DEM bit of the figure.  
        from geopy import distance
        image_size = {}
        image_size['x'] = int(distance.distance((displacement_r3['lats_mg'][-1,0], displacement_r3['lons_mg'][-1,0]),                       # bottom left corner  
                                                (displacement_r3['lats_mg'][-1,-1], displacement_r3['lons_mg'][-1,-1])).meters / 1000)      #  to bottom right, and convert to integere kms
        image_size['y'] = int(distance.distance((displacement_r3['lats_mg'][-1,0], displacement_r3['lons_mg'][-1,0]),                       # bottom left 
                                            (displacement_r3['lats_mg'][0,0], displacement_r3['lons_mg'][0,0])).meters / 1000)              # to to top left, and conver to integer kms
        
        ax_dem.text(-0.5 * displacement_r3['dem'].shape[1], -0.75 * displacement_r3['dem'].shape[0], f"WxH (km): {image_size['x']} x {image_size['y']}\n"              # add these in these labels in the space above the DEM.  
                        f"DEM (m): {int(np.nanmin(displacement_r3['dem'])), int(np.nanmax(displacement_r3['dem']))}", fontsize = 6 )
    
    # 8: Possible save output
    if (figure_type == 'png') or (figure_type == 'both'):
        filename = "_".join(figtitle.split(" "))                                            # figtitle has spaces, but filename must use underscores instead.  
        fig1.savefig(figure_out_dir / f"{filename}.png")# , bbox_inches='tight')            # bbox causes an error "*** AttributeError: 'NoneType' object has no attribute '_get_renderer'
        
    # check if we need to reset the backend (to interactive window)
    if plt.get_backend() != 'Qt5Agg':                                                           
        plt.switch_backend('Qt5Agg')


#%% licsalert_results_explorer()


def licsalert_results_explorer(licsalert_out_dir, fig_width = 18):
    """ The interactive figure for exploring LiCSAlert results.
    
    Inputs:
        licsalert_out_dir | pathlib Path | path to directory of LiCSAlert results.   e.g. Path("./../001_campi_flegrei_example")
        fig_width | int | figure width in inches.  Height is set relative to this and optimised for display on two FHD screens (i.e. 2 x 1920 x 1080)
        
        
    Functions used:
        
        replot_on_ic_select(label)              - just seems to be a string of the name of the IC that is selected.  
        replot_on_pixel_select(event, rax, check)# event is the click, rax is the axes for the radio button, check is the radio button
            
        
        
        plot_ts()
        plot_original_reconstruction()
        
        
        # things to keep track of:
            - ICs
            - pixel
    
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    import matplotlib.gridspec as gridspec
    
    from licsalert.licsalert import reconstruct_ts
    from licsalert.data_importing import open_aux_data, open_tcs, determine_last_licsalert_date
    from licsalert.data_importing import crop_licsalert_results_in_time
    from licsalert.plotting import truncate_colormap, xticks_every_nmonths
    from licsalert.aux import r2_to_r3, moving_average, col_to_ma, determine_abs_max_pixel
    
    
    def plot_ics_ctcs(f, grid, cbar_ax, icasar_sources, mask, sources_tcs, acq_dates, baselines_cumulative):
        """ Plot the ICs and their cumulative time courses.  
        """
        ax_ics = []
        ## Plot the ICs 
        vmin = np.min(icasar_sources)
        vmax = np.max(icasar_sources)
        for source_n, source in enumerate(icasar_sources):
            # plot spatial pattern
            ax_ic = plt.Subplot(f, grid[source_n, 0])                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
            ax_ics.append(ax_ic)
            ic = ax_ic.matshow(col_to_ma(source, mask), vmin = vmin, vmax = vmax)                         # Plot the raw data last cumulative ifg.  
            ax_ic.set_title(f'IC {source_n}: On')
            ax_ic.set_xticks([])
            ax_ic.set_yticks([])
            f.add_subplot(ax_ic)
            
            # plot the cumulative tc
            ax_ctc = plt.Subplot(f, grid[source_n, 1:5])                                                
            ctc = sources_tcs[source_n]['cumulative_tc']
            ctc_smooth, valid = moving_average(ctc)                                                           
            ax_ctc.scatter(baselines_cumulative, ctc, alpha = 0.4, marker = '.', s = 2) 
            ax_ctc.plot(baselines_cumulative, ctc_smooth)                                   
            ax_ctc.axhline(0, c = 'k')
            ax_ctc.grid(True)
            ax_ctc.yaxis.tick_right()
            ax_ctc.tick_params(axis='both', which='both', labelsize=8)
            f.add_subplot(ax_ctc)
            
            if source_n == (n_sources - 1):
                include_tick_labels = True
            else:
                include_tick_labels = False
            xticks_every_nmonths(ax_ctc, acq_dates[0], baselines_cumulative, include_tick_labels = include_tick_labels, 
                                  major_ticks_n_months = 12, minor_ticks_n_months = 1)
            
        f.colorbar(ic, cax = cbar_ax, orientation = 'horizontal')                                    # ICs colorbar, tick only 0 and the max (and check max is not a nan)
        cax_ics.tick_params(axis='both', which='major', labelsize=8, rotation = 315)                                               #
        cax_ics.set_xlabel('IC')            
        return ax_ics
    
    
    def plot_original_reconstruction_dem_resid(fig, ax_orig, ax_reco, ax_resid,
                                               ax_dem, cumulative_r2, 
                                               cumulative_reco_r2, displacement_r2,
                                               pixel, initialise = False):
        """ Make the three plots that show the raw signal (from the input time series) and the reconstruction using various 
        IC components, and the DEM.  
        
        Inputs:
            ax_orig | matplotlib axes | axes to plot the orignal data on
            ax_reco | matplotlib axes | axes to plot the reconstruction data on
            ax_dem  | matplotlib axes | axes to plot the DEM on
            
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
        
        axs = [ax_dem, ax_orig, ax_reco, ax_resid]
        names = ['DEM', 'Original', 'Reconstruction', 'Residual']
        
        for ax, name in zip(axs, names):
            ax.clear()
            # these two share a colorar so work on min and max for both.  
            if  name in ['Original', 'Reconstruction']:
                vmin = np.min(np.concatenate([cumulative_r2[-1,], cumulative_reco_r2[-1]]))
                vmax = np.max(np.concatenate([cumulative_r2[-1,], cumulative_reco_r2[-1]]))
                
                # make a colormap where 0 is grey.  
                cmap_mid = 1 - vmax/(vmax + abs(vmin))
                orig_reco_colours = remappedColorMap(plt.get_cmap('coolwarm'), 
                                                     start=0.0, midpoint=cmap_mid,
                                                     stop=1, name='orig_reco_colours')                    

                if name == 'Original':
                    data_plotted = ax.matshow(col_to_ma(cumulative_r2[-1,], 
                                                        displacement_r2['mask']),
                                              vmin = vmin, vmax = vmax, 
                                              cmap = orig_reco_colours)             
                else:
                    data_plotted = ax.matshow(col_to_ma(cumulative_reco_r2[-1,],
                                                        displacement_r2['mask']),
                                              vmin = vmin, vmax = vmax,
                                              cmap = orig_reco_colours)             
            elif name == 'DEM':
                terrain_cmap = plt.get_cmap('terrain')                                                                                  
                # get rid of the water colours at the bottom.  
                terrain_cmap = truncate_colormap(terrain_cmap, 0.2, 1)                                                                  
                data_plotted = ax_dem.matshow(displacement_r2["dem"], 
                                              cmap = terrain_cmap)                                                
            elif name == 'Residual':
                
                # make a colormap where 0 is grey.  
                residual_r1 = cumulative_r2[-1,] - cumulative_reco_r2[-1,]
                vmax = np.max(residual_r1)
                vmin = np.min(residual_r1)                
                cmap_mid = 1 - vmax/(vmax + abs(vmin))
                resid_colours = remappedColorMap(plt.get_cmap('coolwarm'), 
                                                 start=0.0, midpoint=cmap_mid,
                                                 stop=1, name='resid_colours')                    
                
                # plot the reisudal
                data_plotted = ax.matshow(col_to_ma(residual_r1,
                                                    displacement_r2['mask']), 
                                          cmap = resid_colours)             
            else:
                raise Exception(f"An error has occured when iterating through "
                                f"the four types of data that are plotted.  ")

            ax.scatter(pixel['x'], pixel['y'], s = 20, c = 'r', marker = 'x')

            # cbar
            cax = ax.inset_axes([0.25, 0.0, 0.5 ,0.05]) 
            fig.colorbar(data_plotted, cax=cax, orientation='horizontal')
            cax.xaxis.set_ticks_position('top')
            # text in top left, can get messy if redrawn so only do once.  
            if initialise:
                # update the name to be explicit about how residual calculated
                if name == 'Residual':
                    name = "Residual (Orig. - Reco.)"
                f.text(0.01, 0.99, name, ha='left', va='top',
                       transform = ax.transAxes, color  = 'tab:orange')
            
            ax.set_xticks([])
            ax.set_yticks([])
            fig.add_subplot(ax)
    

    
    
    def plot_ts(fig, ax_ts, cumulative_r2, cumulative_reco_r2, mask, 
                tbaseline_info, pixel):
        """
        Plot the time series for a pixel in two different datasets.  
        
        Inputs:
            fig | matplotlib figure | figure axes is in.  
            ax_ts | matplotlib axes | axes to draw in.  
            cumulative_r2 | r2 array | cumulative ifgs as row vectors.  
            cumulative_reco_r2 | r2 array | reconstrcution of cumultaive data.  
            mask | r2 boolean | true where maskd.  
            tbaseline_info | dict | must contain baselines_cumulative (i.e. 0,6,
                                                                       12,18 for 6 day acquisitions)
            pixel | dict | if None, pixel with maximum deformation is plotted.
                            Or can be dict contiing 'x' and 'y' to choose pixel.  
        Returns:
            Plot in axes
        History:
            2023 | MEG | Written
        """
        from licsalert.aux import r2_to_r3, moving_average
        #from licsalert.plotting import xticks_every_nmonths
        
        data_colours = ['tab:purple', 'tab:orange']
        data_names = ['Original', 'Reconstruction']
        
        cumulative_r3 = r2_to_r3(cumulative_r2, mask)
        cumulative_reco_r3 = r2_to_r3(cumulative_reco_r2, mask)
        

        # iterate through the original and reconstructed data        
        for i, data in enumerate([cumulative_r3, cumulative_reco_r3]):                                      
            # get the pixel to be plotted
            ts = data[:, pixel['y'], pixel['x']]                                                            
            ts_smooth, valid = moving_average(ts)                                                           
            ax_ts.scatter(tbaseline_info['baselines_cumulative'], ts, alpha = 0.4,
                           marker = '.', s = 4,  label = data_names[i], c = data_colours[i])                                       
            ax_ts.plot(tbaseline_info['baselines_cumulative'], ts_smooth,
                       c = data_colours[i])                                   
        
        
        ax_ts.axhline(0, c = 'k')
        ax_ts.grid(True)
        ax_ts.yaxis.tick_right()
        fig.add_subplot(ax_ts)                                                                   
        ax_ts.yaxis.set_label_position("right")
        ax_ts.set_ylabel("LOS displacemnt (m)")
        
        xticks_every_nmonths(ax_ts, tbaseline_info['acq_dates'][0], 
                             tbaseline_info['baselines_cumulative'], 
                             include_tick_labels = True, 
                              major_ticks_n_months = 12, minor_ticks_n_months = 1)
    
        ax_ts.legend(markerscale=5)
        
    
    def replot_on_click(event): 
        """ When the selected pixel changes, replot the time series for that 
        point and the three images with a point showing where was clicked.  
        """
        # if the click is in one of the images axes, that changes the point plotted
        if (event.inaxes is ax_reco) or (event.inaxes is ax_cum) or (event.inaxes is ax_dem) or (event.inaxes is ax_resid):                          
            figure_status['pixel']['x'] = int(event.xdata)                                                                               
            figure_status['pixel']['y'] = int(event.ydata)
            replot()
            
        # if the click is in one of the ic axes, that turns them on and off
        for ic_n, ax in enumerate(figure_status['ic_axs']):
            if event.inaxes is ax:
                # reverse that status of the ic clicked
                figure_status['ic_status'][ic_n] = not figure_status['ic_status'][ic_n]                     # switch from on to off, or off to on.  
                if figure_status['ic_status'][ic_n]:
                    ax.set_title(f"IC {ic_n}: On")
                else:
                    ax.set_title(f"IC {ic_n}: Off")
                    
                # remake the time series using the ICs selected
                _, cumulative_reco_r2 = reconstruct_ts(figure_status['ic_status'], sources_tcs, aux_data, displacement_r2)                               
                figure_status['cumulative_reco_r2'] = cumulative_reco_r2                                                            # store in the mutable object
                replot()
                
                
    
    def replot():
        """ When an update to the figure is required, replot the images and the time series for a pixel.  
        """
        for ax in ax_cum, ax_reco, ax_resid, ax_dem:
            ax.clear()
        
        # plot original, DEM, reconstruction, residual
        plot_original_reconstruction_dem_resid(f, ax_cum, ax_reco, ax_resid, 
                                               ax_dem,  cumulative_r2, 
                                               figure_status['cumulative_reco_r2'], 
                                               displacement_r2, figure_status['pixel'])                            
        # clear the time series plot ready for new plot
        ax_ts.clear()                                                                                                          
        
        # replot the time series using the new reconstruction.  
        plot_ts(f, ax_ts, cumulative_r2, figure_status['cumulative_reco_r2'],    
                displacement_r2['mask'], tbaseline_info, figure_status['pixel'])                                                                 
    
        plt.draw()
    
    print(f"Starting the LiCSAlert interactive time series explorer.  ")
    
    # this figure only works interactively, so ensure backend is set for this
    if plt.get_backend() != 'Qt5Agg':                                                               
        print(f"Updating the backend to Qt5Agg as this figure is interactive.  ")
        plt.switch_backend('Qt5Agg')                                                           
      
    # 1: Open data and some simple processing 
    displacement_r3, tbaseline_info, aux_data = open_aux_data(licsalert_out_dir)
    final_date_dir = determine_last_licsalert_date(licsalert_out_dir)
    sources_tcs, residual_tcs = open_tcs(final_date_dir)    
    
    print(f"A LiCSAlert directory for {final_date_dir.parts[-1]} was found "
          f"so the interactive figure will be created up to this date.  ")
    
    crop = crop_licsalert_results_in_time(
        final_date_dir.parts[-1], tbaseline_info['acq_dates'],
        sources_tcs, residual_tcs,
        None, None, displacement_r3, tbaseline_info
        )
    
    sources_tcs, residual_tcs, _, _, displacement_r3, tbaseline_info = crop; del crop
    
    n_sources = len(sources_tcs)
    n_pixels = np.size(displacement_r2['incremental'], axis = 1)
    
    # convert incremenal to cumulative displacements, 0s on first acquisition. 
    cumulative_r2 = np.concatenate((np.zeros((1, n_pixels)), 
                                    np.cumsum(displacement_r2['incremental'], axis = 0)), axis = 0)
    cumulative_r3 = r2_to_r3(cumulative_r2, displacement_r2['mask'])                                                                                                         # conver to rank 3
    
    # reconstruct cumulative data
    _, cumulative_reco_r2 = reconstruct_ts([1 for i in range(n_sources)], 
                                           sources_tcs, aux_data, displacement_r2)

    x, y, col = determine_abs_max_pixel(cumulative_r3, cumulative_r2)
    pixel = {'x' : x, 'y' : y}                                                                                                                                              
    del y, x, cumulative_r3, col
    
    # 2: start the figure.      
    f = plt.figure(figsize = (fig_width, fig_width /  (2 * (1920 / 1080))))
    f.canvas.manager.set_window_title('LiCSAlert results visualiser')
    
    # set the number of rows to always be even by adding one if needed
    if n_sources % 2 == 1:                                                                                  # number of rows is the smallest even number greater than or equal to the number of sources.  
        n_rows = n_sources + 1
    else:
        n_rows = n_sources
    grid = gridspec.GridSpec(n_rows, 20, wspace=0.2, hspace=0.2)                                            # 
    
    
    # create all the axes from the grid, top row first
    ax_cum   = plt.Subplot(f, grid[:int(n_rows/2), 10:15])                                                
    ax_dem   = plt.Subplot(f, grid[:int(n_rows/2), 5:10])
    ax_reco  = plt.Subplot(f, grid[int(n_rows/2): , 10:15])
    ax_resid = plt.Subplot(f, grid[int(n_rows/2):, 5:10])
    # axes for the line graps (ax_ts) and the ICs colorbar
    ax_ts = plt.subplot(grid[:5,15:])
    cax_ics = f.add_axes([0.125, 0.11, 0.03, 0.02])
    
    # 3: start plotting.  
    # plot DEM, original, residual, reconstruction.  
    plot_original_reconstruction_dem_resid(f, ax_cum, ax_reco, ax_resid, ax_dem, 
                                           cumulative_r2, cumulative_reco_r2, 
                                           displacement_r2, pixel, True)     

    # plot the time series for the point of interest in the two ways.      
    plot_ts(f, ax_ts, cumulative_r2, cumulative_reco_r2, displacement_r2['mask'], 
            tbaseline_info, pixel)    
    
    ax_ics = plot_ics_ctcs(f, grid, cax_ics, aux_data['icasar_sources'], 
                           displacement_r2['mask'],  sources_tcs, 
                           tbaseline_info['acq_dates'], 
                           tbaseline_info['baselines_cumulative'])
    
    # 3:  start the interactive bit
    # mutable object that is mutated by the interactive functions.  
    # get current pixel that time series is plotted for.  
    figure_status = {'pixel'                : pixel,                                                                               
                     "cumulative_reco_r2"   : cumulative_reco_r2,                                                    
                     "ic_axs"               : ax_ics,
                     "ic_status"            : [True for i in range(n_sources)]}                           
    
    cid = f.canvas.mpl_connect('button_press_event', replot_on_click)                   # click on point to plot it.  



#%% LiCSAlert_epoch_figures

def LiCSAlert_epoch_figures(
        processing_date, displacement_r3_current,  reconstructions, residuals,
        tbaseline_info, figure_type, figure_out_dir):
    """
    Plot last cumulative ifg, last incremetnal ifg, reconstrution for last 
    incremental ifg, and residual for last incremental ifg. 
    Also save the data from these figures as a pickle.  
    
    Note that reconstructions and residuals are of the incremental displacements,
    so are 1 shorter in time than the cum_ma in displacement_r3
    
    
    Inputs:
        processing_date          | licsalert date | processing date currently on
        displacement_r3_current | dict | licsalert dict of ifgs
        reconstrutions          | r3 array | reconstructions of incremental data
        residuals               | r3 array | residuals, Assumed not to be cumulative.     
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
        plt.switch_backend('Agg')                                                           
    else: 
        if plt.get_backend() != 'Qt5Agg':                                                               
            plt.switch_backend('Qt5Agg')                                                           
  
    if processing_date.acq_n == 0:
        # this is just the same date repeated as the first acquisition so no ifg
        inc_ifg_date = (
            f"{tbaseline_info['acq_dates'][0]}_"
            f"{tbaseline_info['acq_dates'][processing_date.acq_n]}"
        )
    else:
        # the date of the acquistiions before this one to the current one
        inc_ifg_date = (
            f"{tbaseline_info['acq_dates'][processing_date.acq_n-1]}_"
            f"{tbaseline_info['acq_dates'][processing_date.acq_n]}"
        )
        
    # the first acquisition date to the current one.  
    cum_ifg_date = (
        f"{tbaseline_info['acq_dates'][0]}_"
        f"{tbaseline_info['acq_dates'][processing_date.acq_n]}"
    )

    
    # 1: do the plots
    #cum_1 = np.sum(displacement_r3_current['inc_ma'][:processing_date.acq_n,], axis = 0)
    cum_1 = displacement_r3_current['cum_ma'].original[processing_date.acq_n,]

    # the first date is a special case as there is no data yet    
    if processing_date.acq_n == 0:
        zeros = np.zeros((1, displacement_r3_current['inc_ma'].shape[1]))
        inc_1 = zeros
        recon_1 = zeros
        residual_1 = zeros
        residual_cum = zeros
    else:
        # note -1s here as no incremental data on epoch 0
        inc_1 = displacement_r3_current['inc_ma'].original[processing_date.acq_n-1,]
        recon_1 = reconstructions[processing_date.acq_n-1,]
        residual_1 = residuals[processing_date.acq_n-1,]
        # sum all the residuals through time to get current cumulative
        # residual
        residual_cum = np.sum(residuals[:processing_date.acq_n-1,], axis = 0)
        
    # output each of the pngs
    plot_1_image(cum_1, f"01_cumulative_{cum_ifg_date}", 
                 figure_type, figure_out_dir, cmap = plt.get_cmap('coolwarm'))    
    
    plot_1_image(inc_1, f"02_incremental_{inc_ifg_date}",
                 figure_type, figure_out_dir, cmap = plt.get_cmap('coolwarm'))
    
    plot_1_image(recon_1,  f"03_incremental_reconstruction_{inc_ifg_date}",
                 figure_type, figure_out_dir, cmap = plt.get_cmap('coolwarm'))
    
    plot_1_image(residual_1, f"04_incremental_residual_{inc_ifg_date}",
                 figure_type, figure_out_dir, cmap = plt.get_cmap('coolwarm'))
    
    plot_1_image(residual_cum,  f"05_cumulative_residual_{cum_ifg_date}",
                 figure_type, figure_out_dir, cmap = plt.get_cmap('coolwarm'))

    #  save the data that is plotted as png
    epoch_images = {
        'cumulative'     : cum_1,
        'incremental'    : inc_1,
        'reconstruction' : recon_1,
        'residual'       : residual_1,
        'mask'           : displacement_r3_current['cum_ma'].original.mask[processing_date.acq_n,],
        'residual_cum'   : residual_cum
    }
    with open(figure_out_dir / 'epoch_images_data.pkl', 'wb') as f:
        pickle.dump(epoch_images, f)
    f.close()
    
    # if plt.get_backend() != 'Qt5Agg':                                                           
    #     plt.switch_backend('Qt5Agg')        
        
        
#%% LiCSAlert_aux_figures

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
    
    from licsalert.aux import col_to_ma
    
    aux_fig_dir = parent_dir / "aux_data_figs"                                                                                                                    # make a directory to save some aux figures.  
    aux_fig_dir.mkdir(parents=True, exist_ok=True)                                   

    plot_1_image(
        dem,
        f"DEM",
        figure_type='png',
        figure_out_dir=aux_fig_dir,
        figsize = (18,9)
    )
    
    for source_n, ic in enumerate(icasar_sources):
        source_r2=col_to_ma(ic, mask)
        plot_1_image(
            source_r2,
            f"IC_{source_n}",
            figure_type='png',
            figure_out_dir=aux_fig_dir,
            figsize = (18,9)
        )
        
    aux_images = {'icasar_sources'  : icasar_sources,
                  'dem'             : dem,
                  'mask'            : mask}
    with open(aux_fig_dir / 'aux_images_data.pkl', 'wb') as f:
        pickle.dump(aux_images, f)
    f.close()    
    
        
#%% plot_1_image

def plot_1_image(im_r2,
                 title,
                 figure_type,
                 figure_out_dir,
                 figsize = (18,9),
                 cmap =  plt.get_cmap('viridis')
     ):
    """Plot a single image (column or row vector) when also given its mask.  
    
    """
    import matplotlib.pyplot as plt
    from licsalert.aux import col_to_ma
    
    
    f, ax = plt.subplots(1,1, figsize = figsize)
    im = ax.matshow(im_r2, cmap = cmap)
    f.colorbar(im, label = 'Displacement (m)')
    f.suptitle(title)
    
    if (figure_type == 'png') or (figure_type == 'both'):
        f.savefig(figure_out_dir / f"{title}.png", bbox_inches='tight')
        
    # if figure_type == 'png':
    #     plt.close(f)
    
    
#%% LiCSAlert_mask_figure


def LiCSAlert_mask_figure(
        mask_icasar, mask_epoch, 
        licsbas_date, 
        current_output_dir,
        figure_type
        ):
    """ Create a .png showing the licsbas mask, the ICASAR mask, and the 
    current combined mask (ie the pixels in both).  
    
    Inputs:
        mask_icasar | r2 array | the mask used by ICASAR
        mask_epoch | r2 array | the mask produced by the last run of LiCSBAS
        mask_combined | r2 array | the mask that removes any pixels that aren't
                                   in both the sources and the ifgs
        licsbas_date | string | the date that LiCSAlert is being run to.  
        current_output_dir | Path | the folder that LiCSALert is currently 
                                    outputting to
    Returns:
        .png figure
    History:
        2020/06/25 | MEG | Written
        2020/07/01 | MEG | Major rewrite to suit directory based structure.  
        2020/07/03 | MEG | continue major rewrite, and write docs.  
        2021_10_20 | MEG Simplify for LiCSAlert 2.0
        2025_06_25 | MEG | Change to per epoch masking for LiCSAlert 4.  
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    # # 0 Check matplotlib backend is set correctly 
    # if figure_type == 'png':
    #     plt.switch_backend('Agg')
    # else: 
    #     if plt.get_backend() != 'Qt5Agg':
    #         plt.switch_backend('Qt5Agg')
    
    n_pixels = mask_icasar.shape[0] * mask_icasar.shape[1]
    
    title = f"LiCSBAS_last_date_{licsbas_date}"

    # 1 Figure showing the masks
    f1,axes = plt.subplots(1,5, figsize = (18,6))
    axes[0].imshow(mask_icasar)
    axes[0].set_title(f'ICA mask\n({n_pixels - np.sum(mask_icasar)} pixels)')
    
    axes[1].imshow(mask_epoch)
    axes[1].set_title(
        f'Current epoch mask\n({n_pixels - np.sum(mask_epoch)} pixels)'
        )
    
    # find the pixels that are not masked in both
    mask_combined = np.invert(
        np.logical_and(
            np.invert(mask_icasar),
            np.invert(mask_epoch)
            )
        )
    
    axes[2].imshow(mask_combined)
    axes[2].set_title(
        f'Combined mask\n({n_pixels - np.sum(mask_combined)} pixels)'
        )
    
    # find the pixels that are valid in ICA, but not the epoch
    missing_epoch = np.invert(
        np.logical_and(
            np.invert(mask_icasar),
            mask_epoch
            )
        )
    axes[3].imshow(missing_epoch)
    axes[3].set_title(
        f'Valid in ICA but not epoch\n({n_pixels-np.sum(missing_epoch)} pixlels)'
        )
    
    
    # find the pixels that are valid in the epoch, but not ICA
    missing_ica = np.invert(
        np.logical_and(
            mask_icasar,
            np.invert(mask_epoch)
            )
        )
    axes[4].imshow(missing_ica)
    axes[4].set_title(
        f'Valid in epoch but not ICA\n({n_pixels-np.sum(missing_ica)} pixlels)'
        )
    
    f1.suptitle(title)
    
    f1.canvas.manager.set_window_title(title)
    if (figure_type == 'png') or (figure_type == 'both'):
        f1.savefig(
            current_output_dir / "mask_status.png", 
            bbox_inches='tight'
            )


#%% xticks_every_nmonths()


def xticks_every_nmonths(ax_to_update, day0_date, time_values, include_tick_labels, 
                         major_ticks_n_months = 3, minor_ticks_n_months = 1):
    """Given an axes, update the xticks so the major ones are the 1st of every n months (e.g. if every 3, would be: jan/april/july/october).  
    
    Inputs:
        ax_to_update | matplotlib axes | the axes to update.  
        day0_date | string | in form yyyymmdd
        time_values | rank 1 array | cumulative temporal baselines, e.g. np.array([6,18, 30, 36, 48])
        include_tick_labels | boolean | if True, tick labels are added to the ticks.  
        n_months | int | x ticks are very n months.  e.g. 2, 3,4,6,12 (yearly)  Funny spacings (e.g. every 5) not tested.  
    Returns:
        updates axes
    History:
        2021_09_27 | MEG | Written
        2022_02_17 | MEG | modify so can be monhtly spacing other than every 3 months.  
    """
    import pdb
    import numpy as np
    import datetime as dt
    import copy
    
    from dateutil.relativedelta import relativedelta                                                    # add 3 months and check not after end
    #import matplotlib.pyplot as plt

    
    def create_tick_labels_every_nmonths(day0_date_dt, dayend_date_dt, n_months = 1):
        """ Given a spacing of every n_months, get the dates and days since the first date for ticks every n_months.  
        e.g. every month, every 6 months.  
        
        Inputs:
            day0_date_dt | datetime | date of x = 0 on axis.  
            dayend_date_dt | datetime | date of last x value. 
            n_months | int | frequency of ticks. 
            
        Returns:
            ticks | dict | contains datetimes : datetimes for each tick
                                    yyyymmdd : strings to use as labels in form yyyy/mm/dd
                                    n_day     : day number of tick.  
        History:
            2022_03_29 | MEG | Written
        """
    
        # 1: find first tick date (the first of the jan/ april/jul /oct)                        
        date_tick0 = copy.deepcopy(day0_date_dt)                                                                 # version that can be modified as we iterate through.  
        while not ( (date_tick0.day) == 1 and (date_tick0.month in (np.arange(0, 12, n_months) + 1))):           # i.e. whilst it's not the 1st of the month, and not jan/apr/jul/oct....
            date_tick0 +=  dt.timedelta(1)                                                                       # then add one day and keep going.  
    
        # 2: get all the other first of the quarters as datetimes (ie keep adding n months until we're gone past the day end date)
        ticks = {'datetimes' : [date_tick0],
                 'yyyymmdd'   : [],
                 'n_day'     : []}
       
        while ticks['datetimes'][-1] < (dayend_date_dt - relativedelta(months=+n_months)):                         # while we haven't gone past the last date (and subtract 3 months to make sure we don't go one 3 month jump too far. )
            ticks['datetimes'].append(ticks['datetimes'][-1] + relativedelta(months=+n_months))                    # append the next date which is n_months more.  
        
        # 3: work out what day number each first of the quarter is.  
        for tick_dt in ticks['datetimes']:                                                                      # loop along the list of datetimes (which are each tick) 
            ticks['yyyymmdd'].append(dt.datetime.strftime(tick_dt, "%Y/%m/%d"))                                 # as a string that can be used for the tick label (hence why include / to make more readable)
            ticks['n_day'].append((tick_dt - day0_date_dt).days)
            
        return ticks

    
    xtick_label_angle = 315                                                                              # this angle will read from top left to bottom right (i.e. at a diagonal)
    
    tick_labels_days = ax_to_update.get_xticks().tolist()                                                # get the current tick labels
    day0_date_dt = dt.datetime.strptime(day0_date, "%Y%m%d")                                             # convert the day0 date (date of day number 0) to a datetime.  
    dayend_date_dt = day0_date_dt +  dt.timedelta(int(time_values[-1]))                                  # the last time value is the number of days we have, so add this to day0 to get the end.  

    ticks_major = create_tick_labels_every_nmonths(day0_date_dt, dayend_date_dt, n_months = major_ticks_n_months)
    ticks_minor = create_tick_labels_every_nmonths(day0_date_dt, dayend_date_dt, n_months = minor_ticks_n_months)                        # these are used as the minor ticks every month.  
    

        
    # 4: Update the figure.  
    ax_to_update.set_xticks(ticks_major['n_day'])                                                                   # apply major tick labels to the figure
    ax_to_update.set_xticks(ticks_minor['n_day'], minor = True)                                                                   # apply major tick labels to the figure

    if include_tick_labels:
        ax_to_update.set_xticklabels(ticks_major['yyyymmdd'], rotation = xtick_label_angle, ha = 'left', size = 8)            # update tick labels, and rotate
        #plt.subplots_adjust(bottom=0.15)
        #ax_to_update.set_xlabel('Date')
    else:
        ax_to_update.set_xticklabels([])                                                                    # remove any tick lables if they aren't to be used.  
    
    # add vertical lines every year.  
    for major_tick_n, datetime_majortick in enumerate(ticks_major['datetimes']):
        if datetime_majortick.month == 1:                                                                       # if it's the january tick (i.e the 1st of the year)
            ax_to_update.axvline(x = ticks_major['n_day'][major_tick_n], color='k', alpha=0.1, linestyle='--')                          
             



#%% make_colormap() and remappedColorMap()
         

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

#%% truncate_colormap()

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """ Take a colorbar and crop it.  Useful for removing blue parts of "terrain""
    """
    import matplotlib.colors as colors
    import numpy as np
    
    new_cmap = colors.LinearSegmentedColormap.from_list(
    'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
    cmap(np.linspace(minval, maxval, n)))
    return new_cmap 


#%% status_fig_one_volc()


def status_fig_one_volc(name_and_index, licsalert_status, day_list, fig_type = 'png',
                        volc_dirs = None, out_dir = None, remap = True):
    """ Given a volcano name and its index (col number in the licsalert matrix),
    plot all times for that volcano.  
    Inputs:
        name_and_index | tuple | name and index 
        licsalert_status | r3 array | n_times x n_volcs x 2
        day_list | list of datetimes | date for each row in licsalert matrix. 
    Returns:
        Figure
    History:
        2023_09_04 | MEG | Written
    """
    
    import numpy as np
    from datetime import datetime as dt
    from pathlib import Path
    from glob import glob
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import matplotlib.gridspec as gridspec
    from adjustText import adjust_text
    
    from licsalert.temporal import day_list_to_baselines
    from licsalert.status import remove_dates_with_no_status
    
    # index the 3d one for volcanoes x times x metric to just times x metric for one volcano
    volcano_status = licsalert_status[:, name_and_index[1], :]                          
    # most days don't have an acquisition.  Remove them.  
    volcano_status_crop, day_list_crop = remove_dates_with_no_status(volcano_status, day_list)      
    tbaselines = day_list_to_baselines(day_list_crop)                                               # get the temporal baselines in days
    
    if remap:
        volcano_status_crop_remap = np.zeros(volcano_status_crop.shape)                                                                 # initalise
        for row_n, row in enumerate(volcano_status_crop):                                                                               # loop through each valid time
            volcano_status_crop_remap[row_n, 0] = remap_sigmas(row[0], name_and_index[0], start = (10,10), end = (100, 11))             # do the rescaling
            volcano_status_crop_remap[row_n, 1] = remap_sigmas(row[1], name_and_index[0], start = (10,10), end = (100, 11))             # ditto
        volcano_status_crop = volcano_status_crop_remap
            
    # one static plot
    if fig_type == 'png':
        f, ax = plt.subplots()
        points = ax.scatter(volcano_status_crop[:,0], volcano_status_crop[:,1], c = tbaselines)
        ax.set_xlabel('Sigma for exising def.')
        ax.set_ylabel('Sigma for new def.')
        ax.set_ylim([0,11])
        ax.set_xlim([0,11])
        
        # colorbar which works in day numbers but displays dates.  
        cbar = f.colorbar(points)
        original_tick_labels = cbar.get_ticks()
        new_tick_labels = []
        for tick_label in original_tick_labels:
            try:
                index = list(tbaselines).index(tick_label)                                      # these should be equivalent
                new_tick_labels.append(dt.strftime(day_list[index], '%Y_%m_%d'))                # 
            except:
                new_tick_labels.append('')
        cbar.set_ticklabels(new_tick_labels)
        cbar.set_label(f"Acquisition date")
        
        
        # label the most extreme points with their dates
        label_ratio = 0.2
        texts = []
        x_threshold = (1 - label_ratio) * np.max(volcano_status_crop[:,0])
        y_threshold = (1 - label_ratio) * np.max(volcano_status_crop[:,1])
        for point_n, status in enumerate(volcano_status_crop):
            if (status[0] > x_threshold) or (status[1] > y_threshold):
                texts.append(ax.annotate(dt.strftime(day_list_crop[point_n], '%Y%m%d'), (status[0], status[1])))
        adjust_text(texts, only_move={'points':'y', 'texts':'y'})                                                                                   # adjust how volcanoes of interest are labelled so that the labels don't overlap
        # save the figure
        if not out_dir.is_dir():
            try:
                out_dir.mkdir(parents=True)
            except Exception as e:
                print(f"Error creating directory for the pngs: {e}")
        f.savefig(out_dir / f"{name_and_index[0]}.png", bbox_inches='tight')            
                
        
    # plots for each time step that can be converted into a gif.  
    elif fig_type == 'gif':
        plt.switch_backend('Agg')                                                           # works when there is no X11 forwarding, and when displaying plots during creation would be annoying.  
        volc_dir = volc_dirs[name_and_index[1]]                                             # get the directory of the volcano we're working with.  
        for date_n, status in enumerate(volcano_status_crop):
        
            
            f = plt.figure(figsize = (20,8))
            grid = gridspec.GridSpec(2, 6, wspace = 0.0, hspace = 0.0)                        # divide into 2 sections, 1/5 for ifgs and 4/5 for components
            ax_2d = plt.Subplot(f, grid[-2,0])                                                                                      # create an axes for the IC (spatial source)
            ax_im = plt.Subplot(f, grid[:,1:])                                                                                      # create an axes for the IC (spatial source)
            ax_2d.set_ylim([0,11])
            ax_2d.set_xlim([0,11])
            ax_2d.set_xlabel('Sigma for exising def.')
            ax_2d.set_ylabel('Sigma for new def.')
        
            # ax_2d.scatter(status[0], status[1], c = tbaselines[date_n])
            ax_2d.scatter(volcano_status_crop[:date_n+1,0], volcano_status_crop[:date_n+1,1], c = tbaselines[:date_n+1])                # python excludes end when indexing so +1 to include
            ax_2d.set_aspect('equal')
            
            # Plot the corresponding licsalert figure
            date_string = f"{dt.strftime(day_list_crop[date_n], '%Y%m%d')}"
            date_dir = Path(volc_dir) / date_string
            licsalert_fig_path = glob(str(date_dir/ 'LiCSAlert_figure_with_*_monitoring_interferograms.png'))
            licsalert_fig = mpimg.imread(licsalert_fig_path[0])
            ax_im.imshow(licsalert_fig)
            ax_im.axis('off')
            
            # tidy up axes, save fig etc.  
            for ax in [ax_2d, ax_im]:
                f.add_subplot(ax)
            f.suptitle(f"{name_and_index[0]} : {date_string}")

            # save the figure
            if not out_dir.is_dir():
                try:
                    out_dir.mkdir(parents=True)
                except Exception as e:
                    print(f"Error creating directory for the pngs: {e}")
            f.savefig(out_dir / f"{date_string}.png", bbox_inches='tight')            
            print(f"Saved figure {date_n} of {volcano_status_crop.shape[0]} ")
    if plt.get_backend() != 'Qt5Agg':                                                           # check if we need to reset the backend (to interactive window)
        plt.switch_backend('Qt5Agg')   
        
        
        
#%% status_fig_all_volcs()


def status_fig_all_volcs(licsalert_status, volc_names, day_list, out_dir, figsize, 
                         xlim = 11, ylim = 11, plot_frequency = 6, label_fs = 12):
    """
     This is a 2D plot, with the new deformation unrest metric on one axis, 
     and the existing deformation unrest metric on the other.
    
    """
    
    import numpy as np
    from copy import deepcopy
    import matplotlib.pyplot as plt
    from datetime import datetime as dt
    from adjustText import adjust_text
    plt.switch_backend('Agg')                                                           #  
    #plt.switch_backend('QtAgg')                                                           #   interactive, used to debug.                    
    
    from licsalert.temporal import day_list_to_baselines
    
    #n_frames = len(day_list)
    
    tbaselines = day_list_to_baselines(day_list)                                               # get the temporal baselines in days
    n_volcs = licsalert_status.shape[1]
    
    #licsalert_status_current = np.zeros((licsalert_status.shape[1], 2))                     # first dim is volcano number, 2nd is 
    licsalert_status_current  = deepcopy(licsalert_status[0,])
    day_counter = np.repeat(tbaselines[0], n_volcs)

    volcs_off_chart = []    
    volcs_off_chart_sigmas = []
    
    for day_n, day in enumerate(day_list):                                                                      # loop through all possible days
        for volc_n in range(n_volcs):                                                                              # on each day, loop through all volcanoes     
            existing_def = licsalert_status[day_n, volc_n, 0]                                                   # get the number of sigmas for existing deformation
            new_def = licsalert_status[day_n, volc_n, 1]                                                    # and new deformation
            if  (existing_def != 0) and (new_def != 0):                                                     # check if we have values (and not a date in which they're just zeros) 
                existing_def = remap_sigmas(existing_def, volc_names[volc_n], start = (10,10), end = (100, 11))
                new_def = remap_sigmas(new_def, volc_names[volc_n], start = (10,10), end = (100, 11))
                licsalert_status_current[volc_n,0]  = existing_def                                              # update the current status
                licsalert_status_current[volc_n,1]  = new_def                                                   # continued
                day_counter[volc_n] = tbaselines[day_n]                                                             # also update the counter that records which day number the data are from.  
            
        
        # start the figure.  
        if day_n % plot_frequency == 0:                                                                             # only make hte figure depending on frequency (% is remainder operator)
        
            f, ax = plt.subplots(1, figsize = figsize)                                                                  # one figure per day.  
            points = ax.scatter(licsalert_status_current[:, 0], licsalert_status_current[:, 1], 
                                c = day_counter, vmin = int(tbaselines[0]), vmax = int(tbaselines[-1]))
            
            # deal with the colorbar
            cbar = f.colorbar(points)
            original_tick_labels = cbar.get_ticks()
            new_tick_labels = []
            for tick_label in original_tick_labels:
                try:
                    index = list(tbaselines).index(tick_label)                                      # these should be equivalent
                    new_tick_labels.append(dt.strftime(day_list[index], '%Y_%m_%d'))                # 
                except:
                    new_tick_labels.append('')
            cbar.set_ticklabels(new_tick_labels)
            cbar.set_label(f"Acquisition date")
            
            # label volcano names if above a threshold
            texts = []
            label_ratio = 0.2                                                                                                                           # the top ratio are labelled.  
            x_threshold = (1 - label_ratio) * np.nanmax(licsalert_status_current[:,0])
            y_threshold = (1 - label_ratio) * np.nanmax(licsalert_status_current[:,1])
            for volc_n, status in enumerate(licsalert_status_current):                                                                                  # loop through every volcano
                if (status[0] > x_threshold) or (status[1] > y_threshold):                                                                              # if one measure is above threshold    
                    texts.append(ax.annotate(volc_names[volc_n], (licsalert_status_current[volc_n, 0], licsalert_status_current[volc_n, 1]),
                                 fontsize = label_fs))           # label the point.  
            
                    
            # tidy up the figure
            ax.set_xlabel('Sigma for exising def.')
            ax.set_ylabel('Sigma for new def.')
            title = dt.strftime(day, '%Y_%m_%d')
            ax.set_title(title)
    
            ax.set_ylim([0,ylim])                                                                                                   # axes size and ticks and labels.  
            ax.set_xlim([0,xlim])
            ticks = np.concatenate((np.arange(0,11, 2), np.array([11])))
            tick_labels = np.concatenate((np.arange(0,11, 2), np.array([100])))
            ax.set_xticks(ticks)
            ax.set_xticklabels(tick_labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(tick_labels)
            plt.grid(alpha = 0.2)
            
            # to be removed as now remap large values
            # if len(volcs_off_chart) > 0:                                                        # if there are volcanoes off the chart here
            #     off_volcs = f"Volcanoes off chart:\n"
            #     for volcs_off_chart_sigma, volc_off_chart in zip(volcs_off_chart_sigmas, volcs_off_chart):
            #         off_volcs = off_volcs + f"{volc_off_chart} (x:{volcs_off_chart_sigma[0]:.2f} y:{volcs_off_chart_sigma[1]:.2f})\n"
            #     ax.annotate(off_volcs, (xlim - 0.02, ylim - 0.1), horizontalalignment='right',verticalalignment='top',)
            
            adjust_text(texts, only_move={'points':'y', 'texts':'y'})                                                                                   # adjust how volcanoes of interest are labelled so that the labels don't overlap
            
            # save the figure
            if not out_dir.is_dir():
                try:
                    out_dir.mkdir(parents=True)
                except Exception as e:
                    print(f"Error creating directory for the pngs: {e}")
            f.savefig(out_dir / f"{title}.png", bbox_inches='tight')            
            print(f"Saved figure {title}")

    if plt.get_backend() != 'Qt5Agg':                                                           # check if we need to reset the backend (to interactive window)
        plt.switch_backend('Qt5Agg')                                                           #  


#%% remap_sigmas()


def remap_sigmas(signal, volc_name, start = (10,10), end = (100, 11)):
    """
    Signal below start[0] are not changed are not changed.  
    Signals in the range set by start [0] and end [0] are remapped into range 
    start[1] and end [1].  
    Above end[0], they are set to end[1].  
    to startValues are considered with input (signal) on the x, and output (signal_remap) on the y.  
    Inputs:
        signal | rank 1 | signal 
        volc_name | string | sigma signal is from this volcano.  
        start | tuple | (signal, signal_remap) of start of remapping.  
        end | tuple | (signal, signal_remap) of end of remapping.  
    Returns:
        signal_remap | rank 1 | signal in new range.  
    """
    import numpy as np
    
    if np.isnan(signal):
        
        return signal
    else:
        threshold = start[0]
        
        if signal <= start[0]:                                                              # if below range, don't change.                                   
            signal_remap = signal
        elif (start[0] < signal) and (signal <= end[0]):                                    # if in range, remap.  
            gradient = (end[1] - start[1]) / (end[0] - start[0])                            # delta y / delta x
            c = start[1] - gradient * start[0]                                              # rearrange y = mx +c, subbing in start point.  
            signal_remap = (gradient * signal)  + c
        elif end[1] < signal:                                                               # if above, set to maximum
            signal_remap = end[1]
            print(f"Warning: {volc_name} has a sigma value of {signal},"
                  f"which is above the maximum remapping value of {end[1]}. "
                  f"It has been set to {end[1]}.  ")
        else:
            raise Exception(f"Failed to determine which interval the signal lies in.  Exiting.  ")
        return signal_remap