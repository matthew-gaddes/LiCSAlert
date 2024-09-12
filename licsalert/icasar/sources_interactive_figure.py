#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 14:35:36 2024

@author: matthew
"""

import pdb


#%%


def plot_2d_interactive_fig(S_pca, S_hists, mask, spatial, sica_tica, 
                            hdbscan_param, tsne_param, 
                            n_converge_bootstrapping, n_converge_no_bootstrapping,
                            inset_axes_side = {'x':0.1, 'y':0.1}, 
                            arrow_length = 10., figsize = (19,10), 
                            labels = None, legend = None, markers = None, 
                            figures = 'window', png_path = './', 
                            fig_filename = '2d_interactive_plot'):
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

        
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from matplotlib.widgets import Slider, Button
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from licsalert.aux import col_to_ma
    
    
        
    # mutable container to store the scatter plot object
    sc_container = {}
    # mutable list to store the 2d array of xy coords. for 
    xy_tsne = [None]
    # mutable sources container (used to plot a point when hovered over)
    S_all_r3 = [None]
    # mutable source axes (for PCs and ICs)
    source_axes = [None]
    # mutable list to store the 2d array of xy coords. 
    xy_tsne = [None]
    # mutable to store the S_ica (ica sources)
    S_ica = [None]
    # mutable to store other less useful outputs
    source_outputs = {}
    
    
    button_clicked = [False]
    
    # get the number of PCs for each run
    n_pca_comps = [i[0].shape[0] for i in  S_hists]
    
    
    
    # Draw the figure
    fig = plt.figure(figsize = figsize)                                                                
    axes1 = fig.add_axes([0.1, 0.25, 0.35, 0.6])                                                        
    
    # Try and add various labels from the labels dict
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
    
    # Add sliders
    axcolor = 'lightgoldenrodyellow'
    horizontal_slider_dx = 0.3
    # Slider for n_pca_comp
    ax_a = plt.axes([0.15, 0.9, horizontal_slider_dx, 0.03], facecolor=axcolor)
    slider_a = Slider(ax_a, '# PCA components', n_pca_comps[0], n_pca_comps[-1], 
                      valinit=int(np.average(n_pca_comps)), valstep = 1)
    # Slider for tsne perplexity
    ax_b = plt.axes([0.15, 0.1, horizontal_slider_dx, 0.03], facecolor=axcolor)
    slider_b = Slider(ax_b, 'tsne: perplexity', 5, 50, 
                      valinit=tsne_param[0], valstep = 1)
    # Slider for tsne early exageration
    ax_c = plt.axes([0.15, 0.15, horizontal_slider_dx, 0.03], facecolor=axcolor)
    slider_c = Slider(ax_c, 'tsne: early exageration', 0, 20, 
                      valinit=tsne_param[1], valstep = 1)
    # Slider for hdbscan min_cluster_size
    ax_d = plt.axes([0.015, 0.25, 0.0225, 0.63], facecolor=axcolor)
    slider_d = Slider(ax_d, 'HDBSCAN: min_cluster_size', 0, 200, 
                      valinit=hdbscan_param[0], valstep = 1, 
                      orientation='vertical')
    slider_d.label.set_rotation(90)
    slider_d.label.set_verticalalignment('center')
    slider_d.label.set_position((0.0, 0.5))

    # Slider for hdbscan min_samples
    ax_e = plt.axes([0.04, 0.25, 0.0225, 0.63], facecolor=axcolor)
    slider_e = Slider(ax_e, 'HDBACAN: min_samples', 0, 50, 
                      valinit=hdbscan_param[1], valstep = 1,
                      orientation='vertical')
    slider_e.label.set_rotation(90)
    slider_e.label.set_verticalalignment('center')
    slider_e.label.set_position((0.08, 0.5))
    
    slider_axes = [ax_a, ax_b, ax_c, ax_d, ax_e]

    # button to advance when user has configured
    button_dx = 0.1
    ax_button = plt.axes([(0.15 + horizontal_slider_dx/2 + button_dx/2 ), 
                          0.025, button_dx, 0.04])
    button_axes = [ax_button]
              
    # call the function that highlights a source if you hover on it.  
    fig.canvas.mpl_connect("motion_notify_event", lambda val: hover(val, axes1,
                                                                    sc_container,
                                                                    fig,
                                                                    xy_tsne, 
                                                                    arrow_length, 
                                                                    inset_axes_side,
                                                                    sica_tica,
                                                                    S_all_r3,
                                                                    slider_axes,
                                                                    source_axes,
                                                                    button_axes))



    # arguemnts required by update function.  
    params = {'fig': fig,
            'axes1': axes1,
            'slider_a': slider_a,
            'slider_b': slider_b,
            'slider_c': slider_c,
            'slider_d': slider_d,
            'slider_e': slider_e,
            'n_pca_comps': n_pca_comps,
            'S_hists': S_hists,
            'mask': mask,
            'spatial': spatial,
            'sica_tica': sica_tica,
            'n_converge_bootstrapping': n_converge_bootstrapping,
            'n_converge_no_bootstrapping': n_converge_no_bootstrapping,
            'S_all_r3': S_all_r3,
            'xy_tsne': xy_tsne,
            'markers': markers,
            'sc_container': sc_container,
            'legend': legend,
            'source_axes': source_axes,
            'S_pca': S_pca,
            'S_ica' : S_ica,
            'source_outputs' : source_outputs}


    # call the slider function to plot the points during 1st plot initialisation
    update(None, update_tsne = True, **params)

    # Call update function when slider value changes
    # pca slider, so update tsne
    slider_a.on_changed(lambda val: update(val, update_tsne = True, **params))
    # tsne sliders, so update tsne
    slider_b.on_changed(lambda val: update(val, update_tsne = True, **params))
    slider_c.on_changed(lambda val: update(val, update_tsne = True, **params))
    #hdbscan sliders, so don't update tsne
    slider_d.on_changed(lambda val: update(val, update_tsne = False, **params))
    slider_e.on_changed(lambda val: update(val, update_tsne = False, **params))

    # Initialize the user_data variable (mutable?)
    chosen_settings = {}

    # Initialize a flag to indicate when the button has been clicked
    button = Button(ax_button, 'Continue')
        
    # Set the callback function for the button
    button.on_clicked(lambda val: on_button_click(val, slider_a, slider_b,
                                                  slider_c, slider_d,
                                                  slider_e, button_clicked))
    
    # Display the interactive figure
    plt.show()#block=False)


    # Wait for user input, if there is a window (i.e. don't if making only 
    # a png)
    if (figures == 'window') or (figures == 'png+window'):
        while not button_clicked[0]:
            plt.pause(1.)
            print("Waiting for the user to select the ICA sources using the "
                  "interactive window.  ")
        print(f"{S_ica[0].shape[0]} ICA sources have been chosen by the user. "
              "Continuing")
    else:
        print(f"{S_ica[0].shape[0]} ICA sources have been chosen automatically."
          "Continuing")

    # only save the figure if a png is selected.      
    if figures == 'window':
        pass
    elif (figures == "png") or (figures == 'png+window'):
        fig.savefig(f"{png_path}/{fig_filename}.png")
        plt.close()
    else:
        raise Exception(f"'figures' was set incorrectly, and should be  "
                        "'window', 'png', or 'png+window', but is {figures} "
                        "Exiting.")
    
    # return the ICA sources from the mutable they are stored in 
    return S_ica[0], source_outputs
    
    

#%%

def hover(event, axes1, sc_container, fig, xy_tsne, arrow_length, 
          inset_axes_side, sica_tica, S_all_r3, slider_axes, source_axes,
          button_axes):
    """ 
    """
    # determine if the mouse is in the axes
    if event.inaxes == axes1:                                                       
        # cont is a boolean of if hovering on point, ind is a dictionary about 
        # the point being hovered over.  Note that two or more points can be in this.  
        cont, ind = sc_container['sc'].contains(event)                                              
        # if on point
        if cont:                                                                    
            # remove the axes and arrow created when hovering on the point 
            # (in case cursor moves from one point to next without going off a point)
            remove_axes2_and_arrow(fig, axes1, slider_axes, source_axes, button_axes)                                             
            
            # get the index of which data point we're hovering on in a simpler form.      
            point_n = ind['ind'][0]                                                 
            
            # record the limits so we can reset them
            xlim_orig = axes1.get_xlim()
            ylim_orig = axes1.get_ylim()
            
            # 1: Add the annotation arrow (from inset axes to data point)
            # calculate the length of the arrow, which depends on which quadrant 
            # we're in (as the arrow always goes away from the plot)
            arrow_lengths = calculate_insetaxes_offset([axes1.get_xlim(), 
                                                        axes1.get_ylim()], 
                                                      [xy_tsne[0][0,point_n], 
                                                       xy_tsne[0][1,point_n]], 
                                                      arrow_length)                               
            # add the arrow.  Notation is all a bit backward as head is 
            # fixed at end, so it has to be drawn backwards.  
            # clip_on makes sure it's visible, even if it goes off the edge of the axes.  
            axes1.arrow(xy_tsne[0][0,point_n] + arrow_lengths[0], 
                        xy_tsne[0][1,point_n] + arrow_lengths[1],                                       
                        -arrow_lengths[0], -arrow_lengths[1], 
                        clip_on = True, zorder = 999)

            # get the start of the arrow (not on the data point)
            x_data = xy_tsne[0][0,point_n] + arrow_lengths[0]
            y_data = xy_tsne[0][1,point_n] + arrow_lengths[1]
            
            # and conver to position in figure in ragne [0 ,1]
            x_fig, y_fig = axes1.transData.transform((x_data, y_data))
            fig_width, fig_height = fig.get_size_inches() * fig.dpi
            x_normalized = x_fig / fig_width
            y_normalized = y_fig / fig_height

            # 2: Add the inset axes                
            # four possible quadrants for the inset axes.  Note they are 
            # shifted so they don't lie on top of the arrow.  
            if arrow_lengths[0] > 0 and arrow_lengths[1] > 0:                                                          
                inset_axes = fig.add_axes([x_normalized, y_normalized,                                                               
                                            inset_axes_side['x'], 
                                            inset_axes_side['y']], anchor = 'SW')               
            elif arrow_lengths[0] < 0 and arrow_lengths[1] > 0:                                                        
                inset_axes = fig.add_axes([x_normalized - inset_axes_side['x'], 
                                           y_normalized,                                        
                                            inset_axes_side['x'], 
                                            inset_axes_side['y']], anchor = 'SE')     
            elif arrow_lengths[0] > 0 and arrow_lengths[1] < 0:                                                        
                inset_axes = fig.add_axes([x_normalized, 
                                           y_normalized - inset_axes_side['y'],                                        
                                            inset_axes_side['x'], 
                                            inset_axes_side['y']], anchor = 'NW')                 
            else:                                                                                                      
                inset_axes = fig.add_axes([x_normalized - inset_axes_side['x'],
                                           y_normalized - inset_axes_side['y'],                 
                                            inset_axes_side['x'], 
                                            inset_axes_side['y']], anchor = 'NE')                
                
            # 3: Plot data on the inset axes
            if sica_tica == 'sica':
                inset_axes.matshow(S_all_r3[0][point_n,])
            else:
                raise Exception("Not tested yet.  ")
                #inset_axes.plot(temporal_data['xvals'], temporal_data['tcs_r2'][point_n,])                            # draw the inset axes time course graph
                inset_axes.axhline(0)
                
            inset_axes.set_xticks([])                                                                                       # and remove ticks (and so labels too) from x
            inset_axes.set_yticks([])                                                                                       # and from y
            
            # and reset the limits
            axes1.set_xlim(xlim_orig)
            axes1.set_ylim(ylim_orig)

            
            fig.canvas.draw_idle()                                                                                          
        else:                                                                       # else not on a point
            remove_axes2_and_arrow(fig, axes1, slider_axes, source_axes, button_axes)                                             # remove the axes and arrow created when hovering on the point                       
    else:                                                                           # else not in the axes
        remove_axes2_and_arrow(fig, axes1, slider_axes, source_axes, button_axes)                                                 # remove the axes and arrow created when hovering on the point (in case cursor moves from one point to next without going off a point)



#%%

def update(val, fig, axes1, slider_a, slider_b, slider_c, slider_d, slider_e, 
           n_pca_comps, S_hists, mask, spatial, sica_tica, 
           n_converge_bootstrapping, n_converge_no_bootstrapping, 
           S_all_r3, xy_tsne, markers, sc_container, legend, source_axes, 
           S_pca, S_ica, source_outputs,  update_tsne = True):
    
    """
    """
    
    import numpy as np
    import matplotlib.pyplot as plt 
    
    from licsalert.icasar.icasar_funcs import tsne_and_cluster
    from licsalert.icasar.icasar_funcs import sources_list_to_r2_r3
    
    # remove points from previous run
    axes1.clear()
    
    # get the slider values
    n_pca = slider_a.val
    hdbscan_param = (slider_d.val, slider_e.val)
    tsne_param = (slider_b.val, slider_c.val)
        
    # determine with multiple ICA runs we're using. 
    S_index = n_pca_comps.index(n_pca)
    outputs = tsne_and_cluster(S_hists[S_index], mask, spatial, sica_tica, 
                               n_pca, hdbscan_param, tsne_param,
                               n_converge_bootstrapping, 
                               n_converge_no_bootstrapping, 
                               update_tsne)
    (sources_all_r3, S_ica[0], labels_hdbscan, xy_tsne_new, marker_dict, 
     legend_dict, labels_colours, Iq_sorted, n_clusters) = outputs
    S_all_r3[0] = sources_all_r3
    # if tsne was not updated, previous function returns None
    if xy_tsne_new is None:
        pass
    else:
        # update the mutable list
        xy_tsne[0] = xy_tsne_new.T

    # do the plotting
    if markers is None:      
        # either plot all the points
        sc_container['sc'] = axes1.scatter(xy_tsne[0][0,:], xy_tsne[0][1,:], 
                                           c=labels_colours, s=100)                                             
    else:                                                                                                                                     
        # or plot them based on the different markers 
        n_markers = len(markers['styles'])                                                                                                    
        for n_marker in range(n_markers):                                                                                                     
            point_args = np.ravel(np.argwhere(markers['labels'] == n_marker))                                                                 
            try:
                sc_container['sc'] = axes1.scatter(xy_tsne[0][0,:][point_args], 
                                   xy_tsne[0][1,:][point_args], 
                                   c=labels_colours[point_args], s=100, 
                                   marker = markers['styles'][n_marker])        
            except:
                pass
        # draw the scatter plot again with all the points (regardless of 
        # marker style), but with invisible markers.  As the last to be drawn, 
        # these are the ones that are hovered over, and indexing works as 
        # all the points are draw this time.  
        sc_container['sc'] = axes1.scatter(xy_tsne[0][0,:], xy_tsne[0][1,:], 
                                           c=labels_colours, s=100,  
                                           alpha = 0.0)                                                                       
    axes1.set_aspect('equal')
    
    # add things to the mutable
    source_outputs['xy'] = xy_tsne[0]
    source_outputs['labels'] = labels_hdbscan
    source_outputs['Iq_sorted'] = Iq_sorted
    source_outputs['n_clusters'] = n_clusters
    # get all the sources for the 2d points (depends on n_pca) as row vectors
    sources_all_r2, _ = sources_list_to_r2_r3(S_hists[S_index], mask)           
    source_outputs['sources_all_r2'] = sources_all_r2


    # Remove any previous PCs and ICs (ready to plot updated ones)
    for ax in source_axes:
        try:
            ax.remove()
        except:
            pass
        
    # Plot the PCs
    plot_signals(fig, source_axes, S_pca, mask, 0.5, 0.95, 0.5, 1., title = 'PC', 
                 legend_dict = None, obscure = n_pca)
    
    # Plot the ICs
    plot_signals(fig, source_axes, S_ica[0], mask, 0.5, 0.95, 0., 0.5, title = 'IC',
                 legend_dict = legend_dict, obscure = None)
    fig.canvas.draw_idle()
    

#%%

def on_button_click(event, slider_a, slider_b, slider_c, slider_d, slider_e,
                    button_clicked):
    """
    """
    
    import matplotlib.pyplot as plt

    # assign the current settings to the mutable in the parent function.  
    chosen_settings = {
        'n_pca': slider_a.val,
        'hdbscan_param' : (slider_d.val, slider_e.val),
        'tsne_param' : (slider_b.val, slider_c.val)
    }
    
    # Set the flag to True
    button_clicked[0] = True           

    plt.close()                     
    
    
#%%

def plot_signals(fig, source_axes, r2_signals, mask, x_start, x_stop, y_start, y_stop,
                 title, legend_dict = None, obscure = None):
    """
    """
    from matplotlib.lines import Line2D                                  
    
    from licsalert.aux import col_to_ma
    
    n = r2_signals.shape[0]
    
    if n <= 0:
        raise ValueError("Number of subplots must be greater than 0.")
    if not (0 <= x_start <= 1) or not (0 <= x_stop <= 1):
        raise ValueError("x_start and x_stop must be between 0 and 1.")
    if not (0 <= y_start <= 1) or not (0 <= y_stop <= 1):
        raise ValueError("y_start and y_stop must be between 0 and 1.")
    if x_start >= x_stop:
        raise ValueError("x_stop must be greater than x_start.")
    if y_start >= y_stop:
        raise ValueError("y_stop must be greater than y_start.")
    
    fig_width = x_stop - x_start  
    fig_height = y_stop - y_start  

    # Number of columns for each row
    ncols = (n + 1) // 2

    # Width and height of each subplot
    width = fig_width / ncols
    height = fig_height / 2

    for i in range(n):
        # determine if it should have low opacity:
        opacity = 1.
        if obscure is not None:
            if i >= obscure:
                opacity = 0.5
        
        row = i // ncols
        col = i % ncols
        left = x_start + col * width
        bottom = y_start + (1 - (row + 1) / 2) * fig_height  

        ax = fig.add_axes([left, bottom, width, height])
        if ax not in source_axes:
            source_axes.append(ax)
            
        #ax.set_title(f'Subplot {i+1}')
        mat = ax.matshow(col_to_ma(r2_signals[i, :], mask),
                         alpha = opacity)
        ax.set_xticks([])
        ax.set_yticks([])
        
        ax.text(0.05, 0.95, f"   {title}{i}", transform=ax.transAxes, 
                fontsize=12, verticalalignment='top', 
                bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', 
                          facecolor='white'))
        
        
        # Possibly add dot showing which cluster source came from
        if legend_dict is not None:
            cluster_col = legend_dict['elements'][i].get_markerfacecolor()
            line = Line2D([4], [3], marker='o', color=cluster_col, 
                          markersize=10, zorder = 999)
            ax.add_line(line)


#%%

def remove_axes2_and_arrow(fig, axes1, slider_axes, source_axes, button_axes):
    """ Given a figure that has a second axes and an annotation arrow due to a 
    point having been hovered on, remove this axes and annotation arrow.  
    Inputs:
        fig | matplotlib figure 
    Returns:
    History:
        2020/09/08 | MEG | Written
    """
    
    import matplotlib
    
    # 1: try and remove any axes except the primary one and the slider ones
    try:
        for ax_n, ax in enumerate(fig.axes):
            if ( (ax_n > 0) and 
                (ax not in slider_axes) and
                (ax not in source_axes) and 
                (ax not in button_axes)):
                ax.remove()                
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
    fig.canvas.draw_idle()                                          


#%%

def calculate_insetaxes_offset(lims, points, offset_length):
    """
    The offsets between the inset axes and the point are different depending 
    on which quadrant of the graph the point is in.  
    Inputs:
        lims | list | length is equal to the number of dimensions.  Filled with tuples of the axes limits.  
        point | list | length is equal to the number of diemsions. Filled with points.  
        offset_length | float | length of the arrow.  
    Returns:
        offsets | list | length is equal to the number of dimensions.  
                        Length of offset for inset axes in each dimension.  
    History:
        2020/09/08 | MEG | Written
    """
    import numpy as np
    offsets = []
    for dim_n in range(len(lims)):                                        # loop through each dimension.  
        dim_centre = np.mean(lims[dim_n])
        # determine which side of centre point is.  
        if points[dim_n] < dim_centre:
            offsets.append(-offset_length)
        else:
            offsets.append(offset_length)
    return offsets