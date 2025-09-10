#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 13:41:27 2025

@author: matthew
"""

import pdb

#%% automatic_pixel_epoch_selection()

def automatic_pixel_epoch_selection(
        displacement_r3, 
        tbaseline_info,
        baseline_end,
        volcano_dir,
        figures,
        interactive=False
        ):
    """
    Given a time series with a time varying mask (i.e. pixels come 
    in and out of coherene), build a time series with a consistent mask
    that uses only some of these acquisitions to build a compromise 
    between temporal resolution and number of pixels
    
    
    Inputs
    
    spatial_ICASAR_data = {'ifgs_dc'       : displacement_r2['mixtures_mc'][:(baseline_end.acq_n+1),],                             
                           'mask'          : displacement_r2['mask'],
                           'lons'          : displacement_r2['lons'],
                           'lats'          : displacement_r2['lats'],
                           
                           
   'ifg_dates_dc'  : tbaseline_info['ifg_dates'][:(baseline_end.acq_n+1)]}                             
    
    volcano_dir | Path | outdir for figures.  
    
    """
    
    import numpy as np
    import numpy.ma as ma
    
    from licsalert.pixel_selection import calculate_optimal_n_epochs
    from licsalert.pixel_selection import consistent_pixels_plot
    from licsalert.aux import r3_to_r2


    # crop the input data in time
    cum_ma_baseline = displacement_r3['cum_ma'].original[:(baseline_end.acq_n+1)]
    acq_dates_baseline = tbaseline_info['acq_dates'][:(baseline_end.acq_n+1)]                                 
    
    # determine the number of pixels for each epoch
    n_pixels, n_pixels_idx, total_pix = calculate_valid_pixels(
        cum_ma_baseline
        )
    
    # and how those change as we add epochs
    n_pix_epoch = intersect_valid_pixels(
        cum_ma_baseline,
        acq_dates_baseline,
        verbose = False
        )
    
    # calculate optimal number of epochs
    epoch_values, optimal_epoch_n = calculate_optimal_n_epochs(
        n_pix_epoch,
        volcano_dir,
        figures,
        )
    # tidy as not needed.  
    del epoch_values
    
    
    # Figure, which can be set to interactive for debugging.  
    consistent_pixels_plot(
        cum_ma_baseline,
        acq_dates_baseline,
        volcano_dir,
        figures,
        optimal_epoch_n,
        interactive,
        )
        
    
    # build the time series that uses this number of epochs
    selected_epochs = sorted(list(n_pixels_idx[:optimal_epoch_n]))    
    cum_ma_ica = cum_ma_baseline[selected_epochs, ]
    acq_dates_ica = [acq_dates_baseline[i] for i in selected_epochs]
    
    # convert to ICA data form
    # debug            
    # import matplotlib.pyplot as plt
    # for im in cum_ma_ica:
    #     f, ax = plt.subplots(1)
    #     ax.matshow(im)
    
    mask_bool = ma.getmaskarray(cum_ma_ica)   # nomask → array(False)
    # 2. Logical OR along the time axis:  True if masked at least once
    mask_2d = np.any(mask_bool, axis=0)
    
    # f, ax = plt.subplots()
    # ax.matshow(mask_2d)
    
    # make the mask that's consistent in time.  
    mask_3d = np.repeat(
        mask_2d[np.newaxis,],
        cum_ma_ica.shape[0],
        axis = 0
        )
    
    # and apply to the data.  
    cum_ma_ica_consistent = cum_ma_ica.copy()
    cum_ma_ica_consistent.mask = mask_3d
    
    # debug plot
    # for im in cum_ma_ica_consistent:
    #     f, ax = plt.subplots()
    #     ax.matshow(im)
    #     plt.pause(0.5)
        
    # flatten to ICA standard (image is a row vector)            
    displacement_r2_ica = r3_to_r2(cum_ma_ica_consistent)

    tbaseline_info_ica = {
        'acq_dates' : acq_dates_ica
        }

    return displacement_r2_ica, tbaseline_info_ica



#%% consistent_pixels_plot()

def consistent_pixels_plot(
        cum_r3,
        acq_dates,
        volcano_dir,
        figures,
        optimal_epoch_n,
        interactive=False,
        ):
    """
    
    Allows a user to select the number of epochs that are included, and shows
    the resulting set of consistent pixels.  The epochs are sorted by the number
    of pixels, so that the first preserves the most pixels etc.  
    
    Left shows Mask, right shows masked cumualtive displacement,
    lower left shows # pixels for each epoch, 
    lower right shows number of consistent pixels for each number of epochs added.  
    
    
    Given:
      - cum_r3 (numpy MaskedArray or boolean 3D array) of shape (T, H, W), where 
        cum_r3[t, :, :] indicates whether each pixel is valid (False mask or True).
      - acq_dates (list of length T), each entry is a yyyymmdd string representing that image's date.
    
    Returns:
      selected_times (list): The sorted order (by valid pixels) of epochs.
      final_consistent_pixels (2D array): The final intersection (boolean) of selected cum_r3.
    """
    import numpy as np
    import numpy.ma as ma
    import matplotlib.pyplot as plt
    from datetime import datetime as dt
    import matplotlib.ticker as ticker
    
    from licsalert.monitoring_functions import calculate_valid_pixels
    from licsalert.monitoring_functions import intersect_valid_pixels
    
    
            
    def draw_plots(num_epochs):
        """
        """
        print(f"Updating the figure to {num_epochs} epochs...", end = '')
        
        
        # get the cum_r3 up to the epoch number selected
        cum_r3_subset = cum_r3[sorted(n_pixels_idx[:num_epochs]), ]
        
        # get the mask of pixels valid at all times
        consistent_pixels = np.logical_or.reduce(cum_r3_subset.mask, axis=0)
        
        # determine which epochs (acquisition dates) are currently selected
        acq_dates_selected = [acq_dates[i] for i in n_pixels_idx[:num_epochs]]
        acq_dates_deselected = [
            i for i in acq_dates if i not in acq_dates_selected
            ]
        
        # Update the subplots:
        ax_cum.clear()
        ax_mask.clear()
        ax_dates.clear()
        ax_pixs.clear()
        cax.clear()
        

        # Top Left: Show the last added image masked by the current consistent mask
        # (This visualizes displacement using the last epoch's data)
        # get the cumulative displacment in the time series that we have
        im = ax_cum.matshow(
            ma.array(cum_r3_subset[-1,], mask=consistent_pixels) - 
            ma.array(cum_r3_subset[0,], mask=consistent_pixels)
            )
        ax_cum.set_title('Masked cumulative displacement')
        ax_cum.tick_params(
            axis='x', bottom=True, top=False, labelbottom=True, labeltop=False
            )
        # and the colorbar for the cumulative displacement
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("Displacement (mm)")
        

        # Top Right: Show the consistent pixel map.
        # (Here we invert consistent_pixels so that True (valid) appears white)
        ax_mask.matshow(~consistent_pixels, cmap='gray', interpolation='nearest')
        ax_mask.set_title(f"Mask after adding epoch {n_pixels_idx[num_epochs-1]}")
        ax_mask.tick_params(
            axis='x', bottom=True, top=False, labelbottom=True, labeltop=False
            )

        # Bottom Left: Plot the dates corresponding to the selected epochs.
        acq_date_status = [1 if acq_date in acq_dates_selected else 0 for acq_date in acq_dates]
        
        ax_dates.scatter([dt.strptime(i, '%Y%m%d') for i in acq_dates], n_pixels, c = acq_date_status)
        ax_dates.set_ylim([0, total_pix])
        ax_dates.yaxis.set_major_formatter(CustomSciFormatter())
        ax_dates.set_title("Selected Dates So Far")
        ax_dates.set_xlim(
            left=dt.strptime(acq_dates[0], '%Y%m%d'),
            right=dt.strptime(acq_dates[-1], '%Y%m%d')
            )
        ax_dates.tick_params(axis='x', labelrotation=90)
        ax_dates.set_ylabel('# pixels')


        # Bottom Right: Plot the valid pixel count (number of unmasked pixels) 
        # at each step.
        ax_pixs.plot(np.arange(0, len(n_pix_epoch)), n_pix_epoch)
        # plot the dot which shows how many epochs are selected
        ax_pixs.scatter(num_epochs, n_pix_epoch[num_epochs], c = 'r')
        ax_pixs.set_xlim(left=1, right=len(acq_dates))
        ax_pixs.set_ylim(bottom=0, top=total_pix)
        ax_pixs.yaxis.set_major_formatter(CustomSciFormatter())
        ax_pixs.set_title("Valid Pixel Count")
        ax_pixs.set_xlabel('# epochs selected')

        fig.canvas.draw_idle()
        print(f" Done.")
    
    
    class CustomSciFormatter(ticker.Formatter):
        def __call__(self, x, pos=None):
            # Format x in scientific notation with no decimals
            s = f"{x:.0e}"
            # Remove the plus sign from the exponent (e.g., turn "3e+05" into "3e5")
            s = s.replace("e+0", "e").replace("e+", "e")
            return s


    # determine the number of pixels for each epoch
    n_pixels, n_pixels_idx, total_pix = calculate_valid_pixels(cum_r3)
    
    # and how those change as we add epochs
    n_pix_epoch = intersect_valid_pixels(cum_r3, acq_dates, verbose = False)

    # Set up the figure and axes
    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(nrows=2, ncols=2, height_ratios=[4, 1])
    ax_cum = fig.add_subplot(gs[0, 1])
    ax_mask = fig.add_subplot(gs[0, 0])
    ax_dates = fig.add_subplot(gs[1, 0])
    ax_pixs = fig.add_subplot(gs[1, 1])
    
    cax = fig.add_axes([0.905, 0.45, 0.01, 0.3])  
    
    # Ensure date lables aren't lost
    plt.subplots_adjust(bottom=0.15)
    
    # either interactive to choose the number of epochs
    if interactive:
        draw_plots(num_epochs = 1)
        def on_click(event):
            if event.inaxes == ax_pixs:
                draw_plots(num_epochs = int(event.xdata))
                
        fig.canvas.mpl_connect("button_press_event", on_click)
        plt.show()
        # keep interactive
        plt.pause(999)
                
    # or draw figure with number of epochs passed as an argument.  
    else:
        draw_plots(optimal_epoch_n)
        
        #and then save or display
        outpath=volcano_dir/'baseline_pixel_selection'
        outpath.mkdir(parents=True, exist_ok=True)
        if figures == 'window':
            pass
        elif figures == "png":
            try:
                fig.savefig(
                    outpath/"baseline_pixels_vs_epochs.png",
                    bbox_inches='tight'
                    )
                plt.close(fig)
            except:
                print(f"Failed to save the figure.  Trying to continue.  ")
        elif figures == 'both':
            try:
                fig.savefig(
                    outpath/"baseline_pixels_vs_epochs.png",
                    bbox_inches='tight'
                    )
            except:
                print(f"Failed to save the figure.  Trying to continue.  ")
    
    

#%% calculate_valid_pixels()

def calculate_valid_pixels(cum_r3):
    """
    
    returns:
        n_pixels | number of pixels on each epoch.  
        n_pixels_idx | | index of epochs when starting with the largest number
                        of pixels first (i.e. index of best ifg, 2nd best ifg etc.)
        total_pix | int | maximum number of pixels in image (includes water)
    
    """
    
    import numpy as np

    # Compute number of valid pixels in each epoch
    # (mask is True where pixel is masked)
    n_pixels = [np.sum(~cum_r2.mask) for cum_r2 in cum_r3]

    # debug plot    
    # f, ax = plt.subplots()
    # ax.scatter([dt.strptime(i, '%Y%m%d') for i in acq_dates], n_pixels)

    # Sort epochs in descending order of valid pixel count 
    # i.e. epoch number with most pixels first
    sorted_pairs = sorted(enumerate(n_pixels), key=lambda x: x[1], reverse=True)
    n_pixels_idx, _ = zip(*sorted_pairs)

    # Total pixel count in a single image
    total_pix = cum_r3.shape[1] * cum_r3.shape[2]
    
    return n_pixels, n_pixels_idx, total_pix


#%% intersect_valid_pixels()

def intersect_valid_pixels(cum_r3, acq_dates, verbose = False):
    """ Calculate the number of pixels that are valid at all epochs for
    the time series as the number of epochs is increased in order of the 
    number of pixels on each date (i.e. add dates with lots of pixels first)
    """
    import numpy as np
    import numpy.ma as ma

    n_epochs, ny, nx = cum_r3.shape

    # get the number of pixels per epoch, and the order
    n_pixels, n_pixels_idx, total_pix = calculate_valid_pixels(cum_r3)
    
    # initialise
    n_pix_epoch = []
    cum_r3_resampled = ma.zeros((0, ny, nx))
    
    # nothing is masked before we start
    mask_r2 = np.zeros((0, ny, nx), dtype=bool)
    
    for progress, epoch_n in enumerate(n_pixels_idx):
        if verbose:
            print(
                f"Step {progress} of {n_epochs}: Adding epoch "
                f"{acq_dates[epoch_n]} to the time series...", end = ''
                )
        
        # get displacement to that epoch
        cum_r2 = cum_r3[epoch_n:epoch_n+1, ]

        # intersect current mask with previous to make new mask        
        # mask is True where masked, so need to be not True at all times
        mask_r2 = ~np.all(~
                 np.concatenate((mask_r2, cum_r2.mask), axis = 0),
                 axis=0
             )[np.newaxis, :, :]
        
        # record the number of pixels we have.  Note, True for masked, so 
        # invert to count non-masked
        n_pix_epoch.append(np.sum(~mask_r2))
       
        if verbose:
            print(f" Done: {n_pix_epoch[-1]} of {total_pix} pixels remain.  ")
            
        
    # apply the mask
    cum_r3_resampled.mask = mask_r2
        
    return n_pix_epoch


#%% calculate_optimal_n_epochs()

def calculate_optimal_n_epochs(
        n_pix_epoch,
        volcano_dir,
        figures,
        ):
    """ 
    Given a list of the number of pixels valid at all times in a time series 
    as this grows in length, calculate the optimal length of time series
    based on the product of the number of pixels and the number of epochs.  
    
    
    Inputs:
        n_pix_epoch | list of ints | number of pixels consistent through 
                                    the time series as we add more and more
                                    epochs.  
                                    
    Returns:
        epoch_metrics \ list of values | measure of quality of the time series
                                        when that number of epochs have been
                                        chosen.  
                                        
        optimal_epoch_n | int | best time series is when this number of epochs
                                has been added.  
                                
    History:
        2025_04_09 | MEG | Written
    """
    import numpy as np
        
    # stores the metric of how good the time series is with that many epochs
    # included.  
    epoch_metrics = []
    for epoch_n, n_pix in enumerate(n_pix_epoch):
        
        # simple metric that takes into account how many times and pixels
        epoch_metrics.append(epoch_n * n_pix)
        
        # this value could also include a factor for ts_length
        # and gap size.  ?
        # instead of the approach I started below that tries to work out which to add
        # here we just add by the pixel order and try to work out when to stop. 
        
    # find the optimal (highest)
    optimal_epoch_n = np.argmax(epoch_metrics)
    
    # plot showing results
    plot_ts_metric_vs_pixels(
        epoch_metrics, 
        n_pix_epoch,
        volcano_dir,
        figures,
        optimal_epoch_n
        )
        
    return epoch_metrics, int(optimal_epoch_n)


#%% plot_ts_metric_vs_pixels()

def plot_ts_metric_vs_pixels(
        epoch_metrics, 
        n_pix_epoch,
        volcano_dir,
        figures, 
        optimal_epoch_n=None
        ):
    """
    Inputs:
        epoch_metrics | list of ints | quality score for the time series 
                                       when a certain number of epochs
                                       have been added.  
        n_pix_epoch | list of ints | number of pixels consistent throughout
                                    time when a certrain number of epochs
                                    have been added.  
        optimal_epoch_n | int | time series is considered optimal when this
                                number of epochs have been added.  
            
    Returns:
        figure
        
    History:
        2025_04_09 | MEG | Written
    
    """
    
    import matplotlib.pyplot as plt
    import numpy as np

    f, ax = plt.subplots(constrained_layout=True)
    
    # **Define colours for each line**
    color1 = 'tab:blue'
    color2 = 'tab:orange'
    
    # Plot the first line on the primary y-axis
    p1, = ax.plot(
        np.arange(len(epoch_metrics)), 
        epoch_metrics,
        color=color1, 
        label='Time series quality metric'
        )
    
    # **Set label and tick colours for the left axis**
    ax.set_ylabel('Time series quality metric', color=color1)
    ax.tick_params(axis='y', labelcolor=color1)
    
    # Create a secondary y-axis
    ax_pixels = ax.twinx()
    p2, = ax_pixels.plot(
        np.arange(len(n_pix_epoch)), 
        n_pix_epoch,
        color=color2, 
        label='Pixel Count'
        )
    
    # **Set label and tick colours for the right axis**
    ax_pixels.set_ylabel('Pixel Count', color=color2)
    ax_pixels.tick_params(axis='y', labelcolor=color2)
    
    # **Combine legends from both axes and add a legend to the primary axis**
    lines = [p1, p2]
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, loc='upper right')
    
    # vertical line on optimal value
    if optimal_epoch_n != None:
        p3 = ax.axvline(
            x=optimal_epoch_n, color='black', linestyle='--', 
            label='Optimal Epoch #'
            )
        
        
    # possibly save and close
    outpath=volcano_dir/'baseline_pixel_selection'
    outpath.mkdir(parents=True, exist_ok=True)
    if figures == 'window':
        pass
    elif figures == "png":
        try:
            f.savefig(
                outpath/"baseline_pixels_vs_epochs_quality_metric.png",
                bbox_inches='tight'
                )
            plt.close(f)
        except:
            print(f"Failed to save the figure.  Trying to continue.  ")
    elif figures == 'both':
        try:
            f.savefig(
                outpath/"baseline_pixels_vs_epochs_quality_metric.png",
                bbox_inches='tight'
                )
        except:
            print(f"Failed to save the figure.  Trying to continue.  ")