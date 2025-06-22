#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 13:41:27 2025

@author: matthew
"""

#%% calculate_valid_pixels()

def calculate_valid_pixels(cum_r3):
    """
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

def calculate_optimal_n_epochs(n_pix_epoch):
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
        
    def plot_ts_metric_vs_pixels(
            epoch_metrics, n_pix_epoch, optimal_epoch_n=None
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

        f, ax = plt.subplots()
        
        # **Define colours for each line**
        color1 = 'tab:blue'
        color2 = 'tab:orange'
        
        # Plot the first line on the primary y-axis
        p1, = ax.plot(
            np.arange(len(epoch_metrics)), epoch_metrics, color=color1, 
            label='Epoch Values'
            )
        
        # **Set label and tick colours for the left axis**
        ax.set_ylabel('Epoch Metric Values', color=color1)
        ax.tick_params(axis='y', labelcolor=color1)
        
        # Create a secondary y-axis
        ax_pixels = ax.twinx()
        p2, = ax_pixels.plot(
            np.arange(len(n_pix_epoch)), n_pix_epoch, color=color2, 
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
                label='Optimal Epoch ~'
                )
        
        

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
    plot_ts_metric_vs_pixels(epoch_metrics, n_pix_epoch, optimal_epoch_n)
        
    return epoch_metrics, optimal_epoch_n
