#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:09:28 2020

@author: matthew
"""

#%%


def record_mask_changes(mask_sources, mask_ifgs, mask_shared):
    """ 
    
    
    
    Given a new mask and the date that the time series until, update the mask_change dictionary that keeps track of this.  
    Inputs:
        current_mask | rank 2 boolean | the current mask, as produced by LiCSBAS
        current_date | string | the most recent data that the time series spans until, in form yyyymmdd
        mask_change | dict or None | dictionary where each key is the date the time series ends up, and the variable the mask
        png_path | string | path to folder of where to save the .png files.  Needs trailing /
    Returns:
        mask_change | dict or None | dictionary where each key is the date the time series ends up, and the variable the mask
    History:
        2020/06/25 | MEG | Written
        2020/07/01 | MEG | Major rewrite to suit directory based structure.  
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import numpy.ma as ma
        
    if mask_change == None:                                                                 # when the function is first called, this won't exist (so is set to None)
        initialising = True
        mask_change = {}                                                                    # initiatalise 
        mask_change['dates'] = [current_date]                                               # one list will have the dates for the mask, initiate with first date
        mask_change['masks'] = [current_mask]                                                       # and one will have the masks, initate with first mask
    else:
        initialising = False                                                                # if mask is not None, we aren't initiliasing as it must have previously run in order to create mask_change
        mask_change['dates'].append(current_date)
        mask_change['masks'].append(current_mask)
    
    #import ipdb; ipdb.set_trace()
    # Save the mask as a png, and the change in mask to the previous one
    f1,axes = plt.subplots(1,2)
    axes[0].imshow(mask_change['masks'][-1])
    axes[0].set_title(mask_change['dates'][-1])
    if initialising:                                                                                        # first run, so won't be a difference between this and past mask
        axes[1].imshow(ma.array(np.ones(current_mask.shape), mask = np.ones(current_mask.shape)))           # an array of ones, but every element is masked as we know it hasn't changed yet.  
        axes[1].set_title("First mask so no change")
    else:
        axes[1].imshow(ma.masked_where(mask_change['masks'][-2] == mask_change['masks'][-1] , np.ones(current_mask.shape)))           # create an array of ones that is maksed everywhere that the two masks are the same, and then plot this.  
        #axes[1].imshow((mask_change['masks'][-2]).astype(int) - (mask_change['masks'][-1]).astype(int))         # difference in the masks
        axes[1].set_title(f"{mask_change['dates'][-2]} - {mask_change['dates'][-1]} ")      # and differnce inthe dates
        
    f1.canvas.set_window_title(current_date)
    f1.savefig(f"{png_path}{mask_change['dates'][-1]}.png", bbox_inches='tight')
    plt.close(f1)

    # Create the line graph of number of pixels
    n_updates = len(mask_change['dates'])
    x_vals = np.arange(n_updates)
    n_pixs = np.zeros((n_updates, 2))                                     # 1st column will be number of non-masked and 2nd number of masked
    for ifg_n in range(n_updates):
        n_pixs[ifg_n,0] = len(np.argwhere(mask_change['masks'][ifg_n] == False))
        n_pixs[ifg_n,1] = len(np.argwhere(mask_change['masks'][ifg_n] == True))
    
    #import ipdb; ipdb.set_trace()
    f2,ax = plt.subplots(1)
    ax.plot(x_vals, n_pixs[:,0], label = 'Non-masked pixels')
    ax.plot(x_vals, n_pixs[:,1], label = 'Masked pixels')
    ax.scatter(x_vals, n_pixs[:,0])
    ax.scatter(x_vals, n_pixs[:,1])
    ax.set_title(current_date)
    f2.canvas.set_window_title(current_date)
    leg = ax.legend()
    ax.set_ylabel('# pixels')
    ax.set_xlabel('Time series end (yyyymmdd)')
    ax.set_xticks(x_vals)
    ax.set_xticklabels(mask_change['dates'])
    ax.set_ylim(bottom = 0)
    #plt.grid()
    f2.savefig(f'{png_path}pixel_history_{current_date}.png', bbox_inches='tight')
    plt.close(f2)
        
    return mask_change


# def record_mask_changes(current_mask, current_date, mask_change = None,
#                         png_path = './mask_history/'):
#     """ Given a new mask and the date that the time series until, update the mask_change dictionary that keeps track of this.  
#     Inputs:
#         current_mask | rank 2 boolean | the current mask, as produced by LiCSBAS
#         current_date | string | the most recent data that the time series spans until, in form yyyymmdd
#         mask_change | dict or None | dictionary where each key is the date the time series ends up, and the variable the mask
#         png_path | string | path to folder of where to save the .png files.  Needs trailing /
#     Returns:
#         mask_change | dict or None | dictionary where each key is the date the time series ends up, and the variable the mask
#     History:
#         2020/06/25 | MEG | Written
#     """
#     import matplotlib.pyplot as plt
#     import numpy as np
#     import numpy.ma as ma
        
#     if mask_change == None:                                                                 # when the function is first called, this won't exist (so is set to None)
#         initialising = True
#         mask_change = {}                                                                    # initiatalise 
#         mask_change['dates'] = [current_date]                                               # one list will have the dates for the mask, initiate with first date
#         mask_change['masks'] = [current_mask]                                                       # and one will have the masks, initate with first mask
#     else:
#         initialising = False                                                                # if mask is not None, we aren't initiliasing as it must have previously run in order to create mask_change
#         mask_change['dates'].append(current_date)
#         mask_change['masks'].append(current_mask)
    
#     #import ipdb; ipdb.set_trace()
#     # Save the mask as a png, and the change in mask to the previous one
#     f1,axes = plt.subplots(1,2)
#     axes[0].imshow(mask_change['masks'][-1])
#     axes[0].set_title(mask_change['dates'][-1])
#     if initialising:                                                                                        # first run, so won't be a difference between this and past mask
#         axes[1].imshow(ma.array(np.ones(current_mask.shape), mask = np.ones(current_mask.shape)))           # an array of ones, but every element is masked as we know it hasn't changed yet.  
#         axes[1].set_title("First mask so no change")
#     else:
#         axes[1].imshow(ma.masked_where(mask_change['masks'][-2] == mask_change['masks'][-1] , np.ones(current_mask.shape)))           # create an array of ones that is maksed everywhere that the two masks are the same, and then plot this.  
#         #axes[1].imshow((mask_change['masks'][-2]).astype(int) - (mask_change['masks'][-1]).astype(int))         # difference in the masks
#         axes[1].set_title(f"{mask_change['dates'][-2]} - {mask_change['dates'][-1]} ")      # and differnce inthe dates
        
#     f1.canvas.set_window_title(current_date)
#     f1.savefig(f"{png_path}{mask_change['dates'][-1]}.png", bbox_inches='tight')
#     plt.close(f1)

#     # Create the line graph of number of pixels
#     n_updates = len(mask_change['dates'])
#     x_vals = np.arange(n_updates)
#     n_pixs = np.zeros((n_updates, 2))                                     # 1st column will be number of non-masked and 2nd number of masked
#     for ifg_n in range(n_updates):
#         n_pixs[ifg_n,0] = len(np.argwhere(mask_change['masks'][ifg_n] == False))
#         n_pixs[ifg_n,1] = len(np.argwhere(mask_change['masks'][ifg_n] == True))
    
#     #import ipdb; ipdb.set_trace()
#     f2,ax = plt.subplots(1)
#     ax.plot(x_vals, n_pixs[:,0], label = 'Non-masked pixels')
#     ax.plot(x_vals, n_pixs[:,1], label = 'Masked pixels')
#     ax.scatter(x_vals, n_pixs[:,0])
#     ax.scatter(x_vals, n_pixs[:,1])
#     ax.set_title(current_date)
#     f2.canvas.set_window_title(current_date)
#     leg = ax.legend()
#     ax.set_ylabel('# pixels')
#     ax.set_xlabel('Time series end (yyyymmdd)')
#     ax.set_xticks(x_vals)
#     ax.set_xticklabels(mask_change['dates'])
#     ax.set_ylim(bottom = 0)
#     #plt.grid()
#     f2.savefig(f'{png_path}pixel_history_{current_date}.png', bbox_inches='tight')
#     plt.close(f2)
        
#     return mask_change



#%%

def update_mask_sources_ifgs(mask_sources, sources, mask_ifgs, ifgs):
    """ Given two masks of pixels, create a mask of pixels that are valid for both.  
    Inputs:
        mask_sources | boolean rank 2| original mask
        sources  | r2 array | sources as row vectors
        mask_ifgs | boolean rank 2| new mask
    Returns (in initiate mode):
        mask_combnied | boolean rank 2| original mask
    History:
        2020/02/19 | MEG |  Written      
        2020/06/26 | MEG | Major rewrite.  
    """
    import numpy as np
    import numpy.ma as ma
    from LiCSAlert_aux_functions import col_to_ma
    
    
    def apply_new_mask(ifgs, mask_old, mask_new):
        """Apply a new mask to a collection of ifgs (or sources) that are stored as row vectors with an accompanying mask.  
        Inputs:
            ifgs | r2 array | ifgs as row vectors
            mask_old | r2 array | mask to convert a row of ifg into a rank 2 masked array
            mask_new | r2 array | the new mask to be applied.  Note that it must not unmask any pixels that are already masked.  
        Returns:
            ifgs_new_mask | r2 array | as per ifgs, but with a new mask.  
        History:
            2020/06/26 | MEG | Written
        """
        n_pixs_new = len(np.argwhere(mask_new == False))                                        
        ifgs_new_mask = np.zeros((ifgs.shape[0], n_pixs_new))                        # initiate an array to store the modified sources as row vectors    
        for ifg_n, ifg in enumerate(ifgs):                                 # Loop through each source
            ifg_r2 = col_to_ma(ifg, mask_old)                             # turn it from a row vector into a rank 2 masked array        
            ifg_r2_new_mask = ma.array(ifg_r2, mask = mask_new)              # apply the new mask   
            ifgs_new_mask[ifg_n, :] = ma.compressed(ifg_r2_new_mask)       # convert to row vector and places in rank 2 array of modified sources
        return ifgs_new_mask
    
    
    mask_both = ~np.logical_and(~mask_sources, ~mask_ifgs)                                       # make a new mask for pixels that are in the sources AND in the current time series
    n_pixs_sources = len(np.argwhere(mask_sources == False))                                  # masked pixels are 1s, so invert with 1- bit so that non-masked are 1s, then sum to get number of pixels
    n_pixs_new = len(np.argwhere(mask_ifgs == False))                                          # ditto for new mask
    n_pixs_both = len(np.argwhere(mask_both == False))                                        # ditto for the mutual mask
    print(f"Updating masks and ICA sources.  Of the {n_pixs_sources} in the sources and {n_pixs_new} in the current LiCSBAS time series, "
          f"{n_pixs_both} are in both and can be used in this iteration of LiCSAlert.  ")
    
    ifgs_new_mask = apply_new_mask(ifgs, mask_ifgs, mask_both)
    sources_new_mask = apply_new_mask(sources, mask_sources, mask_both)
    
    return ifgs_new_mask, sources_new_mask, mask_both
    


#%%
 
def detect_new_ifgs(folder_ifgs, folder_LiCSAlert):
    """ Determine the number of LiCSAR ifgs in a folder, and detect if this changes.  Note that it has different returns, depending on if 
    it is in the simple "initate" mode or not (if not, also returns a flag of if new interferograms have been detected).  
    Inputs:

    Rerturns:


    History:
        2020/02/18 | MEG |  Written        
        2020/06/29 | MEG | Major rewrite to use a folder based structure    

    """
    import os 
    import datetime
    
    # 0: Get the last acquisition that LiCSAR has been run until.  
    LiCSAR_ifgs = sorted([f.name for f in os.scandir(folder_ifgs) if f.is_dir()])     # get names of folders produced by LiCSAR (ie the ifgs), and keep chronological.  
    LiCSAR_last_acq = LiCSAR_ifgs[-1][-8:]                                                  # this is the last date that LiCSAR has processed up to 
    
    # 1: Get the last date that LiCAlert has been run until
    LiCSAlert_dates = sorted([f.name for f in os.scandir(folder_LiCSAlert) if f.is_dir()])     # get names of folders produced by LiCSAR (ie the ifgs), and keep chronological.  
    try:
        LiCSAlert_dates.remove('LiCSBAS')                                                       # note that the LiCSBAS folder also gets caught by this, and needs removing as it's not a date.  
    except:
        pass                                                                                    # however, on the first ever run it doesn't exist.  
    try:
        LiCSAlert_dates.remove('ICASAR_results')                                                       # note that the ICSAR folder also gets caught by this, and needs removing as it's not a date.  
    except:
        pass                                                                                    # however, on the first ever run it doesn't exist.  
        
    
    if len(LiCSAlert_dates) == 0:
        print("It appears that LiCSAlert hasn't been run for this volcano yet.  Setting the 'new_ifgs_flag' to True.  ")
        new_ifgs_flag = True
    else:
        LiCSAlert_last_run = LiCSAlert_dates[-1]                                            # folder are YYYYMMDD so last one is last time it was run until.  
        fmt = '%Y%m%d'                                                                      # tell datetime the format of the date (here year, month, day with no sepeartions)
        LiCSAR_last_acq_dt = datetime.datetime.strptime(LiCSAR_last_acq, fmt)               # convert from string to datetime
        LiCSAlert_last_run_dt = datetime.datetime.strptime(LiCSAlert_last_run, fmt)         # ditto
        if LiCSAR_last_acq_dt > LiCSAlert_last_run_dt:
            new_ifgs_flag = True
        else:
            new_ifgs_flag = False
            
    return new_ifgs_flag, LiCSAR_last_acq
    


#%%

def read_config_file(config_file):
    """Given a .txt file of arguments, read it into dictionaries
    Inputs:
        config_file | string | .txt file to open
    Returns:
        LiCSAR_settings | dict | 
        LiCSBAS_settings | dict | as per used by map_profile_wrapper
        ICASAR_settings | dict | as per used by map_profile_wrapper
       
    History:
        2020/05/29 | MEG | Written
        2020/06/29 | MEG | Modified for use with LiCSAlert
        2020/06/30 | MEG | Add LiCSAlert settings
    """
    import configparser    
   
   
    LiCSAR_settings = {}                                                                       # initiate
    LiCSBAS_settings = {}
    LiCSAlert_settings = {}
    ICASAR_settings = {}
   
    config = configparser.ConfigParser()                                                       # read the config file
    config.read_file(open(config_file))
   
    LiCSAR_settings['frame'] = config.get('LiCSAR', 'frame')                                    # 1:  LiCSAR settings
    
    west = float(config.get('LiCSBAS', 'west'))                                                 # 2: LiCSBAS
    east = float(config.get('LiCSBAS', 'east'))                                                 # Some need to be extracted individually...
    south = float(config.get('LiCSBAS', 'south'))
    north = float(config.get('LiCSBAS', 'north'))
    LiCSBAS_settings['lon_lat'] = [west, east, south, north]                                    # and then merged together into one item in the dictionary (e.g. a list or tuple)
    
    LiCSAlert_settings['downsample_run'] = float(config.get('LiCSAlert', 'downsample_run'))       # 3 LiCSAlert settings
    LiCSAlert_settings['downsample_plot'] = float(config.get('LiCSAlert', 'downsample_plot'))                 
    
    ICASAR_settings['n_comp'] = int(config.get('ICASAR', 'n_comp'))                             # 4: ICASAR settings
    n_bootstrapped =  int(config.get('ICASAR', 'n_bootstrapped'))                 
    n_not_bootstrapped =  int(config.get('ICASAR', 'n_not_bootstrapped'))                 
    ICASAR_settings['bootstrapping_param'] = (n_bootstrapped, n_not_bootstrapped)
    
    HDBSCAN_min_cluster_size =  int(config.get('ICASAR', 'HDBSCAN_min_cluster_size'))                 
    HDBSCAN_min_samples =  int(config.get('ICASAR', 'HDBSCAN_min_samples'))                 
    ICASAR_settings['hdbscan_param'] = (HDBSCAN_min_cluster_size, HDBSCAN_min_samples)
    
    tsne_perplexity = int(config.get('ICASAR', 'tsne_perplexity'))                 
    tsne_early_exaggeration = int(config.get('ICASAR', 'tsne_early_exaggeration'))                 
    ICASAR_settings['tsne_param'] = (tsne_perplexity, tsne_early_exaggeration)
    
    ica_tolerance = float(config.get('ICASAR', 'ica_tolerance'))                 
    ica_max_iterations = int(config.get('ICASAR', 'ica_max_iterations'))                 
    ICASAR_settings['ica_param'] = (ica_tolerance, ica_max_iterations)
    
    return LiCSAR_settings, LiCSBAS_settings, LiCSAlert_settings, ICASAR_settings