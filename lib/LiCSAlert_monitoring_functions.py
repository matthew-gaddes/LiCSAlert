#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:09:28 2020

@author: matthew
"""

#%%



def LiCSAlert_monitoring_mode(volcano, LiCSBAS_bin, LiCSAlert_bin, ICASAR_bin, LiCSAR_frames_dir, LiCSAlert_volcs_dir, n_para=1):
    """
       
    Inputs:
        LiCSBAS_bin | string | Path to folder containing LiCSBAS functions.  
        LiCSAlert_bin | string | Path to folder containing LiCSAlert functions.  
        ICASAR_bin | string | Path to folder containing ICASAR functions.  
        LiCSAR_frames_dir | string | path to the folder containing LiCSAR frames.  Needs trailing /
        LiCSAlert_volcs_dir | string | path to the folder containing each volcano.  Needs trailing /
        n_para | int | Sets number of parallel processes used by LiCSBAS.  
    Returns:
        Directory stucture.  
        
    History:
        2020/06/29 | MEG | Written as a script
        2020/07/03 | MEG | Convert to a funtcion
        2020/11/11 | RR | Add n_para argument
        2020/11/16 | MEG | Pass day0_data info to LiCSAlert figure so that x axis is not in terms of days and is instead in terms of dates.  
                
     """
    # 0 Imports etc.:        

    import sys
    import os
    import pickle
    
    if ICASAR_bin not in sys.path:                                                  # check if already on path
        sys.path.append(ICASAR_bin)                                                 # and if not, add
    
    from LiCSAlert_functions import LiCSBAS_for_LiCSAlert, LiCSBAS_to_LiCSAlert, LiCSAlert_preprocessing, LiCSAlert, LiCSAlert_figure
    from LiCSAlert_monitoring_functions import read_config_file, detect_new_ifgs, update_mask_sources_ifgs, record_mask_changes
    from downsample_ifgs import downsample_ifgs
    from ICASAR_functions import ICASAR
        
    # Python version of Tee used to output print functions to the terminal and a log file.  Taken from stack excange
    class Tee(object):
        def __init__(self, *files):
            self.files = files
        def write(self, obj):
            for f in self.files:
                f.write(obj)
                f.flush() # If you want the output to be visible immediately
        def flush(self) :
            for f in self.files:
                f.flush()
    
  
    
    volcano_dir = f"{LiCSAlert_volcs_dir}{volcano}/"
    LiCSBAS_dir = f"{volcano_dir}LiCSBAS/"
    LiCSAR_settings, LiCSBAS_settings, LiCSAlert_settings, ICASAR_settings = read_config_file(f"{volcano_dir}LiCSAlert_settings.txt")                  # read various settings from the volcanoes config file
                                                                                                                                                        # LiCSAR_settings: frame
                                                                                                                                                        # LiCSBAS_settings: lon_lat  
                                                                                                                                                        # ICSAR_settings: n_comp, bootstrapping_param, hdbscan_param, tsne_param, ica_param
    # Determine if new ifgs are present
    new_ifg_flag, LiCSAR_last_acq = detect_new_ifgs(f"{LiCSAR_frames_dir}{LiCSAR_settings['frame']}/GEOC/", volcano_dir)
    

    if new_ifg_flag:                                                                                                                # if new interferograms have been detected, the LiCSAlert outputs can be updated
        # -0: determine if this is the first run:
        previous_LiCSAlert_dates = sorted([f.name for f in os.scandir(volcano_dir) if f.is_dir()])                   # get names of folders produced by LiCSAlert , and keep chronological.  
        unused_folders = ['LiCSBAS', 'ICASAR_results']                                                      # these folders also exist, but aren't dates so we need to delete them.  
        for unused_folder in unused_folders:
            try:
                previous_LiCSAlert_dates.remove(unused_folder)                                                       # loop through trying to delete
            except:
                pass                                                                                    
        if len(previous_LiCSAlert_dates) == 0:                                                                       # if there are no dates, it must be the first time it's been run
            initialise = True
            print(f"This is the first time that LiCSAlert has been run for {volcano}.  ")
        else:
            initialise = False
    
        # 1: Create a folder (YYYYMMDD) for the outputs.  
        try:
            print(f'A new LiCSAR acquisition has been detected.  Creating a folder to contain the new LiCSAlert results ({LiCSAR_last_acq})... ', end = '')
            if not os.path.exists(f"{volcano_dir}{LiCSAR_last_acq}"):
                os.mkdir(f"{volcano_dir}{LiCSAR_last_acq}")                                                                       
            print('Done!')
        except:
            print('Failed!')
            raise Exception('Unable to create folder {LiCSAR_last_acq} - perhaps it already exists?  Exiting...  ')
          
           
        # 2: Create a log file.  
        f = open(f"{volcano_dir}{LiCSAR_last_acq}/LiCSALert_log.txt", 'w')
        original = sys.stdout
        sys.stdout = Tee(sys.stdout, f)
            
        # 3: Create or update the LiCSBAS time series
        try:
            os.mkdir(LiCSBAS_dir)                                                                                                                  # if it's the first run, a folder will be needed for LiCSBAS
        except:
            pass                                                                                                                                   # assume if we can't make it, the folder already exists from a previous run.  
        print(f"Running LiCSBAS.  See 'LiCSBAS_log.txt' for the status of this.  ")
        LiCSBAS_for_LiCSAlert(LiCSAR_settings['frame'], LiCSAR_frames_dir, LiCSBAS_dir, f"{volcano_dir}{LiCSAR_last_acq}/",                        # run LiCSBAS to either create or extend the time series data.  
                              LiCSBAS_bin, LiCSBAS_settings['lon_lat'], n_para=n_para)                                                             # Logfile is sent to the directory for the current date
        displacement_r2, baseline_info = LiCSBAS_to_LiCSAlert(f"{LiCSBAS_dir}TS_GEOCmldir/cum.h5", figures=False)                                  # open the h5 file produced by LiCSBAS
        displacement_r2 = LiCSAlert_preprocessing(displacement_r2, LiCSAlert_settings['downsample_run'], LiCSAlert_settings['downsample_plot'])    # mean centre, and crate downsampled versions (either for general use to make                                                                                                                            # things faster), or just for plotting (to make LiCSAlert figures faster)                         
        
        
        # 4: Possibly run ICASAR
        if 'ICASAR_results' in [f.name for f in os.scandir(volcano_dir) if f.is_dir()]:                                                         # if an ICASAR folder already exists, simply open the saved file rather than rerunning.  
            with open(f"{volcano_dir}ICASAR_results/ICASAR_results.pkl", 'rb') as f:
                sources = pickle.load(f)   
                mask_sources = pickle.load(f)
                tcs  = pickle.load(f)    
                source_residuals = pickle.load(f)    
                Iq_sorted = pickle.load(f)    
                n_clusters = pickle.load(f)    
            f.close()                                                                                                                           
        else:
            print(f"A folder of ICASAR results has not been found so running this now... ", end = '')                                       # or if not, run it
            spatial_ICASAR_data = {'mixtures_r2' : displacement_r2['incremental'],
                                   'mask'        : displacement_r2['mask']}
            
            sources, tcs, residual, Iq, n_clusters, S_all_info, r2_ifg_means  = ICASAR(spatial_data = spatial_ICASAR_data, 
                                                                                       out_folder = f"{volcano_dir}ICASAR_results/", **ICASAR_settings,
                                                                                       ica_verbose = 'short', figures = 'png')
            mask_sources = displacement_r2['mask']                                                                                                          # rename a copy of the mask
            print('Done! ')
        n_baseline_ifgs = tcs.shape[0]                                                                                                                                              # LiCSAlert needs to know how long the baseline stage from ICSASAR was.  
        
        # 5: Deal with changes to the mask of pixels 
        displacement_r2_combined = {}                                                                                                                                               # a new dictionary to save the interferograms sampled to the combined mask in 
        displacement_r2_combined['incremental'], sources_mask_combined, mask_combined = update_mask_sources_ifgs(mask_sources, sources, 
                                                                                                                 displacement_r2['mask'], displacement_r2['incremental'])           # the new mask overwrites the mask in displacement_r2
        displacement_r2_combined['mask'] = mask_combined                                                                                                                            # also put the combined mask in the dictionary
        displacement_r2_combined["incremental_downsampled"], displacement_r2_combined["mask_downsampled"] = downsample_ifgs(displacement_r2_combined["incremental"], displacement_r2_combined["mask"],
                                                                                                                            LiCSAlert_settings['downsample_plot'], verbose = False)
        if initialise:
            previous_output_dir = None                                                                                                                                  # there is no previous output directory
        else:
            previous_output_dir = f"{volcano_dir}{previous_LiCSAlert_dates[-1]}"
        record_mask_changes(mask_sources, displacement_r2['mask'], mask_combined, LiCSAR_last_acq, f"{volcano_dir}{LiCSAR_last_acq}/", previous_output_dir)             # record any changes in the mask (ie pixels that are now masked due to being incoherent).  
        
        # 6: Run LiCSAlert
        #def LiCSAlert(sources, time_values, ifgs_baseline, ifgs_monitoring = None, t_recalculate = 10, verbose=False):
        
        sources_tcs_baseline, residual_tcs_baseline = LiCSAlert(sources_mask_combined, baseline_info["baselines_cumulative"],                                                               # the LiCSAlert algoirthm, using the sources with the combined mask (sources_mask_combined)
                                                                displacement_r2_combined['incremental'][:n_baseline_ifgs,], displacement_r2_combined['incremental'][n_baseline_ifgs:,],     # baseline ifgs and monitoring ifgs
                                                                t_recalculate=10, verbose=False)                                                                                            # recalculate lines of best fit every 10 acquisitions
        
        LiCSAlert_figure(sources_tcs_baseline, residual_tcs_baseline, sources_mask_combined, displacement_r2_combined, n_baseline_ifgs,                                 # creat the LiCSAlert figure
                         baseline_info["baselines_cumulative"], out_folder = f"{volcano_dir}{LiCSAR_last_acq}", day0_date = baseline_info['imdates'][0])                #
        
        sys.stdout = original                                                                                                                       # return stdout to be normal.  
        f.close()                                                                                                                                   # and close the log file.  
    else:
        print(f'No new LiCSAR acquisitions have been detected for {volcano}, and the LiCSAlert outputs are thought to be up to date.  ')
        
    
    
    




#%%


def record_mask_changes(mask_sources, mask_ifgs, mask_combined, current_date, current_output_dir, previous_output_dir = None):
    """ Record changes to the masks used in LiCSAlert, as this is dependent on the mask provided by LiCSBAS.  Creates a variety of .png images showing the mask,
    how many pixels remain for LiCSAlert to use, and a .pkl so this can be compared to the last time it was run.  
    
    Inputs:
        mask_sources | r2 array | the mask used by ICASAR
        mask_ifgs | r2 array | the mask produced by the last run of LiCSBAS
        mask_combined | r2 array | the mask that removes any pixels that aren't in bothh the sources and the ifgs
        current_date | string | the date that LiCSAlert is being run to.  
        current_output_dir | string | the folder that LiCSALert is currently outputting to
        previous_output_dir | string | If it's not hte first time LiCSAlert was run, this is the folder that LiCSALert previously output to
    Returns:
        2 x png figures
        .pkl of the masks and dates.  
    History:
        2020/06/25 | MEG | Written
        2020/07/01 | MEG | Major rewrite to suit directory based structure.  
        2020/07/03 | MEG | continue major rewrite, and write docs.  
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import numpy.ma as ma
    import pickle
    
    # 0: try to open a .pkl with the mask changes in, or initiate it
    try:
        with open(f"{previous_output_dir}/mask_history.pkl", 'rb') as f:
            dates = pickle.load(f)   
            masks_combined = pickle.load(f)
            masks_ifgs = pickle.load(f)
        f.close()
        initialising = False
    except:                                                                             # if we can't open file, assume it is because it doesn't exist as first run of function.  
        initialising = True                                                             # this flag used to control plotting as it's different for the first one.  
        dates = []
        masks_combined = []
        masks_ifgs = []
    
    # 1: append current masks
    dates.append(current_date)
    masks_combined.append(mask_combined)
    masks_ifgs.append(mask_ifgs)
    
    # 2: Save the file
    with open(f'{current_output_dir}mask_history.pkl', 'wb') as f:
        pickle.dump(dates, f)
        pickle.dump(masks_combined, f)
        pickle.dump(masks_ifgs, f)
    f.close()
    
    
    # 3: Figure showing the masks
    f1,axes = plt.subplots(1,4, figsize = (12,6))
    axes[0].imshow(mask_sources)
    axes[0].set_title('(ICASAR) sources mask')
    axes[1].imshow(masks_ifgs[-1])
    axes[1].set_title('Current LiCSBAS mask')
    axes[2].imshow(masks_combined[-1])
    axes[2].set_title('Current combined mask')
    if initialising:                                                                                        # first run, so won't be a difference between this and past mask
        axes[3].imshow(np.full(mask_combined.shape, np.nan))                                                # just nans so it's blank
        axes[3].set_title("First mask so no change")
    else:
        axes[3].imshow(ma.masked_where(masks_combined[-2] == masks_combined[-1] , np.ones(mask_combined.shape)))           # create an array of ones that is maksed everywhere that the two masks are the same, and then plot this.  
        axes[3].set_title(f"{dates[-2]} - {dates[-1]} ")                                                # and differnce inthe dates
        
    f1.canvas.set_window_title(current_date)
    f1.savefig(f"{current_output_dir}mask_changes.png", bbox_inches='tight')
    plt.close(f1)

    # 4 Figure showing how the number of pixels has changes with time
    n_updates = len(dates)
    x_vals = np.arange(n_updates)
    n_pixs = np.zeros((n_updates, 2))                                     # 1st column will be number of non-masked and 2nd number of masked
    for ifg_n in range(n_updates):
        n_pixs[ifg_n,0] = len(np.argwhere(masks_combined[ifg_n] == False))
        n_pixs[ifg_n,1] = len(np.argwhere(masks_combined[ifg_n] == True))
    
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
    ax.set_xticklabels(dates)
    #ax.set_ylim(bottom = 0)
    #plt.grid()
    f2.savefig(f"{current_output_dir}mask_changes_graph.png", bbox_inches='tight')
    plt.close(f2)
        




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
    
    
    
    def check_LiCSAlert_products(dates):
        """ Given a list of dates in which LiCSAlert has been run, check that the required outputs are present in each folder.  
        Inputs:
            dates | list of strings | dates that LiCSAlert was run until.  In form YYYYMMDD
        Returns:
            dates_incomplete | list of strings | dates that a LiCSAlert folder exisits, but it doesn't have all the ouptuts.  
        History:
            2020/11/13 | MEG | Written
        """
        from pathlib import Path
        import os
        import fnmatch                                                                      # used to compare lists and strings using wildcards
        
        constant_outputs = ['mask_changes_graph.png',                                      # The output files that are expected to exist and never change name
                            'mask_changes.png',
                            'mask_history.pkl']
        variable_outputs = ['LiCSAlert_figure_with_*_monitoring_interferograms.png']        # The output files that are expected to exist and change name.  

        dates_incomplete = []
        for date in dates:
            date_folder = Path(folder_LiCSAlert + date)                                     # join to make a path to the current date folder
            date_files = sorted([f.name for f in os.scandir(date_folder)])                  # get names of the files in this folder (no paths)
            
            # 1: look for the products that don't change name (ie all but the first)
            all_products_complete = True                                                                            # initiate as True
            for constant_output in constant_outputs:                                                                # loop through the outputs that can't change name
                all_products_complete = (all_products_complete) and (constant_output in date_files)                 # update boolean 
            # 2: look for the product that does change name (the main figure)
            for variable_output in variable_outputs:                                                                # loop through the outputs that can change name
                output = fnmatch.filter(date_files, "LiCSAlert_figure_with_*_monitoring_interferograms.png")        # check for file with wildcard for changing name
                all_products_complete = (all_products_complete) and (len(output) > 0)                               # empty list if file doesn't exit, use to update boolean
                
            if all_products_complete is False:
                dates_incomplete.append(date)                                                                       # create a list of dates for which otputs are missing
        return dates_incomplete
        
    
   
    # 0: Get the last acquisition that LiCSAR has been run until.  
    LiCSAR_ifgs = sorted([f.name for f in os.scandir(folder_ifgs) if f.is_dir()])                # get names of folders produced by LiCSAR (ie the ifgs), and keep chronological.  
    # check for empty directory:
    if not LiCSAR_ifgs:                                                                          # RR addition?  To check that the list isn't empty
        print(f"No files found in {folder_ifgs} ... ")
        return False, False                                                                     # return back to parent function (new_ifg_flag, LiCSAR_last_acq)
    LiCSAR_last_acq = LiCSAR_ifgs[-1][-8:]                                                      # this is the last date that LiCSAR has processed up to, in the form YYYYMMDD
    
    # 1: Get the last date that LiCAlert has been run until
    LiCSAlert_dates = sorted([f.name for f in os.scandir(folder_LiCSAlert) if f.is_dir()])      # get names of folders produced by LiCSAR (ie the ifgs), and keep chronological.  
    for unneeded_folder in ['LiCSBAS', 'ICASAR_results']:                                       # these folders get caught in the dates list, but aren't dates so need to be deleted.  
        try:
            LiCSAlert_dates.remove(unneeded_folder)                                             # note that the LiCSBAS folder also gets caught by this, and needs removing as it's not a date.  
        except:
            pass                                                                                # however, on the first ever run these don't exist.  
    if len(LiCSAlert_dates) == 0:                                                                # if the list is of length 0, LiCSAlert has not been run yet.  
        new_ifgs_flag = True                                                                    # if LiCSAlert hasn't been run yet but the LiCSAR_ifgs is not empty, there must be new interferograms.  
        return new_ifgs_flag, LiCSAR_last_acq                                                   # return to parent function.  
    else:
        LiCSAlert_last_run = LiCSAlert_dates[-1]                                            #    folder are YYYYMMDD so last one is last time it was run until.  

    # 2: Compare 0 (the last LiCSAR acquisition) and 1 (the last date LiCSAlert has been run until)
    if len(LiCSAlert_dates) == 0:                                                                       # if there are no LiCSAlert dates, it hasn't been run for this volcano.  
        #print("It appears that LiCSAlert hasn't been run for this volcano yet.  Setting the 'new_ifgs_flag' to True.  ")
        new_ifgs_flag = True
    else:
        fmt = '%Y%m%d'                                                                      # tell datetime the format of the date (here year, month, day with no sepeartions)
        LiCSAR_last_acq_dt = datetime.datetime.strptime(LiCSAR_last_acq, fmt)               # convert from string to datetime
        LiCSAlert_last_run_dt = datetime.datetime.strptime(LiCSAlert_last_run, fmt)         # ditto
        if LiCSAR_last_acq_dt > LiCSAlert_last_run_dt:                                      # if the last licsar ifg is after the last LiCSAlert run.  
            new_ifgs_flag = True
        else:
            new_ifgs_flag = False

    # 3: Check that LiCSAlert has run succesfully before in each folder
    dates_incomplete = check_LiCSAlert_products(LiCSAlert_dates)
    if len(dates_incomplete) > 0:
        print(f"LiCSBAS products have not be found for {dates_incomplete}")
        
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
