#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:09:28 2020

@author: matthew
"""



#%%



def LiCSAlert_monitoring_mode(volcano, LiCSBAS_bin, LiCSAlert_bin, ICASAR_bin, LiCSAR_frames_dir, LiCSAlert_volcs_dir, n_para=1):
    """The main function for running LiCSAlert is monitoring mode.  It was designed to work with LiCSAR interferograms, but could be used with 
    any product that creates interfegorams that LiCSBAS can use (which is used for the time series calculation.  )
    
    One of the key functions called in this function is "run_LiCSAlert_status", which returns LiCSAlert status.  
       
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
    import shutil
    import copy
    
    if ICASAR_bin not in sys.path:                                                  # check if already on path
        sys.path.append(ICASAR_bin)                                                 # and if not, add
    
    from LiCSAlert_functions import LiCSBAS_for_LiCSAlert, LiCSBAS_to_LiCSAlert, LiCSAlert_preprocessing, LiCSAlert, LiCSAlert_figure, shorten_LiCSAlert_data
    from LiCSAlert_monitoring_functions import read_config_file, detect_new_ifgs, update_mask_sources_ifgs, record_mask_changes
    from LiCSAlert_aux_functions import compare_two_dates, Tee, get_baseline_end_ifg_n
    from downsample_ifgs import downsample_ifgs
    from ICASAR_functions import ICASAR
        
    # 0: begin
    volcano_dir = f"{LiCSAlert_volcs_dir}{volcano}/"
    LiCSBAS_dir = f"{volcano_dir}LiCSBAS/"
    LiCSAR_settings, LiCSBAS_settings, LiCSAlert_settings, ICASAR_settings = read_config_file(f"{volcano_dir}LiCSAlert_settings.txt")                      # read various settings from the volcanoes config file
                                                                                                                                                           # LiCSAR_settings: frame | LiCSBAS_settings: lon_lat | ICSAR_settings: n_comp, bootstrapping_param, hdbscan_param, tsne_param, ica_param
    # 1: Determine the status of LiCSAlert, and update the user.      
    LiCSAlert_status = run_LiCSAlert_status(f"{LiCSAR_frames_dir}{LiCSAR_settings['frame']}/GEOC/", volcano_dir, LiCSAlert_settings['baseline_end'],       # Determine the status for LiCSAlert for this volcano
                                            f"{volcano_dir}LiCSAlert_history.txt")                                                                         # note that this logs by appending to a file in the volcano's directory.  
    print(f"This is the current 'LiCSAlert_status' for this folder: {LiCSAlert_status}")                                                                   #  print the contents of the dict which are fairly straightforward.  
    
    if LiCSAlert_status['run_LiCSAlert']:
        if not os.path.exists(f"{volcano_dir}{LiCSAlert_status['LiCSAR_last_acq']}"):
            os.mkdir(f"{volcano_dir}{LiCSAlert_status['LiCSAR_last_acq']}")                                                                       
        f_run_log = open(f"{volcano_dir}{LiCSAlert_status['LiCSAR_last_acq']}/LiCSAlert_log.txt", 'w')                                                                          # 2020/11/16 | MEG | Fix typo in LiCSAlert.  
        original = sys.stdout
        sys.stdout = Tee(sys.stdout, f_run_log)
    
        
        # 2: if required, run LiCSBAS
        if LiCSAlert_status['run_LiCSBAS']:
            try:
                os.mkdir(LiCSBAS_dir)                                                                                                                  # if it's the first run, a folder will be needed for LiCSBAS
            except:
                pass                                                                                                                                   # assume if we can't make it, the folder already exists from a previous run.  
            print(f"Running LiCSBAS.  See 'LiCSBAS_log.txt' for the status of this.  ")
            LiCSBAS_for_LiCSAlert(LiCSAR_settings['frame'], LiCSAR_frames_dir, LiCSBAS_dir, f"{volcano_dir}{LiCSAlert_status['LiCSAR_last_acq']}/",    # run LiCSBAS to either create or extend the time series data.  
                                  LiCSBAS_bin, LiCSBAS_settings['lon_lat'], n_para=n_para)                                                             # Logfile is sent to the directory for the current date
            displacement_r2, temporal_baselines = LiCSBAS_to_LiCSAlert(f"{LiCSBAS_dir}TS_GEOCmldir/cum.h5", figures=False)                             # open the h5 file produced by LiCSBAS, lons and lats are in geocode info and same resolution as the ifgs
            displacement_r2 = LiCSAlert_preprocessing(displacement_r2, LiCSAlert_settings['downsample_run'], LiCSAlert_settings['downsample_plot'])    # mean centre, and create downsampled versions (both for general use and again for plotting)
            
        # 3: If required, run ICASA to find the sources, or just load the sources from a previous run
        if LiCSAlert_status['run_ICASAR']:
            print(f"Running ICASAR... ", end = '')                                       # or if not, run it
            LiCSAlert_settings['baseline_end_ifg_n'] = get_baseline_end_ifg_n(temporal_baselines['acq_dates'], LiCSAlert_settings['baseline_end'])               # if this is e.g. 14, the 14th ifg would not be in the baseline stage
            spatial_ICASAR_data = {'mixtures_r2' : displacement_r2['incremental'][:(LiCSAlert_settings['baseline_end_ifg_n']+1),],                             # only take up to the last 
                                   'mask'        : displacement_r2['mask']}
            
            sources, tcs, residual, Iq, n_clusters, S_all_info, r2_ifg_means  = ICASAR(spatial_data = spatial_ICASAR_data,                                      # Run ICASAR (slow))
                                                                                       out_folder = f"{volcano_dir}ICASAR_results/", **ICASAR_settings,
                                                                                       ica_verbose = 'short', figures = 'png',
                                                                                       lons = displacement_r2['lons'], lats = displacement_r2['lats'])          # ICs are geocoded as kmz, so need the lon and lats here.  
            mask_sources = displacement_r2['mask']                                                                                                              # rename a copy of the mask
            print('Done! ')
        else:
            with open(f"{volcano_dir}ICASAR_results/ICASAR_results.pkl", 'rb') as f_icasar:                                                                     # Or open the products from a previous ICASAR run.  
                sources = pickle.load(f_icasar)   
                mask_sources = pickle.load(f_icasar)
                tcs  = pickle.load(f_icasar)    
                source_residuals = pickle.load(f_icasar)    
                Iq_sorted = pickle.load(f_icasar)    
                n_clusters = pickle.load(f_icasar)    
            f_icasar.close()                                                                                                                           
            LiCSAlert_settings['baseline_end_ifg_n'] = get_baseline_end_ifg_n(temporal_baselines['acq_dates'], LiCSAlert_settings['baseline_end'])            # if this is e.g. 14, the 14th ifg would not be in the baseline stage
    
        # 5: Deal with changes to the mask of pixels 
        #import pdb; pdb.set_trace()
        displacement_r2_combined = {}                                                                                                                                               # a new dictionary to save the interferograms sampled to the combined mask in (ie the same pixels as the sources)
        displacement_r2_combined['incremental'], sources_mask_combined, mask_combined = update_mask_sources_ifgs(mask_sources, sources,                                             # comapres the mask for the ICs and the mask for the current ifgs and finds a set of pixels that are present in both
                                                                                                                 displacement_r2['mask'], displacement_r2['incremental'])           # the new mask overwrites the mask in displacement_r2
        displacement_r2_combined['mask'] = mask_combined                                                                                                                            # also put the combined mask in the dictionary
        displacement_r2_combined["incremental_downsampled"], displacement_r2_combined["mask_downsampled"] = downsample_ifgs(displacement_r2_combined["incremental"], displacement_r2_combined["mask"],      # create the downsampled plotting ones
                                                                                                                            LiCSAlert_settings['downsample_plot'], verbose = False)
        
        # note - what will happen to existing products in the processed_with_errors folders?
        
        # 6: Main loop to run LiCSAlert for each date that is required
        processing_dates = copy.deepcopy(LiCSAlert_status['pending'])
        for processed_with_error in LiCSAlert_status['processed_with_errors']:
            if processed_with_error not in processing_dates:
                processing_dates.append(processed_with_error)
        processing_dates = sorted(processing_dates)
        print(f"LiCSAlert will be run for the following dates: {processing_dates}")
        for processing_date in processing_dates:
            print(f"Running LiCSAlert for {processing_date}")
            ifg_n = temporal_baselines['acq_dates'].index(processing_date)
            
            # 6a: Create a folder (YYYYMMDD) for the outputs.  
            if not os.path.exists(f"{volcano_dir}{processing_date}"):                                   # True if folder exists, so enter if statement if doesn't exist (due to not)
                os.mkdir(f"{volcano_dir}{processing_date}")                                            # if doesn't exist, make it                           
            else:                                                                                       # if the folder does already exist
                if processing_date == LiCSAlert_status['LiCSAR_last_acq']:                              # if the folder exists and was used for the log file:
                    pass              
                else:
                    print(f"The folder {processing_date} appears to exists already.  This is usually due to the date not having all the required LiCSAlert products, and LiCSAlert"
                          f" is now trying to fill this date again.  ")
                    shutil.rmtree(f"{volcano_dir}{processing_date}")                                        # delete the folder and all its contents
                    os.mkdir(f"{volcano_dir}{processing_date}")                                             # and remake the folder
                
            # 6b: Update the mask to create a figure of changes.  Not robustly written
            previous_date = temporal_baselines['acq_dates'][ifg_n-1]                                                                  # find the date one before the one being processed.  
            previous_date_after_baseline = compare_two_dates(LiCSAlert_settings['baseline_end'], previous_date)                     # check that this date is not during the baseline.  
            if previous_date_after_baseline:                                                                                        # if it's not,
                previous_output_dir = f"{volcano_dir}{temporal_baselines['acq_dates'][ifg_n-1]}"                                      # get the previous output directory
            else:
                previous_output_dir = None                                                                                          # if it is, there's no previous output directory
            record_mask_changes(mask_sources, displacement_r2['mask'], mask_combined, processing_date, f"{volcano_dir}{processing_date}/", previous_output_dir)             # record any changes in the mask (ie pixels that are now masked due to being incoherent).  
            
            
            # 6c: LiCSAlert stuff
            displacement_r2_combined_current = shorten_LiCSAlert_data(displacement_r2_combined, n_end=ifg_n+1)                        # get the ifgs available for this loop (ie one more is added each time the loop progresses),  +1 as indexing and want to include this data
            cumulative_baselines_current = temporal_baselines['baselines_cumulative'][:ifg_n+1]                                       # also get current time values.  +1 as indexing and want to include this data
            
            sources_tcs_baseline, residual_tcs_baseline = LiCSAlert(sources_mask_combined, cumulative_baselines_current,                                              # the LiCSAlert algoirthm, using the sources with the combined mask (sources_mask_combined)
                                                                displacement_r2_combined_current['incremental'][:(LiCSAlert_settings['baseline_end_ifg_n']+1),],               # baseline ifgs
                                                                displacement_r2_combined_current['incremental'][(LiCSAlert_settings['baseline_end_ifg_n']+1):,],               # monitoring ifgs
                                                                t_recalculate=10, verbose=False)                                                                      # recalculate lines of best fit every 10 acquisitions
        
            LiCSAlert_figure(sources_tcs_baseline, residual_tcs_baseline, sources_mask_combined, displacement_r2_combined_current, LiCSAlert_settings['baseline_end_ifg_n'],  # creat the LiCSAlert figure
                             cumulative_baselines_current, out_folder = f"{volcano_dir}{processing_date}", day0_date = temporal_baselines['acq_dates'][0])    #
            
            
        sys.stdout = original                                                                                                                       # return stdout to be normal.  
        f_run_log.close()                                                                                                                                   # and close the log file.  


#%%
def LiCSAlert_dates_status(LiCSAlert_required_dates, LiCSAlert_dates, folder_LiCSAlert):
    """ Given a list of dates in which LiCSAlert has been run, check that the required outputs are present in each folder.  
    Inputs:
        dates | list of strings | dates that LiCSAlert was run until.  In form YYYYMMDD
    Returns:
        dates_incomplete | list of strings | dates that a LiCSAlert folder exisits, but it doesn't have all the ouptuts.  
    History:
        2020/11/13 | MEG | Written
        2020_11_17 | MEG | Overhauled ready for version 2
    """
    from pathlib import Path
    import os
    import fnmatch                                                                      # used to compare lists and strings using wildcards
    
    constant_outputs = ['mask_changes_graph.png',                                      # The output files that are expected to exist and never change name
                        'mask_changes.png',
                        'mask_history.pkl']
    variable_outputs = ['LiCSAlert_figure_with_*_monitoring_interferograms.png']        # The output files that are expected to exist and change name.  

    # 0: The dates that still need to be processed
    pending  = []
    for LiCSAlert_required_date in LiCSAlert_required_dates:
        if LiCSAlert_required_date not in LiCSAlert_dates:
            pending.append(LiCSAlert_required_date)
        
    # 1: All the dates that haven been processed are either processed, or processed with errors.  Loop through and decide which.  
    processed = []
    processed_with_errors = []
    for LiCSAlert_date in LiCSAlert_dates:
        LiCSAlert_date_folder = Path(folder_LiCSAlert + LiCSAlert_date)                                     # join to make a path to the current LiCSAlert_date folder
        LiCSAlert_date_files = sorted([f.name for f in os.scandir(LiCSAlert_date_folder)])                  # get names of the files in this folder (no paths)
        
        # 1: look for the products that don't change name (ie all but the first)
        all_products_complete = True                                                                                 # initiate as True
        for constant_output in constant_outputs:                                                                     # loop through the outputs that can't change name
            all_products_complete = (all_products_complete) and (constant_output in LiCSAlert_date_files)            # update boolean 
        # 2: look for the product that does change name (the main figure)
        for variable_output in variable_outputs:                                                                     # loop through the outputs that can change name
            output = fnmatch.filter(LiCSAlert_date_files, "LiCSAlert_figure_with_*_monitoring_interferograms.png")   # check for file with wildcard for changing name
            all_products_complete = (all_products_complete) and (len(output) > 0)                                    # empty list if file doesn't exit, use to update boolean
        
        if all_products_complete:
            processed.append(LiCSAlert_date)
        else:
            processed_with_errors.append(LiCSAlert_date)
            
    return processed, processed_with_errors, pending




def run_LiCSAlert_status(folder_ifgs, folder_LiCSAlert, date_baseline_end, LiCSAlert_history_file):
    """ When 'LiCSAlert_monitoring_mode' is run in a folder, this function determines which steps
    need to be done, ranging from nothing (everything is up to date), through to runnings LiCSBAS, ICASAR, and LiCSAlert.  
    
    Inputs:
        folder_ifgs | path | path to LiCSAR ifgs.  
        folder_LiCSAlert | path | path to where LiCSAlert_monitoring_mode is being run.  
    Rerturns:
        LiCSAlert_status | dict | contains: run_LiCSBAS | Boolean | True if LiCSBAS will be required
                                            run_ICASAR | Boolean | True if ICASAR will be required.  
                                            run_LiCSAlert | Boolean | True if required.  

    History:
        2020/02/18 | MEG |  Written        
        2020/06/29 | MEG | Major rewrite to use a folder based structure    
        2020/11/17 | MEG | Write the docs and add compare_two_dates function.  
        2020/11/24 | MEG | Major update to provide more information on status of volcano being processed.  
        2021_04_15 | MEG | Fix bug in how LiCSAR_last_acq was not always last date in LiCSAR_dates, as LiCSAR_dates was not always in chronological order.  

    """
    import os 
    import datetime
    import sys
    from pathlib import Path
    from LiCSAlert_aux_functions import compare_two_dates, LiCSAR_ifgs_to_s1_acquisitions, Tee
    
    def get_LiCSAlert_required_dates(LiCSAR_dates, date_baseline_end):
        """ Given a list of LICSAR_dates, determine which ones are after the baseline stage ended.  
        Inputs:
            LiCSAR_dates | list | Dates of Sentinel-1 acquisitions used to make LiCSAR ifgs in form YYYYMMDD
            baseline_end | string | Date that LiCSAlert baseline stage ended at, in form YYYYMMDD
        Returns:
            LiCSAlert_required_dates | list | Dates for which there should be LiCSAlert outputs, given the current LiCSAR time series.  
        History:
            2020_11_18 | MEG | Written
        """
        from LiCSAlert_aux_functions import compare_two_dates
        LiCSAlert_required_dates = []
        for LiCSAR_date in LiCSAR_dates:
            after_baseline_end = compare_two_dates(date_baseline_end, LiCSAR_date)               # check if the date was after the end of the baseline stage.  
            if after_baseline_end:
                LiCSAlert_required_dates.append(LiCSAR_date)
        return LiCSAlert_required_dates
    

    history_file = open(LiCSAlert_history_file, 'a')                                             # append to the history file
    original = sys.stdout                                                                        # record the orignal form of stdout, to be used to return it to this after using Tee
    sys.stdout = Tee(sys.stdout, history_file)                                                   # class that sends terminal to file and to terminal.  
    now = datetime.datetime.now()                                                                # get the current time, ready for recording
    print(f"\nLiCSAlert is being run for this volcano at {now.strftime('%d/%m/%Y %H:%M:%S')}")   # record the current time in the LiCSAlert history file

    
    # 0: Get the LiCSAR dates.  
    LiCSAR_ifgs = sorted([f.name for f in os.scandir(folder_ifgs) if f.is_dir()])                # get names of folders produced by LiCSAR (ie the ifgs), and keep chronological.  
    if not LiCSAR_ifgs:                                                                          # RR addition.  To check that the list isn't empty
        print(f"No files found in {folder_ifgs} ... ")
        run_LiCSBAS = run_ICASAR = run_LiCSAlert = False                                        # nothing to run, if there are no interferograms available from LiCSAR.  
        return run_LiCSBAS, run_ICASAR, run_LiCSAlert
    LiCSAR_dates = sorted(LiCSAR_ifgs_to_s1_acquisitions(LiCSAR_ifgs))                           # a list of the unique dates that LiCSAR ifgs span.  Sorted so that oldest if first, newest is last.  
    LiCSAR_last_acq = LiCSAR_dates[-1]                                                           # get the date of the last Sentinel-1 acquisition used by LiCSARu    
    
    # 1: Determine if we can run LiCSAlert yet (ie past the baseline stage)
    LiCSAR_past_baseline = compare_two_dates(date_baseline_end, LiCSAR_last_acq)                 # determine if the last LiCSAR date is after the baseline stage has ended.  
    if not LiCSAR_past_baseline:                                                                 # if we haven't passed the baseline stage yet, we can only wait for more Sentinel acquisitions
        print(f"LiCSAR is up to date until {LiCSAR_last_acq}, but the baseline stage is set to end on {date_baseline_end} "
              f" and, as this hasn't been reached yet, LiCSAlert cannot be run yet.")
        run_LiCSBAS = run_ICASAR = run_LiCSAlert = False                                        # nothing to run
        return run_LiCSBAS, run_ICASAR, run_LiCSAlert
    
    # 3: If we are passed the basline stage, determine what dates LiCSAlert has been run for/which need to be run/ which have errors etc.    
    LiCSAlert_required_dates = get_LiCSAlert_required_dates(LiCSAR_dates, date_baseline_end)                # Determine which dates LiCSAlert should have an output for, which are the dates after the baseline has ended (regardless of if we actually have them)
    LiCSAlert_dates = sorted([f.name for f in os.scandir(folder_LiCSAlert) if f.is_dir()])                  # get names of folders produced by LiCSAR (ie the ifgs), and keep chronological.  
    
    # 4: Determine if we need to run ICASAR
    if 'ICASAR_results' not in LiCSAlert_dates:                                                             # the line above will also catch the ICASAR_results folder, and can be used to check if it exists.  
        run_ICASAR = True                                                                                   # if it doesn't exist, it will need to be run.      
    else:                                                                                                   # However, just because the folder exisits, doesn't mean that it has the crucial files (ie it could have crashed)
        ICASAR_files = [f.name for f in os.scandir(Path(folder_LiCSAlert) / 'ICASAR_results')]              # if the folder exists, see what's in it
        if 'ICASAR_results.pkl' in ICASAR_files:
            run_ICASAR = False                                                                                  # if it exists, it will not need to be run
        else:
            run_ICASAR = True                                                                                # The folder exisits, but the important file doesn't, so ICASAR still needs to be run
        
    # 5: Determine which dates have been processed, processed with errors, or processing is pending.  
    for unneeded_folder in ['LiCSBAS', 'ICASAR_results']:                                                   # these folders get caught in the dates list, but aren't dates so need to be deleted.  
        try:
            LiCSAlert_dates.remove(unneeded_folder)                                                         # note that the LiCSBAS folder also gets caught by this, and needs removing as it's not a date.  
        except:
            pass                                                                                            # however, on the first ever run these don't exist.  
    
    processed, processed_with_errors, pending = LiCSAlert_dates_status(LiCSAlert_required_dates, LiCSAlert_dates, folder_LiCSAlert)     # do the determing.  

    if (len(processed_with_errors) > 0) or (len(pending) > 0):                                              # set boolean flags based on results of which dates exist
        run_LiCSBAS = run_LiCSAlert = True
    else:
        run_LiCSBAS = run_LiCSAlert = False

    # 6: Bundle together all outputs to return as a single dict.  
    LiCSAlert_status = {'run_LiCSBAS'             : run_LiCSBAS,                                            # LiCSBAS is run if some dates with errors exist, or if there are pending dates.  Sane as for LiCSAlert
                        'run_ICASAR'              : run_ICASAR,                                             # ICASAR is run unless the folder containing it is found
                        'run_LiCSAlert'           : run_LiCSAlert,                                          # LiCSAlert is run if some dates with errors exist, or if there are pending dates.  Same as for LICSBAS.  
                        'processed_with_errors'   : processed_with_errors,                                  # any dates that have missing LiCSAlert products.  
                        'pending'                 : pending,                                                # new dates to be processed (ie update LiCSBAS time series, run LiCSAlert)
                        'LiCSAR_last_acq'         : LiCSAR_last_acq}                                        # date of the most recent sentinel-1 acquisition that there is an interferogram for.  

    history_file.close()                                                                                    # close the logging file 
    sys.stdout = original                                                                                   # return stdout to be normal (i.e. just to the terminal)

    return LiCSAlert_status
    




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
    """ Given two masks of pixels, create a mask of pixels that are valid for both.  Also return the two sets of data with the new masks applied.  
    Inputs:
        mask_sources | boolean rank 2| original mask
        sources  | r2 array | sources as row vectors
        mask_ifgs | boolean rank 2| new mask
        ifgs  | r2 array | ifgs as row vectors
    Returns:
        ifgs_new_mask
        sources_new_mask
        mask_both | boolean rank 2| original mask
    History:
        2020/02/19 | MEG |  Written      
        2020/06/26 | MEG | Major rewrite.  
        2021_04_20 | MEG | Add check that sources and ifgs are both rank 2 (use row vectors if only one source, but it must be rank2 and not rank 1)
    """
    import numpy as np
    import numpy.ma as ma
    from auxiliary_functions import col_to_ma

        
    
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
    
    
    # check some inputs.  Not exhuastive!
    if (len(sources.shape) != 2) or (len(ifgs.shape) != 2):
        raise Exception(f"Both 'sources' and 'ifgs' must be rank 2 arrays (even if they are only a single source).  Exiting. ")
    
    mask_both = ~np.logical_and(~mask_sources, ~mask_ifgs)                                       # make a new mask for pixels that are in the sources AND in the current time series
    n_pixs_sources = len(np.argwhere(mask_sources == False))                                  # masked pixels are 1s, so invert with 1- bit so that non-masked are 1s, then sum to get number of pixels
    n_pixs_new = len(np.argwhere(mask_ifgs == False))                                          # ditto for new mask
    n_pixs_both = len(np.argwhere(mask_both == False))                                        # ditto for the mutual mask
    print(f"Updating masks and ICA sources.  Of the {n_pixs_sources} in the sources and {n_pixs_new} in the current LiCSBAS time series, "
          f"{n_pixs_both} are in both and can be used in this iteration of LiCSAlert.  ")
    
    ifgs_new_mask = apply_new_mask(ifgs, mask_ifgs, mask_both)                                  # apply the new mask to the old ifgs and return the non-masked elemts as row vectors.  
    sources_new_mask = apply_new_mask(sources, mask_sources, mask_both)                         # ditto for the sources.  
    
    return ifgs_new_mask, sources_new_mask, mask_both
    



#%%
 
def detect_new_ifgs(folder_ifgs, folder_LiCSAlert):
    """ Determine the number of LiCSAR ifgs in a folder, and detect if this changes.  Note that it has different returns, depending on if 
    it is in the simple "initate" mode or not (if not, also returns a flag of if new interferograms have been detected).  
    Inputs:
        folder_ifgs | path | path to LiCSAR ifgs.  
        folder_LiCSAlert | path | path to where LiCSAlert_monitoring_mode is being run.  
    Rerturns:
        new_ifgs_flag | Boolean | True if there are LiCSAR ifgs after the last LiCSAlert product.  
        LiCSAR_last_acq | string | date as string in format YYYYMMDD for the last LiCSAR acquisition.  

    History:
        2020/02/18 | MEG |  Written        
        2020/06/29 | MEG | Major rewrite to use a folder based structure    
        2020/11/17 | MEG | Write the docs and add compare_two_dates function.  

    """
    import os 
    import datetime
    from LiCSAlert_aux_functions import compare_two_dates, LiCSAR_ifgs_to_s1_acquisitions
    
    
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
    
    s1_acquisitions = LiCSAR_ifgs_to_s1_acquisitions(LiCSAR_ifgs)                               # a list of the unique dates that LiCSAR ifgs span.  
    LiCSAR_last_acq = LiCSAR_ifgs[-1][-8:]                                                      # this is the last date that LiCSAR has processed up to, in the form YYYYMMDD
    
    # 1: Get the last date that LiCAlert has been run until
    LiCSAlert_dates = sorted([f.name for f in os.scandir(folder_LiCSAlert) if f.is_dir()])      # get names of folders produced by LiCSAR (ie the ifgs), and keep chronological.  
    for unneeded_folder in ['LiCSBAS', 'ICASAR_results']:                                       # these folders get caught in the dates list, but aren't dates so need to be deleted.  
        try:
            LiCSAlert_dates.remove(unneeded_folder)                                             # note that the LiCSBAS folder also gets caught by this, and needs removing as it's not a date.  
        except:
            pass                                                                                # however, on the first ever run these don't exist.  
    if len(LiCSAlert_dates) == 0:                                                               # if the list is of length 0, LiCSAlert has not been run yet.  
        new_ifgs_flag = True                                                                    # if LiCSAlert hasn't been run yet but the LiCSAR_ifgs is not empty, there must be new interferograms.  
        return new_ifgs_flag, LiCSAR_last_acq, s1_acquisitions                                  # return to parent function.  
    else:
        LiCSAlert_last_run = LiCSAlert_dates[-1]                                                # folder are YYYYMMDD so last one is last time it was run until.  

    # If there are any LiCSAlert_dates:
    if LiCSAlert_dates:
        # Check that LiCSAlert has run succesfully before in each folder
        dates_incomplete = check_LiCSAlert_products(LiCSAlert_dates)
        if len(dates_incomplete) > 0:
            print(f"LiCSBAS products have not be found for {dates_incomplete}")
        # Remove incomplete dates from LiCSAlert_dates:
        LiCSAlert_dates = [i for i in LiCSAlert_dates
                           if i not in dates_incomplete]    
        # If there are any dates left:
        if LiCSAlert_dates:
            LiCSAlert_last_run = LiCSAlert_dates[-1]                                            # folder are YYYYMMDD so last one is last time it was run until.  


    # 2: Compare 0 (the last LiCSAR acquisition) and 1 (the last date LiCSAlert has been run until)
    if len(LiCSAlert_dates) == 0:                                                                       # if there are no LiCSAlert dates, it hasn't been run for this volcano.  
        #print("It appears that LiCSAlert hasn't been run for this volcano yet.  Setting the 'new_ifgs_flag' to True.  ")
        new_ifgs_flag = True
    else:
        new_ifgs_flag = compare_two_dates(LiCSAlert_last_run, LiCSAR_last_acq)              # check if the last LiCSAR acquisition is after the last date LiCSAlert was run until.  
       
    return new_ifgs_flag, LiCSAR_last_acq, s1_acquisitions
    


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
        2020/11/17 | MEG | Add the argument baseline_end to LiCSAlert_settings
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
    LiCSAlert_settings['baseline_end'] = str(config.get('LiCSAlert', 'baseline_end'))                 
    
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
