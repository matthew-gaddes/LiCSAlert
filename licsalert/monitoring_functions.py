#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:09:28 2020

@author: matthew
"""
import pdb

#%%

def LiCSAlert_monitoring_mode(region, volcano, LiCSAlert_pkg_dir, ICASAR_pkg_dir, licsbas_dir, licsalert_dir,
                              override_settings = None):
    """The main function for running LiCSAlert is monitoring mode.  It was designed to work with LiCSAR interferograms, but could be used with 
    any product that creates interfegorams that LiCSBAS can use (which is used for the time series calculation.  )
    
    One of the key functions called in this function is "run_LiCSAlert_status", which returns LiCSAlert status.  
       
    Inputs:
        LiCSAlert_bin | string | Path to folder containing LiCSAlert functions.  
        ICASAR_bin | string | Path to folder containing ICASAR functions.  
        LiCSAR_frames_dir | string | path to the folder containing LiCSAR frames.  Needs trailing /
        LiCSAlert_volcs_dir | string | path to the folder containing each volcano.  Needs trailing /
        override_settings | dict of dicts| Each volcano directory has a .txt file of settings, but by passing settings here, they are applied to all volcanoes.  Mainly useful for 
                                            changing the downsampling factor so that processing is fast and debudding easy.  Contains a dict of icasar settings, and licsalert settings.  
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
    from datetime import datetime
    
    def create_licsalert_outdir(outdir):
        """This directory may already exist (if processing was done but errors occured).  Make it again if so.  
        """
        if not os.path.exists(outdir):                                       # if folder doesn't exist...
            os.mkdir(outdir)                                                 # if doesn't exist, make it                           
        else:                                                                                       # if the folder does already exist
            print(f"The folder {processing_date} appears to exists already.  This is usually due to the date not having all the required LiCSAlert products, and LiCSAlert"
                      f" is now trying to fill this date again.  ")
            shutil.rmtree(outdir)                                        # delete the folder and all its contents
            os.mkdir(outdir)                                             # and remake the folder
    
    def append_licsbas_date_to_file(licsalert_dir, region, volcano, licsbas_json_creation_time):
        """Append the licsbas .json timestamp to the file that records these for each volcano.  """
    
        with open(licsalert_dir / region / volcano / 'licsbas_dates.txt', "a") as licsbas_dates_file:                      # open the file
                #licsbas_dates_file.write(f"{datetime.strftime(licsbas_json_timestamp, '%Y-%m-%d %H:%M:%S')}\n")            # and write this first dummy date to it,
                licsbas_dates_file.write(f"{licsbas_json_creation_time}\n")            # and write this first dummy date to it,
    
    
    if str(ICASAR_pkg_dir) not in sys.path:                                                  # check if already on path
        sys.path.append(str(ICASAR_pkg_dir))                                                 # and if not, add
    
    from licsalert.licsalert import LiCSAlert_preprocessing, LiCSAlert, LiCSAlert_figure, shorten_LiCSAlert_data
    #from monitoring_functions import read_config_file, detect_new_ifgs, update_mask_sources_ifgs, record_mask_changes
    from licsalert.aux import compare_two_dates, Tee
    from licsalert.downsample_ifgs import downsample_ifgs

    
    # 1: Log all outputs to a file for that volcano:
    volcano_dir = licsalert_dir / region / volcano
    f_run_log = open(volcano_dir / "LiCSAlert_history.txt", 'a')                                                           # append to the single txt file for that volcano                             
    original = sys.stdout
    sys.stdout = Tee(sys.stdout, f_run_log)

    print(f"\n\n\n\n\nLiCSAlert is being run for {volcano} at {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}")                              # record the current time in the LiCSAlert history file
        
    # 2: get settings for this volcano
    LiCSAlert_settings, ICASAR_settings = read_config_file(volcano_dir / "LiCSAlert_settings.txt")                                          # read various settings from the volcano's config file
    if override_settings is not None:                                                                                                       # if ther are any override settings
        if 'licsalert' in override_settings.keys():                                                                                         # loook for licsalert ones
            for arg, arg_setting in override_settings['licsalert'].items():                                                                 # loop through them
                LiCSAlert_settings[arg] = arg_setting                                                                                       # and apply
        if 'icasar' in override_settings.keys():                                                                                            # same for icasar settings
            for arg, arg_setting in override_settings['icasar'].items():
                ICASAR_settings[arg] = arg_setting
            
    # 3: Open the licsbas data and determine licsalert status
    displacement_r2, _, tbaseline_info, ref_xy, licsbas_json_creation_time = LiCSBAS_json_to_LiCSAlert(licsbas_dir / region / f"{volcano}.json")          # open the .json LiCSBAS data for this volcano 
    append_licsbas_date_to_file(licsalert_dir, region, volcano, licsbas_json_creation_time)                                                               # append licsbas .json file date to list of file dates used.  
    del displacement_r2['cumulative'], licsbas_json_creation_time                                                                                         #  this is not needed and is deleted for safety.
    displacement_r2 = LiCSAlert_preprocessing(displacement_r2, LiCSAlert_settings['downsample_run'], LiCSAlert_settings['downsample_plot'])           # mean centre and downsize the data
    LiCSAlert_status = run_LiCSAlert_status(tbaseline_info['acq_dates'], licsalert_dir / region / volcano, LiCSAlert_settings['baseline_end'])        # determine the LiCSAlert status for this volcano (ie do we need to run ICASAR, is the time series up to date etc.  )
    print(f"LiCSAlert status:  1) Run ICASAR: {LiCSAlert_status['run_ICASAR']}   2) Run LiCSAlert: {LiCSAlert_status['run_LiCSAlert']}")
   
    # 4: possible run licsalert
    if LiCSAlert_status['run_LiCSAlert']:                                                                                       # if LiCSAlert will be run...
    
        # 4a: either load or make the ICs (from icasar)
        icasar_sources, icasar_mask, LiCSAlert_settings['baseline_end_ifg_n'] = load_or_create_ICASAR_results(LiCSAlert_status['run_ICASAR'], displacement_r2, tbaseline_info,                       # either load or create the ICASAR sources
                                                                                                              LiCSAlert_settings['baseline_end'],  volcano_dir / "ICASAR_results", ICASAR_settings)
    
        # 4b: Deal with changes to the mask of pixels 
        licsbas_mask = displacement_r2['mask']                                                                                                                              # make a copy of the licsbas mask before it gets overwritten with the new combined mask
        displacement_r2['incremental'], icasar_sources, displacement_r2['mask'] = update_mask_sources_ifgs(icasar_mask, icasar_sources,                                      # comapres the mask for the ICs and the mask for the current ifgs and finds a set of pixels that are present in both
                                                                                                           displacement_r2['mask'], displacement_r2['incremental'])          # the new mask overwrites the mask in displacement_r2
        
        displacement_r2["incremental_downsampled"], displacement_r2["mask_downsampled"] = downsample_ifgs(displacement_r2["incremental"], displacement_r2["mask"],           # create the downsampled plotting ones
                                                                                                          LiCSAlert_settings['downsample_plot'], verbose = False)
        
        # 4c: Main loop to run LiCSAlert for each date that is required
        LiCSAlert_status['combined_processing'] = sorted(list(set(LiCSAlert_status['pending'] + LiCSAlert_status['processed_with_errors'])))             # get a sorted and unique list of all the dates we need LiCSAlert for
        for processing_date in LiCSAlert_status['combined_processing']:                                                        # 
            print(f"Running LiCSAlert for {processing_date}")
            ifg_n = tbaseline_info['acq_dates'].index(processing_date)                                                                  # instead of working in dates, switch this to ifg_n in the sorted interferograms.   
            create_licsalert_outdir(volcano_dir / processing_date)                                                                      # Create a folder (YYYYMMDD) for the outputs.  
            record_mask_changes(icasar_mask, licsbas_mask, displacement_r2['mask'], tbaseline_info['acq_dates'][-1], volcano_dir / processing_date)                # create a figure showing the masks (licsbas, icasar, and combined)
                                                                                                                                                                   # 6c: LiCSAlert 
            displacement_r2_current = shorten_LiCSAlert_data(displacement_r2, n_end=ifg_n+1)                                                                       # get the ifgs available for this loop (ie one more is added each time the loop progresses),  +1 as indexing and want to include this data
            cumulative_baselines_current = tbaseline_info['baselines_cumulative'][:ifg_n+1]                                                                        # also get current time values.  +1 as indexing and want to include this data
                       
            sources_tcs_baseline, residual_tcs_baseline = LiCSAlert(icasar_sources, cumulative_baselines_current,                                                  # the LiCSAlert algoirthm, using the sources with the combined mask (sources_mask_combined)
                                                                    displacement_r2_current['incremental'][:(LiCSAlert_settings['baseline_end_ifg_n']+1),],        # baseline ifgs
                                                                    displacement_r2_current['incremental'][(LiCSAlert_settings['baseline_end_ifg_n']+1):,],        # monitoring ifgs
                                                                    t_recalculate=10, verbose=False)                                                               # recalculate lines of best fit every 10 acquisitions
                                                                                                                                                                   # 6d LiCAlert figure
            LiCSAlert_figure(sources_tcs_baseline, residual_tcs_baseline, icasar_sources, displacement_r2_current, LiCSAlert_settings['baseline_end_ifg_n'],       # creat the LiCSAlert figure
                             cumulative_baselines_current, out_folder = volcano_dir / processing_date, day0_date = tbaseline_info['acq_dates'][0])    
            
            if processing_date == LiCSAlert_status['combined_processing'][-1]:                                                                                      # if it's the last date possible....
                save_licsalert_products(displacement_r2['dem'], sources_tcs_baseline, residual_tcs_baseline, volcano_dir / processing_date,                         # save some useful outputs,
                                        remove_others = True, acq_dates = tbaseline_info['acq_dates'] )                                                             # and clean up by removing this file from other dates (as we only need the last one)
                       
    sys.stdout = original                                                                                                                                      # return stdout to be normal.  
    f_run_log.close()                                                                                                                                          # and close the log file.  


#%%

def create_licsalert_update_list(licsbas_dir, licsalert_dir, template_file_dir, licsbas_json_type = 'unfiltered'):
    """ Given a dir of licsbas products and a dir of licsalert time series, determine which licsalert time series need to be updated.  
    This is done by comparing the modified dates of the licsbas .json files, and the date
    of the licsbas file last used by licsalert which is stored as a .txt file in the directory 
    for each volcano.  
    Inputs:
        licsbas_dir | path | location of directory containing the region directories.  
        licsalert_dir | path | location of directory containing the licsalert region directories.  
        template_file_dir | path | location of template files that are copied in to each volcano's licsalert directory if it doesn't already exist.  
        licsbas_json_type | string | 'filtered', 'unfiltered', or 'both'
        
    Returns:
        volc_names_current | list of dicts | name and region of each volcano that is up to date (ie) curent.  
        volc_names_to_update | list of dicts | names and region of each volcano that need to be updated (ie. LiCSAlert run).  
    History:
        2021_10_01 | MEG | Written
    
    """
    import os
    from os import listdir
    import glob
    from datetime import datetime
    

    def get_required_licsbas_jsons(licsbas_dir, licsbas_json_type = 'unfiltered'):
        """ Given a directory of LiCSBAS .json files (i.e. as produced on Jasmin), return paths
        to all of either the filtered or normal files.  Note that each volcano is expected to be
        in a region (e.g. afrida / pacific island etc.)
        Inputs:
            licsbas_dir | path | location of directory containing the region directories.  
            licsbas_json_type | string | 'filtered', 'unfiltered', or 'both'
        Returns:
            volc_paths | list of strings | either to filtered or normal files
        History:
            2021_10_01 | MEG | Written
            2021_10_25 | MEG | Update to return both filt and non filtered, but no web products.  
        """
        from pathlib import Path
        
        regions = listdir(licsbas_dir)
        volc_path_strs = []                                                                                 # list of paths, but as strings
        volc_paths = []                                                                                     # list of pathlib Paths to be returned

        for region in regions:                                                                              # loop through all regions to get all the .json files.  
            volc_path_strs.extend(glob.glob(str(licsbas_dir / region / '*.json')))                            # contains both _filt.json and .json for web and normal resolution, returns a list which we use to extend the other list
            
        for volc_path_str in volc_path_strs:                                                            # loop through each file
            if '_web' not in volc_path_str:                                                                     # make sure it's not a low resolution _web file
            
                if (licsbas_json_type == 'filtered') or (licsbas_json_type == 'both'):
                    if '_filt.' in volc_path_str:
                        volc_paths.append(Path(volc_path_str))
                        
                if (licsbas_json_type == 'unfiltered') or (licsbas_json_type == 'both'):        
                    if '_filt.' not in volc_path_str:
                        volc_paths.append(Path(volc_path_str))

        return volc_paths

        
    def initiate_licsalert_volcano(licsalert_dir, region, volc, template_file_dir):
        """If a volcano has not been processed with LiCSAlert yet (ie no output directory exists), make the directory and 
        copy in the template files.  
        Inputs:
            licsalert_dir | Path | to licsalert folder
            region | string | volcano region
            volc | string | volcano form, in the licsbas form (ie with track info etc.  )
            template_file_dir | Path | to where the template files are stored.  
        Returns:
            directories and files.  
        History:
            2021_10_22 | MEG | Written
        """
        import os
        import glob
        from shutil import copyfile
        from pathlib import Path
               
        os.makedirs(licsalert_dir / region / volc)                                                              # make the volcanoes directory
        template_files = glob.glob(str(template_file_dir / "*"))                                                # get the paths to each template file
        for template_file in template_files:                                                                    # loop through them
            filename = list(Path(template_file).parts)[-1]                                                      # get jut the name
            copyfile(template_file, licsalert_dir / region / volc / filename)                                   # copy them to the new dir
            
        with open(licsalert_dir / "all_volcs_products" / "licsalert_all_volcs_log.txt", "a") as all_volcs_log:
            message = f"A licsbas .json was found for {volc} but no licsalert directory.  The directory has now been created, and the template files copied in.  "
            print(message)
            all_volcs_log.write(message + '\n')



    # 1: get all the licsbas volcanoes (either filt or not), and loop through seeing how they compare to the files used by licsalert for that volcano    
    licsbas_volc_paths = get_required_licsbas_jsons(licsbas_dir, licsbas_json_type)                                           # get a list of the full paths to all required .json files (filtered and non filtered, but no web)
    
    volc_names_current = []
    volc_names_to_update = []
    
    for licsbas_volc_path in licsbas_volc_paths:                                                            # loop through each volcano to get the times
        region = licsbas_volc_path.parts[-2]
        volc_name = licsbas_volc_path.parts[-1].split('.')[0]                                               # last part of the path is the file name, and drop the .json bit.  

        if os.path.isdir(licsalert_dir / region / volc_name) == False:                                      # if the directory for that volcano doesn't exist
            initiate_licsalert_volcano(licsalert_dir, region, volc_name, template_file_dir)                 # make it and copy in the template files, then update the all_volcs log 

        licsbas_json_time = datetime.fromtimestamp(os.path.getmtime(licsbas_volc_path))                     # and determine when the .json file was modified
        licsbas_json_time = licsbas_json_time.replace(microsecond = 0)                                      # incuding micro seconds just makes it more complicated

        f = open(licsalert_dir / region / volc_name / 'licsbas_dates.txt', "r")                             # open the file
        lines = f.readlines()                                                                               # lines of the file to a list                                                                               
        licsbas_time_used_by_licsalert = datetime.strptime(lines[-1][:-1],  '%Y-%m-%d %H:%M:%S')                   # last line is the date of the .json file last usedd, -1 at the end as \n (linebreak) counts as one character.   
        
        if licsbas_json_time > licsbas_time_used_by_licsalert:                                          # if time from licsbas .json file is after the time of the licsbas .json licsalert recorded as having last used
            volc_names_to_update.append({'name'   : volc_name,                                          # then volcano will need updating
                                         'region' : region})
        else:                                                                                           # else (ie same age, or maybe the licsbas .json is beforet the time recorded as having been used by licsalert)
            volc_names_current.append({'name'   : volc_name,                
                                      'region' : region})

    # 2: Update the user and log file
    message = f"\n\n\n\n\nRunning 'create_licsalert_update_list' at {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}\n    Current volcanoes:\n"
    for volc_name_current in volc_names_current:
        message += f"        {volc_name_current['name']}\n"
    message += f"    Volcanoes to update:\n"
    for volc_name_to_update in volc_names_to_update:
        message += f"        {volc_name_to_update['name']}\n"
    print(message)
    with open(licsalert_dir / "all_volcs_products" / "licsalert_all_volcs_log.txt", "a") as all_volcs_log:
        all_volcs_log.write(message)


    return volc_names_current, volc_names_to_update



#%%

def LiCSBAS_json_to_LiCSAlert(json_file):
    """Given a licsbas .json file (as produced by the processing on Jasmin), extract all the information in it
    in a form that is compatible with LiCSAlert (and ICASAR).  
    Inputs:
        json_file | path | path to file, including extnesion.  
    Returns:
        displacment_r3 | dict | Keys: cumulative, incremental.  Stored as masked arrays.  Mask should be consistent through time/interferograms
                                Also lons and lats, which are the lons and lats of all pixels in the images (ie rank2, and not column or row vectors)    
                                Also Dem, mask
        displacment_r2 | dict | Keys: cumulative, incremental, mask.  Stored as row vectors in arrays.  
                                Also lons and lats, which are the lons and lats of all pixels in the images (ie rank2, and not column or row vectors)    
                                Also Dem, mask
        tbaseline_info | dict| acq_dates : acquisition dates as strings
                              daisy_chain : names of the daisy chain of ifgs, YYYYMMDD_YYYYMMDD
                              baselines : temporal baselines of incremental ifgs
                              baselines_cumluative : cumulative baselines.  
          ref_area | dict | x start x stop etc.  
      History:
          2021_10_05 | MEG | Written
    """
    import json
    import numpy as np
    import numpy.ma as ma
    from datetime import datetime
    import os
    
    def nested_lists_to_numpy(nested_list):
        """when opened, .json files have the data as as nested lists.  This function converts
        those to numpy arrays, and works with any rank of data (but only tested with 1, 2, and 3)
        Inputs:
            nested_list | list of lists etc. | 
        returns:
            data | numpy array | automatically sizes to the corret rank.  Often flipped in the vertical though.  
        History:
            2021_10_04 | MEG | Written
        """
        def dimension_unpacker_recursive(nested_list, dims):
            """Given nested lists, determine how many items are in each list.  Recursive!  
            Inputs:
                nested_list | list of lists (of lists?) \ nested list.  
                dims | empty list | will be filled with number of entres in each dimension
            returns:
                dims | list | number of entries in each dimension.  
            History:
                2021_10_04 | MEG | Written
            """
            dims.append(len(nested_list))                                                           # add the number of entries for this dimension to the list
            if type(nested_list[0]) == list:                                                        # if the 1st item of our list is still a list
                dims = dimension_unpacker_recursive(nested_list[0], dims)                           # then apply the same funtion to this item
            return dims
        
        def data_unpacker_recursive(nested_list, index):
            """Given a nested list of data and the index of the data we want as a tuple, extract it.  
            Inputs:
                nested_list | list in list etc.  |
                index | tuple | e.g. (0,10)  or (5, 16, 45)
            Returns:
                data | single value, could be float or int?  
            History:
                2021_10_04 | MEG | Written.  
            """
            nested_list = nested_list[index[0]]                                             # index our nested list with the first available index
            if type(nested_list) == list:                                                   # if it's still a list (and not just a single number)
                nested_list = data_unpacker_recursive(nested_list, index[1:])               # call the function again on this nested list, but using the next index along 
            return nested_list
        
        dims = []                                                                           # initiate to store the number of entres in each dimension
        dims = dimension_unpacker_recursive(nested_list, dims)                              # find the number of entries in each dimension (and so how many dimensions there are too)
        data = np.zeros(dims)                                                               # initiate to store the data
        for index, _ in np.ndenumerate(data):                                               # loop through each item in the data
            data[index] = data_unpacker_recursive(nested_list, index)                       # and fill it with the correct data item from the nested lists
        return data
    
 
    def rank3_ma_to_rank2(ifgs_r3, consistent_mask = False):
        """A function to take a time series of interferograms stored as a rank 3 array,
        and convert it into the ICA(SAR) friendly format of a rank 2 array with ifgs as
        row vectors, and an associated mask.

        For use with ICA, the mask must be consistent (ie the same pixels are masked throughout the time series).

        Inputs:
            ifgs_r3 | r3 masked array | ifgs in rank 3 format
            consistent_mask | boolean | If True, areas of incoherence are consistent through the whole stack
                                        If false, a consistent mask will be made.  N.b. this step can remove the number of pixels dramatically.
        """

        n_ifgs = ifgs_r3.shape[0]
        # 1: Deal with masking
        mask_coh_water = ifgs_r3.mask                                                                                               #get the mask as a rank 3, still boolean
        if consistent_mask:
            mask_coh_water_consistent = mask_coh_water[0,]                                                                          # if all ifgs are masked in the same way, just grab the first one
        else:
            mask_coh_water_sum = np.sum(mask_coh_water, axis = 0)                                                                   # sum to make an image that shows in how many ifgs each pixel is incoherent
            mask_coh_water_consistent = np.where(mask_coh_water_sum == 0, np.zeros(mask_coh_water_sum.shape),
                                                                          np.ones(mask_coh_water_sum.shape)).astype(bool)           # make a mask of pixels that are never incoherent
        ifgs_r3_consistent = ma.array(ifgs_r3, mask = ma.repeat(mask_coh_water_consistent[np.newaxis,], n_ifgs, axis = 0))          # mask with the new consistent mask

        # 2: Convert from rank 3 to rank 2
        n_pixs = ma.compressed(ifgs_r3_consistent[0,]).shape[0]                                                        # number of non-masked pixels
        ifgs_r2 = np.zeros((n_ifgs, n_pixs))
        for ifg_n, ifg in enumerate(ifgs_r3_consistent):
            ifgs_r2[ifg_n,:] = ma.compressed(ifg)

        return ifgs_r2, mask_coh_water_consistent
    
    def daisy_chain_from_acquisitions(acquisitions):
        """Given a list of acquisiton dates, form the names of the interferograms that would create a simple daisy chain of ifgs.  
        Inputs:
            acquisitions | list | list of acquistiion dates in form YYYYMMDD
        Returns:
            daisy_chain | list | names of daisy chain ifgs, in form YYYYMMDD_YYYYMMDD
        History:
            2020/02/16 | MEG | Written
        """
        daisy_chain = []
        n_acqs = len(acquisitions)
        for i in range(n_acqs-1):
            daisy_chain.append(f"{acquisitions[i]}_{acquisitions[i+1]}")
        return daisy_chain

    
    def baseline_from_names(names_list):
        """Given a list of ifg names in the form YYYYMMDD_YYYYMMDD, find the temporal baselines in days_elapsed
        Inputs:
            names_list | list | in form YYYYMMDD_YYYYMMDD
        Returns:
            baselines | list of ints | baselines in days
        History:
            2020/02/16 | MEG | Documented 
        """
        from datetime import datetime
                
        baselines = []
        for file in names_list:
            master = datetime.strptime(file.split('_')[-2], '%Y%m%d')   
            slave = datetime.strptime(file.split('_')[-1][:8], '%Y%m%d')   
            baselines.append(-1 *(master - slave).days)    
        return baselines


    
    ########## Begin
    
    displacement_r2 = {}
    displacement_r3 = {}
    tbaseline_info = {}
    ref_xy = []
    
    with open(json_file, 'r') as json_data:
        licsbas_data = json.load(json_data)
    
    #import pdb; pdb.set_trace()
    licsbas_json_timestamp = licsbas_data['timestamp']                                                  # this is usually made 10-15 seconds beforet the files final write to disk time
    licsbas_json_creation_time = datetime.fromtimestamp(os.path.getmtime(json_file))                     # which is this time
    licsbas_json_creation_time = licsbas_json_creation_time.replace(microsecond = 0)
    
    print(f"Opening the LiCSBAS .json file with timestamp {licsbas_json_timestamp} and system creation time {licsbas_json_creation_time}.  ")
    
    # 0: Get the reference area.  
    ref_list = licsbas_data['refarea'] 
    ref_xy = {'x_start' : int(ref_list[0]),                                            # convert the correct part of the string to an integer
              'x_stop'   : int(ref_list[1]),
              'y_start'  : int(ref_list[2]),
              'y_stop'   : int(ref_list[3])}
        
    # 1: get the lons and lats of each pixel in the image
    lons_mg, lats_mg = np.meshgrid(licsbas_data['x'], licsbas_data['y'])    
    if lats_mg[0,0] < lats_mg[-1,0]:                                                                            # if top left latitude is less than bottom left latitude
        lats_mg = np.flipud(lats_mg)                                                                            # flip the lats
    displacement_r2['lons'] = lons_mg
    displacement_r2['lats'] = lats_mg
    displacement_r3['lons'] = lons_mg
    displacement_r3['lats'] = lats_mg
    
    # 2: get the mask
    mask_coh_water = 1 - nested_lists_to_numpy(licsbas_data['mask'])                        # flip as vertical always flipped when working with .json files. Also, LiCSAlert uses 1 for area to be masked, this returns opposite (so 1 - to invert)
    
    # 3: get the data, mask, and convert to rank 2 in both cumululative and incremental
    if 'data_raw' in licsbas_data.keys():                                                                   # data can be called one of two things.      
        cumulative_r3 = nested_lists_to_numpy(licsbas_data['data_raw'])                                    # 
    elif 'data_filt' in licsbas_data.keys():
        cumulative_r3 = nested_lists_to_numpy(licsbas_data['data_filt'])                                    # 
    else:
        raise Exception(f"The deformation information is expected to be stored in the .json file as either 'data_raw' or 'data_filt', but neither of these were present so can't continue.  ")
    cumulative_r3 *= 0.001                                                                             # licsbas standard is mm, convert to m
    n_im, length, width = cumulative_r3.shape                                                          # time series size, n_im = n_acquisisions
    mask_coh_water_r3 = np.repeat(mask_coh_water[np.newaxis,], n_im, axis = 0)                         # new version that has expanded to be the same size as the cumulative ifgs
    
    # 3a: Reference the time series
    ifg_offsets = np.nanmean(cumulative_r3[:, ref_xy['y_start']: ref_xy['y_stop'], ref_xy['x_start']: ref_xy['x_stop']], axis = (1,2))                                          # get the offset between the reference pixel/area and 0 for each time
    cumulative_r3 = cumulative_r3 - np.repeat(np.repeat(ifg_offsets[:,np.newaxis, np.newaxis], cumulative_r3.shape[1],  axis = 1), cumulative_r3.shape[2], axis = 2)         # do the correction (first make ifg_offsets teh same size as cumulative).      
    
    # 3b: Deal with masking etc.
    cumulative_r3_ma = ma.array(cumulative_r3, mask=mask_coh_water_r3)                                  # mask the cumulative ifgs (first one should be all 0s), note that this could still have nans in it
    for nan_pixel in np.argwhere(np.isnan(cumulative_r3_ma)):                                           # find any pixels that are nan in this, and iterate over
        mask_coh_water_r3[:, nan_pixel[1], nan_pixel[2]] = 1                                            # and modify the mask so that they are masked for all times
    displacement_r3["cumulative"] = ma.array(cumulative_r3, mask=mask_coh_water_r3)                     # recreate the masked array using the new mask.  

    displacement_r3["incremental"] = np.diff(displacement_r3['cumulative'], axis = 0)                   # difference these to get the incremental ifgs, should be one less in time dimension than previous.  
    if displacement_r3["incremental"].mask.shape == ():                                                 # in the case where no pixels are masked, the mask can disappear
        displacement_r3["incremental"].mask = mask_coh_water_r3[1:]                                        # add a new mask, note that we omit the first one (1:) as we have one less ifg when handling incremental
        
    displacement_r2['cumulative'], displacement_r2['mask'] = rank3_ma_to_rank2(displacement_r3['cumulative'])      # convert from rank 3 to rank 2 and a mask
    displacement_r2['incremental'], _ = rank3_ma_to_rank2(displacement_r3['incremental'])                          # also convert incremental, no need to also get mask as should be same as above
    
    # 4: Get the acquisition dates, then calcualet ifg_names, temporal baselines, and cumulative temporal baselines
    tbaseline_info["acq_dates"] = sorted([''.join(date_hyphen_format.split('-')) for date_hyphen_format in licsbas_data['dates'] ])         # convert from yyy-mm-dd to yyyymmdd
    tbaseline_info["ifg_dates"] = daisy_chain_from_acquisitions(tbaseline_info["acq_dates"])                                                # get teh dates of the incremental ifgs
    tbaseline_info["baselines"] = baseline_from_names(tbaseline_info["ifg_dates"])                                                          # and their temporal baselines
    tbaseline_info["baselines_cumulative"] = np.cumsum(tbaseline_info["baselines"])                                                         # cumulative baslines, e.g. 12 24 36 48 etc
    
    # 5: Try to get the DEM
    try:
        dem = nested_lists_to_numpy(licsbas_data['elev'])                                                 # 
        displacement_r2['dem'] = dem                                                                      # and added to the displacement dict in the same was as the lons and lats
        displacement_r3['dem'] = dem                                                                      # 
    except:
        print(f"Failed to open the DEM from the hgt file for this volcano, but trying to continue anyway.")

    # #################### Debuging Sierra Negra
    # print(f"Lon min: {np.min(lons_mg)} Lon max: {np.max(lons_mg)}")
    # print(f"Lat min: {np.min(lats_mg)} Lon max: {np.max(lats_mg)}")
    # from icasar.aux import r2_arrays_to_googleEarth
    # r2_arrays_to_googleEarth(dem[np.newaxis, :,:], lons_mg, lats_mg, layer_name_prefix = 'layer', kmz_filename = 'ICs', out_folder = './')
    # r2_arrays_to_googleEarth(displacement_r3['incremental'][10:11,], lons_mg, lats_mg, layer_name_prefix = 'layer', kmz_filename = 'ifg', out_folder = './')
    # import pdb; pdb.set_trace()
    return displacement_r2, displacement_r3, tbaseline_info, ref_xy, licsbas_json_creation_time


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
        2021_10_20 | MEG | reduce number of mask figures sought.  
    """
    from pathlib import Path
    import os
    import fnmatch                                                                      # used to compare lists and strings using wildcards
    
    constant_outputs = ['mask_status.png']                                      # The output files that are expected to exist and never change name

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
        LiCSAlert_date_folder = folder_LiCSAlert / LiCSAlert_date                                     # join to make a path to the current LiCSAlert_date folder
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


#%%

def run_LiCSAlert_status(licsbas_dates, volcano_path, date_baseline_end):
    """ When 'LiCSAlert_monitoring_mode' is run in a directory, this function determines which steps
    need to be done, ranging from nothing (everything is up to date), through to runnings LiCSBAS, ICASAR, and LiCSAlert.  
    
    Inputs:
        licsbas_dates | list |  ['YYYYMMDD', 'YYYYMMDD' etc.]  of each Sentinel-1 acquisition.  
        volcano_path | Path | path to the licsalert directory for the current volcano.  
        date_baseline_end | string | YYYYMMDD of when to end the baseline stage.  Needed to determine if licsbas data goes past this (and we can start licsalert)
        
    Rerturns:
        LiCSAlert_status | dict | contains: 
                                            run_ICASAR | Boolean | True if ICASAR will be required.  
                                            run_LiCSAlert | Boolean | True if required.  
                                            'processed_with_errors'   : processed_with_errors,                                  # any dates that have missing LiCSAlert products.  
                                            'pending'                 : pending,                                                # new dates to be processed (ie update LiCSBAS time series, run LiCSAlert)
                                            'licsbas_last_acq'        : licsbas_last_acq}                                        # date of the most recent sentinel-1 acquisition that there is an interferogram for.  

    History:
        2020/02/18 | MEG |  Written        
        2020/06/29 | MEG | Major rewrite to use a folder based structure    
        2020/11/17 | MEG | Write the docs and add compare_two_dates function.  
        2020/11/24 | MEG | Major update to provide more information on status of volcano being processed.  
        2021_04_15 | MEG | Fix bug in how LiCSAR_last_acq was not always last date in LiCSAR_dates, as LiCSAR_dates was not always in chronological order.  
        2021_10_05 | MEG | Update to version 2 to work with .json files.  
        2021_10_13 | MEG | Simplify, comment out logging to a .txt file.  

    """
    import os 
    import datetime
    from datetime import datetime as dt
    import sys
    from pathlib import Path
    from licsalert.aux import compare_two_dates, LiCSAR_ifgs_to_s1_acquisitions, Tee
    
    def get_LiCSAlert_required_dates(LiCSAR_dates, date_baseline_end):
        """ Given a list of dates, determine which ones are after the baseline stage ended.  
        Inputs:
            LiCSAR_dates | list | Dates of Sentinel-1 acquisitions used to make LiCSAR ifgs in form YYYYMMDD
            baseline_end | string | Date that LiCSAlert baseline stage ended at, in form YYYYMMDD
        Returns:
            LiCSAlert_required_dates | list | Dates for which there should be LiCSAlert outputs, given the current LiCSAR time series.  
        History:
            2020_11_18 | MEG | Written
        """
        from licsalert.aux import compare_two_dates
        LiCSAlert_required_dates = []
        for LiCSAR_date in LiCSAR_dates:
            after_baseline_end = compare_two_dates(date_baseline_end, LiCSAR_date)               # check if the date was after the end of the baseline stage.  
            if after_baseline_end:
                LiCSAlert_required_dates.append(LiCSAR_date)
        return LiCSAlert_required_dates
    
    # 1: Determine if we can run LiCSAlert yet (ie past the baseline stage)
    licsbas_last_acq = licsbas_dates[-1]
    licsbas_past_baseline = compare_two_dates(date_baseline_end, licsbas_last_acq)                 # determine if the last LiCSAR date is after the baseline stage has ended.  
    if not licsbas_past_baseline:                                                                 # if we haven't passed the baseline stage yet, we can only wait for more Sentinel acquisitions
        print(f"LiCSBAS is up to date until {licsbas_last_acq}, but the baseline stage is set to end on {date_baseline_end} "
              f" and, as this hasn't been reached yet, LiCSAlert cannot be run yet.")
        run_LiCSBAS = run_ICASAR = run_LiCSAlert = False                                        # nothing to run
        return run_LiCSBAS, run_ICASAR, run_LiCSAlert
    
    # 2: If we are past the basline stage, determine what dates LiCSAlert has been run for/which need to be run/ which have errors etc.    
    LiCSAlert_required_dates = get_LiCSAlert_required_dates(licsbas_dates, date_baseline_end)                # Determine which dates LiCSAlert should have an output for, which are the dates after the baseline has ended (regardless of if we actually have them)
    LiCSAlert_dates = sorted([f.name for f in os.scandir(volcano_path) if f.is_dir()])                       # get names of folders produced by LiCSAR (ie the ifgs), and keep chronological.  
    
    # 4: Determine if we need to run ICASAR
    if 'ICASAR_results' not in LiCSAlert_dates:                                                             # the line above will also catch the ICASAR_results folder, and can be used to check if it exists.  
        run_ICASAR = True                                                                                   # if it doesn't exist, it will need to be run.      
    else:                                                                                                   # However, just because the folder exisits, doesn't mean that it has the crucial files (ie it could have crashed)
        ICASAR_files = [f.name for f in os.scandir(volcano_path / 'ICASAR_results')]              # if the folder exists, see what's in it
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
    processed, processed_with_errors, pending = LiCSAlert_dates_status(LiCSAlert_required_dates, LiCSAlert_dates, volcano_path)     # do the determing.  

    if (len(processed_with_errors) > 0) or (len(pending) > 0):                                              # if there are any dates in either pending or processed with errors
        run_LiCSAlert = True                                                                                # 
    else:
        run_LiCSAlert = False  

    # 6: Bundle together all outputs to return as a single dict.  
    LiCSAlert_status = {'run_ICASAR'              : run_ICASAR,                                             # ICASAR is run unless the folder containing it is found
                        'run_LiCSAlert'           : run_LiCSAlert,                                          # LiCSAlert is run if some dates with errors exist, or if there are pending dates.  Same as for LICSBAS.  
                        'processed_with_errors'   : processed_with_errors,                                  # any dates that have missing LiCSAlert products.  
                        'pending'                 : pending,                                                # new dates to be processed (ie update LiCSBAS time series, run LiCSAlert)
                        'licsbas_last_acq'        : licsbas_last_acq}                                        # date of the most recent sentinel-1 acquisition that there is an interferogram for.  
    
    # 7: Print the information to the terminal (and the sys.out file)
    print(f"LiCSBAS date    Past baseline    LiCSAlert ")
    for date_n, date in enumerate(sorted(licsbas_dates)):
        print(f"{date}          ", end = '')
        if dt.strptime(date, '%Y%m%d') >   dt.strptime(date_baseline_end, '%Y%m%d'):
            print('yes          ', end = '')
        else:
            print('no           ', end = '')
        if date in processed:
            print("processed", end = '')
        elif date in processed_with_errors:
            print("processed with errors", end = '')
        elif date in pending:
            print("pending", end = '')
        print('')

    return LiCSAlert_status




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
        2021/10_04 | MEG | Remove LiCSAlert and LiCBAS settings as not needed in new version.  
        2021_10_05 | MEG |Also get the create_all_ifgs_flag
    """
    import configparser    
   
    LiCSAlert_settings = {}
    ICASAR_settings = {}
   
    config = configparser.ConfigParser()                                                       # read the config file
    config.read_file(open(config_file))
    
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
    
    status_string =  config.get('ICASAR', 'create_all_ifgs_flag')            
    if status_string == 'True':
        ICASAR_settings['create_all_ifgs_flag'] = True
    elif status_string == ' False':
        ICASAR_settings['create_all_ifgs_flag'] = True
    else:
        print(f"create_all_ifgs_flag was not understood, so setting it to True.  ")
        ICASAR_settings['create_all_ifgs_flag'] = True
    
    return LiCSAlert_settings, ICASAR_settings


#%%


def load_or_create_ICASAR_results(run_ICASAR, displacement_r2, tbaseline_info, baseline_end, out_dir, ICASAR_settings):
    """
    ICASAR results are always required by LiCSAlert, and these need either to be computed (usually only once at the start),
    or loaded (more common).  
    
    Inputs:
        run_ICASAR | boolean | whether ICASAR should be run, or the results from a previous run loaded.  
        displacment_r2 | dict | interferograms and associated data that are used by the ICASAR algorithm  
        tbaseline_info | dict | various temporal information, such as ifg_dates
        baseline_end | string | YYYYMMDD that the baseline ends on.  
        acq_dates | list of strings | dates of each Sentinel-1 acquisition.  
    Returns:
        icasar_sources | rank 2 array | sources as row vectors
    History:
        2021_10_15 | MEG | Written.  
    """
    from licsalert.aux import get_baseline_end_ifg_n
    import icasar
    from icasar.icasar_funcs import ICASAR
    import pickle
      
    baseline_end_ifg_n = get_baseline_end_ifg_n(tbaseline_info['acq_dates'], baseline_end)                                                            # if this is e.g. 14, the 14th ifg would not be in the baseline stage
    if run_ICASAR:
        print(f"\nRunning ICASAR.")                                      
        
        spatial_ICASAR_data = {'ifgs_dc' : displacement_r2['incremental'][:(baseline_end_ifg_n+1),],                             # only take up to the last 
                               'mask'        : displacement_r2['mask'],
                               'lons'        : displacement_r2['lons'],
                               'lats'        : displacement_r2['lats'],
                               'ifg_dates_dc'   : tbaseline_info['ifg_dates'][:(baseline_end_ifg_n+1)]}                             # ifg dates (yyyymmdd_yyyymmdd), but only up to the end of the baseline stage
        if 'dem' in displacement_r2.keys():
            spatial_ICASAR_data['dem'] = displacement_r2['dem']
        
        sources, tcs, residual, Iq, n_clusters, S_all_info, r2_ifg_means  = ICASAR(spatial_data = spatial_ICASAR_data,                                      # Run ICASAR (slow))
                                                                                   out_folder = out_dir, **ICASAR_settings,
                                                                                   ica_verbose = 'short', figures = 'png')
        mask_sources = displacement_r2['mask']                                                                                                              # rename a copy of the mask
        
    else:
        with open(out_dir / "ICASAR_results.pkl", 'rb') as f_icasar:                                                                     # Or open the products from a previous ICASAR run.  
            sources = pickle.load(f_icasar)   
            mask_sources = pickle.load(f_icasar)
            tcs  = pickle.load(f_icasar)    
            source_residuals = pickle.load(f_icasar)    
            Iq_sorted = pickle.load(f_icasar)    
            n_clusters = pickle.load(f_icasar)    
        f_icasar.close()                                                                                                                           
    
    return sources, mask_sources,  baseline_end_ifg_n




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
    from licsalert.aux import col_to_ma

        
    
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

def record_mask_changes(icasar_mask, licsbas_mask, mask_combined, licsbas_date, current_output_dir):
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
    f1.savefig(current_output_dir / "mask_status.png", bbox_inches='tight')
    plt.close(f1)


#%%
def save_licsalert_products(dem, sources_tcs_baseline, residual_tcs_baseline, outdir, remove_others = True, acq_dates = None):
    """ Save the useful LiCSAlert outputs in a pickle.  Also include the DEM.  To stop there being a version of this file in every date,
    it also tries to delete the previous pickle (by looking in all possible acq_dates)
    
    Inputs:
        dem | rank 2 array | the dem
        sources_tcs_baselines | list of dicts | information about each source that LiCSAlert is using.  
        residual_tcs_baseline | list of 1 dict | as above, but for the residual.  
        outdir | Path | where to save the output.  
        acq_dates | list of strings | other acqusition dates that we should try and delete old versions of this file from.  
        remove_others | boolean | If True, try to remove old version of this file from other acquision dates.  
    Returns:
    History:
        2020_10_20 | MEG | Written
    """
    import pickle
    import os
    
    # 1: save the outputs
    with open(outdir / 'licsalert_results.pkl', 'wb') as f:
        pickle.dump(dem, f)
        pickle.dump(sources_tcs_baseline, f)
        pickle.dump(residual_tcs_baseline, f)
    f.close()
      
    # 2: delete the previous  licsalert products file (which is done by trying on all possible dates)
    if remove_others and acq_dates is not None:                                                         # can't delete if we don't have the acquision dates.  
        current_date = outdir.parts[-1]
        acq_dates.remove(current_date)
        for date in acq_dates:
            try:
                os.remove(outdir.parents / date / 'licsalert_results.pkl')
            except:
                pass


