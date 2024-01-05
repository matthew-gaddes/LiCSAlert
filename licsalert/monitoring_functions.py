#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:09:28 2020

@author: matthew
"""
import pdb

#%%

def LiCSAlert_monitoring_mode(outdir, region, volcano,                                              # for naming out directory.  
                              licsbas_dir = None, licsbas_jasmin_dir = None, data_as_arg = None,    # 3 possible ways to pass data to the function.      
                              licsalert_settings = None, icasar_settings = None):                   # dicts of settings.  
    """The main function for running LiCSAlert is monitoring mode.  It was designed to work with LiCSAR interferograms, but could be used with 
    any product that creates interfegorams that LiCSBAS can use (which is used for the time series calculation.  )
    
    One of the key functions called in this function is "run_LiCSAlert_status", which returns LiCSAlert status.  
       
    Inputs:
        outdir | pathlib Path | parent out directory.  
        region | string or None | When processing large numbers of volcanoes, it can be ueful to group them into regions.  
                                  Outputs would then be outdir / region / volcano.  
        volcano | string | Name of volcano, but only used to name out directories so could be anything.  

        
        Three options to pass data to the function.  Choose one of the three methods
            licsbas_dir  | pathlib Path | If you have used  LiCSBAS, simply provide the directory of the LiCSBAS outputs.                          
            licsbas_jasmin_dir | pathlib Path | If using a LiCSBAS time series that was automatically created on Jasmin for volcano monitoring, 
                                                simply give the directory of the .json files.  
            data_as_arg | dict of dicts | displacement_r2: dict_keys(['dem', 'cumulative', 'mask', 'incremental', 'lons', 'lats']) 
                                          tbaseline_info: dict_keys(['acq_dates', 'ifg_dates', 'baselines', 'baselines_cumulative'])
                                          If there are 327 acq dates (as per the example), there are 327 ifg_dates (as these are the short
                                          temporal baseline ifgs joining the acquisitions, )
        
        licsalert_settings | dict | (['baseline_end', 'figure_intermediate', 'figure_type', 'downsample_run', 'downsample_plot', 'residual_type'])
        icasar_settings | dict | (['n_comp', 'bootstrapping_param', 'tsne_param', 'ica_param', 'hdbscan_param', 'sica_tica', 'ifgs_format']) 

    Returns:
        Directory stucture.  
        
    History:
        2020/06/29 | MEG | Written as a script
        2020/07/03 | MEG | Convert to a funtcion
        2020/11/11 | RR | Add n_para argument
        2020/11/16 | MEG | Pass day0_data info to LiCSAlert figure so that x axis is not in terms of days and is instead in terms of dates.  
        2023_11_15 | MEG | Update input arguments.  
     """
    import sys
    import os
    import pickle
    import shutil
    import copy
    from datetime import datetime
    import numpy.ma as ma
    
    
    def create_licsalert_outdir(outdir):
        """This directory may already exist (if processing was done but errors occured).  Make it again if so.  
        """
        if not os.path.exists(outdir):                                       # if folder doesn't exist...
            os.mkdir(outdir)                                                 # if doesn't exist, make it                           
        else:                                                                                       # if the folder does already exist
            print(f"The folder {outdir} appears to exists already.  This is usually "
                  f"due to the date not having all the required LiCSAlert products, "
                  f"and LiCSAlert is now trying to fill this date again.  ")
            shutil.rmtree(outdir)                                        # delete the folder and all its contents
            os.mkdir(outdir)                                             # and remake the folder
    
    def append_licsbas_date_to_file(outdir, region, volcano, licsbas_json_creation_time):
        """Append the licsbas .json timestamp to the file that records these for each volcano.  """
    
        with open(outdir / region / volcano / 'licsbas_dates.txt', "a") as licsbas_dates_file:                      # open the file
                #licsbas_dates_file.write(f"{datetime.strftime(licsbas_json_timestamp, '%Y-%m-%d %H:%M:%S')}\n")            # and write this first dummy date to it,
                licsbas_dates_file.write(f"{licsbas_json_creation_time}\n")            # and write this first dummy date to it,
    

    import licsalert
    from licsalert.data_importing import LiCSBAS_to_LiCSAlert, LiCSBAS_json_to_LiCSAlert
    #from licsalert.data_importing import crop_licsalert_results_in_time
    from licsalert.data_exporting import save_licsalert_aux_data
    from licsalert.licsalert import LiCSAlert_preprocessing, LiCSAlert, shorten_LiCSAlert_data, write_volcano_status, load_or_create_ICASAR_results
    from licsalert.licsalert import licsalert_date_obj
    from licsalert.aux import compare_two_dates, Tee, find_nearest_date
    from licsalert.downsample_ifgs import downsample_ifgs
    from licsalert.plotting import LiCSAlert_figure, LiCSAlert_epoch_figures, LiCSAlert_aux_figures, LiCSAlert_mask_figure
    
    # 1: Log all outputs to a file for that volcano:
    if region == None:
        volcano_dir = outdir / volcano
    else:
        volcano_dir = outdir / region / volcano
    volcano_dir.mkdir(parents=True, exist_ok=True)                                   
        
    f_run_log = open(volcano_dir / "LiCSAlert_history.txt", 'a')                                                                           # append to the single txt file for that volcano                             
    original = sys.stdout
    sys.stdout = Tee(sys.stdout, f_run_log)

    print(f"\n\n\n\n\nLiCSAlert is being run for {volcano} at {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}")                              # record the current time in the LiCSAlert history file
        
    # 2: get settings for this volcano
    try:
        LiCSAlert_settings, ICASAR_settings = read_config_file(volcano_dir / "LiCSAlert_settings.txt")                                          # read various settings from the volcano's config file
    except:
        print(f"Failed to open the config file for this volcano.  Trying to populate the LiCSAlert_settings and ICASAR_settings dictionaries with from the "
              "function input arguments ")
        try:
            LiCSAlert_settings = {}
            for arg, arg_setting in licsalert_settings.items():                                                                 # loop through them
                LiCSAlert_settings[arg] = arg_setting                                                                                       # and apply
            del licsalert_settings
                
            ICASAR_settings = {}
            for arg, arg_setting in icasar_settings.items():
                ICASAR_settings[arg] = arg_setting
            del icasar_settings
            ICASAR_settings['figures'] = LiCSAlert_settings['figure_type']                                                                          # alwyas the same as for LiCSAlert
        except:
            raise Exception(f"Failed to load the settings from arguments passed to the function.  Try checking licsalert_settings and icasar_settings.  Exiting.  ")

    # 3: Open the the input data, which can be in various formats (3 so far), and deal with mean centering (either in time or space)
    if licsbas_jasmin_dir is not None:
        print(f"LiCSAlert is opening a JASMIN COMET Volcano Portal timeseries json file.  ")
        displacement_r2, _, tbaseline_info, ref_xy, licsbas_json_creation_time = LiCSBAS_json_to_LiCSAlert(licsbas_jasmin_dir / region / f"{volcano}.json")          # open the .json LiCSBAS data for this volcano 
        append_licsbas_date_to_file(outdir, region, volcano, licsbas_json_creation_time)                                                        # append licsbas .json file date to list of file dates used (in the text file for each volano)
        del licsbas_json_creation_time                                                                                                                  # delete for safety
        
    elif licsbas_dir is not None:
        print(f"LiCSAlert is opening a LiCSBAS directory.  ")
        displacement_r2, tbaseline_info, _ = LiCSBAS_to_LiCSAlert(licsbas_dir, figures=True,  filtered = False, n_cols=5,                              # open various LiCSBAS products, spatial ones in displacement_r2, temporal ones in tbaseline_info
                                                                  crop_pixels = None, return_r3 = False, ref_area = True, mask_type = 'dem')
    else:
        print(f"LiCSAlert is using data that was passed to the function as an argument  ")
        displacement_r2 = data_as_arg['displacement_r2']
        tbaseline_info = data_as_arg["tbaseline_info"]
        
    try:
        del displacement_r2['cumulative']                                                                                                                 #  this is not needed and is deleted for safety.
        print(f"LiCSAlert removed 'cumulative' from 'displacement_r2' as it expects only the incremental displacements.  ")
    except:
        pass


    displacement_r2 = LiCSAlert_preprocessing(displacement_r2, tbaseline_info, ICASAR_settings['sica_tica'],                                                    #  mean centre in space (sica) or time (tica) and downsize the data.  
                                              LiCSAlert_settings['downsample_run'], LiCSAlert_settings['downsample_plot'])                                      # mixtures_mc and means contains the daisy chain of ifgs mean centered either in space or time, depending on whther sica or tica
    
    # 3b determine if the date provided happens to be an acquisition date (the easy case)
    updated_date = find_nearest_date(LiCSAlert_settings['baseline_end'], tbaseline_info['acq_dates'])
    if updated_date != LiCSAlert_settings['baseline_end']:
        print(f"\nAs the baseline_end date did not lie on an acquisition date, "
              f"it has been updated so that it does.  Previously, it was "
              f"{LiCSAlert_settings['baseline_end']}, but it is now {updated_date}")
    LiCSAlert_settings['baseline_end'] = licsalert_date_obj(updated_date, tbaseline_info['acq_dates'])

    # 4: possible run licsalert
    LiCSAlert_status = run_LiCSAlert_status(tbaseline_info['acq_dates'], volcano_dir, LiCSAlert_settings['baseline_end'].date, LiCSAlert_settings['figure_intermediate'])                       # NOTE volcano_dir used to be licsalert_dir / region / volcano.   determine the LiCSAlert status for this volcano (ie do we need to run ICASAR, is the time series up to date etc.  )
    print(f"LiCSAlert status:  1) Run ICASAR: {LiCSAlert_status['run_ICASAR']}   2) Run LiCSAlert: {LiCSAlert_status['run_LiCSAlert']}")
    

    if LiCSAlert_status['run_LiCSAlert']:                                                                                       # if LiCSAlert will be run...
    
        icasar_sources, icasar_mask, ics_labels = load_or_create_ICASAR_results(LiCSAlert_status['run_ICASAR'], displacement_r2, tbaseline_info,                             # either load or create the ICASAR sources.  
                                                                                LiCSAlert_settings['baseline_end'],  volcano_dir / "ICASAR_results", ICASAR_settings)        # Note that this uses mixtures_mc, which are mean centered in space or time already depending on if sica or tica is being used.  

        # 4b: Deal with changes to the mask of pixels 
        licsbas_mask = displacement_r2['mask']                                                                                                                              # make a copy of the licsbas mask before it gets overwritten with the new combined mask

        displacement_r2['incremental'], icasar_sources, displacement_r2['mask'] = update_mask_sources_ifgs(icasar_mask, icasar_sources,                                      # comapres the mask for the ICs and the mask for the current ifgs and finds a set of pixels that are present in both
                                                                                                           displacement_r2['mask'], displacement_r2['incremental'])          # the new mask overwrites the mask in displacement_r2
        
        displacement_r2["incremental_downsampled"], displacement_r2["mask_downsampled"] = downsample_ifgs(displacement_r2["incremental"], displacement_r2["mask"],           # create the downsampled plotting ones.  Note that these are also downsampled
                                                                                                          LiCSAlert_settings['downsample_plot'], verbose = False)
        
        # 4c: Loop through each monitoring date to run LiCSAlert
        print(f"SHORTENING THE DATA TO PROCESS")
        processing_dates = []
        # processing_dates.extend(LiCSAlert_status['dates_monitoring'].pending)
        # processing_dates.extend(LiCSAlert_status['dates_monitoring'].processed_with_errors)
        processing_dates.extend(LiCSAlert_status['dates_baseline'].pending)
        processing_dates.extend(LiCSAlert_status['dates_baseline'].processed_with_errors)
        
                            
        for processing_date in processing_dates:

            processing_date = licsalert_date_obj(processing_date, tbaseline_info['acq_dates'])
            # -1 from length to make more intuative for user (eg if 10 dates, 9 of 9 is last)
            print(f"Running LiCSAlert for {processing_date.date} (acquisition # {processing_date.acq_n} "
                  f"of {len(tbaseline_info['acq_dates']) - 1})", end = '')                                                 
            if processing_date.dt > LiCSAlert_settings['baseline_end'].dt:
                print(f" : monitoring date")
            else:
                print(f" : baseline date")
            
            create_licsalert_outdir(volcano_dir / processing_date.date)                                                                                                 # Create a folder (YYYYMMDD) for the outputs.  
            LiCSAlert_mask_figure(icasar_mask, licsbas_mask, displacement_r2['mask'], tbaseline_info['acq_dates'][-1], 
                                  volcano_dir / processing_date.date, figure_type = LiCSAlert_settings['figure_type'])                                                                                   # create a figure showing the masks (licsbas, icasar, and combined)
            
            
            # compare datetimes to see if we are passed baseline stage.  
            if processing_date.dt > LiCSAlert_settings['baseline_end'].dt:
                # if not in baseline, run LiCSAlert
                licsalert_result = LiCSAlert(icasar_sources, tbaseline_info['baselines_cumulative'][:processing_date.acq_n+1],                                               
                                             ifgs_baseline = displacement_r2['incremental_mc_space'][:(LiCSAlert_settings['baseline_end'].acq_n),],                             
                                             ifgs_monitoring = displacement_r2['incremental_mc_space'][(LiCSAlert_settings['baseline_end'].acq_n):processing_date.acq_n,], 
                                             mask = displacement_r2['mask'],
                                             t_recalculate=LiCSAlert_settings['t_recalculate'], verbose=False)                                                                                 # recalculate lines of best fit every 10 acquisitions
                sources_tcs, residual_tcs, reconstructions, residuals  = licsalert_result; del licsalert_result
            else:
                # if in baseline, check if LiCSAlert has been run, and run if not  
                if 'sources_tcs' not in locals().keys():
                    licsalert_result = LiCSAlert(icasar_sources, tbaseline_info['baselines_cumulative'],                                               
                                                 ifgs_baseline = displacement_r2['incremental_mc_space'][:(LiCSAlert_settings['baseline_end'].acq_n),],                             
                                                 ifgs_monitoring = displacement_r2['incremental_mc_space'][(LiCSAlert_settings['baseline_end'].acq_n):,], 
                                                 mask = displacement_r2['mask'],
                                                 t_recalculate=LiCSAlert_settings['t_recalculate'], verbose=False)                                                                                 # recalculate lines of best fit every 10 acquisitions
                    sources_tcs, residual_tcs, reconstructions, residuals  = licsalert_result; del licsalert_result
                    

                    
            dayend_date = licsalert_date_obj(tbaseline_info['acq_dates'][-1], tbaseline_info['acq_dates'])                         # the highest x value (time) of the figure, although data doesn't necesarily plot to here

            LiCSAlert_figure(sources_tcs, residual_tcs, icasar_sources, displacement_r2,                                                  # the main licsalert figure
                             figure_date = processing_date, acq_dates = tbaseline_info['acq_dates'], baselines_cs = tbaseline_info['baselines_cumulative'],
                             baseline_end_date = LiCSAlert_settings['baseline_end'],  dayend_date = dayend_date,
                             #baseline_end_date = LiCSAlert_settings['baseline_end'],  dayend_date = None,
                             
                             figure_type = LiCSAlert_settings['figure_type'], figure_out_dir = volcano_dir / processing_date.date,
                             ifg_xpos_scaler = LiCSAlert_settings['inset_ifgs_scaling'])    
            
            
            LiCSAlert_epoch_figures(processing_date, displacement_r2, reconstructions, 
                                    residuals, tbaseline_info, figure_type = LiCSAlert_settings['figure_type'], 
                                    figure_out_dir = volcano_dir / processing_date.date)                                                                              
            save_epoch_data(sources_tcs, residual_tcs, volcano_dir / processing_date.date)                                                             # save info about the time courses
            write_volcano_status(sources_tcs, residual_tcs, ics_labels, volcano_dir / processing_date.date)                                           # write the change in def  / new def text file.  

            
        LiCSAlert_aux_figures(volcano_dir, icasar_sources, displacement_r2['dem'], displacement_r2['mask'])                                                   # also plot the ICS and the DEM
        save_licsalert_aux_data(volcano_dir, displacement_r2, tbaseline_info)

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

def run_LiCSAlert_status(licsbas_dates, volcano_path, date_baseline_end, figure_intermediate):
    """ When 'LiCSAlert_monitoring_mode' is run in a directory, this function determines which steps
    need to be done, ranging from nothing (everything is up to date), through to running, ICASAR, and LiCSAlert.  
    
    Inputs:
        licsbas_dates | list |  ['YYYYMMDD', 'YYYYMMDD' etc.]  of each Sentinel-1 acquisition.  
        volcano_path | Path | path to the licsalert directory for the current volcano.  
        date_baseline_end | string | YYYYMMDD of when to end the baseline stage.  Needed to determine if licsbas data goes past this (and we can start licsalert)
        figure_intermediate | boolean | If true, only the last LiCSAlert date will be processed.  
        
    Rerturns:
        LiCSAlert_status | dict | contains: 
                                            run_ICASAR | Boolean | True if ICASAR will be required.  
                                            run_LiCSAlert | Boolean | True if required.  
                                            'processed_with_errors'   : processed_with_errors,                                  # any dates that have missing LiCSAlert products.  
                                            'not_processed'                 : pending,                                                # not processed
                                            'to process'
                                            'licsbas_last_acq'        : licsbas_last_acq}                                        # date of the most recent sentinel-1 acquisition that there is an interferogram for.  

    History:
        2020/02/18 | MEG |  Written        
        2020/06/29 | MEG | Major rewrite to use a folder based structure    
        2020/11/17 | MEG | Write the docs and add compare_two_dates function.  
        2020/11/24 | MEG | Major update to provide more information on status of volcano being processed.  
        2021_04_15 | MEG | Fix bug in how LiCSAR_last_acq was not always last date in LiCSAR_dates, as LiCSAR_dates was not always in chronological order.  
        2021_10_05 | MEG | Update to version 2 to work with .json files.  
        2021_10_13 | MEG | Simplify, comment out logging to a .txt file.  
        2023_04_03 | MEG | Update with option to only process last date.  

    """
    import os 
    import datetime
    from datetime import datetime as dt
    import sys
    from pathlib import Path
    from licsalert.aux import compare_two_dates, LiCSAR_ifgs_to_s1_acquisitions, Tee
    
    
    
    def determine_baseline_or_monitoring(LiCSAR_dates, date_baseline_end):
        """ Given a list of dates, determine which ones are in the baseline stage
        (before date_baseline_end) and which ones are in monitoring (after date_baseline_end).  
                                                                     .  
        Inputs:
            LiCSAR_dates | list | Dates of Sentinel-1 acquisitions used to make LiCSAR ifgs in form YYYYMMDD
            baseline_end | string | Date that LiCSAlert baseline stage ended at, in form YYYYMMDD
        Returns:
            baseline_monitoring_dates | dict | baseline and monitoring lists.  
        History:
            2020_11_18 | MEG | Written
            2023_11_31 | MEG | Update to divide into baseline and monitoring.  
        """
        from licsalert.aux import compare_two_dates
        baseline_monitoring_dates = {'baseline'   : [],
                                     'monitoring' : []}

        for LiCSAR_date in LiCSAR_dates:
            after_baseline_end = compare_two_dates(date_baseline_end, LiCSAR_date)               # check if the date was after the end of the baseline stage.  
            if after_baseline_end:
                baseline_monitoring_dates['monitoring'].append(LiCSAR_date)
            else:
                baseline_monitoring_dates['baseline'].append(LiCSAR_date)
        return baseline_monitoring_dates
    
    

    def get_existing_licsalert_dates(volcano_path):
        """ Determine which dates at least have directories made by LiCSAlert
        Also check if the ICASAR directory exists and has the pickle of 
        ICASAR results 
        """
        
        dirs = sorted([f.name for f in os.scandir(volcano_path) if f.is_dir()])                       # get names of folders produced by LiCSAR (ie the ifgs), and keep chronological.  
        
        # 1: Determine if we need to run ICASAR
        if 'ICASAR_results' not in dirs:                                                             # the line above will also catch the ICASAR_results folder, and can be used to check if it exists.  
            run_ICASAR = True                                                                                   # if it doesn't exist, it will need to be run.      
        else:                                                                                                   # However, just because the folder exisits, doesn't mean that it has the crucial files (ie it could have crashed)
            ICASAR_files = [f.name for f in os.scandir(volcano_path / 'ICASAR_results')]              # if the folder exists, see what's in it
            if 'ICASAR_results.pkl' in ICASAR_files:
                run_ICASAR = False                                                                                  # if it exists, it will not need to be run
            else:
                run_ICASAR = True                                                                                # The folder exisits, but the important file doesn't, so ICASAR still needs to be run
            
        # 2: Determine which dates have been processed, processed with errors, or processing is pending.  
        for unneeded_folder in ['LiCSBAS', 'ICASAR_results', 'aux_data_figs']:                                    # these folders get caught in the dates list, but aren't dates so need to be deleted.  
            try:
                dirs.remove(unneeded_folder)                                                         # note that the LiCSBAS folder also gets caught by this, and needs removing as it's not a date.  
            except:
                pass                                                                                            # however, on the first ever run these don't exist.  
        
        return dirs, run_ICASAR
    
        
    class licsalert_dates:
        def __init__(self, existing_dates, required_dates, licsalert_dir, baseline_dates = True,):
            self.existing_dates = existing_dates
            self.required_dates = required_dates
            self.baseline_dates = baseline_dates
            self.licsalert_dir = licsalert_dir
            self.figure_intermediate = figure_intermediate
            self.determine_pending_or_processed()                       # if the date has been processed or not                   
            if len(self.processed) > 0:
                self.determine_processed_correctly()                        # if it has been processed, do all products exist.  
            else:
                self.processed_succesfully = []                     # if nothing has been processed, nothing can be succesful
                self.processed_with_errors = []                     # if nothing has been processed, there can't be any errors.  
            
            
        def determine_pending_or_processed(self):
            """ Iterate through the required dates and see if 
            they exist, creating a list of those that do (processed), 
            and those that still need to be (pending)
            """
            self.pending  = []
            self.processed = []
            for required_date in self.required_dates:
                if required_date in self.existing_dates:                        # if the date exists
                    self.processed.append(required_date)                        # lets say it's been processed
                else:
                    self.pending.append(required_date)                          # else procesing is pending.  
                    
        def determine_processed_correctly(self):
            """ For the processed dates, check they have been processed correctly 
            (no missing files)
            """
            import fnmatch

            # these files don't change in name depending on the date
            constant_outputs = ['epoch_images_data.pkl',
                                'mask_status.png',
                                'time_course_info.pkl',
                                'volcano_status.txt',]
            # these fiels do change in name depending on the date
            variable_outputs = ['LiCSAlert_figure_on*.png',
                                '01_cumulative*.png',
                                '02_incremental*.png',
                                '03_reconstruction*.png',
                                '04_residual*.png']

            
            self.processed_with_errors = []
            self.processed_succesfully = []
            
            for processed_date in self.processed:
                date_dir = self.licsalert_dir / processed_date                                     # join to make a path to the current LiCSAlert_date folder
                date_dir_files = sorted([f.name for f in os.scandir(date_dir)])                  # get names of the files in this folder (no paths)
                all_products_complete = True                                                                                 # initiate as True
                
                for constant_output in constant_outputs:                                                                     # loop through the outputs that can't change name
                    all_products_complete = (all_products_complete) and (constant_output in date_dir_files)            # update boolean 
                
                for variable_output in variable_outputs:                                                                     # loop through the outputs that can change name
                    output = fnmatch.filter(date_dir_files, variable_output)                                                # check for file with wildcard for changing name
                    all_products_complete = (all_products_complete) and (len(output) > 0)                                    # empty list if file doesn't exit, use to update boolean
                    
                if all_products_complete:
                    self.processed_succesfully.append(processed_date)
                else:
                    self.processed_with_errors.append(processed_date)

    


    # 1: Determine if we can run LiCSAlert yet (ie past the baseline stage)
    licsbas_last_acq = licsbas_dates[-1]
    licsbas_past_baseline = compare_two_dates(date_baseline_end, licsbas_last_acq)                 # determine if the last LiCSAR date is after the baseline stage has ended.  
    if not licsbas_past_baseline:                                                                 # if we haven't passed the baseline stage yet, we can only wait for more Sentinel acquisitions
        print(f"LiCSBAS is up to date until {licsbas_last_acq}, but the baseline stage is set to end on {date_baseline_end} "
              f" and, as this hasn't been reached yet, LiCSAlert cannot be run yet.")
        
        LiCSAlert_status = {'run_ICASAR'              : False,                                             # ICASAR is run unless the folder containing it is found
                            'run_LiCSAlert'           : False,                                          # LiCSAlert is run if some dates with errors exist, or if there are pending dates.  Same as for LICSBAS.  
                            'processed_with_errors'   : None,                                  # any dates that have missing LiCSAlert products.  
                            'not_processed'           : None,                                                # new dates to be processed (ie update LiCSBAS time series, run LiCSAlert)
                            'create_figure'           : None,
                            'licsbas_last_acq'        : licsbas_last_acq}                                        # date of the most recent sentinel-1 acquisition that there is an interferogram for.  
        return LiCSAlert_status
        
    existing_licsalert_dirs, run_ICASAR = get_existing_licsalert_dates(volcano_path)    
    baseline_monitoring_dates = determine_baseline_or_monitoring(licsbas_dates, date_baseline_end)                # Determine which dates LiCSAlert should have an output for, which are the dates after the baseline has ended (regardless of if we actually have them)
    
    if not figure_intermediate:
        baseline_monitoring_dates = {'baseline'   : [],                                                         # no figures required for the baseline dates
                                     'monitoring' : [sorted(baseline_monitoring_dates['monitoring'])[-1]]}        # figure and processing only required on last date.  
        
    dates_baseline  = licsalert_dates(existing_licsalert_dirs , baseline_monitoring_dates['baseline'],
                                         volcano_path, baseline_dates = True)                                       # get the status for all the dates in the baseline stage (i.e. do we need to make a figure for them)


    
    dates_monitoring  = licsalert_dates(existing_licsalert_dirs , baseline_monitoring_dates['monitoring'],
                                         volcano_path, baseline_dates = False)                                      # get the status for all the dates in the monitoring stage (i.e. do we need to run LiCSALert and make a figure for them)

    # determine if there are any dates that need licsalert or the licsalert figure
    if ((len(dates_monitoring.processed_with_errors) > 0) or 
        (len(dates_monitoring.pending) > 0) or 
        (len(dates_baseline.processed_with_errors) > 0) or 
        (len(dates_baseline.pending) > 0)):                    
        run_LiCSAlert = True                                                                                        
    else:
        run_LiCSAlert = False  


    # 7: Bundle together all outputs to return as a single dict.  
    LiCSAlert_status = {'run_ICASAR'              : run_ICASAR,                                             # ICASAR is run unless the folder containing it is found
                        'run_LiCSAlert'           : run_LiCSAlert,                                          # LiCSAlert is run if some dates with errors exist, or if there are pending dates.  Same as for LICSBAS.  
                        'licsbas_last_acq'        : licsbas_last_acq,                                        # date of the most recent sentinel-1 acquisition that there is an interferogram for.  
                        'dates_baseline'          : dates_baseline,
                        'dates_monitoring'        : dates_monitoring}
    
    # 8: Print the information to the terminal (and the sys.out file)
    print(f"LiCSBAS date    Past baseline    LiCSAlert      ")
    for date_n, date in enumerate(sorted(licsbas_dates)):                                               # iterate over all acquistisions
        print(f"{date}          ", end = '')
        if dt.strptime(date, '%Y%m%d') >   dt.strptime(date_baseline_end, '%Y%m%d'):                     # check if after end of baseline   
            print('yes             ', end = '')
        else:
            print('no              ', end = '')
        
        if (date in dates_baseline.processed_succesfully) or (date in dates_monitoring.processed_succesfully):                                                                           # if it has been processed
            print("processed succesfuly", end = '')
        elif (date in dates_baseline.processed_with_errors) or (date in dates_monitoring.processed_with_errors):                                                             # if it has been processed with errors
            print("processed with errors           ", end = '')
        elif (date in dates_baseline.pending):
            print("pending (figure only)", end = '')
        elif (date in dates_monitoring.pending):                                                                           # if it is pending.  
             print("pending ", end = '')
        else:
             print("not requested", end = '')
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
    status_string =  config.get('LiCSAlert', 'figure_intermediate')            
    if status_string == 'True':
        LiCSAlert_settings['figure_intermediate'] = True
    elif status_string == ' False':
        LiCSAlert_settings['figure_intermediate'] = True
    else:
        print(f"The LiCSAlert setting figure_intermediate was not understood in the config file, so setting it to True.  ")
        LiCSAlert_settings['figure_intermediate'] = True
    LiCSAlert_settings['figure_type'] =  config.get('LiCSAlert', 'figure_type')            
    
    
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
    
    ICASAR_settings['ifgs_format'] =  config.get('ICASAR', 'ifgs_format')            
    

    
    return LiCSAlert_settings, ICASAR_settings


#%%



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
def save_epoch_data(sources_tcs, residual_tcs, outdir):
    """ Save the useful LiCSAlert outputs in a pickle.  Also include the DEM.  
    
    Inputs:
        dem | rank 2 array | the dem
        sources_tcs | list of dicts | information about each source that LiCSAlert is using.  
        residual_tcs | list of 1 dict | as above, but for the residual.  
        outdir | Path | where to save the output.  
        acq_dates | list of strings | other acqusition dates that we should try and delete old versions of this file from.  
        remove_others | boolean | If True, try to remove old version of this file from other acquision dates.  
    Returns:
        
    History:
        2020_10_20 | MEG | Written
        2023_10_19 | MEG | Remove part of function that deletes old file from previous dates.  
        2023_10_25 | MEG | Update names of variables. 
    """
    import pickle
    import os
    
    # 1: save the outputs
    with open(outdir / 'time_course_info.pkl', 'wb') as f:
        pickle.dump(sources_tcs, f)                                        # list with dict for each IC.  each dict contains dict_keys(['cumulative_tc', 'gradient', 'lines', 'sigma', 'distances', 't_recalculate'])
        pickle.dump(residual_tcs, f)                                       # as above, but list is length 1 as there is only one residual.   
    f.close()
      
    # # 2: delete the previous  licsalert products file (which is done by trying on all possible dates)
    # if remove_others and acq_dates is not None:                                                         # can't delete if we don't have the acquision dates.  
    #     current_date = outdir.parts[-1]
    #     acq_dates.remove(current_date)
    #     for date in acq_dates:
    #         try:
    #             os.remove(outdir.parents / date / 'licsalert_results.pkl')
    #         except:
    #             pass


