#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 09:59:48 2023

@author: matthew
"""

import pdb

#%%

def import_insar_data(
        volcano, volcano_dir, region, licsalert_settings, icasar_settings, 
        licsbas_settings,
        licsbas_jasmin_dir = None, licsbas_dir = None, alignsar_dc = None,
        data_as_arg = None
        ):
    """
    Open InSAR data for LiCSAlert.  Data can either be a:
        licsbas_jasmin_dir - the path a a COMET volcano portal .json
        licsbas_dir  - the path to a LiCSBAS time series 
        alignsar_dc - an AlignSAR data cube.  
        data_as_arg - data processed in a different way.  
        
    Inputs:

        volcano | str | name of volcano frame. 
        volcano_dir | Path | outdir for volcano
        region | string | region that volcano lies in.  Optional?  
        
        licsalert_settings | dict
        icasar_settings  | dict 
        licsbas_settings | dict 
        
        licsbas_jasmin_dir | path |  COMET volcano portal
        licsbas_dir | path | 
        alignsar_dc
        data_as_arg | dict |
    
    Returns:
        displacement_r2 | dict, things like displacement, mask, DEM, lons etc.  
        
    History:
        2025_02_21 | MEG | Written.  
    
    """
    
    from licsalert.monitoring_functions import read_config_file

    
    # 3: Open the the input data, which can be in various formats (3 so far), 
    #    And and input settings for that type of data
    
    # 3.1: a JASMIN dir and associated text file of settings.  
    if licsbas_jasmin_dir is not None:
        # read various settings from the volcano's config file

        
        
        print(f"LiCSAlert is opening a JASMIN COMET Volcano Portal timeseries"
              " json file.  ")
        products = LiCSBAS_json_to_LiCSAlert(
            licsbas_jasmin_dir / region / f"{volcano}.json",
            licsbas_settings['crop_side_length'], 
            licsbas_settings['mask_type']
            )          
        
        (displacement_r2, _, tbaseline_info, ref_xy, 
          licsbas_json_creation_time) = products
        
        
        # #pdb.set_trace()
        # ############################# Begin debug of imported data
        # from licsalert.aux import col_to_ma
        # import matplotlib.pyplot as plt   
        # plt.switch_backend('qt5agg')
        # f, ax = plt.subplots(1); ax.matshow(displacement_r2['mask'])
        # f, ax = plt.subplots(1); ax.matshow(
        #     col_to_ma(displacement_r2['cumulative'][-1,],
        #                 displacement_r2['mask']))
        
        
        # # write the data to a file
        # import pickle
        # with open('extracted_data.pkl', "wb") as f:
        #     pickle.dump(displacement_r2, f)
        #     pickle.dump(tbaseline_info, f)
        # ############################# end debug
        
        
        # append licsbas .json file date to list of file dates used 
        # (in the text file for each volano)
        # append_licsbas_date_to_file(outdir, region, volcano, 
        #                             licsbas_json_creation_time)                                                        
        with open(volcano_dir / 'licsbas_dates.txt', "a") as f:
                f.write(f"{licsbas_json_creation_time}\n")     
        
        
        # delete for safety
        del licsbas_json_creation_time         

    # remaining two ways to pass data to function.  
    else:
        # check user has provided inputs.  
        check_required_args(
            licsalert_settings, ['t_recalculate', 'baseline_end', 
                                  'figure_intermediate',  'figure_type', 
                                  'downsample_run',  'downsample_plot', 
                                  'inset_ifgs_scaling'], 'licsalert_settings'
            )    
    
        check_required_args(
            icasar_settings, ['n_pca_comp_start', 'n_pca_comp_stop'],
            'icasar_settings'
            )

        # 3.2: As a LiCSBAS direcotry
        if licsbas_dir is not None:
            # check licsbas_settings are provided.  
            check_required_args(licsbas_settings, ['mask_type', 'filtered'],
                                'licsbas_settings')
            print(f"LiCSAlert is opening a LiCSBAS directory.  ")
            # if there are no licsbas settings, set them to the default.  
            if licsbas_settings == None:
                licsbas_settings = {"filtered"              : False,
                                    "date_start"            : None,
                                    "date_end"              : None,
                                    'mask_type'             : 'licsbas',
                                    'crop_pixels'           : None}
            ts = LiCSBAS_to_LiCSAlert(licsbas_dir, figures=True, n_cols=5,                              
                                      **licsbas_settings)
            displacement_r2, tbaseline_info = ts
            del ts

        elif alignsar_dc is not None:
            print(f"LiCSAlert is opening a an AlignSAR data cube.  ")
            
            ts = AlignSAR_to_LiCSAlert(alignsar_dc, )
            displacement_r2, tbaseline_info = ts
            

        # 3.3: Data processed with users own approach/software.  
        else:
            print(f"LiCSAlert is using data that was passed to the function as an argument  ")
            displacement_r2 = data_as_arg['displacement_r2']
            tbaseline_info = data_as_arg["tbaseline_info"]
            if licsbas_settings is not None:
                print(f"'licsbas_settings' can only be provided if a "
                      f"licsbas_dir is provided as the input data.  As data "
                      f"is being passed to LiCSAlert in a different way, "
                      f"licsbas_settings will be removed.  ")
                del licsbas_settings
                
    # however data is used, these two arguments must agree.  
    icasar_settings['figures'] = licsalert_settings['figure_type']   
        

    # Check data that has been ingested.  
    displacement_r2, tbaseline_info = check_input_data(
        displacement_r2, tbaseline_info
        )
    
    
    return displacement_r2, tbaseline_info
    

#%%

    

#%%

def check_input_data(displacement_r2, tbaseline_info):
    """ A function to perform any checks on the input data (e.g. the presence
    of nans).  
    
    Inputs: 
        displacements_r2 | dict | standard LiCSAlert dict 
        tbaseline_info | dict | standard LiCSAlert dict 
        
    Returns:
        displacements_r2 | dict | standard LiCSAlert dict 
        tbaseline_info | dict | standard LiCSAlert dict 
        
    History:
        2025_01_22 | MEG | Written.  
        
    """
    import numpy as np
    
    # remove cumualtive from the data as that shouldn't be included
    try:
        del displacement_r2['cumulative']                                                                                                                 #  this is not needed and is deleted for safety.
        print(f"LiCSAlert removed 'cumulative' from 'displacement_r2' as it expects "
              f"only the incremental displacements.  ")
    except:
        pass
    
    # check for nans 
    if np.isnan(displacement_r2['incremental']).any():
        raise Exception(
            "'displacement_r2['incremental']' (i.e. the displacments "
            "between each time step) contains nans.  All nans should be "
            " masked.  Exiting.  "
            )
        
    return displacement_r2, tbaseline_info

#%%

def AlignSAR_to_LiCSAlert(alignsar_dc):
    """ Open an AlignSAR datacube and format the data ready for use with 
    LiCSAlert. 
    """

    import numpy as np    
    import numpy.ma as ma
    import h5py as h5
    from datetime import datetime, timedelta
    
    from licsalert.aux import r3_to_r2
    from licsalert.temporal import daisy_chain_from_acquisitions
    
    
    
    # initialise
    displacement_r2 = {}
    tbaseline_info = {}
    
    # Open the .nc (data cube) file 
    cumh5 = h5.File(alignsar_dc, "r")
    #print(cumh5.keys())
    
    # 1: open the cumualtive displacement
    cumulative = cumh5['cum_displacement'][:]
    # debug plot
    # import matplotlib.pyplot as plt
    # f, ax = plt.subplots()
    # ax.matshow(cumulative[-1])
    
    # lats need reversing so N up
    cumulative = np.flip(cumulative, axis = 1)
    # convert from mm to m
    cumulative *= 0.001
    (_, length, width) = cumulative.shape

    # 2: reference time seres - N/A (referenced allready)
    
    # 3: Open the mask and DEM
    # Alignsar is 1 for valid data, 0 for masked
    # switch to LiCSAlert of 1 for masked, 0 for not.  
    mask_AS = 1 - np.flipud(cumh5['mask'][:])
    
    # also mask any pixels that ever go to nan (1 for masked)
    mask_cum = np.isnan(cumulative)                                                                       
    mask_cum = (np.sum(mask_cum, axis = 0) > 0)                                                             

    # combine the masks
    mask = np.logical_or(mask_AS, mask_cum)

    
    # debug plot
    # import matplotlib.pyplot as plt
    # f, ax = plt.subplots(3); ax[0].matshow(mask_AS); ax[1].matshow(mask_cum)
    # ax[2].matshow(mask)
    
    # apply the mask to the DEM
    dem = np.flipud(cumh5['DEM'][:])
    dem_ma = ma.array(dem, mask = mask)
    displacement_r2['dem'] = dem_ma
    
    # 3a apply the mask to the data
    cum_r3 = ma.array(
        cumulative, mask=np.repeat(mask[np.newaxis,], cumulative.shape[0], 0)
        )
    
    # convert the data to rank 2
    r2_data = r3_to_r2(cum_r3)
    displacement_r2['cumulative'] = r2_data['ifgs']
    displacement_r2['mask'] = r2_data['mask']
    del r2_data
    
    # also calculate incremental displacements
    displacement_r2['incremental'] = np.diff(
        displacement_r2['cumulative'], axis = 0
        )
    
    # 5: get the lons and lats of each pixel in the ifgs
    # lon and lat are stored as row vectors so expand to LiCSAlert meshgrids.  
    lons, lats = np.meshgrid(cumh5['lon'][:], cumh5['lat'][:])
    lats = np.flipud(lats)
    displacement_r2['lats'] = lats
    displacement_r2['lons'] = lons
    
    # 6 Get ENU look vector info N/A - not available.  
    
    # 7 Get the time info
    # open as array, then convert to list of ints
    tbaseline_info['baselines_cumulative'] = cumh5['time'][:].astype(float)
     
    # calculate the shortest temporal baselines between these.  
    tbaseline_info['baselines'] = np.diff(
        tbaseline_info['baselines_cumulative']
        ).tolist()
    
    # get the date of the first acquisition (as a string and a datetime) 
    acq0_date_str = str(cumh5['time'].attrs['units']).split(' ')[2]
    acq0_dt = datetime.strptime(acq0_date_str,'%Y-%m-%d')
    
    # and the dates of all subsequent acquisitions (as a list of strings)
    tbaseline_info['acq_dates'] = [
        (acq0_dt + timedelta(i)).strftime('%Y%m%d') for
        i in tbaseline_info['baselines_cumulative']
        ]

    # and also the names of the daisy chain of ifgs.  
    tbaseline_info['ifg_dates'] = daisy_chain_from_acquisitions(
        tbaseline_info['acq_dates']
        )
    
    # Debug: test output sizes
    # for key, value in tbaseline_info.items():
    #     print(f"{key} : {len(value)}")
       
    
    return displacement_r2, tbaseline_info

#%%

def check_required_args(settings_dict, required_inputs, settings_name):
    """ Check that input dictionary contains the required keys.  
    Inputs:
        settings_dict | dict | contains keys and values of the setting.  
        required_inputs | list | keys that must be in dict. 
        settings_name | str | name of dict, to make error messages clearer
    Returns:
        exception if required is missing.  
    History:
        2024_01_26 | MEG | Written.  
    """
    
   
    if settings_dict is None:
        raise Exception(f"'{settings_name}' is None, but should be a dictionary "
                        f"of settings (see the examples).  Exiting.  ")
    
    for required_input in required_inputs:
        if not (required_input in settings_dict.keys()):
            raise Exception(f"'{required_input}' was not found in "
                            f"'{settings_name}', and it is not optional.  "
                            f"Exiting.  ")
    print(f"All the required arguments were found in the dictionary.  ")

#%%


def crop_ts_data_in_time(date_start, date_end,
                         displacement_r2, tbaseline_info):
    """ Crop time series of data (e.g. licsbas time series) in time.  
    Inputs:
        date_start | yyyymmdd | string, inclusive.  
        date_end | yyyymmdd | string, inclusive.  
        displacement_r2 | dict | time series info
        tbaseline_info | dict | time series time info.  
    Returns:
        displacement_r2 | dict | time series info, cropped
        tbaseline_info | dict | time series time info, cropped
    History:
        2024_01_11 | MEG  | Written
    """
    
    from copy import deepcopy
    from licsalert.licsalert import licsalert_date_obj
    
    # convert strings 
    date_start = licsalert_date_obj(date_start, tbaseline_info['acq_dates'])
    date_end = licsalert_date_obj(date_end, tbaseline_info['acq_dates'])
    
    displacement_r2_crop = deepcopy(displacement_r2)
    tbaseline_info_crop = deepcopy(tbaseline_info)
    
    # crop the time series information in time, for products with n_acq -1 entries
    for key in ['incremental', 'incremental_downsampled', 'incremental_mc_space', 'means_space', 
                'incremental_mc_time', 'means_time', 'mixtures_mc', 'means']:
        if key in displacement_r2.keys():
            displacement_r2_crop[key] = displacement_r2_crop[key][date_start.acq_n : date_end.acq_n, ]
            
    # crop the time series information in time, for products with n_acq entries (so +1) 
    for key in ['cumulative']:
        if key in displacement_r2.keys():
            displacement_r2_crop[key] = displacement_r2_crop[key][date_start.acq_n : date_end.acq_n+1, ]
        
    # crop time info, for products with n_acq -1 entries
    for key in ['ifg_dates', 'baselines']:
        tbaseline_info_crop[key] = tbaseline_info[key][date_start.acq_n : date_end.acq_n]
    
    # crop time info, for products with n_acq entries (so +1)
    for key in ['acq_dates', 'baselines_cumulative']:
        tbaseline_info_crop[key] = tbaseline_info[key][date_start.acq_n : date_end.acq_n+1]
        
    return displacement_r2_crop, tbaseline_info_crop



#%%


def crop_licsalert_results_in_time(processing_date,
                                   acq_dates, sources_tcs, residual_tcs,
                                   reconstructions, residuals, 
                                   displacement_r2, tbaseline_info):
    """ Crop some licsalert products in time.  
    Inputs:
        start_date | string YYYYMMDD | Date to crop to. Must be an acquisition date.  
                                        Included!
        end_date | string YYYYMMDD | As above, also included.  
        acq_dates | list of strings YYYYMMDD | Dates of acqusitions.  
    Returns:
        cropped in time deep copies.  
    History:
        2023_12_08 | MEG | Written
    """
    from copy import deepcopy
    from licsalert.licsalert import licsalert_date_obj
    
    # 
    processing_date = licsalert_date_obj(processing_date, acq_dates)
    
    data = crop_ts_data_in_time(acq_dates[0], processing_date.date, displacement_r2, 
                                tbaseline_info)
    displacement_r2_crop, tbaseline_info_crop = data
    
    # make deep copies that will become the cropped outputs.  
    sources_tcs_crop = deepcopy(sources_tcs)
    residual_tcs_crop = deepcopy(residual_tcs)
    reconstructions_crop = deepcopy(reconstructions)
    residuals_crop = deepcopy(residuals)
    
    # crop the time course information in time, +1 to make inclusive (as has 
    # as many entries as there are acquisitions
    for source_tc in sources_tcs_crop:
        for key in ['cumulative_tc', 'distances']:
            source_tc[key] = source_tc[key][:processing_date.acq_n+1, ]
        # as this is a list, need indexing without trailling ,
        source_tc['lines'] = source_tc['lines'][:processing_date.acq_n+1]

    # crop the residual information in time
    for residual_tc in residual_tcs_crop:
        for key in ['cumulative_tc', 'distances']:
            residual_tc[key] = residual_tc[key][:processing_date.acq_n, ]
        # as this is a list, need indexing without trailling ,
        residual_tc['lines'] = residual_tc['lines'][:processing_date.acq_n]
            

    # optional ones to crop. 
    if reconstructions_crop is not None:
        reconstructions_crop = reconstructions_crop[:processing_date.acq_n]                    
    if reconstructions_crop is not None:    
        residuals_crop = residuals_crop[:processing_date.acq_n, ]
    
    return sources_tcs_crop, residual_tcs_crop, reconstructions_crop, residuals_crop, displacement_r2_crop, tbaseline_info_crop


#%%


def open_aux_data(licsalert_dir):
    """Open all the data stored in the pickle files in aux_images_data
    Inputs:
        licsalert_dir | pathlib Path | output directory when LiCSAlert was run.  
    Returns:
        displacement_r2 | dict | contains (['dem', 'mask', 'incremental', 'lons', 'lats', 'E', 'N', 'U', 'incremental_downsampled', 'mask_downsampled'])
        aux_data | dict | ['icasar_sources', 'dem', 'mask'])
    History:
        2023_10_25 | MEG | Written. 
    """
    
    import pickle
    
    with open(licsalert_dir / "aux_data_figs" / 'original_ts_data.pkl', 'rb') as f:
        displacement_r2 = pickle.load(f)
        tbaseline_info = pickle.load(f)
    f.close()
    
    with open(licsalert_dir / "aux_data_figs" / 'aux_images_data.pkl', 'rb') as f:
        aux_data = pickle.load(f)
    f.close()
    
    return displacement_r2, tbaseline_info, aux_data
    

#%% 

def determine_last_licsalert_date(licsalert_dir):
    """ Given a directory of LiCSAlert outputs, find the youngest date.  
    Inputs:
        licsalert_dir | pathlib Path | output directory when LiCSAlert was run.  
    Returns:
        final_date_dir | pathlib Path | directory of youngest licsalert output.  
    History:
        20??_??_?? | MEG | Written
        2024_12_04 | MEG | New and more robust function to get date dirs.  
        
    """
    from glob import glob
    from pathlib import Path
    
    
    def filter_valid_dates(directories):
        """
        Filters a list of directory names, keeping only those that are valid 
        dates in the yyyymmdd format.
        Parameters:
            directories (list): List of directory names (strings).
            
        Returns:
            list: A list of directory names that are valid dates.
        """
        from datetime import datetime
        valid_dates = []
        for directory in directories:
            try:
                # Attempt to parse the directory name as a date (yyyymmdd)
                datetime.strptime(Path(directory).parts[-1], "%Y%m%d")
                valid_dates.append(directory)
            except ValueError:
                # If parsing fails, it's not a valid date
                continue
        return valid_dates
    
    # get all the directories and files in the licsalert dir
    licsalert_items = sorted(glob(str(licsalert_dir / '*')))
    # get only those that are dates.  
    licsalert_items = filter_valid_dates(licsalert_items)
        
    if len(licsalert_items) == 0:
        raise Exception(
            "There are no LiCSAlert date directories (of the form YYYYMMDD), "
            "so there are no LiCSAlert results to open.  Exiting."
            )
   
    #get the last date
    final_date_dir = Path(sorted(licsalert_items)[-1])
    
    return final_date_dir
    
    
#%%

    
def open_tcs(final_date_dir):
    """ Open the time course data.  
    Inputs:
        
    Returns:
        sources_tcs | list of dicts | One item in list for each source, 
                                      each item contains: 
                                          ['cumulative_tc', 
                                           'gradient', 
                                           'lines', 
                                           'sigma', 
                                           'distances', 
                                           't_recalculate'])
    History:
        2023_10_25 | MEG | Written
        2024_01_04 | MEG | move functionality that finds youngest direcotyr 
                            to new function.  
    """
    
    import pickle

    try:
        with open(final_date_dir / 'time_course_info.pkl', 'rb') as f:
            sources_tcs = pickle.load(f)
            residual_tcs = pickle.load(f)
        f.close()
    except:
        raise Exception(
            "Unable to open the time course information for this date.  "
            "Perhaps it's a baseline date so doesn't have enough information "
            "to create the licsalert figure? Exiting  "
            )

    return sources_tcs, residual_tcs


#%%



class ifg_timeseries():
    def __init__(self, mixtures, ifg_dates):
        self.mixtures = mixtures
        self.ifg_dates = ifg_dates
        self.print_timeseries_info()
        self.mean_centre_in_space()
        self.mean_centre_in_time()
        self.baselines_from_names()
                    
    def print_timeseries_info(self):
        print(f"This interferogram timeseries has {self.mixtures.shape[0]} times and {self.mixtures.shape[1]} pixels.  ")
        
    def mean_centre_in_space(self):
        import numpy as np
        self.means_space = np.mean(self.mixtures, axis = 1)
        self.mixtures_mc_space = self.mixtures - self.means_space[:, np.newaxis]
        
    def mean_centre_in_time(self):
        import numpy as np
        self.means_time = np.mean(self.mixtures, axis = 0)
        self.mixtures_mc_time = self.mixtures - self.means_time[np.newaxis, :]
        
        
    def baselines_from_names(self):
        from datetime import datetime, timedelta
        baselines = []
        for file in self.ifg_dates:
            master = datetime.strptime(file.split('_')[-2], '%Y%m%d')
            slave = datetime.strptime(file.split('_')[-1][:8], '%Y%m%d')
            baselines.append(-1 *(master - slave).days)
        self.t_baselines = baselines

#%%



def LiCSBAS_to_LiCSAlert(LiCSBAS_out_folder, filtered = False, figures = False,
                         n_cols=5, crop_pixels = None, mask_type = 'dem', 
                         date_start = None, date_end = None,
                         draw_manual_mask = None):
    """ A function to prepare the outputs of LiCSBAS for use with LiCSALERT. Note that this includes the step of referencing the time series to the reference area selected by LiCSBAS.  
    LiCSBAS uses nans for masked areas - here these are converted to masked arrays.   Can also create three figures: 1) The Full LiCSBAS ifg, and the area
    that it has been cropped to 2) The cumulative displacement 3) The incremental displacement.  

    Inputs:
        h5_file | string | path to h5 file.  e.g. cum_filt.h5
        figures | boolean | if True, make figures
        n_cols  | int | number of columns for figures.  May want to lower if 
                        plotting a long time series
        crop_pixels | tuple | coords to crop images to.  x then y, 00 is top left.  e.g. (10, 500, 600, 900).  
                                x_start, x_stop, y_start, y_stop, No checking that inputted values make sense.  

        mask_type | string | 'dem' or 'licsbas'  If dem, only the pixels masked in the DEM are masked (i.e. pretty much only water.  ).  Note that any pixels that are incoherent in the time series are also masked (ie if incoherent in one ifg, will be masked for all ifgs.  )
                                                If licsbas, the pixels that licsbas thinks should be masked are masked (ie water + incoherent).  Note that the dem masked is added to the licsbas mask to ensure things are consistent, but there should be no change here (every pixel masked in the DEM is also masked in the licsbas mask)
        draw_manual_mask | None | Not currently used here.  

    Outputs:
        displacment_r2 | dict | Keys: cumulative, incremental, mask.  Stored as row vectors in arrays.  
                                Also lons and lats, which are the lons and lats of all pixels in the images (ie rank2, and not column or row vectors)    
                                Also Dem, mask, and  E N U (look vector components in east north up diretcion)
        tbaseline_info | dict| acq_dates : acquisition dates as strings
                                e.g. ['20141231', 20150105',  ]
                               ifg_dates : strings of the daisy chain of ifgs
                               e.g. ['20141231_20150105', ]
                              baselines : temporal baselines of incremental ifgs
                              e.g [6, 6, 6]
                              baselines_cumulatve : relative to  the first 
                              acquistion.  i.e. [0, 6, 12, 18)]
                                      NOTE - not a list, a 1D array of floats!

    2019/12/03 | MEG | Written
    2020/01/13 | MEG | Update depreciated use of dataset.value to dataset[()] when working with h5py files from LiCSBAS
    2020/02/16 | MEG | Add argument to crop images based on pixel, and return baselines etc
    2020/11/24 | MEG | Add option to get lons and lats of pixels.  
    2021/04/15 | MEG | Update lons and lats to be packaged into displacement_r2 and displacement_r3
    2021_04_16 | MEG | Add option to also open the DEM that is in the .hgt file.  
    2021_05_07 | MEG | Change the name of baseline_info to tbaseline_info to be consistent with LiCSAlert
    2021_09_22 | MEG | Add functionality to extract the look vector componenets (ENU files)
    2021_09_23 | MEG | Add option to extract where the LiCSBAS reference area is.  
    2021_09_28 | MEG | Fix cropping option.  
    2021_11_15 | MEG | Use LiCSBAS reference pixel/area information to reference time series.  
    2021_11_17 | MEG | Add funtcionality to work with LiCSBAS bytes/string issue in reference area.  
    2022_02_02 | MEG | Iimprove masking, and add option to use the LiCSbas mask (i.e. pixels licsbas deems are incohere/poorly unwrapped etc.  )
    """

    import h5py as h5
    import numpy as np
    import numpy.ma as ma
    import matplotlib.pyplot as plt
    import os
    import re
    import pathlib
    #from pathlib import Path
    import pdb
    
    from licsalert.aux import col_to_ma, r3_to_r2
    from licsalert.temporal import daisy_chain_from_acquisitions
    from licsalert.temporal import baseline_from_names
    from licsalert.icasar.aux import add_square_plot
    from licsalert.aux import find_nearest_date

    
    def create_lon_lat_meshgrids(corner_lon, corner_lat, post_lon, post_lat, ifg):
        """ Return a mesh grid of the longitudes and latitues for each pixels.  Not tested!
        I think Corner is the top left, but not sure this is always the case
        """
        ny, nx = ifg.shape
        x = corner_lon +  (post_lon * np.arange(nx))
        y = corner_lat +  (post_lat * np.arange(ny))
        xx, yy = np.meshgrid(x,y)
        geocode_info = {'lons_mg' : xx,
                        'lats_mg' : yy}
        return geocode_info

    def get_param_par(mlipar, field):
        """
        Get parameter from mli.par or dem_par file. Examples of fields are;
         - range_samples
         - azimuth_lines
         - range_looks
         - azimuth_looks
         - range_pixel_spacing (m)
         - azimuth_pixel_spacing (m)
         - radar_frequency  (Hz)
        """
        import subprocess as subp        


    
        
    def read_img(file, length, width, dtype=np.float32, endian='little'):
        """
        Read image data into numpy array.
        endian: 'little' or 'big' (not 'little' is regarded as 'big')
        """
        if endian == 'little':
            data = np.fromfile(file, dtype=dtype).reshape((length, width))
        else:
            data = np.fromfile(file, dtype=dtype).byteswap().reshape((length, width))
        return data


    # -1: Check for common argument errors:
    if not isinstance(LiCSBAS_out_folder, pathlib.PurePath):
        raise Exception(f"'LiCSBAS_out_folder' must be a pathlib Path, but instead is a {type(LiCSBAS_out_folder)}. Exiting.  ")
    
    # 0: Work out the names of LiCSBAS folders - not tested exhaustively! 
    LiCSBAS_folders = {}
    LiCSBAS_folders['all'] = os.listdir(LiCSBAS_out_folder)

    # 1: Loop though looking for the TS direcotry    
    for LiCSBAS_folder in LiCSBAS_folders['all']:                                                                   
        if bool(re.match(re.compile('TS_.'), LiCSBAS_folder)):
            LiCSBAS_folders['TS_'] = LiCSBAS_folder
    # and warn if failed to find it
    if ('TS_' not in LiCSBAS_folders):
        raise Exception(f"Unable to find the TS_* directory that "
                        "contain the LiCSBAS results.  Perhaps the LiCSBAS "
                        "directories have unusual names?  Exiting.  ")
    
    # # 2a: Loop though looking for the ifgs directory, 
    # for LiCSBAS_folder in LiCSBAS_folders['all']:                                                                   
    #     if re.match(re.compile('GEOCml.+clip'), LiCSBAS_folder):                                                    
    #             LiCSBAS_folders['ifgs'] = LiCSBAS_folder

    # # 2b: If we haven't found it already
    # if 'ifgs' not in LiCSBAS_folders.keys():                                                                        
    #     for LiCSBAS_folder in LiCSBAS_folders['all']:                                                               
    #         if re.match(re.compile('GEOCml.+'), LiCSBAS_folder):
    #             LiCSBAS_folders['ifgs'] = LiCSBAS_folder

    # # 2c if we haven't found it already
    # if 'ifgs' not in LiCSBAS_folders.keys():                                                                        
    #     for LiCSBAS_folder in LiCSBAS_folders['all']:                                                               
    #         if re.match(re.compile('GEOC'), LiCSBAS_folder):
    #             LiCSBAS_folders['ifgs'] = LiCSBAS_folder
    
    # error message if the directories can't be found
    # elif ('ifgs' not in LiCSBAS_folders):
    #     raise Exception(f"Unable to find the ifgs directory that "
    #                     "contain the LiCSBAS results.  Perhaps the LiCSBAS "
    #                     "directories have unusual names?  Exiting.  ")

    # 1: Open the h5 file with the incremental deformation in.  
    # displacement_r3 = {}                                                                                        # here each image will 1 x width x height stacked along first axis
    displacement_r2 = {}                                                                                        # here each image will be a row vector 1 x pixels stacked along first axis
    tbaseline_info = {}

    if filtered:
        print(f"Opening the LiCSBAS filtered results.  ")
        cumh5 = h5.File(
            LiCSBAS_out_folder / LiCSBAS_folders['TS_'] / 'cum_filt.h5' ,'r'
            )
    else:
        print(f"Opening the LiCSBAS unfiltered results.  ")
        cumh5 = h5.File(
            LiCSBAS_out_folder / LiCSBAS_folders['TS_'] / 'cum.h5' ,'r'
            )
    
    cumulative = cumh5['cum'][()]                                                                               # get cumulative displacements as a rank3 numpy array
    cumulative *= 0.001                                                                                         # LiCSBAS default is mm, convert to m
    (_, length, width) = cumulative.shape

    # pdb.set_trace()

    # 2: Reference the time series
    ref_str = cumh5['refarea'][()] 
    if not isinstance(ref_str, str):                                                                           # ref_str is sometimes a string, sometimes not (dependent on LiCSBAS version perhaps? )
        ref_str = ref_str.decode()                                                                             # assume that if not a string, a bytes object that can be decoded.        
        
    ref_xy = {'x_start' : int(ref_str.split('/')[0].split(':')[0]),                                            # convert the correct part of the string to an integer
              'x_stop' : int(ref_str.split('/')[0].split(':')[1]),
              'y_start' : int(ref_str.split('/')[1].split(':')[0]),
              'y_stop' : int(ref_str.split('/')[1].split(':')[1])}
    
    try:                                                                                                                                                             # reference the time series
        ifg_offsets = np.nanmean(cumulative[:, ref_xy['y_start']: ref_xy['y_stop'], ref_xy['x_start']: ref_xy['x_stop']], axis = (1,2))                              # get the offset between the reference pixel/area and 0 for each time
        cumulative = cumulative - np.repeat(np.repeat(ifg_offsets[:,np.newaxis, np.newaxis], cumulative.shape[1],  axis = 1), cumulative.shape[2], axis = 2)         # do the correction (first make ifg_offsets teh same size as cumulative).      
        print(f"Succesfully referenced the LiCSBAS time series using the pixel/area selected by LiCSBAS.  ")
    except:
        print(f"Failed to reference the LiCSBAS time series - use with caution!  ")
    
    
    #3: Open the mask and the DEM
    mask_licsbas = read_img(LiCSBAS_out_folder / LiCSBAS_folders['TS_'] / 'results' /  'mask', length, width)                   # this is 1 for coherenct pixels, 0 for non-coherent/water, and masked for water
    mask_licsbas = np.logical_and(mask_licsbas, np.invert(np.isnan(mask_licsbas)))                                          # add any nans to the mask (nans become 0)
    mask_licsbas = np.invert(mask_licsbas)                                                                              # invert so that land is 0 (not masked), and water and incoherent are 1 (masked)
    
    
    # v2
    dem = read_img(
        LiCSBAS_out_folder / LiCSBAS_folders['TS_'] / 'results' / 'hgt', 
        length, width
        )
    
    mask_dem = np.isnan(dem)
    
    mask_cum = np.isnan(cumulative)                                                                       # get where masked
    # sum all the pixels in time (dim 0), and if ever bigger than 0 then must
    # be masked at some point so mask.  
    mask_cum = (np.sum(mask_cum, axis = 0) > 0)                                                             
    
    if mask_type == 'dem':
        mask = np.logical_or(mask_dem, mask_cum)                                                        # if dem, mask water (from DEM), and anything that's nan in cumulative (mask_cum)
    elif mask_type == 'licsbas':
        mask = np.logical_or(mask_licsbas, np.logical_or(mask_dem, mask_cum))
    else:
        raise Exception(f"'mask_type' can be either 'dem' or 'licsbas', but "
                        f"not '{mask_type}'.  Exiting.  ")
        
    # mask_r3 = np.repeat(mask[np.newaxis,], cumulative.shape[0], 0)
    dem_ma = ma.array(dem, mask = mask)
    displacement_r2['dem'] = dem_ma                                                                      # and added to the displacement dict in the same was as the lons and lats
    # displacement_r3['dem'] = dem_ma                                                                      # 
    
    
    # 3: Mask the data  
    # displacement_r3["cumulative"] = ma.array(cumulative, mask=mask_r3)                                   # rank 3 masked array of the cumulative displacement
    # displacement_r3["incremental"] = np.diff(displacement_r3['cumulative'], axis = 0)                           # displacement between each acquisition - ie incremental
    # if displacement_r3["incremental"].mask.shape == ():                                                         # in the case where no pixels are masked, the diff operation on the mask collapses it to nothing.  
    #     displacement_r3["incremental"].mask = mask_r3[1:]                                                # in which case, we can recreate the mask from the rank3 mask, but dropping one from the first dimension as incremental is always one smaller than cumulative.  
    # n_im, length, width = displacement_r3["cumulative"].shape                                   

    # if figures:                                                 
    #     ts_quick_plot(displacement_r3["cumulative"], title = 'Cumulative displacements')
    #     ts_quick_plot(displacement_r3["incremental"], title = 'Incremental displacements')

    # mask the images, consistent mask through time.  
    cum_r3 = ma.array(cumulative, mask=np.repeat(mask[np.newaxis,], cumulative.shape[0], 0))

    r2_data = r3_to_r2(cum_r3)
    displacement_r2['cumulative'] = r2_data['ifgs']
    displacement_r2['mask'] = r2_data['mask']

    displacement_r2['incremental'] = np.diff(
        displacement_r2['cumulative'], axis = 0
        )

    # 4: work with the acquisiton dates to produces names of daisy chain ifgs, and baselines
    tbaseline_info["acq_dates"] = cumh5['imdates'][()].astype(str).tolist()                                     # get the acquisition dates
    tbaseline_info["ifg_dates"] = daisy_chain_from_acquisitions(
        tbaseline_info["acq_dates"]
        )
    tbaseline_info["baselines"] = baseline_from_names(
        tbaseline_info["ifg_dates"]
        )
    # calculate the cumulative baselines from the baselines, ensure 0 at start.  
    tbaseline_info["baselines_cumulative"] = np.concatenate(
        (np.zeros((1)), np.cumsum(tbaseline_info["baselines"])), axis = 0
        )
    
    # 5: get the lons and lats of each pixel in the ifgs
    geocode_info = create_lon_lat_meshgrids(cumh5['corner_lon'][()], cumh5['corner_lat'][()], 
                                            cumh5['post_lon'][()], cumh5['post_lat'][()], 
                                            displacement_r2['mask'])
    displacement_r2['lons'] = geocode_info['lons_mg']                                                                                        # add to the displacement dict
    displacement_r2['lats'] = geocode_info['lats_mg']
    # displacement_r3['lons'] = geocode_info['lons_mg']                                                                                        # add to the displacement dict (rank 3 one)
    # displacement_r3['lats'] = geocode_info['lats_mg']

    # 6: Get the E N U files (these are the components of the ground to satellite look vector in east north up directions.  )   
    try:
        for component in ['E', 'N', 'U']:
            # old version
            # look_vector_component = read_img(LiCSBAS_out_folder / LiCSBAS_folders['ifgs'] / f"{component}.geo", length, width)
            # new version
            look_vector_component = read_img(
                (LiCSBAS_out_folder / LiCSBAS_folders['TS_'] / 'results' /
                 f"{component}.geo", length, width)
                )
            displacement_r2[component] = look_vector_component
            # displacement_r3[component] = look_vector_component
    except:
        print(f"Failed to open the E N U files (look vector components), but trying to continue anyway.")
        
    if crop_pixels is not None:
        print(f"Cropping the images in x from {crop_pixels[0]} to {crop_pixels[1]} "
              f"and in y from {crop_pixels[2]} to {crop_pixels[3]} (NB matrix notation - 0,0 is top left).  ")
        
        if figures:
            ifg_n_plot = 1                                                                                      # which number ifg to plot.  Shouldn't need to change.  
            title = f'Cropped region, ifg {ifg_n_plot}'
            fig_crop, ax = plt.subplots()
            fig_crop.canvas.manager.set_window_title(title)
            ax.set_title(title)
            ax.imshow(col_to_ma(displacement_r2['incremental'][ifg_n_plot,:], displacement_r2['mask']),
                                interpolation='none', aspect='auto')                                            # plot the uncropped ifg
        
        # # loop through and crop, determing if r2 or r3 array.  
        # for product in displacement_r3:
        #     if len(displacement_r3[product].shape) == 2:                                                                                  
        #         resized_r2 = displacement_r3[product][crop_pixels[2]:crop_pixels[3], crop_pixels[0]:crop_pixels[1]]               
        #         displacement_r2[product] = resized_r2
        #         displacement_r3[product] = resized_r2
        #     elif len(displacement_r3[product].shape) == 3:                                                                                
        #         resized_r3 = displacement_r3[product][:, crop_pixels[2]:crop_pixels[3], crop_pixels[0]:crop_pixels[1]]            
        #         displacement_r3[product] = resized_r3
        #         displacement_r2[product], displacement_r2['mask'] = rank3_ma_to_rank2(resized_r3)      
        #     else:
        #         pass

        if figures:
            add_square_plot(crop_pixels[0], crop_pixels[1], crop_pixels[2], crop_pixels[3], ax)                 # draw a box showing the cropped region    

    # 7: Crop in time
    # if no values have been provided, set so that doesn't crop in time
    if date_start is None:
        date_start = tbaseline_info['acq_dates'][0]
    else:
        updated_date = find_nearest_date(
            date_start, tbaseline_info['acq_dates']
            )
        if updated_date != date_start:
            print(
                f"The date to crop the LiCSBAS time series to didn't lie on " 
                f"an acquisition date so has been updated from {date_start} "
                f"to {updated_date}"
                )
            date_start = updated_date; del updated_date
        print(
            f"Cropping the LiCSBAS time series to start on {date_start} "
            " (inclusive).  "
            )
        
    if date_end is None:
        date_end = tbaseline_info['acq_dates'][-1]
    else:
        updated_date = find_nearest_date(
            date_end, tbaseline_info['acq_dates']
            )
        if updated_date != date_end:
            print(
                f"The date to crop the LiCSBAS time series to didn't lie on " 
                f"an acquisition date so has been updated from {date_end} "
                f"to {updated_date}"
                )
            date_end = updated_date; del updated_date
        print(
            f"Cropping the LiCSBAS time series to start on {date_start} "
            " (inclusive).  "
            )

        print(f"Cropping the LiCSBAS time series to end on {date_end} (inclusive).  ")

    # do the cropping
    data = crop_ts_data_in_time(date_start, date_end, displacement_r2, 
                                tbaseline_info)
    
    displacement_r2, tbaseline_info = data
          
    return displacement_r2, tbaseline_info

            


#%%
def LiCSBAS_json_to_LiCSAlert(json_file, crop_side_length, mask_type):
    """Given a licsbas .json file (as produced by the processing on Jasmin), 
    extract all the information in it
    in a form that is compatible with LiCSAlert (and ICASAR).  
    Inputs:
        json_file | path | path to file, including extnesion.  
        crop_side_length | None or int | Possibly crop imakes to square of this
                                        side length in km.  
        mask_type | str | either :
                         'licsbas' to mask pixels that are in the licsbas 
                         mask, and those that ever go to nan.  
                         'nan_once' to mask only pixels that are ever nan, 
                         'nan_variable' to just mask the pixels that are nan
                         in the licsbas data (so the mask varies with time)

    Returns:
        displacment_r3 | dict | Keys: cumulative, incremental.  Stored as 
        masked arrays.  Mask should be consistent through time/interferograms
                                Also lons and lats, which are the lons and lats 
                                of all pixels in the images (ie rank2, and not 
                                                             column or row 
                                                             vectors)    
                                Also Dem, mask
        displacment_r2 | dict | Keys: cumulative, incremental, mask.  Stored as
                                row vectors in arrays.  
                                Also lons and lats, which are the lons and lats 
                                of all pixels in the images (ie rank2, and not
                                                             column or row 
                                                             vectors)    
                                Also Dem, mask
        tbaseline_info | dict| acq_dates : acquisition dates as strings
                              daisy_chain : names of the daisy chain of ifgs, 
                              YYYYMMDD_YYYYMMDD
                              baselines : temporal baselines of incremental 
                              ifgs
                              baselines_cumluative : cumulative baselines.  
          ref_area | dict | x start x stop etc.  
      History:
          2021_10_05 | MEG | Written
          2024_11_01 | MEG | Add support for web (downsampled) .json files.  
    """
    import json
    import numpy as np
    import numpy.ma as ma
    from datetime import datetime
    import os
    from copy import deepcopy

    from licsalert.temporal import daisy_chain_from_acquisitions
    from licsalert.temporal import baseline_from_names
    from licsalert.aux import r3_to_r2
    
    def nested_lists_to_numpy(nested_list):
        """when opened, .json files have the data as as nested lists.  This 
        function convertsthose to numpy arrays
        Inputs:
            nested_list | list of lists etc. | 
            
        returns:
            data | numpy array | automatically sizes to the corret rank. 
                                    Often flipped in the vertical though.  
        History:
            2021_10_04 | MEG | Written
            2024_11_01 | MEG | Complete re-write to be faster and handle
                                downsampled (_web_) .json files.  
        """
        
        def replace_nulls(data):
            """ If a list of lists uses 'null' for no data, np.array will
            not be able to handle this, so it must be converted.  
            """
            if isinstance(data, list):
                return [replace_nulls(item) for item in data]
            elif data == 'null' or data is None:
                return np.nan
            else:
                return data
            
        # convert any 'null' entries to something easier for numpy to use
        cleaned_data = replace_nulls(nested_list)
        # nested list to numpy array.  
        data = np.array(cleaned_data,  dtype=np.float64)
        
        return data
    
    
    def create_mask(cumulative_r3, licsbas_mask, mask_type):
        """ Create the mask
        """
        
        
        # 3.1: get the mask from LiCSBAS
        # returns -9e18 for water, 0 for incoherent pixels, and 1 for coherent
        mask_licsbas = nested_lists_to_numpy(licsbas_mask)
        # ######## debug plot
        # import matplotlib.pyplot as plt
        # test = np.where(
        #     mask_licsbas < -1, -1 * np.ones(mask_licsbas.shape), mask_licsbas 
        #                 )
        # f, ax = plt.subplots(1); ax.matshow(test)
        # ####### end debug
        
        # get the below 0 (usually -9e18) pixels (1 for where water and masked)
        mask_water = (mask_licsbas < 0)
        #f, ax = plt.subplots(1); ax.matshow(mask_water)
        
        # get incoherent pixels (output 1 for masked as incoherent)
        # (also includes water pixels)
        mask_coherence = (mask_licsbas == 0)
        
        # 3.2 also mask any pixels that are ever nan:
        # initialise (1 for masked, 0 for not)
        #mask_nan = np.zeros(cumulative_r3.shape[1:])
        mask_nan_r3 = np.zeros(cumulative_r3.shape)
        for epoch_n, ifg in enumerate(cumulative_r3):
            epoch_mask = np.isnan(ifg)
            mask_nan_r3[epoch_n,] = epoch_mask
            #mask_nan = np.logical_or(mask_nan, epoch_mask)
        # reduce to a single mask for all epochs
        mask_nan = np.any(mask_nan_r3, axis = 0)
       
        # 3.3 combine, masked if True in either.  
        if mask_type == 'nan_once':
            mask = np.logical_or.reduce([mask_water, mask_nan])
            mask_r3 = np.repeat(mask[np.newaxis,], n_im, axis = 0)
            print("Masking the LiCSBAS time series based on pixels that go to nan.")
        elif mask_type == 'licsbas':
            mask = np.logical_or.reduce([mask_water, mask_coherence, mask_nan])
            mask_r3 = np.repeat(mask[np.newaxis,], n_im, axis = 0)
            print("Masking the LiCSBAS time series based on the LiCSBAS mask and "
                  "on pixels that go to nan.")
            
        elif mask_type == 'nan_variable':
            print("Masking the LiCSBAS time series so that any pixels that are "
                  "nan on a certain date remain nan on that date (i.e the mask "
                  "changes for each time step).  No other masks (such as water "
                  "or LiCSBAS are used).  ")
            mask = mask_nan_r3
            mask_r3 = mask
            
        else:
            raise Exception(f"'mask_type' must be either 'nan_once', 'licsbas, or "
                            f"'nan_variable, but is {mask_type}.  Exiting.")
            
        ########### debug
        # import matplotlib.pyplot as plt
        # f, axes = plt.subplots(1,4)
        # axes[0].matshow(mask_water)
        # axes[1].matshow(mask_coherence)
        # axes[2].matshow(mask_nan)
        # axes[3].matshow(mask)
        # for ax, title in zip(axes, ['mask_water', 'mask_coherence', 'mask_nans', 
        #                             'mask']):
        #     ax.set_title(title)
        ############# end 
            
        return mask, mask_r3
    
    def calculate_ref_in_time_offset(cumulative_r3, ref_xy):
        """
        """
    
        def find_center(ref_xy):
            """  Find the centre of a reference region (i.e. are to pixel)"""
            # Calculate the center in x and y dimensions
            x_center = round((ref_xy['x_start'] + ref_xy['x_stop']) / 2)
            y_center = round((ref_xy['y_start'] + ref_xy['y_stop']) / 2)
            
            return {'x_center': x_center, 'y_center': y_center}
        
        def nearest_non_nan_index(x_center, y_center, non_nan_indices):
            """  Find which pixel in non_nan_idices is nearest to the centre """
            # Get the center of the specified region
            center = np.array([y_center, x_center])
            
            # Calculate distances to the center for each non-NaN index
            distances = np.linalg.norm(non_nan_indices - center, axis=1)
            
            # Find the index of the nearest non-NaN pixel
            nearest_index = non_nan_indices[np.argmin(distances)]
            
            
            ref_xy = {'x_start' : nearest_index[1],
                      'x_stop' : nearest_index[1] + 1,
                      'y_start' : nearest_index[0],
                      'y_stop' : nearest_index[0] + 1}
            
            return ref_xy
    
        
        # get the offset between the reference pixel/area and 0 for each time
        ifg_offsets = np.nanmean(
            cumulative_r3[:, ref_xy['y_start']: ref_xy['y_stop'],
                          ref_xy['x_start']: ref_xy['x_stop']], axis = (1,2))
        
        
        
        if np.isnan(ifg_offsets).any():
            print("Their appears to be an error in the reference pixel as it is "
                  "nan on one or more acquisition dates.  Trying to find a nearby "
                  "pixel to use.    ")
            
            # Check if there are any NaNs along the time dimension (axis=0)
            nan_mask = np.isnan(cumulative_r3).any(axis=0)
            
            # Find indices where `nan_mask` is False (pixels with no NaNs across time)
            no_nan_indices = np.argwhere(~nan_mask)
            
            if no_nan_indices.shape[0] == 0:
                raise Exception(f"No pixels were found that never go to nan, "
                                "so referencing cannot be done.  Exiting. ")
                        
            # find the centre of the refernce area (or just return the ref pixel)
            ref_center = find_center(ref_xy)

            # find the nearest non-nan pixel to the reference area
            ref_xy_new = nearest_non_nan_index(
                ref_center['x_center'], ref_center['y_center'], no_nan_indices
                )
            
            # update the reference pixel
            print(f"The refence pixel/area has been updated from: \n"
                  f"{ref_xy}"
                  "to: \n"
                  f"{ref_xy_new}")
            ref_xy = ref_xy_new

            # get the offset between the reference pixel/area and 0 for each time
            ifg_offsets = np.nanmean(
                cumulative_r3[:, ref_xy['y_start']: ref_xy['y_stop'],
                              ref_xy['x_start']: ref_xy['x_stop']], axis = (1,2))                                          
            
        return ifg_offsets
    
    # def accumulate_nan_mask(images):
    #     """
    #     Accumulates a boolean mask of where NaN values are in a 4D array of 
    #     images.
    
    #     Parameters:
    #     images (numpy.ndarray): A 4D numpy array of shape (time, x, y, 
    #                                                         channels).
    
    #     Returns:
    #     numpy.ndarray: A 2D boolean mask of shape (x, y) where True indicates 
    #     the presence of a NaN in any of the images.
    #     """
    #     import numpy as np
    #     # Initialize the mask with False values (no NaNs found yet)
    #     nan_mask = np.zeros(images.shape[1:], dtype=bool)
        
    #     # initilise to count nans
    #     nan_count = np.zeros(images.shape[1:])
    
    #     # Iterate over each image and accumulate the NaN mask
    #     for t in range(images.shape[0]):
    #         # union of sets
    #         nan_mask = nan_mask | np.isnan(images[t,])
            
    #         nan_count[np.isnan(images[t,])] += 1
    #     return nan_mask, nan_count

    
    # initiliase
    displacement_r3 = {}
    tbaseline_info = {}
    ref_xy = []
    
    with open(json_file, 'r') as json_data:
        licsbas_data = json.load(json_data)
    
    # this is usually made 10-15 seconds beforet the files final write to disk time
    licsbas_json_timestamp = licsbas_data['timestamp']                                                  
    # which is this time
    licsbas_json_creation_time = datetime.fromtimestamp(os.path.getmtime(json_file))                     
    # remove the miliseconds part
    licsbas_json_creation_time = licsbas_json_creation_time.replace(microsecond = 0)
    
    print(f"Opening the LiCSBAS .json file with timestamp {licsbas_json_timestamp} "
          f"and system creation time {licsbas_json_creation_time}.  ")
    
    # 0: Get the reference area.  
    ref_list = licsbas_data['refarea'] 
    ref_xy = {'x_start' : int(ref_list[0]),                                          
              'x_stop'   : int(ref_list[1]),
              'y_start'  : int(ref_list[2]),
              'y_stop'   : int(ref_list[3])}
        
    # 1: get the lons and lats of each pixel in the image
    lons_mg, lats_mg = np.meshgrid(licsbas_data['x'], licsbas_data['y'])    
    # if top left latitude is less than bottom left latitude
    if lats_mg[0,0] < lats_mg[-1,0]:                                                                            
        # flip the lats
        lats_mg = np.flipud(lats_mg)                                                                            

    displacement_r3['lons'] = lons_mg
    displacement_r3['lats'] = lats_mg
    
   
    
    # 2: get incremental and cumulative, mask, and convert to rank 2 in both 
    # data can be called one of two things.      
    if 'data_raw' in licsbas_data.keys():                                                                   
        cumulative_r3 = nested_lists_to_numpy(licsbas_data['data_raw'])
    elif 'data_filt' in licsbas_data.keys():
        cumulative_r3 = nested_lists_to_numpy(licsbas_data['data_filt'])
    else:
        raise Exception(f"The deformation information is expected to be stored"
                        "in the .json file as either 'data_raw' or 'data_filt'"
                        ", but neither of these were present so can't"
                        "continue.  ")
   
        
    # licsbas standard is mm, convert to m
    cumulative_r3 *= 0.001                                                                             
    n_im, length, width = cumulative_r3.shape                                                          
    

    # calculate the offsets in time (as time series is not refernced in time)
    ifg_offsets = calculate_ref_in_time_offset(cumulative_r3, ref_xy)
    # do the correction 
    cumulative_r3_referenced = cumulative_r3 - np.repeat(
        np.repeat(ifg_offsets[:,np.newaxis, np.newaxis], 
                  cumulative_r3.shape[1],  axis = 1), 
        cumulative_r3.shape[2], axis = 2)         

    
    ######### debug check referencing
    # pdb.set_trace()
    # import matplotlib.pyplot as plt
    # f, ax = plt.subplots(1, 2)
    # ax[0].matshow(cumulative_r3[-1,])
    # ax[1].matshow(cumulative_r3_referenced[-1,])
    ###########
    
    # 4: build the masks for the data (but do not apply it yet)
    mask, mask_r3 = create_mask(cumulative_r3, licsbas_data['mask'], mask_type)
    
    
    # 3.4 and mask the data, same mask at all times
    displacement_r3["cumulative"]  = ma.array(cumulative_r3, mask=mask_r3)
    
    # difference these to get the incremental ifgs, 
    # should be one less in time dimension than previous.  
    displacement_r3["incremental"] = np.diff(displacement_r3['cumulative'], 
                                             axis = 0)                   
    # in the case where no pixels are masked, the mask can disappear
    if displacement_r3["incremental"].mask.shape == ():                                                 
        # add a new mask, note that we omit the first one (1:)
        # as we have one less ifg when handling incremental
        displacement_r3["incremental"].mask = mask_r3[1:]                                        
        
    ######### debug
    # pdb.set_trace()
    # import matplotlib.pyplot as plt
    # f, ax = plt.subplots(1)
    # ax.matshow(displacement_r3['cumulative'][-1,])
    ###########
    

    # 4: Get the acquisition dates, then calcualet ifg_names, temporal baselines, and cumulative temporal baselines
    tbaseline_info["acq_dates"] = sorted([''.join(date_hyphen_format.split('-')) for date_hyphen_format in licsbas_data['dates'] ])         # convert from yyy-mm-dd to yyyymmdd
    tbaseline_info["ifg_dates"] = daisy_chain_from_acquisitions(tbaseline_info["acq_dates"])                                                # get teh dates of the incremental ifgs
    tbaseline_info["baselines"] = baseline_from_names(tbaseline_info["ifg_dates"])                                                          # and their temporal baselines
    # make cumulative baselines, and ensure that 0 at the start.  
    tbaseline_info["baselines_cumulative"] = np.concatenate((np.zeros((1)), np.cumsum(tbaseline_info["baselines"])), axis = 0)                                                            
    
    # 5: Try to get the DEM (simple numpy array, uses 1e-20 for water (April 25))
    try:
        dem = nested_lists_to_numpy(licsbas_data['elev'])                                                 # 
        displacement_r3['dem'] = dem                                                                      # 
    except:
        print(f"Failed to open the DEM from the hgt file for this volcano, but trying to continue anyway.")
        
    # 6: Possibly crop in space
    if crop_side_length != None:
        displacement_r3 =square_crop_r3_data_in_space (displacement_r3, 
                                                       crop_side_length)

    # 7: make r2 data from r3 data.          
    displacement_r2 = {}
    displacement_r2['lons'] = displacement_r3['lons']
    displacement_r2['lats'] = displacement_r3['lats']
    displacement_r2['dem'] = displacement_r3['dem']
    
    # we can only make the r2 data if the mask doesn't vary with time
    if mask_type != 'nan_variable':
        
        # returns a dict of an array of ifgs as row (ifgs) and the mask (mask)
        r2_data = r3_to_r2(displacement_r3['cumulative'])
 
        # unpack
        displacement_r2['cumulative'] = r2_data['ifgs']
        displacement_r2['mask'] = r2_data['mask']
        del r2_data
        
        # same for the incremental, but mask is the same as cumualtive so omit
        r2_data = r3_to_r2(displacement_r3['incremental'])                          
        displacement_r2['incremental'] = r2_data['ifgs']
        del r2_data
        
        
        return displacement_r2, displacement_r3, tbaseline_info, ref_xy, licsbas_json_creation_time

    # if the mask does, just return r3 data
    else:
        return None, displacement_r3, tbaseline_info, ref_xy, licsbas_json_creation_time
    
    
    
    

#%%

def square_crop_r3_data_in_space(displacement_r3, crop_side_length = 20000):
    """ Given r3 data, crop it to be square with a given side length.  
    Not tested if the crop_side_length is larger than the image.  
    
    Inputs:
        displacement_r3 | dict | lons, lats, dem, cumualtive, incremental
        crop_side_length | int | new image size in m
    Returns:
        displacement_r3_crop | dict | as above, but cropped in space
    History:
        2024_02_01 | MEG | Written
        
    """   ############################### debug
    # #%%
    # import matplotlib.pyplot as plt    
    # f, ax = plt.subplots(1); ax.matshow(mask_coh_water)
    
    # plot 10 if the ifgs (gives an error at the end)
    # for i in np.linspace(0, cumulative_r3.shape[0], 10):
    #     f, ax = plt.subplots(1); ax.matshow(cumulative_r3[int(i),])
    #     plt.pause(3)
    #     plt.close()
    
    
    # def accumulate_nan_mask(images):
    #     """
    #     Accumulates a boolean mask of where NaN values are in a 4D array of 
    #     images.
    
    #     Parameters:
    #     images (numpy.ndarray): A 4D numpy array of shape (time, x, y, 
    #                                                         channels).
    
    #     Returns:
    #     numpy.ndarray: A 2D boolean mask of shape (x, y) where True indicates 
    #     the presence of a NaN in any of the images.
    #     """
    #     import numpy as np
    #     # Initialize the mask with False values (no NaNs found yet)
    #     nan_mask = np.zeros(images.shape[1:], dtype=bool)
        
    #     # initilise to count nans
    #     nan_count = np.zeros(images.shape[1:])
    
    #     # Iterate over each image and accumulate the NaN mask
    #     for t in range(images.shape[0]):
    #         # union of sets
    #         nan_mask = nan_mask | np.isnan(images[t,])
            
    #         nan_count[np.isnan(images[t,])] += 1
    #     return nan_mask, nan_count

    # nan_mask, nan_count = accumulate_nan_mask(cumulative_r3)
    # f, ax = plt.subplots(1); ax.matshow(nan_mask)
    # f, ax = plt.subplots(1); im = ax.matshow(nan_count); f.colorbar(im, ax=ax)
    # f, ax = plt.subplots(1); im = ax.matshow(mask_coh_water)
    # #%%
    ############################### end debug
    import numpy as np
    from geopy import distance
    
    print(f"Cropping the LiCSBAS .json data to be a square of side length "
          f"{crop_side_length} m.  ")
    
    # size of current image in pixels
    ny, nx = displacement_r3['lons'].shape
    
    # get the size of the image in metres.  
    image_size = {}
    # size in metres from bottom left to bottom right
    image_size['x'] = int(distance.distance((displacement_r3['lats'][-1,0], 
                                             displacement_r3['lons'][-1,0]),
                                            (displacement_r3['lats'][-1,-1],
                                             displacement_r3['lons'][-1,-1])).meters )      
    
    # and from bottom left to top left
    image_size['y'] = int(distance.distance((displacement_r3['lats'][-1,0],
                                             displacement_r3['lons'][-1,0]),
                                            (displacement_r3['lats'][0,0], 
                                             displacement_r3['lons'][0,0])).meters)              
    
    if (((image_size['x']) < crop_side_length) or 
        (((image_size['y']) < crop_side_length))):
        
         print(f"The crop_side_length ({crop_side_length}) is larger than  "
               f"the image size in either x or y.  Setting it to the smaller "
               f" of the two to continue.  ")
         crop_side_length = np.min([image_size['x'], image_size['y']])
                
    # get the size of a pixel in metres.  
    pixel_size = {}
    pixel_size['x'] = image_size['x'] / displacement_r3['lons'].shape[1]
    pixel_size['y'] = image_size['y'] / displacement_r3['lats'].shape[0]
    
    # determine how many pixels new image will be.          
    ny_new = int((crop_side_length) / pixel_size['y'])
    nx_new = int((crop_side_length) / pixel_size['x'])
    
    # find what the start and stop pixels will be to do the cropping
    y_start = int((ny/2) - (ny_new/2))
    y_stop = int((ny/2) + (ny_new/2))
    
    x_start = int((nx/2) - (nx_new/2))
    x_stop = int((nx/2) + (nx_new/2))

    # crop each item in space
    displacement_r3_crop = {}
    displacement_r3_crop['lons'] = displacement_r3['lons'][y_start:y_stop,
                                                           x_start:x_stop]
    displacement_r3_crop['lats'] = displacement_r3['lats'][y_start:y_stop,
                                                           x_start:x_stop]
    displacement_r3_crop['dem'] = displacement_r3['dem'][y_start:y_stop,
                                                           x_start:x_stop]
    displacement_r3_crop['cumulative'] = displacement_r3['cumulative'][:, y_start:y_stop,
                                                                           x_start:x_stop]
    displacement_r3_crop['incremental'] = displacement_r3['incremental'][:, y_start:y_stop,
                                                                            x_start:x_stop]
    
    return displacement_r3_crop

#%%

# def append_licsbas_date_to_file(outdir, region, volcano, licsbas_json_creation_time):
#     """Append the licsbas .json timestamp to the file that records these for each volcano.  """

#     with open(outdir / region / volcano / 'licsbas_dates.txt', "a") as licsbas_dates_file:                      # open the file
#             #licsbas_dates_file.write(f"{datetime.strftime(licsbas_json_timestamp, '%Y-%m-%d %H:%M:%S')}\n")            # and write this first dummy date to it,
#             licsbas_dates_file.write(f"{licsbas_json_creation_time}\n")            # and write this first dummy date to it,