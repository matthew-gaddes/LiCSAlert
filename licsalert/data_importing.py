#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 09:59:48 2023

@author: matthew
"""

import pdb

#%%


def crop_licsalert_results_in_time(processing_date, acq_dates, sources_tcs, residual_tcs,
                                   reconstructions, residuals, 
                                   displacement_r2, tbaseline_info):
    """ Crop some licsalert products in time.  
    Inputs:
        processing_date | string YYYYMMDD | Date to crop to. Must be an acquisition date.  
        acq_dates | list of strings YYYYMMDD | Dates of acqusitions.  
    Returns:
        cropped in time deep copies.  
    History:
        2023_12_08 | MEG | Written
    """
    from copy import deepcopy
    from licsalert.licsalert import licsalert_date_obj
    
    # make deep copies that will become the cropped outputs.  
    sources_tcs_crop = deepcopy(sources_tcs)
    residual_tcs_crop = deepcopy(residual_tcs)
    displacement_r2_crop = deepcopy(displacement_r2)
    tbaseline_info_crop = deepcopy(tbaseline_info)
    reconstructions_crop = deepcopy(reconstructions)
    residuals_crop = deepcopy(residuals)
    
    # crop the time course information in time, +1 to make inclusive
    processing_date = licsalert_date_obj(processing_date, acq_dates)
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
            
    # crop the time series information in time
    for key in ['incremental', 'incremental_downsampled', 'incremental_mc_space', 'means_space', 
                'incremental_mc_time', 'means_time', 'mixtures_mc', 'means']:
        displacement_r2_crop[key] = displacement_r2_crop[key][:processing_date.acq_n, ]
        
    # croppping in time, note +1 to make upper index inclusive.  
    for key in ['acq_dates', 'baselines_cumulative']:
        tbaseline_info_crop[key] = tbaseline_info[key][:(processing_date.acq_n + 1)]
    # but no +1 here as tehre is no baselines of 0 for ifg acq0_acq0
    for key in ['ifg_dates', 'baselines']:
        tbaseline_info_crop[key] = tbaseline_info[key][:(processing_date.acq_n)]
        
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
        2024_12_04 | MEG | Written
    """
    from glob import glob
    from pathlib import Path
    
    
    licsalert_items = sorted(glob(str(licsalert_dir / '*')))
    
    # remove any items that are not a licsalert data directory.  
    delete_args = []
    for item_n, licsalert_item in enumerate(licsalert_items):
        item_name = Path(licsalert_item).parts[-1]
        if item_name in ["ICASAR_results", "LiCSAlert_history.txt", "aux_data_figs"]:
            delete_args.append(item_n)
    
    for delete_arg in delete_args[::-1]:
        del licsalert_items[delete_arg]
        
    if len(licsalert_items) == 0:
        raise Exception(f"There are no LiCSAlert date directories (of the form YYYYMMDD), so there "
                        f"are no LiCSAlert results to open.  Exiting.")
    
    final_date_dir = Path(sorted(licsalert_items)[-1])
    
    return final_date_dir
    
    
#%%

    
def open_tcs(final_date_dir):
    """ Open the time course data.  
    Inputs:
        
    Returns:
        sources_tcs | list of dicts | One item in list for each source, each item contains ['cumulative_tc', 'gradient', 'lines', 'sigma', 'distances', 't_recalculate'])
    History:
        2023_10_25 | MEG | Written
        2024_01_04 | MEG | move functionality that finds youngest direcotyr to new function.  
    """
    
    import pickle

    try:
        with open(final_date_dir / 'time_course_info.pkl', 'rb') as f:
            sources_tcs = pickle.load(f)
            residual_tcs = pickle.load(f)
        f.close()
    except:
        raise Exception(f"Unable to open the time course information for this date.  Perhaps it's "
                        f"a baseline date so doesn't have enough information to create the licsalert figure? Exiting  ")

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



def LiCSBAS_to_LiCSAlert(LiCSBAS_out_folder, filtered = False, figures = False, n_cols=5, crop_pixels = None, return_r3 = False, 
                      ref_area = True, mask_type = 'dem'):
    """ A function to prepare the outputs of LiCSBAS for use with LiCSALERT. Note that this includes the step of referencing the time series to the reference area selected by LiCSBAS.  
    LiCSBAS uses nans for masked areas - here these are converted to masked arrays.   Can also create three figures: 1) The Full LiCSBAS ifg, and the area
    that it has been cropped to 2) The cumulative displacement 3) The incremental displacement.  

    Inputs:
        h5_file | string | path to h5 file.  e.g. cum_filt.h5
        figures | boolean | if True, make figures
        n_cols  | int | number of columns for figures.  May want to lower if plotting a long time series
        crop_pixels | tuple | coords to crop images to.  x then y, 00 is top left.  e.g. (10, 500, 600, 900).  
                                x_start, x_stop, y_start, y_stop, No checking that inputted values make sense.  
                                Note, generally better to have cropped (cliped in LiCSBAS language) to the correct area in LiCSBAS_for_LiCSAlert
        return_r3 | boolean | if True, the rank 3 data is also returns (n_ifgs x height x width).  Not used by ICASAR, so default is False
        ref_area | boolean | If True, the reference area (in pixels, x then y) used by LiCSBAS is returned to the user.  
                            ##### Regardless of how this is set, the reference area is always extractd to reference the time series.  #######
        mask_type | string | 'dem' or 'licsbas'  If dem, only the pixels masked in the DEM are masked (i.e. pretty much only water.  ).  Note that any pixels that are incoherent in the time series are also masked (ie if incoherent in one ifg, will be masked for all ifgs.  )
                                                If licsbas, the pixels that licsbas thinks should be masked are masked (ie water + incoherent).  Note that the dem masked is added to the licsbas mask to ensure things are consistent, but there should be no change here (every pixel masked in the DEM is also masked in the licsbas mask)

    Outputs:
        displacment_r3 | dict | Keys: cumulative, incremental.  Stored as masked arrays.  Mask should be consistent through time/interferograms
                                Also lons and lats, which are the lons and lats of all pixels in the images (ie rank2, and not column or row vectors)    
                                Also Dem, mask, and  E N U (look vector components in east north up diretcion)
        displacment_r2 | dict | Keys: cumulative, incremental, mask.  Stored as row vectors in arrays.  
                                Also lons and lats, which are the lons and lats of all pixels in the images (ie rank2, and not column or row vectors)    
                                Also Dem, mask, and  E N U (look vector components in east north up diretcion)
        tbaseline_info | dict| imdates : acquisition dates as strings
                              daisy_chain : names of the daisy chain of ifgs, YYYYMMDD_YYYYMMDD
                              baselines : temporal baselines of incremental ifgs

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
    
    from licsalert.aux import col_to_ma
    from licsalert.icasar.aux import add_square_plot
    
    

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


    def ts_quick_plot(ifgs_r3, title):
        """
        A quick function to plot a rank 3 array of ifgs.
        Inputs:
            title | string | title
        """
        n_ifgs = ifgs_r3.shape[0]
        n_rows = int(np.ceil(n_ifgs / n_cols))
        fig1, axes = plt.subplots(n_rows,n_cols)
        fig1.suptitle(title)
        for n_ifg in range(n_ifgs):
            ax=np.ravel(axes)[n_ifg]                                                                            # get axes on it own
            matrixPlt = ax.imshow(ifgs_r3[n_ifg,],interpolation='none', aspect='equal')                         # plot the ifg
            ax.set_xticks([])
            ax.set_yticks([])
            fig1.colorbar(matrixPlt,ax=ax)                                                                       
            ax.set_title(f'Ifg: {n_ifg}')
        for axe in np.ravel(axes)[(n_ifgs):]:                                                                   # delete any unused axes
            axe.set_visible(False)

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
        value = subp.check_output(['grep', field,mlipar]).decode().split()[1].strip()
        return value
    
        
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

    for LiCSBAS_folder in LiCSBAS_folders['all']:                                                                   # 1: Loop though looking for the TS direcotry
        if bool(re.match(re.compile('TS_.'), LiCSBAS_folder)):                                                      # the timeseries output, which is named depending on mutlitlooking and clipping.  
            LiCSBAS_folders['TS_'] = LiCSBAS_folder
    
    for LiCSBAS_folder in LiCSBAS_folders['all']:                                                                   # 2a: Loop though looking for the ifgs directory, which depends on lots of things.  
        if re.match(re.compile('GEOCml.+clip'), LiCSBAS_folder):                                                    # multilooked and clipped
                LiCSBAS_folders['ifgs'] = LiCSBAS_folder

    if 'ifgs' not in LiCSBAS_folders.keys():                                                                        # 2b: If we haven't found it already
        for LiCSBAS_folder in LiCSBAS_folders['all']:                                                               # Loop though looking for the other way the ifgs directory can be called
            if re.match(re.compile('GEOCml.+'), LiCSBAS_folder):                                                    # or just multilooked 
                LiCSBAS_folders['ifgs'] = LiCSBAS_folder

    if 'ifgs' not in LiCSBAS_folders.keys():                                                                        # 2c if we haven't found it already
        for LiCSBAS_folder in LiCSBAS_folders['all']:                                                               # loop through
            if re.match(re.compile('GEOC'), LiCSBAS_folder):                                                        # neither multilooked or clipped
                LiCSBAS_folders['ifgs'] = LiCSBAS_folder
    
    if ('TS_' not in LiCSBAS_folders) or ('ifgs' not in LiCSBAS_folders):
        raise Exception(f"Unable to find the TS_* and ifgs  directories that contain the LiCSBAS results.  Perhaps the LiCSBAS directories have unusual names?  Exiting.  ")


    # 1: Open the h5 file with the incremental deformation in.  
    displacement_r3 = {}                                                                                        # here each image will 1 x width x height stacked along first axis
    displacement_r2 = {}                                                                                        # here each image will be a row vector 1 x pixels stacked along first axis
    tbaseline_info = {}

    if filtered:
        cumh5 = h5.File(LiCSBAS_out_folder / LiCSBAS_folders['TS_'] / 'cum_filt.h5' ,'r')                       # either open the filtered file from LiCSBAS
    else:
        cumh5 = h5.File(LiCSBAS_out_folder / LiCSBAS_folders['TS_'] / 'cum.h5' ,'r')                            # or the non filtered file from LiCSBAS
    tbaseline_info["acq_dates"] = cumh5['imdates'][()].astype(str).tolist()                                     # get the acquisition dates
    cumulative = cumh5['cum'][()]                                                                               # get cumulative displacements as a rank3 numpy array
    cumulative *= 0.001                                                                                         # LiCSBAS default is mm, convert to m


    # 2: Open the parameter file to get the number of pixels in width and height (though this should agree with above)   
    try:
        width = int(get_param_par(LiCSBAS_out_folder / LiCSBAS_folders['ifgs'] / 'slc.mli.par', 'range_samples'))
        length = int(get_param_par(LiCSBAS_out_folder / LiCSBAS_folders['ifgs'] / 'slc.mli.par', 'azimuth_lines'))
    except:
        print(f"Failed to open the 'slc.mli.par' file, so taking the width and length of the image from the h5 file and trying to continue.  ")
        (_, length, width) = cumulative.shape

    # 3: Reference the time series
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
        
    #######
    # print(f"DEBUG")
    # f, ax = plt.subplots()
    # ax.plot(cumh5['cum'][:,76, 122 ])
    # ax.set_title('After referencing')
    # pdb.set_trace()
    ##############
    
    
    
    #3: Open the mask and the DEM
    mask_licsbas = read_img(LiCSBAS_out_folder / LiCSBAS_folders['TS_'] / 'results' /  'mask', length, width)                   # this is 1 for coherenct pixels, 0 for non-coherent/water, and masked for water
    mask_licsbas = np.logical_and(mask_licsbas, np.invert(np.isnan(mask_licsbas)))                                          # add any nans to the mask (nans become 0)
    mask_licsbas = np.invert(mask_licsbas)                                                                              # invert so that land is 0 (not masked), and water and incoherent are 1 (masked)
    
    dem = read_img(LiCSBAS_out_folder / LiCSBAS_folders['ifgs'] / 'hgt', length, width)
    mask_dem = np.isnan(dem)
    
    mask_cum = np.isnan(cumulative)                                                                       # get where masked
    mask_cum = (np.sum(mask_cum, axis = 0) > 0)                                                             # sum all the pixels in time (dim 0), and if ever bigger than 0 then must be masked at some point so mask.  
    
    if mask_type == 'dem':
        mask = np.logical_or(mask_dem, mask_cum)                                                        # if dem, mask water (from DEM), and anything that's nan in cumulative (mask_cum)
    elif mask_type == 'mask_licsbas':
        mask = np.logical_or(mask_licsbas, np.logical_or(mask_dem, mask_cum))
        
    mask_r3 = np.repeat(mask[np.newaxis,], cumulative.shape[0], 0)
    dem_ma = ma.array(dem, mask = mask)
    displacement_r2['dem'] = dem_ma                                                                      # and added to the displacement dict in the same was as the lons and lats
    displacement_r3['dem'] = dem_ma                                                                      # 
    
    
    # 3: Mask the data  
    displacement_r3["cumulative"] = ma.array(cumulative, mask=mask_r3)                                   # rank 3 masked array of the cumulative displacement
    displacement_r3["incremental"] = np.diff(displacement_r3['cumulative'], axis = 0)                           # displacement between each acquisition - ie incremental
    if displacement_r3["incremental"].mask.shape == ():                                                         # in the case where no pixels are masked, the diff operation on the mask collapses it to nothing.  
        displacement_r3["incremental"].mask = mask_r3[1:]                                                # in which case, we can recreate the mask from the rank3 mask, but dropping one from the first dimension as incremental is always one smaller than cumulative.  
    n_im, length, width = displacement_r3["cumulative"].shape                                   

    # if figures:                                                 
    #     ts_quick_plot(displacement_r3["cumulative"], title = 'Cumulative displacements')
    #     ts_quick_plot(displacement_r3["incremental"], title = 'Incremental displacements')

    displacement_r2['cumulative'], displacement_r2['mask'] = rank3_ma_to_rank2(displacement_r3['cumulative'])      # convert from rank 3 to rank 2 and a mask
    displacement_r2['incremental'], _ = rank3_ma_to_rank2(displacement_r3['incremental'])                          # also convert incremental, no need to also get mask as should be same as above

    # 4: work with the acquisiton dates to produces names of daisy chain ifgs, and baselines
    tbaseline_info["ifg_dates"] = daisy_chain_from_acquisitions(tbaseline_info["acq_dates"])
    tbaseline_info["baselines"] = baseline_from_names(tbaseline_info["ifg_dates"])
    # calculate the cumulative baselines from the baselines, ensure 0 at start.  
    tbaseline_info["baselines_cumulative"] = np.concatenate((np.zeros((1)), np.cumsum(tbaseline_info["baselines"])), axis = 0)                                                            
    
    # 5: get the lons and lats of each pixel in the ifgs
    geocode_info = create_lon_lat_meshgrids(cumh5['corner_lon'][()], cumh5['corner_lat'][()], 
                                            cumh5['post_lon'][()], cumh5['post_lat'][()], displacement_r3['incremental'][0,:,:])             # create meshgrids of the lons and lats for each pixel
    displacement_r2['lons'] = geocode_info['lons_mg']                                                                                        # add to the displacement dict
    displacement_r2['lats'] = geocode_info['lats_mg']
    displacement_r3['lons'] = geocode_info['lons_mg']                                                                                        # add to the displacement dict (rank 3 one)
    displacement_r3['lats'] = geocode_info['lats_mg']

    # 6: Get the E N U files (these are the components of the ground to satellite look vector in east north up directions.  )   
    try:
        for component in ['E', 'N', 'U']:
            look_vector_component = read_img(LiCSBAS_out_folder / LiCSBAS_folders['ifgs'] / f"{component}.geo", length, width)
            displacement_r2[component] = look_vector_component
            displacement_r3[component] = look_vector_component
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
        
        # loop through and crop, determing if r2 or r3 array.  
        for product in displacement_r3:
            if len(displacement_r3[product].shape) == 2:                                                                                  
                resized_r2 = displacement_r3[product][crop_pixels[2]:crop_pixels[3], crop_pixels[0]:crop_pixels[1]]               
                displacement_r2[product] = resized_r2
                displacement_r3[product] = resized_r2
            elif len(displacement_r3[product].shape) == 3:                                                                                
                resized_r3 = displacement_r3[product][:, crop_pixels[2]:crop_pixels[3], crop_pixels[0]:crop_pixels[1]]            
                displacement_r3[product] = resized_r3
                displacement_r2[product], displacement_r2['mask'] = rank3_ma_to_rank2(resized_r3)      
            else:
                pass
        
    

        if figures:
            add_square_plot(crop_pixels[0], crop_pixels[1], crop_pixels[2], crop_pixels[3], ax)                 # draw a box showing the cropped region    
  
   

    if return_r3:
        if ref_area:
            return displacement_r3, displacement_r2, tbaseline_info, ref_xy
        else:
            return displacement_r3, displacement_r2, tbaseline_info
    else:
        if ref_area:
            return displacement_r2, tbaseline_info, ref_xy
        else:
            return displacement_r2, tbaseline_info


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
    # make cumulative baselines, and ensure that 0 at the start.  
    tbaseline_info["baselines_cumulative"] = np.concatenate((np.zeros((1,1)), np.cumsum(tbaseline_info["baselines"])), axis = 0)                                                            
    
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
