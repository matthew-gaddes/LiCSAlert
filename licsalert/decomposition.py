#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:52:44 2024

@author: matthew
"""
import pdb

#%%


def plot_decomposition_timeseries(interpolated_ts, all_dates, frame_names, 
                                  title = ''):
    """
    Plots columns of a 2D array against a list of dates, so can be used to 
    visualise the displacement of a pixel through time in multiple frames.  
    
    Inputs:
        interpolated_ts | r2 array | n_acqs x n_frames
        all_dates | list | date for each row in interpolated_ts
        frame_names | list | name for each frame (column of interpolated_ts)
        title | string | Used to name the figure and the window.

    Returns:
        Figure
        
    2024_08_29 | MEG | Written
    """
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime
    
    # Convert date strings to datetime objects
    dates = [datetime.strptime(date, "%Y%m%d") for date in all_dates]
    
    # Check dimensions
    if interpolated_ts.shape[0] != len(all_dates):
        raise ValueError("Number of date strings must match the number of rows in the data array")

    # Create the plot
    f, ax = plt.subplots(figsize=(10, 6))
    f.suptitle(title)
    f.canvas.manager.set_window_title(title)
    
    # Plot each column
    for i in range(interpolated_ts.shape[1]):
        ax.scatter(dates, interpolated_ts[:, i], 
                   label=f'{frame_names[i]}', marker = '.')
        

    
    # Format the x-axis to show dates
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval = 6))
    
    # Rotate and align the date labels
    f.autofmt_xdate()
    
    # Add labels and title
    ax.set_xlabel('Date')
    ax.set_ylabel('LOS disp. (m)')
    #f.suptitle('Line Graph of Each Column vs Dates')
    f.legend()
    
    # Show the plot
    plt.show()



#%%
class licsbas_frame():
    """A class used to handle multiple LiCSBAS frames when decomposing them 
    into UN  and E
    
    """
    
    def __init__(self, data_dir):
        """ path to LiCSBAS data
        """
        self.data_dir = data_dir

        # frame name is assumed to be in the standard form and the first
        # part of the directory.  e.g. 044A_04913_071213
        self.frame_name = data_dir.parts[-1][:17]

    def import_data(self):
        """ use the path to import the data.  Requires LiCSAlert to be on
        the path.  
        """
        from licsalert.data_importing import LiCSBAS_to_LiCSAlert    
        from licsalert.data_importing import LiCSBAS_json_to_LiCSAlert
        
        def has_single_json_file(directory):
            """ """
            import os
            # Get a list of all .json files in the directory
            json_files = [f for f in os.listdir(directory) if f.endswith('.json')]
            
            # Check if there is exactly one .json file
            json_status = (len(json_files) == 1)
            
            if json_status:
                json_file = json_files[0]
            else:
                json_file = None
            return json_status, json_file
        
        
        # determine the type of intput data
        json_status, json_file = has_single_json_file(self.data_dir)
        
        # import data if it's in a .json file
        if json_status:        
            print(f"A single .json file was found in the directory, so "
                  "LiCSAlert is assuming that the data is in a COMET "
                  "Volcano Portal .json file")
            
            op = LiCSBAS_json_to_LiCSAlert(
                    self.data_dir / json_file, crop_side_length = 4e3,  
                    mask_type = 'nan_once'
                    )
            
            (self.displacement_r2, _, self.tbaseline_info, _, _) = op; del op
            
        # or if its just a LiCSBAS directory.  
        else:
            op = LiCSBAS_to_LiCSAlert(
                self.data_dir, filtered = False,  figures = True,  n_cols=5,
                crop_pixels = None, mask_type = 'dem',  date_start = None,
                date_end = None
                )
    
            (self.displacement_r2, self.tbaseline_info) = op; del op
        
        
    def downsample_data(self, ds_factor):
        """ Downsample the data in space (useful for fast exploratory runs).  
        """
        from licsalert.licsalert import LiCSAlert_preprocessing
        
        displacement_r2_ds = LiCSAlert_preprocessing(
            self.displacement_r2, self.tbaseline_info, 'sica', 
            ds_factor, 1.0
            ) 
        
        to_remove = ['incremental_downsampled', 'mask_downsampled', 
                     'incremental_mc_space', 'means_space', 
                     'incremental_mc_time', 'means_time', 'mixtures_mc', 
                     'means']
        
        for item in to_remove:
            del displacement_r2_ds[item]
        
        self.displacement_r2 = displacement_r2_ds
        # number of pixels in this frame.  
        self.n_pixels = displacement_r2_ds['cumulative'].shape[1]
        
        
#%%

def resample_in_space(licsbas_frames, pixel_spacing_m = 100):
    """ Sample multiple LiCSBAS outputs to a common grid.  
    
    Inputs:
        licsbas_frames | list of licsbas_frame objects | 
        pixel_spacing_m  | int | New grid to sample to.  Larger values = 
                                 coarser.
                                 
     Returns:
         licsbas_frames_resampled | list of licsbas_frame objects | as above,
                                     now sampled to a common grid.  
     
    History:
        2024_08_28 | MEG Written
            
    """
    import numpy as np
    import numpy.ma as ma 
    from copy import deepcopy
    from licsalert.aux import r2_to_r3, r3_to_r2
    
    def find_lon_lat_limits(licsbas_frames):
        """
        """
        # Initialize min and max with the first coordinate values
        min_lon = np.min(licsbas_frames[0].displacement_r2['lons'])
        max_lon = np.max(licsbas_frames[0].displacement_r2['lons'])
        min_lat = np.min(licsbas_frames[0].displacement_r2['lats'])
        max_lat = np.max(licsbas_frames[0].displacement_r2['lats'])
        
        for lbf in licsbas_frames:
            min_lon_current = np.min(lbf.displacement_r2['lons'])
            max_lon_current = np.max(lbf.displacement_r2['lons'])
            
            min_lat_current = np.min(lbf.displacement_r2['lats'])
            max_lat_current = np.max(lbf.displacement_r2['lats'])
    
            if min_lon_current < min_lon:
                min_lon = min_lon_current
            if max_lon_current > max_lon:
                max_lon = max_lon_current
            if min_lat_current < min_lat:
                min_lat = min_lat_current
            if max_lat_current > max_lat:
                max_lat = max_lat_current
                
        return min_lon, max_lon, min_lat, max_lat
    
    def find_nearest_pixels(lons1, lats1, lons2, lats2):
        """
        Finds the indices of the nearest pixels in the second meshgrid (lons2, lats2) 
        for each pixel in the first meshgrid (lons1, lats1).
    
        Parameters:
        lons1, lats1 : 2D arrays
            Meshgrids of longitudes and latitudes for the first set of pixels.
        lons2, lats2 : 2D arrays
            Meshgrids of longitudes and latitudes for the second set of pixels.
    
        Returns:
        nearest_indices : 3D array
            Array of shape (n, m, 2) where each element contains the (x, y) index 
            of the nearest pixel in the second meshgrid corresponding to the pixel 
            in the first meshgrid.
        """
        
    
        # Get the shape of the first meshgrid
        shape = lons1.shape
    
        # Flatten the meshgrids
        lons1_flat = lons1.flatten()
        lats1_flat = lats1.flatten()
        lons2_flat = lons2.flatten()
        lats2_flat = lats2.flatten()
    
        # Initialize array to store the (x, y) indices of the nearest pixels
        nearest_indices = np.zeros((lons1_flat.size, 2), dtype=int)
    
        # Iterate over each pixel in the first meshgrid
        for i, (lon1, lat1) in enumerate(zip(lons1_flat, lats1_flat)):
            # Calculate the Euclidean distance to all pixels in the second meshgrid
            distances = np.sqrt((lons2_flat - lon1)**2 + (lats2_flat - lat1)**2)
            
            # Find the index of the minimum distance
            min_idx = np.argmin(distances)
    
            # Convert the flat index back to 2D (x, y) index
            nearest_indices[i] = np.unravel_index(min_idx, lons2.shape)
    
        # Reshape the nearest_indices array back to the shape of the original meshgrids with the last dimension being (x, y)
        nearest_indices = nearest_indices.reshape(shape + (2,))
    
        return nearest_indices
    
    # make new grid that all data is contained within        
    min_lon, max_lon, min_lat, max_lat = find_lon_lat_limits(licsbas_frames)
    
    # convert pixel spacing
    pixel_spacing_deg = pixel_spacing_m / 111320
    #1d vector for lons and lats
    lons = np.arange(min_lon, max_lon, pixel_spacing_deg)
    lats = np.arange(min_lat, max_lat, pixel_spacing_deg)
    # make 2d, check lats are correct
    lons_mg, lats_mg = np.meshgrid(lons, lats)
    if lats_mg[0,0] < lats_mg[-1,0]:
        lats_mg = np.flipud(lats_mg)
        
    print(f"The grid for resampling is size: {lons_mg.shape} with pixels " 
          f"spaced every {pixel_spacing_m} m, or every {pixel_spacing_deg} "
          "deg.")
    
    # # debug test
    # f, ax = plt.subplots(1, 2)
    # ax[0].matshow(lons_mg)
    # ax[0].set_title('lons_mg')
    # ax[1].matshow(lats_mg)
    # ax[1].set_title('lats_mg')

        
    # do resampling
    licsbas_frames_resampled = deepcopy(licsbas_frames)
    
    for frame_n, lbf in enumerate(licsbas_frames_resampled):
        
        # for the new meshgrid, find the indices of the nearest original pixel
        nearest_indices = find_nearest_pixels(lons_mg, lats_mg,
            lbf.displacement_r2['lons'], lbf.displacement_r2['lats'],
            )
        
        #convert to rank 3 ready to resample in space
        lbf.displacement_r2['cumulative'] = r2_to_r3(
            lbf.displacement_r2['cumulative'], lbf.displacement_r2['mask'])
        
        # remove these as not needed 
        del lbf.displacement_r2['incremental'], lbf.displacement_r2['mask']

        # iterate through and sample old pixels to new grid
        print(f"\nFor frame {lbf.frame_name}, resampling: ")
        for item in ['dem', 'E', 'N', 'U']:
            print(f"    {item}", end = '')
            # open the current 2D data
            data = licsbas_frames[frame_n].displacement_r2[item]
            # initiliase empty for the new data at new sampling.  
            new_data = np.zeros(lons_mg.shape)
            # iterate through pixels of new gid.  
            for nx in range(lons_mg.shape[1]):
                for ny in range(lons_mg.shape[0]):
                    
                    # find the old value at the location
                    pixel_value = data[nearest_indices[ny, nx, 0], 
                                       nearest_indices[ny, nx, 1]]
                    # and copy to the new one
                    new_data[ny, nx] = pixel_value
            print(f"  (Originally: {data.shape} Now: {new_data.shape})")
            # write to the licsbas frame
            lbf.displacement_r2[item] = new_data
        
        # time series of data and mask are handled differently.  
        print(f"    cumulative and mask", end = '')
        n_acqs = lbf.displacement_r2['cumulative'].shape[0]
        data = lbf.displacement_r2['cumulative']
        new_data = ma.zeros((n_acqs, lons_mg.shape[0], lons_mg.shape[1]))
        for nx in range(lons_mg.shape[1]):
            for ny in range(lons_mg.shape[0]):
                
                # find the old value at the location
                pixel_value = (data[:, nearest_indices[ny, nx, 0], 
                                    nearest_indices[ny, nx, 1]])
                # and copy to the new one
                new_data[:, ny, nx] = pixel_value
                    
        # convert r3 back to r2
        r2_data = r3_to_r2(new_data)
        lbf.displacement_r2['cumulative'] = r2_data['ifgs']
        lbf.displacement_r2['mask'] = r2_data['mask']
        print(f"  (Originally: {data.shape[1:]} Now: {new_data.shape[1:]})")

        # update the long lat info to new grid  
        lbf.displacement_r2['lons'] = lons_mg
        lbf.displacement_r2['lats'] = lats_mg            

        # # debug plot
        # f, ax = plt.subplots(1,2)
        # ax[0].matshow(licsbas_frames[0].displacement_r2['dem'])
        # ax[0].set_title('Original frame')
        # ax[1].matshow(licsbas_frames_resampled[0].displacement_r2['dem'])
        # ax[1].set_title('Sampled to new grid')

        # update the number of pixels        
        lbf.n_pixels = r2_data['ifgs'].shape[1]
        
        print("\n")
        
        
    # pdb.set_trace()
    # for lbf in licsbas_frames_resampled:
    #     print(lbf.n_pixels)

    #     f, ax = plt.subplots()
    #     ax.matshow(lbf.displacement_r2['dem'])
    return licsbas_frames_resampled
   
#%%

def consistent_mask_across_frames(licsbas_frames):
    """  Different frames have different areas of incoherence.  
    Despite being on the same grid, they may have a different set of 
    valid pixels.  
    
    This is designed for location with only two tracks.  With more frames,
    a pixel doesn't have to be coherent in all frames to still have multiple
    look angles that could be decomposed.  
    
    Masked arrays this acts on: cumulative, DEM
    Normal arrays it doesn't: lons, alts, E, N, U
    
    Inputs:
        licsbas_frames | list of licsbas_frame objects | 
        
    Returns:
        licsbas_frames_masked | list of licsbas_frame objects | as above, but
                                all with the same mask
    History:
        2024_08_28 | MEG Written
    """
    
    from copy import deepcopy
    import numpy as np
    import numpy.ma as ma
    from licsalert.aux import r2_to_r3, r3_to_r2
    
    licsbas_frames_masked = deepcopy(licsbas_frames)
    
    # mask is 1 where masked (water etc).  
    masks = [lbf.displacement_r2['mask'] for lbf in licsbas_frames]

    # debug 
    # print(licsbas_frames[0].displacement_r2['lons'].shape[0] * 
    #       licsbas_frames[0].displacement_r2['lons'].shape[1])
    # for mask in masks:
    #     print(len(np.argwhere(mask == False)))
    #     print(len(np.argwhere(mask == True)))

    mask = np.logical_or.reduce(masks)
    
    # iterate through the frames applying the new mask
    for lbf in licsbas_frames_masked:
        
        print(f"LiCSAR frame {lbf.frame_name} originally had {lbf.n_pixels} "
              f" pixels", end = '' )

        # dem uses nans for no data.  Make anywhere mask is True also nan        
        lbf.displacement_r2['dem'][mask] = np.nan
        
        cumulative_r3 = r2_to_r3(lbf.displacement_r2['cumulative'],
                                 lbf.displacement_r2['mask'])

        # make the mask rank 3
        mask_r3 = test = ma.repeat(mask[np.newaxis,], cumulative_r3.shape[0], 
                                   axis = 0)
        # apply the new mask
        cumulative_r3.mask = mask_r3
        
        # convert back to r2
        r2_data = r3_to_r2(cumulative_r3)
        lbf.displacement_r2['cumulative'] = r2_data['ifgs']
        lbf.displacement_r2['mask'] = r2_data['mask']
        
        # update the number of pixels.  Should always be constant (as shared
        # mask)
        lbf.n_pixels = r2_data['ifgs'].shape[1]
        print(f", but now has {lbf.n_pixels}")
        
    return licsbas_frames_masked


#%%



def resample_to_daily(licsbas_frames):
    """
    Each LiCSBAS time series has measurements to each acquisition.  
    An array is created with a row for each day, and the measurments
    for each frame are inserted into columns of this.  The majority of entries
    are nans as acqusitions are every 6/12 days etc.  
    
    Inputs:
        licsbas_frames | list of licsbas_frame objects | Measurments
                            only every e.. 6/12/24 days.  
    
    Returns:
        disp_daily
        all_dates
        max_pixels
        
    History:
        2024_08_29 | MEG | Written.  
        
    """
    import numpy as np
    from datetime import datetime, timedelta
    from copy import deepcopy
    
    def find_date_extremes(licsbas_frames):
        """
        """
        from datetime import datetime
        # Initialize empty list to hold all dates
        all_dates = []
        
        # Iterate through each dictionary and extract 'acq_dates'
        for frame in licsbas_frames:
            all_dates.extend(frame.tbaseline_info['acq_dates'])
        
        # Convert all dates to datetime objects
        all_dates = [datetime.strptime(date, "%Y%m%d") for date in all_dates]
        
        # Find the first (earliest) and last (latest) dates
        first_date = min(all_dates)
        last_date = max(all_dates)
        
        # Convert back to string format if needed
        first_date_str = first_date.strftime("%Y%m%d")
        last_date_str = last_date.strftime("%Y%m%d")
        
        return first_date_str, last_date_str
    
    
    # get the number of frames
    n_frames = len(licsbas_frames)
    print(f"{n_frames} frames have been detected.  ")

    # get the number of pixels for each frame.      
    n_pixels = [frame.displacement_r2['cumulative'].shape[1] 
                for frame in licsbas_frames]
    
    # get the total number of days 
    first_date_str, last_date_str = find_date_extremes(licsbas_frames)
    start_date = datetime.strptime(first_date_str, "%Y%m%d")
    end_date = datetime.strptime(last_date_str, "%Y%m%d")
    all_dates = [(start_date + timedelta(days=x)).strftime('%Y%m%d')
                 for x in range((end_date - start_date).days + 1)]
    n_days = len(all_dates)
    
    # get the maximum displacement pixel for each frame (column #).  
    max_pixels = [  
         np.unravel_index(np.argmax(frame.displacement_r2['cumulative']), 
         frame.displacement_r2['cumulative'].shape)[1] 
         for frame in licsbas_frames
         ]
    
    # get the number of pixels for each frame.  
    n_pixels = [lbf.n_pixels for lbf in licsbas_frames]

    # initialise array to store cumulative disps
    disp_daily = np.nan * np.ones((n_days, n_frames, np.max(n_pixels)))
    

    # copy the data to the daily displacements    
    for frame_n, lbf in enumerate(licsbas_frames):
        disp_cumulative = lbf.displacement_r2['cumulative']
        n_pix = disp_cumulative.shape[1]
        for date_n, date in enumerate(lbf.tbaseline_info['acq_dates']):
            day_idx = all_dates.index(date)
            disp_daily[day_idx, frame_n, :n_pix] = disp_cumulative[date_n, :]
            
    return disp_daily, all_dates, max_pixels


#%%

def interpolate_missing_times(data):
    """
    Interpolates missing rows (times, rows of NaNs) in a 2D array.

    Inputs:
        data | r2 array | rows are times, columns are variables.  Rows of nans
                            where variables are missing at certain times.  

    Returns: 
        data | r2 array | as above, but with nans filled.  

    History:
        2024_08_29 | MEG | Written.  
            
    """
    import numpy as np
    
    # # Copy the data to avoid modifying the original array
    # data = np.array(data, dtype=float)
    
    # Get the indices of non-NaN rows
    not_nan_indices = np.where(~np.isnan(data[:, 0]))[0]
    
    # Iterate over columns
    for col in range(data.shape[1]):
        print(f"    Interpolating pixel {col} of {data.shape[1]} pixels.  ")
        if np.all(np.isnan(data[:, col])):
            # this column has no data and can't be interpolated.  
            #print("        Pixel contains nans at all times.  Skipping.  ")
            pass
        else:
            # Iterate over the missing rows
            for i in range(len(not_nan_indices) - 1):
                start_index = not_nan_indices[i]
                end_index = not_nan_indices[i + 1]
                
                if end_index - start_index > 1:
                    # Interpolate between the start and end index
                    x = np.arange(start_index, end_index)
                    y = data[start_index:end_index + 1, col]
                    
                    # Only interpolate where y values are NaN
                    if np.isnan(y).any():
                        valid_x = np.arange(start_index, end_index + 1)[~np.isnan(y)]
                        valid_y = y[~np.isnan(y)]
                        
                        # Linear interpolation
                        interpolated_values = np.interp(x, valid_x, valid_y)
                        
                        # Assign the interpolated values to the array
                        data[x, col] = interpolated_values

    return data

#%%


def sample_new_times(cumulative_daily_interp, all_dates, desired_dates,
                     licsbas_frames):
    """
    Given measurements of displacement that are daily (usually interpolating
    between satellite acquisitions), sample them on new days.  Usually used
    to sample interpolated ascending to descending (or vice versa).  
    
    Inputs:
        cumulative_daily_interp | r3 array | (n_times, n_frames, n_pixels)
        all_dates | list | date for each row in cumulative_daily_interp
        desired | list | dates that data will be sampled to.  
        licsbas_frames | list of licsbas_frame | used to get the number of
                                                pixels in each frame.  
        
    Returns:
        cumulative_t_resampled | r3 array | (n_times, n_frames, n_pixels)
                                            Now with far fewer rows as sampled
                                            to the new_dates dates.  
        new_dates | list | new dates.  Note that this is not necessaily the 
                            same as the one input, as if a date has nans
                            it is removed.  
                                        
    History:
        2024_08_29 | MEG | Written
        
        
    """
    import numpy as np  
    from copy import deepcopy
    
    # make a copy to avoid changes
    new_dates = deepcopy(desired_dates)
    
    _, n_frames, n_pixels_max = cumulative_daily_interp.shape

    # get the number of pixels in each frame
    n_frames_pixels = [i.n_pixels for i in licsbas_frames]

    n_times_new = len(new_dates)

    # initilaise to store output
    cumulative_t_resampled = np.zeros((n_times_new, n_frames, n_pixels_max))
    
    nan_days = []

    # loop through all new dates and find cumulative displacement on that date
    for day_n, new_date in enumerate(new_dates):
        # find which daily solution number the current date is 
        day_index = all_dates.index(new_date)
        
        # if there are nans on that date in any frame:
        nans_present = False
        for frame_n in range(n_frames):
            # get the cumulative displacement to that date (1d)
            cum_on_date = cumulative_daily_interp[day_index, frame_n, 
                                                  :n_frames_pixels[frame_n]]
            # update the nans present flag if there are nans for this frame
            nans_present = nans_present or np.isnan(cum_on_date).any()
            
        # either record or pass, depending on nans status.  
        if nans_present:
            # record that as a date to remove late
            nan_days.append(new_date)
        else:
            cumulative_t_resampled[day_n, ] = cumulative_daily_interp[day_index, :, :]
            
            
    # remove any dates we don't have data for in all frames 
    # usually this is as we can't interpolate before the start or after the 
    # end of the time series.  
    
    if len(nan_days) > 0:
        
        for nan_day in sorted(nan_days)[::-1]:
            day_index = new_dates.index(nan_day)
            del new_dates[day_index]
            cumulative_t_resampled = np.delete(
                cumulative_t_resampled, day_index, axis=0
                )
        
    return cumulative_t_resampled, new_dates
    

#%% 

def apply_consistent_zero(cumulative_t_resampled):
    """    Different frames have all zeros on their first acquisition, which
    is now not necesarily the same date.  This function zeros the data
    to the first date in the array.  
    
    Inputs:
        cumulative_t_resampled | r3 array | (n_times, n_frames, n_pixels)
    
    Returns:
        cumulative_t_resampled_0 | r3 array | (n_times, n_frames, n_pixels),
                                              first row is now all zeros.  
                                              
    History:
        2024_08_29 | MEG | Written.  
    """
    
    import numpy as np 
    
    cumulative_t_resampled_0 = np.zeros(cumulative_t_resampled.shape)
    
    n_acqs = cumulative_t_resampled.shape[0] 
    offset = cumulative_t_resampled[0,]
       
    for acq_n in range(n_acqs):
        cumulative_t_resampled_0[acq_n] = (cumulative_t_resampled[acq_n] -
                                           offset)

    return cumulative_t_resampled_0


#%%


def extract_los_info(licsbas_frames):
    """ For a list of licsbas_frame objects, extract the components of the 
    ground to satellite look vectors for each pixel.  
    
    Inputs:
        licsbas_frames | list of licsbas_frame objects | 
        
    Returns:
        E | r2 array | n_pixels x n_frames 
        N | r2 array | n_pixels x n_frames 
        U | r2 array | n_pixels x n_frames 
        
    History:
        2024_08_29 | MEG | Written.  
    """
    
    import numpy as np
    import numpy.ma as ma
    
    # get the number of pixels
    n_pixels = licsbas_frames[0].n_pixels
    
    # get the number of frames
    n_frames = len(licsbas_frames)
    
    # initialise
    comps = {}
    
    # get the mask, same for all frame in my implementation
    mask = licsbas_frames[0].displacement_r2['mask']
    
    for comp in ['E', 'N', 'U']:
        
        # initialise
        comps[comp] = np.zeros((n_pixels, n_frames))
        
        for frame_n, lbf in enumerate(licsbas_frames):
            # convert the components into a masked array.  
            comp_masked = ma.array(lbf.displacement_r2[comp], mask = mask)
            
            # get values only whre the mask is valid.  
            comps[comp][:, frame_n] = ma.compressed(comp_masked)
            
    return comps['E'], comps['N'], comps['U']
            

#%%


def decompose_timeseries(cumulative_resampled, comp_e, comp_n, comp_u):
    """ Given a time series of images, with each image viewed from different
    look angles, decompose each image in UN (Up North) and E (East) components.  
    
    Inputs:
        cumulative_resampled | r3 array | n_times x n_frames x n_pixels
        comp_e | r2 array | n_times x n_frames, east component of look vector
                            for each pixel in each frame.  
        comp_n | r2 array | n_times x n_frames, nort component of look vector
                            for each pixel in each frame.  
        comp_u | r2 array | n_times x n_frames, up component of look vector
                            for each pixel in each frame.  
                            
    Returns:
        cumulative_un | r2 array | n_times x n_pixels, Up North displacement
                                    for each pixel at each time.  
                                    
        cumulative_e | r2 array |  n_times x n_pixels, East displacement
                                    for each pixel at each time.  
                                    
    History:
        2024_08_29 | MEG | Written.  

    """
    import numpy as np
    
    n_acqs, n_frames, n_pixels = cumulative_resampled.shape
    
    # initialise
    cumulative_un = np.zeros((n_acqs, n_pixels))
    cumulative_e = np.zeros((n_acqs, n_pixels))
        
    
    for n_acq in range(n_acqs):
        
        # d1 is n_pixels x n_frames
        d1 = cumulative_resampled[n_acq, ].T

        # decompose that time step.  
        m_un, m_e, var_un, var_e = decompose_image(d1, comp_e, comp_n, comp_u)
    
        # debug plot
        # f, axes = plt.subplots(1,4)
        # axes[0].matshow(col_to_ma(d1[:,0], licsbas_frames[0].displacement_r2['mask']))
        # axes[1].matshow(col_to_ma(d1[:,1], licsbas_frames[0].displacement_r2['mask']))
        # axes[2].matshow(col_to_ma(m_un, licsbas_frames[0].displacement_r2['mask']))
        # axes[3].matshow(col_to_ma(m_e, licsbas_frames[0].displacement_r2['mask']))
        # for ax, title in zip(axes, ['Frame 1', 'Frame 2', 'UN', 'E']):
        #     ax.set_title(title)
        
        cumulative_un[n_acq, ] = m_un
        cumulative_e[n_acq, ] = m_e
        
        print(f"Succesfully decomposed image {n_acq} of {n_acqs}")
        
    return cumulative_un, cumulative_e


#%%



def decompose_image(image, compE, compN, compU):
    """ Given a set of images that cove the same time span but are from 
    different look angles (e.g. interferograms spaning the same time period
    but from different tracks, or a map of velocities from different tracks),
    decompose these into the Up-North (UN) direction, and East (e) direction. 
    
    Code mofified from Andrew Watson's "Decompose InSAR Velocities"
    
    Inputs:
        image | r2 array | n_pixles x n_frames
        compE | r2 array | n_pixles x n_frames, east component of look vector
        compN | r2 array | n_pixles x n_frames, north component of look vector
        compU | r2 array | n_pixles x n_frames, up component of look vector
        
    Returns:
        m_un | r2 array |  1 x n_pixels, Up North component of displacement
        m_e | r2 array |  1 x n_pixels, East component of displacement
        var_un | r2 array |  1 x n_pixels, Variance for each pixel, currently
                                            all equal.  
        var_e | r2 array |  1 x n_pixels, Variance for each pixel, currently
                                            all equal.  
    History:
        2024_08_28 | MEG | Written.  
    """
    
    import numpy as np
    
    n_pixels, n_frames = image.shape
    
    #% calculate UN component vector, first by calculating the incidence angle
    #% and heading. Incidence angle is measured from the vertical, and azimuth
    #% is measured negatively counterclockwise from north. This is to match the
    #% definitions in Qi Ou's work.
    # Convert compU to 'inc'
    inc = 90 - np.degrees(np.arcsin(compU))
    # Calculate 'az'
    az = np.degrees(np.arccos(compE / np.sin(np.radians(inc)))) - 180
    # Calculate 'compUN' (n_pixels x n_frames)
    compUN = np.sqrt(1 - np.sin(np.radians(inc))**2 * np.cos(np.radians(az))**2)

   
    #initiliase to store Up North (UN) and East (E) displacement for each pixel
    m_un = np.nan * np.ones((1, n_pixels))
    m_e = np.nan * np.ones((1, n_pixels))
    var_un = np.nan * np.ones((1, n_pixels))
    var_e = np.nan * np.ones((1, n_pixels))
    
   
    # iterate over all the pixels (these are the ones that aren't masked)
    for n_pixel in range(n_pixels):
        
        # this should be the standard deviation of that pixel in each frame.  
        # replace with 1s for simplicity.  
        Qd = np.eye(n_frames)          
        
        # G is n_frames x 2
        G = np.concatenate((compUN[n_pixel : n_pixel+1, :].T,
                            compE[n_pixel : n_pixel+1 , :].T), axis = 1)
        
        # d is n_frames x 1 (observed displacement in each LOS)
        d = image[n_pixel : n_pixel + 1, :].T
        

        # Solve        
        # Calculate W as the inverse of Qd
        W = np.linalg.inv(Qd)
        # Calculate m
        m = np.linalg.inv(G.T @ W @ G) @ (G.T @ W @ d)
        # Calculate Qm
        Qm = np.linalg.inv(G.T @ W @ G)
        
        # record the two components of displacement and variance. 
        m_un[0, n_pixel] = m[0,0]
        m_e[0, n_pixel] = m[1,0]
        var_un[0, n_pixel] = Qm[0,0]
        var_e[0, n_pixel] = Qm[1,1]
        
    return m_un, m_e, var_un, var_e
        