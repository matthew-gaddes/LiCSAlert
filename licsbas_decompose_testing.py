#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 10:57:24 2024

@author: matthew
"""

print(f"Started")

import sys
import pickle
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pdb
from copy import deepcopy

import licsalert
from licsalert.aux import col_to_ma

# from licsalert.licsalert import reconstruct_ts_from_dir

#%% Visualise ascending and descnding frame look vectors.  


# Example usage:

# Define angles for the two vectors
# Each tuple is (theta_deg, phi_deg_from_north)
# S1 inc: 29.1 - 46, mean = 37.55
# s1 heading: 350 (asc), 190 (desc)
# S1 LOS is +90' to heading, so 80 (asc), 280
# vector1_angles = (37.55, 260)  # inc, and from north, ascending
# vector2_angles = (37.55, 100) # descending
# plot_vectors_and_plane(vector1_angles, vector2_angles)


#%%

def plot_timeseries(interpolated_ts, all_dates, frame_names, title = ''):
    """
    Plots each column of a 2D array against a list of dates.

    Parameters:
    - data (numpy.ndarray): 2D array with shape (n_times, n_columns).
    - date_strings (list of str): List of date strings corresponding to each row in the array.

    Returns:
    - None
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

from licsalert.decomposition import licsbas_frame

from licsalert.decomposition import resample_in_space
from licsalert.decomposition import consistent_mask_across_frames
from licsalert.decomposition import resample_to_daily
from licsalert.decomposition import interpolate_missing_times, sample_new_times
from licsalert.decomposition import apply_consistent_zero

from licsalert.decomposition import extract_los_info, decompose_timeseries
            

#%% Set the data directory

#downsample = 0.5
downsample = 0.2
# originally 100m pixels, so attempt to preserve.  
pixel_spacing_m = 100 / downsample


licsbas_dir_d = Path(
    "/home/matthew/university_work/data/00_LiCSBAS_time_series/"
    "022D_04826_121209_rationalized")                                         

licsbas_dir_a = Path(
    "/home/matthew/university_work/data/00_LiCSBAS_time_series/"
    "044A_04913_071213_campi_flegrei_rationalized")              

licsbas_frames = [
    licsbas_frame(licsbas_dir_a),
    licsbas_frame(licsbas_dir_d),
    ]




#%% Open the data and downsample 

for frame in licsbas_frames:
    frame.import_data()
    frame.downsample_data(downsample)


#%% Resample all frames to a common grid.  

licsbas_frames = resample_in_space(licsbas_frames, pixel_spacing_m)  


#%% Ensure the coherence mask is consistent across frames 

licsbas_frames = consistent_mask_across_frames(licsbas_frames)

#%% resample to the data to all days, with nans for days with no data

# cumulative_daily is (n_days x n_frames x n_pixels), but most of the days
# are nans.  

cumulative_daily, all_dates, max_pixels = resample_to_daily(licsbas_frames)

plot_timeseries(cumulative_daily[:, [0,1], max_pixels], all_dates,
                [f.frame_name for f in licsbas_frames], title = 'Step 01: '
                'Resampled to daily ')



#%% Interpolate to fill days with no data (the nans).  

print("Creating a copy of the data to interpolate in time.  This can "
      "be slow.  ")


cumulative_daily_interp = deepcopy(cumulative_daily)

# interpolate each frame in time to fill the nans.  
for frame_n in range(len(licsbas_frames)):
    cumulative_daily_interp[:, frame_n, :] = interpolate_missing_times(
                                            cumulative_daily[:, frame_n, :])

# still (n_times x n_frames x n_pixels), but nans are filled.  
plot_timeseries(cumulative_daily_interp[:, [0,1], max_pixels], all_dates,
                [f.frame_name for f in licsbas_frames], title = 'Step 02: '
                'Interpolate to fill missing days')
  

# debug plot
# f, ax = plt.subplots(1, len(licsbas_frames))
# f.suptitle('Daily and interpolated')
# for frame_n in range(len(licsbas_frames)):
#     ax[frame_n].matshow(cumulative_daily_interp[:, frame_n,:])

#%% Sample to chosen days

# 1st frame are the ascending dates
cumulative_t_resampled, new_dates = sample_new_times(
    cumulative_daily_interp, all_dates, 
    licsbas_frames[0].tbaseline_info['acq_dates'], licsbas_frames
    )

plot_timeseries(cumulative_t_resampled[:, [0,1], max_pixels], new_dates,
                [f.frame_name for f in licsbas_frames], title = 'Step 03: '
                'Sample to consistent days')

# debug plot
# f, ax = plt.subplots(1, len(licsbas_frames))
# for frame_n in range(len(licsbas_frames)):
#     ax[frame_n].matshow(cumulative_t_resampled[:, frame_n,:])

   

#%% ensure 0 at start of time series

cumulative_t_resampled = apply_consistent_zero(cumulative_t_resampled)
    
plot_timeseries(cumulative_t_resampled[:, [0,1], max_pixels], new_dates,
                [f.frame_name for f in licsbas_frames], title = 'Step 04:'
                'Ensure 0 on first epoch date')
    
#%% Decompose

comp_e, comp_u, comp_n = extract_los_info(licsbas_frames)
    
m_un, m_e = decompose_timeseries(cumulative_t_resampled, comp_e, comp_u, 
                                 comp_n)





