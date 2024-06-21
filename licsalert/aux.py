#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 18:10:21 2020

@author: matthew
"""

import pdb

#%%

def get_licsalert_date_dirs(parent_dir):
    """Get the paths to only the LiCSAlert date directories
    (i.e. not the ones like ICASAR_results)
    
    Inputs:
        parent_dir | str (or maybe Path?) | licsalert frame directory
    Returns:
        licsalert_date_dirs | list of strings | names of date dirs.  Not full paths.  
    History:
        2024_06_14 | MEG | Written.  
    """
    import re
    from glob import glob
    import os
    from pathlib import Path
    
    # Regular expression for eight digits
    pattern = re.compile(r'^\d{8}$')

    licsalert_date_dirs = []

    # get all the dirs
    child_dirs = [d for d in glob(str(Path(parent_dir) / '*')) if os.path.isdir(d)]

    # loop through and check if they fit the YYYYMMDD form
    for child_dir in child_dirs:
        # Extract the directory name
        dir_name = os.path.basename(child_dir)
        if pattern.match(dir_name):
            licsalert_date_dirs.append(dir_name)
    licsalert_date_dirs = sorted(licsalert_date_dirs)
        
    return licsalert_date_dirs



#%%

def determine_abs_max_pixel(cumulative_r3, cumulative_r2):
    """ Given a time series as a rank 3 tensor, find the pixel
     with the largest absolution deformation (ie could be negative 
     or positive)
     Inputs:
         cumulative_r3 | r3 masked array | time x y x x
     Returns:
         y | int | y pixel of max abs deformation
         x | int | x pixel of max abs deformation
     History:
         2024_01_05 | MEG | Written
    """
    import numpy as np
    import numpy.ma as ma    
    
    def_max = ma.max(cumulative_r3)
    def_min = ma.min(cumulative_r3)
    
    if def_max > np.abs(def_min):
        t, y, x = np.unravel_index(ma.argmax(cumulative_r3), cumulative_r3.shape)
        col_r2 = np.argwhere(cumulative_r2 == def_max)[0,-1]
    else:
        t, y, x = np.unravel_index(ma.argmin(cumulative_r3), cumulative_r3.shape)
        col_r2 = np.argwhere(cumulative_r2 == def_min)[0,-1]
    return x, y, col_r2

#%%


def moving_average(ts, window = 11):
    """ A simple moving average function that reduces the window size at either edge of the 1d data.  
        Therefore the smoothed data is the same size as the original.  
        Can be tested with something like this:
            #smooth, valid =  moving_average(np.arange(20), window = 7)                                # testing fuctnion
    Inputs:
        ts | rank 1 data | the 1d time series to be smoothed.  
        window | int | odd number that sets the window size.  If 3, one data point either side of the data will be used for the smoothing.  
    
    Returns:
        ts_smooth | rank 1 data | smoothed version of the time series.  
        valid | rank 1 data | 0 if average for that point has edge effects (i.e. the ful window could not be used for averaging), 1 if full window was used.  
    History:
        2022_01_?? | MEG | Written
        2022_02_17 | MEG | add "valid" return.  
    """
    import numpy as np
    
    if window % 2 == 0:
        raise Exception(f"'window' must be odd, but is even.  Exiting.")
    half_window = int((window-1)/2)
    
    n_points = ts.shape[0]
    ts_smooth = np.zeros(n_points)
    valid = np.ones(n_points)
    
    for point_n in range(n_points):
        window_stop_adjusted = False                                                    # set (or reset)
#        pdb.set_trace()
        window_start = point_n - half_window
        if window_start < 0:
            window_start = 0
            window_stop = point_n + (point_n - window_start) + 1
            window_stop_adjusted = True
            valid[point_n] = 0
        
        if not window_stop_adjusted:
            window_stop = point_n + half_window +1
            if window_stop > n_points:
                window_stop = n_points                                          # if there are 20 points, this will also be 20 so can be indexced up to
                window_start = point_n - (window_stop - point_n) + 1
                valid[point_n] = 0

        ts_smooth[point_n] = np.mean(ts[window_start : window_stop])

        #print(f"for {ts[point_n]}: {ts[window_start : window_stop]}")               # debug/test
    return ts_smooth, valid




#%%



def r2_to_r3(ifgs_r2, mask):
    """ Given a rank2 of ifgs as row vectors, convert it to a rank3.   Copied from insar_tools.general to avoid making insar_tools a dependency.  
    Inputs:
        ifgs_r2 | rank 2 array | ifgs as row vectors 
        mask | rank 2 array | to convert a row vector ifg into a rank 2 masked array        
    returns:
        phUnw | rank 3 array | n_ifgs x height x width
    History:
        2020/06/10 | MEG  | Written
    """
    import numpy as np
    import numpy.ma as ma
        
    n_ifgs = ifgs_r2.shape[0]
    ny, nx = col_to_ma(ifgs_r2[0,], mask).shape                                   # determine the size of an ifg when it is converter from being a row vector
    
    ifgs_r3 = np.zeros((n_ifgs, ny, nx))                                                # initate to store new ifgs
    for ifg_n, ifg_row in enumerate(ifgs_r2):                                           # loop through all ifgs
        ifgs_r3[ifg_n,] = col_to_ma(ifg_row, mask)                                  
    
    mask_r3 = np.repeat(mask[np.newaxis,], n_ifgs, axis = 0)                            # expand the mask from r2 to r3
    ifgs_r3_ma = ma.array(ifgs_r3, mask = mask_r3)                                      # and make a masked array    
    return ifgs_r3_ma


#%%

# def get_baseline_end_ifg_n(LiCSBAS_imdates, baseline_end):
#     """ Given a list of the dates that there are steps in the LICSBAS time series for (i.e. when there was a Sentinel-1 acquisition),
#     find which number is the last before the baseline stage ends.  Note the baseline stage can end on any date (i.e. not one when there's an acquisition)
                                                                                                                
#     The number returned is the imdate number that is the last in the baseline, starting counting at 0.  If you wanted to index imdates and include this,
#     you would have to do date_n + 1.  
    
#     Inputs:
#         LiCSBAS_imdates | list of strings | dates of Sentinel-1 acquisitions, in form YYYYMMDD
#         baseline_end | string | dates of end of baseline stage, in form YYYYMMDD
#     Returns:
#         date_n | int | the number of the last date_n that is in the baseline stage, starting counting at 0.  
#     History:
#         2020/11/25 | MEG | Written
    
#     """
#     acq_after_baseline = False                                                                  # initiate as False
#     date_n = 0                                                                                  # initate the counter
#     while acq_after_baseline is False:                                                          # while loop
#         acq_after_baseline = compare_two_dates(baseline_end, LiCSBAS_imdates[date_n])           # check if imdate is after baseline_end
#         if acq_after_baseline is False:
#             date_n += 1                                                                         # if not, update the counter
#         else:
#             date_n -= 1                                                                         # but if it is, go back one to get the ifg number of the last ifg in the baseline stage (ie before it switched to True).  
#     return date_n

#%%

# Python version of Tee used to output print functions to the terminal and a log file.  Taken from stack exchange.  
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

#%%

def LiCSAR_ifgs_to_s1_acquisitions(LiCSAR_ifgs):
    """ Given a list of LiCSAR ifgs, determine the Sentinel-1 acquisition dates.  
    Inputs:
        LiCSAR_ifgs | list of strings | e.g. [20200306_20200312, 20200306_20200318]
    Returns:
        s1_acquisitions | list of strings 
    History:
        2020/11/17 | MEG | Written
    """
    s1_acquisitions = []
    for LiCSAR_ifg in LiCSAR_ifgs:
        date1 = LiCSAR_ifg[:8]
        date2 = LiCSAR_ifg[9:]
        for date in [date1, date2]:
            if date not in s1_acquisitions:
                s1_acquisitions.append(date)
    return s1_acquisitions

#%%

def baselines_from_ifgnames(names_list):
    """Given a list of ifg names in the form YYYYMMDD_YYYYMMDD, find the temporal baselines in days_elapsed (e.g. 12, 6, 12, 24, 6 etc.  )
    Inputs:
        names_list | list | in form YYYYMMDD_YYYYMMDD
    Returns:
        baselines | list of ints | baselines in days
    History:
        2020/02/16 | MEG | Documented
    """

    from datetime import datetime, timedelta

    baselines = []
    for file in names_list:

        master = datetime.strptime(file.split('_')[-2], '%Y%m%d')
        slave = datetime.strptime(file.split('_')[-1][:8], '%Y%m%d')
        baselines.append(-1 *(master - slave).days)
    return baselines

#%%

def compare_two_dates(date1, date2, fmt ='%Y%m%d'):
        """ Given two dates as strings (by default in YYYYMMDD format), determine if date 2 is after date 1.  
        Inputs:
            date1 | string | date, usually in format YYYYMMDD, but another format could be used if the fmt option is used.  
            date2 | string | date, usually in format YYYYMMDD, but another format could be used if the fmt option is used.  
            fmt | string | format used by date1 and date 2.  For formats, try: https://www.journaldev.com/23365/python-string-to-datetime-strptime
        Returns:
            date2_after | boolean | True is date2 is after date1
        History:
            2020/11/17 | MEG | Written from an existing script.  
        """
        import datetime

        date1_dt = datetime.datetime.strptime(date1, fmt)               # convert from string to datetime
        date2_dt = datetime.datetime.strptime(date2, fmt)               # ditto
        if date2_dt > date1_dt:                                         # see if date2 is later in time than date1
            date2_after = True
        else:
            date2_after = False
        return date2_after


#%%

def create_folder(folder):
    """ Try to create a folder to save function outputs.  If folder already exists,
    funtion will try to delete it and its contents.  
    
    Inputs: 
        folder | string |path to new folder.  e.g. './my_folder' 
    Returns:
        new folder
    History:
        2020/06/25 | MEG | Written
    """
    import shutil
    import os
    try:
        print(f"Trying to remove the existing outputs folder ({folder})... ", end = '')
        shutil.rmtree(folder)                                                                       # try to remove folder
        print('Done!')
    except:
        print("Failed!")                                                                                # 
    try:
        print(f"Trying to create a new outputs folder ({folder})... ", end = '')                                    # try to make a new folder
        os.mkdir(folder)                                                                       
        print('Done!')
    except:
        print("Failed!") 
        
        
        #%%


def col_to_ma(col, pixel_mask):
    """ A function to take a column vector and a 2d pixel mask and reshape the column into a masked array.  
    Useful when converting between vectors used by BSS methods results that are to be plotted
    
    Inputs:
        col | rank 1 array | 
        pixel_mask | array mask (rank 2)
        
    Outputs:
        source | rank 2 masked array | colun as a masked 2d array
    
    2017/10/04 | collected from various functions and placed here.  
    
    """
    import numpy.ma as ma 
    import numpy as np
    
    source = ma.array(np.zeros(pixel_mask.shape), mask = pixel_mask )
    source.unshare_mask()
    source[~source.mask] = col.ravel()   
    return source


def add_square_plot(x_start, x_stop, y_start, y_stop, ax, colour = 'k'):
    """Draw localization square around an area of interest, x_start etc are in pixels, so (0,0) is top left.  
    Inputs:
        x_start | int | start of box
        x_stop | int | etc. 
        y_start | int |
        y_ stop | int |
        ax | axes object | axes on which to draw
        colour | string | colour of bounding box.  Useful to change when plotting labels, and predictions from a model.  
    
    Returns:
        box on figure
        
    History:
        2019/??/?? | MEG | Written
        2020/04/20 | MEG | Document, copy to from small_plot_functions to LiCSAlert_aux_functions
    """
        
    ax.plot((x_start, x_start), (y_start, y_stop), c= colour)           # left hand side
    ax.plot((x_start, x_stop), (y_stop, y_stop), c= colour)             # bottom
    ax.plot((x_stop, x_stop), (y_stop, y_start), c= colour)             # righ hand side
    ax.plot((x_stop, x_start), (y_start, y_start), c= colour)             # top
    
    
    

#%%

def find_nearest_date(given_date, date_list):
    """  Given a list of dates in the form yyyymmdd and a date in the form
    yyyymmdd, find the date nearest to the given one.  
    """
    from datetime import datetime
    given_datetime = datetime.strptime(given_date, '%Y%m%d')
    date_list = [datetime.strptime(date, '%Y%m%d') for date in date_list]
    time_diffs = [abs(given_datetime - date) for date in date_list]
    nearest_date_index = time_diffs.index(min(time_diffs))

    return date_list[nearest_date_index].strftime('%Y%m%d')