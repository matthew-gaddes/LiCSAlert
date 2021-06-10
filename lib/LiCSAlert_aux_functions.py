#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 18:10:21 2020

@author: matthew
"""

#%%

def get_baseline_end_ifg_n(LiCSBAS_imdates, baseline_end):
    """ Given a list of the dates that there are steps in the LICSBAS time series for (i.e. when there was a Sentinel-1 acquisition),
    find which number is the last before the baseline stage ends.  Note the baseline stage can end on any date (i.e. not one when there's an acquisition)
                                                                                                                
    The number returned is the imdate number that is the last in the baseline, starting counting at 0.  If you wanted to index imdates and include this,
    you would have to do date_n + 1.  
    
    Inputs:
        LiCSBAS_imdates | list of strings | dates of Sentinel-1 acquisitions, in form YYYYMMDD
        baseline_end | string | dates of end of baseline stage, in form YYYYMMDD
    Returns:
        date_n | int | the number of the last date_n that is in the baseline stage, starting counting at 0.  
    History:
        2020/11/25 | MEG | Written
    
    """
    acq_after_baseline = False                                                                  # initiate as False
    date_n = 0                                                                                  # initate the counter
    while acq_after_baseline is False:                                                          # while loop
        acq_after_baseline = compare_two_dates(baseline_end, LiCSBAS_imdates[date_n])           # check if imdate is after baseline_end
        if acq_after_baseline is False:
            date_n += 1                                                                         # if not, update the counter
        else:
            date_n -= 1                                                                         # but if it is, go back one to get the ifg number of the last ifg in the baseline stage (ie before it switched to True).  
    return date_n

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