#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:34:45 2024

@author: matthew
"""
import pdb

#%%

def daisy_chain_from_acquisitions(acquisitions):
    """Given a list of acquisiton dates, form the names of the interferograms 
    that would create a simple daisy chain of ifgs.  
    Inputs:
        acquisitions | list | list of acquistiion dates in form YYYYMMDD
    Returns:
        daisy_chain | list | names of daisy chain ifgs, in form 
        YYYYMMDD_YYYYMMDD
    History:
        2020/02/16 | MEG | Written
    """
    daisy_chain = []
    n_acqs = len(acquisitions)
    for i in range(n_acqs-1):
        daisy_chain.append(f"{acquisitions[i]}_{acquisitions[i+1]}")
    return daisy_chain

#%%

def create_day_list(d_start, d_stop):
    """
    """
    from datetime import datetime, timedelta
    import numpy as np
    
    dstart = datetime.strptime(d_start, '%Y%m%d')
    dstop = datetime.strptime(d_stop, '%Y%m%d')
    
    acq_dates = [dstart]
    dcurrent = acq_dates[-1]

    while dcurrent < dstop:
        dnext = dcurrent + timedelta(days = 1)                           # add the temp baseline to find the new date   
        if dnext < dstop:                                                        # check we haven't gone past the end date 
            acq_dates.append(dnext)                                              # if we haven't, record            
            dcurrent = acq_dates[-1]                                             # record the current date
        else:
            break                                                                # remember to exit the while if we have got to the last date

    return acq_dates



#%%


def day_list_to_baselines(day_list_crop):
    """ Given a list of datetimes, get the temporal baselines relative to the 
    first one.  
    Inputs:
        day_list_crop | list of datetimes | 
    Returns:
        tbaselines | r1 array | baselines relatie to first date.  
    History:
        2023_09_04 | MEG | Written
    """
    import numpy as np 
    
    tbaselines = np.zeros((len(day_list_crop)))
    for day_n, day in enumerate(day_list_crop):
        tbaselines[day_n] = (day - day_list_crop[0]).days
    return tbaselines


#%%

def baseline_from_names(names_list):
    """Given a list of ifg names in the form YYYYMMDD_YYYYMMDD, find the 
    temporal baselines in days_elapsed
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