#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:34:45 2024

@author: matthew
"""
import pdb

#%%

def create_day_list(d_start, d_stop):
    """
    """
    from datetime import datetime, timedelta
    import numpy as np
    
    dstart = datetime.strptime(d_start, '%Y%m%d')                              # 
    dstop = datetime.strptime(d_stop, '%Y%m%d')                                # 
    
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
    """ Given a list of datetimes, get the temporal baselines relative to the first one.  
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

