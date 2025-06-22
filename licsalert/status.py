#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:15:16 2024

@author: matthew
"""

import pdb



#%% get_volc_names_fron_dir_of_frames()

def get_volc_names_fron_dir_of_frames(licsalert_dir, regions = False):
    """ Given a directory of lics volcano frames (possibly organised 
    into regions), get the names of the volcanoes present.  
    e.g. if there's volcano_156D_12345_123456 and volcano_149A_12345_123456',
    it will return ['volcano']
    
    Inputs:
        licsalert_dir | pathlib Path | location of licbas (or licsalert) 
                                        volcano frames.  
        regions | boolean | If true, the above directories are organised into 
                            regions.  
            
    Returns:
        volc_names | list | human readable volcano names.  e.g.:
                            ['corcovado', 'huequi', 'callaqui']
                            
    History:
        2024_11_15 | MEG | Written.  
    """
    
    volc_frame_paths, volc_frame_names = get_all_volcano_dirs(
        licsalert_dir,  regions = regions
        )
    
    
    volc_names = []
    
    for volc_frame_name in volc_frame_names:
        camel_name = volc_frame_name[:-18]
        # remove underscores
        human_name = camel_name.replace('_', ' ')
        volc_names.append(human_name)
        
    # remove duplicates
    volc_names = list(set(volc_names))
    
    return volc_names
        


#%% extract_licsalert_status()



def extract_licsalert_status(volcs, day_list, combined_status_method='window'):
    """ For every volcano that we do licsalert for, extract the 2 status 
    values for all possible times.  
    
    Also calculate the combined status from these values.  
    
    Also calculate the combined unrest metric (which is currently just a 
    the maximum of the new / existing deformation metrics for each frame)
    
    Inputs:
        volc_dirs | list of paths as strings | 
        day_list | list of datetimes | consecutive days.  
        
    Returns:
        licsalert_status | r3 array | n_times x n_volcs x 2
    History:
        2023_09_04 | MEG | Adapted from a sript.  
    """
    import numpy as np
    from glob import glob
    from pathlib import Path
    import os
    from datetime import datetime
    from copy import deepcopy
    from licsalert.aux import get_licsalert_date_dirs

    # loop through each volcano    
    for volc_n, volc in enumerate(volcs):
        
        # initiliase dict to store the status for each frame at each time
        volc.status = {}
        
        # only try to make a licsalert status if there are licsalert resultsq
        if volc.processing_status.licsalert_result == True:
    
            # and the frames of that volcano (only the ones with data though)
            for frame_n, frame_path in enumerate(volc.frame_status):
                           
                # check that there is data for that frame
                if frame_path != 'NA':
                    print(f"Creating the LiCSAlert status for frame {volc.frames[frame_n]}")
                    # create to store output for this frame
                    frame_status = {'dates' : [],
                                    'day_ns' : [],
                                    'existing_defs' : [],
                                    'new_defs'      :[]}
                    # loop through all licsalert dates
                    for date_dir in get_licsalert_date_dirs(frame_path):
                        # and get the licsalert status for that date
                        try:
                            licsalert_date = datetime.strptime(Path(date_dir).parts[-1], '%Y%m%d')
                            day_n = day_list.index(licsalert_date)
                            # read the status from the .txt file.  
                            f = open(Path(frame_path)/date_dir/"volcano_status.txt" , "r")
                            existing_def = float(f.readline()[:-1])
                            new_def = float(f.readline()[:-1])
                            
                            # if no error at this point, assume OK to record values
                            frame_status['existing_defs'].append(existing_def)
                            frame_status['new_defs'].append(new_def)
                            
                            frame_status['dates'].append(licsalert_date)
                            frame_status['day_ns'].append(day_n)
                            
                            f.close()
                        except:
                            pass
                    
                    # add the status for this frame to the dict for all frames
                    volc.status[f"{volc.frames[frame_n][-17:]}"] = frame_status
                    
            # after we have processed all frames for that volcano, 
            # combine the statuses to make the combined status.  
            # (that is a a new def and existing def metric for each time)
            # 2 x times for each volcano
            volc.status_combined = calculate_status_combined(
                volc, day_list, combined_status_method
                )
            
            # also calculate the overall status
            # (that is the maximum of new and existing for each volcano)
            # 1 x times 
            calculate_status_overall(volc)
            
            # calcaulte the integral of the overall status    
            calculate_status_cumulative(volc)
       

def calculate_status_combined(
        volc, day_list, combined_status_method = 'window' 
        ):
    """ For a single volcano, iterate along each day and calculate the 
    combined status from all frames for that day.  This can be done using one
    of two methods (previous, where we take the maximum of the previous value
    from all frames, which if the time series stops on a high value can be 
    months before, or window, where we take the maximum in a window.
    
    """

    # initiliase as empty
    status_combined = {'dates' : [],
                       'day_ns' : [],
                       'existing_defs' : [],
                       'new_defs'      :[]}
    # calculate for all days
    for day in day_list:
        # add the time info
        status_combined['dates'].append(day)
        day_n = day_list.index(day)
        status_combined['day_ns'].append(day_n)
        
        # collate across all frames for that volcano
        day_existing_defs = []
        day_new_defs = []
        
        # v1: take maximum of nearest value for each frame
        # which could many days or months prior 
        # value is the dict of the status for that volcano
        # (dates, day_ns, existing_defs, new_defs)
        if combined_status_method == 'previous':
            day_existing_def, day_new_def = combined_status_previous(
                volc, day, day_n
                )
        # v2: take the maximum of a sliding window
        elif combined_status_method == 'window':
            
            day_existing_def, day_new_def = combined_status_window(
                volc, day, day_n, days_prior=12
                    )
        else:
            raise Exception(
                f"'combined_status_method' should either be 'window' or "
                f"'previous', but is {combined_status_method}.  Exiting.  ")
        
            
        # record the status, whichever way it was calculated.  
        status_combined['existing_defs'].append(day_existing_def)
        status_combined['new_defs'].append(day_new_def)
            
    # rerurn the dict that contains the combined status for that volcano (i.e.
    # across all frames) for all days.  
    return status_combined

   
def calculate_status_overall(volc):
    """ Calculate the overall status, which is just the maximum
    of the new deformation unrest and existing deformation unrest
    for each time for the volcano.  
    """
    from copy import deepcopy
    
    # make a deepcopy to edit
    volc.status_overall = deepcopy(volc.status_combined)
    
    # initialise
    volc.status_overall['existing_defs']
    
    
    # Take the maximum value from each pair of values to make a single
    # unrest value for all times 
    volc.status_overall['unrest_metric'] = [
        max(a, b) for a, b in zip(
            volc.status_overall['new_defs'], 
            volc.status_overall['existing_defs'], 
            )
        ]
    
    # remove these from the overall status now the unrest metric
    # has been calculated
    del volc.status_overall['new_defs']
    del volc.status_overall['existing_defs']
    

def combined_status_previous(volc, day, day_n):
    """ Calculate the combined status metric as the maximum of
    the previous values for each frame.  If a frame has a high 
    value and then the time series stops, this value will
    continue to be used indefinately, even though newer values
    are available in other frames.  
    """
    
    # collate across all frames for that volcano
    day_existing_defs = []
    day_new_defs = []
    
    for frame, value in volc.status.items():

        # find the nearest LiCSAlert status before the current day
        # first find the day number for the nearest previous day
        closest_previous_day = nearest_smaller_or_equal(
            value['day_ns'], day_n
            )
        
        # if no data yet (as time series for the frame hasn't started)
        if closest_previous_day is None:
            day_existing_defs.append(0)
            day_new_defs.append(0)
            
        # if no data as the time series for that frame has stopped
        elif day > value['dates'][-1]:
            # pdb.set_trace()
            day_existing_defs.append(0)
            day_new_defs.append(0)
            
        else:
            # closest_previous_day is done relative to the start of 
            # day_list.  Convert this to be day number relative to 
            # the start of the time series for that frame (the frame
            # may not go back to the start of day_list)
            closest_day_index = value['day_ns'].index(
                closest_previous_day
                )
            
            # record the two unrest metrics for that dat
            day_existing_defs.append(
                value['existing_defs'][closest_day_index]
                )
            day_new_defs.append(
                value['new_defs'][closest_day_index]
                )
            
    # these are the statuses for each frame on that day, take max across all frames
    return max(day_existing_defs), max(day_new_defs)




def combined_status_window(volc, day, day_n, days_prior=12):
    """
    Extracts existing_defs and new_defs from value_dict within a given interval.

    Parameters:
        value_dict (dict): Dictionary with 'dates', 'existing_defs', and 'new_defs'.
        target_date (str): End date in 'yyyymmdd' format.
        days_prior (int): Number of days prior to the target_date.

    Returns:
        dict: Dictionary with extracted 'existing_defs' and 'new_defs'.
    """
    from datetime import datetime, timedelta
    import numpy as np
    
    def extract_defs_in_interval(value_dict, target_date, days_prior):
        """
        Extracts existing_defs and new_defs from value_dict within a given interval.
    
        Parameters:
            value_dict (dict): Dictionary with 'dates', 'existing_defs', and 'new_defs'.
            target_date (str): End date in 'yyyymmdd' format.
            days_prior (int): Number of days prior to the target_date.
    
        Returns:
            dict: Dictionary with extracted 'existing_defs' and 'new_defs'.
        """
        dates = value_dict['dates']
        target_datetime = target_date
        start_datetime = target_datetime - timedelta(days=days_prior)
    
        existing_defs_interval = []
        new_defs_interval = []
    
        for idx, current_date in enumerate(dates):
            if start_datetime <= current_date <= target_datetime:
                existing_defs_interval.append(value_dict['existing_defs'][idx])
                new_defs_interval.append(value_dict['new_defs'][idx])
    
        return existing_defs_interval, new_defs_interval
    
    # collate across all frames for that volcano in the window
    day_existing_defs = []
    day_new_defs = []
    
    for frame, value in volc.status.items():
        # get the new and existing metrics for each frame
        existing_interval_1frame, new_interval_1frame= extract_defs_in_interval(
            value, day, days_prior
            )
        # and record
        day_existing_defs.extend(existing_interval_1frame)
        day_new_defs.extend(new_interval_1frame)
        

    # these are the statuses for each frame on that day, take max across all frames
    # (return nan if empty)
    return (
        max(day_existing_defs) if day_existing_defs else np.nan,
        max(day_new_defs) if day_new_defs else np.nan
    )



def calculate_status_cumulative(volc):
    """ Calculate the cumulative unrest metric
    (which is just the integral of the unrest_metric)
    """
    from copy import deepcopy
    import numpy as np
    
    # make a deepcopy to edit
    volc.status_cumulative = deepcopy(volc.status_overall)

            
    volc.status_cumulative['unrest_cum'] = np.cumsum(
        volc.status_cumulative['unrest_metric']
        )
   
    # remove these from the copied version as now redundant
    del volc.status_cumulative['unrest_metric']
    
                

        
def nearest_smaller_or_equal(lst, target):
    """ Find the closest item in the list that is smaller or equal.  
    E.g. looking for day 12, options are [1, 7, 13], 7 will be returned    

    """
    smaller_numbers = [num for num in lst if num <= target]
    if smaller_numbers:
        return max(smaller_numbers)
    return None  # If there are no smaller numbers
    
   
# def extract_defs_in_interval(value_dict, target_dt, days_prior):
#     """
#     Extracts existing_defs and new_defs from value_dict within a given interval.

#     Parameters:
#         value_dict (dict): Dictionary with 'dates', 'existing_defs', and 'new_defs'.
#         target_date (str): End date in 'yyyymmdd' format.
#         days_prior (int): Number of days prior to the target_date.

#     Returns:
#         dict: Dictionary with extracted 'existing_defs' and 'new_defs'.
#     """
#     dates = [datetime.strptime(str(d), "%Y%m%d") for d in value_dict['dates']]
#     #target_dt = datetime.strptime(target_date, "%Y%m%d")
#     start_dt = target_dt - timedelta(days=days_prior)

#     existing_defs_interval = []
#     new_defs_interval = []

#     for idx, current_date in enumerate(dates):
#         if start_dt <= current_date <= target_dt:
#             existing_defs_interval.append(value_dict['existing_defs'][idx])
#             new_defs_interval.append(value_dict['new_defs'][idx])

#     return {
#         'existing_defs': existing_defs_interval,
#         'new_defs': new_defs_interval
#     }
        

#%%


def remove_dates_with_no_status(volcano_status, day_list):
    """ Given the licsalert status for one volcano, remove any days that don't have an entry.  
    As S1 acquisitoins are every 6/12/24 days etc, this is most days.  
    Inputs:
        volcano_status | r2 array | n_times x 2, where first column is sigma for existing deformtaion, and 2nd is sigma for new deformaiton. 
        day_list | list of datetimes | 
    Returns:
        as above, but with dates with no status removed.  
    History:
        2023_09_04 | MEG | Written
    """
    
    import numpy as np
    from copy import deepcopy
    
    volcano_status_crop = deepcopy(volcano_status)
    day_list_crop = deepcopy(day_list)
    
    
    # get the indexes of ones to delete (that are just 0s)
    day_del_indexes = []
    for day_n in range(len(day_list)):
        if np.sum(volcano_status[day_n, :]) == 0:
            #print(f"Noting for {volc_dirs[volc_n]}")
            day_del_indexes.append(day_n)
            
    # do the deleting    
    volcano_status_crop = np.delete(volcano_status_crop, day_del_indexes, axis = 0)                            # each volcano is an item in dimension 1 (columns)
    for day_del_index in sorted(day_del_indexes, reverse = True):
        del day_list_crop[day_del_index]
    
    return volcano_status_crop, day_list_crop
       
 

#%%



def get_all_volcano_dirs(licsalert_dir, omits = None, regions = True):
    """ Get the paths to all the licsalert volcano dirs that are split across regions.  
    Inputs:
        licsalert_dir | pathlib Path | directory LiCSAlert outputs, 
                                        possibly containing a subdirectory of regions.  
        omit | list or None | List of volcanoes to remove.  
        regions | boolean | If True, volcanoes are separated into region 
                            directories, as per the COMET Volcano portal
    Returns:
        volc_dirs | list of pathlib Paths | list of each volcano frame's licsalert directory
        volc_names | list | each volcano frames name
    """
    from glob import glob
    from pathlib import Path
    import os
    
    
    def contains_licsalert_date_dir(directories):
        """
        """
        
        import re
        
        # Regular expression for eight digits
        pattern = re.compile(r'^\d{8}$')
        
        for directory in directories:
            # Extract the directory name
            dir_name = os.path.basename(directory)
            if pattern.match(dir_name):
                return True
        return False
    
    
    # get the paths to each frame licsalert resykts
    volc_frame_paths = []
    # make a list to iterate through, which depends on if there are regions.  
    if regions:
        licsalert_dirs = sorted(glob(str(licsalert_dir / '*')))
    else:
        licsalert_dirs = [licsalert_dir]
        
        
    
    # iterate through each region
    for dir1 in licsalert_dirs:
        
        # return a list of volcano frames
        volc_dirs = sorted(glob(str(Path(dir1) / '*')))
        
        # iterate through each volcano 
        for volc_dir in volc_dirs:
            licsalert_products = sorted(glob(str(Path(volc_dir) / '*')))
            # get only the directories
            licsalert_dirs = [item for item in licsalert_products if os.path.isdir(item)]
            
            # see if there's a licsalert date dir (so it ran succesfully)
            licsalert_success = contains_licsalert_date_dir(licsalert_dirs)

            # if licsalert ran for this volcano frame, add it to the list            
            if licsalert_success:
                volc_frame_paths.append(volc_dir)
                
    
    if len(volc_frame_paths) == 0:
        raise Exception(f"No licsalert directories were found for any volcanoes.  "
                        f"Perhaps 'licsalert_dir' is set incorrectly?  It is "
                        f"currently set to: {licsalert_dir} \n"
                        f"Exiting.  ")
        
    # possibly remove some volcanoes
    if omits is not None:
        for omit in omits:
            for volc_frame in volc_frame_paths:
                # get the frame name from the full path
                if Path(volc_frame).parts[-1] == omit:
                    # and remove it in the omit list
                    volc_frame_paths.remove(volc_frame)
        
    # get just the volcano name
    volc_frame_names = []
    for volc_frame in volc_frame_paths:
        volc_frame_names.append(Path(volc_frame).parts[-1])
        
    return volc_frame_paths, volc_frame_names


    


#%%

def find_volcano(volc_frames, volc_name):
    """ Given a volcano name, find corresponding LiCSAR frames for that name.  
    Inputs:
        volc_frames | list of strings | e.g.  ['tungnafellsjokull_147A_02466_191712','tungnafellsjokull_147A_02488_131211']
        volc_name | string | volcano to search for.  E.g. sierra_negra*.  Wildcards are allowed!
    Returns:
        name_indexes | list of tuples | LiCSAR grame and corresponding index.  
    History:
        2023_09_04 | MEG | Written
    """
    import fnmatch
    name_and_indexes = []
    possible_matches = fnmatch.filter(volc_frames, volc_name)
    for possible_match in possible_matches:
        name_and_indexes.append((possible_match, volc_frames.index(possible_match)))
    return name_and_indexes

