#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:15:16 2024

@author: matthew
"""

import pdb

#%%



#%%

def extract_licsalert_status(volcs, day_list):
    """ For every volcano that we do licsalert for, extract the 2 status values for all possible times.  
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
    
    from licsalert.aux import get_licsalert_date_dirs
    
            
    def nearest_smaller(lst, target):
        smaller_numbers = [num for num in lst if num < target]
        if smaller_numbers:
            return max(smaller_numbers)
        return None  # If there are no smaller numbers
    
    # def remove_empty_licsalert_volcs(licsalert_status, volc_names, volc_dirs):
    #     """  Some ovlcanoes do not have any licsalert outputs so just have a column of 0s.  
    #     Remove these.  
    #     Inputs:
    #         licsalert_status | r2 array | n_times x n_volcs x 2
    #         volc_names | list of strings |
    #         volc_dirs | list of paths as strings |
    #     Returns:
    #         As above, but with some removed.  
    #     History:
    #         2023_09_04 | MEG | Written
    #     """
    #     from copy import deepcopy
        
    #     licsalert_status_crop = deepcopy(licsalert_status)
    #     volc_names_crop = deepcopy(volc_names)
    #     volc_dirs_crop = deepcopy(volc_dirs)
        
    #     # get the indexes of ones to delete (that are just 0s)
    #     volc_del_indexes = []
    #     for volc_n in range(len(volc_names)):
    #         if np.sum(licsalert_status[:, volc_n]) == 0:
    #             #print(f"Noting for {volc_dirs[volc_n]}")
    #             volc_del_indexes.append(volc_n)
                
    #     # do the deleting    
    #     licsalert_status_crop = np.delete(licsalert_status_crop, volc_del_indexes, axis = 1)                            # each volcano is an item in dimension 1 (columns)
    #     for volc_n in sorted(volc_del_indexes, reverse = True):
    #         del volc_names_crop[volc_n]
    #         del volc_dirs_crop[volc_n]
                    
    #     return licsalert_status_crop, volc_names_crop, volc_dirs_crop
    
    # loop through each volcano    
    for volc_n, volc in enumerate(volcs):
        # and the frames of that volcano (only the ones with data though)
        volc.status = {}
        for frame_n, frame_path in enumerate(volc.frame_status):
            # check that there is data for that frame
            if frame_path != 'NA':
                print(f"Creating the LiCSAlert status for frame {volc.frames[frame_n]}")
                # create to store output
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
                volc.status[f"{volc.frames[frame_n][-17:]}"] = frame_status
                
    
        # combine the statuses
        combined_status = {'dates' : [],
                           'day_ns' : [],
                           'existing_defs' : [],
                           'new_defs'      :[]}
        for day in day_list:
            combined_status['dates'].append(day)
            day_n = day_list.index(day)
            combined_status['day_ns'].append(day_n)
            
            # collate across all frames for that volcano
            day_existing_defs = []
            day_new_defs = []
            for frame, value in volc.status.items():
                closest_previous_day = nearest_smaller(value['day_ns'], day_n)
                # if no data yet
                if closest_previous_day is None:
                    day_existing_defs.append(0)
                    day_new_defs.append(0)
                else:
                    closest_day_index = value['day_ns'].index(closest_previous_day)
                    day_existing_defs.append(value['existing_defs'][closest_day_index])
                    day_new_defs.append(value['new_defs'][closest_day_index])
        
            combined_status['existing_defs'].append(max(day_existing_defs))
            combined_status['new_defs'].append(max(day_new_defs))
            
        volc.combined_status = combined_status
                

            
            

            
            
            
                
        
        
    

    

                



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
    volc_frames = []
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
                volc_frames.append(volc_dir)
                
    
    if len(volc_frames) == 0:
        raise Exception(f"No licsalert directories were found for any volcanoes.  "
                        f"Perhaps 'licsalert_dir' is set incorrectly?  It is "
                        f"currently set to: {licsalert_dir} \n"
                        f"Exiting.  ")
        
    # possibly remove some volcanoes
    if omits is not None:
        for omit in omits:
            for volc_frame in volc_frames:
                # get the frame name from the full path
                if Path(volc_frame).parts[-1] == omit:
                    # and remove it in the omit list
                    volc_frames.remove(volc_frame)
        
    # get just the volcano name
    volc_names = []
    for volc_frame in volc_frames:
        volc_names.append(Path(volc_frame).parts[-1])
        
    return volc_frames, volc_names


    


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

