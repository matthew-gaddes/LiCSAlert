#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:15:16 2024

@author: matthew
"""

import pdb

#%%

def extract_licsalert_status(volc_dirs, volc_names, day_list):
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
    
    def remove_empty_licsalert_volcs(licsalert_status, volc_names, volc_dirs):
        """  Some ovlcanoes do not have any licsalert outputs so just have a column of 0s.  
        Remove these.  
        Inputs:
            licsalert_status | r2 array | n_times x n_volcs x 2
            volc_names | list of strings |
            volc_dirs | list of paths as strings |
        Returns:
            As above, but with some removed.  
        History:
            2023_09_04 | MEG | Written
        """
        from copy import deepcopy
        
        licsalert_status_crop = deepcopy(licsalert_status)
        volc_names_crop = deepcopy(volc_names)
        volc_dirs_crop = deepcopy(volc_dirs)
        
        # get the indexes of ones to delete (that are just 0s)
        volc_del_indexes = []
        for volc_n in range(len(volc_names)):
            if np.sum(licsalert_status[:, volc_n]) == 0:
                #print(f"Noting for {volc_dirs[volc_n]}")
                volc_del_indexes.append(volc_n)
                
        # do the deleting    
        licsalert_status_crop = np.delete(licsalert_status_crop, volc_del_indexes, axis = 1)                            # each volcano is an item in dimension 1 (columns)
        for volc_n in sorted(volc_del_indexes, reverse = True):
            del volc_names_crop[volc_n]
            del volc_dirs_crop[volc_n]
                    
        return licsalert_status_crop, volc_names_crop, volc_dirs_crop
            
        
    
    # initiliase as empty    
    licsalert_status = np.zeros((len(day_list), len(volc_dirs), 2))                            
    
    for volc_n, volc_dir in enumerate(volc_dirs):
        date_dirs = sorted([d for d in glob(str(Path(volc_dir) / '*')) if os.path.isdir(d)])                    # get just the directories for that volcano (could be licsalert dates, could be ICASAR etc.)
        for date_dir in date_dirs:                                                                                  # loop through assuming that all are dates
            if (Path(date_dir).parts[-1] == 'ICASAR_results') or (Path(date_dir).parts[-1] == 'aux_figures'):           # but skip the occasional one that isn't a licsalert date directory
                pass
            else:
                try:
                    licsalert_date =  datetime.strptime(Path(date_dir).parts[-1], '%Y%m%d')                              #  get the date of the directory that we're in
                    day_n = day_list.index(licsalert_date)                                                              # find which row number this will be in the big outut array.   
                    
                    f = open(Path(date_dir) / "volcano_status.txt" , "r")                                               # get the volcano status
                    licsalert_status[day_n, volc_n, 0] = float(f.readline()[:-1])                                       # sigmas for existing deformation
                    licsalert_status[day_n, volc_n, 1] = float(f.readline()[:-1])                                       # sigmas for new deformaiton
                    f.close()
                    print(f"Found a licsalert status in {Path(date_dir).parts[-2:]} ")
                except:
                     print(f"Failed to find a licsalert status in {Path(date_dir).parts[-2:]} ")
                     
    # remove any volcanes that there are no licsalert products for.                       
    output = remove_empty_licsalert_volcs(licsalert_status, volc_names, volc_dirs)                      
    (licsalert_status_crop, volc_names_crop, volc_dirs_crop) = output; del output
    
          
    return licsalert_status_crop, volc_names_crop, volc_dirs_crop
 


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
    
    # get hte paths to each frame licsalert resykts
    volc_dirs = []
    if regions:
        region_dirs = sorted(glob(str(licsalert_dir / '*')))
        for region_dir in region_dirs:
            volc_dirs.extend(sorted(glob(str(Path(region_dir) / '*'))))
    else:
        volc_dirs = sorted(glob(str(licsalert_dir / '*')))
    
    if len(volc_dirs) == 0:
        raise Exception(f"No licsalert directories were found for any volcanoes.  "
                        f"Perhaps 'licsalert_dir' is set incorrectly?  It is "
                        f"currently set to: {licsalert_dir} \n"
                        f"Exiting.  ")
        
    # possibly remove some volcanoes
    if omits is not None:
        for omit in omits:
            for volc_dir in volc_dirs:
                # get the frame name from the full path
                if Path(volc_dir).parts[-1] == omit:
                    # and remove it in the omit list
                    volc_dirs.remove(volc_dir)
        
    # get just the volcano name
    volc_names = []
    for volc_dir in volc_dirs:
        volc_names.append(Path(volc_dir).parts[-1])
        
    

    return volc_dirs, volc_names


    


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

