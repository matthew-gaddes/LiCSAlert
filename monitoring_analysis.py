#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:28:07 2023

@author: matthew

# initialise something to store:


  dim 1: day
  dim 2; volcano
  dim 3: which value

  plus: 
    list of days
    list of volcano names.  

  # need to get number of days 
  # need to get number of volcanoes.  


# get the region directories.  
# loop through these.  
  # get the volcano directories.  
  # loop through these.  


"""

import numpy as np
import matplotlib.pyplot as plt 
from pathlib import Path
import pdb
import os
from glob import glob

licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/01_test/")
#licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/02_all/")




#%% debug scripts 


def matrix_show(matrix, title=None, ax=None, fig=None, save_path = None, vmin0 = False):
    """Visualise a matrix
    Inputs:
        matrix | r2 array or masked array
        title | string
        ax | matplotlib axes
        save_path | string or None | if a string, save as a .png in this location.  
        vmin0 | boolean | 
        

    2017/10/18 | update so can be passed an axes and plotted in an existing figure
    2017/11/13 | fix bug in how colorbars are plotted.
    2017/12/01 | fix bug if fig is not None
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if ax is None:
        fig, ax = plt.subplots()
    matrix = np.atleast_2d(matrix)                   # make at least 2d so can plot column/row vectors

    if isinstance(matrix[0,0], np.bool_):           # boolean arrays will plot, but mess up the colourbar
        matrix = matrix.astype(int)                 # so convert

    if vmin0:
        matrixPlt = ax.imshow(matrix,interpolation='none', aspect='auto', vmin = 0)
    else:
        matrixPlt = ax.imshow(matrix,interpolation='none', aspect='auto')
    fig.colorbar(matrixPlt,ax=ax)
    if title is not None:
        ax.set_title(title)
        fig.canvas.manager.set_window_title(f"{title}")

    if save_path is not None:                                                   # possibly save the figure
        if title is None:                                                       # if not title is supplied, save with a default name
            fig.savefig(f"{save_path}/matrix_show_output.png")
        else:
            fig.savefig(f"{save_path}/{title}.png")                             # or with the title, if it's supplied 
            
    
    plt.pause(1)                                                                    # to force it to be shown when usig ipdb



#%%


def get_all_volcano_dirs(licsalert_dir):
    """ Get the paths to all the licsalert volcano dirs that are split across regions.  
    """
    volc_dirs = []
    #region_dirs = os.listdir(licsalert_path)
    region_dirs = sorted(glob(str(licsalert_dir / '*')))
    for region_dir in region_dirs:
    # for region_dir in region_dirs[:2]:
        # print(f"JUST LOOKING IN 2nd DIR")
        volc_dirs.extend(sorted(glob(str(Path(region_dir) / '*'))))
        
    volc_names = []
    for volc_dir in volc_dirs:
        volc_names.append(Path(volc_dir).parts[-1])
        
    return volc_dirs, volc_names
    



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
            del volc_dirs[volc_n]
                    
        return licsalert_status_crop, volc_names_crop, volc_dirs_crop
            
        
    
    from datetime import datetime
    licsalert_status = np.zeros((len(day_list), len(volc_dirs), 2))                            # initiliase as empty
    
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
                     
    licsalert_status_crop, volc_names_crop, volc_dirs_crop = remove_empty_licsalert_volcs(licsalert_status, volc_names, volc_dirs)                      # remove any volcanes that there are no licsalert products for.  
                     
    return licsalert_status_crop, volc_names_crop, volc_dirs_crop
    


def find_volcano(volc_names, volc_name):
    """ Given a volcano name, find corresponding LiCSAR frames for that name.  
    Inputs:
        volc_names | list of strings | e.g.  ['tungnafellsjokull_147A_02466_191712','tungnafellsjokull_147A_02488_131211']
        volc_name | string | volcano to search for.  E.g. sierra_negra*.  Wildcards are allowed!
    Returns:
        name_indexes | list of tuples | LiCSAR grame and corresponding index.  
    History:
        2023_09_04 | MEG | Written
    """
    import fnmatch
    name_and_indexes = []
    possible_matches = fnmatch.filter(volc_names, volc_name)
    for possible_match in possible_matches:
        name_and_indexes.append((possible_match, volc_names.index(possible_match)))
    return name_and_indexes




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
    



def licsalert_ts_one_volc(name_and_index, licsalert_status, day_list):
    """ Given a volcano name and its index (col number in the licsalert matrix),
    plot all times for that volcano.  
    Inputs:
        name_and_index | tuple | name and index 
        licsalert_status | r3 array | n_times x n_volcs x 2
        day_list | list of datetimes | date for each row in licsalert matrix. 
    Returns:
        Figure
    History:
        2023_09_04 | MEG | Written
    """
    
    def day_list_to_baselines(day_list_crop):
        """ Given a list of datetimes, get the temporal baselines relative to the first one.  
        Inputs:
            day_list_crop | list of datetimes | 
        Returns:
            tbaselines | r1 array | baselines relatie to first date.  
        History:
            2023_09_04 | MEG | Written
        """
        tbaselines = np.zeros((len(day_list_crop)))
        for day_n, day in enumerate(day_list_crop):
            tbaselines[day_n] = (day - day_list_crop[0]).days
        return tbaselines
    
    volcano_status = licsalert_status[:, name_and_index[1], :]                          # index the 3d one for volcanoes x times x metric to just times x metric for one volcano
    volcano_status_crop, day_list_crop = remove_dates_with_no_status(volcano_status, day_list)      # most days don't have an acquisition.  Remove them.  
    tbaselines = day_list_to_baselines(day_list_crop)                                               # get the temporal baselines in days
    
    
    # one static plot
    from datetime import datetime as dt
    f, ax = plt.subplots()
    ax.scatter(volcano_status_crop[:,0], volcano_status_crop[:,1], c = tbaselines)
    ax.set_xlabel('Sigma for exising def.')
    ax.set_ylabel('Sigma for new def.')
    ax.set_ylim([0,5])
    ax.set_xlim([0,5])
    
    label_ratio = 0.2
    x_threshold = (1 - label_ratio) * np.max(volcano_status_crop[:,0])
    y_threshold = (1 - label_ratio) * np.max(volcano_status_crop[:,1])
    for point_n, status in enumerate(volcano_status_crop):
        if (status[0] > x_threshold) or (status[1] > y_threshold):
            ax.annotate(dt.strftime(day_list_crop[point_n], '%Y%m%d'), (status[0], status[1]))
    
    

#%% Initialise


volc_dirs, volc_names = get_all_volcano_dirs(licsalert_dir)                             # get path to outputs for each volcanoes, and their names (snake case)
    
day_list = create_day_list("20200101", "20230901")                                      # make a datetime for each day in region of interest.  

licsalert_status, volc_names_crop, volc_dirs_crop = extract_licsalert_status(volc_dirs, volc_names, day_list)           

#%% figure for one volcano.  

        
# name_and_indexes = find_volcano(volc_names_crop, 'sierra_negra*')                         # couldn't find any.  
name_and_indexes = find_volcano(volc_names_crop, 'ale_bagu*')
     
        
    
licsalert_ts_one_volc(name_and_indexes[-1], licsalert_status, day_list)    


#%%

# def 2d_gif_for_all_volcs(licsalert_status, volc_names_crop, day_list):
#     """
#     """
    
    
#     n_frames = len(day_list)



    
    

    
    # 2nd pane is licalert figure? read correct one and display.  
    # make gif.  
    
    # animation for all times.  
    # f, ax = plt.subplots()
    # ax.set_ylim([0,5])
    # ax.set_xlim([0,5])
    # for date_n, status in enumerate(volcano_status_crop):
    #     ax.scatter(status[0], status[1], c = tbaselines[date_n])
    #     plt.show()
    #     plt.pause(1)
        
        
    

                            

        
#%%



    

    
    
    
    
   
    
    
    
    
matrix_show(licsalert_status[:,:,0], title = 'change to def')    

matrix_show(licsalert_status[:,:,1], title = 'new  def')    
    
    
    
    
    
    
    
    
    
    
    
    
    





