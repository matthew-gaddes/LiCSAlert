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

#licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/01_test/")
licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/02_all/")




#%% from other modules


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




def pngs_to_gif(input_folder, output_file, image_duration = 0.5):
    """A function to conbime a folder of .pngs into one gif.
    Inputs:
        input_folder | pathlib Path | the location of the pngs.  E.g. "output_data".  Must be a pathlib Path
        output_file  | string | e.g. "output.gif"
        image_duration | float | the time each png is displayed for in seconds.  e.g. 0.5
    Returns:
        gif
    History:
        2018/07/23 | MEG: adapted from a script.
        2018/12/05 | paths must now be absolute.
        2020/08/18 | MEG | Update to handle paths formatted using the pathlib module.  
        2021_05_04 | MEG | Update to handle paths better.  
        2023_09_04 | MEG | Add funtion to make sure all images are the same size.  
    """

    import imageio
    import os
    from PIL import Image
    import re

    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [ atoi(c) for c in re.split('(\d+)', text) ]
    
    def determine_minimim_image_size(images):
        """ Iterate along the list of images and find the size of the smallest image
        """
        for image_n, image in enumerate(images):
            ny, nx, _ = image.shape
            if image_n == 0:
                min_y = ny
                min_x = nx
            else:
                if ny < min_y:
                    min_y = ny
                if nx < min_x:
                    min_x = nx
        return min_y, min_x

           
    def crop_if_larger(images, min_y, min_x):
        """ If an image is larger than the minimum, crop it.  
        """
        for image_n, image in enumerate(images):
            ny, nx, _ = image.shape
            if (ny > min_y) or (nx > min_x):
                print(f"Resizing image from {ny} to {min_y} in y, and {nx} to {min_x} in x ")
                images[image_n] = image[:min_y, :min_x, :]
        return images

    # sort the files into a human order
    image_names = os.listdir(input_folder)              # get 
    image_names.sort(key=natural_keys)                  # sort into the correct (human) order

    # write the images to the .gif
    images = []
    for counter, filename in enumerate(image_names):
        images.append(imageio.imread(input_folder / filename))                                  # input_folder must be a pathlib Path to allow simple extension with /
        print(f"Processed {counter} of {len(image_names)} images.  ")
    
    min_y, min_x = determine_minimim_image_size(images)                             # fnid the size of the minimum image
    images = crop_if_larger(images, min_y, min_x)                                   # make sure there are none bigger than that.  

            
    
    print('Writing the .gif file...', end = "")
    imageio.mimsave(output_file, images, duration = image_duration, loop = 1)
    print('Done!')



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


def licsalert_ts_one_volc(name_and_index, licsalert_status, day_list, fig_type = 'png',
                          volc_dirs = None, out_dir = None):
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
    
    from datetime import datetime as dt
    import matplotlib.image as mpimg
    import matplotlib.gridspec as gridspec
    
    volcano_status = licsalert_status[:, name_and_index[1], :]                          # index the 3d one for volcanoes x times x metric to just times x metric for one volcano
    volcano_status_crop, day_list_crop = remove_dates_with_no_status(volcano_status, day_list)      # most days don't have an acquisition.  Remove them.  
    tbaselines = day_list_to_baselines(day_list_crop)                                               # get the temporal baselines in days
    
    
    # one static plot
    if fig_type == 'png':
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
    
    elif fig_type == 'gif':
        plt.switch_backend('Agg')                                                           #  works when there is no X11 forwarding, and when displaying plots during creation would be annoying.  
        volc_dir = volc_dirs[name_and_index[1]]                         # get the directory of the volcano we're working with.  
        for date_n, status in enumerate(volcano_status_crop):
        
            f = plt.figure(figsize = (25,8))
            grid = gridspec.GridSpec(2, 6, wspace = 0., hspace = 0.)                        # divide into 2 sections, 1/5 for ifgs and 4/5 for components
            ax_2d = plt.Subplot(f, grid[-2,0])                                                                                      # create an axes for the IC (spatial source)
            ax_im = plt.Subplot(f, grid[:,1:])                                                                                      # create an axes for the IC (spatial source)
            ax_2d.set_ylim([0,5])
            ax_2d.set_xlim([0,5])
            ax_2d.set_xlabel('Sigma for exising def.')
            ax_2d.set_ylabel('Sigma for new def.')
        
            # ax_2d.scatter(status[0], status[1], c = tbaselines[date_n])
            ax_2d.scatter(volcano_status_crop[:date_n+1,0], volcano_status_crop[:date_n+1,1], c = tbaselines[:date_n+1])                # python excludes end when indexing so +1 to include
            ax_2d.set_aspect('equal')
            
            # Plot the corresponding licsalert figure
            date_string = f"{dt.strftime(day_list_crop[date_n], '%Y%m%d')}"
            date_dir = Path(volc_dir) / date_string
            licsalert_fig_path = glob(str(date_dir/ 'LiCSAlert_figure_with_*_monitoring_interferograms.png'))
            licsalert_fig = mpimg.imread(licsalert_fig_path[0])
            ax_im.imshow(licsalert_fig)
            ax_im.axis('off')
            
            # tidy up axes, save fig etc.  
            for ax in [ax_2d, ax_im]:
                f.add_subplot(ax)
            f.suptitle(date_string)
            
            f.savefig(out_dir / f"{date_string}.png", bbox_inches='tight')            
            print(f"Saved figure {date_n} of {volcano_status_crop.shape[0]} ")
    if plt.get_backend() != 'Qt5Agg':                                                           # check if we need to reset the backend (to interactive window)
        plt.switch_backend('Qt5Agg')                                                           #  

            
        

            


#%% Initialise


day_list = create_day_list("20200101", "20230901")                                      # make a datetime for each day in region of interest.  

volc_dirs, volc_names = get_all_volcano_dirs(licsalert_dir)                                               # get path to outputs for each volcanoes, and their names (snake case)
licsalert_status, volc_names, volc_dirs = extract_licsalert_status(volc_dirs, volc_names, day_list)           # get the licsalert status where valid, and shorten the volc_dirs and volc_names to only have them for where licsalert is valid.  

#%% figure for one volcano at all times.  

# name_and_indexes = find_volcano(volc_names_crop, 'sierra_negra*')                         # couldn't find any.  
name_and_indexes = find_volcano(volc_names, 'ale_bagu*')
     
licsalert_ts_one_volc(name_and_indexes[-1], licsalert_status, day_list, fig_type = 'png')                                                                           # 1 figure for all times
#licsalert_ts_one_volc(name_and_indexes[-1], licsalert_status, day_list, fig_type = 'gif', volc_dirs = volc_dirs, out_dir = Path('monitoring_steps_1volc')    )      # figures for all times, with corresponding licsaelrt
#pngs_to_gif(Path('monitoring_steps_1volc'), Path('monitoring_steps_1volc') / "animation.gif", image_duration = 0.5)                                                 # convert all times into a gif

#%% 

def licsalert_ts_all_volcs(licsalert_status, volc_names, day_list, out_dir, figsize, xlim = 10, ylim = 10):
    """
    
    One plot per day?  
    """
    
    from copy import deepcopy
    import matplotlib.pyplot as plt
    from datetime import datetime as dt
    from adjustText import adjust_text
    plt.switch_backend('Agg')                                                           #  
    # plt.switch_backend('QtAgg')                                                           #   interactive, used to debug.                    
    
    
    #n_frames = len(day_list)
    
    tbaselines = day_list_to_baselines(day_list)                                               # get the temporal baselines in days
    n_volcs = licsalert_status.shape[1]
    
    #licsalert_status_current = np.zeros((licsalert_status.shape[1], 2))                     # first dim is volcano number, 2nd is 
    licsalert_status_current  = deepcopy(licsalert_status[0,])
    day_counter = np.repeat(tbaselines[0], n_volcs)

    volcs_off_chart = []    
    volcs_off_chart_sigmas = []
    
    for day_n, day in enumerate(day_list):                                                                      # loop through all possible days
    
        
        # determine if some sigmas are larger than 5

        for volc_n in range(n_volcs):                                                                              # on each day, loop through all volcanoes     
            existing_def = licsalert_status[day_n, volc_n, 0]                                                   # get the number of sigmas for existing deformation
            new_def = licsalert_status[day_n, volc_n, 1]                                                    # and new deformation
            if  (existing_def != 0) and (new_def != 0):                                                     # check if we have values (and not a date in which they're just zeros) 
            
                if volc_names[volc_n] in volcs_off_chart:                                               # if we are dealing with a volcano that was previously off the plot
                    list_index = volcs_off_chart.index(volc_names[volc_n])                              # find where int he list it is
                    del volcs_off_chart[list_index], volcs_off_chart_sigmas[list_index]                 # so we can delete it as its location is to be updated
                    

                if (existing_def > xlim) and (new_def > ylim):                                                   # Option 1 of 4: too large in both
                    volcs_off_chart.append(volc_names[volc_n])
                    volcs_off_chart_sigmas.append((existing_def, new_def))
                    existing_def = xlim
                    new_def = ylim
                elif (existing_def > xlim):                                                                      # Option 2 of 4: too large in one
                    volcs_off_chart.append(volc_names[volc_n])
                    volcs_off_chart_sigmas.append((existing_def, new_def))
                    existing_def = xlim
                    new_def = new_def
                elif (new_def > ylim):                                                                           # Option 3 of 4: too large in other
                    volcs_off_chart.append(volc_names[volc_n])
                    volcs_off_chart_sigmas.append((existing_def, new_def))
                    new_def = ylim
                    existing_def = existing_def
                else:                                                                                           # Option 4 of 4: neither too large.  
                    new_def = new_def
                    existing_def = existing_def

                licsalert_status_current[volc_n,0]  = existing_def                                              # update the current status
                licsalert_status_current[volc_n,1]  = new_def                                                   # continued
                day_counter[volc_n] = tbaselines[day_n]                                                             # also update the counter that records which day number the data are from.  
            
        
        # start the figure.  
        f, ax = plt.subplots(1, figsize = figsize)                                                                  # one figure per day.  

        points = ax.scatter(licsalert_status_current[:, 0], licsalert_status_current[:, 1], 
                            c = day_counter, vmin = int(tbaselines[0]), vmax = int(tbaselines[-1]))
        
        # deal with the colorbar
        cbar = f.colorbar(points)
        original_tick_labels = cbar.get_ticks()
        new_tick_labels = []
        for tick_label in original_tick_labels:
            try:
                index = list(tbaselines).index(tick_label)                                      # these should be equivalent
                new_tick_labels.append(dt.strftime(day_list[index], '%Y_%m_%d'))                # 
            except:
                new_tick_labels.append('')
        cbar.set_ticklabels(new_tick_labels)
        cbar.set_label(f"Acquisition date")
        
        # label volcano names if above a threshold
        texts = []
        label_ratio = 0.2                                                                                                                           # the top ratio are labelled.  
        x_threshold = (1 - label_ratio) * np.nanmax(licsalert_status_current[:,0])
        y_threshold = (1 - label_ratio) * np.nanmax(licsalert_status_current[:,1])
        for volc_n, status in enumerate(licsalert_status_current):                                                                                  # loop through every volcano
            if (status[0] > x_threshold) or (status[1] > y_threshold):                                                                              # if one measure is above threshold    
                texts.append(ax.annotate(volc_names[volc_n], (licsalert_status_current[volc_n, 0], licsalert_status_current[volc_n, 1])))           # label the point.  
        adjust_text(texts, only_move={'points':'y', 'texts':'y'})                                                                                   # adjust how volcanoes of interest are labelled so that the labels don't overlap
                
        # tidy up the figure
        ax.set_xlabel('Sigma for exising def.')
        ax.set_ylabel('Sigma for new def.')
        ax.set_ylim([0,ylim])
        ax.set_xlim([0,xlim])
        title = dt.strftime(day, '%Y_%m_%d')
        ax.set_title(title)
        
        
        
        if len(volcs_off_chart) > 0:                                                        # if there are volcanoes off the chart here
            off_volcs = f"Volcanoes off chart:\n"
            for volcs_off_chart_sigma, volc_off_chart in zip(volcs_off_chart_sigmas, volcs_off_chart):
                off_volcs = off_volcs + f"{volc_off_chart} (x:{volcs_off_chart_sigma[0]:.2f} y:{volcs_off_chart_sigma[1]:.2f})\n"
            ax.annotate(off_volcs, (xlim - 0.02, ylim - 0.1), horizontalalignment='right',verticalalignment='top',)
        
        f.savefig(out_dir / f"{title}.png", bbox_inches='tight')            
        print(f"Saved figure {title}")

    if plt.get_backend() != 'Qt5Agg':                                                           # check if we need to reset the backend (to interactive window)
        plt.switch_backend('Qt5Agg')                                                           #  
        
        
    
    

licsalert_ts_all_volcs(licsalert_status, volc_names, day_list,  out_dir = Path('monitoring_steps_all_volcs'), figsize = (20, 12)    )


#%%
           

f, ax = plt.subplots(1)
xs = np.arange(0,100)
ys = np.log10(xs)
ax.plot(xs, ys)

        
#%%



    

    
    
    
    
   
    
    
    
    
#matrix_show(licsalert_status[:,:,0], title = 'change to def')    

#matrix_show(licsalert_status[:,:,1], title = 'new  def')    
    
    
    
    
    
    
    
    
    
    
    
    
    





