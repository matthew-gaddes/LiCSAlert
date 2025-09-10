#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:28:07 2023

@author: matthew

"""

print("Succesfully started the script.  ")

import numpy as np
import matplotlib.pyplot as plt 
from pathlib import Path
import pdb
import os
from glob import glob
import sys
from datetime import datetime
import pandas as pd

# needed to save user-defined objects
import cloudpickle 




#%% plot_unrest_metric_all_frames_one_volc()

def plot_unrest_metric_all_frames_one_volc(
        volc, ticks_labels_every_n = 6, scatter_points = False):
    """
    Plot the time series for new deformation and existing deformation 
    for each frame at a volcano, and the combined metric.  
    
    """
    import matplotlib.dates as mdates
    from matplotlib.lines import Line2D
    from datetime import datetime
    
    # remove no data region?
    
    def determine_figure_y_lims(volc):
        """
        """
        n_frames = len(volc.frame_status)
        frames = volc.frame_status

        # set initial bounds that shouldn't be exceeded for a while
        y_low = 0
        y_high = 1
        
        # iterate through and plot each frame.  
        for n_frame in range(n_frames):
    
            # get just the frame name with no volcano name.  
            frame_name = Path(frames[n_frame]).parts[-1][-17:]
            
            # if there's no status for a frame, it will return as NA
            if frame_name != 'NA':
                
                # make a 2 x times array of unrest values            
                unrest_metrics = np.array(
                    [volc.status[frame_name]['new_defs'],
                     volc.status[frame_name]['existing_defs']]
                         )
            
                # keep track of figure x limits
                if np.max(unrest_metrics) > y_high:
                    y_high = np.max(unrest_metrics)
                
        # add a nudge factor so points don't touch edge.  
        y_high +=1 
                
        return y_low, y_high
    
    def determine_figure_xlims(volc):
        """
        """
        from datetime import datetime
        
        n_frames = len(volc.frame_status)
        frames = volc.frame_status

        # set initial bounds that shouldn't be exceeded for a while
        fig_start_dt = datetime(9999, 1, 1, 00, 00, 00)
        fig_end_dt = datetime(1, 1, 1, 00, 00, 00)
        
        
        # iterate through and plot each frame.  
        for n_frame in range(n_frames):
    
            # get just the frame name with no volcano name.  
            frame_name = Path(frames[n_frame]).parts[-1][-17:]
            
            # if there's no status for a frame, it will return as NA
            if frame_name != 'NA':
                
                # keep track of figure x limits
                if volc.status[frame_name]['dates'][0] < fig_start_dt:
                    fig_start_dt = volc.status[frame_name]['dates'][0]
                if volc.status[frame_name]['dates'][-1] > fig_end_dt:
                    fig_end_dt = volc.status[frame_name]['dates'][-1]
                
        return fig_start_dt, fig_end_dt
    


    def plot_unrest_metric_one_frame(
            ax, status, d_start_dt, d_end_dt, y_low, y_high, y_title,
            show_xtick_labels = False, 
            ):
        """
        """
        
        dates = status['dates']
        new_def = status['new_defs']
        existing_def = status['existing_defs']
        
        
        # Plot left data on primary y-axis
        line1, = ax.plot(dates, new_def, color=left_colour)
        if scatter_points:
             ax.scatter(dates, new_def, color=left_colour, marker = '.')
             
        if show_xtick_labels:
            ax.set_xlabel('Date')
        #ax.set_ylabel("New deformation", color=left_colour)
        ax.set_ylabel(y_title, color='k', fontsize = 8)
        ax.tick_params(axis='y', labelcolor=left_colour)
        ax.set_ylim(y_low, y_high)
        
        # Create a second y-axis sharing the same x-axis
        ax2 = ax.twinx()
        line2, = ax2.plot(dates, existing_def, color=right_colour)
        if scatter_points:
            ax2.scatter(dates, existing_def, color=right_colour, marker = '.')
        #ax2.set_ylabel('Existing Deformation', color=right_colour)
        ax2.tick_params(axis='y', labelcolor=right_colour)
        ax2.set_ylim(y_low, y_high)
               
        ax.set_xlim(left = d_start_dt, right = d_end_dt)
        # Set up the major locator: every ticks_labels_every_n months
        major_locator = mdates.MonthLocator(interval=ticks_labels_every_n)
        ax.xaxis.set_major_locator(major_locator)
        
        if show_xtick_labels:
            # Use a formatter to display dates in YYYY_MM_DD format
            major_formatter = mdates.DateFormatter("%Y_%m_%d")
        else:
            # Turn off tick labels using a NullFormatter
            from matplotlib.ticker import NullFormatter
            major_formatter = NullFormatter()
        ax.xaxis.set_major_formatter(major_formatter)
        
        if show_xtick_labels:
            # Rotate the date labels 45 degrees for better readability
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)

        
        # Set up the minor locator: every month (no labels by default)
        minor_locator = mdates.MonthLocator(interval=1)
        ax.xaxis.set_minor_locator(minor_locator)
        
    
    x_low_dt, x_high_dt = determine_figure_xlims(volc)
    y_low, y_high = determine_figure_y_lims(volc)

    left_colour = 'tab:blue'
    right_colour = 'tab:orange'
    
    # determine the number of frames for this volcano
    n_frames = len(volc.frame_status)
    frames = volc.frame_status
    

    
    # dates = volc.status_combined['dates']
    # new_def = volc.status_combined['new_defs']
    # existing_def = volc.status_combined['existing_defs']
    
    # Create a new figure and primary axis
    fig, axes = plt.subplots(n_frames+1, 1, figsize=(10, 5))
    
    # iterate through and plot each frame.  
    for n_frame in range(n_frames):

        # get just the frame name with no volcano name.  
        frame_name = Path(frames[n_frame]).parts[-1][-17:]
   
        # print(
        #     f"{volc.status[frame_name]['dates'][0]} - " 
        #     f"{volc.status[frame_name]['dates'][-1]}" 
        #     )
        if frame_name != 'NA':
            plot_unrest_metric_one_frame(
                axes[n_frame], volc.status[frame_name], x_low_dt, x_high_dt,
                y_low, y_high, frame_name
                )
        
    # bottom row is combined status
    plot_unrest_metric_one_frame(
        axes[-1], volc.status_combined, x_low_dt, x_high_dt, 
        y_low, y_high, 'Combined', show_xtick_labels = True
        )
        
    # add the legend
    # Create two dummy handles that don't draw anything.
    dummy_handles = [
        Line2D([], [], linestyle='None'),  # blank handle for first word
        Line2D([], [], linestyle='None')   # blank handle for second word
    ]

    # Create the legend with no visible handles.
    legend = fig.legend(
        dummy_handles, ['New deformation', 'Existing deformation'],
        handlelength=0,  # no line or marker shown
        frameon=True, loc = 'lower left')   

    # Set the colors of the legend texts.
    for text, color in zip(legend.get_texts(), [left_colour, right_colour]):
        text.set_color(color)
        
    # -999 puts it behind tick labels
    legend.set_zorder(0)

    plt.title(f"{volc.name}")
    fig.tight_layout()  # Adjust layout to fit labels and title
    plt.show()

#%% plot_unrest_metric_all_volcs()


def plot_unrest_metric_all_volcs(volcs):
    """
    Plots unrest metrics for all volcanoes. When hovering over a line,
    an annotation box near the cursor displays the volcano's name and the line is highlighted.
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    lines = []
    # Dictionary to store each line's original linewidth.
    default_linewidth = {}
    
    # Plot each volcano's data and store the line objects.
    for volc in volcs:
        line, = ax.plot(
            volc.status_overall['dates'], volc.status_overall['unrest_metric'],
            label=volc.name
        )
        lines.append(line)
        default_linewidth[line] = line.get_linewidth()  # Save original linewidth.
    
    # Create the legend.
    legend = fig.legend(frameon=True, loc='lower left')
    
    # Create an annotation that will act as the tooltip.
    annot = ax.annotate(
        "", xy=(0, 0), xytext=(20, 20),
        textcoords="offset points",
        bbox=dict(boxstyle="round", fc="w"),
        arrowprops=dict(arrowstyle="->")
    )
    annot.set_visible(False)
    
    def on_move(event):
        # Only process events within the axes.
        if event.inaxes != ax:
            annot.set_visible(False)
            # Reset all lines to their default linewidth.
            for line in lines:
                line.set_linewidth(default_linewidth[line])
            fig.canvas.draw_idle()
            return
        
        # Flag to indicate whether any line is hovered.
        found = False
        for line in lines:
            contains, _ = line.contains(event)
            if contains:
                # Update and show the annotation.
                annot.xy = (event.xdata, event.ydata)
                annot.set_text(line.get_label())
                annot.set_visible(True)
                
                # Highlight the hovered line (increase linewidth).
                for l in lines:
                    if l == line:
                        l.set_linewidth(default_linewidth[l] + 2)
                    else:
                        l.set_linewidth(default_linewidth[l])
                
                fig.canvas.draw_idle()
                found = True
                break
        if not found:
            # Hide the annotation and reset line styles.
            if annot.get_visible():
                annot.set_visible(False)
                for l in lines:
                    l.set_linewidth(default_linewidth[l])
                fig.canvas.draw_idle()
    
    # Connect the mouse motion event to the handler.
    fig.canvas.mpl_connect("motion_notify_event", on_move)
    
    plt.show()



#%% order_volcs_by_cumulative_unrest()



def order_volcs_by_cumulative_unrest(volcs):
    """ Each volcano in volc has a cumulative unrest metric (the integral
    of the unrest metric).  This re-orders the volcs list so that the volcano
    with the highest cumulaive unrest metric comes first.  
    """
    cumulative_unrests = []
    
    for volc in volcs:
        cumulative_unrests.append(
            volc.status_cumulative['unrest_cum'][-1]
            )
        
    # Zip indices with data.
    combined = list(enumerate(cumulative_unrests))
    # Sort by the data value.
    sorted_combined = sorted(combined, key=lambda x: x[1], reverse=True)
    # Unzip the sorted pairs.
    sorted_indices = [index for index, _ in sorted_combined]
    sorted_data = [item for _, item in sorted_combined]
    
    volcs_sorted = [volcs[i] for i in sorted_indices]
    
    return volcs_sorted


#%% find_volc_indexes()

def find_volc_indexes(volcs, volc_names):
    """ Given a list of volcs, find which number a volcano is when specified
    by name
    """
    import difflib
    
    volc_indexes = []
    volc_names_found = []
    
    # extract the volcano names
    volcs_extracted_names = [volc.name for volc in volcs]
    
    for volc_name in volc_names:
    
        matches = difflib.get_close_matches(
            volc_name, volcs_extracted_names, n=1
            )
        
        if len(matches) > 0:
            volc_indexes.append(volcs_extracted_names.index(matches[0]))
            volc_names_found.append(volc_name)
        else:
            print(f"Unable to find a volcano for {volc_name}")
            
    return volc_indexes, volc_names_found
            



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




#%%
# def find_directories_with_more_than_n_subdirectories(root_dir, n=2):
#     # List to store directories with more than n subdirectories
#     result_directories = []

#     # Walk through the directory tree
#     for root, dirs, files in os.walk(root_dir):
#         # Count the number of subdirectories in the current directory
#         num_subdirectories = len(dirs)
        
#         # If the number of subdirectories exceeds n, add the directory to the result list
#         if num_subdirectories > n:
#             result_directories.append(root)

#     return result_directories

# # Example usage:
# directory_to_search = '/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/05_jasmin_clones/02_2024_01_19_cov_only/'

# # Call the function
# directories_with_more_than_2_subdirectories = find_directories_with_more_than_n_subdirectories(directory_to_search, n=2)

# directories_with_more_than_2_subdirectories = sorted(directories_with_more_than_2_subdirectories)

# # Write the result to a text file
# with open('directories_with_more_than_2_subdirectories.txt', 'w') as f:
#     for directory in directories_with_more_than_2_subdirectories:
#         f.write(directory + '\n')
    


#%% Imports

from licsalert.jasmin_tools import open_comet_frame_files
from licsalert.jasmin_tools import write_jasmin_download_shell_script
from licsalert.jasmin_tools import update_volcs_with_data
from licsalert.jasmin_tools import volc_names_to_licsalert_volcs
from licsalert.jasmin_tools import comet_db_to_licsalert_volcs
#from licsalert.jasmin_tools import get_lon_lat_of_volcs_from_ts_data
from licsalert.jasmin_tools import get_lon_lat_of_volcs_from_db


from licsalert.status import extract_licsalert_status, get_all_volcano_dirs
from licsalert.status import find_volcano, get_volc_names_fron_dir_of_frames


from licsalert.temporal import create_day_list
from licsalert.volcano_portal_tools import get_portal_public_volcanoes
from licsalert.plotting import offset_volc_lls


from licsalert.plotting import status_fig_one_volc, status_fig_all_volcs
        
#%% Things to set

#%%  ################# jasmin outputs 2023 (for fringe)
# #licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/01_test/")
# licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/06_LiCSAlert/09_jasmin_clone_2023_08_30/02_all/")
# out_dir = Path("./monitoring_jasmin/one_volc/ale_bagu")
# d_start = "20200101"
# d_stop = "20230901"
# test_volcano = 'ale_bagu*'
# regions = True
################## end jasmin outputs 2023 (for fringe)

#%% ################## jasmin outputs 2024 (for Cities on Volcanoes)

# option 1: volcanoes on public portal  
# cov_volc_names = ['Antisana',
#                  'Cerro Overo',
#                  'Cordon del Azufre',
#                  'Cotopaxi',
#                  'Domuyo',
#                  'Fernandina',
#                  #'Guagua Pichinca',                    # volcano portal spelling
#                  'Guagua Pichincha',                     # licsbas spelling
#                  'Laguna del Maule',
#                  'Masaya',
#                  'Nevados de Chillan',
#                  'Planchón-Peteroa',                    # 
#                  'Sabancaya',
#                  'San Miguel',
#                  'San Salvador',
#                  'Sangay',
#                  'Santa Ana',
#                  'Sierra Negra',
#                  'Tungurahua',
#                  'Turrialba',
#                  'Villarrica']

# # option 2: Camila deforming volcanoes list
# cov_volc_names = [ # camilla deforming list
#                   "Cordon Caulle",
#                   "Wolf",
#                   "Sabancaya",
#                   "Laguna del Maule",
#                   "Cordon del Azufre",
#                   "Sangay",
#                   "Domuyo",
#                   "Alcedo",
#                   "Mevados de Chillan",
#                   "Cerro Azul",
#                   "Fernandina",
#                   "Sierra Negra",
#                   # below are Visible on portal, def observed, central america
#                   "Masaya",                 
#                   "San Miguel",
#                   "Santa Ana",
#                   "Turrialba"]



# # ?
# # licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/"
# #                      "06_LiCSAlert/05_jasmin_clones/02_2024_01_19_cov_only")

# # current clone of Jasmin data
# licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/"
#                       "06_LiCSAlert/05_jasmin_clones/2024_06_11_comet_talk")

# # print(f"USING THE TEMPORARY TEST DIRECTORY")
# # licsalert_dir = Path("/home/matthew/university_work/03_automatic_detection_algorithm/"
# #                       "06_LiCSAlert/05_jasmin_clones/test")


# out_dir = Path("./status_outputs/03_comet")
# # d_start = "20200101"                # LiCSAlert baseline ends at end of 2019?  
# d_start = "20180101"                # LiCSAlert baseline ends at end of 2019?  
# d_stop = "20240621"
# test_volcano = 'sabancaya*'
# regions = True

# # volcs_omit = ["guagua_pichincha_121A_08908_131313",  # no coherence
# #               "planchon-peteroa_018A_12472_131313",  # single pixels of large value
# #               ]
# volcs_omit = []

################## end jasmin outputs 2024 (for Cities on Volcanoes)
#%% ################## jasmin outputs 2024 Chile 

# # current clone of Jasmin data
# licsalert_dir = Path(
#     "/home/matthew/university_work/03_automatic_detection_algorithm/"
#     "06_LiCSAlert/05_jasmin_clones/2024_11_14_chile_v2"
#     )

# # names of all comet frames.  
# comet_volcano_frame_index_dir = Path('./comet_volcano_frames/')

# out_dir = Path("./status_outputs/04_chile_november_2024")
# d_start = "20210101"                
# d_stop = "20241101"
# #test_volcano = 'sabancaya*'
# regions = False
# volcs_omit = []
# generate_bash_dl_file = True
# ################## end jasmin outputs 2024 for Chile

#%% ################## 10 years of Sentinel paper (2025_02_12)


# # current clone of Jasmin data - bug with baseline data
# # licsalert_dir = Path(
# #     "/home/matthew/university_work/03_automatic_detection_algorithm/"
# #     "06_LiCSAlert/05_jasmin_clones/2025_02_13_s1_10yr_licsalert/"
# #     "2025_02_13_licsalert/licsalert_sync"
# #     )

# # #v2, 
# licsalert_dir = Path(
#     "/home/matthew/university_work/03_automatic_detection_algorithm/"
#     "06_LiCSAlert/05_jasmin_clones/2025_02_28_s1_10yr/licsalert_sync"
#     )

# # # #v3, only status files, but for all volcs
# # licsalert_dir = Path(
# #     "/home/matthew/university_work/03_automatic_detection_algorithm/"
# #     "06_LiCSAlert/05_jasmin_clones/2025_03_06_s1_10yr_status_only/licsalert"
# #     )

# # volc_names = [  
# #     "Campi Flegrei",
# #     "Fernandina",
# #     "Darwin",
# #     "Sierra Negra",
# #     "Cerro Azul",
# #     "Kilauea",
# #     "Wolf",
# #     "Laguna del Maule",
# #     "Ecuador",
# #     "Reykjanes",
# #     "Alcedo",
# #     "Domuyo",
# #     "Krysuvik",
# #     "Fogo",
# #     "Fujisan",
# #     # bottom left of figure
# #     "Mauna Loa",
# #     "Mauna Kea",
# #     "Etna",
# #     "Askja",
# #     "Iztaccihuatl",
# #     "Pico",
# #     "Bardarbunga",
# #     "Esjufjoll",
# #     "Kverkfjoll"
# # ]


# # names of all comet frames.  
# comet_volcano_frame_index_dir = Path('./comet_volcano_frames/')
# comet_volcano_frame_index_dir = Path('./comet_volcano_frames_fogo_edit/')

# #out_dir = Path("./status_outputs/04_chile_november_2024")
# d_start = "20140101"                
# d_stop = "20250219"
# #test_volcano = 'sabancaya*'
# regions = True
# volcs_omit = []
# generate_bash_dl_file = True
# ################## end jasmin outputs 2024 for Chile

# combined_status_method = 'previous'
# combined_status_method = 'window'


#%% COMET summer 25


# ## Generate the volc_names list
# volc_priorities_used = ["A1"]
# lists = comet_db_to_licsalert_volcs(
#     "./comet_volcano_database.pkl", 
#     volc_priorities_used
#     )

# volc_names_a1 = []
# for sublist in lists:
#     volc_names_a1.extend(sublist)
# del lists



# # #v3, only status files, but for all volcs
# print(f"USING THE TEST DIR OF THE DATA CLONED FROM JASMIN")
# licsalert_dir = Path(
#     "/home/matthew/university_work/03_automatic_detection_algorithm/"
#     "06_LiCSAlert/05_jasmin_clones/2025_05_22_comet_2025/licsalert_sync_test/"
#     )


# # names of all comet frames.  
# #comet_volcano_frame_index_dir = Path('./comet_volcano_frames/')
# comet_volcano_frame_index_dir = Path('./comet_volcano_frames_fogo_edit/')


# d_start = "20140101"                
# d_stop = "20250219"
# #test_volcano = 'sabancaya*'
# regions = True
# volcs_omit = []
# generate_bash_dl_file = True


# #combined_status_method = 'previous'
# combined_status_method = 'window'               # old method uses 'previous'





# #%%  #################  local test volcs (processed with LiCSBAS by me?  )
# # all_volcs_figs_dir = Path('monitoring_steps_all_volcs_test')
# # licsalert_dir = Path("/home/matthew/university_work/31_from_bright/2023_09_08/01a_LiCSAlert_batch_mode_2/")
# # regions = False

# # # Option for the 1 volc we plot in detail.  
# # # option 1
# # out_dir = Path("./monitoring_test/one_volc/sierra_negra")
# # d_start = "20170101"
# # d_stop = "20230901"
# # test_volcano = '*sierra*'

# # # option 2
# # # out_dir = Path("./monitoring_test/one_volc/erta_ale")
# # # d_start = "20160601"
# # # d_stop = "20230901"
# # # test_volcano = '*erta*'

# ################# end local test volcs            

# #%% Step 00: volcano names:
   
# # volc_names might be passed as a list of strings.  
# if "volc_names" not in locals():
        
#     # if it's not, try and generate it.  
#     try:
#         volc_names_dir = get_volc_names_fron_dir_of_frames(
#             licsalert_dir, regions = regions
#             )
        
#         print(
#             "Succesfully created a list of volcanoes (volc_names) from a "
#             "licsalert directory.  Continuing.  "
#             )
#     except:
#         print(
#             "Failed to get the volcano names from a licsalert directory.  Perhaps "
#             "this is the first time this has been run?  You could manually "
#             "provide a list of volc_names, but for now we'll try to generate "
#             "one from the volcanoes that are visible on the COMET volcano portal"
#             )
    
#         try:
#             # get the names of the COMET portal public volcanoes
#             volc_names = get_portal_public_volcanoes()
#             print(
#                 "Succesfully generated a list of volcanoes that are visible "
#                 "publicly on the COMET volcano portal.  Continuing.  "
#                 )
#         except:
#             print(
#                 "Failed to generate a list of volcanoes from the COMET volcano "
#                 "portal.  Exiting.  ")




# # choose which volc_names to use.  
# volc_names = volc_names_a1
# # volc_names = volc_names_dir
# # print("SELECTING ONLY A SUBSET OF THE VOLC_NAMES")
# # volc_names = sorted(volc_names)[165:168]

# # open the info on COMET volcano frames and what region they're in.  
# # region is a key, and each value is a list of comet frames in that region
# comet_volcano_frame_index = open_comet_frame_files(comet_volcano_frame_index_dir)

# # convert volc_names to volcs, which is a list of comet_volcano objects
# # output is jasmin_sync_script.sh
# volcs = volcano_name_to_comet_frames(volc_names, comet_volcano_frame_index)    

# # tidy up.  
# del volc_names


# # for i, v in enumerate(volcs):
# #     print(f"{i} {v.name}")


# #%% Step 01: Generate bash file to copy selected data from Jasmin



# # if generate_bash_dl_file:
# #     jasmin_dir = Path(
# #         "mgaddes@xfer-vm-03.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/"
# #         "projects/LiCS/volc-portal/processing_output/licsalert"
# #         )
# #     local_dir = Path('./licsalert_sync/')
# #     bash_sync_fle = "jasmin_sync_script.sh"
    

# #     # write a shell script to download the LiCSAlert data from jasmin 
# #     write_jasmin_download_shell_script(
# #         jasmin_dir, local_dir, bash_sync_fle, volcs, exclude_json_gz = True,
# #         exclude_original_ts = False, exclude_ICASAR_data = True,
# #         exclude_epoch_images = True, exclude_epoch_data = True)






# #%% Step 02:
    
# # # run jasmin sync file on FOE-linux (to copy data from jasmin to school server)
# # # copy from FOE-linux to local.  

# #%% Step 02a: possibly move directories that don't have LiCSAlert results
# # # to a different directory

# # from licsalert.jasmin_tools import move_small_directories

# # move_small_directories(
# #     licsalert_dir, licsalert_dir.parent / "to_delete", n = 6)

# # print("\n\n\nSELECTING ONLY A SUBSET OF THE DATA\n\n\n")
# # volcs = volcs[0:5]


# #%% Step 03: Compile licalert status for all frames at all times.  

# # get path to licsalert frames for each volcano, and their names (snake case)
# volc_frames, volc_frame_names = get_all_volcano_dirs(
#     licsalert_dir, volcs_omit,regions
#     )




# # update volcs so it only contains volcanoes we have data for.  
# # also  convert to a list with an entry for each volcano, with frame
# # info accessed from the objects attributes.  
# volcs = update_volcs_with_data(
#     volcs, volc_frames, volc_frame_names,verbose=False, remove_all_nans=False
#     )


# #  add the lon lat info for each volcano
# # Method 1: from data 
# # volcs = get_lon_lat_of_volcs_from_ts_data(volcs)
# # method 2: from comet database
# get_lon_lat_of_volcs_from_db(
#     volcs,
#     pd.read_pickle("comet_volcano_database.pkl")
#     )


# # adjust the lon lats of the volcanoes so when plotted they don't overlap
# offset_volc_lls(volcs,  threshold = 2., offset = 2.1, attempts = 200)

# # simple list of datetime for each status day
# day_list = create_day_list(d_start, d_stop)


# # get the licsalert status for each volcano frame (update volcs)
# extract_licsalert_status(volcs, day_list, combined_status_method)           

# # tidy up
# del volc_frames, volc_frame_names


# #%%


# #%% Step 04: figure for one volcano at all times.  
# # (the one with the LiCSAlert inset bit?)

# # name_and_indexes = find_volcano(volc_frame_names, test_volcano)



# # out_dir_1_volc = out_dir / name_and_indexes[0][0]

# # status_fig_one_volc(name_and_indexes[-1], licsalert_status, day_list, 
# #                 fig_type = 'png', out_dir = out_dir_1_volc)

# # status_fig_one_volc(name_and_indexes[-1], licsalert_status, day_list, 
# #                     fig_type = 'gif', volc_dirs = volc_frame_dirs, 
# #                     out_dir = out_dir_1_volc/ "gif") 
# # #pngs_to_gif(out_dir / "gif", out_dir  / "gif" / "animation.gif", image_duration = 0.5)






# #%% Step 05: 2D plot figure for all volcs at all times.  


# # status_fig_all_volcs(licsalert_status, volc_frame_names, day_list,  
# #                      out_dir = out_dir / "all_volcs", 
# #                      figsize = (20  , 12), plot_frequency = 6, label_fs = 16)

# # # # convert the directory of pngs to a single gif
# # # pngs_to_gif(out_dir / "all_volcs", out_dir / "all_volcs" / "animation.gif",
# # #             image_duration = 250)                                                 



# #%% Step 06: New worldmap figure


# """
# This needs editing to show volcanoes in grey if they have no licsalert
# status, and then colour by status.  

# """

# def licsalert_status_map(
#         volcs, outdir,  sigma_min = 0., sigma_max = 10., 
#         d_start = None, d_end = None,
#         plot_frequency = "monthly",  figure_type = 'window',
#         no_licsbas_c='tab:red',
#         no_licsalert_c='tab:pink'
#         ):
#     """ Plot multiple volcano statuses on the worldmap.  Colour indicates
#     status.  
    
#     Inputs:
#         volcs | list of comet_volcanos | contains info such as lon_lat, 
#                                          lon_lat_offset (as above, but shifted
#                                          so points don't lie on top of each other)
#          outdir | Path | ouput png files, if backed is 'agg'
#          d_start | str | start date for period to plot daily / monthly yearly
#          d_end | str | end date for period to plot daily / monthly yearly

#          sigma_max | int or float | maximum number of sigmas for colourscale 
#                                  i.e. if set to 10, any signal that is 10 sigmas
#                                      or more will plot as maximum colour (yellow)
#          plot_frequency | string | daily / monthly / yearly 
#          backend | string | 'agg' if exporting pngs, 'qt5agg' for interactive. 
#      Returns:
#          Figure
         
#      History:
#          2024_06_21 | MEG | Written.  
         
#     """
    
#     import matplotlib.pyplot as plt
#     import matplotlib.colors as colors
#     import matplotlib.cm as cm
#     import cartopy.crs as ccrs
#     import datetime as dt
#     import math
#     import warnings
#     from copy import deepcopy
    

#     def plot_licsalert_result(volc, dayn_index):
#         """ 
#         """
        
#         existing_def = volc.status_combined['existing_defs'][dayn_index]
#         new_def = volc.status_combined['new_defs'][dayn_index]
        
        
#         # determine the maximum from the two unrest metrics
#         combined_def = max(existing_def, new_def)
        
#         # only plot if larger than minimum threshold
#         if combined_def > sigma_min:
        
#             # if we do, plot the max of existing or new deformation
#             sc = ax.scatter(
#                 volc.lon_lat_offset[0], volc.lon_lat_offset[1], 
#                 c=combined_def, s=50, transform=ccrs.PlateCarree(),
#                 vmin = sigma_min, vmax = sigma_max
#                 )
#             scatter_objects.append((sc, volc.name))
            
#             # also plot the lines from the shifted points to their true point.  
#             ax.plot(
#                 [volc.lon_lat_offset[0], volc.lon_lat[0]],
#                 [volc.lon_lat_offset[1], volc.lon_lat[1]],
#                 c = 'k',  transform=ccrs.PlateCarree()
#                 )
            
#     def plot_volc_manual_colour(volc, c):
#         """
#         """
        
#         sc = ax.scatter(
#             volc.lon_lat_offset[0], volc.lon_lat_offset[1], 
#             c=c, s=50, transform=ccrs.PlateCarree(),
#             )
        
#         # also plot the lines from the shifted points to their true point.  
#         ax.plot(
#             [volc.lon_lat_offset[0], volc.lon_lat[0]],
#             [volc.lon_lat_offset[1], volc.lon_lat[1]],
#             c = 'k',  transform=ccrs.PlateCarree()
#             )
    

#     def add_status_legend(
#             ax,
#             *,
#             colour_nolicsbas,
#             colour_nolicsalert,
#             bbox=(0.02, 0.05, 0.3, 0.08)   # (x0, y0, width, height) in axes fraction
#         ):
#         from matplotlib.lines import Line2D
        
#         """
#         Add a two-row legend with coloured dots.
    
#         Parameters
#         ----------
#         ax : matplotlib Axes
#             The map axes you plotted on.
#         colour_nolicsbas : str
#             Dot colour for “No LiCSBAS ts.”    (default 'tab:red')
#         colour_nolicsalert : str
#             Dot colour for “No LiCSAlert result” (default 'tab:pink')
#         bbox : tuple
#             (x0, y0, width, height) in **axes-fraction units** (0-1).
#             Move/resize this to taste.
    
#         Returns
#         -------
#         legend : matplotlib Legend
#             The legend handle (in case you want to tweak later).
#         """
#         handles = [
#             Line2D([0], [0], marker='o', color='w',
#                    markerfacecolor=colour_nolicsbas, markeredgecolor='k',
#                    markersize=10, label="No LiCSBAS ts."),
#             Line2D([0], [0], marker='o', color='w',
#                    markerfacecolor=colour_nolicsalert, markeredgecolor='k',
#                    markersize=10, label="No LiCSAlert result"),
#         ]
    
#         legend = ax.legend(
#             handles=handles,
#             loc="lower left",
#             bbox_to_anchor=bbox,
#             borderpad=0.8,
#             frameon=True,
#             fancybox=True
#         )
#         return legend


#     # if plt.get_backend() != backend:   
#     #     print(f"The matplotlib backend was previously {plt.get_backend()}, "
#     #           f"but has been switched to ", end = '')
#     #     plt.switch_backend(backend)
#     #     print(f"{plt.get_backend()}")
    
#     # Check matplotlib backend is set correctly 
#     if figure_type == 'png':
#         plt.switch_backend('Agg')                                                           
#     else: 
#         if plt.get_backend() != 'Qt5Agg':                                                               
#             plt.switch_backend('Qt5Agg')                                                           

#     # 0: small steps to handle dates:    
#     # all volcanoes should shre the same day list (list of datetimes, daily),
#     # but it won't exist for some volcs that don't have a licsalert result
#     day_list = next(
#         (v.status_combined['dates'] for v in volcs 
#          if hasattr(v, "status_combined")), None
#         )
    
#     d_start_dt = dt.datetime.strptime(d_start, "%Y%m%d")
#     d_end_dt = dt.datetime.strptime(d_end, "%Y%m%d")


#     # 1 figure for each day.
#     for day_n, day in enumerate(day_list):

        
#         if (day > d_start_dt) and (day < d_end_dt):
#             # check if we should plot this one
#             if plot_frequency == "daily":
#                 plot_today = True
#             elif plot_frequency == 'monthly':
#                 plot_today = day.day == 1
#             elif plot_frequency == 'yearly':
#                 plot_today = (day.day == 1) and (day.month == 1)
#         else:
#             plot_today = False
        
#         if plot_today:
            
#             # get figure date as string
#             day_str = dt.datetime.strftime(day, "%Y%m%d")
#             print(
#                 f"Creating a LiCSAlert status figure for {day_str}.  "
#                 )

#             fig = plt.figure(figsize=(20, 11))
#             ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
            
#             # Make the map global
#             ax.set_global()
#             ax.stock_img()
#             ax.coastlines()
        
#             scatter_objects = []
#             # loop through each volcano to plot it as a point on the map.
#             for volc_n, volc in enumerate(volcs):
                
#                 # three options here.  Either no LiCSBAS ts, no licsalert result
#                 # or a licsalert result
        
#                 if volc.processing_status.licsbas_ts == False:
#                     plot_volc_manual_colour(
#                         volc, no_licsbas_c
#                         )
#                 if volc.processing_status.licsalert_result == False:    
#                     plot_volc_manual_colour(
#                         volc, no_licsalert_c
#                         )
#                 else:
#                     dayn_index = volc.status_combined['dates'].index(day)
#                     plot_licsalert_result(volc, dayn_index)
        
                
            
#             #colorbar (possibly with no sc created by ax.scatter)
#             cbar_ax = fig.add_axes([0.35, 0.15, 0.3, 0.01])  
#             norm = colors.Normalize(vmin=sigma_min, vmax=sigma_max)
#             cmap = cm.viridis  
#             sm = cm.ScalarMappable(norm=norm, cmap=cmap)
#             cbar = plt.colorbar(sm, cax=cbar_ax, orientation = 'horizontal')
#             cbar.set_label("# sigma from background")
            
#             # title, first make the date a string.  
#             date_str = dt.datetime.strftime(
#                 next(
#                     (v.status_combined['dates'][dayn_index] for v in volcs 
#                      if hasattr(v, "status_combined")), None
#                     ),
#                 "%Y/%m/%d")
#             fig.suptitle(date_str)

    
#             ## save as a png and close the figure.                     
#             # Suppress warnings for tight_layout
#             with warnings.catch_warnings():
#                 warnings.simplefilter("ignore")
#                 plt.tight_layout()
            
#             # start the interactive part of the figure
#             if (figure_type == 'window') or (figure_type == 'both'):
#                 # Create an annotation object
#                 annot = ax.annotate("", xy=(0,0), xytext=(10,10),
#                                     textcoords="offset points",
#                                     bbox=dict(boxstyle="round", fc="w"),
#                                     arrowprops=dict(arrowstyle="->"),
#                                     transform=ccrs.PlateCarree())
#                 annot.set_visible(False)
                
#                 # Function to update the annotation
#                 def update_annot(ind, scatter_obj, label):
#                     pos = scatter_obj.get_offsets()[ind["ind"][0]]
#                     annot.xy = pos
#                     annot.set_text(label)
#                     annot.get_bbox_patch().set_alpha(0.4)
                
#                 # Function to check if mouse is over a scatter point
#                 def hover(event):
#                     vis = annot.get_visible()
#                     if event.inaxes == ax:
#                         for scatter_obj, label in scatter_objects:
#                             cont, ind = scatter_obj.contains(event)
#                             if cont:
#                                 update_annot(ind, scatter_obj, label)
#                                 annot.set_visible(True)
#                                 fig.canvas.draw_idle()
#                                 return
#                     if vis:
#                         annot.set_visible(False)
#                         fig.canvas.draw_idle()
                
#                 # Connect the hover event to the hover function
#                 fig.canvas.mpl_connect("motion_notify_event", hover)

#             # 8: Possible save output
#             if (figure_type == 'png') or (figure_type == 'both'):
#                 fig.savefig(outdir / f"{day_str}.png")
#                 print(f"Saved the LiCSAlert status map for {str}")
            
#         else:
#             pass
        
    
    

#     add_status_legend(
#         ax,
#         colour_nolicsbas=no_licsbas_c,       # must match your plot_volc_manual_colour
#         colour_nolicsalert=no_licsalert_c,
#         bbox=(0.01, 0.01, 0.35, 0.1)      # tweak position/size as you like
#     )

    
# licsalert_status_map(
#     volcs, sigma_min = 0., sigma_max = 10.,
#     #d_start = '20151231', d_end = '20170101',
#     d_start = '20241231', d_end = '20251231',
#     outdir = Path("./status_outputs/03_comet/"),
#     plot_frequency = 'yearly', figure_type = 'window',
#     no_licsbas_c='tab:red',
#     no_licsalert_c='tab:pink'
#     )

# # possibly create an animation from multiple frames
# # pngs_to_gif(out_dir / "short_status_map", out_dir / "status_map" / "status_map_animation.gif",
# #             image_duration = 750)                                                 
    
    


#%% Load / save before plotting




# print(f"\n\n Saving a pickle of 'volcs' and 'day_list' \n\n ")
# with open("monitoring_analysis_locals_RSE_volcs.pkl", 'wb') as f:
#     cloudpickle.dump(volcs, f)
#     cloudpickle.dump(day_list, f)




print(f"\n\n OPENING A PICKLE OF 'volcs' AND 'day_list' \n\n ")
# with open("monitoring_analysis_locals.pkl", 'rb') as f:    
with open("monitoring_analysis_locals_RSE_volcs.pkl", 'rb') as f:    
    volcs = cloudpickle.load(f)
    day_list = cloudpickle.load(f)



#%%  Plot a figure for each volcano showing the status for each frame at all times

fig_indexes, _ = find_volc_indexes(volcs, ['reykjanes'])
fig_indexes, _ = find_volc_indexes(volcs, ['wolf'])
fig_indexes, _ = find_volc_indexes(volcs, ['la palma'])             # not in data?
fig_indexes, _ = find_volc_indexes(volcs, ['laguna del maule'])             # 

for volc in volcs[fig_indexes[0]:fig_indexes[0]+1]:
    print(f"Plotting the unrest metrics for all frames at {volc.name}")
    plot_unrest_metric_all_frames_one_volc(
        volc, ticks_labels_every_n = 6, scatter_points=True
        )


"""
# Campi Flegrei?
# Sierra Negra?
# La Palma?
Can we make one subplot for with eruption and one without?

"""

#%% Order by cumulative unrest

volcs = order_volcs_by_cumulative_unrest(volcs)


# interactive figure 
#plot_unrest_metric_all_volcs(volcs[8:20])

volcs_to_remove = [
    'chagulak',
    ]

# v1
# user defined as interesting
# volcs_to_plot = [
#     'reykjanes',
#     'mauna_loa',
#     'wolf',
#     'alcedo',
#     'ecuador',
#     'darwin',
#     'askja',
#     'domuyo',
#     'campi_flegrei'
#     'kilauea',
#     'mauna_kea'
#     ]

# v2 3 with eruption, 3 without
volcs_to_plot = [
    'reykjanes',
    'mauna_loa',
    'wolf',
    'Askja',
    'Campi Flegrei',
    'Sierra Negra',
    'Fogo',
    'Fujisan'
    ]


fig_indexes, _ = find_volc_indexes(volcs, volcs_to_plot)

    

#%%


def unrest_metric_all_volcs_publication_figure(
        volcs, n_rows, n_per_row, x_lims, baseline_end,
        labels_dict1 = None, labels_dict2 = None, 
        y_lims=None, arrow_length = 0.1, text_angle=33, figsize = (7.84, 7.84)):
    """
    Plots unrest metrics for volcanoes on multiple subplots. Each subplot displays a subset of volcanoes.
    Additionally, annotations are added outside (above) the first and second subplots based on the provided
    dictionaries of labels.
    
    Parameters:
      volcs: list of volcano objects (each must have attributes:
             - status_overall: dict with keys 'dates' (list of datetime objects) and 'unrest_metric' (list)
             - name: volcano name)
      n_rows: number of subplot rows (expected to be 2 for this version)
      n_per_row: number of volcanoes per row
      x_lims: tuple or list with left and right limits for x-axis (e.g. [datetime(2021,1,1), datetime(2025,1,1)])
      labels_dict1: dictionary for annotations on the top of the first subplot.
                    Keys are label texts; values are date strings in the format "YYYY_MM_DD"
      labels_dict2: dictionary for annotations on the top of the second subplot.
                    Same format as labels_dict1.
      y_lims: optional, tuple or list with bottom and top y-limits.
    """
    
    
    from datetime import datetime
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    
    def color_extender(colors, n_volcs):
        """ Extend a list of colours so that it's longer than n_volcs """
        if len(colors) < n_volcs:
            colors = colors + colors
        else:
            return colors
        # check if it needs extending again, recursive
        colors = color_extender(colors, n_volcs)
        return colors
    
    def plot_for_subset(volcs_subset, ax):
        """Plot each volcano's data on ax and add a subplot legend."""
        for volc in volcs_subset:
            line, = ax.plot(
                volc.status_overall['dates'],
                volc.status_overall['unrest_metric'],
                label=volc.name,
                c=colors[color_counter[0]]
            )
            lines.append(line)
            color_counter[0] += 1
            
        # add the legend for that subplot
        ax.legend(loc='upper left', ncol=1, frameon=True, prop={'size': 10})
    
    # Get a list of colors, extended if necessary.
    colors = list(mcolors.TABLEAU_COLORS)
    colors = color_extender(colors, len(volcs))
    
    # Mutable counter for colors.
    color_counter = [0]
    
    # List to store all line objects
    lines = []
    
    #
    y_label = 'LiCSAlert Unrest metric'
    
    # Create subplots
    fig, axes = plt.subplots(n_rows, 1, figsize=figsize, )
    plt.subplots_adjust(hspace = 0.4)
    
    for row in range(n_rows):
        subset = volcs[(row * n_per_row): ((row + 1) * n_per_row)]
        plot_for_subset(subset, axes[row])
    
    # Set the x and y limits and add grid lines to all axes.
    for ax in axes:
        ax.set_xlim(left=x_lims[0], right=x_lims[1])
        if y_lims is not None:
            ax.set_ylim(bottom=y_lims[0], top=y_lims[1])
        ax.grid(
            True, which='both', linestyle='--', linewidth=0.5, color='gray', 
            alpha=0.5
            )
        ax.set_ylabel(y_label)
    
    # Only leave x tick labels on the bottom subplot.
    for ax in axes[:-1]:
        ax.set_xticklabels([])
    
    
    # shade the baseline stage:
    for ax in axes:
        ax.axvspan(x_lims[0], baseline_end, color='grey', alpha=0.1)
    
    # x label
    axes[-1].set_xlabel('Date')
    
    # For first subplot using labels_dict1
    if labels_dict1 is not None:
        y_bottom, y_top = axes[0].get_ylim()
        offset = arrow_length * (y_top - y_bottom)
        for label_text, date_str in labels_dict1.items():
            date_obj = datetime.strptime(date_str[0], "%Y_%m_%d")
            axes[0].annotate(
                label_text,
                xy=(date_obj, y_top),                # arrow tip at the top edge of subplot 0
                xytext=(date_obj, y_top + offset),     # text placed above the axes
                textcoords='data',
                arrowprops=dict(arrowstyle="->"),
                rotation=text_angle,                           # rotated text
                ha='center',
                va='bottom',
                clip_on=False,
                color = date_str[1]
            )
        
    # For second subplot using labels_dict2
    if labels_dict2 is not None:
        y_bottom2, y_top2 = axes[1].get_ylim()
        offset2 = arrow_length * (y_top2 - y_bottom2)
        for label_text, date_str in labels_dict2.items():
            date_obj = datetime.strptime(date_str[0], "%Y_%m_%d")
            axes[1].annotate(
                label_text,
                xy=(date_obj, y_top2),                # arrow tip at the top edge of subplot 1
                xytext=(date_obj, y_top2 + offset2),     # text placed above the axes
                textcoords='data',
                arrowprops=dict(arrowstyle="->"),
                rotation=0,
                ha='center',
                va='bottom',
                clip_on=False,
                color = date_str[1]
        )
    
    plt.show()
    
    fig.savefig('RSE_fig.png', bbox_inches='tight')
    fig.savefig('RSE_fig.pdf', format='pdf', dpi=1200, bbox_inches='tight')
    
    pdb.set_trace()

# Example dictionaries for annotations.
labels_eruptive = {
    'Mauna Loa': ['2022_11_27', 'tab:orange'],
    'Wolf': ['2022_01_07', 'tab:green'],
    'Fagradalsfjall 1': ['2021_03_19', 'tab:blue'],
    'Fagradalsfjall 2': ['2021_12_22', 'tab:blue'],
    'Fagradalsfjall 3': ['2022_08_03', 'tab:blue'],
    'Fagradalsfjall 4': ['2023_07_10', 'tab:blue'],
    'Svartsengi': ['2023_12_18', 'tab:blue'],
}

labels_no_erupt = {
    'Askja inflation': ['2021_01_01', 'tab:red'],
    #'Campi Flegrei acc. 1': '2021_02_01',
    'Campi Flegrei\nacceleration': ['2024_01_01', 'tab:purple'],
    'Sierra Negra\nacceleration': ['2022_01_01', 'tab:brown']
}

# figure
unrest_metric_all_volcs_publication_figure(
    [volcs[i] for i in fig_indexes], 3, 3,
    x_lims=[datetime(2019, 1, 1), datetime(2025, 1, 1)], 
    baseline_end=datetime(2021, 1, 1),
    labels_dict1=labels_eruptive,
    labels_dict2=labels_no_erupt,
    y_lims=None, figsize=(7.84, 9)
)


#%% SI figure showing all the other volcs

fig_indexes_si = [i for i in range(len(volcs)) if i not in fig_indexes]

unrest_metric_all_volcs_publication_figure(
    [volcs[i] for i in fig_indexes_si], 5, 3,
    x_lims=[datetime(2019, 1, 1), datetime(2025, 1, 1)],
    baseline_end=datetime(2021, 1, 1),
    y_lims=None
)





