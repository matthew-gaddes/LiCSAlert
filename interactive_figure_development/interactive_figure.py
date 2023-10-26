#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 14:49:10 2023

@author: matthew

Thoughts:
    
    reconstructed on left, cumulative on right, 
    
    ICs - get these from aux figures.  
    
    time courses?
    
    radio buttons to select ICs.  
    
    


"""

from numpy import pi, sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

from pathlib import Path
import pickle
import sys
from glob import glob

#%%

print("Running")

#%%

#def licsalert_ic_combiner(licsalert_dir):
    
    
    
licsalert_dir = Path("./../campi_flegrei_example")
    
    # open the licsbas data.  
    
sys.path.append("./..")
    
import licsalert
from licsalert.aux import col_to_ma
    
def open_aux_data(licsalert_dir):
    """Open all the data stored in the pickle files in aux_images_data
    Inputs:
        licsalert_dir | pathlib Path | output directory when LiCSAlert was run.  
    Returns:
        displacement_r2 | dict | contains (['dem', 'mask', 'incremental', 'lons', 'lats', 'E', 'N', 'U', 'incremental_downsampled', 'mask_downsampled'])
        aux_data | dict | ['icasar_sources', 'dem', 'mask'])
    History:
        2023_10_25 | MEG | Written. 
    """
    
    with open(licsalert_dir / "aux_data_figs" / 'original_ts_data.pkl', 'rb') as f:
        displacement_r2 = pickle.load(f)
    f.close()
    
    with open(licsalert_dir / "aux_data_figs" / 'aux_images_data.pkl', 'rb') as f:
        aux_data = pickle.load(f)
    f.close()
    
    return displacement_r2, aux_data
    
    
    
def open_tcs(licsalert_dir):
    """ Open the time course data.  
    Inputs:
        licsalert_dir | pathlib Path | output directory when LiCSAlert was run.  
    Returns:
        sources_tcs | list of dicts | One item in list for each source, each item contains ['cumulative_tc', 'gradient', 'lines', 'sigma', 'distances', 't_recalculate'])
    History:
        2023_10_25 | Written | MEG
    """
    
    licsalert_items = sorted(glob(str(licsalert_dir / '*')))
    
    # remove any items that are not a licsalert data directory.  
    delete_args = []
    for item_n, licsalert_item in enumerate(licsalert_items):
        item_name = Path(licsalert_item).parts[-1]
        if item_name in ["ICASAR_results", "LiCSAlert_history.txt", "aux_data_figs"]:
            delete_args.append(item_n)
    
    for delete_arg in delete_args[::-1]:
        del licsalert_items[delete_arg]
    
    final_date_dir = Path(sorted(licsalert_items)[-1])
    
    with open(final_date_dir / 'time_course_info.pkl', 'rb') as f:
        sources_tcs = pickle.load(f)
        residual_tcs = pickle.load(f)
    f.close()

    return sources_tcs


    
    
displacement_r2, aux_data = open_aux_data(licsalert_dir)
sources_tcs = open_tcs(licsalert_dir)    

f, axes = plt.subplots(1,2)

cum_ifg_r1 = np.sum(displacement_r2['incremental'], axis = 0)
axes[1].matshow(col_to_ma(cum_ifg_r1, displacement_r2['mask']))
    
    
# switch to a gridspec.  
# plot ICS and cumulative tcs.  
# get the means. 
# function to reconstruct data
# slection button (can select multiople)
# redo reconstrutcion on button change.  
# tica vs sica?  

    
    
sys.exit()
    
    


def signal(amp, freq):
    return amp * sin(2 * pi * freq * t)

axis_color = 'lightgoldenrodyellow'

fig, axes = plt.subplots(1,2)

ax = axes[0]



# Adjust the subplots region to leave some space for the sliders and buttons
fig.subplots_adjust(left=0.25, bottom=0.25)

t = np.arange(0.0, 1.0, 0.001)
amp_0 = 5
freq_0 = 3

# Draw the initial plot
# The 'line' variable is used for modifying the line later
[line] = ax.plot(t, signal(amp_0, freq_0), linewidth=2, color='red')
ax.set_xlim([0, 1])
ax.set_ylim([-10, 10])

# Add two sliders for tweaking the parameters

# Define an axes area and draw a slider in it
amp_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axis_color)
amp_slider = Slider(amp_slider_ax, 'Amp', 0.1, 10.0, valinit=amp_0)

# Draw another slider
freq_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axis_color)
freq_slider = Slider(freq_slider_ax, 'Freq', 0.1, 30.0, valinit=freq_0)

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    line.set_ydata(signal(amp_slider.val, freq_slider.val))
    fig.canvas.draw_idle()
amp_slider.on_changed(sliders_on_changed)
freq_slider.on_changed(sliders_on_changed)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    freq_slider.reset()
    amp_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

# Add a set of radio buttons for changing color
color_radios_ax = fig.add_axes([0.025, 0.5, 0.15, 0.15], facecolor=axis_color)
color_radios = RadioButtons(color_radios_ax, ('red', 'blue', 'green'), active=0)
def color_radios_on_clicked(label):
    line.set_color(label)
    fig.canvas.draw_idle()
color_radios.on_clicked(color_radios_on_clicked)

plt.show()