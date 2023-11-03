#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 14:49:10 2023

@author: matthew

Thoughts:

    
    plot ICS and cumulative tcs

    add the dem
    
"""

from numpy import pi, sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

from pathlib import Path
import pickle
import sys
from glob import glob
import pdb


debug_scripts = "/home/matthew/university_work/python_stuff/python_scripts"

if debug_scripts not in sys.path:                                                                             # check if already on path
    sys.path.append(debug_scripts)
from small_plot_functions import matrix_show, quick_linegraph

#%%

print("Running")

#%%

#def licsalert_ic_combiner(licsalert_dir):
    
    
    
licsalert_dir = Path("./../001_campi_flegrei_example")
    
    # open the licsbas data.  
    
sys.path.append("./..")
  

import licsalert
from licsalert.plotting import licsalert_results_explorer

licsalert_results_explorer(licsalert_dir, fig_width = 18)

#%% 


    


    
#%%
    