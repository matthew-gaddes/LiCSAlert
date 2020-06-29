#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 18:10:21 2020

@author: matthew
"""


#%%

#%%

def create_folder(folder):
    """ Try to create a folder to save function outputs.  If folder already exists,
    funtion will try to delete it and its contents.  
    
    Inputs: 
        folder | string |path to new folder.  e.g. './my_folder' 
    Returns:
        new folder
    History:
        2020/06/25 | MEG | Written
    """
    import shutil
    import os
    try:
        print(f"Trying to remove the existing outputs folder ({folder})... ", end = '')
        shutil.rmtree(folder)                                                                       # try to remove folder
        print('Done!')
    except:
        print("Failed!")                                                                                # 
    try:
        print(f"Trying to create a new outputs folder ({folder})... ", end = '')                                    # try to make a new folder
        os.mkdir(folder)                                                                       
        print('Done!')
    except:
        print("Failed!") 
        
        
        #%%


def col_to_ma(col, pixel_mask):
    """ A function to take a column vector and a 2d pixel mask and reshape the column into a masked array.  
    Useful when converting between vectors used by BSS methods results that are to be plotted
    
    Inputs:
        col | rank 1 array | 
        pixel_mask | array mask (rank 2)
        
    Outputs:
        source | rank 2 masked array | colun as a masked 2d array
    
    2017/10/04 | collected from various functions and placed here.  
    
    """
    import numpy.ma as ma 
    import numpy as np
    
    source = ma.array(np.zeros(pixel_mask.shape), mask = pixel_mask )
    source.unshare_mask()
    source[~source.mask] = col.ravel()   
    return source


def add_square_plot(x_start, x_stop, y_start, y_stop, ax, colour = 'k'):
    """Draw localization square around an area of interest, x_start etc are in pixels, so (0,0) is top left.  
    Inputs:
        x_start | int | start of box
        x_stop | int | etc. 
        y_start | int |
        y_ stop | int |
        ax | axes object | axes on which to draw
        colour | string | colour of bounding box.  Useful to change when plotting labels, and predictions from a model.  
    
    Returns:
        box on figure
        
    History:
        2019/??/?? | MEG | Written
        2020/04/20 | MEG | Document, copy to from small_plot_functions to LiCSAlert_aux_functions
    """
        
    ax.plot((x_start, x_start), (y_start, y_stop), c= colour)           # left hand side
    ax.plot((x_start, x_stop), (y_stop, y_stop), c= colour)             # bottom
    ax.plot((x_stop, x_stop), (y_stop, y_start), c= colour)             # righ hand side
    ax.plot((x_stop, x_start), (y_start, y_start), c= colour)             # top