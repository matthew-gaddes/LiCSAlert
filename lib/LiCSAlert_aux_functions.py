#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 18:10:21 2020

@author: matthew
"""

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