#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 17:55:25 2020

@author: matthew
"""


def downsample_ifgs(ifgs, mask, scale = 0.1, verbose = True):
    """ A function to take ifgs as row vectors (and their associated mask) and return them downsampled (for fast plotting)
    Inputs:
        ifgs | rank 2 array | ifgs as rows
        mask | rank 2 mask | to convert a row interferogram into a rank 2 masked array
        scale | flt | <1 and downsample, >1 might make it upsample/interpolate?  Not tested
    Outputs:
        ifgs_ds | rank 2 array | downsampled ifgs as rows
        mask_ds | rank 2 mask | for converting ifgs_ds row vectors into rank 2 masked arrays 
        
    2018/03/?? | MEG | written
    2018/07/09 | MEG | update skimage.transform.rescale arguments to supress warnings.  
    2020/03/08 | MEG | Major rewrite to deal with smearing/interpolating of the masks when using integer instead of boolean values.  
    """
    
    import numpy as np
    from skimage.transform import rescale
    import numpy.ma as ma
    from licsalert.aux import col_to_ma
    
    # 1: Check inputs
    if np.array_equal(mask, mask.astype(bool)):                                                                         # force the user to use a boolean mask
        pass
    else:
        raise Exception(f"The 'mask' must contain boolean values.  I.e., not 0s and 1s.  Exiting....")
    
    # 2: initate some items  
    n_pixels = np.sum(np.logical_not(mask))                                                                             # number of pixels is sum of False (not masked) pixels, and convert False to True with Not
    n_ifgs = np.size(ifgs, axis = 0)                                                                                    # get no. of ifgs
       
    # 3: Downsample the mask, and make sure it stays boolean
    mask_ds = rescale(mask, scale, multichannel = False, anti_aliasing = False).astype(bool)
    n_pixels_ds = np.sum(np.logical_not(mask_ds))                                                                       # get number of pixels that are not masked
    
    if verbose:
        print(f'Interferograms are being downsampled from {n_pixels} pixels to {n_pixels_ds} pixels.  ')
    
    # 4: Downsample the ifgs by looping through them
    ifgs_ds = np.zeros((n_ifgs, n_pixels_ds))                                                                           # initiate array to store rows 
    for i, single_ifg in enumerate(ifgs):                                                                               # loop through each ifg (which is a row)
        ifg_ma = col_to_ma(single_ifg, mask)                                                                            # make into a rank 2 masked array
        ifg_rescale = rescale(ifg_ma, scale, multichannel = False, anti_aliasing = False)                               # rescale, no longer a ma
        ifg_rescale_ma = ma.array(ifg_rescale, mask = mask_ds)                                                          # convert back to ma
        ifgs_ds[i,:] = ma.compressed(ifg_rescale_ma)                                                                    # back to being a row vector
        
    return ifgs_ds, mask_ds
