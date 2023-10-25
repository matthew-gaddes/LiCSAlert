#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:09:28 2020

@author: matthew
"""
import pdb



#%%

def save_licsalert_aux_data(out_dir, displacement_r2):
    """ Save the displacment_r2 dict containing all the time series 
    data that was used by LiCSAlert.
    
    Inputs:
        out_dir \ pathlib Path \ directory where LiCSAlert outputs are written to. 
        displacement_r2 | dict | dict that contains dict_keys(['dem', 'mask', 'incremental', 'lons', 'lats', 'E', 'N', 'U', 'incremental_downsampled', 'mask_downsampled'])
    Returns:
        pickle file
    History:
        2023_10_25 | MEG | Written
    """                
    import pickle
    
    with open(out_dir / "aux_data_figs" / 'original_ts_data.pkl', 'wb') as f:
        pickle.dump(displacement_r2, f)
    f.close()    
    
