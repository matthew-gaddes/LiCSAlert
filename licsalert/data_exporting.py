#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:09:28 2020

@author: matthew
"""
import pdb



#%%

def save_licsalert_aux_data(
        out_dir,
        displacement_r3,
        tbaseline_info,
        sica_tica,
        ):
    """ Save the displacment_r2 dict containing all the time series 
    data that was used by LiCSAlert.
    
    Inputs:
        out_dir \ pathlib Path \ directory where LiCSAlert outputs are written to. 
        displacement_r3 | dict | dict that contains dict_keys
                        (['dem', 'mask', 'inc_ma', 'lons', 'lats', 'E', 'N',
                          'U', 'inc_ma_downsampled', 'cum_ma', 'cum_ma_downsampled'])
        tbaseline_info | dict | contains dict_keys(['acq_dates', 'ifg_dates', 'baselines', 'baselines_cumulative'])
        sica_tica | str | record how ICA was performed 
    Returns:
        pickle file
    History:
        2023_10_25 | MEG | Written
        2025_07_03 | MEG | Update to use rank 3 arrays.  
    """                
    # v1: don't rely on cloudpickle, but can't save custom object
    # import pickle
    # import copy
    
    # displacement_r3_out = copy.deepcopy(displacement_r3)
    
    # # remove the MeanCeneteredArrays (custom LiCSAlert objects)
    # del displacement_r3_out['cum_ma'], displacement_r3_out['inc_ma']
    
    # # replace them with a simple masked array of the data
    # displacement_r3_out['cum_ma'] = displacement_r3['cum_ma'].original
    # displacement_r3_out['inc_ma'] = displacement_r3['inc_ma'].original
    
        
    # with open(out_dir / "aux_data_figs" / 'original_ts_data.pkl', 'wb') as f:
    #     pickle.dump(displacement_r3_out, f)
    #     pickle.dump(tbaseline_info, f)
    # f.close()    
    displacement_r3['sica_tica']=sica_tica
    
    import cloudpickle
    
    with open(out_dir / "aux_data_figs" / 'original_ts_data.pkl', 'wb') as f:
        cloudpickle.dump(displacement_r3, f)
        cloudpickle.dump(tbaseline_info, f)
    f.close()    
    
