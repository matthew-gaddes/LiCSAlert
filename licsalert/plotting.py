#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 15:57:28 2023

@author: matthew
"""

import pdb

#%%

def LiCSAlert_aux_figures(displacement_r2_current, reconstructions, residuals, tbaseline_info,
                          figure_type, figure_out_dir):
    """
    Plot last cumulative ifg, last incremetnal ifg, reconstrution for last incremental ifg, and residual for last incremental ifg
    
    Inputs:
        displacement_r2_current | dict | licsalert dict of ifgs
        reconstrutions          | r2 array | reconstructions as row vectors
        residuals               | r2 array | reconstructions as row vectors
        tbasline_info           | dict | licsalert dict of tbaseline info
        figure_type             | string | png / window / both
        figure_out_dir          | Path   | out dir path.  
        
    Returns:
        figure
        
    History: 
        2023_04_03 | MEG | Written
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from licsalert.aux import col_to_ma
    
    
    # 0 Check matplotlib backend is set correctly 
    if figure_type == 'png':
        plt.switch_backend('Agg')                                                           #  works when there is no X11 forwarding, and when displaying plots during creation would be annoying.  
    else: 
        if plt.get_backend() != 'Qt5Agg':                                                               # check what the backend is 
            plt.switch_backend('Qt5Agg')                                                           #  and switch to interactive if it wasn't already.  
        
    ifg_n = displacement_r2_current['incremental'].shape[0]
    inc_ifg_date = tbaseline_info['ifg_dates'][ifg_n-1]
    cum_ifg_date = f"{tbaseline_info['ifg_dates'][0][:8]}_{tbaseline_info['ifg_dates'][ifg_n-1][-8:]}"
    
    # 1: do the plots
    plot_1_image(np.sum(displacement_r2_current['incremental'], axis = 0), displacement_r2_current['mask'], f"01_cumulative_{cum_ifg_date}", 
                 figure_type, figure_out_dir)
    plot_1_image(displacement_r2_current['incremental'][-1, :], displacement_r2_current['mask'], f"02_incremental_{inc_ifg_date}",
                 figure_type, figure_out_dir)
    plot_1_image(reconstructions[:,-1], displacement_r2_current['mask'], f"03_reconstruction_{inc_ifg_date}",
                 figure_type, figure_out_dir)
    plot_1_image(residuals[:, -1], displacement_r2_current['mask'], f"04_residual_{inc_ifg_date}",
                 figure_type, figure_out_dir)
    
    if plt.get_backend() != 'Qt5Agg':                                                           # check if we need to reset the backend (to interactive window)
        plt.switch_backend('Qt5Agg')        
        
        
#%%

def plot_1_image(im_r1, mask, title, figure_type, figure_out_dir, figsize = (18,9)):
    """
    """
    import matplotlib.pyplot as plt
    from licsalert.aux import col_to_ma
    
    
    f, ax = plt.subplots(1,1, figsize = figsize)
    im = ax.matshow(col_to_ma(im_r1, mask))
    f.colorbar(im, label = 'Displacement (m)')
    f.suptitle(title)
    
    if (figure_type == 'png') or (figure_type == 'both'):
        f.savefig(figure_out_dir / f"{title}.png", bbox_inches='tight')
        
    if figure_type == 'png':
        plt.close(f)