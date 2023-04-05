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
    
    # if plt.get_backend() != 'Qt5Agg':                                                           # check if we need to reset the backend (to interactive window)
    #     plt.switch_backend('Qt5Agg')        
        
        
#%%

def plot_1_image(im_r1, mask, title, figure_type, figure_out_dir, figsize = (18,9)):
    """Plot a single image (column or row vector) when also given its mask.  
    
    """
    import matplotlib.pyplot as plt
    from licsalert.aux import col_to_ma
    
    
    f, ax = plt.subplots(1,1, figsize = figsize)
    im = ax.matshow(col_to_ma(im_r1, mask))
    f.colorbar(im, label = 'Displacement (m)')
    f.suptitle(title)
    
    if (figure_type == 'png') or (figure_type == 'both'):
        f.savefig(figure_out_dir / f"{title}.png", bbox_inches='tight')
        
    # if figure_type == 'png':
    #     plt.close(f)
    
    
#%%


def plot_mask_changes(icasar_mask, licsbas_mask, mask_combined, licsbas_date, current_output_dir, figure_type):
    """ Create a .png showing the licsbas mask, the ICASAR mask, and the current combined mask (ie the pixels in both).  
    
    Inputs:
        icasar_mask | r2 array | the mask used by ICASAR
        licsbas_mask | r2 array | the mask produced by the last run of LiCSBAS
        mask_combined | r2 array | the mask that removes any pixels that aren't in bothh the sources and the ifgs
        licsbas_date | string | the date that LiCSAlert is being run to.  
        current_output_dir | Path | the folder that LiCSALert is currently outputting to
    Returns:
        .png figure
    History:
        2020/06/25 | MEG | Written
        2020/07/01 | MEG | Major rewrite to suit directory based structure.  
        2020/07/03 | MEG | continue major rewrite, and write docs.  
        2021_10_20 | MEG Simplify for LiCSAlert 2.0
    """
    import matplotlib.pyplot as plt
    
    # 0 Check matplotlib backend is set correctly 
    if figure_type == 'png':
        plt.switch_backend('Agg')                                                           #  works when there is no X11 forwarding, and when displaying plots during creation would be annoying.  
    else: 
        if plt.get_backend() != 'Qt5Agg':                                                               # check what the backend is 
            plt.switch_backend('Qt5Agg')                                                           #  and switch to interactive if it wasn't already.  
    
    title = f"LiCSBAS_last_date_{licsbas_date}"

    # 1 Figure showing the masks
    f1,axes = plt.subplots(1,3, figsize = (12,6))
    axes[0].imshow(icasar_mask)
    axes[0].set_title('(ICASAR) sources mask')
    axes[1].imshow(licsbas_mask)
    axes[1].set_title('LiCSBAS mask')
    axes[2].imshow(mask_combined)
    axes[2].set_title('Current combined mask')
    f1.suptitle(title)
    
    f1.canvas.manager.set_window_title(title)
    if (figure_type == 'png') or (figure_type == 'both'):
        f1.savefig(current_output_dir / "mask_status.png", bbox_inches='tight')
