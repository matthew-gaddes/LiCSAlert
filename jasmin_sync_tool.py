#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 11:44:57 2024

@author: matthew
"""

print(f"Succesfully started the script.  ")

from pathlib import Path
import pdb
import numpy as np
from glob import glob 
import sys
import subprocess

from licsalert.jasmin_tools import open_comet_frame_files, volcano_name_to_comet_frames
from licsalert.jasmin_tools import write_jasmin_download_shell_script



#%% Things to set


comet_volcano_frame_index_dir = Path('./comet_volcano_frames/')

cov_volc_names = ['Antisana',
                 'Cerro Overo',
                 'Cordon del Azufre',
                 'Cotopaxi',
                 'Domuyo',
                 'Fernandina',
                 #'Guagua Pichinca',                    # volcano portal spelling
                 'Guagua Pichincha',                     # licsbas spelling
                 'Laguna del Maule',
                 'Masaya',
                 'Nevados de Chillan',
                 'Planchón-Peteroa',                    # accent o?
                 'Sabancaya',
                 'San Miguel',
                 'San Salvador',
                 'Sangay',
                 'Santa Ana',
                 'Sierra Negra',
                 'Tungurahua',
                 'Turrialba',
                 'Villarrica']

jasmin_dir = Path('mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert')
local_dir = Path('./licsalert_sync/')
bash_sync_fle = "sync_script.sh"


#%%

# open the info on COMET volcano frames and what region they're in.  
comet_volcano_frame_index = open_comet_frame_files(comet_volcano_frame_index_dir)

# get the frames for the volcanoes of interest.  
cov_volcs = volcano_name_to_comet_frames(cov_volc_names, comet_volcano_frame_index)    
        
# write a shell script to download them.  
write_jasmin_download_shell_script(jasmin_dir, local_dir, bash_sync_fle,
                                   cov_volcs, omit_json = True)