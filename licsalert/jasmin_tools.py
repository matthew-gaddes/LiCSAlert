#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:35:43 2024

@author: matthew
"""
import pdb
from glob import glob
from pathlib import Path

def open_comet_frame_files(comet_volcano_frame_index_dir):
    """ Given a directory that contains files of the name of each COMET
    volcano frame and which region it's in, turn these into a dictionary of 
    lists, where each region is a key and its value is a list of frames in that
    region.  
    
    Inputs:
        comet_volcano_frame_index_dir | Path | direcotry containing text files
                                              of comet frame names
    Returns:
        comet_volcano_frame_index | dict | each region is a key, and its value
                                            is a list of the frames in that region.  
                                            
    History:
        2024_01_18 | MEG | Written.  
    """
    
    comet_volcano_frame_index = {}

    volcano_frame_files = sorted(glob(str(comet_volcano_frame_index_dir / '*')))
    
    for vff in volcano_frame_files:

        # get the region name from the path        
        region = Path(vff).parts[-1].split('.')[0]
        
        # Open the file and read its lines into a list
        with open(vff, 'r') as file:
            lines = file.readlines()
        
        # Remove newline characters from each line (if needed)
        lines = [line.strip() for line in lines]
        
        comet_volcano_frame_index[region] = lines
        
    return comet_volcano_frame_index

#%%

def volcano_name_to_comet_frames(volc_names, comet_volcano_frame_index):
    """ Given a list of volcanoes of interest, find the frames that name them
    from the dictionary of comet volcano frames
    
    Inputs:
        volc_names | list | volcano names, can have accents, capital, and spaces. 
        comet_volcano_frame_index | dict | each region is a key, and its value
                                            is a list of the frames in that region.  
                                            
    Returns:
        volc_frames | dict | each volcano name is a key, and the value is a list
                                of the comet frames for that volcano.  
        
    History:
        2024_01_18 | MEG | Written
    """
    
    
        
    class comet_volcano:
        def __init__(self, name):
            """ 
            """
            self.name = name
            self.computer_name = self.create_simple_name()
            self.frames = []
            
            
        def create_simple_name(self):
            """
            """
            from unidecode import unidecode
            # convert to snake case and all lowercase as this is format used by comet 
            volc_name_converted = volc_name.lower().replace(' ', '_')
            # swap any accents to nearest ASCII 
            volc_name_converted = unidecode(volc_name_converted)
            return volc_name_converted
        
        
        
    volc_frames = []
    for volc_name in sorted(volc_names):
        volc = comet_volcano(volc_name)
        
        # loop through each region
        for region in list(comet_volcano_frame_index.keys()):
            
            # get all the frames for that region
            frames = comet_volcano_frame_index[region]
            
            # loop through the frames to see if the volcano name is in the frame name
            for frame in frames:
                if volc.computer_name in frame:
                    volc.region = region
                    volc.frames.append(frame)
        volc_frames.append(volc)
        
    return volc_frames



#%%
def write_jasmin_download_shell_script(jasmin_dir, local_dir, out_file,
                                       volcs, exclude_json_gz =True,
                                       exclude_original_ts = True,
                                       exclude_fastica = True):
    """ Given a dictionary of volcano frames, create a shell script to download
    the LiCSAlert directories for each one.  
    Inputs:
        jasmin_dir | Path | path to licsalert products on Jasmin, including server.  
                            e.g. 'mgaddes@xfer1.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/
                            projects/LiCS/volc-portal/processing_output/licsalert'
                            
        local_dir | Path | directory that shell script will sync into
        out_file | Path | path to outfile (i.e. includes name)
        volc | list | list of comet_volcano objects.  One item for each volcano,
                        .frames returns all the frames for that volcano.  
        exclude_json_gz | Boolean | if True, all the .json.gz files that the 
                                    web portal uses are not added to the download
                                    script.  
        exclude_original_ts | Boolean | if True, aux_data_figs/original_ts_data.pkl
                                        is not added to the download script for
                                        each frame.  
        exclude_fastica | Boolean | if True, ICASAR_results/FastICA_results.pkl
                                    is not added to the download script for each
                                    frame
        
    Returns:
        shell script
    History:
        2023_01_19 | MEG | Written.  
    """

    # this exluces the directory of json files for each licsalert date, or doesn't

    
    if exclude_json_gz:
        json_gz = "--exclude '*json.gz'"
    else:
        json_gz = ""
    
    if exclude_original_ts:
        original_ts = "--exclude 'original_ts_data.pkl'"
    else:
        original_ts = ""
        
    if exclude_fastica:
        fastica = "--exclude 'FastICA_results.pkl'"
    else:
        fastica = ""
    
    json_section = f"{json_gz} {original_ts} {fastica}"
    
    with open(out_file, "w") as script_file:
        # Write shebang (optional, depends on your use case)
        script_file.write("#!/bin/bash\n\n")
        # Add a warning to the user about how to use the shell script
        script_file.write(f"# Note - this will only work if you can login to Jasmin\n")
        script_file.write(f"# which requires running a command such as:\n")
        script_file.write(f"# eval $(ssh-agent -s); ssh-add ~/.ssh/id_rsa_jasmin\n")
        script_file.write(f"# after you have setup a key pair with Jasmin.\n\n\n")
                          

        # iterate over each volcano
        for volc in volcs:
            # and the grames for that volcano
            for frame in volc.frames:
                # build the command for that frame        
                command = (f"rsync -av {json_section} "
                           f"{jasmin_dir / volc.region / frame} "
                           f"{local_dir / volc.region }")
                script_file.write(command + "\n")
