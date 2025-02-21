#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:35:43 2024

@author: matthew
"""
import pdb
from glob import glob
from pathlib import Path


#%%




def move_small_directories(indir, outdir, n, regions = True):
    """
    Move subdirectories with fewer than `n` items to `outdir`.
    
    This is useful when working with Jasmin licsalert dirs that will only 
    have a few products if LiCSAlert failed, compared to one with hundreds
    if it ran succesfully.  
    
    Parameters:
        indir (str): The path to the main directory containing subdirectories.
        outdir (str): The path to the output directory.
        n (int): The minimum number of items a subdirectory must have to remain.
        
    Returns:
        nan
        
    History:
        2024_??_?? | MEG | Written
        2025_02_19 | MEG | Update to include regions. 
    """
    import os
    import shutil
    
    # Ensure the output directory exists
    os.makedirs(outdir, exist_ok=True)

    # if we have regions, get the path to each region directory
    if regions:
        region_names = sorted(os.listdir(indir))
        region_dirs = [indir / region_name for region_name in region_names]
    
    # if we don't, just say the indir is the region dir.  
    else:
        region_dirs = [indir]
        
    
    for region_dir in region_dirs:
        volc_names = sorted(os.listdir(region_dir))
        
        for volc_name in volc_names:
            # get the path to current volcano
            volc_dir = sorted(os.listdir(region_dir / volc_name))
            
            # Check if it's a directory
            if os.path.isdir(str(volc_dir)):
                
                # Count items in the subdirectory
                item_count = len(os.listdir(volc_dir))
                
                # Move the directory if it has fewer than `n` items
                if item_count < n:
                    print(f"Moving {volc_dir} to {outdir}")
                    shutil.move(volc_dir, outdir)
            
            




#%%


def get_lon_lat_of_volcs(volcs):
    """ Given a list of volcs, add the lon and lat of the centre of the 
    LiCSBAS data (averaged across all the frames)
    
    Inputs:
        volcs | list of comet volcs |
    Returns:
        volcs | list of comet volcs | now with lon_lat attribute
    History:
        2024_06_12 | MEG | Written
        
    """
    import numpy as np
    
    from licsalert.data_importing import open_aux_data
    
    print(f"Finding the lon and lat of the centre of each volcano.  ")
    
    for volc in volcs:
        lons = []
        lats = []
        for frame_n, frame in enumerate(volc.frame_status):
            print(
                f"Opening the lon and lat for {volc.frames[frame_n]}...", 
                end = ''
                )
            try:
                # open the time series 
                displacement_r2, tbaseline_info, aux_data = open_aux_data(
                    Path(frame)
                    )
                # get the shape of the images
                ny, nx = displacement_r2['lons'].shape
                # get the lon and lat of the centre pixel
                lons.append(displacement_r2['lons'][int(ny/2), int(nx/2)])
                lats.append(displacement_r2['lats'][int(ny/2), int(nx/2)])
                print("Succeeded")
            except:
                print(f"Failed")
                    
                
        # average across the frames for that volcano
        if (len(lons) > 0) and (len(lats) > 0):
            volc.lon_lat = (sum(lons)/ len(lons), sum(lats) / len(lats))
        else:
            volc.lon_lat = (np.nan, np.nan)
            # warn user if we failed.  
            print(
                f"Failed to extract the lon and lat for any of the frames "
                f"of {volc.name} (volc.region)"
                )
        print("\n")
    
    return volcs
    

#%%


def update_volcs_with_data(volcs, volc_frame_dirs, volc_frame_names,
                           verbose = False):
    """ Given the LiCSAlert volcs list (list of comet volcano items),
    update to contain only ones that actually have data.  
    
    Inputs:
        volcs | list | one item for each volcano
        volc_frame_dirs | list | each item is path to directory of licsalert 
                                results
        volc_frame_names | list | names of frames in above list.  Must be the 
                                same length
                                
    Returns:
        volcs | list | one item for each volcano, volcanoes with no data 
                        removed.  frame_status attribute now has path to data,
                        if it exists and NA if it doesn't 
                        
    History:
        2024_06_12 | MEG | Written
    
    """
    
    def contains_only_na(lst):
        for item in lst:
            if item != "NA":
                return False
        return True
    
    del_indexes = []
    # loop through each volcano
    for volc_n, volc in enumerate(volcs):
        
        frame_status = []
        # and each frame that that volcanoe has
        for frame in volc.frames:
            
            # to see if there is licsalert data for that frame
            if  frame in volc_frame_names:
                # get the index of the frame
                index = volc_frame_names.index(frame)
                frame_status.append(volc_frame_dirs[index])
            else:
                frame_status.append('NA')
        volc.frame_status = frame_status
        
        # check if we should remove the volcano as there's no data
        if contains_only_na(frame_status):
            del_indexes.append(volc_n)

    # Sort the indexes in descending order to avoid shifting issues, then delete
    for index in sorted(del_indexes, reverse=True):
        del volcs[index]
        
    # possibly update user on terminal
    if verbose:
        for volc in volcs:
            print(volc.name)
            print(volc.frames)
            print(volc.frame_status)
            print('\n')
        
        
    return volcs
    

#%%

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
        2024_06_12 | MEG | Fix bug that if volcan name was within other name, 
                           that would be returned (e.g. if volcano is abc, abcdef
                           would also be returned. )     
    """
    
    
        
    class comet_volcano:
        def __init__(self, name):
            """ 
            """
            self.name = name
            self.computer_name = self.create_computer_name()
            self.frames = []
            
            
        def create_computer_name(self):
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
                # the last 18 digits of a frame are the track / burst etc.  
                # before that is the volcano name in snake_case
                if volc.computer_name == frame[:-18]:
                    volc.region = region
                    volc.frames.append(frame)
        volc_frames.append(volc)
        
    return volc_frames



#%%
def write_jasmin_download_shell_script(
        jasmin_dir, local_dir, out_file,volcs, 
        exclude_json_gz =True, exclude_original_ts = True,
        exclude_ICASAR_data = True,
        exclude_epoch_images = True, exclude_epoch_data = True, 
        ):
    """ Given a dictionary of volcano frames, create a shell script to download
    the LiCSAlert directories for each one.  
    Inputs:
        jasmin_dir | Path | path to licsalert products on Jasmin, including server.  
                            e.g. 'mgaddes@xfer-vm-01.jasmin.ac.uk:/gws/nopw/j04/nceo_geohazards_vol1/
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
        exclude_ICASAR_data | Boolean | if True, ICASAR_results/FastICA_results.pkl
                                        and ICASAR_results.pkl and pca_results.pkl
                                        are not transferred.  
        exclude_epoch_iamges | Boolean | If True, the 01.. 02.. images
                                        for each epoch are not transferred.  
        exclude_epoch_data | Boolean | If True, the two .pkl files for 
                                        each epoch are not transferred.  
        
    Returns:
        shell script
    History:
        2023_01_19 | MEG | Written.  
        2025_02_13 | MEG | Add ability to exclude more files. 
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
        
    if exclude_ICASAR_data:
        icasar_data = ("--exclude 'FastICA_results.pkl' "
                       "--exclude 'ICASAR_results.pkl' "
                       "--exclude 'pca_results.pkl'"                       
                       )
    else:
        icasar_data = ""
        
    if exclude_epoch_images:
        epoch_images = "--exclude '0*_*.png' --exclude 'mask_status.png'"
    else:
        epoch_images = ""
        
    if exclude_epoch_data:
        epoch_data = ("--exclude 'epoch_images_data.pkl'"  
                     " --exclude 'time_course_info.pkl'")
    else:
        epoch_data = ""
        
    
    
    json_section = (
        f"{json_gz} {original_ts} {icasar_data} {epoch_images} {epoch_data}"
        )
    
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
