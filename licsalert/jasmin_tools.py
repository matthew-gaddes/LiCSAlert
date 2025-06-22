#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:35:43 2024

@author: matthew
"""
import pdb
from glob import glob
from pathlib import Path
from types import SimpleNamespace

#%%


class licsalert_volc:
    def __init__(self, name):
        """ 
        """
        self.name = name
        # self.computer_name = self.create_computer_name()
        self.processing_status = SimpleNamespace(
            name=True,
            licsbas_ts=False,
            licsalert_result=False
            )
        
        
    def create_computer_name(self):
        """
        """
        from unidecode import unidecode
        # convert to snake case and all lowercase as this is format used by comet 
        volc_name_converted = self.name.lower().replace(' ', '_')
        # swap any accents to nearest ASCII 
        volc_name_converted = unidecode(volc_name_converted)
        return volc_name_converted
    
    
    
    
    def determine_region_and_frames_from_computer_name(
            self, comet_volcano_frame_index
            ):
        """ Given a comet volcano computer_name, find the frames and region from the 
        comet volcano frame index (list of all frames split but region)
        """
        
        # print(f"    {self.name}")
        
        # stores the frames for a certain volcano
        volc_regions = []
        volc_frames = []
        
        succesful_locate=False
        
        # loop through each COMET region
        for region in list(comet_volcano_frame_index.keys()):
            
            # get all the frames for that region 
            # (i.e. multiple for each volcano)
            frames = comet_volcano_frame_index[region]
            
            # loop through the frames to see if the volcano name is in the frame name
            for frame in frames:
                # the last 18 digits of a frame are the track / burst etc.  
                # before that is the volcano name in snake_case
                if self.computer_name == frame[:-18]:
                    
                    succesful_locate=True
                    # record, regions should awlays agree
                    volc_regions.append(region)        
                    volc_frames.append(frame)
                    
                    
        if succesful_locate:
            # record all the frames for that volcano
            self.frames = volc_frames
            # check that the region awlays agrees
            unique_regions = list(set(volc_regions))
            if len(unique_regions) > 1:
                raise Exception(
                    "    Multiple regions have been found for one volcano, "
                    "exiting"
                    )
            else:
                self.region = unique_regions[0]
            self.processing_status.licsbas_ts=True
                
        else:
            print(f" Unable to find a COMET frame")
            
            
    def determine_frames_from_computer_name_and_region(
            self, comet_volcano_frame_index
            ):
        """ Given a comet volcano computer_name, find the frames and region from the 
        comet volcano frame index (list of all frames split but region)
        """
        
        print(f"    {self.name}", end='')
        
        # stores the frames for a certain volcano
        volc_frames = []

        # intialise to record if we find frames for that volcano        
        succesful_locate=False
        
        # get all the frames for that region 
        # (i.e. multiple for each volcano)
        if self.region != None:
            frames = comet_volcano_frame_index[self.region]
            print(f" | Region: {self.region}", end = '')
        
        # if there is no region, we can try and find the frames based only on
        # the name and then exit this method wth return
        else:
            print(f" | Region: currently unknown, searching by name", end='')
            self.determine_region_and_frames_from_computer_name(
                comet_volcano_frame_index
                )

            return
        
        # loop through the frames to see if the volcano name is in the frame name
        for frame in frames:
            # the last 18 digits of a frame are the track / burst etc.  
            # before that is the volcano name in snake_case
            if self.computer_name == frame[:-18]:
                succesful_locate=True
                volc_frames.append(frame)
                    
        if succesful_locate:
            # record all the frames for that volcano
            self.frames = volc_frames
            self.processing_status.licsbas_ts=True
            # print(f" | Frames: {self.frames}")
            print(f" | #frames: {len(self.frames)}")
                
        else:
            print(f" | #frames: 0    <-------- warning")


#%%



def comet_db_to_licsalert_volcs(pkl_path, priorities):
    """
    Load a volcano DataFrame from a pickle and return a nested list of names,
    one sub-list per requested priority, in the same order as `priorities`.

    Parameters
    ----------
    pkl_path : str
        Path to the .pkl file that contains the DataFrame.
        The DataFrame must have at least the columns 'name' and 'priority'.
    priorities : List[str]
        Priority classes you want (e.g. ["A1", "A2"]).

    Returns6
    -------
    List[List[str]]
        Outer list follows the order of `priorities`.
        Each inner list contains the volcano names that match that priority.
        If a priority is absent in the DataFrame you get an empty list.
    """
    import pandas as pd
    
    # 1. Load the DataFrame
    df = pd.read_pickle(pkl_path)

    # 2. Normalise priority column to *strings*; strip spaces, replace NaNs with "None"
    df = df.copy()
    df["priority"] = (
        df["priority"]
          .astype(str)            # turn None/NaN into "nan"
          .where(~df["priority"].isin(["nan", "None"]), "None")  # unify null label
          .str.strip()
    )

    # 3. Build the volcs list
    # make a new dataframe that only has volcs of the required priority.  
    df_subset = df[df["priority"].isin(priorities)]

    # initialise
    volcs = []
    for row in df_subset.itertuples(index=False):
        # initialise using custom class
        volc = licsalert_volc(row.name)
        #add some attributes
        volc.computer_name = row.vportal_name
        volc.region = row.vportal_area
        volcs.append(volc)
    
    return volcs


#%% volcano_name_to_comet_frames()

def volc_names_to_licsalert_volcs(volc_names, comet_volcano_frame_index):
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
        
    volcs = []
    # iterate over each volcano
    for volc_n, volc_name in enumerate(sorted(volc_names)):
        
        # initialise as custom object , only deal with name
        volc = comet_volcano(volc_name)

        # find the region, and any comet frames associated with that vol
        # if volc_n == 125:
        #     pdb.set_trace()
        volc.determine_region_and_frames(comet_volcano_frame_index)

        # add to list of all volcs        
        volcs.append(volc)

    return volcs



    


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
            
            


#%% get_lon_lat_of_volcs_from_db()


def get_lon_lat_of_volcs_from_db(volcs, df_source):
    """
    For every volcano object in `volcs`, look up its coordinates in a
    DataFrame and add an attribute `lon_dat = (lon, lat)`.

    Parameters
    ----------
    volcs : list
        List of user-defined volcano objects.  Each must have a `.name` attribute.
    df_source : str | pandas.DataFrame
        Path to a pickle/CSV/Parquet file **or** an already-loaded DataFrame
        containing at least the columns `name_col`, `lon_col`, `lat_col`.
    name_col, lon_col, lat_col : str
        Column names in the DataFrame.
    case_insensitive : bool, default True
        Match names ignoring case and surrounding whitespace.

    Returns
    -------
    None
        Objects are modified in place; nothing is returned.
    """
    import pandas as pd
    from typing import List, Union

    

    # 1. Load or accept the DataFrame
    if isinstance(df_source, pd.DataFrame):
        df = df_source.copy()
    else:
        # pick the loader you need (pickle / csv / parquet)
        df = pd.read_pickle(df_source)

    # 2. Build a fast name → (lon, lat) lookup dict
    norm = lambda s: str(s).strip().lower()


    # build the dict (key is name, value is lon lat tuples)
    coords_lookup = {
        norm(row['name']): (row['lon'], row['lat'])
        for _, row in df[['name', 'lon', 'lat']].iterrows()
    }


    # 3. Add lon lat to each volcano
    for v in volcs:
        key = norm(v.name)
        coords = coords_lookup.get(key)

        if coords is None:
            # You can choose to raise instead, or log:
            print(f"Warning: '{v.name}' not found in DataFrame")
            continue

        # Attach the new attribute
        setattr(v, "lon_lat", coords)





#%% get_lon_lat_of_volcs_from_ts_data()


def get_lon_lat_of_volcs_from_ts_data(volcs):
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
    

#%% update_volcs_with_data()


def update_volcs_with_data(volcs, volc_frame_dirs, volc_frame_names,
                           verbose = False, remove_all_nans = False):
    """ Given the LiCSAlert volcs list (list of comet volcano items),
    update to contain only ones that actually have data.  
    
    Inputs:
        volcs | list | one item for each volcano
        volc_frame_dirs | list | each item is path to directory of licsalert 
                                results
        volc_frame_names | list | names of frames in above list.  Must be the 
                                same length
                                
        remove_all_nans | boolean | if a volcano has only nans for the frame status
                                    i.e. there's no data, possibly remove it 
                                    from the volcs list.  
                                
    Returns:
        volcs | list | one item for each volcano, volcanoes with no data 
                        removed.  frame_status attribute now has path to data,
                        if it exists and NA if it doesn't 
                        
    History:
        2024_06_12 | MEG | Written
        2025_05_23 | MEG | Update so that all nans are not removed.  
    
    """
    
    def contains_only_na(lst):
        for item in lst:
            if item != "NA":
                return False
        return True
    
    del_indexes = []
    # iterate over all volcanoes
    for volc_n, volc in enumerate(volcs):
        
        # if there's no licsbas time series, there's can't be any licsalert
        if volc.processing_status.licsbas_ts:
            
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
            
            if contains_only_na(frame_status):
                # check if we should remove the volcano as there's no data
                if remove_all_nans:
                    del_indexes.append(volc_n)
                    # always update the user when a volcano is being removed.  
                    print(
                        f"Removing volcano {volc.name} as it does not have "
                        "any frames with LiCSAlert results.  "
                        )
            # if there is licsalert data, update the object to record this.  
            else:
                volc.processing_status.licsalert_result=True

            
    # Sort the indexes in descending order to avoid shifting issues, then delete
    # (only deletes if del_indexes exits)
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
    

#%% open_comet_frame_files()

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



#%% write_jasmin_download_shell_script

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
        script_file.write("# Note - this will only work if you can login to Jasmin\n")
        script_file.write("# which requires running a command such as:\n")
        script_file.write("# eval $(ssh-agent -s); ssh-add ~/.ssh/id_rsa_jasmin\n")
        script_file.write("# after you have setup a key pair with Jasmin.\n\n\n")

        # write a SIGINT so Ctrl+C exits the script , rather than just the 
        # current rsync command.  
        script_file.write("# Trap SIGINT (Ctrl+C)\n")
        script_file.write(
            "trap 'echo \"Script interrupted, exiting...\"; exit 1' INT\n\n"
            )

        # determine how many commands there will be
        frames_total = 0
        for volc in volcs:
            # some volcs may not have any frames.  
            if hasattr(volc, 'frames'):
                for frame in volc.frames:
                    frames_total += 1
        

        frame_counter = 0
        # iterate over each volcano
        for volc in volcs:
            if hasattr(volc, 'frames'):
                # and the frames for that volcano
                for frame in volc.frames:
                    # simple echo to say which frame we're downloading
                    script_file.write(
                        f"echo 'Attempting to rsync {frame} ({frame_counter} of {frames_total})'\n"
                        )
                    
                    # build the command for that frame        
                    command = (f"rsync -av {json_section} "
                               f"{jasmin_dir / volc.region / frame} "
                               f"{local_dir / volc.region }")
                    script_file.write(command + "\n\n")
                    
                    # update the counter
                    frame_counter += 1
