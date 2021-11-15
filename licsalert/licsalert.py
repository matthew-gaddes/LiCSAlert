# -*- coding: utf-8 -*-
"""
A selection of functions used by LiCSAlert

@author: Matthew Gaddes
"""

#%%

def LiCSAlert_batch_mode(displacement_r2, n_baseline_end, out_folder, 
                         ICASAR_settings, run_ICASAR = True, ICASAR_path = None, ic_classifying_model = None,
                         intermediate_figures = False, downsample_run = 1.0, downsample_plot = 0.5, t_recalculate = 10):
    """ A function to run the LiCSAlert algorithm on a preprocssed time series.  To run on a time series that is being 
    updated, use LiCSAlert_monitoring_mode.  
    
    Inputs:
        displacement_r2  | dict | Required: 
                                     incremental | rank 2 array | row vectors of the ifgs, incremental (ie not cumulative)
                                     mask  | rank 2 array | mask to conver the row vectors to rank 2 masked arrays.  
                                     ifg_dates | list of strings | YYYYMMDD_YYYYMMDD of interferograms.  
                                  Optional:
                                     lons | rank 2 array | lons of each pixel in the image.  Changed to rank 2 in version 2.0, from rank 1 in version 1.0  .  If supplied, ICs will be geocoded as kmz.  
                                     lats | rank 2 array | lats of each pixel in the image. Changed to rank 2 in version 2.0, from rank 1 in version 1.0
                                     dem | rank 2 array | height in metres of each pixel in the image.  If supplied, IC vs dem plots will be produced.  
                                                           
        n_baseline_end | int | the interferogram number which is the last in the baseline stage.  
        out_folder | path or string | name of folder in which to save ouputs.  
        ICASAR_settings | dict | contains all the settings for the ICASAR algorithm.  See ICASAR for details.  
        run_ICASAR | boolean | If false, the resutls from a previous run of ICASAR are used, if True it is run again (which can be time consuming)
        ICASAR_path | path or string | location of ICASAR package.  Note, important to remember the /lib/ part at teh end.  
        ic_classifying_model | path or string | location of Keras model for classification and localisation of deformation in a single ifg.  
        intermediate_figures | boolean | if True, figures for all time steps in the monitoring phase are created (which is slow).  If False, only the last figure is created.  
        downsample_run | float | data can be downsampled to speed things up
        downsample_plot | float | and a 2nd time for fast plotting.  Note this is applied to the restuls of the first downsampling, so is compound
        t_recalculate | int | rolling lines of best fit are recalcaluted every X times (nb done in number of data points, not time between them)
    Returns:
        out_folder with various items.  
    History:
        2020/09/16 | MEG | Created from various scripts.        
        2021_08_09 | MEG | add option to set r_recalculate, rather than hard coding it.  
        
    Stack overview:
        LiCSAlert_preprocessing
        ICASAR
        LiCSAlert
        LiCSAlert_figure
    """
    import numpy as np
    from pathlib import Path
    import glob
    import os
    import sys
    import pickle
    import shutil
    
    from licsalert.licsalert import LiCSAlert, LiCSAlert_figure, save_pickle, shorten_LiCSAlert_data, LiCSAlert_preprocessing
    from licsalert.downsample_ifgs import downsample_ifgs
    from licsalert.aux import LiCSAR_ifgs_to_s1_acquisitions, baselines_from_ifgnames
    
    # 0: Check some inputs and raise Exceptions if they're not suitable.  
    if ICASAR_path is None:
        raise Exception(f"LiCSAlert requires the path to the local copy of the ICASAR package, which can be downloaded from https://github.com/matthew-gaddes/ICASAR  Exiting...."  )
    else:
        sys.path.append(str(ICASAR_path))                  # location of ICASAR functions
        import icasar
        from icasar.icasar_funcs import ICASAR
    
    if 'out_folder' in ICASAR_settings:
        print(f"An 'out_folder' can't be set when running ICASAR within LiCSAlert.  Deleting this item from the 'ICASAR_settings' dictionary and continuing.  ")
        del ICASAR_settings['out_folder']
        
    for required_to_be_r2 in ['lons', 'lats']:                                                                          # check that these aren't rank 1 (should be rank 2 - ie a value for every pixel)
        if required_to_be_r2 in displacement_r2:                                                                        # though we also need to check they are 
            if len(displacement_r2[required_to_be_r2].shape) != 2:
                raise Exception(f"{required_to_be_r2} is not rank 2 (i.e. a value for every pixel) so exiting.   ")
                
    for required_array in ['incremental', 'mask', 'ifg_dates']:                                                             # some keys are always required in displacement_r2.  Check they exist.  
        if required_array not in displacement_r2.keys():                                                                    # loop through checkin.  
            raise Exception(f"{required_array} is not in displacement_r2, but is one of the required three ('incremental', 'mask', and 'ifg_dates').  Exiting.  ")
                

    # 1: Sort out the ouput folder, which depends on if ICASAR will be run.  
    out_folder = Path(out_folder)
    if os.path.exists(out_folder):                                                                                                     # the directory containing both LiCSAlert products and the 'ICASAR_outputs' directory exists
        print(f"The out_folder ({out_folder}) already exists, deleting the LiCSAlert outputs within it.   ")
        files = glob.glob(str(out_folder / 'LiCSAlert*'))                                                                              # get all the files in the folder that have LiCSALert in the name (ie not the ICASAR directory)
        for f in files:
            try:
                os.remove(f)                                                                                                           # try to delete assuming it's a file
            except:
                shutil.rmtree(f)                                                                                                       # if that fails, assume it's a directory and try to delete a different way.  
        if run_ICASAR:                                                                                                                 # if ICASAR is being run, the directory may or may not need to be deleted.  
            if ICASAR_settings['load_fastICA_results'] == True:                                                                        # if we will load the results of an existing run, the directory must remain 
                print(f"As run_ICASAR is selected but'load_fastICA_results' is also True, leaving the ICASAR_outputs folder.  ")
            else:
                print(f"As run_ICASAR is selected but'load_fastICA_results' is False, deleting the ICASAR_outputs folder.  ")           # if we won't load the results of an existing run, the directory can be deleted.  
                try:
                    shutil.rmtree(out_folder / 'ICASAR_outputs')                                                                            # delete the directory
                except:
                    print(f"Failed to remove the ICASAR folder ({str(out_folder / 'ICASAR_outputs')}) which is probably "
                          f"as it doesn't exist, so trying to continue anyway.  ")
        else:                                                                                                                           # if we're not running ICASAR, leave the directory alone.  
            print(f"As run_ICASAR is False, leaving the ICASAR_outputs folder")
    else:
        os.mkdir(out_folder)                                                                                                            # if nothing exits (i.e. first run), just make the folder                                  

    # 2: Prepare a dictionary to store information about the temporal information.  
    tbaseline_info = {'acq_dates'            : LiCSAR_ifgs_to_s1_acquisitions(displacement_r2['ifg_dates']),                    # get the unique acquisition (epoch) dates
                      'baselines_cumulative' : np.cumsum(baselines_from_ifgnames(displacement_r2['ifg_dates']))}                # and the cumulative temporal baselines (in days)
    
   
    # 3: Either run ICASAR to find latent spatial sources in baseline data, or load the results from a previous run.  
    displacement_r2 = LiCSAlert_preprocessing(displacement_r2, downsample_run, downsample_plot)                         # mean centre and downsize the data
    
    if run_ICASAR:                                                                                                      # set to True if we need to run ICASAR
        baseline_data = {'mixtures_r2' : displacement_r2['incremental'][:n_baseline_end],                               # prepare a dictionary of data for ICASAR
                         'mask'        : displacement_r2['mask'],
                         'ifg_dates'   : displacement_r2['ifg_dates'][:n_baseline_end]}                                 # ifg dates are stored in displacement_r2,    

        for optional_data in ['dem', 'lons', 'lats']:                                                                  # these are not guaranteed to be supplied to LiCSAlert, to check if they have been, loop through each one.  
            if optional_data in displacement_r2:                                                                       # if it's in the main displacement_r2 dict
                baseline_data[optional_data] = displacement_r2[optional_data]                                          # add it to the baseline data.  
        
        sources, tcs, residual, Iq, n_clusters, S_all_info, means = ICASAR(spatial_data = baseline_data, out_folder = out_folder/"ICASAR_outputs", **ICASAR_settings)       # run ICASAR to recover the latent sources from the baseline stage
        sources_downsampled, _ = downsample_ifgs(sources, displacement_r2["mask"], downsample_plot)                                                                         # downsample for recovered sources for plots
    else:
        try:
            with open(str(out_folder / "ICASAR_outputs/ICASAR_results.pkl"), 'rb') as f:
                sources = pickle.load(f)    
                tcs  = pickle.load(f)    
                source_residuals = pickle.load(f)    
                Iq_sorted = pickle.load(f)    
                n_clusters = pickle.load(f)    
            del tcs, source_residuals, Iq_sorted, n_clusters                                                                                      # these ICASAR products are not needed by LiCSAlert
            sources_downsampled, _ = downsample_ifgs(sources, displacement_r2["mask"], downsample_plot)                     # downsample the sources as this can speed up plotting
        except:
            raise Exception(f"Unable to open the results of ICASAR (which are usually stored in 'ICASAR_results') "
                            f"Try re-running and enabling ICASAR with 'run_ICASAR' set to 'True' in the 'LiCSAlert_settings' dictionary.  ")
    

    # 4: Possible use the VUDL-net-21 model to detect and locate deformation in the ICs.      
    if ic_classifying_model != None:
        print(f"Using the supplied neural network to determine which ICs are deformation and which are atmosphere.  ",
              f"This will require keras to be imported.  ")
        from LiCSAlert_neural_network_functions import sources_though_cnn                                                                        # don't import if we don't need to as this needs Keras and may not be installed
        sources_labels = {'defo_sources' : ['Dyke', 'Sill/Point', 'Atmosphere' ]}                                                                # one hot encodings to labels for Vudl-net}
        sources_labels['Y_class'], sources_labels['Y_loc'] = sources_though_cnn(sources, tcs, displacement_r2['mask'], ic_classifying_model)
    else:
        sources_labels = None
        

    # 5a: Either do LiCSAlert and the LiCSAlert figure for all time steps, 
    if intermediate_figures:                                                                                                # controls if we enter the intermediate ifgs loop.  
        for ifg_n in np.arange(n_baseline_end+1, displacement_r2["incremental"].shape[0]+1):
            
            displacement_r2_current = shorten_LiCSAlert_data(displacement_r2, n_end=ifg_n)                                  # get the ifgs available for this loop (ie one more is added each time the loop progresses)
            baselines_cumulative_current = tbaseline_info['baselines_cumulative'][:ifg_n]                                   # also get current time values
        
        
            sources_tcs_monitor, residual_monitor = LiCSAlert(sources, baselines_cumulative_current, displacement_r2_current["incremental"][:n_baseline_end],               # do LiCSAlert
                                                              displacement_r2_current["incremental"][n_baseline_end:], t_recalculate=t_recalculate, 
                                                              out_file = out_folder / 'LiCSAlert_results.pkl')    
        
            LiCSAlert_figure(sources_tcs_monitor, residual_monitor, sources, displacement_r2_current, n_baseline_end, 
                              baselines_cumulative_current, time_value_end=tbaseline_info['baselines_cumulative'][-1], out_folder = out_folder,
                              day0_date = tbaseline_info['acq_dates'][0], sources_labels = sources_labels)                                                                                 # main LiCSAlert figure, note that we use downsampled sources to speed things up

    # 5b: Or just do LiCSAlet and the LiCSAlert figure for the final time step (much quicker, but only one LiCSAlert figure is created)
    else:
        sources_tcs_monitor, residual_monitor = LiCSAlert(sources, tbaseline_info['baselines_cumulative'], displacement_r2["incremental"][:n_baseline_end],                       # Run LiCSAlert once, on the whole time series.  
                                                          displacement_r2["incremental"][n_baseline_end:], t_recalculate=t_recalculate, 
                                                          out_file = out_folder / 'LiCSAlert_results.pkl')    
        
        LiCSAlert_figure(sources_tcs_monitor, residual_monitor, sources, displacement_r2, n_baseline_end,                                                       # and only make the plot once
                          tbaseline_info['baselines_cumulative'], time_value_end=tbaseline_info['baselines_cumulative'][-1], day0_date = tbaseline_info['acq_dates'][0], 
                          out_folder = out_folder, sources_labels = sources_labels)                 
 

#%%

def LiCSAlert(sources, time_values, ifgs_baseline, ifgs_monitoring = None, t_recalculate = 10, verbose=False, out_file=None):
    """ Main LiCSAlert algorithm for a daisy-chain timeseries of interferograms.  
    
    Inputs:
        sources | r2 array | sources (from ICASAR) as row vectors, as per ICA, that can be turned back to interferograms with a rank 2 boolean mask of which pixels are masked and the col_to_ma function.  
        time_values | r1 array | time values for each point in the time series, commonly (12,24,36) for Sentinel-1 data.  Could also be described as the cumulative temporal baselines.  
        ifgs_baseline | r2 array | ifgs used in baseline stage as row vectors
        ifgs_monitoring | r2 array or None | ifgs used in monitoring stage as row vectors
        t_recalculate | int | rolling lines of best fit are recalcaluted every X times (nb done in number of data points, not time between them)
        verbose | boolean | if True, various information is printed to screen.  
        out_file | None or string | If not None, name of file that is used to save the LiCSALert outputs.  
        
    Outputs
        sources_tcs_monitor | list of dicts | list, with item for each time course.  Each dictionary contains:
        cumualtive_tc | cumulative time course : the cumulative sum of how that IC is used to fit each interferogram.  Column vector.  
        gradient      | gradient of the cumulative time course in the baseline stage.  float.  
        lines         | n_times x n_times.  The lines of best fit used to calculate the line-of-best-fit to point distance at each time.  
                        Each column is a line of best fit, but many values are zeros as the lines of best fit are redrawn periodically.  
        sigma         | sigma (standard deviation) of the line-to-point distances in the baseline stage
        distances     | number of standard deviations each point is from it's line of best fit.  
        t_recalcaulte | integer, how often the lines of best fit are redrawn.  
                                                                
        residual_tcs_monitor | list of dicts | As per above, but only for the cumulative residual (i.e. list is length 1)
    
    History:
        2019/12/XX | MEG |  Written from existing script.  
        2020/02/16 | MEG |  Update to work with no monitoring interferograms
        2020/12/14 | MEG | Improve docs and add option to save outputs.  
    """
    import numpy as np
    import pickle
    from licsalert.licsalert import bss_components_inversion, residual_for_pixels, tcs_baseline, tcs_monitoring  
        
    # -1: Check common input errors (number of pixels first, number of time second)
    if sources.shape[1] != ifgs_baseline.shape[1]:                                                                                              # 2nd dimension ([1] bit) is the number of pixels, which we compare.  
        raise Exception(f"The sources don't have the same number of pixels ({sources.shape[1]}) as the interferograms "
                        f"({ifgs_baseline.shape[1]}), so can't be used to fit them.  This is usually due to changing "
                        f"the cropped region but not re-running ICASAR.  Exiting...")
    else:
        pass
    
    if ifgs_monitoring is None:                                                                                                                     # in the case that there are no monitoring ifgs
        if time_values.shape[0] != ifgs_baseline.shape[0]:                                                                                          # check the number of times agree
            raise Exception(f"There appears to be a mismatch between the number of times ({time_values.shape[0]}, set by time_values), "            # error if they don't agree
                            f"and the number of interferograms ({ifgs_baseline.shape[0]}, in ifgs_baseline).  Exiting...")
        else:
            print(f"There are {time_values.shape[0]} times (set by time_values), which agrees with the  "                                           # or update user that all ok if they don
                  f"{ifgs_baseline.shape[0]} interferograms in ifgs_baseline.  " )
    else:
        if time_values.shape[0] != (ifgs_baseline.shape[0] + ifgs_monitoring.shape[0]):                                                             # check in the case that there are monitoring ifgs
            raise Exception(f"There appears to be a mismatch between the number of times ({time_values.shape[0]}, set by time_values), "            # error if they don't agree
                            f"and the number of interferograms ({ifgs_baseline.shape[0]}, in ifgs_baseline; and {ifgs_monitoring.shape[0]}, in ifgs_monitoring)).  Exiting..." )
        else:
            print(f"There are {time_values.shape[0]} times (set by time_values), which agress with the  "                                            # or update user that all ok if they do
                  f"{ifgs_baseline.shape[0]} interferograms in ifgs_baseline and {ifgs_monitoring.shape[0]} in ifgs_monitoring.  " )
        
    
    # 0: Ensure we can still run LiCSAlert in the case that we have no monitoring interferograms (yet)
    n_times_baseline = ifgs_baseline.shape[0]
    if ifgs_monitoring is None:
        ifgs_all = np.copy(ifgs_baseline)                                                                # if there are no monitoring ifgs, ifgs_all is just the set of baseline ifgs
        n_times_monitoring = 0                                                                           # there are no monitoring ifgs
    else:
        ifgs_all = np.vstack((ifgs_baseline, ifgs_monitoring))                                           # ifgs are row vectors, so stack vertically
        n_times_monitoring = ifgs_monitoring.shape[0]
    #print(f"LiCSAlert with {n_times_baseline} baseline interferograms and {n_times_monitoring} monitoring interferogram(s).  ")    
        
    # 1: calculating time courses/distances etc for the baseline data
    tcs_c, _ = bss_components_inversion(sources, ifgs_baseline, cumulative=True)                         # compute cumulative time courses for baseline interferograms (ie simple inversion to fit each ifg in ifgs_baseline using sources.  )
    sources_tcs = tcs_baseline(tcs_c, time_values[:n_times_baseline], t_recalculate)                     # lines, gradients, etc for time courses.  Note done by tcs_baseline, which contrasts with tcs_monitoring (which is used lower down)
    _, residual_cb = residual_for_pixels(sources, sources_tcs, ifgs_baseline)                            # get the cumulative residual for the baseline interferograms
    residual_tcs = tcs_baseline(residual_cb, time_values[:n_times_baseline], t_recalculate)              # lines, gradients. etc for residual 
    del tcs_c, residual_cb
    
    #2: Calculate time courses/distances etc for the monitoring data
    if ifgs_monitoring is not None:
        tcs_c, _ = bss_components_inversion(sources, ifgs_monitoring, cumulative=True)                      # compute cumulative time courses for monitoring interferograms (ie simple inversion to fit each ifg in ifgs_baseline using sources.  )
        sources_tcs_monitor = tcs_monitoring(tcs_c, sources_tcs, time_values)                               # update lines, gradients, etc for time courses  Note done by tcs_monitoring, which contrasts with tcs_baseline (which is used before)
    
        #3: and update the residual stuff                                                                            # which is handled slightly differently as must be recalcualted for baseline and monitoring data
        _, residual_c_bm = residual_for_pixels(sources, sources_tcs_monitor, ifgs_all)                               # get the cumulative residual for baseline and monitoring (hence _cb)    
        residual_tcs_monitor = tcs_monitoring(residual_c_bm, residual_tcs, time_values, residual=True)               # lines, gradients. etc for residual 
    

    if ifgs_monitoring is None:
        if out_file is not None:
            with open(out_file, 'wb') as f:
                pickle.dump(sources_tcs, f)
                pickle.dump(residual_tcs, f)
            f.close()
        return sources_tcs, residual_tcs
   
    else:
        if out_file is not None:
            with open(out_file, 'wb') as f:
                pickle.dump(sources_tcs_monitor, f)
                pickle.dump(residual_tcs_monitor, f)
            f.close()
        return sources_tcs_monitor, residual_tcs_monitor
    


#%%

def residual_for_pixels(sources, sources_tcs, ifgs, n_skip=None):
    """
    Given spatial sources and their time courses, reconstruct the entire time series and calcualte:
        - RMS of the residual between each reconstructed and real ifg
        - RMS of the cumulative residual between each reconstucted and real ifg
    N.b. if you only want to use some of the time courses (ie time course is 100 ifgs long, but only using the last 50 ifgs),
         this first n_skip ifgs can be set.  

    Inputs:
        sources | r2 array | sources as row vectors
        sources_tcs | list of dicts | As per LiCSAlert, a list with an item for each sources, and each item is a dictionary of various itmes
        ifgs | r2 array | interferograms as row vectors
        n_skip | None or int | if an int, the first n_skip values of the timecourses will be skipped.  

    Outputs:
        residual_ts | r2 array | Column vector of the RMS residual between that ifg, and its reconstruction
        residual_cs | r2 array | Column vector of the RMS residual between an ifg and the cumulative residual (for each pixel)
                                 N.b. the point is that if we have a strong atmosphere, it then reverses in the next ifg
                                 so the cumulative for each pixel goes back to zero

    2019/01/XX | MEG | Written, in discussion with AH
    2019/12/06 | MEG | Comment and documentation
    2020/01/02 | MEG | Update to use new LiCSAlert list of dictionaries
    2020/02/06 | MEG | Fix bug as had forgotten to convert cumulative time courses to be incremental
    """

    import numpy as np

    def list_dict_to_r2(sources_tcs):
        """ Extract the cumulative time courses from the sources_tcs list of dicts, and convert them to incremental
        time courses.  
        """
        
        n_sources = len(sources_tcs)
        n_ifgs = sources_tcs[0]["cumulative_tc"].shape[0]                               # as many rows as time steps
        tcs_r2 = np.zeros((n_ifgs, n_sources))                                          # initiate
        for n_source, source_tc in enumerate(sources_tcs):                              # loop through each source
            tc_c = source_tc["cumulative_tc"]                                           # and copy the cumulative time course out
            tc = np.diff(np.vstack((np.array([[0]]), tc_c)), axis = 0)                    # convert to incremental time course, add 0 at start to ensure that we get something of the same length.  
            tcs_r2[:,n_source:n_source+1] = tc                                          # store as column vector
        return tcs_r2
                    

    (n_sources, n_pixs) = sources.shape                                         # number of sources and number of pixels
    tcs = list_dict_to_r2(sources_tcs)                                          # get the incremental time courses as a rank 2 array
    if n_skip is not None:                                                      # crop/remove the first ifgs
        tcs = tcs[n_skip:,]                                                        # usually the baseline ifgs when used with monitoring data
    
    data_model_residual = ifgs - (tcs @ sources)                                # residual for each pixel at each time
    data_model_residual_cs = np.cumsum(data_model_residual, axis = 0)           # summing the residual for each pixel cumulatively through time   
    residual_ts = np.zeros((data_model_residual.shape[0], 1))                                   # initiate, n_ifgs x 1 array
    residual_cs = np.zeros((data_model_residual.shape[0], 1))                                   # initiate, n_ifgs x 1 array
    for row_n in range(data_model_residual.shape[0]):                                           # loop through each ifg
        residual_ts[row_n, 0] = np.sqrt(np.sum(data_model_residual[row_n,:]**2)/n_pixs)         # RMS of residual for each ifg
        residual_cs[row_n, 0] = np.sqrt(np.sum(data_model_residual_cs[row_n,:]**2)/n_pixs)      # RMS of residual for cumulative ifgs.  i.e. if atmosphere reverses, should cancel in 2nd cumulative ifg.  
    
    return residual_ts, residual_cs

#%%

def tcs_baseline(tcs_c, time_values, t_recalculate):
    """
    Given cumulative time courses (tsc_c), the time values for each entry (time_values), and the recalculation time 
    (t_recalculate), create a list with a dictionary about each time course.  The dictionary contains the cumulative 
    time courses, their gradient, and the ditsances each point is from a line of best fit at that gradient, and the 
    lines of best fit, redrwwn every t_recalculate.  
    
    Inputs:
        tcs_c | rank2 array | Cumulative Time CourseS, as column vectors
        time_values  | r1 array | time values for each point in the time course.  for Sentinel-1, commonly (12,24,36 etc)
        t_recalculate | int | the number of ifgs used to calculate the rolling lines of best fit
    
    Outputs:
        sources_tcs | list of dicts | as per description.  
        
    History:
        2020/01/02 | MEG | Written
    """
    import numpy as np
    
    n_times, n_tcs = tcs_c.shape                                                                    # there will be as many time courses (tcs) as there are sources, which are rows
    sources_tcs = []                                                                           # each time course will be a list
    for n_tc in range(n_tcs):
        tc_dict = {}                                                                           # initialise
        # 1: add the cumulative time course
        tc_dict["cumulative_tc"] = tcs_c[:, n_tc:n_tc+1]                                       # keep cumulative time course as a column vector
        
        # 2: add the gradient 
        gradient, y_intercept = np.polyfit(time_values, tc_dict["cumulative_tc"], 1)            # gradients 1st, y intercept second        
        tc_dict["gradient"] = gradient[0]                                                       # add to dictionary
        
        # 3: Lines of best fit
        line_yvals = np.polyval((tc_dict["gradient"], y_intercept), time_values)                # line of best fit, using the calcaulted y value
        lines_yvals = np.nan * np.ones((n_times, n_times))                                      # initiate as nans
        for time_step in range(n_times):                                                        # loop through each time step, to calculate it
            start_time = time_step - (t_recalculate-1)
            if start_time < 0:                                                                  # can't have negative time
                start_time = 0
            lines_yvals[start_time:time_step+1, time_step] = line_yvals[start_time:time_step+1]              # copy the y values for that little bit of line
        tc_dict["lines"] = lines_yvals                                                          # add to dictionary
        
        # 4: line to point distances (which are stored in the dict in terms of how many sigmas they are)
        line_point_distances = tc_dict["cumulative_tc"] - line_yvals[:,np.newaxis]
        tc_dict["sigma"] = np.std(line_point_distances)
        tc_dict["distances"] = np.abs(line_point_distances / tc_dict["sigma"])                  # ie the number of standard deviations a point is from the line of best fit
        
        # 5: Record the t_recalculate parameter
        tc_dict["t_recalculate"] = t_recalculate
        
        sources_tcs.append(tc_dict)                                                            # append to list, where each item is a time_course dictionary
    return sources_tcs
    
#%%    
    
def tcs_monitoring(tcs_c, sources_tcs, time_values, residual=False):
    """
    Given an extension to a time series (i.e. not in monitoring mode) and given the time courses for the baseline
    data (sources_tcs), create a new list of dicts for the complete data (i.e. rolling lines of best fit, line to point
    distances).  
    
    Inputs:
        tcs_cs | rank 2 array | time courses as column vectors
        sources_tcs | list of dicts | 
        time_values  | r1 array | time values for each point in the time course.  for Sentinel-1, commonly (12,24,36 etc)
        
    Outputs:
        sources_tcs | list of dicts | as per description.  
        
    History:
        2020/01/02 | MEG | Written
    """
    import numpy as np
    import copy                                                                                   # needed for the deep copy of the list of dictionaries
    
    # 1: Small initial steps
    sources_tcs_monitor = copy.deepcopy(sources_tcs)                                              # copy so that we can change without affecting version calcaulted only from baseline data
    t_recalculate = sources_tcs_monitor[0]["t_recalculate"]                                       # get the recalculation time used during the baseline stage
    n_times_baseline = sources_tcs_monitor[0]["cumulative_tc"].shape[0]
    n_times_tc = tcs_c.shape[0]                                                              # number of time steps is just the number of interferograms
    
    if residual:
        n_times_monitor = n_times_tc - n_times_baseline
    else:
        n_times_monitor = n_times_tc
    n_times_total = n_times_baseline + n_times_monitor       
    
    for source_n, source_tc in enumerate(sources_tcs_monitor):
        # 1: Update the cumulative time course 
        if residual is False:                                                                                   # residual timecourse doesn't start from 0, so don't need to add value to it
            tcs_c[:, source_n] += source_tc["cumulative_tc"][-1]                                                # add last value of cumulative tc, so that we continue from that value (and don't reset ot zero0)
            source_tc["cumulative_tc"] = np.vstack((source_tc["cumulative_tc"], tcs_c[:, source_n][:, np.newaxis]))        # extend the cumulative time course
        else:
            source_tc["cumulative_tc"] = tcs_c[:, source_n][:, np.newaxis]                                    # don't need to extend as have full time course
        
        # 2: Lines of best fit, and line to point distances for the monitoring data
        source_tc["lines"] = np.pad(source_tc["lines"], [(0,n_times_monitor),(0,n_times_monitor)], "constant", constant_values=(np.nan))        # resize, keeping original values in top left corner
        source_tc["distances"] = np.pad(source_tc["distances"], [(0, n_times_monitor), (0,0)], "constant", constant_values=(0))                 # lengthen to incorporate the monitoring data
        for n_ifg in np.arange(n_times_baseline, n_times_total):                                                                                # loop through each monitoring ifg
            line_y_intercept = (np.mean(source_tc["cumulative_tc"][n_ifg-t_recalculate: n_ifg]) 
                             - (source_tc["gradient"]*np.mean(time_values[n_ifg-t_recalculate: n_ifg])))                   # find the y-intercept of a the line 
            line_yvals = (time_values[n_ifg-t_recalculate: n_ifg+1] * source_tc["gradient"]) + line_y_intercept                                                                      # predict the y values given the gradient and y-intercept of the line, note that also predicting y value for next point
            source_tc["lines"][n_ifg-t_recalculate: n_ifg+1, n_ifg] = line_yvals
            source_tc["distances"][n_ifg,] = (np.abs(source_tc["cumulative_tc"][n_ifg,] - line_yvals[-1]))/source_tc["sigma"]      
    return sources_tcs_monitor
    
#%%

def LiCSAlert_figure(sources_tcs, residual, sources, displacement_r2, n_baseline_end, time_values, day0_date=None,
                     time_value_end=None, out_folder=None, ifg_xpos_scaler = 15, sources_labels = None):
    """
    The main fucntion to draw the LiCSAlert figure.  
    
    Inputs:
        sorces_tcs | list of dicts | Each source is an item in the list and has its own dictionary, containing information such as the lines of best, the graident learned
                                     in the baseline stage, the lines-of-best-fit to points distances etc.  
         residual | list of dicts | Same structure as above, but as there is only one residual, list is of length 1.  Shuld contain: cumulative timecourse 
                                    (cumualtive_tc), gradient, lines, sigma, distances, and t_recalculate.  
        sources | r2 array or None | sorces (recoverd by ICASAR) as row vectors.  N.b. must be the same size as the downsampled mask in displacement_r2
                                     If set to None, the full resolution source will be plotted        
        displacement_r2 | dict | contains ifgs as row vectors in "incremental" and their mask ("mask"), and also downsampled versions for faster figures,
                                  and their mask.  Downsampled ones used in plotting!  
        mask | r2 array | to convert an ifg (or source) as a row vector into a rank 2 masked array
        n_baseline_end | int | number of ifgs at which we switch from baseline to monitoring
        time_values | r1 array | time values to end date of each ifg.  e.g. [12,24,36] etc.  Also could be called cumulative baselines
        day0_date | string or None |  date of start of time series / first acquisition.  Used along with time values to make x tick labels as dates

        time_value_end | int or None | if an int, axes are extended to this time value, even if it's larrger than the last value of time_values
                                        N.b. time is in days, so e.g. set it to 96, or 192
        out_folder | pathlib Path | If not None, output pngs are saved to this location and the matplotib figures closed
        ifg_xpos_scaler | int | To be positioned correctly in the x direction, the ifgs that are plotted on the upper row must not be taller
                                than the axis they lie within.  Increasing this value makes the ifgs smaller, and therefore fit.  
        WIP sources_labels | dict |  keys and variables:
                                    'defo_sources' : ['dyke', 'sill', 'no_def' ]
                                    'Y_class' : n_sources x 3, labels as one hot encoding.  
                                    'Y_loc' : n_sources x 4, location of deformation.  [0,0,0,0] if none present.  
        
    Returns:
        figure
        
    History:
        2020/01/XX | MEG | Written
        2020/01/10 | MEG | update to add "upper_time_values"
        2020/02/16 | MEG | add ifg_xpos_scaler to make sure ifgs are plotted with the correct x value.  
        2020/03/08 | MEG | Change plotting of ifgs and sources to awlays be the downsampled ones.  
        2020/04/20 | MEG | Update so that x tick labels are dates and not numbers since time series started.  
        2020/06/23 | MEG | Write documentation for dates argument.  
        2020/12/15 | MEG | Determine whether sources are downsampled automatically, and also raise exception if the number of pixels doesn't agree.  
        2021_09_28 | MEG | Update various plotting featuers (add cumulative ifgs, add DEM with lons and lats.  )
        2021_09_29 | MEG | Add ifg date
    
    """
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib import ticker
    import matplotlib as mpl
    from matplotlib.ticker import MultipleLocator
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes            
    plt.switch_backend('Agg')                                                           #  works when there is no X11 forwarding, and when displaying plots during creation would be annoying.  
    #plt.switch_backend('Qt5Agg')                                                           #  for debugging
    import datetime as dt 
    
    from licsalert.aux import col_to_ma

    
    def calcualte_line_args(n_ifgs, t_recalculate):
        """Lines of best fit are calculated for eahc time step, but we don't want
        to plot them all (as they lie on top of each other, mostly).  
        Therefore, calcaulte the numbers of which to plot so that they don't overlap.  '
        """
        line_args = []    
        for k in range(n_ifgs):
            if k % t_recalculate == 0 and k != 0:                                   # pick which ones, but have to exclue the 0th one
                line_args.append(k-1)
            if k == (n_ifgs-1) and k not in line_args:                              # and always pick the last one
                line_args.append(k)
        return line_args
    

    def plot_ifgs(ifgs, pixel_mask, figure, gridspec_area, time_values, xlim, ylabel = '', day0_date = None, cumulative = False):
        """ Plot all the ifgs (baseline and monitoring) within the grispec_area.  
        """
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes            
        # 1: Create a single wide axes for the inset axes to be plotted on
        ax_ifgs = plt.Subplot(figure, gridspec_area)                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
        fig1.add_subplot(ax_ifgs)                                                                   # add to figure
        ax_ifgs.set_yticks([])                                                                      # no y ticks
        ax_ifgs.set_ylim(bottom = 0, top = 1)
        ax_ifgs.set_xlim(left = 0, right = xlim)                                                    # set x axis upper limit to be the number of acquisitions in the time series
        if day0_date is not None:       
            xticks_every_3months(ax_ifgs, day0_date, time_values, include_tick_labels = False)      # update the ticks (but not labels) to be the same as the time course and residual axis

        for ifg_n, source in enumerate(ifgs):                                                                                # ifgs are rows, loop through
            iax = ax_ifgs.inset_axes([time_values[ifg_n], 0., (xlim/ifg_xpos_scaler), 1.], transform=ax_ifgs.transData)      # [xpos, ypos, xwidth, ywidth], note the x_pos_scaler that reduces the size of the inset axes to make sure it remains in tehe right place
            ifg_plot = iax.imshow(col_to_ma(source, pixel_mask), cmap = plt.get_cmap('coolwarm'))                                       # plot on the axes
            iax.set_xticks([])                                                                                               # images so get rid of x and y ticks (nb this is just for the inset axes)
            iax.set_yticks([])
            # print(time_values[ifg_n])                                                                       # for debugging
            # plt.pause(1)                                                                                    # "
            if ifg_n == (ifgs.shape[0] -1):                                                                                             # if it's the last interferogram
                cbar_ax = inset_axes(iax, width="7%", height="40%",   loc='lower left',  bbox_to_anchor=(1.05, 0.1, 1, 1),              # Colorbar: isnet axes just to left of the main axis
                                     bbox_transform=iax.transAxes,borderpad=0)
                
                #cbar = fig1.colorbar(ifg_plot, cax = cbar_ax, ticks = [np.nanmin(source), np.nanmax(source)])                            # colorbar, tick only 0 and the max (and check max is not a nan)
                cbar = fig1.colorbar(ifg_plot, cax = cbar_ax)                           
                cbar.set_ticks([np.nanmin(source), np.nanmax(source)])
                cbar.set_ticklabels([f"{np.nanmin(source):.3} m", f"{np.nanmax(source):.3} m"])
                cbar_ax.tick_params(labelsize=6)                                                                                         #                                    
                if day0_date is None:                                                                                                   # if we don't have date information,
                    pass                                                                                                                # we can't add it via the colorbar title
                    #cbar_ax.set_title('LOS disp. (m)', fontsize = 6, loc = 'left')
                else:                                                                                                                    # if we do ahve date information
                    day0_date_dt = dt.datetime.strptime(day0_date, "%Y%m%d")                                                             # label the ifg with date interval, start by converting day0 date
                    if cumulative:                                                                                                       # cumulative ifg spans from start to end
                        ifg_start_date = day0_date_dt
                    else:
                        ifg_start_date = day0_date_dt + dt.timedelta(int(time_values[-2]))                                               # but incremental is from previous date
                    ifg_end_date = day0_date_dt + dt.timedelta(int(time_values[-1]))                                                     # to final date
                    cbar_ax.set_title(f"{dt.datetime.strftime(ifg_start_date, '%Y%m%d')}\n"
                                      f"{dt.datetime.strftime(ifg_end_date, '%Y%m%d')}", fontsize = 6, loc = 'left')
        ax_ifgs.set_ylabel(ylabel, fontsize = 7, labelpad = -1)
        

    def colourbar_for_sources(icasar_sources):
        """ Creat a colourbar for the ICA sources that is centered on 0, and cropped so that each side is equal
        (i.e. if data lies in range [-1 10], will only go slightly blue, but up to max red, with grey at 0)
        """
        
        ics_min = np.min(icasar_sources)                                                       # 
        ics_max = np.max(icasar_sources)
        ic_colours = plt.get_cmap('coolwarm')
        cmap_mid = 1 - ics_max/(ics_max + abs(ics_min))                                     # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
        ic_colours_cent = remappedColorMap(ic_colours, start=0.0, midpoint=cmap_mid, stop=1, name='ic_colours_cent')                    # make the colours for plotting the ICs
        return ic_colours_cent
    
    def sigma_bar_plotter(ax_tc, xvals, yvals, sigma_cmap):
        """
        Add the bar plots onto an axes to show how many sigmas from the line each point is.  
        Inputs:
            ax_tcs | axes object | axes on which to plot
            xvals | x values of bars - usually time
            yvals | height of bars - the number of sigmas from the mean that point is
        """
        ax_tc2 = ax_tc.twinx()                                                       # instantiate a second axes that shares the same x-axis
        for point_n, xval in enumerate(xvals):                                       # loop along each bar
            ax_tc2.bar(xval, yvals[point_n], width=10, alpha = 0.3, color = sigma_cmap(yvals[point_n]/5))   # and plot in required colour
        ax_tc2.set_yticklabels([])                                                      # turn off y labels
        ax_tc2.set_ylim(top = 10)                                                       # set so in range 0 to 10

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        import matplotlib.colors as colors
        new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
        return new_cmap 
    
    def xticks_every_3months(ax_to_update, day0_date, time_values, include_tick_labels):
        """Given an axes, update the xtcisk so the major ones are the 1st of jan/april/july/october, and the minor ones are the first of the 
        other months.  
        Inputs:
            ax_to_update | matplotlib axes | the axes to update.  
            day0_date | string | in form yyyymmdd
            time_values | rank 1 array | cumulative temporal baselines, e.g. np.array([6,18, 30, 36, 48])
        Returns:
            updates axes
        History:
            2021_09_27 | MEG | Written
        """
        
        from dateutil.relativedelta import relativedelta                                                    # add 3 months and check not after end
        from matplotlib.ticker import AutoMinorLocator      
        
        xtick_label_angle = 315
        
        tick_labels_days = ax_to_update.get_xticks().tolist()                                                # get the current tick labels
        day0_date_dt = dt.datetime.strptime(day0_date, "%Y%m%d")                                            
        dayend_date_dt = day0_date_dt +  dt.timedelta(int(time_values[-1]))                                 # the last time value is the number of days we have, so add this to day0 to get the end.  
    
        # 1: find first tick date (the first of the jan/ april/jul /oct)                        
        date_tick0 = day0_date_dt                                                                           
        while not ( (date_tick0.day) == 1 and (date_tick0.month == 1  or date_tick0.month == 4 or date_tick0.month == 7 or date_tick0.month == 10 )):
            date_tick0 +=  dt.timedelta(1)
            
        # 2: get all the other first of the quarters
        ticks = {'datetimes' : [date_tick0],
                 'yyyymmdd'   : [],
                 'n_day'     : []}
       
        while ticks['datetimes'][-1] < (dayend_date_dt - relativedelta(months=+3)):                         # subtract 3 months to make sure we don't go one 3 month jump too far. 
            ticks['datetimes'].append(ticks['datetimes'][-1] + relativedelta(months=+3))
        
        # 3: work out what day number each first of the quarter is.  
        for tick_dt in ticks['datetimes']:                                                                   # find the day nubmers from this.             
            ticks['yyyymmdd'].append(dt.datetime.strftime(tick_dt, "%Y/%m/%d"))
            ticks['n_day'].append((tick_dt - day0_date_dt).days)
            
        # 4: Update the figure.  
        ax_to_update.set_xticks(ticks['n_day'])                                                                   # apply major tick labels to the figure
        minor_locator = AutoMinorLocator(3)                                                                       # there are three months in each quarter, so a minor tick every month
        ax_to_update.xaxis.set_minor_locator(minor_locator)                                                       # add to figure.  
        if include_tick_labels:
            ax_to_update.set_xticklabels(ticks['yyyymmdd'], rotation = xtick_label_angle, ha = 'left')            # update tick labels, and rotate
            plt.subplots_adjust(bottom=0.15)
            ax_to_update.set_xlabel('Date')
        else:
            ax_to_update.set_xticklabels([])                                                                    # remove any tick lables if they aren't to be used.  
        
        # add vertical lines every year.  
        for major_tick_n, datetime_majortick in enumerate(ticks['datetimes']):
            if datetime_majortick.month == 1:
                ax_to_update.axvline(x = ticks['n_day'][major_tick_n], color='k', alpha=0.1, linestyle='--')                          
                   





    # -1: Check that the sizes of the sources and the interferograms agree.  Raise error if not.  
    if sources.shape[1] == displacement_r2['incremental'].shape[1]:
        sources_downsampled = False
    elif sources.shape[1] == displacement_r2['incremental_downsampled'].shape[1]:
        sources_downsampled = True
    else:
        raise Exception(f"There appears to be a mismatch in the number of pixels contained within the sources ({sources.shape[1]} pixels) "
                        f"and the interferograms ({displacement_r2['incremental'].shape[1]} pixels) or the downsampled interferograms "
                        f"({displacement_r2['incremental_downsampled'].shape[1]} pixels).  The sources must have the same number of pixels "
                        f"as one of these so that their mask can be used to turn the sources from row vectors into images.  ")
        
    
    # 0: Start, some definitions that shouldn't need changing (ie hard coded variables)
    #line_best_fit_alpha = 0.7
    dot_marker_size = 12
    if time_value_end is None:
        t_end = time_values[-1]                                                         # last time value - i.e. the right hand x value of all the axes
    else:
        t_end = time_value_end                                                          # or if it was provided, just use that value
    
    
    # 1 set some preliminary stuff
    t_recalculate = sources_tcs[0]["t_recalculate"]
    n_ics = len(sources_tcs)
    n_times = time_values.shape[0] 
    n_ifgs = displacement_r2["incremental"].shape[0]
    line_args= calcualte_line_args(n_times, t_recalculate)                      # which lines of best fit to plot
    c = mpl.colors.ColorConverter().to_rgb                                      
    cmap_discrete = make_colormap(  [c('black'), c('orange'), 0.33, c('orange'), c('yellow'), 0.66, c('yellow'), c('red')])     # custom colorbar for number of sigmas from line
    cmap_sources = colourbar_for_sources(sources)
    figtitle = f'LiCSAlert figure with {(n_ifgs-n_baseline_end):03d} monitoring interferograms'

    # 2 Initiate the figure    
    fig1 = plt.figure(figsize=(14,8))
    fig1.canvas.set_window_title(figtitle)
    grid = gridspec.GridSpec((n_ics + 3), 11, wspace=0.3, hspace=0.1)                        # divide into 2 sections, 1/5 for ifgs and 4/5 for components

    # 3: Plot the ifgs along the top
    plot_ifgs(np.cumsum(displacement_r2["incremental_downsampled"], axis = 0), displacement_r2["mask_downsampled"], fig1, grid[0,1:], 
              time_values, t_end, ylabel = 'Cumulative', day0_date = day0_date, cumulative = True)                                                                # cumulative ifgs
    plot_ifgs(displacement_r2["incremental_downsampled"], displacement_r2["mask_downsampled"], fig1, grid[1,1:], 
              time_values, t_end, ylabel = ' Incremental', day0_date = day0_date, cumulative = False)                                                               # incremental ifgs

    # 4: Plot each source and its time course 
    try:
        baseline_monitor_change = np.mean([time_values[n_baseline_end-1], time_values[n_baseline_end]])                                     # Vertical line will be drawn at this time value to show that we switch from baseline to monitoring
    except:
        baseline_monitor_change = np.mean([time_values[n_baseline_end-1], time_values[n_baseline_end-1] + 12])                              # But the above won't work if there are no monitoring ifgs, so just guess next ifg will be after 12 days and draw line as if that were true (ie 6 days after last point)
    for row_n, source_tc in enumerate(sources_tcs):
        # 4a: Plot the source
        ax_source = plt.Subplot(fig1, grid[row_n+2,0])                                                                                      # create an axes for the IC (spatial source)
        if sources_downsampled:
            im = ax_source.imshow(col_to_ma(sources[row_n], displacement_r2["mask_downsampled"]), cmap = cmap_sources, vmin = np.min(sources), vmax = np.max(sources))   # plot the downsampled source
        else:
            im = ax_source.imshow(col_to_ma(sources[row_n], displacement_r2["mask"]), cmap = cmap_sources, vmin = np.min(sources), vmax = np.max(sources))                # or plot the full resolution source
        ax_source.set_xticks([])
        ax_source.set_yticks([])
        ax_source.set_ylabel(f"IC {row_n}")
        
        # 4b: Possible annotate it with VUDL-net-21 output
        ######################## start WIP
        if sources_labels != None:
            if row_n == 0:
                ax_source.set_title('VUDL-Net-21\nprediction', fontsize = 8, color = 'tab:orange')
            from LiCSAlert_neural_network_functions import centre_to_box, add_square_plot
            #start_stop_locs_pred = centre_to_box(locs_predicted[plot_args[n_plot]])                                         # covert from centre width notation to start stop notation, # [x_start, x_stop, Y_start, Y_stop]
            start_stop_locs = centre_to_box(sources_labels['Y_loc'][row_n])                                         # covert from centre width notation to start stop notation, # [x_start, x_stop, Y_start, Y_stop]
            add_square_plot(start_stop_locs[0], start_stop_locs[1], 
                            start_stop_locs[2], start_stop_locs[3], ax_source, colour='tab:orange')                           # box around deformation
        
            ax_orange = ax_source.twinx()
            ax_orange.set_ylabel(f"{sources_labels['defo_sources'][np.argmax(sources_labels['Y_class'][row_n,])]}", fontsize = 8, color = 'tab:orange')
            #ax_orange.yaxis.label.set_color("tab:orange")
            ax_orange.set_yticks([])
            
        fig1.add_subplot(ax_source)
        ######################## end WIP
        
        
        # 4c: plot the time courses for that IC, and the rolling lines of best fit
        ax_tc = plt.Subplot(fig1, grid[row_n+2,1:])
        ax_tc.scatter(time_values, source_tc["cumulative_tc"], c = source_tc["distances"], marker='o', s = dot_marker_size, cmap = cmap_discrete, vmin = 0, vmax = 5, )                        # 
        for line_arg in line_args:                                                                                          # line args sets which lines of best fit to plot (there is a line of best fit for each point, but it's too busy if we plot them all)
            ax_tc.plot(time_values, source_tc["lines"][:,line_arg], c = 'k')                                                # ie each column is a line of best fit
    
        # 4d: tidy up some stuff on the axes
        ax_tc.axhline(y=0, color='k', alpha=0.3)  
        ax_tc.axvline(x = baseline_monitor_change, color='k', alpha=0.3)                          #line the splits between baseline and monitoring ifgs
        ax_tc.set_xlim(left = 0, right = t_end)
        if day0_date is not None:
            xticks_every_3months(ax_tc, day0_date, time_values, include_tick_labels = False)
        fig1.add_subplot(ax_tc)
        sigma_bar_plotter(ax_tc, time_values, source_tc["distances"], cmap_discrete)                # draw the bar graph showing sigma values
        ax_tc.yaxis.tick_right()                                                                    # has to be called after sigma_bar_plotter
        
                                                                
    # 5: Plot the residual
    ax_residual = plt.Subplot(fig1, grid[-1,1:])                                                                    # plot on the last row
    ax_residual.scatter(time_values, residual[0]["cumulative_tc"], marker='o', s = dot_marker_size, cmap = cmap_discrete, vmin = 0, vmax = 5, c = residual[0]["distances"])         # 
    for line_arg in line_args:                                                                                      # plot the rolling line of best fit, but not all of them (only those in line_args)
        ax_residual.plot(time_values, residual[0]["lines"][:,line_arg], c = 'k')                                    # each column is a line of best fit
    ax_residual.axhline(y=0, color='k', alpha=0.3)
    ax_residual.axvline(x = baseline_monitor_change, color='k', alpha=0.3)                          #line the splits between baseline and monitoring ifgs
    ax_residual.set_xlim(left = 0, right = t_end)                    # and finaly tidy up axis and labels etc.  
    ax_residual.yaxis.tick_right()
    ax_residual.yaxis.set_label_position("right")
    ax_residual.set_ylabel('RMS\nresidual')
    fig1.add_subplot(ax_residual)
    sigma_bar_plotter(ax_residual, time_values, residual[0]["distances"], cmap_discrete)                    # draw the bar graph showing sigma values
    ax_residual.yaxis.tick_right()                                                                        # has to be called after sigma_bar_plotter
    if day0_date is not None:
        xticks_every_3months(ax_residual, day0_date, time_values, include_tick_labels = True)                   # create the ticks and labels on the 1st of the quater.  
    
    ## 6: add the two colorbars
    cax = fig1.add_axes([0.12, 0.08, 0.005, 0.1])                                      # source strength
    ics_cbar = fig1.colorbar(im, cax=cax, orientation='vertical')
    tick_locator = ticker.MaxNLocator(nbins=4)
    ics_cbar.locator = tick_locator
    ics_cbar.update_ticks()
    ics_cbar.set_label('IC (m)')
    ics_cbar.ax.yaxis.set_label_position('left')
    
    cax2 = fig1.add_axes([0.17, 0.08, 0.005, 0.1])                                            # number of sigmas from the mean
    norm = mpl.colors.Normalize(vmin=0, vmax=5)
    std_cbar = mpl.colorbar.ColorbarBase(cax2, cmap=cmap_discrete, norm=norm,orientation='vertical')
    tick_locator2 = ticker.MaxNLocator(nbins=5)
    std_cbar.locator = tick_locator2
    std_cbar.update_ticks()
    std_cbar.set_label(r'$\sigma$ from trend line')
    std_cbar.ax.yaxis.set_label_position('left')
    
    # 7: Possibly add the DEM
    if 'dem' in displacement_r2.keys():                                                                                         # DEM is not alway included.  
        ax_dem = plt.Subplot(fig1, grid[1,0])                                                                                   # create an axes for the IC (spatial source)
        terrain_cmap = plt.get_cmap('terrain')                                                                                  # appropriate colours for a dem
        terrain_cmap = truncate_colormap(terrain_cmap, 0.2, 1)                                                                  # but crop (truncate) the blue parts as we are only interested in land
        dem_plot = ax_dem.imshow(displacement_r2["dem"], cmap = terrain_cmap)                                                   # plot the DEM
        ax_dem.xaxis.tick_top()                                                                                                 #
        ax_dem.tick_params(axis='both', which='major', labelsize=7)                                                             # adjust fontsize
        ax_dem.tick_params(axis='both', which='minor', labelsize=7)
        ax_dem.set_xticks([0, displacement_r2['dem'].shape[1]])                                                                 # tick only the min and max in each direction
        ax_dem.set_yticks([0, displacement_r2['dem'].shape[0]])

        ax_dem.xaxis.set_label_position('top')
        if ('lons' in displacement_r2.keys()) and ('lats' in displacement_r2.keys()):                                           # if we have lons and lats, we can update the tick lables to be lons and lats.  
            
            ax_dem.set_xticklabels([str(round(displacement_r2['lons'][-1,0], 2)) + "$^\circ$", str(round(displacement_r2['lons'][-1,-1], 2)) + "$^\circ$" ])
            ax_dem.set_yticklabels([str(round(displacement_r2['lats'][0,0], 2))  + "$^\circ$", str(round(displacement_r2['lats'][-1,0], 2)) + "$^\circ$"])
        # ax_dem.set_ylabel('Latitude ($^\circ$)', fontsize = 6)
        # ax_dem.set_xlabel("Longitude ($^\circ$)", fontsize = 6)
            
        #colorbar for the DEM, just gets in the way.              
        # axins = inset_axes(ax_dem, width="7%", height="50%",   loc='lower left',  bbox_to_anchor=(1.05, 0., 1, 1),              # isnet axes just to left of the main axix for a colorbar
        #                     bbox_transform=ax_dem.transAxes,borderpad=0)
        # #fig1.colorbar(dem_plot, cax = axins, ticks = [0, np.nanmax(displacement_r2['dem'])])                                    # colorbar, tick only 0 and the max (and check max is not a nan)
        # fig1.colorbar(dem_plot, cax = axins, ticks = [])                                    # colorbar, tick only 0 and the max (and check max is not a nan)
        # #axins.tick_params(axis='both', which='major', labelsize=6, rotation = 90)                                               #
        fig1.add_subplot(ax_dem)
        
        # work out the size of the ICs/ DEM and add to the DEM bit of the figure.  
        from geopy import distance
        image_size = {}
        image_size['x'] = int(distance.distance((displacement_r2['lats'][-1,0], displacement_r2['lons'][-1,0]),                       # bottom left corner  
                                                (displacement_r2['lats'][-1,-1], displacement_r2['lons'][-1,-1])).meters / 1000)      #  to bottom right, and convert to integere kms
        image_size['y'] = int(distance.distance((displacement_r2['lats'][-1,0], displacement_r2['lons'][-1,0]),                       # bottom left 
                                            (displacement_r2['lats'][0,0], displacement_r2['lons'][0,0])).meters / 1000)              # to to top left, and conver to integer kms
        
        ax_dem.text(-0.5 * displacement_r2['dem'].shape[1], -0.75 * displacement_r2['dem'].shape[0], f"WxH (km): {image_size['x']} x {image_size['y']}\n"              # add these in these labels in the space above the DEM.  
                        f"DEM (m): {int(np.nanmin(displacement_r2['dem'])), int(np.nanmax(displacement_r2['dem']))}", fontsize = 6 )
    
    # 8: Possible save output
    if out_folder is not None:
        filename = "_".join(figtitle.split(" "))                                            # figtitle has spaces, but filename must use underscores instead.  
        fig1.savefig(out_folder / f"{filename}.png", bbox_inches='tight')
        plt.close(fig1)



#%%
        
        
def LiCSBAS_for_LiCSAlert(LiCSAR_frame, LiCSAR_frames_dir, LiCSBAS_out_dir, logfile_dir, LiCSBAS_bin, lon_lat = None, downsampling = 1, n_para=1):
    """ Call this to either create a LiCSBAS timeseries from LiCSAR products, or to update one when new products become available.  
    Not all LiCSBAS features are supported! 
    
    Inputs:
        LiCSAR_frame | string | name of frame being processed.  
        LiCSAR_frames_dir | string | The path to the LiCSAR frame which contains that volcano.  Needs trailing /
        LiCSBAS_out_dir | string | path to where LiCSBAS products will be stored.  Needs trailing /
        LiCSBAS_bin | string | The LiCSBAS functions must have been added to your shell's path.  This is used to check this has been done correctly.  
        logfile | string | path to directory where logfile will be appended to. Needs trailing /
        lon_lat | list | west east south north to be clipped to, or None.  
        downsampling | int | >=1, sets the downsampling used in LiCSBAS (mulitlooking in both range and azimuth?)
        
    Returns:
        All products described in the LiCSBAS documentation.  
        Of importance for use with LiCSAlert are:
            cum.h5
            
    History:
        2020/02/15 | MEG | Written
        2020/06/22 | MEG | Add option to set processing directory (LiCSAR_dir)
        2020/06/24 | MEG | Add option to call step_05 and clip to geographic region.  
        2020/06/30 | MEG | Simplify inputs
        2020/07/02 | MEG | Add logfile_dir, and change from os.system to subprocess.call so output can be appended to a logfile (and still be displayed to a terminal)
        2020/11/11 | RR | Add n_para argument for new version of LiCSBAS
        2020/11/13 | MEG | Add LiCSBAS_bin argument to check that path is set correctly.  
        
    """

    # import sys
    import os
    import subprocess
    import sys

    # Get the user's PATH:
    user_path = os.environ['PATH'].split(':') 
    if LiCSBAS_bin.rstrip('/') not in user_path:                                                  # check if already on path
        raise Exception(f"Error - the LiCSBAS scripts don't appear to be on your path.  As these functions are called from the command line, "
                        f"the path can't be updated from within Python.  This can usually be rectified by adding a line such as this to your ~/.bashrc file: "
                        f"source <your_LiCSBAS_path>/LiCSBAS/bashrc_LiCSBAS.sh \n The LiCSBAS documentation may also be useful: "
                        f"https://github.com/yumorishita/LiCSBAS/wiki/1_Installation Exiting.  ")
        
    # Inputs args - probably a better way to change these (rather than hard-coding)    
    p11_unw_thre = 0.5
    p11_coh_thre = 0.1
    p12_loop_thre = 1.5                 # in rads
    
    #Rarely changed
    p13_inv_alg = "LS"              	# LS (default) or WLS
    p13_mem_size = 4000	                # default: 4000 (MB)
    p13_gamma = 0.0001              	# default: 0.0001
    p13_n_unw_r_thre = 1	            # default: 1
    p13_keep_incfile = "n"	            # y/n. default: n

    # make directory names in the style used by LiCSBAS.  
    GEOCdir = f"{LiCSAR_frames_dir}{LiCSAR_frame}/GEOC"                      # GEOC dir, where LiCSAR ifgs are stored
    GEOCmldir = f"{LiCSBAS_out_dir}GEOCml{downsampling}"                     # multilooked directory, where LiCSBAS products are stored
    TSdir = f"{LiCSBAS_out_dir}TS_GEOCmldir"                                 # time series directory, where LiCSBAS products are stored
    GEOCmldirclip = f"{LiCSBAS_out_dir}GEOCmldirclip"                        # clipped products, produced by step_05
       
    # Convert format (LiCSBAS02)  NB: This will automatically skip files that have already been converted.  
    subprocess.call(f"LiCSBAS02_ml_prep.py -i {GEOCdir} -o {GEOCmldir} -n {downsampling} --n_para {n_para}" + f" >&1 | tee -a {logfile_dir}LiCSBAS_log.txt", shell=True)                 # This creates the files in GEOCmlXXX, including the png preview of unw, note that 1 is stdout, -a to append

    # LiCSBAS03 - GACOS
    # LiCSBAS04 - mask    

    # LiCSBAS05 - clip to region of interest (using lat and long, but can also use pixels)
    if lon_lat is not None:
        LiCSBAS_lon_lat_string = f"{lon_lat[0]}/{lon_lat[1]}/{lon_lat[2]}/{lon_lat[3]}"                                     # conver to a string which includes / (and so python does not see them as four numbers divided!)
        subprocess.call(f"LiCSBAS05op_clip_unw.py -i {GEOCmldir} -o {GEOCmldirclip} -g {LiCSBAS_lon_lat_string} --n_para {n_para}" + f" >&1 | tee -a {logfile_dir}LiCSBAS_log.txt", shell=True)                 # # do the clipping, -g of form west/east/south/north   N.b.!  As above, careful with / being treated as divide by Python!
        GEOCmldir = GEOCmldirclip                                                                                           # update so now using the clipped products
    
    # LiCSBAS11 - check unwrapping, based on coherence
    subprocess.call(f"LiCSBAS11_check_unw.py -d {GEOCmldir} -t {TSdir} -c {p11_coh_thre} -u {p11_unw_thre}" + f" >&1 | tee -a {logfile_dir}LiCSBAS_log.txt", shell=True)   

    # LiCSBAS12 - check unwrapping, based on loop closure
    subprocess.call(f"LiCSBAS12_loop_closure.py -d {GEOCmldir} -t {TSdir} -l {p12_loop_thre} --n_para {n_para}" + f" >&1 | tee -a {logfile_dir}LiCSBAS_log.txt", shell=True)   

    # LiCSBAS13 - SB inversion
    subprocess.call(f"LiCSBAS13_sb_inv.py -d {GEOCmldir} -t {TSdir} --inv_alg {p13_inv_alg} --mem_size {p13_mem_size} --gamma {p13_gamma} --n_para {n_para} --n_unw_r_thre {p13_n_unw_r_thre} --keep_incfile {p13_keep_incfile} " + f" >&1 | tee -a {logfile_dir}LiCSBAS_log.txt", shell=True)   

    # LiCSBAS 14 - velocity standard dev
    # LiCSBAS 15 - mask using noise indicies
    # LiCSBAS 16 - fiter
    


   

#%%
    
def LiCSAlert_preprocessing(displacement_r2, downsample_run=1.0, downsample_plot=0.5, verbose=True, mean_centre = True):
    """A function to downsample the data at two scales (one for general working [ie to speed things up], and one 
    for faster plotting.  )  Also, data are mean centered, which is required for ICASAR and LiCSAlert.  
    Note that the downsamples are applied consecutively, so are compound (e.g. if both are 0.5, 
    the plotted data will be at 0.25 the resolution of the original data).  
    
    Inputs:
        displacement_r2 | dict | input data stored in a dict as row vectors with a mask    
                                 Also lons and lats as rank 2 arrays.  
        downsample_run | float | in range [0 1], and used to downsample the "incremental" data
        downsample_plot | float | in range [0 1] and used to downsample the data again for the "incremental_downsample" data
        verbose | boolean | if True, informatin returned to terminal
        mean_centre | boolean | add option to control if mean centered.  
        
    Outputs:
        displacement_r2 | dict | input data stored in a dict as row vectors with a mask
                                 updated so that "incremental" is downsampled (and its mask), 
                                 and a new key is created, called "incremental_downsampled" 
                                 that is downsamled further for fast plotting                                 
    History:
        2020/01/13 | MEG | Written
        2020/12/15 | MEG | Update to also downsample the lons and lats in the ICASAR geocoding information.  
        2021_04_14 | MEG | Update so handle rank2 arrays of lons and lats properly.  
        2021_05_05 | MEG | Add check that lats are always the right way up, and fix bug in lons.  
        2021_10_13 | MEG | Add function to also downsample ENU grids.  
        2021_11_15 | MEG | Add option to control mean centering, and warning that it is happening.  
    """
    import numpy as np
    from licsalert.downsample_ifgs import downsample_ifgs
    from licsalert.aux import col_to_ma
    from skimage.transform import rescale

    
    n_pixs_start = displacement_r2["incremental"].shape[1]                                          # as ifgs are row vectors, this is the number of pixels we start with
    shape_start = displacement_r2["mask"].shape                                                     # and the shape of the ifgs (ny, nx)
    
    # 0: Mean centre the interferograms (ie. break any connection to a reference pixel/region that the interferogram was set to)
    if mean_centre:
        if verbose:
            print(f"LiCSAlert_preprocessing: mean centering the interferograms before downsampling")
        displacement_r2["incremental"] = displacement_r2["incremental"] - np.mean(displacement_r2["incremental"], axis = 1)[:,np.newaxis]                            # mean centre the data (along rows) 

    # 1: Downsample the ifgs for use in all following functions.  
    if downsample_run != 1.0:                                                                                       # if we're not actually downsampling, skip for speed
        displacement_r2["incremental"], displacement_r2["mask"] = downsample_ifgs(displacement_r2["incremental"], displacement_r2["mask"],
                                                                                  downsample_run, verbose = False)

    # 2: Downsample for plotting
    displacement_r2["incremental_downsampled"], displacement_r2["mask_downsampled"] = downsample_ifgs(displacement_r2["incremental"], displacement_r2["mask"],
                                                                                                      downsample_plot, verbose = False)


    # 3: Downsample the geocode info (if provided) in the same was as step 1:
    if ('lons' in displacement_r2) and ('lats' in displacement_r2):                                             # check if we have lon lat data as not alway strictly necessary.  
        ifg1 = col_to_ma(displacement_r2['incremental'][0,:], pixel_mask = displacement_r2['mask'])             # get the size of a new ifg (ie convert the row vector to be a rank 2 masked array.  )
        lons = np.linspace(displacement_r2['lons'][0,0], displacement_r2['lons'][-1,-1], ifg1.shape[1])          # remake the correct number of lons (to suit the number of pixels in the new ifgs, first as 1d
        displacement_r2['lons'] = np.repeat(lons[np.newaxis, :], ifg1.shape[0], axis = 0)                       # then as 2D
        lats = np.linspace(displacement_r2['lats'][0,0], displacement_r2['lats'][-1,0], ifg1.shape[0])          # remake the correct number of lats (to suit the number of pixels in the new ifgs, first as 1d
        displacement_r2['lats'] = np.repeat(lats[:, np.newaxis], ifg1.shape[1], axis = 1)                       # then as 2D
        if displacement_r2['lats'][0,0] < displacement_r2['lats'][-1,0]:                                        # if the lat in the top row is less than the lat in the bottom row
            displacement_r2['lats'] = np.flipud(displacement_r2['lats'])                                        # something is reveresed, so flip up down so that the highest lats are in the top row.  
        
    # 4: and also downsample other simple data if it's included:
    for product in ['dem', 'E', 'N', 'U']: 
        if product in displacement_r2.keys():
            displacement_r2[product] = rescale(displacement_r2[product], downsample_run, multichannel = False, anti_aliasing = False)                               # do the rescaling
        
        
    print(f"Interferogram were originally {shape_start} ({n_pixs_start} unmasked pixels), "
          f"but have been downsampled to {displacement_r2['mask'].shape} ({displacement_r2['incremental'].shape[1]} unmasked pixels) for use with LiCSAlert, "
          f"and have been downsampled to {displacement_r2['mask_downsampled'].shape} ({displacement_r2['incremental_downsampled'].shape[1]} unmasked pixels) for figures.  ")

    return displacement_r2


#%%
def bss_components_inversion(sources, interferograms, cumulative = True):
    """
    A function to fit an interferogram using components learned by BSS, and return how strongly
    each component is required to reconstruct that interferogramm, and the

    Inputs:
        sources | n_sources x pixels | ie architecture I.  Mean centered
        interferogram | n_ifgs x pixels | Doesn't have to be mean centered, ifgs are rows
        cumulative | Boolean | if true, m and residual (mean_l2_norm) are returned as cumulative sums.

    Outputs:
        m | rank 1 array | the strengths with which to use each source to reconstruct the ifg.
        mean_l2norm | float | the misfit between the ifg and the ifg reconstructed from sources

    2019/12/30 | MEG | Update so handles time series (and not single ifgs), and can return cumulative values
    """
    import numpy as np

    interferograms -= np.mean(interferograms)                     # mean centre
    (n_ifgs, n_pixels) = interferograms.shape

    d = interferograms.T                                                 # a column vector (p x 1)
    g = sources.T                                                       # a matrix of ICA sources and each is a column (p x n_sources)
    m = np.linalg.inv(g.T @ g) @ g.T @ d                                # m (n_sources x n_ifgs)
    d_hat = g@m                                                         # reconstructed ifgs, as column vectors
    d_resid = d - d_hat                                                 # residual between each ifg and its reconstruction

    m = m.T                                                             # make these column vectors
    residual = np.zeros((n_ifgs,1))                                     # residuals, as column vectors
    for i in range(n_ifgs):
        residual[i,] = np.sqrt(np.sum(d_resid[:,i]**2))/n_pixels         # the mean l2 norm for each ifg

    if cumulative:
        m = np.cumsum(m, axis=0)
        residual = np.cumsum(residual, axis=0)
    return m, residual


#%%
def time_course_rescaler(timecourses, temp_baselines):
    """A script to normalise timecourses so that ones that span long temporal baselines are normalised
    
    Neer finished?
    
    Inputs:
        timecourses | 2d array | as column (e.g. 77x5)
        temp_baselines | 2d array | as a column

    Outputs:
        timecourses | 2d array | as column (e.g. 77x5)
    """

    # timecourses_scaled = 12 * np.divide(timecourses, temp_baselines)                    # x12 so that no change for normal S1 12 day baseline
    # return timecourses_scaled




#%%
def filenames_to_baselines(phUnw_files):
    """
    Given a list of LiCSAR filenames, works out the temporal baselines from these.
    Inputs:
        phUnw_files | lst of strings | filename of each interferogram

    Outputs:
        baselines | column vector of temporal baselines
        baselines_cumulative | column vector of the cumulative sum of the baselines

    2017/12/1? | written as a script
    2017/12/19 | converte to a function
    2018/07/31 | update how dates are removed from filename to work with geocoded data (.geo.unw)
    2018/08/21 | change in format of how dates are removed.

    """
    import datetime
    import nump as np

    baselines = np.zeros((len(phUnw_files),1))

    i = 0
    for filename in phUnw_files:
        master_date = filename[:8]
        slave_date = filename[9:17]

        fmt = '%Y%m%d'                                                                      # tell datetime the format of the date (here year, month, day with no sepeartions)
        master_date_dt = datetime.datetime.strptime(master_date, fmt)
        master_date_tt = master_date_dt.timetuple()
        slave_date_dt = datetime.datetime.strptime(slave_date, fmt)
        slave_date_tt = slave_date_dt.timetuple()

        if master_date_tt.tm_year == slave_date_tt.tm_year:                                 # most ifgs don't span new year so are easy to deal with
            temporal_baseline = slave_date_tt.tm_yday - master_date_tt.tm_yday
        else:                                                                               # some do
            if master_date_tt.tm_year == 2016:                                              # and some could span leap years
                n_year_days = 366
            else:
                n_year_days = 365
            temporal_baseline = slave_date_tt.tm_yday + (n_year_days - master_date_tt.tm_yday)

        baselines[i,:] = temporal_baseline
        i += 1

    baselines_cumulative = np.cumsum(baselines)
    return baselines, baselines_cumulative


#%%
    
def save_pickle(fileout_name, *argv):
    """Save any number of items to a pickle file.  
    Inputs:
        fileout_name | str | filename
        *argv | items to be saved in pickle
    Outputs:
        fileout_name.pkl
    """
    import pickle

    with open(f"{fileout_name}.pkl", 'wb') as f:
        for arg in argv:                                        # loop through all inputs and save
            pickle.dump(arg, f)

  
    
#%%
            
def shorten_LiCSAlert_data(displacement_r2, n_end, n_start=0, verbose=False):
    """ Given a dictionary of ifgs for use with LiCSAlert, crop temporally (ie. the fist/vertical axis).  
    Inputs:
        displacement_r2 | dict | displacement is stored as row vectors in this
        n_end | int | ifg number to stop cropping at
        n_start | int | ifg number to start croppting at
        verbose | boolean | 
        
    Returns:
        displacement_r2_short | dict | as per input, but temporally cropped
        
    History:
        2020/01/10 | MEG  | Written
    """    
    import copy                                                                                                # needed to deepcopy dict
    displacement_r2_short = copy.deepcopy(displacement_r2)
    keys_to_shorten = ["incremental", "incremental_downsampled"]                                                    # only these items in the dict will be cropped
    for key_to_shorten in keys_to_shorten:
        try:
            displacement_r2_short[key_to_shorten] = displacement_r2_short[key_to_shorten][n_start:n_end,]
            if verbose:
                print(f"Succesfuly shortened {key_to_shorten} in the dictionary of interferograms")
        except:
            pass
            if verbose:
                print(f"{key_to_shorten} was not found in the dictionary of interferograms")
    return displacement_r2_short



#%% Small functions used by multiple function in this file
    
                     

def make_colormap(seq):
    """
    Taken from Stackechange - https://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
    Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).

    e.g. useage:  rvb = make_colormap(  [c('black'), c('orange'), 0.33, c('orange'), c('yellow'), 0.66, c('yellow'), c('red')])
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors


    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)



def remappedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the median value of a colormap, and scale the
    remaining color range (i.e. truncate the colormap so that it isn't
    compressed on the shorter side) . Useful for data with a negative minimum and
    positive maximum where you want the middle of the colormap's dynamic
    range to be at zero.
    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and 0.5; if your dataset mean is negative you should leave
          this at 0.0, otherwise to (vmax-abs(vmin))/(2*vmax)
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0; usually the
          optimal value is abs(vmin)/(vmax+abs(vmin))
          Only got this to work with:
              1 - vmin/(vmax + abs(vmin))
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          0.5 and 1.0; if your dataset mean is positive you should leave
          this at 1.0, otherwise to (abs(vmin)-vmax)/(2*abs(vmin))

      2017/??/?? | taken from stack exchange
      2017/10/11 | update so that crops shorter side of colorbar (so if data are in range [-1 100],
                   100 will be dark red, and -1 slightly blue (and not dark blue))
      '''
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt

    if midpoint > 0.5:                                      # crop the top or bottom of the colourscale so it's not asymetric.
        stop=(0.5 + (1-midpoint))
    else:
        start=(0.5 - midpoint)


    cdict = { 'red': [], 'green': [], 'blue': [], 'alpha': []  }
    # regular index to compute the colors
    reg_index = np.hstack([np.linspace(start, 0.5, 128, endpoint=False),  np.linspace(0.5, stop, 129)])

    # shifted index to match the data
    shift_index = np.hstack([ np.linspace(0.0, midpoint, 128, endpoint=False), np.linspace(midpoint, 1.0, 129)])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)
    return newcmap
