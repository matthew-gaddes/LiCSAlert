# -*- coding: utf-8 -*-
"""
Created on Tue May 29 11:23:10 2018

@author: eemeg
"""

import pdb

#%%

def ICASAR(n_pca_comp_start, n_pca_comp_stop, 
           spatial_data = None, temporal_data = None, figures = "window", 
           sica_tica = 'sica', ifgs_format = 'cum', max_n_all_ifgs = 1000,                                                     # this row of arguments are only needed with spatial data.  
           bootstrapping_param = (200,0), ica_param = (1e-4, 150), 
           tsne_param = (30,12), hdbscan_param = (35,10),
           out_folder = './ICASAR_results/', ica_verbose = 'long', 
           inset_axes_side = {'x':0.1, 'y':0.1}, 
           load_fastICA_results = False, label_sources = False):
    """
    Perform ICASAR, which is a robust way of applying sICA to data.  As PCA is also performed as part of this,
    the sources and time courses found by PCA are also returned.  Note that this can be run with eitehr 1d data (e.g. time series for a GPS station),
    or on 2d data (e.g. a time series of interferograms) by providing a 'mask', that is used to convert 1d row vectors to 2d masked arrays.  
    
    A note on reference areas/pixels:
        ICASAR requires each interferogram to be mean centered (ie the mean of all the pixels for a single interferogram is 0).  
        Therefore, when the time series is reconstructed using the result of ICASAR (i.e. tcs * sources), these will produce 
        the mean centered time series.  If you wish to work 
    
    Inputs:
        n_pca_comp_start | int | PCA is run for dimensionality reduction, start
                                 with this number of components.  
        n_pca_comp_stop | int | and stop running it for every number of components
                                until this value is reached.  
        spatial_data | dict or None | Required: 
                                         ifgs_dc | rank 2 array | row vectors of the daisy chain (i.e incremental) ifgs
                                         mask  | rank 2 array | mask to conver the row vectors to rank 2 masked arrays.  
                                         ifg_dates_dc | list | dates of the interferograms in the form YYYYMMDD_YYYYMMDD.  If supplied, IC strength vs temporal baseline plots will be produced.  
                                      Optional (ie don't have to exist in the dictionary):
                                         lons | rank 2 array | lons of each pixel in the image.  Changed to rank 2 in version 2.0, from rank 1 in version 1.0  .  If supplied, ICs will be geocoded as kmz.  
                                         lats | rank 2 array | lats of each pixel in the image. Changed to rank 2 in version 2.0, from rank 1 in version 1.0
                                         dem | rank 2 array | height in metres of each pixel in the image.  If supplied, IC vs dem plots will be produced.  
                                         
                                         
        temporal_data | dict or None | contains 'ifgs_inc' as time signals as row vectors and 'xvals' which are the times for each item in the time signals.   
        sica_tica | string | If not provided, spatial ICA (sICA) is performed, that is the spatial patterns are independent, rather than the time courses.  
        ifgs_format | string | 'all' / 'inc' / 'cum'
                                Spatial_data contains incremental ifgs (i.e. the daisy chain), these can be recombined to create interferograms so giving three options:
                                        1) incremental ifgs - just the daisy chain
                                        2) cumulative ifgs - relative to a single master, chosen as teh first acquisition.  
                                        3) all ifgs - all possible combinations between all acquisitions.Encompases 1 and 2
                                            e.g. for 3 interferogams between 4 acquisitions: a1__i1__a2__i2__a3__i3__a4
                                           This option would also make: a1__i4__a3, a1__i5__a4, a2__i6__a4
        max_n_all_ifgs | If after creating all the ifgs there are more than this number, select only this many at random.  Useful as the number of ifgs created grows with the square of the number of ifgs.  
        figures | string,  "window" / "png" / "none" / "png+window" | controls if figures are produced, noet none is the string none, not the NoneType None
        
        bootstrapping_param | tuple | (number of ICA runs with bootstrap, number of ICA runs without bootstrapping )  e.g. (100,10)
        ica_param | tuple | Used to control ICA, (ica_tol, ica_maxit)
        hdbscan_param  | tuple | Used to control the clustering (min_cluster_size, min_samples)
        tsne_param     | tuple | Used to control the 2d manifold learning  (perplexity, early_exaggeration)
        
        

        out_folder | string | if desired, can set the name of the folder results are saved to.  Should end with a /
        
        ica_verbose | 'long' or 'short' | if long, full details of ICA runs are given.  If short, only the overall progress 
        inset_axes_side | dict | inset axes side length as a fraction of the full figure, in x and y direction in the 2d figure of clustering results.  
        

        load_fastICA_results | boolean | The multiple runs of FastICA are slow, so if now paramters are being changed here, previous runs can be reloaded.  

    Outputs:
        S_best | rank 2 array | the recovered sources as row vectors (e.g. 5 x 1230)
        mask | rank 2 boolean | Same as inputs, but useful to save.  mask to convert the ifgs as rows into rank 2 masked arrays.  Used for figure outputs, an
        tcs    | rank 2 array | the time courses for the recoered sources (e.g. 17 x 5)
        source_residuals    | ? | the residual when each input mixture is reconstructed using the sources and time courses 
        Iq_sorted              | ?|     the cluster quality index for each centrotype
        n_clusters  | int | the number of clusters found.  Doens't include noise, so is almost always 1 less than the length of Iq
        S_all_info | dictionary| useful for custom plotting. Sources: all the sources in a rank 3 array (e.g. 500x500 x1200 for 6 sources recovered 200 times)
                                                            labels: label for each soure
                                                            xy: x and y coordinats for 2d representaion of all sources
        phUnw_mean | r2 array | the mean for each interfeorram.  subtract from (tcs * sources) to get back original ifgs.  
        
    History:
        2018/06/?? | MEG | Written
        2019/11/?? | MEG | Rewrite to be more robust and readable
        2020/06/03 | MEG | Update figure outputs.  
        2020/06/09 | MEG | Add a raise Exception so that data cannot have nans in it.  
        2020/06/12 | MEG | Add option to name outfolder where things are saved, and save results there as a pickle.  
        2020/06/24 | MEG | Add the ica_verbose option so that ICASAR can be run without generating too many terminal outputs.  
        2020/09/09 | MEG | Major update to now handle temporal data (as well as spatial data)
        2020/09/11 | MEG | Small update to allow an argument to be passed to plot_2d_interactive_fig to set the size of the inset axes.  
        2020/09/16 | MEG | Update to clarify the names of whether variables contain mixtures or soruces.  
        2021/04/13 | MEG | Update so that lons and lats are now rank2 tensors (ie matrices with a lon or lat for each pixel)
        2021/04/13 | MEG | Add option to create_all_ifgs_from_incremental
        2021_10_07 | MEG | Add option to limit the number of ifgs created from incremental. (e.g. if 5000 are generated but default value of 1000 is used, 1000 will be randomly chosen from the 5000)
        2021_10_20 | MEG | Also save the 2d position of each source, and its HDBSSCAN label in the .pickle file.  
        2022_01_18 | MEG | Add option to use cumulative (i.e. single master) interferograms.  
    
    Stack overview:
        PCA_meg2                                        # do PCA
        maps_tcs_rescale                                # rescale spatial maps from PCA so they have the same range, then rescale time courses so no change.  makes comparison easier.  )
        pca_variance_line                               # plot of variance for each PC direction
        component_plot                                  with PCA sources
        bootstrap_ICA                                   with bootstrapping
        bootstrap_ICA                                   and without bootstrapping
        bootstrapped_sources_to_centrotypes             # run HDBSCAN (clustering), TSNE (2d manifold) and make figure showing this.  Choose source most representative of each cluster (centrotype).  
        plot_2d_interactive_fig                         # interactive figure showing clustering and 2d manifold representaiton.  
        bss_components_inversion                        # inversion to get time coures for each centrotype.  
        component_plot                                  # with ICASAR sources
        r2_arrays_to_googleEarth                        # geocode spatial sources and make a .kmz for use with Google Earth.  
    """


    import matplotlib.pyplot as plt
    plt.switch_backend('Qt5Agg')
    import sys
    import pdb
    # print(f"\n\n\nDEBUG IMPORT\n\n\n")
    # sys.path.append("/home/matthew/university_work/python_stuff/python_scripts")
    # from small_plot_functions import matrix_show

    # external functions
    import numpy as np
    import numpy.ma as ma  
    import matplotlib.pyplot as plt
    import shutil                                                                # used to make/remove folders etc
    import os                                                                    # ditto
    import pickle                                                                # to save outputs.  
    from pathlib import Path
    import pdb
    
    # internal functions
    from licsalert.data_importing import ifg_timeseries
    from licsalert.icasar.blind_signal_separation import PCA_meg2
    from licsalert.icasar.aux import  bss_components_inversion, maps_tcs_rescale
    from licsalert.icasar.aux import r2_to_r3, create_all_ifgs, create_cumulative_ifgs
    from licsalert.icasar.plotting import r2_arrays_to_googleEarth
    from licsalert.icasar.plotting import plot_pca_variance_line, plot_temporal_signals
    from licsalert.icasar.plotting import two_spatial_signals_plot
    
    
    
    from licsalert.icasar.sources_interactive_figure import plot_2d_interactive_fig
    
    
    
    # -5 Check inputs, unpack either spatial or temporal data, and check for nans
    if temporal_data is None and spatial_data is None:                                                                  # check inputs
        raise Exception("One of either spatial or temporal data must be supplied.  Exiting.  ")
    elif temporal_data is not None and spatial_data is not None:
        raise Exception("Only either spatial or temporal data can be supplied, but not both.  Exiting.  ")
    elif temporal_data is None and spatial_data is not None:        
        print(f"From the data passed to ICASAR, I understand it to be spatial data (e.g. time series of images)")
        spatial = True
    elif temporal_data is not None and spatial_data is None:        
        print(f"From the data passed to ICASAR, I understand it to be temporal data (e.g. time series of sound recordings)")
        spatial = False
    
    # -4: main section for checking inputs
    if spatial:                                                                                                                 # if we have spatial data
        if np.max(np.isnan(spatial_data['ifgs_dc'])):
            raise Exception("Unable to proceed as the data ('spatial_data['dc_ifgs']') contains Nans.  ")
            
        if isinstance(spatial_data['ifgs_dc'], ma.MaskedArray):
            raise Exception("The interferograms in spatial_data ('ifgs_dc') "
                            "is a masked array, but must be a numpy array.  "
                            "Exiting.   ")

        
        mask = spatial_data['mask']                                                                                         # the mask that converts row vector mixtures into 2d (rank 2) arrays.  
        n_ifgs = spatial_data['ifgs_dc'].shape[0]                                                                           # get the number of incremental ifgs
        if n_ifgs != len(spatial_data['ifg_dates_dc']):                                                                     # and check it's equal to the list of ifg dates (YYYYMMDD_YYYYMMDD)
            raise Exception(f"There should be an equal number of incremental interferogram and dates (in the form YYYYMMDD_YYYYMMDD), but they appear to be different.  Exiting...")
            
        # check the arrays are the right size
        spatial_data_r2_arrays = ['mask', 'dem', 'lons', 'lats']                                                             # we need to check the spatial data is the correct resolution (ie all the same)
        spatial_data_r2_arrays_present = list(spatial_data.keys())                                                           # we alse need to determine which of these spatial data we actually have.  
        spatial_data_r2_arrays = [i for i in spatial_data_r2_arrays if i in spatial_data_r2_arrays_present]                  # remove any from the check list incase they're not provided.  
        for spatial_data_r2_array1 in spatial_data_r2_arrays:                                                                # first loop through each spatial data
            for spatial_data_r2_array2 in spatial_data_r2_arrays:                                                            # second loo through each spatial data
                if spatial_data[spatial_data_r2_array1].shape != spatial_data[spatial_data_r2_array2].shape:                 # check the size is equal
                    raise Exception(f"All the spatial data should be the same size, but {spatial_data_r2_array1} is of shape {spatial_data[spatial_data_r2_array1].shape}, "
                                    f"and {spatial_data_r2_array2} is of shape {spatial_data[spatial_data_r2_array2].shape}.  Exiting.")
        if 'dem' not in spatial_data_r2_arrays_present:                                                                  #  the dem is not compulsory
            spatial_data['dem'] = None                                                                                    # so set it to None if not available.  
            
        # check that the outputs can be geocoded.      
        if ('lons' in spatial_data) and ('lats' in spatial_data):                                                       # 
            print(f"As 'lons' and 'lats' have been provided, the ICs will be geocoded.  ")
            if (len(spatial_data['lons'].shape) != 2) or (len(spatial_data['lats'].shape) != 2):
                raise Exception(f"'lons' and 'lats' should be rank 2 tensors (i.e. matrices with a lon or lat for each pixel in the interferogram.  Exiting...  ")
            ge_kmz = True
        elif ('lons' in spatial_data) and ('lats' not in spatial_data):
            raise Exception(f"Either both or neither of 'lons' and 'lats' should be provided, but only 'lons' was.  Exiting...  ")
        elif ('lons' not in spatial_data) and ('lats' in spatial_data):
            raise Exception(f"Either both or neither of 'lons' and 'lats' should be provided, but only 'lats' was.  Exiting...  ")
        else:
            ge_kmz = False
        
        if sica_tica == "sica":                                                                                                 # spatial sources can be recoverd one of two ways.  sICA
            print(f"Spatial patterns will be statitically independent and ", end = '')
            if ifgs_format == 'all':
                print(f"All possible interferograms will be made between the acquisitions (up to max_n_all_ifgs).  ")
            elif ifgs_format == 'inc':
                print(f"The incremental (daisy chain) interferograms will be used.  ")
            elif ifgs_format == 'cum':
                print(f"The cumulative (single master) interferograms will be used.  ")
            else:
                raise Exception(f"'ifgs_format' must be either 'all', 'inc' or 'cum, but is {ifgs_format}.  Exiting.  ")
        elif sica_tica == "tica":
            if ifgs_format != 'cum':
                print(f"With tica (tICA),  only the cumulative ifgs can be used.  Updating the ifgs_format setting to reflect this.  ")
                ifgs_format = 'cum'
            print(f"For tICA, the cumulative interferograms are required.  These will be calculated from the incremental (daisy chain) interferograms.  ")
                
    else:                                                                                                                           # or we could do temporal data.  
        xvals = temporal_data['xvals']
        print(f"Deleting the 3 arguments that only apply to spatial data (sica_tica, create_all_ifgs_flag, max_n_all_ifgs) ")
        del sica_tica, ifgs_format, max_n_all_ifgs
        if np.max(np.isnan(temporal_data['mixtures_r2'])):
            raise Exception("Unable to proceed as the data ('spatial_data['mixtures_r2']') contains Nans.  ")
   
                       
    #-3:  sort out various things for figures, and check input is of the correct form
    if type(out_folder) == str:
        print(f"Trying to convert the 'out_folder' arg which is a string to a pathlib Path.  ")
        out_folder = Path(out_folder)
    fig_kwargs = {"figures" : figures}
    if figures == "png" or figures == "png+window":                                                         # if figures will be png, make 
        fig_kwargs['png_path'] = out_folder                                                                  # this will be passed to various figure plotting functions
    elif figures == 'window' or figures == 'none':
        pass
    else:
        raise ValueError("'figures' should be 'window', 'png', 'png+window', or 'None'.  Exiting...")
    if figures == 'png':
        plt.switch_backend('agg')                                                                       # with this backend, no windows are created during figure creation.  
    
    # -2: Check argument
    if ica_verbose == 'long':
        fastica_verbose = True
    elif ica_verbose == 'short':
        fastica_verbose = False
    else:
        print(f"'ica_verbose should be either 'long' or 'short'.  Setting to 'short' and continuing.  ")
        ica_verbose = 'short'
        fastica_verbose = False


    # -1: create a folder that will be used for outputs
    if os.path.exists(out_folder):                                                                      # see if the folder we'll write to exists.  
        if load_fastICA_results:                                                                        # we will need the .pkl of results from a previous run, so can't just delete the folder.  
            existing_files = os.listdir(out_folder)                                                     # get all the ICASAR outputs.  
            print(f"As 'load_fastICA' is set to True, all but the FastICA_results.pkl file will be deleted.   ")
            for existing_file in existing_files:
                if existing_file == 'FastICA_results.pkl':                                              # if it's the results from the time consuming FastICA runs...
                    pass                                                                                # ignore it    
                else:
                    os.remove(out_folder / existing_file)                                               # but if not, delete it.  
        else:
            print("Removing the existing outputs directory and creating a new empty one... ", end = '')                                 # if we don't care about the FastICA results file, just delete the folder and then make a new one.  
            shutil.rmtree(out_folder)                                                                   # try to remove folder
            os.mkdir(out_folder)
            print("Done.")
    else:
        os.mkdir(out_folder)                                                                            # if it never existed, make it.  
                                                                                           
    
    n_converge_bootstrapping = bootstrapping_param[0]                 # unpack input tuples
    n_converge_no_bootstrapping = bootstrapping_param[1]  
    

    # -0:  Create all interferograms, create the arary of mixtures (X), and mean centre
    if spatial:
        print(f"Creating all possible variations of the time series (incremental/daisy chain, cumulative, and all possible).  ")        # we allrady have the daisy chain ifgs.  
        ifgs_all_r2, ifg_dates_all = create_all_ifgs(spatial_data['ifgs_dc'], spatial_data['ifg_dates_dc'], max_n_all_ifgs)             # create all ifgs, even if we don't use them.  
        ifgs_cum_r2, ifg_dates_cum = create_cumulative_ifgs(spatial_data['ifgs_dc'], spatial_data['ifg_dates_dc'])                      # create the cumulative ifgs, even if we don't use them
        ifgs_dc = ifg_timeseries(spatial_data['ifgs_dc'], spatial_data['ifg_dates_dc'])                                                 # create a class (an ifg_timeseries) using the daisy chain ifgs
        ifgs_all = ifg_timeseries(ifgs_all_r2, ifg_dates_all)                                                                           # create a class (an ifg_timeseries) using all possible ifgs
        ifgs_cum = ifg_timeseries(ifgs_cum_r2, ifg_dates_cum)                                                                           # create a class (an ifg_timeseries) using the cumualtive ifgs.  
        del ifgs_all_r2, ifg_dates_all, ifgs_cum_r2, ifg_dates_cum
        
        if sica_tica == 'sica':
            X_mean = ifgs_dc.means_space                                                                                                # daisy chain ifgs are supplied, so only interested in returning daisy chain means.  
            # get the mean centered signals, depending on which of the 3 ifg formats
            if ifgs_format == 'all':
                X_mc = ifgs_all.mixtures_mc_space                                                                                       # the mixtures (used by PCA and ICA) can be all possible ifgs...
            elif ifgs_format == 'inc':
                X_mc = ifgs_dc.mixtures_mc_space                                                                                        # or just the incremental (daisy chain) ifgs.....
            elif ifgs_format == 'cum':
                X_mc = ifgs_cum.mixtures_mc_space                                                                                       # or the cumulative (single master) ifgs.  
        elif sica_tica == 'tica':                                                                                                       # if we're doing temporal ica with spatial data, the mixtures need to be the transpose
            X_mc = ifgs_cum.mixtures_mc_time.T                                                                                          # as cumulative and transpose, effectively the time series for each point.  
            X_mean = ifgs_cum.means_time
    else:
        X = temporal_data['mixtures_r2']
        X_mean = np.mean(X, axis = 1)[:,np.newaxis]                                                                                     # get the mean for each ifg (ie along rows.  )
        X_mc = X - X_mean                                                                                                               # mean centre the data (along rows)
        
        
    # 1: do sPCA once (and possibly create a figure of the PCA sources)
    print('Performing PCA to whiten the data....', end = "")
    success = False
    count = 0
    while (success == False) and (count < 10):
        
        try:
            outputs = PCA_meg2(X_mc, verbose = False)                    
            (PC_vecs, PC_vals, PC_whiten_mat, PC_dewhiten_mat, x_mc,
             x_decorrelate, x_white) = outputs; del outputs
            success = True
        except:
            success = False
            count += 1 
    if not success:
        raise Exception(f"PCA failed after {count} attempts.  This is most "
                        "common with tICA when the number of observations are "
                        "lower than the number of variables and the compact "
                        "trick is being used.  ")
    A_pca = PC_vecs                                                                                                                     
    S_pca = x_decorrelate     

    # rescale the signals 
    if spatial:
        if sica_tica == 'sica':
            S_pca, A_pca = maps_tcs_rescale(S_pca, A_pca)                                                        
        elif sica_tica == 'tica':
            # spatial patterns are in A, time courses are in S
            A_pca, S_pca_cum = maps_tcs_rescale(A_pca.T, S_pca.T)                                               
            # these are the cumulative time courses, should still be mean centered (have checked this)  
            S_pca_cum = S_pca_cum.T                                                                                                     
            A_pca = A_pca.T
            del S_pca
    else:
        S_pca = S_pca                                                                                                       # truncate to desired number of components
        A_pca =  A_pca
        
    
    # check that these have been set within limits
    if n_pca_comp_start > n_pca_comp_stop:
        raise Exception(f"The starting number of PCA components must be "
                        "larger than the stopping number, but they are "
                        f"{n_pca_comp_start} and {n_pca_comp_stop}.  Exiting.  ")

    # plot the PCA figures 
    if fig_kwargs['figures'] != "none":
        plot_pca_variance_line(PC_vals, title = '01_PCA_variance_line', **fig_kwargs)
        
            
        if spatial:
            if sica_tica == 'sica':
                S_pca = S_pca[:n_pca_comp_stop,:]
                A_pca = A_pca[:,:n_pca_comp_stop]
                inversion_results = bss_components_inversion(S_pca, 
                                                             [ifgs_dc.mixtures_mc_space,
                                                              ifgs_all.mixtures_mc_space])            
                source_residuals = inversion_results[0]['residual']                                                                     
                A_pca_dc = inversion_results[0]['tcs'].T                                                                                
                A_pca_all = inversion_results[1]['tcs'].T                                                                               
                two_spatial_signals_plot(S_pca, spatial_data['mask'], 
                                         spatial_data['dem'], 
                                         A_pca_dc, A_pca_all, 
                                         ifgs_dc.t_baselines, 
                                         ifgs_all.t_baselines,
                                         f"02_PCA_{n_pca_comp_stop}_components", 
                                         spatial_data['ifg_dates_dc'], 
                                         fig_kwargs)
            elif sica_tica == 'tica':
                A_pca = A_pca[:, :n_pca_comp_stop]
                S_pca_cum = S_pca_cum[:n_pca_comp_stop, :]
                S_pca_dc = np.diff(S_pca_cum, axis = 1, prepend = 0)                                                                    
                two_spatial_signals_plot(A_pca.T, spatial_data['mask'], 
                                         spatial_data['dem'], 
                                         S_pca_dc.T, S_pca_cum.T, 
                                         ifgs_dc.t_baselines, 
                                         ifgs_cum.t_baselines,                            
                                         f"02_PCA_{n_pca_comp_stop}_components", 
                                         spatial_data['ifg_dates_dc'], 
                                         fig_kwargs)                    
        else:
            
            S_pca = S_pca[:n_pca_comp_stop, ]
            plot_temporal_signals(S_pca, '02_PCA_sources', **fig_kwargs)
    
    # 2: Make or load the results of the multiple ICA runs.  
    if load_fastICA_results:
        print(f"Loading the results of multiple FastICA runs.  ")
        try:
            with open(out_folder / 'FastICA_results.pkl', 'rb') as f:
                S_hists = pickle.load(f)   
                A_hists = pickle.load(f)
        except:
            print(f"Failed to open the results from the previous runs of FastICA. "
                  f"Switching 'load_fastICA_results' to False and trying to "
                  "continue anyway.  ")
            load_fastICA_results = False
    if not load_fastICA_results:
       print(f"No results were found for the multiple ICA runs, so these will now be performed.  ")
       # initiliase to store
       S_hists = []
       A_hists = []
       #  +1 to make upper limit inclusive
       for n_comp in np.arange(n_pca_comp_start, n_pca_comp_stop + 1):
           print(f"\n\nStarting the multiple ICA runs with {n_comp} "
                 "principal components")
           S_hist, A_hist = perform_multiple_ICA_runs(n_comp, X_mc,
                                                      bootstrapping_param, 
                                                      ica_param, x_white, 
                                                      PC_dewhiten_mat, ica_verbose)                                                    
           S_hists.append(S_hist)
           A_hists.append(A_hist)
           
           # # debug (should be n_pca_comp_start to n_pca_comp_stop)
           # for S_hist in S_hists:
           #     print(f"Number of dimensions used: {S_hist[0].shape[0]}")
           
       with open(out_folder / 'FastICA_results.pkl', 'wb') as f:
           pickle.dump(S_hists, f)
           pickle.dump(A_hists, f)
            
    if spatial:                                                                                                                        
        if sica_tica == 'sica':
            # if it's spatial, images are rows in S
            n_pixels_loaded = S_hists[0][0].shape[1]
        elif sica_tica == 'tica':
            # if it's temporal, images are columns in A
            n_pixels_loaded = A_hist[0].shape[0]                                                                
        # check that the number of pixels in an image is correct
        # (often get errors here as data could be loaded after settings were chnaged
        if n_pixels_loaded != np.sum(1-spatial_data['mask']):
            raise Exception(f"There are {S_hist[0].shape[1]} pixels in the ICASAR sources that have been loaded, but"
                            f" {np.sum(1-spatial_data['mask'])} pixels in the current mask.  This normally happens when the"
                            f" FastICA results that are being loaded are from a different set of data.  If not, something "
                            f" is inconsitent with the mask and the coherent pixels.  Exiting.  ")

    plot_2d_labels = {'title' : '04_clustering_and_manifold_results',
                      'xlabel' : 'TSNE dimension 1',
                      'ylabel' : 'TSNE dimension 2'}
        
    if spatial:
        
        
        
        # # Save each variable individually to the pickle file
        # with open('figure_inputs.pkl', 'wb') as f:
        #     pickle.dump(S_pca, f)
        #     pickle.dump(S_hists, f)
        #     pickle.dump(mask, f)
        #     pickle.dump(spatial, f)
        #     pickle.dump(sica_tica, f)
        #     pickle.dump(hdbscan_param, f)
        #     pickle.dump(tsne_param, f)
        #     pickle.dump(n_converge_bootstrapping, f)
        #     pickle.dump(n_converge_no_bootstrapping, f)
        #     pickle.dump(inset_axes_side, f)
            
            
        
        S_ica, source_outputs = plot_2d_interactive_fig(S_pca, S_hists, mask, spatial, sica_tica, 
                                                        hdbscan_param, tsne_param,
                                                        n_converge_bootstrapping, 
                                                        n_converge_no_bootstrapping,
                                                        inset_axes_side = inset_axes_side,
                                                        fig_filename = plot_2d_labels['title'], 
                                                        **fig_kwargs)
        
        # unpack to use the previous naming convention.  
        labels_hdbscan = source_outputs['labels']
        xy_tsne = source_outputs['xy']
        n_clusters = source_outputs['n_clusters']
        Iq_sorted = source_outputs['Iq_sorted']
        sources_all_r2 = source_outputs['sources_all_r2']
        

    else:
        raise Exception(f"Fucntion removed. ")
        # temporal_data_S_all = {'tcs_r2' : sources_all_r2,
        #                        'xvals'  : temporal_data['xvals'] }                                                                               # make a dictionary of the sources recovered from each run
        # plot_2d_interactive_fig(xy_tsne.T, colours = labels_colours, 
        #                         temporal_data = temporal_data_S_all,                                        # make the 2d interactive plot
        #                         labels = plot_2d_labels, legend = legend_dict, 
        #                         markers = marker_dict, 
        #                         inset_axes_side = inset_axes_side,
        #                         fig_filename = plot_2d_labels['title'], 
        #                         **fig_kwargs)


    
    # 5: Make time courses using centrotypes (i.e. S_ica, the spatial patterns found by ICA), or viceversa if tICA 
    if spatial: 
        if sica_tica == 'sica':
            inversion_results = bss_components_inversion(S_ica, [ifgs_dc.mixtures_mc_space, ifgs_all.mixtures_mc_space])                                           # invert to fit the incremetal ifgs and all ifgs using spatial patterns (that are in Sica)
            source_residuals = inversion_results[0]['residual']                                                                                                    # also get how well each daisy chain (mean cenetered) interferogram is fit                 
            A_ica_dc = inversion_results[0]['tcs'].T                                                                                                               # in sICA, time courses are in A
            A_ica_all = inversion_results[1]['tcs'].T                                                                                                              # time courses to fit all possible interferograms.                                   
            if fig_kwargs['figures'] != "none":                
                dem_to_sources_comparisons, tcs_to_tempbaselines_comparisons =   two_spatial_signals_plot(S_ica, spatial_data['mask'], spatial_data['dem'], 
                                                                                                          A_ica_dc, A_ica_all, ifgs_dc.t_baselines, ifgs_all.t_baselines,
                                                                                                          "03_ICA_sources", spatial_data['ifg_dates_dc'], fig_kwargs)
        elif sica_tica == 'tica':
            S_ica_cum = S_ica                                                                                                                                     # if temporal, sources are time courses, and are for the cumulative ifgs (as the transpose of these was given to the ICA function).  Note that these each have a mean of 0
            S_ica_dc = np.diff(S_ica_cum, axis = 1, prepend = 0)                                                                                                  # the diff of the cumluative time courses is the incremnetal (daisy chain) time course.  Prepend a 0 to make it thesame size as the original diays chain (ie. the capture the difference between 0 and first value).  Note that these are no longer each mean centered
            del S_ica           
            inversion_results = bss_components_inversion(S_ica_cum, [ifgs_cum.mixtures_mc_time.T])                                                                # inversion to fit the time series for each pixel (why it's the transpose), using the sources which are time courses (as its tICA).  Note that eveything here is mean centered in time
            A_ica = inversion_results[0]['tcs']                                                                                                                   # in tICA, spatial sources are row vector (ie these are the images)
            source_residuals = inversion_results[0]['residual']                                                                                                   # how well the cumulative time courses are fit?  Not clear.  
            if fig_kwargs['figures'] != "none":
                dem_to_sources_comparisons, tcs_to_tempbaselines_comparisons = two_spatial_signals_plot(A_ica, spatial_data['mask'], spatial_data['dem'],                           # here A are the images, simply the result of the inversion and no consideration given to their mean
                                                                                                        S_ica_dc.T, S_ica_cum.T, ifgs_dc.t_baselines, ifgs_cum.t_baselines,          # 
                                                                                                        "03_ICA_sources", spatial_data['ifg_dates_dc'], fig_kwargs)                    
    else:
        inversion_results = bss_components_inversion(S_ica, [X_mc])                                                 # invert to fit the mean centered mixture.    
        source_residuals = inversion_results[0]['residual']                                                         # how well we fit those
        A_ica = inversion_results[0]['tcs'].T                                                                       # and the time coruses to remake them.                  
        plot_temporal_signals(S_ica, '04_ICASAR_sources', **fig_kwargs)
 
     
    # 7: Possibly geocode the recovered sources and make a Google Earth file.     
    if spatial:
        if ge_kmz:
            print('Creating a Google Earth .kmz of the geocoded independent components... ', end = '')
            if sica_tica == 'sica':
                S_ica_r3 = r2_to_r3(S_ica, mask)                                                                                                        # if sica, sources (S) are images
            elif sica_tica == 'tica':
                S_ica_r3 = r2_to_r3(A_ica, mask)                                                                                                            # but if tica, time courses (A) are images.  
            r2_arrays_to_googleEarth(S_ica_r3, spatial_data['lons'], spatial_data['lats'], 'IC', out_folder = out_folder)                              # note that lons and lats should be rank 2 (ie an entry for each pixel in the ifgs)
            print('Done!')
    
    # 8: Try to determine which sources relate to which processes (e.g. deformation, topo. correlated APS, or turbulent APS.  ) 
    if label_sources:
        if fig_kwargs['figures'] == "none":
            print(f"'label_sources' rquires figures to be made (as the correlations are calculated when making figures.  Continuing, but 'label_sources' is being set to False ")
            label_sources == False
        else:
            if spatial:
                print(f"\n")
                if sica_tica == 'sica':
                    n_sources = S_ica.shape[0]                                                                                          # if spatial, sources are images and rows.  
                elif sica_tica == 'tica':
                    n_sources = S_ica_cum.shape[0]                                                                                      # if temporal, sources are (cumulative) time courses and rows.  
                label_sources_output = {'source_names' : ['deformation', 'topo_cor_APS', 'turbulent_APS'],
                                        'labels'       : np.zeros((n_sources, 3))   }
                label_sources_output['labels'][np.argmax(np.abs(tcs_to_tempbaselines_comparisons['cor_coefs'])), 0] = 1                 # teh deformation labels (column 0) - most correalted in time 
                label_sources_output['labels'][np.argmax(np.abs(dem_to_sources_comparisons['cor_coefs'])), 1] = 1                       # the topo. correalted label (column 1) - most correalted with DEM
                for source_n in range(n_sources):                                                         # loop through all    
                    if np.max(label_sources_output['labels'][source_n, :]) != 1:                                                        # if don't have a label yet...
                        label_sources_output['labels'][source_n, 2] = 1                                                                 # label as turbulent (column 2)
                    print(f"Source {source_n}: {label_sources_output['source_names'][np.argmax(label_sources_output['labels'][source_n])]}")
                print(f"\n")
                
            else:
                print(f"'label_sources' is only supported with spatial data.  Setting this to False and trying to conitnue.  ")           # labelling of sources is intended for use with time series of intererograms and not temporal data.  
                label_sources == False
                
    if plt.get_backend() == 'agg':
        plt.switch_backend('Qt5Agg')

    # 9: Save the results: 
    S_all_info = {'sources' : sources_all_r2,                                                                # package into a dict to return
                  'labels' : labels_hdbscan,
                  'xy' : xy_tsne       }
    print('Saving the key results as a .pkl file... ', end = '')                                            
    if spatial:
        with open(out_folder / 'ICASAR_results.pkl', 'wb') as f:                            # first save spatial patterns, then time courses (incremental, not cumulative)
            if sica_tica == 'sica':
                pickle.dump(S_ica, f)                                                       # spatial patterns are row vectors
            elif sica_tica == 'tica':
                pickle.dump(A_ica, f)                                                       # spatial images are in A, neither mean centered in time or space
            pickle.dump(mask, f)
            if sica_tica == 'sica':                             
                pickle.dump(A_ica_dc, f)                                                    # time courses
            elif sica_tica == 'tica':
                pickle.dump(S_ica_dc.T, f)                                                  # time courses are sources and in S.  Note that because these are incremental, they are no longer mean centered.  

            pickle.dump(source_residuals, f)
            pickle.dump(n_clusters, f)
            pickle.dump(xy_tsne, f)
            pickle.dump(labels_hdbscan, f)
            if label_sources:
                pickle.dump(label_sources_output, f)
        f.close()
        
        
        # save the PCA outputs.  
        with open(out_folder / 'pca_results.pkl', 'wb') as f:
            if sica_tica == 'sica':
                # spatial maps as rows
                pickle.dump(S_pca, f)                                                      
            elif sica_tica == 'tica':
                # spatial sources are columns of A
                pickle.dump(A_pca.T, f)                                                       
            if sica_tica == 'sica':                             
                # time courses as columns
                pickle.dump(A_pca, f)                                                    
            elif sica_tica == 'tica':
                # time courses for daisy chain ifgs
                pickle.dump(S_pca_dc.T, f)                                                  
                # time courses for cumulative ifgs
                pickle.dump(S_pca_cum.T, f)                                                  
        f.close()
        
        print("Done!")
        


        if sica_tica == 'sica':
            if label_sources:
                return S_ica, A_ica_dc, source_residuals, Iq_sorted, n_clusters, S_all_info, X_mean, label_sources_output
            else:
                return S_ica, A_ica_dc, source_residuals, Iq_sorted, n_clusters, S_all_info, X_mean
        elif sica_tica == 'tica':
            if label_sources:
                return A_ica, S_ica_dc.T, source_residuals, Iq_sorted, n_clusters, S_all_info, X_mean, label_sources_output
            else:
                return A_ica, S_ica_dc.T, source_residuals, Iq_sorted, n_clusters, S_all_info, X_mean
        
    else:                                                                       # if temporal data, no mask to save
        with open(out_folder / 'ICASAR_results.pkl', 'wb') as f:
            pickle.dump(S_ica, f)
            pickle.dump(A_ica, f)
            pickle.dump(source_residuals, f)
            pickle.dump(Iq_sorted, f)
            pickle.dump(n_clusters, f)
            pickle.dump(xy_tsne, f)
            pickle.dump(labels_hdbscan, f)
        f.close()
        
        with open(out_folder / 'pca_results.pkl', 'wb') as f:
            pickle.dump(S_pca, f)
            pickle.dump(A_pca, f)
        f.close()
        
        print("Done!")
        return S_ica, A_ica, source_residuals, Iq_sorted, n_clusters, S_all_info, X_mean


    

#%%

def tsne_and_cluster(S_hist, mask, spatial, sica_tica, n_comp,
                     hdbscan_param, tsne_param,
                     n_converge_bootstrapping, n_converge_no_bootstrapping,
                     update_tsne = True):
    """ Given the results of multiple ICA runs, do the representation in 2D 
    space, and the clustering. 
    
    
    Inputs:
        S_hist | list | results of multiople ICA runs (one item in list for each run.  )
                        each item would be e.g. (6,1264).  
        mask | r2 array | to convert a row vector to rank 2 (image)
        spatial | boolean | True if data is spatial (e.g. ifgs)
        sica_tics | str | either sica or tica, for spatial or temporal. 
        n_comp | int | number of PCA components retained.  
        hdbscan_param | tuple | min cluster size, min samples
        tsne_param | tuple | perplexity, early exageration.  
        n_converge_bootstrapping | int | number of runs with bootstrapping
        n_converge_no_bootstrapping | int | number of runs without bootstrapping
        update_tsne | boolean | if True, TSNE is rerun.  
        
    Returns: (for e.g. 200 runs with 6 PCs)
        sources_all_r3 | r3 |x all sources recovered during runs (e.g. 1200, 200, 150)
        S_ica | r2 | ICA sources as row vectors (e.g. (6, 1264))
        labels_hdbscan | r1 | cluster label for each source (e.g. (1200,))
        xy_tsne | r2  | xy coords for each point in new 2D space (e.g. 1200 x 2)
        marker_dict  | dict | labels and styles, (0 for 'o', 1 for 'x')
        legend_dict | dict | elemetns and labels, colours for points (and their name)
        labels_colours | list | one colour as string for each point.  e.g. 1200 long,
                                array(['#d62728', '#ff7f0e', '#2ca02c'], dtype='<U7')
                                
    Histoy:
        2024_07_17 | MEG | Written
    """
    
    import numpy as np
    from licsalert.icasar.plotting import prepare_point_colours_for_2d
    from licsalert.icasar.plotting import prepare_legends_for_2d
    

    # 1: Convert the sources from lists from each run to a single matrix.  
    if spatial:
        # if its spatial dat and sica, sources are images
        if sica_tica == 'sica':                                                               
            # r2 one is (n_components x n_runs) x n_pixels, r3 one is (n_components x n_runs) x ny x nx,
            sources_all_r2, sources_all_r3 = sources_list_to_r2_r3(S_hist, mask)                            
        elif sica_tica == 'tica':
            #  that these are 1d displacments for a pixel through time.  
            sources_all_r2 = S_hist[0]                                                                      
            for S_hist_one in S_hist[1:]:                                                                   
                sources_all_r2 = np.vstack((sources_all_r2, S_hist_one))                                    
    else:             
        # time signals                                                                                   
        sources_all_r2 = S_hist[0]                                                                      
        for S_hist_one in S_hist[1:]:                                                                   
            sources_all_r2 = np.vstack((sources_all_r2, S_hist_one))                                    
                        
   
    # 2: Do clustering and 2d manifold representation, plus get centrotypes of 
    # clusters, and make an interactive plot.   
    # do the clustering and project to a 2d plane.  clusters_by_max_Iq_no_noise 
    # is an array of which cluster number is best (ie has the highest Iq)
    outputs  = bootstrapped_sources_to_centrotypes(sources_all_r2, hdbscan_param, 
                                                   tsne_param, update_tsne)        
    (S_ica, labels_hdbscan, xy_tsne, clusters_by_max_Iq_no_noise, Iq) = outputs
    Iq_sorted = np.sort(Iq)[::-1]               
    n_clusters = S_ica.shape[0]                                                                     
    # make a list of colours so that each point with the same label has the 
    # same colour, and all noise points are grey
    # e.g. array(['#d62728', '#ff7f0e', '#2ca02c', ..., '#d62728', '#ff7f0e', '#9467bd']
    labels_colours = prepare_point_colours_for_2d(labels_hdbscan, 
                                                  clusters_by_max_Iq_no_noise)                                                                
       
    # legend with colours for each cluster, and names for each cluster.  
    legend_dict = prepare_legends_for_2d(clusters_by_max_Iq_no_noise, Iq)    
    # boostrapped are labelled as 0, and non bootstrapped as 1
    marker_dict = {'labels' : np.ravel(np.hstack((np.zeros((1, n_comp*n_converge_bootstrapping)), 
                                                  np.ones((1, n_comp*n_converge_no_bootstrapping)))))}        
    marker_dict['styles'] = ['o', 'x']
    
    return sources_all_r3, S_ica, labels_hdbscan, xy_tsne, marker_dict, legend_dict, labels_colours, Iq_sorted, n_clusters
   



#%%
def update_mask_sources_ifgs(mask_sources, sources, mask_ifgs, ifgs):
    """ Given two masks of pixels, create a mask of pixels that are valid for both.  Also return the two sets of data with the new masks applied.  
    Inputs:
        mask_sources | boolean rank 2| original mask
        sources  | r2 array | sources as row vectors
        mask_ifgs | boolean rank 2| new mask
        ifgs  | r2 array | ifgs as row vectors
    Returns:
        ifgs_new_mask
        sources_new_mask
        mask_both | boolean rank 2| original mask
    History:
        2020/02/19 | MEG |  Written      
        2020/06/26 | MEG | Major rewrite.  
        2021_04_20 | MEG | Add check that sources and ifgs are both rank 2 (use row vectors if only one source, but it must be rank2 and not rank 1)
    """
    import numpy as np
    import numpy.ma as ma
    from icasar.aux import col_to_ma

        
    
    def apply_new_mask(ifgs, mask_old, mask_new):
        """Apply a new mask to a collection of ifgs (or sources) that are stored as row vectors with an accompanying mask.  
        Inputs:
            ifgs | r2 array | ifgs as row vectors
            mask_old | r2 array | mask to convert a row of ifg into a rank 2 masked array
            mask_new | r2 array | the new mask to be applied.  Note that it must not unmask any pixels that are already masked.  
        Returns:
            ifgs_new_mask | r2 array | as per ifgs, but with a new mask.  
        History:
            2020/06/26 | MEG | Written
        """
        n_pixs_new = len(np.argwhere(mask_new == False))                                        
        ifgs_new_mask = np.zeros((ifgs.shape[0], n_pixs_new))                        # initiate an array to store the modified sources as row vectors    
        for ifg_n, ifg in enumerate(ifgs):                                 # Loop through each source
            ifg_r2 = col_to_ma(ifg, mask_old)                             # turn it from a row vector into a rank 2 masked array        
            ifg_r2_new_mask = ma.array(ifg_r2, mask = mask_new)              # apply the new mask   
            ifgs_new_mask[ifg_n, :] = ma.compressed(ifg_r2_new_mask)       # convert to row vector and places in rank 2 array of modified sources
        return ifgs_new_mask
    
    
    # check some inputs.  Not exhuastive!
    if (len(sources.shape) != 2) or (len(ifgs.shape) != 2):
        raise Exception(f"Both 'sources' and 'ifgs' must be rank 2 arrays (even if they are only a single source).  Exiting. ")
    
    mask_both = ~np.logical_and(~mask_sources, ~mask_ifgs)                                       # make a new mask for pixels that are in the sources AND in the current time series
    n_pixs_sources = len(np.argwhere(mask_sources == False))                                  # masked pixels are 1s, so invert with 1- bit so that non-masked are 1s, then sum to get number of pixels
    n_pixs_new = len(np.argwhere(mask_ifgs == False))                                          # ditto for new mask
    n_pixs_both = len(np.argwhere(mask_both == False))                                        # ditto for the mutual mask
    print(f"Updating masks and ICA sources.  Of the {n_pixs_sources} in the sources and {n_pixs_new} in the current LiCSBAS time series, "
          f"{n_pixs_both} are in both and can be used in this iteration of LiCSAlert.  ")
    
    ifgs_new_mask = apply_new_mask(ifgs, mask_ifgs, mask_both)                                  # apply the new mask to the old ifgs and return the non-masked elemts as row vectors.  
    sources_new_mask = apply_new_mask(sources, mask_sources, mask_both)                         # ditto for the sources.  
    
    return ifgs_new_mask, sources_new_mask, mask_both
    




#%%

def bootstrapped_sources_to_centrotypes(sources_r2, hdbscan_param, tsne_param,
                                        update_tsne):
    """ Given the products of the bootstrapping, run the 2d manifold and
    clustering algorithms to create centrotypes.  
    Inputs:
        sources_r2      | rank 2 array | all the sources recovered after 
                                        bootstrapping.  If 5 components and 100
                                        bootstrapped runs, this will be 500 x n_pixels
                                        (or n_times)
        hdbscan_param  | tuple | Used to control the clustering (min_cluster_size, min_samples)
        tsne_param     | tuple | Used to control the 2d manifold learning  
                                (perplexity, early_exaggeration)
    Returns:
        S_best | rank 2 array | the recovered sources as row vectors (e.g. 5 x 1230)
        labels_hdbscan | rank 2 array | the cluster number for each of the 
                                        sources in sources_all_r2 e.g 1000,
        xy_tsne | rank 2 array | the x and y coordinates of where each space is 
                                 in the 2D space.  e.g. 1000x2
        clusters_by_max_Iq_no_noise | rank 1 array | clusters ranked by quality
                                                      index (Iq).  e.g. 3,0,1,4,2
        Iq | list | cluster quality index for each cluster.  Entry 0 is Iq 
                    (cluster quality index) for the first cluster

    History:
        2020/08/26 | MEG | Created from a script.  
        2021_04_16 | MEG | Remove unused figure arguments.  
    """
    import numpy as np
    import hdbscan                                                               
    from sklearn.manifold import TSNE                                            

    # unpack tuples
    perplexity = tsne_param[0]                                                   
    early_exaggeration = tsne_param[1]
    min_cluster_size = hdbscan_param[0]                                              
    min_samples = hdbscan_param[1] 

    # 1: Create the pairwise comparison matrix
    print('\nStarting to compute the pairwise distance matrices....', end = '')
    # each row of sources_r2 is a source, is column is a pixel.  There are
    # n_comp * n_bootstrapping sources, D are distances, S are similiarities 
    D, S = pairwise_comparison(sources_r2)                                              
    print('Done!')                                                                    


    #  2: Clustering with all the recovered sources   
    print('Starting to cluster the sources using HDBSCAN....', end = "")  
    clusterer_precom = hdbscan.HDBSCAN(metric = 'precomputed', min_cluster_size = min_cluster_size, 
                                       min_samples = min_samples, cluster_selection_method = 'leaf')
    labels_hdbscan = clusterer_precom.fit_predict(D)                                                                  # D is n_samples x n_samples, then returns a rank 1 which is the cluster number (ie label) for each source
    Iq = cluster_quality_index(labels_hdbscan, S)                                                                     # calculate the cluster quality index, using S (n_samples x n_samples), and the label for each one
                                                                                                                      # note that Iq is ordered by cluster, so the first value is the cluster quality index for 1st cluster (which is usually labelled -1 and the noise points)
    if np.min(labels_hdbscan) == (-1):                                                                                # if HDBSCAN has identified noise
        Iq = Iq[1:]                                                                                                   # delete the first entry, as this is the Iq of the noise (which isn't a cluster)
    clusters_by_max_Iq_no_noise = np.argsort(Iq)[::-1]                                                                # clusters by best Iqfirst (ie cluster)
    print('Done!')


    # 3:  2d manifold with all the recovered source (if required)
    if update_tsne:
        print('Starting to calculate the 2D manifold representation....', end = "")
        manifold_tsne = TSNE(n_components = 2, metric = 'precomputed', perplexity = perplexity, early_exaggeration = early_exaggeration,
                             init = 'random', learning_rate = 200.0)# , square_distances=True)                                                                               # default will change to pca in 1.2.  May be worth experimenting with.  
        xy_tsne = manifold_tsne.fit(D).embedding_
        print('Done!' )
    else:
        xy_tsne = None
    
# testing.      
#     def test_tsne(n_comp, D, perplexity, early_exaggeration):
#         """
#         """
#         from sklearn.manifold import TSNE                                            # t-distributed stochastic neighbour embedding
#         import matplotlib.pyplot as plt
# #        plt.switch_backend('Qt5Agg')
#         import sys
#         sys.path.append("/home/matthew/university_work/python_stuff/python_scripts")
#         from small_plot_functions import matrix_show
        
        
#         manifold_tsne = TSNE(n_components = n_comp, metric = 'precomputed', perplexity = perplexity, early_exaggeration = early_exaggeration,
#                              init = 'random', learning_rate = 200.0, square_distances=True, method = 'exact')                                                                               # default will change to pca in 1.2.  May be worth experimenting with.  
#         xy_tsne_2 = manifold_tsne.fit(D).embedding_
#         D_tsne_2, _= pairwise_comparison(xy_tsne_2)                                              # each row is a source, is column is a pixel.  There ar n_comp * n_bootstrapping sources.  
#         matrix_show([D_tsne_2], title = f'{n_comp} dims') 
        

#     matrix_show([D], title = 'D - high dim data')
#     test_tsne(100, D, perplexity, early_exaggeration)
#     test_tsne(3, D, perplexity, early_exaggeration)
    
    
    
    # 4: Determine the number of clusters from HDBSCAN
    if np.min(labels_hdbscan) == (-1):                                      # if we have noise (which is labelled as -1 byt HDBSCAN), 
        n_clusters = np.size(np.unique(labels_hdbscan)) - 1                 # noise doesn't count as a cluster so we -1 from number of clusters
    else:
        n_clusters = np.size(np.unique(labels_hdbscan))                     # but if no noise, number of clusters is just number of different labels
        
    if n_clusters == 0:
        print("No clusters have been found.  Often, this is caused by running the FastICA algorithm too few times, or setting"
              "the hdbscan_param 'min_cluster_size' too high.  ")
        return None, labels_hdbscan, xy_tsne, clusters_by_max_Iq_no_noise, Iq
    
    else:
        # 4:  Centrotypes (object that is most similar to all others in the cluster)
        print('Calculating the centrotypes and associated time courses...', end = '')
        S_best_args = np.zeros((n_clusters, 1)).astype(int)                         
        for i, clust_number in enumerate(clusters_by_max_Iq_no_noise):                              # loop through each cluster in order of how good they are (i.e. highest Iq first)
            source_index = np.ravel(np.argwhere(labels_hdbscan == clust_number))                    # get the indexes of sources in this cluster
            S_this_cluster = np.copy(S[source_index, :][:, source_index])                           # similarities for just this cluster
            in_cluster_arg = np.argmax(np.sum(S_this_cluster, axis = 1))                            # the sum of a column of S_this... is the similarity between 1 source and all the others.  Look for the column that's the maximum
            S_best_args[i,0] = source_index[in_cluster_arg]                                         # conver the number in the cluster to the number overall (ie 2nd in cluster is actually 120th source)     
        S_best = np.copy(sources_r2[np.ravel(S_best_args),:])                                       # these are the centrotype sources
        print('Done!' )
    
        return S_best, labels_hdbscan, xy_tsne, clusters_by_max_Iq_no_noise, Iq


#%%

def perform_multiple_ICA_runs(n_comp, mixtures_mc, bootstrapping_param, ica_param,
                              mixtures_white = None, dewhiten_matrix = None, ica_verbose = 'long'):
    """
    ICASAR requires ICA to be run many times, wither with or without bootstrapping.  This function performs this.  
    Inputs:
        n_comp | int | the number of souces we aim to recover.  
        mixutres_mc | rank 2 array | mixtures as rows, mean centered along rows.  I.e. of size n_varaibles x n_observations.  
        bootstrapping_param | tuple | (number of ICA runs with bootstrap, number of ICA runs without bootstrapping )  e.g. (100,10)
        ica_param | tuple | Used to control ICA, (ica_tol, ica_maxit)
        mixtures_white | rank 2 | mean centered and decorellated and unit variance in each dimension (ie whitened).  As per mixtures, row vectors.  
        dewhiten_matrix | rank 2 | n_comp x n_comp.  mixtures_mc = dewhiten_matrix @ mixtures_white
        ica_verbose | 'long' or 'short' | if long, full details of ICA runs are given.  If short, only the overall progress 
    Returns:
        S_best | list of rank 2 arrays | the sources from each run of the FastICA algorithm, n_comp x n_pixels.  Bootstrapped ones first, non-bootstrapped second.  
        A_hist | list of rank 2 arrays | the time courses from each run of the FastICA algorithm.  n_ifgs x n_comp.  Bootstrapped ones first, non-bootstrapped second.  
    History:
        2021_04_23 | MEG | Written
    """
    
    # 1: unpack a tuple and check a few inputs.  
    n_converge_bootstrapping = bootstrapping_param[0]                 # unpack input tuples
    n_converge_no_bootstrapping = bootstrapping_param[1]          
    
    if (n_converge_no_bootstrapping > 0) and ((mixtures_white is None) or (dewhiten_matrix is None)):
        raise Exception(f"If runs without bootstrapping are to be performed, the whitened data and the dewhitening matrix must be provided, yet one "
                        f"or more of these are 'None'.  This is as PCA is performed to whiten the data, yet if bootstrapping is not being used "
                        f"the data don't change, so PCA doesn't need to be run (and it can be computationally expensive).  Exiting.  ")
    
    # 2: do ICA multiple times   
    # First with bootstrapping
    A_hist_BS = []                                                                                      # ditto but with bootstrapping
    S_hist_BS = []
    n_ica_converge = 0
    n_ica_fail = 0
    if ica_verbose == 'short' and n_converge_bootstrapping > 0:                                 # if we're only doing short version of verbose, and will be doing bootstrapping
        print(f"FastICA progress with bootstrapping: ", end = '')
    while n_ica_converge < n_converge_bootstrapping:
        # this can be slow as PCA is run each time
        S, A, ica_converged = bootstrap_ICA(mixtures_mc, n_comp, bootstrap = True, 
                                            ica_param = ica_param, verbose = ica_verbose)                
        if ica_converged:
            n_ica_converge += 1
            # record results
            A_hist_BS.append(A)                                     
            S_hist_BS.append(S)                                     
        else:
            n_ica_fail += 1
        if ica_verbose == 'long':
            print(f"sICA with bootstrapping has converged {n_ica_converge} of {n_converge_bootstrapping} times.   \n")              # longer (more info) update to terminal
        else:
            print(f"{int(100*(n_ica_converge/n_converge_bootstrapping))}% ", end = '')                                              # short update to terminal

    # and without bootstrapping
    A_hist_no_BS = []                                                                        # initiate to store time courses without bootstrapping
    S_hist_no_BS = []                                                                        # and recovered sources        
    n_ica_converge = 0                                                                       # reset the counters for the second lot of ica
    n_ica_fail = 0
    if ica_verbose == 'short' and n_converge_no_bootstrapping > 0:                           # if we're only doing short version of verbose, and are actually doing ICA with no bootstrapping
        print(f"FastICA progress without bootstrapping: ", end = '')
    while n_ica_converge < n_converge_no_bootstrapping:
        # fast as PCA is not run each time.  
        S, A, ica_converged = bootstrap_ICA(mixtures_mc, n_comp, bootstrap = False,
                                            ica_param = ica_param,
                                            X_whitened = mixtures_white, 
                                            dewhiten_matrix = dewhiten_matrix, 
                                            verbose = ica_verbose)               
        
        if ica_converged:
            n_ica_converge += 1
            # record results
            A_hist_no_BS.append(A)                                     
            S_hist_no_BS.append(S)                                     
        else:
            n_ica_fail += 1
        if ica_verbose == 'long':
            print(f"sICA without bootstrapping has converged {n_ica_converge} of {n_converge_no_bootstrapping} times.   \n",)
        else:
            print(f"{int(100*(n_ica_converge/n_converge_no_bootstrapping))}% ", end = '')
       
    # 3: change data structure for sources, and compute similarities and distances between them.  
    A_hist = A_hist_BS + A_hist_no_BS                                                                   # list containing the time courses from each run.  i.e. each is: times x n_components
    S_hist = S_hist_BS + S_hist_no_BS                                                                   # list containing the soures from each run.  i.e.: each os n_components x n_pixels
    
    return S_hist, A_hist


     
#%%

def bootstrap_ICA(X, n_comp, bootstrap = True, ica_param = (1e-4, 150), 
                  X_whitened = None, dewhiten_matrix = None, verbose = True):
    """  A function to perform ICA either with or without boostrapping.  
    If not performing bootstrapping, performance can be imporoved by passing the whitened data and the dewhitening matrix
    (so that PCA does not have to be peroformed).  
    
    Inputs:
        X | rank2 array | data as row vectors (ie n_variables x n_samples)
        n_comp | int | number of sources to recover 
        X_whitened | rank2 array | data as row vectors (ie n_variables x n_samples), but whitened (useful if not bootstrapping)
        ica_param | tuple | Used to control ICA, (ica_tol, ica_maxit)
        X_whitened | rank2 array | data as row vectors (e.g. 10 x 20,000 for 10 ifgs of 20000 pixels), but whitened.  Useful to pass to function if not bootstapping as 
                                    this can then be calculated only once.  
        dewhiten_matrix | rank2 array | Converts the time courses recovered when using whitened data back to unwhiteend.  
                                        size is n_ifgs x n_sources.  
                                        X_white = A x S
                                        X = dewhiten x A x S
                                        Needed if not bootstrapping and don't want to  do PCA each time (as above)
        verbose | boolean | If True, the FastICA algorithm returns how many times it took to converge (or if it didn't converge)
    
    Returns:
        S | rank2 array | sources as row vectors (ie n_sources x n_samples)
        A | rank 2 array | time courses as columns (ie n_ifgs x n_sources)
        ica_success | boolean | True is the FastICA algorithm does converge.  
        
    History:
        2020/06/05 | MEG | Written
        2020/06/09 | MEG | Update to able to hand the case in which PCA fails (normally to do with finding the inverse of a matrix)
    
    """
    import numpy as np
    from licsalert.icasar.blind_signal_separation import PCA_meg2, fastica_MEG
    from licsalert.icasar.aux import  maps_tcs_rescale
    
    n_loop_max = 1000                                                               # when trying to make bootstrapped samples, if one can't be found after this many attempts, raise an error.  Best left high.  
    
    n_ifgs = X.shape[0]
    
    # 0: do the bootstrapping and determine if we need to do PCA
    if bootstrap:
        pca_needed = True                                                                                            # PCA will always be needed if bootstrapping 
        input_ifg_args = np.arange(n_comp-1)                                                                          # initiate as a crude way to get into the loop
        n_loop = 0                                                                                                   # to count how many goes it takes to generate a good bootstrap sample
        while len(np.unique(input_ifg_args)) < n_comp and n_loop < 100:                                              # try making a list of samples to bootstrap with providing it has enough unique items for subsequent pca to work
            input_ifg_args = np.random.randint(0, n_ifgs, n_ifgs)                                                    # generate indexes of samples to select for bootstrapping 
            n_loop += 1
        if n_loop == n_loop_max:                                                                                            # if we exited beacuse we were stuck in a loop, error message and stop
            raise Exception(f'Unable to bootstrap the data as the number of training data must be sufficently'
                            f' bigger than "n_components" sought that there are "n_components" unique items in'
                            f' a bootsrapped sample.  ')                                                             # error message
        X = X[input_ifg_args, :]                                                                              # bootstrapped smaple
    else:                                                                                                           # if we're not bootstrapping, need to work out if we actually need to do PCA
        if X_whitened is not None and dewhiten_matrix is not None:
            pca_needed = False
        else:
            pca_needed = True
            print(f"Even though bootstrapping is not being used, PCA is being performed.  "
                  f"This step could be sped up significantly by running PCA beforehand and "
                  f"computing 'X_whiten' and 'dewhiten_matrix' only once.  ")
        
    # 1 get whitened data using PCA, if we need to (ie if X_whitened and dewhiten_matrix aren't provided)
    if pca_needed:
        try:
            pca_vecs, _, _, dewhiten_matrix, _, _, X_whitened = PCA_meg2(X, verbose = False)                               # pca on bootstrapped data
            pca_success = True
        except:
            pca_success = False
    else:
        pca_success = True
    
    # If PCA was a success, do ICA (note, if not neeed, success is set to True)        
    if pca_success:                                                                                         
        X_whitened = X_whitened[:n_comp,]                                                                                                       # reduce dimensionality ready for ICA
        try:                                                                                                                                    # try ICA, as it can occasionaly fail (nans etc)
            W, S, A_white, _, _, ica_success = fastica_MEG(X_whitened, n_comp=n_comp,  
                                                           algorithm="parallel",        
                                                           whiten=False, maxit=ica_param[1],
                                                           tol = ica_param[0], 
                                                           verbose = verbose)         
            A = dewhiten_matrix[:,0:n_comp] @ A_white                                                                                       # turn ICA mixing matrix back into a time courses (ie dewhiten/ undo dimensonality reduction)
            S, A = maps_tcs_rescale(S, A)                                                                                                   # rescale so spatial maps have a range or 1 (so easy to compare)
        except:
            print(f"A FastICA run has failed, continuing anyway.  ")
            ica_success = False
        
    if pca_success and ica_success:
        return S, A, ica_success
    else:                                                                                                   # or if not a success, say that 
        ica_success = False
        return None, None, ica_success



#%%


def pairwise_comparison(sources_r2):
    """ Compte the pairwise distances and similarities for ICA sources.  
        Note that this uses the absolute value of the similarities, so is invariant to sign flips of the data.  
    Inputs:
        sources_r2 | rank 2 array | sources as row vectors
    """
    import numpy as np
    
    S = np.corrcoef(sources_r2)                                                  # Similarity matrix, Return Pearson product-moment correlation coefficients
    S = np.abs(S)                                                                # covariance of 1 and -1 are equivalent for our case
    D = 1 - S                                                                   # convert to dissimilarity    
    return D, S
    
 


#%%

def sources_list_to_r2_r3(sources, mask = None):
    """A function to convert a list of the outputs of multiple ICA runs (which are lists) into rank 2 and rank 3 arrays.  
    Inputs:
        sources | list | list of runs of ica (e.g. 10, or 20 etc.), each item would be n_sources x n_pixels
        mask | boolean | Only needed for two_d.  Converts row vector back to masked array.  
    Outputs:
        sources_r2 | rank 2 array | each source as a row vector (e.g. n_sources_total x n_pixels)
        sources_r3 | rank 3 masked array | each source as a rank 2 image. (e.g. n_souces_total x source_height x source_width )
    History:
        2018_06_29 | MEG | Written
        2020/08/27 | MEG | Update to handle both 1d and 2d signals.  
        2020/09/11 | MEG | Change sources_r3 so that it's now a masked array (sources_r3_ma)
       
    
    """
    import numpy as np
    import numpy.ma as ma
    from licsalert.aux import col_to_ma
    
    
    n_converge_needed = len(sources)
    n_comp = np.size(sources[0], axis = 0)
    n_pixels = np.size(sources[0], axis = 1)
    
    sources_r2 = np.zeros(((n_converge_needed * n_comp), n_pixels))                # convert from list to one big array
    for i in range(n_converge_needed):    
        sources_r2[i*n_comp:((i*n_comp) + n_comp), :] = sources[i]
    n_sources_total = np.size(sources_r2, axis = 0)

    if mask is not None:
        sources_r3 = ma.zeros((col_to_ma(sources_r2[0,:], mask).shape))[np.newaxis, :, :]               # get the size of one image (so rank 2)
        sources_r3 = ma.repeat(sources_r3, n_sources_total, axis = 0)                                    # and then extend to make rank 3
        for i in range(n_sources_total):
            sources_r3[i,:,:] = col_to_ma(sources_r2[i,:], mask)   
    else:
        sources_r3 = None
    return sources_r2, sources_r3



#%%

def cluster_quality_index(labels, S):
    """
    A function to calculate the cluster quality index (Iq).  If a cluster has only one element in it,
    the cluster quality index is set to nan (np.nan)
    Inputs:
        labels | rank 1 array | label number for each data point
        S | rank 2 array | similiarit between each data point 
    Returns:
        Iq | list | cluster quality index 
    2018_05_28 | written
    2018_05_30 | if clusters have only one point in them, set Iq to 0
    """
    import numpy as np
   
    Iq = []                                                                                         # initiate cluster quality index
    for i in np.unique(labels):                                                                     # loop through each label (there will be as many loops here as there are clusters)
        labels_1cluster = np.ravel(np.argwhere(labels == i))
        if np.size(labels_1cluster) < 2:                                                            # check if cluster has only one point in it
            Iq_temp = np.nan
        else:
            S_intra = np.copy(S[labels_1cluster, :][:,labels_1cluster])                                 # The similarties between the items in the cluster
            S_intra = np.where(np.eye(np.size(S_intra, axis = 0)) == 1, np.nan, S_intra)                # change the diagonals to nans
            S_inter = np.copy(S[labels_1cluster, :])                                                    # The similarties between the items in the cluster and those out of the cluster
            S_inter = np.delete(S_inter, labels_1cluster, axis = 1)                                     # horizontal axis remove similarities with itself
            Iq_temp = np.nanmean(S_intra) - np.mean(S_inter)                                         # Iq is the difference between the mean of the distances inside the cluster, and the mean distance between items in the cluster and out of the cluster
        Iq.append(Iq_temp)                                                                          # append whichever value of Iq (np.nan or a numeric value)
    return Iq

  
