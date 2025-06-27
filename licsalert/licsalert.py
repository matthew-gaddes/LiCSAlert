# -*- coding: utf-8 -*-
"""
A selection of functions used by LiCSAlert

@author: Matthew Gaddes
"""
import pdb
import matplotlib.pyplot as plt
 

#%% LiCSAlert()

def LiCSAlert(sources,
              mask_sources,
              inc_mc_space,
              baselines_cs,
              baseline_end,
              t_recalculate = 10,
              verbose=False, 
              n_pix_window = 20,
              residual_type = 'cumulative'
              ):
    
    """ Main LiCSAlert algorithm for a daisy-chain timeseries of interferograms.  
    
    Inputs:
        sources | r2 array | sources (from ICASAR) as row vectors, as per ICA, that can be turned back to interferograms with a rank 2 boolean mask of which pixels are masked and the col_to_ma function.  
        baselines_cs | r1 array | time values for each point in the time series, commonly (12,24,36) for Sentinel-1 data.  Could also be described as the cumulative temporal baselines.  
        ifgs_baseline | r2 array | ifgs used in baseline stage as row vectors
        mask | r2 array boolean | mask to convert a row vector (e.g. a row in sources or ifgs) back to a rank2 array.  
        ifgs_monitoring | r2 array or None | ifgs used in monitoring stage as row vectors
        t_recalculate | int | rolling lines of best fit are recalcaluted every X times (nb done in number of data points, not time between them)
        verbose | boolean | if True, various information is printed to screen.  
        n_pix_window | int | side length of square windows used for the window residual calculation.  
        residual_type | str | 'cumulative' or 'window'.  If cumulative, residual is the mean of the cumulative (ie summer in time for each pixel, then averaged in space across each ifg.  )
                                                        If window, the residual is the ratio of the max in a window (e.g. 20x20 pixels) over the mean for all windows for that time.  Value at each point is the cumulative.  
        
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
        
        d_hat            | rank 2 array | reconstructions (of incremrental ifgs)   
        d_resid         | rank 2 array | residual (d - d_hat)
    
    History:
        2019/12/XX | MEG |  Written from existing script.  
        2020/02/16 | MEG |  Update to work with no monitoring interferograms
        2020/12/14 | MEG | Improve docs and add option to save outputs.  
        2021_11_30 | MEG | Add new residual method (max residual of spatial window / mean of residual of all spatial windows for each ifg)
        2023_04_03 | MEG | also return residual and reconsructions
        2023_10_19 | MEG | remove option to save some products.  
    """
    import numpy as np
    import numpy.ma as ma
    import pickle
    from licsalert.licsalert import bss_components_inversion, tcs_baseline
    from licsalert.licsalert import tcs_monitoring, residual_for_pixels
    from licsalert.aux import update_mask
        
    # def check_input_sizes(baselines_cs, ifgs_all, sources):
    #     """ Check that there are the correct number of interferograms
    #     and the correct number of temporal baselines. 
    #     Also check that the number of pixels is correct between the sources
    #     and the interfergraoms.  
    #     """
    
    #     # check number of pixels
    #     if sources.shape[1] != ifgs_baseline.shape[1]:                                                                                              # 2nd dimension ([1] bit) is the number of pixels, which we compare.  
    #         raise Exception(f"The sources don't have the same number of pixels ({sources.shape[1]}) as the interferograms "
    #                         f"({ifgs_baseline.shape[1]}), so can't be used to fit them.  This is usually due to changing "
    #                         f"the cropped region but not re-running ICASAR.  Exiting...")
        
    #     # print(f"There are {baselines_cs.shape[0]} cumulative temporal baselines (set by baselines_cs), "
    #     #       f"and {ifgs_all.shape[0]} interferograms", end = '')
    #     if not baselines_cs.shape[0] != (ifgs_baseline.shape[0] + ifgs_monitoring.shape[0] + 1):                                                             # check in the case that there are monitoring ifgs
    #         # print(f", which is correct (as there are no daisy chain intefergraom on the "
    #         #       f"acquisition date.  ")
    #         pass
    #     else:
    #         print(f".  ")
    #         raise Exception(
    #             f"There appears to be a mismatch between the number of "
    #             f"cumulative baselines and the number of interferograms "
    #             f"As there is a cumulative temporal baselines on the first " 
    #             f"acquisition (value = 0) but no daisy chain ifg on that date, " +
    #             f"the baselines should be one longer than the number "
    #             f"of interferograms.  "
    #             )
        
        
    # function was designed to work with incremental so convert cumulative:
    # inc_ma = ma.diff(cum_ma, axis=0)
        
    # 0: Ensure we can still run LiCSAlert in the case that we have no 
    # monitoring interferograms (yet)
    # not needed now?
    
    # +1 still needed?  
    # n_times_baseline = baseline_end.acq_n + 1
    # if ifgs_monitoring is None:
    #     ifgs_all = np.copy(ifgs_baseline)                                                                # if there are no monitoring ifgs, ifgs_all is just the set of baseline ifgs
    #     n_times_monitoring = 0                                                                           # there are no monitoring ifgs
    # else:
    #     ifgs_all = np.vstack((ifgs_baseline, ifgs_monitoring))                                           # ifgs are row vectors, so stack vertically
    #     n_times_monitoring = ifgs_monitoring.shape[0]
    # check_input_sizes(baselines_cs, ifgs_all, sources)
        
    # do an inversion to fit the incremental ifgs using the sources, then
    # make the time courses (inversion results) cumulative, starting at 0
    # d_hat is reconsturction, d_resid is data - reconstrutcion.  
    
    
    
    # tcs)_c : time courses as column vectors (i.e. n_times x n_ics)
    # d_hat, d_resid.  Same as cum_ma
    # NOTE that these are not cumulative (as working with incremental data
    # and cumualtive is set to False.  
    # also note that we only use the baseline data here.  
    # time courses (tcs) should be the same length as incremental 
    # input ifgs inc_ma[...] as cumulative is set to false
    tcs, d_hat, d_resid = bss_components_inversion_per_epoch(
        sources,
        mask_sources,
        inc_mc_space[:baseline_end.acq_n, ],
        cumulative=False
    )
    
    # create cumulative time courses (i.e. integrate in time):
    # note that tcs_c_baseline is now 1 longer than incremental input ifgs.  
    tcs_c_baseline = np.concatenate(
        (np.zeros(
            (1, tcs.shape[1])
            ),
         np.cumsum(
             tcs, axis=0)
         ),
        axis = 0
    )         
    del tcs
    
    # f, ax = plt.subplots()
    # ax.matshow(d_hat[10,])
    # # debug plot
    # f, axes = plt.subplots(tcs_cs.shape[1], 1)
    # for ax_n, tc in enumerate(tcs_cs.T):
    #     axes[ax_n].plot(tc)
    # # end 
        
    # calcaulte line of best fit, gradient, point to line distances etc. for 
    # cumulative time courses.  
    sources_tcs = tcs_baseline(
        tcs_c_baseline, 
        baselines_cs[:baseline_end.acq_n+1],
        t_recalculate,
    )                                   
   
    
    pdb.set_trace()
    
    # calculate the residual, which is summed for each pixel, then averaged 
    # across an image at any time.  
    
    # this needs to be similar to how it was caldulated in bss_inversion?  
    # ie per epoch?  Can we just use this output to then calculte Andy's window approach etc? 
    
    residual = residual_for_pixels(
        sources, sources_tcs, ifgs_baseline, mask, n_pix_window, residual_type
        )             
    
    # calculate line of best fit, gradient, line to point distances etc for 
    # cumulative residual.  
    residual_tcs = tcs_baseline(
        residual, baselines_cs[:baseline_end.acq_n], t_recalculate
        )                               
    del tcs_c_baseline, residual

    #2: Calculate time courses/distances etc for the monitoring data
    if ifgs_monitoring is not None:

        # do the inversion to fit the monitoring ifgs with the sources    
        # compute cumulative time courses for monitoring interferograms (ie
        # simple inversion to fit each ifg in ifgs_baseline using sources.  )
        products = bss_components_inversion(
            sources, ifgs_monitoring,  cumulative=True, mask = mask
            )                      
        (tcs_c_monitoring, d_hat_monitoring, d_resid_monitoring) = products
        del products
        
        # remove the 0s on the first row as this isn't the basleine stage.  
        tcs_c_monitoring = np.delete(tcs_c_monitoring, 0,0)
        sources_tcs_monitor = tcs_monitoring(
            tcs_c_monitoring, sources_tcs, baselines_cs
            )                               
    
        #3: and update the residual stuff                                                                            
        # get the cumulative residual for baseline and monitoring (hence _cb)    
        residual_bm = residual_for_pixels(
            sources, sources_tcs_monitor, ifgs_all, mask, n_pix_window, 
            residual_type
            )                               
        residual_tcs_monitor = tcs_monitoring(
            residual_bm, residual_tcs, baselines_cs, residual=True
            )               
    
    # 3: combine the reconsuction and residual for all ifgs
    d_hat = (np.hstack((d_hat_baseline, d_hat_monitoring))).T                                                                           # reconstruction, n_times x n_pixels 
    d_resid = (np.hstack((d_resid_baseline, d_resid_monitoring))).T                                                                     # residual

    if ifgs_monitoring is None:
        return sources_tcs, residual_tcs
    else:
        return sources_tcs_monitor, residual_tcs_monitor, d_hat, d_resid
    

#%% LiCSAlert() backup 2025_06_25


# def LiCSAlert(sources,
#               baselines_cs,
#               ifgs_baseline,
#               mask,
#               ifgs_monitoring = None,
#               t_recalculate = 10,
#               verbose=False, 
#               n_pix_window = 20,
#               residual_type = 'cumulative'
#               ):
    
#     """ Main LiCSAlert algorithm for a daisy-chain timeseries of interferograms.  
    
#     Inputs:
#         sources | r2 array | sources (from ICASAR) as row vectors, as per ICA, that can be turned back to interferograms with a rank 2 boolean mask of which pixels are masked and the col_to_ma function.  
#         baselines_cs | r1 array | time values for each point in the time series, commonly (12,24,36) for Sentinel-1 data.  Could also be described as the cumulative temporal baselines.  
#         ifgs_baseline | r2 array | ifgs used in baseline stage as row vectors
#         mask | r2 array boolean | mask to convert a row vector (e.g. a row in sources or ifgs) back to a rank2 array.  
#         ifgs_monitoring | r2 array or None | ifgs used in monitoring stage as row vectors
#         t_recalculate | int | rolling lines of best fit are recalcaluted every X times (nb done in number of data points, not time between them)
#         verbose | boolean | if True, various information is printed to screen.  
#         n_pix_window | int | side length of square windows used for the window residual calculation.  
#         residual_type | str | 'cumulative' or 'window'.  If cumulative, residual is the mean of the cumulative (ie summer in time for each pixel, then averaged in space across each ifg.  )
#                                                         If window, the residual is the ratio of the max in a window (e.g. 20x20 pixels) over the mean for all windows for that time.  Value at each point is the cumulative.  
        
#     Outputs
#         sources_tcs_monitor | list of dicts | list, with item for each time course.  Each dictionary contains:
#         cumualtive_tc | cumulative time course : the cumulative sum of how that IC is used to fit each interferogram.  Column vector.  
#         gradient      | gradient of the cumulative time course in the baseline stage.  float.  
#         lines         | n_times x n_times.  The lines of best fit used to calculate the line-of-best-fit to point distance at each time.  
#                         Each column is a line of best fit, but many values are zeros as the lines of best fit are redrawn periodically.  
#         sigma         | sigma (standard deviation) of the line-to-point distances in the baseline stage
#         distances     | number of standard deviations each point is from it's line of best fit.  
#         t_recalcaulte | integer, how often the lines of best fit are redrawn.  
                                                                
#         residual_tcs_monitor | list of dicts | As per above, but only for the cumulative residual (i.e. list is length 1)
        
#         d_hat            | rank 2 array | reconstructions (of incremrental ifgs)   
#         d_resid         | rank 2 array | residual (d - d_hat)
    
#     History:
#         2019/12/XX | MEG |  Written from existing script.  
#         2020/02/16 | MEG |  Update to work with no monitoring interferograms
#         2020/12/14 | MEG | Improve docs and add option to save outputs.  
#         2021_11_30 | MEG | Add new residual method (max residual of spatial window / mean of residual of all spatial windows for each ifg)
#         2023_04_03 | MEG | also return residual and reconsructions
#         2023_10_19 | MEG | remove option to save some products.  
#     """
#     import numpy as np
#     import pickle
#     from licsalert.licsalert import bss_components_inversion, residual_for_pixels, tcs_baseline, tcs_monitoring  
        
#     def check_input_sizes(baselines_cs, ifgs_all, sources):
#         """ Check that there are the correct number of interferograms
#         and the correct number of temporal baselines. 
#         Also check that the number of pixels is correct between the sources
#         and the interfergraoms.  
#         """
    
#         # check number of pixels
#         if sources.shape[1] != ifgs_baseline.shape[1]:                                                                                              # 2nd dimension ([1] bit) is the number of pixels, which we compare.  
#             raise Exception(f"The sources don't have the same number of pixels ({sources.shape[1]}) as the interferograms "
#                             f"({ifgs_baseline.shape[1]}), so can't be used to fit them.  This is usually due to changing "
#                             f"the cropped region but not re-running ICASAR.  Exiting...")
        
#         # print(f"There are {baselines_cs.shape[0]} cumulative temporal baselines (set by baselines_cs), "
#         #       f"and {ifgs_all.shape[0]} interferograms", end = '')
#         if not baselines_cs.shape[0] != (ifgs_baseline.shape[0] + ifgs_monitoring.shape[0] + 1):                                                             # check in the case that there are monitoring ifgs
#             # print(f", which is correct (as there are no daisy chain intefergraom on the "
#             #       f"acquisition date.  ")
#             pass
#         else:
#             print(f".  ")
#             raise Exception(
#                 f"There appears to be a mismatch between the number of "
#                 f"cumulative baselines and the number of interferograms "
#                 f"As there is a cumulative temporal baselines on the first " 
#                 f"acquisition (value = 0) but no daisy chain ifg on that date, " +
#                 f"the baselines should be one longer than the number "
#                 f"of interferograms.  "
#                 )
        
#     # 0: Ensure we can still run LiCSAlert in the case that we have no 
#     # monitoring interferograms (yet)
#     n_times_baseline = ifgs_baseline.shape[0] + 1                                                   
#     if ifgs_monitoring is None:
#         ifgs_all = np.copy(ifgs_baseline)                                                                # if there are no monitoring ifgs, ifgs_all is just the set of baseline ifgs
#         n_times_monitoring = 0                                                                           # there are no monitoring ifgs
#     else:
#         ifgs_all = np.vstack((ifgs_baseline, ifgs_monitoring))                                           # ifgs are row vectors, so stack vertically
#         n_times_monitoring = ifgs_monitoring.shape[0]
#     check_input_sizes(baselines_cs, ifgs_all, sources)

#     # # pdb.set_trace()
#     # # start debug
#     # plt.switch_backend('qt5agg')
#     # # for i in range(4):
#     # #     f, ax = plt.subplots()
#     # #     ax.plot(tcs_c[:,i])
#     # #     plt.pause(2)
#     # from licsalert.aux import col_to_ma
#     # f, axes = plt.subplots(1,2)
#     # axes[0].matshow(col_to_ma(ifgs_baseline[0,:], mask))
#     # axes[1].matshow(col_to_ma(sources[3,:], mask))
#     # # end
        
#     # do an inversion to fit the incremental ifgs using the sources, then
#     # make the time courses (inversion results) cumulative, starting at 0
#     # d_hat is reconsturction, d_resid is data - reconstrutcion.  
#     products = bss_components_inversion(
#         sources, ifgs_baseline, cumulative=True, mask = mask
#         )     
#     tcs_c_baseline, d_hat_baseline, d_resid_baseline = products; del products
    
#     # calcaulte line of best fit, gradient, point to line distances etc. for 
#     # cumulative time courses.  
#     sources_tcs = tcs_baseline(
#         tcs_c_baseline, baselines_cs[:n_times_baseline], t_recalculate
#         )                                   
   
#     # calculate the residual, which is summed for each pixel, then averaged 
#     # across an image at any time.  
#     residual = residual_for_pixels(
#         sources, sources_tcs, ifgs_baseline, mask, n_pix_window, residual_type
#         )             
    
#     # calculate line of best fit, gradient, line to point distances etc for 
#     # cumulative residual.  
#     residual_tcs = tcs_baseline(
#         residual, baselines_cs[:n_times_baseline], t_recalculate
#         )                               
#     del tcs_c_baseline, residual

#     #2: Calculate time courses/distances etc for the monitoring data
#     if ifgs_monitoring is not None:

#         # do the inversion to fit the monitoring ifgs with the sources    
#         # compute cumulative time courses for monitoring interferograms (ie
#         # simple inversion to fit each ifg in ifgs_baseline using sources.  )
#         products = bss_components_inversion(
#             sources, ifgs_monitoring,  cumulative=True, mask = mask
#             )                      
#         (tcs_c_monitoring, d_hat_monitoring, d_resid_monitoring) = products
#         del products
        
#         # remove the 0s on the first row as this isn't the basleine stage.  
#         tcs_c_monitoring = np.delete(tcs_c_monitoring, 0,0)
#         sources_tcs_monitor = tcs_monitoring(
#             tcs_c_monitoring, sources_tcs, baselines_cs
#             )                               
    
#         #3: and update the residual stuff                                                                            
#         # get the cumulative residual for baseline and monitoring (hence _cb)    
#         residual_bm = residual_for_pixels(
#             sources, sources_tcs_monitor, ifgs_all, mask, n_pix_window, 
#             residual_type
#             )                               
#         residual_tcs_monitor = tcs_monitoring(
#             residual_bm, residual_tcs, baselines_cs, residual=True
#             )               
    
#     # 3: combine the reconsuction and residual for all ifgs
#     d_hat = (np.hstack((d_hat_baseline, d_hat_monitoring))).T                                                                           # reconstruction, n_times x n_pixels 
#     d_resid = (np.hstack((d_resid_baseline, d_resid_monitoring))).T                                                                     # residual

#     if ifgs_monitoring is None:
#         return sources_tcs, residual_tcs
#     else:
#         return sources_tcs_monitor, residual_tcs_monitor, d_hat, d_resid


#%% construct_baseline_ts()

def construct_baseline_ts(sica_tica, 
                          displacement_r3, 
                          tbaseline_info,
                          baseline_end,
                          interactive=False
):
    """
    A function to prepare the baseline time series for ICA.  First step
    is dropping some epochs that cause many pixels to be lost, followed by
    mean centering 
    """
    import numpy as np
    
    from licsalert.pixel_selection import automatic_pixel_epoch_selection
    from licsalert.temporal import daisy_chain_from_acquisitions
    from licsalert.data_importing import ifg_timeseries
    
    # 1: select a subset of the epochs to create a compromise between lots of 
    # pixels and lots of epochs (i.e. drop the epochs that cause lots of 
    # pixels to be lost when a consistent mask through time is sought)
    displacement_r2_ica, tbaseline_info_ica = automatic_pixel_epoch_selection(
        displacement_r3,
        tbaseline_info,
        baseline_end,
        interactive
        )
    
    # also copy some auxilliary info to the new array
    for key in ['dem', 'lons_mg', 'lats_mg']:
        displacement_r2_ica[key] = displacement_r3[key]
    
    # rename the output to be more readable
    displacement_r2_ica['cumulative'] = displacement_r2_ica['ifgs']
    
    # create the incremental displacments
    displacement_r2_ica['incremental'] = np.diff(
        displacement_r2_ica['cumulative'],
        axis = 0
        )
    
    # calculate the dates of the incremental ifgs
    tbaseline_info_ica['ifg_dates'] = daisy_chain_from_acquisitions(
        tbaseline_info_ica['acq_dates']
    )
        
    
    # 2: mean centre in time or space, according to if sica or tica 
    # create a class (an ifg_timeseries) using the incremtnal measurements, 
    # and handle all mean centering (in time and space)
    
    ifg_ts = ifg_timeseries(
        displacement_r2_ica["incremental"],
        tbaseline_info_ica['ifg_dates']
    )
    
    displacement_r2_ica['incremental_mc_space'] = ifg_ts.mixtures_mc_space
    displacement_r2_ica['means_space'] = ifg_ts.means_space
    
    displacement_r2_ica['incremental_mc_time'] = ifg_ts.mixtures_mc_time
    displacement_r2_ica['means_time'] = ifg_ts.means_time

    if sica_tica == 'sica':
        displacement_r2_ica['mixtures_mc'] = ifg_ts.mixtures_mc_space
        displacement_r2_ica['means'] = ifg_ts.means_space
    elif sica_tica == 'tica':
        displacement_r2_ica['mixtures_mc'] = ifg_ts.mixtures_mc_time
        displacement_r2_ica['means'] = ifg_ts.means_time
    else:
        raise Exception(
            f"'sica_tica' can be either 'sica' or 'tica', but not "
            f"{sica_tica}.  Exiting.  "
            )
        
    # check sizes of things
    n_acq = len(tbaseline_info_ica['acq_dates'])
    
    print(
        f"\nChecking the size of the baseline data that has been created for"
        f"use with ICA:  ")
    
    for key, value in tbaseline_info_ica.items():
        print(f"    {key} : {len(value)}")
        
    for key, value in displacement_r2_ica.items():
        try:
            print(f"    {key} : {value.shape}")
        except:
            pass
   
    print(f"\n\n")
    
    return displacement_r2_ica, tbaseline_info_ica
    
    
    


#%% residual_for_pixels()

def residual_for_pixels(sources, sources_tcs, ifgs, mask, n_pix_window = 20, 
                        residual_type = 'cumulative', n_skip=None):
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
        mask | r2 array of boolean | to convert a row vector back to a rank 2 masked array.  
        n_skip | None or int | if an int, the first n_skip values of the timecourses will be skipped.  
        n_pix_window | int | side length of square windows used for the window residual calculation.  
        residual_type | str | 'cumulative' or 'window'.  If cumulative, residual is the mean of the cumulative (ie summer in time for each pixel, then averaged in space across each ifg.  )
                                                        If window, the residual is the ratio of the max in a window (e.g. 20x20 pixels) over the mean for all windows for that time.  Value at each point is the cumulative.  

    Outputs:
        residual_ts | r2 array | Column vector of the RMS residual between that ifg, and its reconstruction
        residual_cs | r2 array | Column vector of the RMS residual between an ifg and the cumulative residual (for each pixel)
                                 N.b. the point is that if we have a strong atmosphere, it then reverses in the next ifg
                                 so the cumulative for each pixel goes back to zero

    2019/01/XX | MEG | Written, in discussion with AH
    2019/12/06 | MEG | Comment and documentation
    2020/01/02 | MEG | Update to use new LiCSAlert list of dictionaries
    2020/02/06 | MEG | Fix bug as had forgotten to convert cumulative time courses to be incremental
    2021_11_30 | MEG | Add new residual method (ratio of maximum in window for mean of all windows)
    """

    import numpy as np
    import numpy.ma as ma
    from licsalert.aux import r2_to_r3

    def list_dict_to_r2(sources_tcs):
        """ Extract the cumulative time courses from the sources_tcs list of dicts, and convert them to incremental
        time courses.  
        """
        
        n_sources = len(sources_tcs)
        tcs_c_length = sources_tcs[0]["cumulative_tc"].shape[0]                               # as many rows as time steps
        tcs = np.zeros((tcs_c_length - 1, n_sources))                                          # initiate
        for n_source, source_tc in enumerate(sources_tcs):                              # loop through each source
            # convert the cumulative tc to the incremental.  
            tcs[:,n_source:n_source+1] = np.diff(source_tc["cumulative_tc"], axis = 0)                    # convert to incremental time course, add 0 at start to ensure that we get something of the same length.  
        return tcs
                    

    (n_sources, n_pixs) = sources.shape                                         # number of sources and number of pixels
    tcs = list_dict_to_r2(sources_tcs)                                          # get the incremental time courses as a rank 2 array (ie. as column vectors, so something like 55x5, if there are 55 times and 5 sources)
    if n_skip is not None:                                                      # crop/remove the first ifgs
        tcs = tcs[n_skip:,]                                                        # usually the baseline ifgs when used with monitoring data
    

    data_model_residual = ifgs - (tcs @ sources)                                                # residual for each pixel at each time (ie the data - the reconstrutcion)
    data_model_residual_cs = np.cumsum(data_model_residual, axis = 0)                           # summing the residual for each pixel cumulatively through time, same size as above.     
    residual_ts = np.zeros((data_model_residual.shape[0], 1))                                   # initiate, n_ifgs x 1 array
    residual_cs = np.zeros((data_model_residual.shape[0], 1))                                   # initiate, n_ifgs x 1 array
    for row_n in range(data_model_residual.shape[0]):                                           # loop through each ifg
        residual_ts[row_n, 0] = np.sqrt(np.sum(data_model_residual[row_n,:]**2)/n_pixs)         # RMS of residual for each ifg
        residual_cs[row_n, 0] = np.sqrt(np.sum(data_model_residual_cs[row_n,:]**2)/n_pixs)      # RMS of residual for cumulative ifgs.  i.e. if atmosphere reverses, should cancel in 2nd cumulative ifg.  
        
    # method 2: calculate the residual in spatial windows.  
    residuals_cs_r3 = r2_to_r3(data_model_residual_cs, mask)                                        # convert from row vectors to rank 2 images.  n_times x ny x nx
    (n_residuals, ny, nx) = residuals_cs_r3.shape
    n_windows_x = int(np.floor(nx / n_pix_window))                                                  # the number of whole windows we can it in the x direction.  
    rem_x = nx % (n_pix_window * n_windows_x)                                                       # the number of pixels remaining if we do that.  
    border_x = int(rem_x / 2)                                                                       # what half the border will be
    n_windows_y = int(np.floor(ny / n_pix_window))                                                  # as above but for y
    rem_y = ny % (n_pix_window * n_windows_y)
    border_y = int(rem_y / 2)
    
    residual_cs_window = np.zeros((data_model_residual.shape[0], 1))                                   # initiate, n_ifgs x 1 array, to store the ratio of the max residual to the mean residual for each ifg. 
    residual_cs_windows_r3 = np.zeros((n_residuals, n_windows_y, n_windows_x))                                   # initiate, n_ifgs x 1 array
    for n_resid, residual_cs_r3 in enumerate(residuals_cs_r3):
        for row_n in range(n_windows_y):
            for col_n in range(n_windows_x):
                window_region = residuals_cs_r3[n_resid, (border_y + row_n*n_pix_window):(border_y + (row_n+1)*n_pix_window), 
                                                         (border_x + col_n*n_pix_window):(border_x + (col_n+1)*n_pix_window) ]                 # extract the window region
                n_coherent_pixs = (n_pix_window**2) - np.sum(window_region.mask)                                                               # get the number of coherent pixels in the windo (maksed pixels are 1)
                if False not in window_region.mask:                                                                                             # if there are no unmasked pixels (unmasked = False)
                    residual_cs_windows_r3[n_resid, row_n, col_n] = np.nan                                                                      # then the residual for that window is just nan
                else:
                    residual_cs_windows_r3[n_resid, row_n, col_n] = np.sqrt(ma.sum(window_region**2)/n_coherent_pixs)                           # ge the RMS error for the window, and assign to the array that stores them.  
        
        residual_cs_window[n_resid] = np.nanmax(residual_cs_windows_r3[n_resid,]) / np.nanmean(residual_cs_windows_r3[n_resid,])                # calculate the ratio between the largest and mean residual 
        #residual_cs_window[n_resid] = np.nanmax(residual_cs_windows_r3[n_resid,]) / np.nanmin(residual_cs_windows_r3[n_resid,])                # calculate the ratio between the largest and minimum residual
    # end method 2
    

    # debugging/ examining plot
    # import numpy.ma as ma
    # import matplotlib.pyplot as plt
    # n_plot = 15
    # f, axes = plt.subplots(2,n_plot, figsize = (30,7))
    # residuals = residuals_cs_r3[:n_plot,]
    # window_residuals = residual_cs_windows_r3[:n_plot,]
    # for i in range(n_plot):
    #     #axes[0,i].imshow(residuals[i,], vmin = ma.min(residuals), vmax = ma.max(residuals))                                        # share colorbar
    #     axes[0,i].imshow(residuals[i,])                                                                                             # dont 
    #     axes[1,i].imshow(window_residuals[i,], vmin = np.nanmin(window_residuals), vmax = np.nanmax(window_residuals))              # share colorbar
    
    # determine which residual to return, and add 0 to first row (as 0 residual on first acquisition date)
    if residual_type == 'cumulative':
        residual = residual_cs
    elif residual_type == 'window':
        residual = residual_cs_window
    else:
        raise Exception(f"'residual_type' must be either 'cumulative' or 'window'.  Exiting as it's {residual_type}")
        
    residual = np.concatenate((np.zeros((1, residual_cs.shape[1])), residual_cs), axis = 0)
        
    return residual


#%%

def tcs_baseline(tcs_c, baselines_cs, t_recalculate):
    """
    Given cumulative time courses (tsc_c), the time values for each entry (baselines_cs), and the recalculation time 
    (t_recalculate), create a list with a dictionary about each time course.  The dictionary contains the cumulative 
    time courses, their gradient, and the ditsances each point is from a line of best fit at that gradient, and the 
    lines of best fit, redrwwn every t_recalculate.  
    
    Inputs:
        tcs_c | rank2 array | Cumulative Time CourseS, as column vectors
        baselines_cs  | r1 array | time values for each point in the time course.  for Sentinel-1, commonly (12,24,36 etc)
        t_recalculate | int | the number of ifgs used to calculate the rolling lines of best fit
    
    Outputs:
        sources_tcs | list of dicts | as per description.  
        
    History:
        2020/01/02 | MEG | Written
    """
    import numpy as np
    
    n_acqs, n_tcs = tcs_c.shape                                                                    # there will be as many time courses (tcs) as there are sources, which are rows
    sources_tcs = []                                                                           # each time course will be a list
    for n_tc in range(n_tcs):
        tc_dict = {}                                                                           # initialise
        # 1: add the cumulative time course
        tc_dict["cumulative_tc"] = tcs_c[:, n_tc:n_tc+1]                                       # keep cumulative time course as a column vector
        
        # 2: add the gradient 
        gradient, y_intercept = np.polyfit(baselines_cs, tc_dict["cumulative_tc"], 1)            # gradients 1st, y intercept second        
        tc_dict["gradient"] = gradient[0]                                                       # add to dictionary
        
        # 3: generate line of best fit for whole baseline stage. 
        line_yvals = np.polyval((tc_dict["gradient"], y_intercept[0]), baselines_cs)                

        tc_dict["lines"] = []
        for n_acq in range(n_acqs):                                                        # loop through each time step, to calculate it
            # get the start acq of the line of best fit, can't be less than 0
            start_acq = np.max((0, (n_acq+1) - (t_recalculate)))
            # get the y values for the line of best fit for this acquisition (it and previous ones)
            tc_dict["lines"].append(line_yvals[start_acq:n_acq+1])              # copy the y values for that little bit of line
                
        # 4: line to point distances (which are stored in the dict in terms of how many sigmas they are)
        line_point_distances = tc_dict["cumulative_tc"] - line_yvals[:,np.newaxis]
        tc_dict["sigma"] = np.std(line_point_distances)
        tc_dict["distances"] = np.abs(line_point_distances / tc_dict["sigma"])                  # ie the number of standard deviations a point is from the line of best fit
        
        # 5: Record the t_recalculate parameter
        tc_dict["t_recalculate"] = t_recalculate
        
        sources_tcs.append(tc_dict)                                                            # append to list, where each item is a time_course dictionary
    return sources_tcs
    
#%%    
    
def tcs_monitoring(tcs_c, sources_tcs, baselines_cs, residual=False):
    """
    Given an extension to a time series (i.e. in monitoring mode) and given the time courses for the baseline
    data (sources_tcs), create a new list of dicts for the complete data (i.e. rolling lines of best fit, 
    line to point distances).  
    
    Inputs:
        tcs_cs | rank 2 array | time courses as column vectors
        sources_tcs | list of dicts | 
        baselines_cs  | r1 array | time values for each point in the time course.  
                                    for Sentinel-1, commonly (12,24,36 etc if 12 day acquisitions.  )
        residual | boolean | 
        
    Outputs:
        sources_tcs | list of dicts | as per description.  
        
    History:
        2020/01/02 | MEG | Written
    """
    import numpy as np
    import copy                                                                                   # needed for the deep copy of the list of dictionaries
    
    
    # cumulative baselines originally started without 0.  Reformat back to starting
    # without 0.  
    #baselines_cs = baselines_cs[1:]
    
    # 1: Small initial steps
    sources_tcs_monitor = copy.deepcopy(sources_tcs)                                              # copy so that we can change without affecting version calcaulted only from baseline data
    # get from the baseline stage
    t_recalculate = sources_tcs_monitor[0]["t_recalculate"]                                       
    # get the number of monitoring ifgs.  
    n_acqs_baseline = sources_tcs_monitor[0]["cumulative_tc"].shape[0] 
    
    n_acqs = len(baselines_cs)
    
    if residual:
        n_acqs_monitor = n_acqs - n_acqs_baseline
    else:
        n_acqs_monitor = tcs_c.shape[0]                                                              # number of time steps is just the number of interferograms
        if n_acqs != (n_acqs_baseline + n_acqs_monitor):
            raise Exception(
                f"There is a mismatch in the number of acquisitions.  "
                f"The number of acquistions in the baseline stage "
                f"plus the number of acquisitions in the monitoring stage "
                f"should equal the number of acquisitions in the cumulative "
                f"temporal baselines, but doesn't.  Exiting.  "
            )
    
    
    
    
    
    # why do we do this? 
    # if residual:
    #     n_times_monitor = n_times_tc - n_times_baseline
    # else:
    #     n_times_monitor = n_times_tc
    # n_times_total = n_times_baseline + n_times_monitor       
    
    for source_n, source_tc in enumerate(sources_tcs_monitor):
        # 1: Update the cumulative time course 
        if residual is False:                                                                                   # residual timecourse doesn't start from 0, so don't need to add value to it
            # make the monitoring cumulative time course start at the value the basline time course ends at 
            tcs_c[:, source_n] += source_tc["cumulative_tc"][-1]                                                
            # and join them
            source_tc["cumulative_tc"] = np.vstack((source_tc["cumulative_tc"], tcs_c[:, source_n][:, np.newaxis]))        
        else:
            source_tc["cumulative_tc"] = tcs_c[:, source_n][:, np.newaxis]                                    # don't need to extend as have full time course
    
        
        # 2: Lines of best fit, and line to point distances for the monitoring data
        # extend with 0s to length of data now.  
        source_tc["distances"] = np.pad(source_tc["distances"], [(0, n_acqs_monitor), (0,0)], "constant", constant_values=(0))                 
        
        # loop through the monitoring ifgs, 
        for n_acq in np.arange(n_acqs_baseline, n_acqs):                                                                                
            # get the start acquisition for the line of best fit, don't let be negative for first few lines (i.e. line 10 long for acq 6)
            start_acq = np.max((0, n_acq-t_recalculate))
            # get the y interecept for the current line of best fit.  
            line_y_intercept = (np.mean(source_tc["cumulative_tc"][start_acq : n_acq]) 
                             - (source_tc["gradient"]*np.mean(baselines_cs[start_acq : n_acq])))                  

    
            # shift the line of best fit one to the right 
            if start_acq == 0:
                start_acq2 = start_acq
            else:
                start_acq2 = start_acq + 1
                
            line_yvals = (baselines_cs[start_acq2: n_acq+1] * source_tc["gradient"]) + line_y_intercept                                                                      # predict the y values given the gradient and y-intercept of the line, note that also predicting y value for next point
            # and record the line of best fit

            source_tc["lines"].append(line_yvals)
            # and the distance between the line and the point (expressed in number of sigmas).
            source_tc["distances"][n_acq,] = (np.abs(source_tc["cumulative_tc"][n_acq,] - line_yvals[-1]))/source_tc["sigma"]      
    return sources_tcs_monitor
    
# maybe this could all be replaced?  
# for basline, calcualte gradient, and line to point distances.
# for monitoring, just use the same function, but give it the gradient and apply to whole time series?  why the complexity of a second funciton? 

#%%
        
        
def LiCSBAS_for_LiCSAlert(
        LiCSAR_frame,
        LiCSAR_frames_dir,
        LiCSBAS_out_dir,
        logfile_dir,
        LiCSBAS_bin,
        lon_lat = None,
        downsampling = 1,
        n_para=1
    ):
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
    

#%% LiCSAlert_preprocessing()

def LiCSAlert_preprocessing(
        displacement_r3,
        tbaseline_info, 
        sica_tica,
        downsample_run=1.0,
        downsample_plot=0.5,
        verbose=True,
        order = 1,
        anti_aliasing=False
    ):
    """A function to downsample the data at two scales (one for general working [ie to speed things up], and one 
    for faster plotting.  )  Also, data are mean centered, which is required for ICASAR and LiCSAlert.  
    Note that the downsamples are applied consecutively, so are compound (e.g. if both are 0.5, 
    the plotted data will be at 0.25 the resolution of the original data).  
    
    Inputs:
        displacement_r2 | dict | input data stored in a dict as row vectors with a mask    
                                 Also lons and lats as rank 2 arrays.  E N U are components of the look vector for East North Up for each pixel.  
                                 dict_keys(['dem', 'mask', 'incremental', 'lons', 'lats', 'E', 'N', 'U'])
        tbaseline_info | dict |  dict_keys(['acq_dates', 'ifg_dates', 'baselines', 'baselines_cumulative'])
        sica_tica       | string | if sica, spatial ICA, if tica, temporal ica
        downsample_run | float | in range [0 1], and used to downsample the "incremental" data
        downsample_plot | float | in range [0 1] and used to downsample the data again for the "incremental_downsample" data
        verbose | boolean | if True, informatin returned to terminal
        
        order : int, default 1 | Interpolation order passed to `rescale` 
                                    (1 = bilinear, 0 = nearest, etc.).
        anti_aliasing : bool, default False Forwarded to `rescale`.  Set True 
                                    for images, False for scientific rasters.
        
    Outputs:
        displacement_r2 | dict | input data stored in a dict as row vectors with a mask
                                 updated so that "incremental" is downsampled (and its mask), 
                                 and a new key is created, called "incremental_downsampled" 
                                 that is downsamled further for fast plotting                                 
                                 
                                 'dem' - the DEM, downsampled.  
                                 'mask' - mask, downsampled.  
                                 'incremental' - daisy chain of short temporal baseline ifgs, downsampled but not mean centered.  
                                 'lons' - longitude for each pixel, downsampled.  
                                 'lats' - latitude for each pixel, downsmapled. 
                                 'E' - east component of look vector, downsmapled.  
                                 'N' - north component of look vector, downsmapled.
                                 'U' - up componennt of look vector, downsampled.  
                                 'incremental_downsampled' - as above, but downsampled again (so downsampled_run * downsample_plot)
                                 'mask_downsampled' - mask for double downsamled data.  
                                 
                                 'incremental_mc_space' - as for incremental, but mean cenetered in space.  
                                 'means_space' - to undo mean centering in space.  
                                 'incremental_mc_time' - as for incremental, but mean centered in time. 
                                 'means_time' - to undo mean centering in time.  
                                 
                                 'mixtures_mc' - mean centered in space or time, depending on how sica_tica was set.  
                                 'means' - to undo above mean centering.  
                                 
                                 
    History:
        2020/01/13 | MEG | Written
        2020/12/15 | MEG | Update to also downsample the lons and lats in the ICASAR geocoding information.  
        2021_04_14 | MEG | Update so handle rank2 arrays of lons and lats properly.  
        2021_05_05 | MEG | Add check that lats are always the right way up, and fix bug in lons.  
        2021_10_13 | MEG | Add function to also downsample ENU grids.  
        2021_11_15 | MEG | Add option to control mean centering, and warning that it is happening.  
        2023_10_26 | MEG | use ifg_timeseries class to handle mean centering more carefully.  
        
        2025_06_22 | MEG | WIP - had to remove mean centerering stuff 
    """
    import numpy as np
    import numpy.ma as ma
    # from licsalert.downsample_ifgs import downsample_ifgs
    from licsalert.aux import col_to_ma
    from licsalert.data_importing import ifg_timeseries
    
    from skimage.transform import rescale

    # class MeanCenteredArray:
    #     def __init__(self, masked_array):
    #         if not isinstance(masked_array, np.ma.MaskedArray):
    #             raise TypeError("Input must be a masked array")
    #         # interal/private attriubute
    #         self._arr = masked_array
    #         self.mean_centered = _MeanCenteredAccessor(self._arr)
    
    #     @property
    #     def original(self):
    #         return self._arr
    
    
    # class _MeanCenteredAccessor:
    #     def __init__(self, arr):
    #         self._original = arr
    
    #         # Precompute mean-centered versions
    #         # Mean-centre in space (each image mean = 0)
    #         image_means = arr.mean(axis=(1, 2))  # shape (time,)
    #         self._space_centered = arr - image_means[:, np.newaxis, np.newaxis]
    
    #         # Mean-centre in time (each pixel's time series mean = 0)
    #         pixel_means = arr.mean(axis=0)  # shape (y, x)
    #         self._time_centered = arr - pixel_means[np.newaxis, :, :]
    
    #     @property
    #     def space(self):
    #         """Return image-wise mean-centred array (mean over y,x = 0 for each time slice)."""
    #         return self._space_centered
    
    #     @property
    #     def time(self):
    #         """Return time-wise mean-centred array (mean over time = 0 for each pixel)."""
    #         return self._time_centered

    import numpy as np
    
    class MeanCenteredArray:
        def __init__(self, masked_array):
            if not isinstance(masked_array, np.ma.MaskedArray):
                raise TypeError("Input must be a masked array")
            
            self._arr = masked_array
    
            # Precompute means
            self._image_means = self._arr.mean(axis=(1, 2))  # shape (time,)
            self._pixel_means = self._arr.mean(axis=0)       # shape (y, x)
    
            # Accessors
            self.mean_centered = _MeanCenteredAccessor(self._arr, self._image_means, self._pixel_means)
            self.means = _MeanAccessor(self._image_means, self._pixel_means)
    
        @property
        def original(self):
            return self._arr
    
    
    class _MeanCenteredAccessor:
        def __init__(self, arr, image_means, pixel_means):
            self._space_centered = arr - image_means[:, np.newaxis, np.newaxis]
            self._time_centered = arr - pixel_means[np.newaxis, :, :]
    
        @property
        def space(self):
            """Mean-centre in space: each 2D slice has zero mean."""
            return self._space_centered
    
        @property
        def time(self):
            """Mean-centre in time: each pixel's time series has zero mean."""
            return self._time_centered
    
    
    class _MeanAccessor:
        def __init__(self, image_means, pixel_means):
            self._image_means = image_means
            self._pixel_means = pixel_means
    
        @property
        def space(self):
            """Return 1D array of means over space (y,x) for each time step."""
            return self._image_means
    
        @property
        def time(self):
            """Return 2D array of mean time series per pixel (shape: y, x)."""
            return self._pixel_means

    
    
    def rescale_masked(arr_ma, scale, *, order=1, anti_aliasing=False,
                       **kw):
        """
        Down-/up-scale a 3-D MaskedArray with skimage.rescale
        while preserving the mask.
    
        Parameters
        ----------
        arr_ma : np.ma.MaskedArray   (t, y, x)
            Data cube; mask may be nomask.
        scale  : float | tuple
            Scale tuple to feed straight into skimage.rescale
            (for your use-case: (1.0, f, f)).
        order, anti_aliasing, **kw
            Forwarded to skimage.transform.rescale.  `preserve_range`
            is forced to True so numeric range/dtype survive.
    
        Returns
        -------
        np.ma.MaskedArray
            Same rank as input, spatially rescaled, with a rebuilt mask.
        """
        
        import numpy.ma as ma
    
        # --- 1. replace masked values with NaN so they don’t bias interpolation
        data_in = arr_ma.filled(np.nan)
    
        data_out = rescale(
            data_in,
            scale=scale,
            order=order,
            anti_aliasing=anti_aliasing,
            preserve_range=True,      # keep original numeric scale
            **kw
        )
    
        # --- 2. new mask: every NaN generated above becomes masked
        mask_out = np.isnan(data_out)
    
        # --- 3. pack back into a MaskedArray and restore original dtype
        return ma.MaskedArray(
            data_out.astype(arr_ma.dtype, copy=False),
            mask=mask_out
        )
        

    # n_pixs_start = displacement_r2["incremental"].shape[1]               
    # shape_start = displacement_r2["mask"].shape                          
    
    # 1: Downsample the ifgs for use in all following functions.  
    # if we're not actually downsampling, skip for speed
    if downsample_run != 1.0:
        
        # downsample the cumulative masked array
        print(
            f"The cumulative time series was size "
            f"{displacement_r3['cum_ma'].shape} ",end=''
            )
        
        # note that the rescale in time (first dimension) is held at 1.  
        displacement_r3['cum_ma'] = rescale_masked(
            displacement_r3['cum_ma'],
            scale=(1.0, downsample_run, downsample_run),
            order=order,
            anti_aliasing=anti_aliasing,
        )
        
        print(
            f"and has been downsampled in space to "
            f"{displacement_r3['cum_ma'].shape} "
            )
        
        # remake lons and lats at new resolution
        # check if we have lon lat data as not alway strictly necessary.  
        if ('lons_mg' in displacement_r3) and ('lats_mg' in displacement_r3):

            lons = np.linspace(
                displacement_r3['lons_mg'][0,0], 
                displacement_r3['lons_mg'][-1,-1], 
                displacement_r3['cum_ma'].shape[2]
            )
            displacement_r3['lons_mg'] = np.repeat(
                lons[np.newaxis, :], 
                displacement_r3['cum_ma'].shape[1],
                axis = 0
            )                       
            
            lats = np.linspace(
                displacement_r3['lats_mg'][0,0], 
                displacement_r3['lats_mg'][-1,-1], 
                displacement_r3['cum_ma'].shape[1]
            )
            displacement_r3['lats_mg'] = np.repeat(
                lats[np.newaxis, :], 
                displacement_r3['cum_ma'].shape[1],
                axis = 0
            )                       
            
            # poor quality check to ensure that lats aren't upside down            
            if displacement_r3['lats_mg'][0,0] < displacement_r3['lats_mg'][-1,0]:                                        
                displacement_r3['lats_mg'] = np.flipud(displacement_r3['lats_mg'])                                        
            
        # also downsample other simple data if it's included:
        for product in ['dem', 'E', 'N', 'U']: 
            if product in displacement_r3.keys():
                displacement_r3[product] = rescale(
                    displacement_r3[product], 
                    downsample_run, 
                    anti_aliasing = anti_aliasing
                )                               
            
    # 2: Downsample further for plotting.  
    # note that the rescale in time (first dimension) is held at 1.  
    displacement_r3['cum_ma_downsampled'] = rescale_masked(
        displacement_r3['cum_ma'],
        scale=(1.0, downsample_run, downsample_run),
        order=order,
        anti_aliasing=anti_aliasing,
    )
    
    
    # 3: also compute the incremental displacements
    displacement_r3['inc_ma'] = ma.diff(
        displacement_r3['cum_ma'],
        axis=0
        )
    
    displacement_r3['inc_ma_downsampled'] = ma.diff(
        displacement_r3['cum_ma_downsampled'],
        axis=0
        )
    
    # deal with mean centering (which could be in space or time)
    displacement_r3['cum_ma'] = MeanCenteredArray(displacement_r3['cum_ma'])
    displacement_r3['inc_ma'] = MeanCenteredArray(displacement_r3['inc_ma'])
   
    return displacement_r3

#%% bss_components_inversion_per_epoch()

def bss_components_inversion_per_epoch(
        sources,
        mask_sources,
        ifgs,
        cumulative=True
        ):
    """
    
    Returns:
        
        tcs_c | n_times x n_ics
        d_hat | n_pixels x n_times |
        d_resid | n_pixles x n_times
    
    """
    import numpy as np
    import numpy.ma as ma
    from licsalert.aux import col_to_ma, update_mask
    
    n_ics=sources.shape[0]
    n_epochs=ifgs.shape[0]
    
    # initialise to store outputs
    tcs_c=np.zeros((n_epochs, n_ics))
    d_hat=ma.zeros(ifgs.shape)
    d_resid=ma.zeros(ifgs.shape)
    
    

    
    # iterate over each epoch
    for epoch_n, epoch in enumerate(ifgs):
        
        
        # find pixels that are valid in ICs and the current epoch
        mask_epoch = epoch.mask
        mask_sources_epoch = np.invert(
                np.logical_and(
                    np.invert(mask_sources),
                    np.invert(mask_epoch)
                    )
                )
        
        # mask the ICs with the new mask
        sources_masked = update_mask(
            {'incremental' : sources,
             'mask' : mask_sources},
            mask_sources_epoch
            )
        
        # make the epoch with the new mask
        epoch_masked = ma.array(epoch, mask=mask_sources_epoch)
        
        # and convert it to a row vector (but still rank 2 as time is dim 1
        # and this is a time series of length 1)
        epoch_masked_r1 = ma.compressed(epoch_masked)[np.newaxis, :]
        
        
        # invert to fit the epoch using the ICs
        # tcs_c are m, d_hat is reconstruction, d_resid is epoch-reconstruction
        tcs_c_cur, d_hat_cur, d_resid_cur = bss_components_inversion(
            sources_masked['incremental'], 
            epoch_masked_r1, 
            cumulative=cumulative,
            mask=mask_sources_epoch
            )
        
        # convert rank1 vectors back to rank 2 matrices (ie. images)
        d_hat_cur_r2 = col_to_ma(d_hat_cur, mask_sources_epoch)
        d_resid_cur_r2 = col_to_ma(d_resid_cur, mask_sources_epoch)
        
        # and record
        tcs_c[epoch_n,] = tcs_c_cur
        d_hat[epoch_n,] = d_hat_cur_r2
        d_resid[epoch_n,] = d_resid_cur_r2
        
        
        
        
        
    return tcs_c, d_hat, d_resid


#%% bss_components_inversion()
def bss_components_inversion(
        sources,
        interferograms,
        cumulative=True,
        mask = None
        ):
    """
    A function to fit an interferogram using components learned by BSS, and return how strongly
    each component is required to reconstruct that interferogramm, and the

    Inputs:
        sources | n_sources x pixels | ie architecture I.  Mean centered
        interferogram | n_ifgs x pixels | Doesn't have to be mean centered, ifgs are rows
        cumulative | Boolean | if true, m and residual (mean_l2_norm) are returned as cumulative sums.

    Outputs:
        m | rank 1 array | the strengths with which to use each source to reconstruct the ifg.
        d_hat  | r2 array | reconstruction of the time series.  
        d_resid | r2 array | residual (d - d_hat)

    2019/12/30 | MEG | Update so handles time series (and not single ifgs), and can return cumulative values
    2023_04_03 | MEG | Also return the reconstruction (d_hat), and the residual
    2023_12_15 | MEG | 
    """
    import numpy as np

    (n_ifgs, n_pixels) = interferograms.shape

    # a column vector (p x n_ifgs)
    d = interferograms.T                                                 
    
    # a matrix of ICA sources and each is a column (p x n_sources)
    g = sources.T                                                       
    
    #invert to fit,  m (n_sources x n_ifgs)
    m = np.linalg.inv(g.T @ g) @ g.T @ d                                
    
    # reconstructed ifgs, as column vectors
    d_hat = g@m
    
    # residual between each ifg and its reconstruction (for all pixels at all times)
    d_resid = d - d_hat
    
    # make these column vectors
    m = m.T                                                             

    # ##%% debug figure: show how the sources are used to reconstruct an ifg.  
    # pdb.set_trace()
    # plt.switch_backend('qt5agg')
    # from licsalert.aux import col_to_ma
    # n_sources = sources.shape[0]
    # if mask is not None:
    #     for ifg_n, ifg in enumerate(interferograms[:3,]):

    #         sources_rescale = sources * np.repeat(m[ifg_n:ifg_n + 1,:].T, axis = 1, repeats = interferograms.shape[1])
    #         all_data = np.concatenate((ifg[np.newaxis, :], sources_rescale, d_hat[:, ifg_n : ifg_n+1].T, d_resid[:, ifg_n : ifg_n+1].T), axis = 0)
    #         vmin = np.min(all_data)
    #         vmax = np.max(all_data)

    #         f, axes = plt.subplots(4, n_sources + 1)
    #         axes[0,0].matshow(col_to_ma(ifg, mask), vmin = vmin, vmax = vmax)
    #         axes[0,0].set_title(f'ifg # {ifg_n}')

    #         for source_n in range(n_sources):
    #             axes[1, source_n + 1].matshow(col_to_ma(sources[source_n, :], mask))
    #             axes[1, source_n + 1].set_title(f"IC # (unscaled)")            
            
            
    #         axes[2,0].matshow(col_to_ma(d_hat[:, ifg_n : ifg_n+1], mask), vmin = vmin, vmax = vmax)
    #         axes[2,0].set_title(f'reconstruction')
    #         for source_n in range(n_sources):
    #             axes[2, source_n + 1].matshow(col_to_ma(sources_rescale[source_n, :], mask), vmin = vmin, vmax = vmax)
    #             axes[2, source_n + 1].set_title(f"{m[ifg_n:ifg_n + 1,source_n]}")
                
    #         axes[3,0].matshow(col_to_ma(d_resid[:, ifg_n : ifg_n+1], mask), vmin = vmin, vmax = vmax)
    #         axes[3,0].set_title(f'residual')
            
    #         for ax in np.ravel(axes):
    #             ax.set_xticks([])
    #             ax.set_yticks([])
    # ##%%

    # sum through time, and make 0 at the start
    if cumulative:
        m_cs = np.cumsum(m, axis=0)
        m_cs = np.concatenate((np.zeros((1, m.shape[1])), m_cs), axis = 0)          
        return m_cs, d_hat, d_resid                       
    else:
        return m, d_hat, d_resid
        




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




#%%



def write_volcano_status(
        processing_date, sources_tcs_baseline, residual_tcs_baseline, 
        ics_labels, txt_out_dir
        ):
    """ For each time step, write a 2 line text file that contains the number 
    of sigmas from the line of best fit that the deformation source is 
    (change in existing deformation) and that the residual is (new 
                                                               deformation).  
    Inputs:
        processing_date | LiCSAlert date | 
        sources_tcs_baseline | list of dicts | from the licsalert function
        residual_tcs_baseline | list of dicts | from the licsalert function
        ics_labels | dict | something like this: {'source_names':
                                                  ['deformation', 
                                                   'topo_cor_APS', 
                                                   'turbulent_APS'], 
                                                  'labels': 
                                                      array([[1., 0., 0.], 
                                                             [0., 1., 0.], 
                                                             [0., 0., 1.]])}
        txt_out_dir | Path | pathlib Path to write file to.  
        
    Returns:
        text file.  First row is change in def, second row is new def.  
        
    History:
        2023_04_04 | MEG | Written
        2025_02_21 | MEG | Fix bug that baseline used final value from 
                            monitoring
    """
    import numpy as np

    # extract which acquisition number we are working with.      
    acq_n = processing_date.acq_n
    
    # calculate change in def status
    # get which column of the one hot encoding related to deformation
    def_col_n = ics_labels['source_names'].index('deformation')                                     
    
    # look in that column to find which source number (row) has 1 (ie 
    # is positive)
    def_ic_n = np.argwhere(ics_labels['labels'][:, def_col_n] == 1)[0][0]                           
    
    # extract the n_sigmas for existing deformation (as we know def_ic_n)
    def_n_sigmas = sources_tcs_baseline[def_ic_n]['distances'][acq_n] 
    
    # extract the n_sigmas for new deformation (which comes from residual)
    new_n_sigmas = residual_tcs_baseline[0]['distances'][acq_n]                                      
    
    # write output, existing deformation first, new second.  
    with open(txt_out_dir / 'volcano_status.txt', 'w') as f:
        f.write(f"{def_n_sigmas[0]}\n")                                             
        f.write(f"{new_n_sigmas[0]}\n")                                             


#%%


def reconstruct_ts_from_dir(ics_one_hot, licsalert_out_dir):
    """
    Simple wrapper to open the data required to reconstruct the cumulative time series.  
    Inputs:
        ics_one_hot | list | one hot encoding for which components to use in the reconstruction.  Must be tthe same length as the number of sources. 
                            If set to 1, that IC is used in the reconstruction.  
                            If set to 0, that IC is not used in the reconstruction.  
                            E.g. [1,0,0] to use onlyl IC0 and to ignore IC1 and IC2
        licsalert_out_dir | pathlib Path | path to directory containing the results of LiCSAlert run.  
    Returns:
        X_inc_r3 | rank 3 masked array | (n_acqs - 1) x ny x nx.  Reconstruction 
                                       of short temporal baseline / daisy chain /
                                       incremental ifgs.  Mean centering removed 
                                       (i.e. back to original reference pixel)
        X_cum_r3 | rank 3 masked array | (n_acqs) x ny x nx.  Reconstruction 
                                       of cumulative ifgs . Mean centering removed 
                                       (i.e. back to original reference pixel)

    History:
        2023_11_22 | MEG | Written.  
        2024_01_06 | MEG | Return both incremental and cumulative.  
    """
    
    from licsalert.licsalert import reconstruct_ts
    from licsalert.data_importing import open_aux_data, open_tcs
    from licsalert.data_importing import determine_last_licsalert_date
    from licsalert.aux import r2_to_r3
    
    displacement_r2, tbaseline_info, aux_data = open_aux_data(licsalert_out_dir)
    final_date_dir = determine_last_licsalert_date(licsalert_out_dir)
    sources_tcs, residual_tcs = open_tcs(final_date_dir)    

    if len(sources_tcs) != len(ics_one_hot):
        raise Exception(f"There are {len(sources_tcs)} time courses, but ics_one_hot is {len(ics_one_hot)} items long.  These must match to continue.  Exiting.")

    X_inc_r2, X_cum_r2 = reconstruct_ts(ics_one_hot, sources_tcs, aux_data, displacement_r2)                             # this returns a rank 2, n_times x n_pixels, of the cumulative data.  
    X_inc_r3 = r2_to_r3(X_inc_r2, displacement_r2['mask'])                                                          # conver to a rank 3 masked array (i.e. n_times x ny x nx)
    X_cum_r3 = r2_to_r3(X_cum_r2, displacement_r2['mask'])                                                          # conver to a rank 3 masked array (i.e. n_times x ny x nx)
    
    return X_inc_r3, X_cum_r3

#%%

def reconstruct_ts(ics_one_hot, sources_tcs, aux_data, displacement_r2):
    """ Reconstruct a LiCSAlert cumulative time series form the ICs and 
    cumulative time course using a choice of components.
    Automatically detects if sICA or tICA was run and handles mean centering accordingly.  
    
    Inputs:
        ics_one_hot | list | 1 if IC to be used, 0 if not.  Must be same length
                            as number of ICS  e.g. [1,0,0,0]
        sources_tcs | list of dicts | one item in list for each IC.  Contains 
                                        cumulative time course and associated data.  
        aux_data | dict | dict_keys(['icasar_sources', 'dem', 'mask'])
        displacement_r2 | dict | dict_keys(['dem', 'mask', 'lons', 'lats', 'E',
                                            'N', 'U', 'incremental', 'means', 
                                            'incremental_downsampled', 
                                            'mask_downsampled'])
    Returns:
        X_inc | r2 array | incremental (short temporal baseline) ifgs as rows,
        mean centering has been removed.  
        X_icum | r2 array | cumulative ifgs as rows, m
    History:
        2023_10_26 | MEG | Written
        2024_01_06 | MEG | Return both cumualtive and incremental displacements.  
        
    """
    import numpy as np
    
    if len(ics_one_hot) != len(sources_tcs):
        raise Exception(f"'sources_tcs' is of length {len(sources_tcs)} so contains {len(sources_tcs)} sources.  "
                        f"However, 'ics_one_hot' is {len(ics_one_hot)}, which doesn't agree.  Exiting.  ")
        
    n_sources = len(sources_tcs)
    n_acqs = sources_tcs[0]['cumulative_tc'].shape[0]
    n_pixels = displacement_r2['incremental'].shape[1]
    
    # make the cumulative time courses
    A_cum = np.zeros((n_acqs, n_sources))                                                                                                      
    for n_source in range(n_sources):
        # multiply by 1 or 0 to include or not, and then put as column in matrix.  
        A_cum[:, n_source] = np.ravel(ics_one_hot[n_source] *  sources_tcs[n_source]['cumulative_tc'])                                          
        
    # the daisy chain ifgs are just the different between each succesive cumulative ifgs.  
    A_inc = np.diff(A_cum, axis = 0)    
    S = aux_data['icasar_sources']
    
    # mean for each incremental ifg (n_acq - 1) means sICA was run
    if displacement_r2['means'].shape[0] == n_acqs-1:                                   
        means_r2 = np.repeat(displacement_r2['means'][:,np.newaxis], n_pixels, axis = 1 )
    
    # mean for each pixel means tICA was run
    elif displacement_r2['means'].shape[0] == n_pixels:                                   
        means_r2 = np.repeat(displacement_r2['means'][np.newaxis, :], n_acqs-1, axis = 0)
    else:
        raise Exception(f"The number of acquisitions in the time courses ({n_acqs}) "
                        f"does not match the number of times in the interferogams "
                        f"({displacement_r2['means'].shape[0]}).  Exiting.  ")
    
    # reconstruct the daisy chain of short temporal baseline ifgs.  
    X_inc = A_inc@S + means_r2
    
    # reconstruct the cumulative ifgs, ensure that all zeros on first acquisition
    X_cum = np.concatenate((np.zeros((1, n_pixels)), np.cumsum(X_inc, axis = 0)))
    
    return X_inc, X_cum


#%% load_or_create_ICASAR_results()




def load_or_create_ICASAR_results(
        run_ICASAR,
        displacement_r2,
        tbaseline_info, 
        baseline_end,
        out_dir,
        ICASAR_settings
        ):
    """
    ICASAR results are always required by LiCSAlert, and these need either to be computed (usually only once at the start),
    or loaded (more common).  
    
    Inputs:
        run_ICASAR | boolean | whether ICASAR should be run, or the results from a previous run loaded.  
        displacment_r2 | dict | interferograms and associated data that are used by the ICASAR algorithm  
                                mixtures_mc (mean centered) are extracted from this.  
        tbaseline_info | dict | various temporal information, such as ifg_dates
        # baseline_end | string | YYYYMMDD that the baseline ends on.  
    Returns:
        icasar_sources | rank 2 array | sources as row vectors
        mask_sources
        label_sources_output | ? | Dict of one hot encoding of labels (defo / topo correlated atmosphere / turbulent atmsophere)
    History:
        2021_10_15 | MEG | Written.  
        2023_04_03 | MEG | Add return of ICASAR label_sources_output
        2023_12_20 | MEG | Remove functions to determine acq_n of baseline end date.  
    """

    from licsalert.icasar.icasar_funcs import ICASAR
    import pickle
    import numpy as np
    
    def check_means(sources, tcs):
        """ Print some information about the means of the data returned by the ICASAR function.  
        Inputs:
            sources | r2 array | Images as rows.  n_images x n_pixels.  
            tcs | r2 array | Incremental time courses as columns.  n_times x n_images
        History:
            2023_11_08 | Written.  
        """
        import numpy as np 
        source_means = np.mean(sources, axis = 1)
        print(f"Means for each row of sources (these should be images): {source_means}")
        tcs_means = np.mean(tcs, axis = 0)
        print(f"Means for each column of the tcs (these should be time courses): {tcs_means}")
      
    
    if run_ICASAR:
        print(f"\nRunning ICASAR.")                                      
        
        # index the data, note +1 so that the chosen ifg_n is included.  
        # mean centered data has already been mean centered in time or space, as required.  
        spatial_ICASAR_data = {
            'ifgs_dc'       : displacement_r2['mixtures_mc'][:(baseline_end.acq_n+1),],                             
            'mask'          : displacement_r2['mask'],
            'lons'          : displacement_r2['lons_mg'],
            'lats'          : displacement_r2['lats_mg'],
            'ifg_dates_dc'  : tbaseline_info['ifg_dates'][:(baseline_end.acq_n+1)]
            }
        if 'dem' in displacement_r2.keys():
            spatial_ICASAR_data['dem'] = displacement_r2['dem']
        
        if ICASAR_settings['figures'] == 'both':
            ICASAR_settings['figures'] = 'png+window'                                                                                  # update licsalert name to ICASAR name.  
            
        # Run ICASAR (slow). tcs are incremental (i.e. not cumulative)
        outputs = ICASAR(spatial_data = spatial_ICASAR_data,
                        out_folder = out_dir, **ICASAR_settings,
                        ica_verbose = 'short', label_sources = True)                     
        
        (sources, tcs, residual, Iq, n_clusters, S_all_info, r2_ifg_means, 
        ics_labels ) = outputs; del outputs
        
        mask_sources = displacement_r2['mask']                                                                                                              # rename a copy of the mask

        if ICASAR_settings['sica_tica'] == 'tica':                                                                                                         # possibly deal with mean centering which is problematic with tica 
            # check_means(sources, tcs)                                                                                                                     # as the sources are cumulative time courses, and only these are guaranteed  
            sources -= np.repeat(np.mean(sources, axis = 1)[:,np.newaxis], axis = 1, repeats = sources.shape[1])                                            # to be mean centered.  images (sources) and incremental time courses (tcs) won't be.  
            tcs -= np.repeat(np.mean(tcs, axis = 0)[np.newaxis, : ], axis = 0, repeats = tcs.shape[0])
            check_means(sources, tcs)
            
    
    # load from a previous run    
    else:
        with open(out_dir / "ICASAR_results.pkl", 'rb') as f_icasar:                                                                     # Or open the products from a previous ICASAR run.  
            sources = pickle.load(f_icasar)   
            mask_sources = pickle.load(f_icasar)
            tcs  = pickle.load(f_icasar)    
            source_residuals = pickle.load(f_icasar)    
            #Iq_sorted = pickle.load(f_icasar)    
            n_clusters = pickle.load(f_icasar)    
            xy_tsne = pickle.load(f_icasar)                             # not used
            labels_hdbscan = pickle.load(f_icasar)                      # not used
            ics_labels = pickle.load(f_icasar)
        f_icasar.close()                                                                                                                           
    
    return sources, mask_sources,  ics_labels


#%%


class licsalert_date_obj:
    def __init__(self, date, acq_dates):
        """ Express date as a string, as a datetime, as the days since
        the time series started, and in acqusition number since the time
        series started.  
        """
        from datetime import datetime as dt
        self.date = date
        self.dt = dt.strptime(self.date, "%Y%m%d")                                            
        self.day_n = (self.dt - dt.strptime(acq_dates[0], "%Y%m%d")).days
        try:
            self.acq_n = acq_dates.index(self.date)
        except:
            pass


#%%

