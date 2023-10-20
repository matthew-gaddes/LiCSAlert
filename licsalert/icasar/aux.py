#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 18:19:37 2019

@author: matthew
"""



#%%


def maps_tcs_rescale(maps, tcs):
    """
    A function to rescale spaital maps to have unit range and rescale each's time cource (tc)
    so that there is no change to the product of the two matrices
    
    input:
        maps | array | spatial maps as rows (e.g. 2x1775)
        tcs  | array | time courses as columns (e.g. 15x2)
    Output:
        maps_scaled | array | spatial maps as rows with each row having unit range
        tcs_scaled | array | TCs scaled so that new maps x new tcs equals maps x tcs
    
    2017/05/15 | written
    
    """
    import numpy as np
    
    def rescale_unit_range(signals):
        """
        rescale a matrix of row vectors so that each row has a range of 1.  
        Also record the scaling factor required to rescale each row vector
        signals are rows
        Input:
            signals: array with each signal as a new row 
        Output:
            signals_rescale | array | each row has a range of 1
            signals_factor | array | factor that each row is dvided by to adjust range
        """
        import numpy as np
        
        signals_rescale = np.ones(signals.shape)                                    # initiate rescaled array
        signals_factor = np.ones((np.size(signals, axis=0) , 1))                    # initiate array to record scaling factor for each row
           
        for i in np.arange(np.size(signals, axis = 0)):
            signals_factor[i,0] = (np.max(signals[i,:])-np.min(signals[i,:]))
            signals_rescale[i,:] = signals[i,:]/signals_factor[i,0]
            
        return signals_rescale, signals_factor

    
    maps_scaled , scaling = rescale_unit_range(maps)
    tcs_scaled = tcs * np.ravel(scaling)
    return maps_scaled, tcs_scaled




#%%


def bss_components_inversion(sources, interferograms):
    """
    A function to fit an interferogram using components learned by BSS, and return how strongly
    each component is required to reconstruct that interferogramm, and the 
    
    Inputs:
        sources | n_sources x pixels | ie architecture I.  Mean centered
        interferogram | list of (n_ifgs x pixels) | Doesn't have to be mean centered, multiple interferograms to be fit can be fit by making the list as long as required.  
        
    Outputs:
        m | rank 1 array | the strengths with which to use each source to reconstruct the ifg.  
        mean_l2norm | float | the misfit between the ifg and the ifg reconstructed from sources
    """
    import numpy as np
    
    inversion_results = []
    for interferogram in interferograms:
        interferogram -= np.mean(interferogram)                     # mean centre
        n_pixels = np.size(interferogram)
    
    
        d = interferogram.T                                        # now n_pixels x n_ifgs
        g = sources.T                                              # a matrix of ICA sources and each is a column (n_pixels x n_sources)
        
        ### Begin different types of inversions.  
        m = np.linalg.inv(g.T @ g) @ g.T @ d                       # m (n_sources x 1), least squares
        #m = g.T @ np.linalg.inv(g @ g.T) @ d                      # m (n_sources x 1), least squares with minimum norm condition.     COULDN'T GET TO WORK.  
        #m = np.linalg.pinv(g) @ d                                   # Moore-Penrose inverse of G for a simple inversion.  
        # u = 1e0                                                                 # bigger value favours a smoother m, which in turn can lead to a worse fit of the data.  1e3 gives smooth but bad fit, 1e1 is a compromise, 1e0 is rough but good fit.  
        # m = np.linalg.inv(g.T @ g + u*np.eye(g.shape[1])) @ g.T @ d;                # Tikhonov solution  
        
        
        ### end different types of inversion
        d_hat = g@m
        d_resid = d - d_hat
        mean_l2norm = np.sqrt(np.sum(d_resid**2))/n_pixels          # misfit between ifg and ifg reconstructed from sources
        
        inversion_results.append({'tcs'      : m,
                                  'model'    : d_hat,
                                  'residual' : d_resid,
                                  'l2_norm'  : mean_l2norm})
    
    return inversion_results


#%%

def signals_to_master_signal_comparison(signals, master_signal, density = False):
        """ Given an array of signals (as row vectors), compare it to a single signal and compute a kernel
        density estimate, and calculate a line of best fit through the points with R**2 value.  
        Inputs:
            signals | rank 2 | signals as rows.  Even if there's only 1 signal, it still needs to be rank 2
            master_signal | rank 2 | signal as a row, but has to be rank 2
            density | boolean | If True, gaussian kernel density estimate for the points.  Can be slow.  
            
        Returns:
            signal_to_msignal_comparison | dict | keys:
                                                    xyzs | list of rank 2 arrays | entry in the list for each signal, xyz are rows.  
                                                    line_xys | list of rank 2 arrays | entry in the list for each signal, xy are points to plot for the lines of best fit
                                                    cor_coefs | list | correlation coefficients between each signal and the master signal.  
            
        History:
            2021_04_22 | MEG | Written.  
            2021_04_26 | MEG | Add check that the signals are of the same length.  
        """
        import numpy as np
        from scipy.stats import gaussian_kde
        import numpy.polynomial.polynomial as poly                                             # used for lines of best fit through dem/source plots
        
        n_signals, n_pixels = signals.shape                                                    # each signal is a row, observations of that are columns.  
        if n_pixels != master_signal.shape[1]:
            raise Exception(f"The signals aren't of the same length (2nd dimension), as 'signals' is {n_pixels} long, but 'master_signal' is {master_signal.shape[1]} long.  Exiting.  ")
        xyzs = []                                                                              # initiate
        line_xys = []
        cor_coefs = []
        print(f"    Starting to calculate the 2D kernel density estimates for the signals.  Completed ", end = '')
        for signal_n, signal in enumerate(signals):                                            # signal is a row of signals, and loop through them.  
            
            # 1: Do the kernel density estimate
            xy = np.vstack((master_signal, signal[np.newaxis,:]))                              # master signal will be on X and be the top row.  
            x = xy[:1,:]                                    
            y = xy[1:2,:]
            if density:
                z = gaussian_kde(xy)(xy)
                idx = z.argsort()                                                               # need to be sorted so that when plotted, those with the highest z value go on top.                                  
                x, y, z = x[0,idx], y[0,idx], z[idx]  
                xyzs.append(np.vstack((x,y,z)))                                                 # 3 rows, for each of x,y and z
            else:
                xyzs.append(np.vstack((x,y,np.zeros(n_pixels))))                                # if we're not doing the kernel density estimate, z value is just zeros.  
                    
            # 2: calculate the lines of best fit
            line_coefs = poly.polyfit(x, y, 1)                                                  # polynomial of order 1 (i.e. a line of best fit)
            line_yvals = (poly.polyval(x, line_coefs))                                          # calculate the lines yvalues for all x points
            line_xys.append(np.vstack((x, line_yvals)))                                         # x vals are first row, y vals are 2nd row
            
            # 3: And the correlation coefficient
            # import pdb; pdb.set_trace()
            cor_coefs.append(np.corrcoef(x, y)[1,0])                                            # which is a 2x2 matrix, but we just want the off diagonal (as thats the correlation coefficient between the signals)
            
            print(f"{signal_n} ", end = '')
        print('\n')
        
        signal_to_msignal_comparison = {'xyzs' : xyzs,
                                        'line_xys' : line_xys,
                                        'cor_coefs' : cor_coefs}
        
        return signal_to_msignal_comparison



#%%

def create_all_ifgs(ifgs_r2, ifg_dates, max_n_all_ifgs = 1000):
    """Given a rank 2 of incremental ifgs, calculate all the possible ifgs that still step forward in time (i.e. if deformation is positive in all incremental ifgs, 
    it remains positive in all the returned ifgs.)  If acquisition dates are provided, the tmeporal baselines of all the possible ifgs can also be found.  
    Inputs:
        ifgs_r2 | rank 2 array | Interferograms as row vectors.  
        ifg_dates | list of strings | dates in the form YYYYMMDD_YYYYMMDD.  As the interferograms are incremental, this should be the same length as the number of ifgs
    Returns:
        ifgs_r2 | rank 2 array | Only the ones that are non-zero (the diagonal in ifgs_r3) and in the lower left corner (so deformation isn't reversed.  )
    History:
        2021_04_13 | MEG | Written
        2021_04_19 | MEG | add funcionality to calculate the temporal baselines of all possible ifgs.  
        2021_04_29 | MEG | Add functionality to handle networks with breaks in them.  
    """
    import numpy as np
    import datetime as dt
    from datetime import datetime, timedelta
    import random
    from licsalert.icasar.aux import acquisitions_from_ifg_dates
    
    def split_network_to_lists(ifgs_r2, ifg_dates):
        """Given a r2 of ifgs (ie row vectors) and a list of their dates (YYYYMMDD_YYYYMMDD), break these into lists that are 
        continuous (i.e. no gaps)
        """
    
        ifg_dates_continuous = []                                                                   # list of the list the dates for a continuous network
        ifgs_r2_continuous = []                                                                     # and the incremental interferograms in that network.  
        start_continuous_run = 0
        for ifg_n in range(n_ifgs-1):                                                               # iterate                                 
            if (ifg_dates[ifg_n][9:] != ifg_dates[ifg_n+1][:8]):                                   # if the dates don't agree (ie. the slave of one isn't the master of the other), we've got to the end of a network
                ifg_dates_continuous.append(ifg_dates[start_continuous_run:ifg_n+1])                # +1 as want to include the last date in the selection
                ifgs_r2_continuous.append(ifgs_r2[start_continuous_run:ifg_n+1,])
                start_continuous_run = ifg_n+1                                                      # +1 so that when we index the next time, it doesn't include ifg_n
            if ifg_n == n_ifgs -2:                                                                  #  if we've got to the end of the list.              
                ifg_dates_continuous.append(ifg_dates[start_continuous_run:])                       # select to the end.  
                ifgs_r2_continuous.append(ifgs_r2[start_continuous_run:,])
                
        return ifg_dates_continuous, ifgs_r2_continuous
    
    def create_all_possible_ifgs(networks):
        """
        """
        n_networks = len(networks)
        ifgs_all_r2 = []
        dates_all_r1 = []
        for n_network, network in enumerate(networks):                                                         # loop through each network.  
            ifgs_r2_temp = network['ifgs_r2']
            ifg_dates_temp = network['ifg_dates']
            n_acq = ifgs_r2_temp.shape[0] + 1
            
            # 2a: convert from daisy chain of incremental to a relative to a single master at the start of the time series.  
            acq1_def = np.zeros((1, n_pixs))                                                     # deformation is 0 at the first acquisition
            ifgs_cs = np.cumsum(ifgs_r2_temp, axis = 0)                                          # convert from incremental to cumulative.  
            ifgs_cs = np.vstack((acq1_def, ifgs_cs))                                             # add the 0 at first time ifg to the other cumulative ones. 
            
            # 2b: create all possible ifgs
            ifgs_cube = np.zeros((n_acq, n_acq, n_pixs))                                    # cube to store all possible ifgs in
            for i in range(n_acq):                                                          # used to loop through each column
                ifgs_cube[:,i,] = ifgs_cs - ifgs_cs[i,]                                     # make one column (ie all the rows) by taking all the ifgs and subtracting one time from it
               
            # 2c: Get only the positive ones (ie the lower left quadrant)    
            lower_left_indexes = triange_lower_left_indexes(n_acq)                              # get the indexes of the ifgs in the lower left corner (ie. non 0, and with unreveresed deformation.  )
            ifgs_all_r2.append(ifgs_cube[lower_left_indexes[:,0], lower_left_indexes[:,1], :])        # get those ifgs and store as row vectors.  
            
            # 2d: Calculate the dates that the new ifgs run between.  
            acq_dates = acquisitions_from_ifg_dates(ifg_dates_temp)                                                         # get the acquisitions from the ifg dates.  
            ifg_dates_all_r2 = np.empty([n_acq, n_acq], dtype='U17')                                                        # initate an array that can hold unicode strings.  
            for row_n, date1 in enumerate(acq_dates):                                                                       # loop through rows
                for col_n, date2 in enumerate(acq_dates):                                                                   # loop through columns
                    ifg_dates_all_r2[row_n, col_n] = f"{date2}_{date1}"
            ifg_dates_all_r1 = list(ifg_dates_all_r2[lower_left_indexes[:,0], lower_left_indexes[:,1]])             # just get the lower left corner (like for the ifgs)
    
            dates_all_r1.append(ifg_dates_all_r1)
    
        # 3: convert lists back to a single matrix of all interferograms.  
        ifgs_all_r2 = np.vstack(ifgs_all_r2)                                                                            # now one big array of n_ifgs x n_pixels
        dates_all_r1 = [item for sublist in dates_all_r1 for item in sublist]                                           # dates_all_r1 is a list (one for each connected network) of lists (each ifg date).  The turns them to a singe list.  
        
        return ifgs_all_r2, dates_all_r1
    
    
    def triange_lower_left_indexes(side_length):
        """ For a square matrix of size side_length, get the index of all the values that are in the lower
        left quadrant (i.e. all to the lower left of the diagonals).  
        Inputs:
            side_length | int | side length of the square.  e.g. 5 for a 5x5
        Returns:
            lower_left_indexes | rank 2 array | indexes of all elements below the diagonal.  
        History:
            2021_04_13 | MEG | Written.  
        """
        import numpy as np
        zeros_array = np.ones((side_length, side_length))                                                               # initate as ones so none will be selected.  
        zeros_array = np.triu(zeros_array)                                                                              # set the lower left to 0s
        lower_left_indexes = np.argwhere(zeros_array == 0)                                                              # select only the lower lefts
        return lower_left_indexes

    
    ########### begin main function
    n_ifgs, n_pixs = ifgs_r2.shape
        
    
    
    # 1: Determine if the network is continuous, and if not split it into lists
    ifg_dates_continuous, ifgs_r2_continuous = split_network_to_lists(ifgs_r2, ifg_dates)                               # If there is a break (gap) in the network, this splits the ifgs_dates and ifgs based on this.  
    n_networks = len(ifg_dates_continuous)                                                      # get the number of connected networks.  
    
    # 2: build a better set of information for each network.  
    n_ifgs_all_total = 0                                                                                            # the total number of ifgs that can be made in the network
    networks = []                                                                           # each item is a dict, containing ifg_dates, length_dats, n_ifgs, and n_ifgs_all (the number of unique interferograms that can be made)
    for n_network in range(n_networks):  
        network_dict = {'ifgs_r2'       : ifgs_r2_continuous[n_network],
                        'ifg_dates'     : ifg_dates_continuous[n_network]}
        
        network_start_date = dt.datetime.strptime(ifg_dates_continuous[n_network][0][:8], '%Y%m%d')
        network_end_date = dt.datetime.strptime(ifg_dates_continuous[n_network][-1][9:], '%Y%m%d')
        network_dict['length_days'] = (network_end_date - network_start_date).days
        network_dict['n_ifgs'] = network_dict['ifgs_r2'].shape[0]
        network_dict['n_ifgs_all'] = int(((network_dict['n_ifgs']+1)**2 - (network_dict['n_ifgs']+1)) / 2)                  # +1 as there is 1 more acquisition than incremental interferogram.  
        networks.append(network_dict)
        n_ifgs_all_total += network_dict['n_ifgs_all']                                                          # add to the running total of the total number of ifgs
           
    #3: determine if we can make all ifgs, or if we need to make only some of them.  
    if n_ifgs_all_total < max_n_all_ifgs:                                                               # if we can just make all ifgs, 
        ifgs_all_r2, dates_all_r1 = create_all_possible_ifgs(networks)
        
    else:    
        ifgs_all_r2 = []                                                                                                # list with entry for each entwork
        dates_all_r1 = []                                                                                               # list for all networks
        for network in networks:
            n_ifgs_from_network = int(max_n_all_ifgs * (network['n_ifgs_all'] / n_ifgs_all_total))                    # get the fraction of ifgs that each network will provide ot the total.  e.g. network 1 = 800, network 2 = 800, but 1000 in total, 500 from each network

            # 0: get some info about the ifgs
            n_pixs = network['ifgs_r2'].shape[1]
            n_acq = network['n_ifgs'] + 1

            # 1: Make the cumulative ifgs
            acq1_def = np.zeros((1, n_pixs))                                                     # deformation is 0 at the first acquisition
            ifgs_cs = np.cumsum(network['ifgs_r2'], axis = 0)                                          # convert from incremental to cumulative.  
            ifgs_cs = np.vstack((acq1_def, ifgs_cs))                                             # add the 0 at first time ifg to the other cumulative ones.     
                   
            # 2: get the acquisition dates:
            acq_dates = acquisitions_from_ifg_dates(network['ifg_dates'])                                                         # get the acquisitions from the ifg dates.  

            # 3: create the list of acquisiton dates that is the correct (for this network) length
            ifg_acq_start_acq_ends = []                                                                                  # list of tuples of acquisitions each ifg will be between
            while len(ifg_acq_start_acq_ends) < n_ifgs_from_network:                                                     # until we have enough...
                acq_start = np.random.randint(0, n_acq)                                                                 # random start acq
                acq_end = np.random.randint(0, n_acq)                                                                   # random end acq
                if (acq_start < acq_end) and ((acq_start, acq_end) not in ifg_acq_start_acq_ends):                       # if start is before end (ie only in lower left of all acquisitions square) and not already in list of pairs
                    ifg_acq_start_acq_ends.append((acq_start, acq_end))                                                  # add to pair
            
            # 4: make the ifgs for those acquisition dates.  
            for ifg_acq_start_acq_end in ifg_acq_start_acq_ends:                                                        # iterate through all the pairs we need to get 
                ifgs_all_r2.append(ifgs_cs[ifg_acq_start_acq_end[1],] - ifgs_cs[ifg_acq_start_acq_end[0],])               # subtract one ifg from the other (end - start)
                dates_all_r1.append(f"{acq_dates[ifg_acq_start_acq_end[0]]}_{acq_dates[ifg_acq_start_acq_end[1]]}")         # get teh dates that that ifg spans
                
        # 5: convert lists back to a single matrix of all interferograms.  
        ifgs_all_r2 = np.vstack(ifgs_all_r2)                                                                            # now one big array of n_ifgs x n_pixels
                    
    print(f"When creating all interferograms, {ifgs_r2.shape[0]} were passed to the function, and these were found to make {n_network+1} connected networks.  "
          f"From these, {ifgs_all_r2.shape[0]} interferograms were created.  ")
    
    return ifgs_all_r2, dates_all_r1    
    


#%%

def create_cumulative_ifgs(ifgs_dc, ifg_dates_dc):
    """ Given a time series of incremental (daisy chain) interferograms, calculate the cumulative interferograms  relative to the first acquisition.  
    Inputs:
        ifgs_dc | rank 2 array | ifgs as row vectors.  
        ifg_dates_dc | list | list of YYYYMMDD_YYYYMMDD 
    Returns:
        ifgs_cum | rank 2 array | ifgs as row vectors, cumulative and relative to the first acquisition.  
        ifg_dates_cum | list | list of YYYYMMDD_YYYYMMDD 
    History:
        2021_11_29 | MEG | Written.  
        2021_11_30 | MEG | Update so that doesn't create the acquisition 0 - acquisition 0 interferogram of 0 displacement.  

    """
    import numpy as np
    from licsalert.icasar.aux import baseline_from_names
    
    # 0: First make the ifgs, v1 that uses that has acquisition 0 to acquisition 0 as a row of 0 displacements at the start
    # ifgs_cum_0 = np.zeros((1, ifgs_dc.shape[1]))                            # 0 displacement on first acquistion
    # ifgs_cum = np.cumsum(ifgs_dc, axis = 0)                                 # displacement to last date of each daisy chain interferogram.  
    # ifgs_cum = np.vstack((ifgs_cum_0, ifgs_cum))                            # combine so first acuqisiton has 0 displacmenent.  
    # 0b: or ignores a0 to a0:
    ifgs_cum = np.cumsum(ifgs_dc, axis = 0)                                 # displacement to last date of each daisy chain interferogram.  

    
    # 1: then make the ifg dates.  
    acq_0 = ifg_dates_dc[0][:8]
    #ifg_dates_cum = [f"{acq_0}_{acq_0}"]
    ifg_dates_cum = []
    for ifg_date_dc in ifg_dates_dc:
        ifg_dates_cum.append(f"{acq_0}_{ifg_date_dc[9:]}")
    
    return ifgs_cum, ifg_dates_cum
        


#%%


def add_square_plot(x_start, x_stop, y_start, y_stop, ax, colour = 'k'):
    """Draw localization square around an area of interest, x_start etc are in pixels, so (0,0) is top left.  
    Inputs:
        x_start | int | start of box
        x_stop | int | etc. 
        y_start | int |
        y_ stop | int |
        ax | axes object | axes on which to draw
        colour | string | colour of bounding box.  Useful to change when plotting labels, and predictions from a model.  
    
    Returns:
        box on figure
        
    History:
        2019/??/?? | MEG | Written
        2020/04/20 | MEG | Document, copy to from small_plot_functions to LiCSAlert_aux_functions
    """
        
    ax.plot((x_start, x_start), (y_start, y_stop), c= colour)           # left hand side
    ax.plot((x_start, x_stop), (y_stop, y_stop), c= colour)             # bottom
    ax.plot((x_stop, x_stop), (y_stop, y_start), c= colour)             # righ hand side
    ax.plot((x_stop, x_start), (y_start, y_start), c= colour)             # top
    

#%% From InSAR tools

def acquisitions_from_ifg_dates(ifg_dates):
    """Given a list of ifg dates in the form YYYYMMDD_YYYYMMDD, get the unique YYYYMMDDs that acquisitions were made on.  
    Inputs:
        ifg_dates | list | of strings in form YYYYMMDD_YYYYMMDD.  Called imdates in LiCSBAS nomenclature.  
    Returns:
        acq_dates | list | of strings in form YYYYMMDD
    History:
        2021_04_12 | MEG | Written        
    """
    acq_dates = []
    for ifg_date in ifg_dates:                                          # loop through the dates for each ifg
        dates = ifg_date.split('_')                                     # split into two YYYYMMDDs
        for date in dates:                                              # loop through each of these
            if date not in acq_dates:                                   # if it't not already in the list...
                acq_dates.append(date)                                  # add to it
    return acq_dates


def baseline_from_names(names_list):
    """Given a list of ifg names in the form YYYYMMDD_YYYYMMDD, find the temporal baselines in days_elapsed (e.g. 12, 6, 12, 24, 6 etc.  )
    Inputs:
        names_list | list | in form YYYYMMDD_YYYYMMDD
    Returns:
        baselines | list of ints | baselines in days
    History:
        2020/02/16 | MEG | Documented
    """

    from datetime import datetime, timedelta

    baselines = []
    for file in names_list:

        master = datetime.strptime(file.split('_')[-2], '%Y%m%d')
        slave = datetime.strptime(file.split('_')[-1][:8], '%Y%m%d')
        baselines.append(-1 *(master - slave).days)
    return baselines



#%%


def r2_to_r3(ifgs_r2, mask):
    """ Given a rank2 of ifgs as row vectors, convert it to a rank3. 
    Inputs:
        ifgs_r2 | rank 2 array | ifgs as row vectors 
        mask | rank 2 array | to convert a row vector ifg into a rank 2 masked array        
    returns:
        phUnw | rank 3 array | n_ifgs x height x width
    History:
        2020/06/10 | MEG  | Written
    """
    import numpy as np
    import numpy.ma as ma
    
    from licsalert.aux import col_to_ma
    
    n_ifgs = ifgs_r2.shape[0]
    ny, nx = col_to_ma(ifgs_r2[0,], mask).shape                                   # determine the size of an ifg when it is converter from being a row vector
    
    ifgs_r3 = np.zeros((n_ifgs, ny, nx))                                                # initate to store new ifgs
    for ifg_n, ifg_row in enumerate(ifgs_r2):                                           # loop through all ifgs
        ifgs_r3[ifg_n,] = col_to_ma(ifg_row, mask)                                  
    
    mask_r3 = np.repeat(mask[np.newaxis,], n_ifgs, axis = 0)                            # expand the mask from r2 to r3
    ifgs_r3_ma = ma.array(ifgs_r3, mask = mask_r3)                                      # and make a masked array    
    return ifgs_r3_ma
    