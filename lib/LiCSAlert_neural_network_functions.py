#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 16:43:03 2021

@author: matthew
"""



def sources_though_cnn(sources_r2, mask, vudl_net_21_path):
    """Given a matrix of interferograms (or ICs) as row vectors and their mask, convert them to the correcto form for use with VUDL-net-21, 
    and then detect and locate deformation using the model.  
    Inputs:
        sources_r2 | rank 2 array | sources as row vectors.  
        mask | rank 2 array of boolean | to convert a row vector to a rank 2 masked array.  
        vudl_net_21_path | Path (or string) path to Keras model.  
        
    Returns:
        Y_class | n_sources x 3 | classes as softmax ouputs (ie probabilities of each)
        Y_loc | n_soruces x 3 | location of deformaiton in ifg (0, 0, 0, 0 if None)
        
    Hitory:
        2021_06_09 | MEG | Written
    """
    import numpy as np
    import numpy.ma as ma
    from skimage.transform import resize as sk_resize
    import keras
       
    # 1: resize and rescale the data to be suitable for use with VUDL-net-21
    ifgs_r3 = r2_to_r3(sources_r2, mask)                                                                                    # row vectors to rank 3 masked array
    x_scale = ifgs_r3.shape[2] / 224                                                                                        # determine how much we are going to resize the images by in x
    y_scale = ifgs_r3.shape[1] / 224                                                                                        # and in y
    ifgs_r4 = ma.repeat(ifgs_r3[:,:,:,np.newaxis], 3, axis = 3)                                                             # rank 4 masked array with 3 channels (data just repeated across them).  
    ifgs_r4_rescale = custom_range_for_CNN(ifgs_r4, 2, mean_centre = True)                                                  # data could have values not suitable for CNN (e.g. [-10 5], so rescale) 
                                                                                                                            # Note that we rescale so that signals are mean centered, and either the min or max touches -1 or 1
    ifgs_r4_rescale = ma.filled(ifgs_r4_rescale, fill_value = 0)                                                            # convert to normal array, filling with 0s
    ifgs_r4_rescale_resize = sk_resize(ifgs_r4_rescale, (ifgs_r4_rescale.shape[0], 224, 224, ifgs_r4_rescale.shape[3]))     # resize so the correct number of pixels for use with CNN
    
    # 2: Load the model
    classification_model = keras.models.load_model(vudl_net_21_path)
    Y_class, Y_loc = classification_model.predict(ifgs_r4_rescale_resize, verbose = True)                                   # pass through the CNN
    
    Y_loc[:,0] *= x_scale                                                                                                   # rescale x centre
    Y_loc[:,2] *= x_scale                                                                                                   # and half width
    Y_loc[:,1] *= y_scale                                                                                                   # same for y
    Y_loc[:,3] *= y_scale
    
    return Y_class, Y_loc


#%%

def custom_range_for_CNN(r4_array, custom_range, mean_centre = False):
    """ Rescale a rank 4 array so that each channels image lies in custom range
    e.g. input with range of [-5 15] is rescaled to [-125 125] or [-1 1] for use with VGG16 
    
    The function used in the original (paper) implementation of VUDL-NET-21
    
    Inputs:
        r4_array | r4 masked array | works with masked arrays?  
        custom_range | int or float | the data will lie between -custom_range/2 and custom_range/2
    
    2019/03/20 | now includes mean centering so doesn't stretch data to custom range.  
                Instead only stretches until either min or max touches, whilst mean is kept at 0
    """
    import numpy as np
    
    def expand_to_r4(r2_array, shape = (224,224)):
        """
        Calcaulte something for every image and channel in rank 4 data (e.g. 100x224x224x3 to get get 100x3)
        Expand new rank 2 to size of original rank 4 for elemtiwise operations
        """
        import numpy as np
        
        r4_array = r2_array[:, np.newaxis, np.newaxis, :]
        r4_array = np.repeat(r4_array, shape[0], axis = 1)
        r4_array = np.repeat(r4_array, shape[1], axis = 2)
        return r4_array

    
    if mean_centre:
        im_channel_means = np.mean(r4_array, axis = (1,2))                                                  # get the average for each image (in all thre channels)
        im_channel_means = expand_to_r4(im_channel_means, r4_array[0,:,:,0].shape)                          # expand to r4 so we can do elementwise manipulation
        r4_array -= im_channel_means                                                                        # do mean centering    

    im_channel_abs_max = np.max(np.abs(r4_array), axis = (1,2))                                             # get the maximum for each image and each channel, in either the positive or negative direction.  
    im_channel_abs_max = expand_to_r4(im_channel_abs_max, r4_array[0,:,:,0].shape)                          # expand to rank 4 (ie repeat the values for each image and channel across each pixel in the image)
    r4_array = (custom_range/2)* (r4_array/im_channel_abs_max)                                              # rescale so that either the largest or smallest value in the image is 1, then multiply by the range/2, so that value is now range/2
                                                                                                            # 
    return r4_array       


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
    from small_plot_functions import col_to_ma
    
    n_ifgs = ifgs_r2.shape[0]
    ny, nx = col_to_ma(ifgs_r2[0,], mask).shape                                   # determine the size of an ifg when it is converter from being a row vector
    
    ifgs_r3 = np.zeros((n_ifgs, ny, nx))                                                # initate to store new ifgs
    for ifg_n, ifg_row in enumerate(ifgs_r2):                                           # loop through all ifgs
        ifgs_r3[ifg_n,] = col_to_ma(ifg_row, mask)                                  
    
    mask_r3 = np.repeat(mask[np.newaxis,], n_ifgs, axis = 0)                            # expand the mask from r2 to r3
    ifgs_r3_ma = ma.array(ifgs_r3, mask = mask_r3)                                      # and make a masked array    
    return ifgs_r3_ma

#%%

def centre_to_box(centre_width):
    """ A function to convert from centre and width notation to start and stop notation (in both x and y directions).   
    Inputs:
        centre_wdith
    Returns:
    History:
        2020_10_28 | MEG | Wrote the docs.  
        """
    x_start = centre_width[0] - centre_width[2]
    x_stop = centre_width[0] + centre_width[2]
    y_start = centre_width[1] - centre_width[3]
    y_stop = centre_width[1] + centre_width[3]
    return [x_start, x_stop, y_start, y_stop]


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
   
#%%


     