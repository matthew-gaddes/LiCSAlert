# Introduction

LiCSAlert is an algorithm to detect changes in the signals at sub-aerial volcanoes that are imaged by the Sentinel-1 satellites.  The algorithm is detailed fully in [Gaddes et al., 2019](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JB017519) but, briefly:
- The time series is divided into baseline interferograms and monitoring interferograms (e.g. the first third of a time series and the remainder).  
- The [ICASAR](https://github.com/matthew-gaddes/ICASAR) algorithm is run on the baseline interferograms to find a set of spatially independent sources that are thought to represent the latent deformation and atmosphere signals present at a volcano.  
- The temporal behaviour of these signals is calculated through performing a simple inversion to fit each baseline interferogram with them.  
- During monitoring, the same inversion is performed to fit new interferograms.  
    - If the way in which the spatial sources is used changes significantly, the volcano is flagged as having entered a period of unrest (e.g. acceleration of uplift).    
    - If the ability of the spatial sources to fit the interferograms (i.e. the residual) changes significantly, the the volcano is flagged as having entered a period of unrest (e.g. a previously unseen deformation signal appears).  

The figure below shows the results when the algorithm is applied to a time series of Sentinel-1 data that images Sierra Negra volcano (Galapagos Archipelago) before the 2018 eruption.  Note that the caldera floor uplift is capture in the first component, and its acceleration flagged (the red bars).  Additionally, the residual (lowest plot) increases before the eruption, suggesting the spatial pattern of the deformation has changed slightly.  

![sierra_negra_resized](https://user-images.githubusercontent.com/10498635/76204672-c007bf80-61f0-11ea-8ba7-86dc91a5e0f0.gif)

# Installation
A suitable python environment can be created using the conda command:<br>
<code>conda env create --file LiCSAlert.yml</code>

The ICASAR package will also be needed, and the <code>ICASAR_path</code> argument updating to point to it.  


# Batch mode usage
Batch mode usage is simpler than monitoring mode as it does not automatically updated the time series using LiCSBAS when new LiCSAlert products are available.  Simply prepare unwrapped incremental interferograms in your software of choice, mask pixels that you do not wish to include (e.g. water bodies, incoherent areas etc.), and flatten each interferogram to a 1D vector that only contains values for the unmasked pixels.  

There are three groups of inputs:

1) <code>ICASAR_path</code>  | the path to a local copy of ICASAR

2) <code>LiCSAlert_settings</code> 
    - <code>n_baseline_end</code>  |  The number of interferograms to use in the baseline stage.  E.g. 30, to use the first 30, and any remaining will be used for monitoring.  
    - <code>out_folder</code>  |  Where to store the outputs of LiCSAlert
    - <code>run_ICSAR</code>  |  True or False.  If it has been run before, setting this to False will try to load the previous results.  
    - <code>intermediate_figures</code>  |  If True, a figure is made for each time step, but if False, a single figure is made for the whole time series.  Intermediate figures can be useful for making .gif animations.  See the example for the differences in the outputs (and runtime!).    
    - <code>downsample_run</code>  |  Downsampling the data can speed up runs.  
    - <code>downsample_plot</code>   |  Downsampling the data for plotting can speed up making figures.  Note that this is applied after the downsample_run command, so is compound (i.e. 0.5 for downsample_run and 0.5 for downsample_plot produces a final downsampling of 0.25 for the plotted signals).  

3) <code> ICASAR_settings</code>
  These are explained in the [ICASAR wiki](https://github.com/matthew-gaddes/ICASAR/wiki/03-Inputs-and-Tunable-parameters).  


# Monitoring mode usage

It uses [LiCSBAS](https://github.com/yumorishita/LiCSBAS) to create time series, which in turn uses the interefrograms that are automatically created by [LiCSAR](https://comet.nerc.ac.uk/comet-lics-portal/). A simple example is outside the scope of this repository.  
