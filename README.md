LiCSAlert is an algorithm to detect changes in the signals at sub-aerial volcanoes that are imaged by the Sentinel-1 satellites.  

It uses [LiCSBAS](https://github.com/yumorishita/LiCSBAS) to create time series, which in turn uses the interefrograms that are automatically created by [LiCSAR](https://comet.nerc.ac.uk/comet-lics-portal/).  Time series are then divided into a 'baseline' and 'monitoring' stage.  During the baseline stage, [ICASAR](https://github.com/matthew-gaddes/ICASAR) is used to determine which signals are normally present at a volcano.  During the monitoring stage, changes in these signals are tracked, to determine if either a) the behaviour of a signal has changed or b) a new signal has entered the time series.  

The algorithm is detailed fully in [Gaddes et al., 2019](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JB017519).  

The figure below shows the results when the algorithm is applied to a time series of Sentinel-1 data that images Sierra Negra volcano (Galapagos Archipelago) before the 2018 eruption.  Note that the caldera floor uplift is capture in the first component, and its acceleration flagged (the red bars):

![sierra_negra_resized](https://user-images.githubusercontent.com/10498635/76204672-c007bf80-61f0-11ea-8ba7-86dc91a5e0f0.gif)

