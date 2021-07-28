# 21.Hu.JClim
Public Code for "Approximating the Internal Variability of Bias-Corrected Global Temperature Projections with Spatial Stochastic Generators" by Wenjing Hu and Stefano Castruccio. More details can be found at https://doi.org/10.1175/JCLI-D-21-0083.1

**Spatial_model.m**<br />
This code fits a spatial linear model with considering the spatial correlation among regions and monthly variance, and then estimate the sparse correlation matrix with graphical lasso.

**bias_correction.m**<br />
This code fits the spatial model above for Large Ensemble and MERRA2, and then perform bias correction with estimated mean trend and sparse correlation matrix.

**Data**<br />
Large ensemble is downloaded from https://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html and aggregated into 47 regions at monthly level.<br />
MERRA2 is downloaded from https://disc.gsfc.nasa.gov/datasets?project=MERRA-2 and aggregated into 47 regions at monthly level.<br />
ERA5 is downloaded from https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5 and aggregated into 47 regions at monthly level.


