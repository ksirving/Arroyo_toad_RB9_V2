# Arroyo_toad_RB9_V2

# Arroyo_toad_RB9
modelling arroyo toad vulnerability in response to altered flow

Data from:

# Files for Arroyo Toad Distribution Modeling work - folder: *AT DistributionModeling_Final*

*Mike Treglia - mtreglia@gmail.com/michael.treglia@tnc.org*

**Not Intended for Broad Distribution. Please reach out to Mike before sharing**

This folder contains R code, Data, and Results files associated models the
distribution of the Arroyo Toad presented and described in:
* Treglia ML, Fisher RN, Fitzgerald LA (2015) Integrating Multiple Distribution
  Models to Guide Conservation Efforts of an Endangered Toad. PLoS ONE 10(6):
  e0131628. https://doi.org/10.1371/journal.pone.0131628; and
* Treglia, M. L. (2014). Multi-scale Conservation in an Altered Landscape: The
  Case of the Endangered Arroyo Toad in southern California (Doctoral
  dissertation).

The Folders `Current` and `Historic` contain files associated with developing
the *Current* and *Potential* models of arroyo toad distribution. Of note - at
the highest level within each of those folders is a `.R` file that contains the
main code used to develop the models. Much of the code was adapted from code
shared in a workshop by
[Jeff Evans](https://evansmurphy.wixsite.com/evansspatial/about_us) at the
US-IALE meeting in 2013, and the respective work this builds on is referenced in
the manuscript. There is code near the end of each of those files to write
output - basically, the full set of code was run 10 times, with that code for
writing output adjusted to write each run to a respectively numbered folder. The
results of each run, for the Current and Potential models, were averaged for the
final results. Note - some of the code leveraged may now be incorporated into R
packages developed by Jeff suchs as the
[`rfUtilities` package](https://github.com/jeffreyevans/rfUtilities).

The results from each run are in
`Current/randomForests/PARTITIONING/Cur_10x_Final_Feb2014` and
`Historic/randomForests/PARTITIONING/Hist_20x_Final_1981`, respectively. Within
those folders are the folders for each run, files with assessment information
about model performance and results of averaged rasters from each run.

Within these `Current` and `Historic` folders are subfolders with input data
(`Current/randomForests/PARTITIONING/DATA3` and
`Historic/randomForests/PARTITIONING/DATA2_1981`) - the presence (Potential
model) and pesence/absence (Current model) data were aggregated to centroids of
200 m pixels that were the result of rasterizing National Hydrographic Database
flowlines, and derived from myriad sources described in the manuscript and
dissertation. The `.grd` and `gri` files are for raster data (used in the R
`raster` package), with the environmental data used reduced to principal
component layers.

The `FullData` folder contains the full set of data for the 200 m pixels, prior
to PCA transformation.
`FullData/200mCells_FullData_PresAbs_Complete_ThinnedCols.csv` has the original
data, and the shapefile (and associated files)
`FullData/200mCells_PresAbs_2005_Final.shp` is a vector representation of the
pixels. (Note - the shapefile was has PCA values for the Current model;
`FullData/CurrentGridFeb14.grd` is the raster representation of the
PCA-transformed data). The CSV file can be joined to the shapefile based on the
common values of *ID2* in both files. The code used to get from the full set of
variables to the PCA-transformed data is in `FullData/PCA_Dataset_Final.R`.
There are columns in the `.csv` file that do not have explicit names (they were
*Final[ #]*, named by the software in which I did the overlay of the landsat
imagery with these grid cells). It would take some digging/experimentation to
figure out what all of the columns refer to, but the ones used in the model were
ultimately selected in that `.R` file and named appropriately.
