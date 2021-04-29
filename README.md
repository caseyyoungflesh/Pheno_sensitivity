# Pheno_sensitivity


[![Zenodo DOI](https://zenodo.org/badge/290607887.svg)](https://zenodo.org/badge/latestdoi/290607887)


Code for Youngflesh et al. 2021 *Nature Ecology and Evolution*


This repository contains code to model bird arrival as a function of green-up and green-up as a function of time. The companion repository, which contains code to quantify bird arrival can be found [here](https://github.com/phenomismatch/Bird_Phenology).


**Associated publications:**

Youngflesh, C., Socolar, J., Amaral, B.R., Arab, A., Guralnick, R.P., Hurlbert, A.H., LaFrance, R., Mayor, S.J., Miller, D.A.W., Tingley, M.W. 2021. [Migratory strategy drives species-level variation in bird sensitivity to vegetation green-up.](https://www.nature.com/articles/s41559-021-01442-y) __*Nature Ecology and Evolution*__



**Repository structure:**

* `Data/` (ignored) - Data relevant for project
  * `arrival_master_<DATE>` - output from `5-extract-arr-dates.R` in companion repo
  * `environent`
    * `RAW` - raw green-up data
      * `MCD12Q1` - MCD12Q1 product data
      * `MCD12Q2` - MCD12Q2 product data
    * `processed` - processed green-up date
  * `hex_grid_crop` - study area
  
* `Scripts/` - scripts to run analyses
  * `1-process-gr/`
    * `1a-landcover.R` - create land cover mask
    * `1b-greenup.R` - extract green-up data over study area
  * `2-gr-SVC.R` - script to run green-up as a function of time
  * `3-arr-gr-SVC-sens.R` - script to run bird arrival as a function of green-up

* `Results/` (ignored)
