# R-package Representative Days (ReprDays)

R-functions and corresponding script file `script.r` for interactive use, and data files.

Software companion for journal article:

**Martin Densing, Yi Wan (2021).** Low-dimensional scenario generation method of solar and wind availability for representative days in energy modeling, *Applied Energy*, in press. https://doi.org/10.1016/j.apenergy.2021.118075


## Documentation of package

See `representative-days.pdf`, first section: `ReprDays`.


## Use of package

The files are arranged in the folder structure of an R package, which allows to use the package in two ways.

1. *Manually* (strongly recommended). Just load `representative-days.r` in R, which contains only functions. Then execute interactively, i.e. line-by-line, `script.r`. In this script-file, the input-data file-paths have to correspond to the path to the CSV-files, which are provided in the folder `inst\extdata`. 
2. *Build the R package*. 

    1. Download and leave unchanged the folder structure of the package sources (files are: the previously files + package description file `DESCRIPTION`). 
	2. Set the R working directory to the root of the soure files, and execute just three commands:
	
	   1. `library(devtools)` # This sounds trivial, but you may have to install first `rtools`, which can be arduous on MS Windows
	   2. `document()` # create the documentation .Rd-files (accessible  with `?` in R) in an automatically generated folder `man` 
	   3. `build(manual = T)` # builds the zipped R package, including the documentation
	   4. After installing the package in R, the path to the external CSV-files can be queried in R with `system.file("extdata", package = "reprDays")`, which gives also the path to `script.r`, which then can be interactively opened and executed.

## Data

Generation and load data is from the transparency platform of ENTSO-E (https://transparency.entsoe.eu/); data can be used freely if the platform is cited.
