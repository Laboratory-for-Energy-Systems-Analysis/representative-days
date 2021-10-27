# Representative days

R-function and script.r file for paper.


The files are arranged in the folder structure of a proper R package.
Hence, you can use the package in two ways.
1. *Manually* (strongly suggested). Load `representative-days.r`, which contains only functions, in R and execute interactively, i.e. line-by-line, `script.r`, by setting the input-data file-paths in `script.r` to the path to the data in `inst\extdata`. documentation is in `representative-days.pdf`.
2. *Build the R package*. 
    1. Download the whole source file structure of the package (which has a files the previously files + package description file `DESCRIPTION`. 
	2. Set the R working directory to root of soure files, and execute just three commands:
	   1. `library(devtools)` # This sounds trivial, but you may have to install first `rtools`, which can be a little bit a pain on MS Windows
	   2. `document()` # builds the documentation (accessible  with `?` in R) in an automatically generated folder `man` 
	   3. `build(manual = T)` # builds a zipped R package, including documentation
	   4. After installing the package in R.
	   5. The path to the data can be queried in R with `system.file("extdata", package = "reprDays")`, which gives also the path to `script.r`, which then can be interactively executed.

