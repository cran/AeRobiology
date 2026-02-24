# AeRobiology version 2.0

### Changes in AeRobiology version 2.0.2

* Authors Information updated.
* CITATION file updated.
* Removed non-standard things in the check directory

### Changes in AeRobiology version 2.0.1

* Bug solved in interpollen(method="tseries")
* CITATION file added with the manuscript information: *Rojo, J., Picornell, A., Oteros, J. (2019) AeRobiology: the computational tool for biological data in the air. Methods in Ecology and Evolution. DOI.: 10.1111/2041-210X.13203*
* Faster functions plot\_hour() and plot\_heathour(), when the argument *locations = FALSE*

## Major changes in version 2.0

* Two new functions for hourly data processing added: **plot\_hour()** and **plot\_heathour()**.
* A new data set about 3-hourly concentrations from an automatic pollen system (BAA500) have been added: **data("POMO\_pollen")**.
* Added a tutorial in RPubs:  **http://rpubs.com/Picornell/AeRobiology** and a vignette.
* **Interpollen()** method = "neighbour" colnames homogenized.
* Homogenized the name of the functions: **ps\_plot -> plot\_ps; trend\_plot -> analyse\_trend; trends -> plot\_trend; qualitycontrol -> quality\_control; hour\_plot -> plot\_hour; summary\_plot -> plot\_summary; summary\_normplot -> plot\_normsummary**
* Homogenized the name of the data sets: **munich** -> **munich\_pollen**; **POMO\_pollen**
* Additional **examples** added.
* Function **trend\_plot()** updated changing the visualization of the slopes of the trends in end.date and start.date.
* Deleted **global assignment** (<<-) in all functions.
* Eliminated import of **depreciated packages**: *reshape2* and *plyr*. Substituted by active packages *dplyr*, *tidyr* and *purr*.
* **export.result** and **export.plot** arguments have been eliminated in functions when unnecessary (e.g. interpollen()) and set to *FALSE* by default in all functions.
* Name of arguments have been homogenized among functions (e.g. **plot -> export.plot**)
* Plots displayed by default, but not exported.
* Eliminated **prints** of progress messages when unnecessary (e.g. functions interpollen() and plot\_trend()). In some functions such as calculate\_ps(), this functionality is very relevant.
* New argument **"result"** added to some functions (e.g. interpollen(), calculate\_ps(), plot\_hour(), analyse\_trend()). This argument allows the user to decide which result want to be produced by the function e.g. (result = "table" or result = "plot").
* Input data are always formatted to data.frame (**data <- data.frame(data)**).
* Fixed bug in function interpollen()



# AeRobiology version 1.0

### Changes in AeRobiology version 1.0.3

* Import the packages: *purrr*, *colorspace* and *httpuv*
* Function trend\_plot() updated with extra information and better visualization of "significant" argument.
* Change the "maintainer" to only one of the authors, maintainer role will rotate among the authors.

### Changes in AeRobiology version 1.0.2

The changes in specific functions are entitled as updates.function()

#### updates.general

* Three typos in the documentation (pg. 27, 29, 33, 37).

#### updates.iplot\_abundance()

* added a new argument: *exclude*. *exclude* is a *character* string vector with the names of the pollen types to be excluded from the plot.

#### updates.calculate\_ps()

* Updated documentation with default values of the arguments.
* Updated documentation about "natural year". Natural year is also known as "Calendar year", i.e. the period from 1.January to 31.December.

#### updates.iplot\_pollen()

* Fixed bug of date recognition.

#### updates.iplot\_pollen()

* Fixed bug of date recognition.
