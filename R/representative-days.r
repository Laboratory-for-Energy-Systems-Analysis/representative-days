
#' reprDays: A package for computating representative days
#'
#' Package for article \emph{Low-dimensional scenario generation method of solar and wind availability for representative days in energy modeling}, Applied Energy, in press. \url{https://doi.org/X}
#' 
#' The package provides different categories of functions:
#'  
#' @section Import:
#' \itemize{
#'   \item Read hourly solar and wind generation data, and optionally load data
#'   \itemize{
#'     \item Current data is from ENTSO-E's transparency platform. Data from countries DE and AT is converted to hourly. Data should all be in CET time-zone.
#'    \item Provided data is of year 2017-2019 (State: Dec 2021)
#'    \item Absolute wind and solar is converted to availabilities by using an auxiliary capacity data file
#'    \item Absolute load data is normalized by dividing with high quantiles. 
#' }
#' \item Plot time-series of solar, wind (and load). Plot daily profiles per season and per country
#' \item Correlations of wind, solar, and load
#' }
#'
#' @section PCA factor model:
#' \itemize{
#' \item Execute the PCA on the data vector
#' \item Make a principal compoment factor analysis
#' \item Generate representaive-days scenarios using the factor models,and plot the scenarios
#'}
#' @section Analysis of accuracy of scenarios:
#' Compare the scnearios with empirical data and incraeasing scenario number:
#' \itemize{
#' \item Compare correlations
#' \item Compare ECDF with empirical data
#' \item Compare availability per season
#' }
#' @section Hierarchical clustering using medoids:
#' Generate representative days by clustering. Given a season, single and cross-regional scenarios can be generated
#'
#'
#' @section Cross-regional scenarios:
#' \itemize{
#' \item Analysis of daily wind and solar correlations across regions (countries)
#' \item Estimation of joint extreme solar and wind availabilities
#' \item Fitting guassian copulas and t-copulas to daily cross-regional wind and solar data
#' \item Random sampling of copulas
#' \item Given a season, join daily cross-regional random samples with (aggregated from hourly to daily) regional scenarios to provide consistent hourly scenarios across regions 
#' }
#' 
#' @section Auxiliary function (selection):
#' \itemize{
#'   \item  Aggregate availabilities to 96 (24*4) load-period availabilities and write to "avlBEM.csv" for use in BEM
#' }
#' 
#' @docType package
#' @name reprDays
NULL


# ESS in Emacs specifics:
#
# Turn off package namespace evaluation (throws errors instead of warnings):
# (setq ess-r-package-auto-activate nil)
# (setq ess-r-package-auto-set-evaluation-env nil)
#
# Interactive:
# 1. R-windows: Turn off minor-mode package
# 1. R-script: `C-c C-t C-s`


library(xts)
library(spatstat)     # for plots of ecdf() of single-country scenarios
library(fitdistrplus) # for fitting Weibull distribution to first PCA factor
library(corrplot)     # for correlation plots in cross-country analysis
library(copula)       # for cross-regional analysis






### Helper functions and global variabes






#' Set global parameters
#'
#' Initialize \code{ctyNames}, \code{seasNames}, \code{seasons}
#'
#' @param cty Vector of country names
#' @param seas Vector of season names
#' @param month.idx Index of months of the seasons {1,\ldots,12} 
#'
#' @return Global parameters:
#' \describe{
#'  \item{\code{ctyNames}}{\code{c("AT", "CH", "DE", "FR", "IT")}}
#'  \item{\code{seasNames}}{\code{c("WI", "SP", "SU", "FA")}}
#'  \item{\code{seasons}}{Matrix: columns =
#'     "WI, "SP", "SU", "FA"; rows = month-index (1,\ldots,12)}
#' }
set.cty.seas <- function(cty = c("AT", "CH", "DE", "FR", "IT")
                        ,seas =  c("WI", "SP", "SU", "FA")
                        ,month.idx = c(c(12,1,2), 3:5, 6:8, 9:11)
                         )
{
    ctyNames <<- cty
    seasNames <<- seas
    seasLength = 12 %/% length(seasNames)
    seasons <<- matrix(nrow = seasLength
                      ,ncol = length(seasNames)
                      ,month.idx
                      ,dimnames=list(paste0("m",seq(seasLength)),
                                     seasNames))
}


#' Extract hours and months
#' 
#' @param x Time-date, can be a vector of dates (format: POSIXct)
#'
#' @return Hour of day (= 0-23), or index of month (= 1-12)
#'
#' @examples
#' hour(as.POSIXct("2012-10-19"))
#'
hour <- function(x) as.numeric(format(time(x), "%H"))
#' @rdname hour
month <- function(x) as.numeric(format(time(x), "%m"))



#' Vector of load-period tags in a season
#'
#' @param seas Season
#' @param n Number of load periods in a day (default: 24)
#' 
#' @return
#' \code{c("WI-D-01",\ldots,"FA-D-24")}
#' 
load.24 <- function(seas, n=24) paste0(seas,"-D-",sprintf("%02d", 1:n))







### Load data






#' Load hourly wind, solar data (optionally also load)
#'
#' \itemize{
#'  \item The function works currently with only country DE having also offshore wind
#'  \item The wind and solar files must be named like "DE_Generation_2019 (ENTSO)HOURLY.csv"
#'  \item The load files must be named like "DE_Load_2019 (ENTSO)HOURLY.csv"
#' }
#' 
#' @param  fn.path  Path to ENTSO country files
#' @return x Time series (xts class) containing hourly data with column names
#' \itemize{
#' \item "ATsolar", "ATwind", \ldots, "ITwind", "DEOFFwind", "DEONOFFwind"
#' \item If also \emph{load} is read:  "ATsolar", "ATwind", "ATload",\ldots
#' }
#' 
x.read <- function(
                   fn.path = "Data/"
                  ,yearRead = c(2017,2018,2019)
                  ,load = FALSE
                  ,colSol = "Solar."
                  ,colWiOn = "Wind.Onshore."
                  ,colWiDEOff = "DEWind.Offshore."
                  ,colLoad = "load" 
                   ) {
    for (yr in yearRead) {
        for (cty in ctyNames) {
            fn  = paste0(fn.path, cty, "_Generation_", yr, " (ENTSO)HOURLY.csv")
            y = read.csv(fn, stringsAsFactor=FALSE)
            colnames(y) = paste0(cty, colnames(y))
            if (load) {
                fn2 = paste0(fn.path, cty, "_Load_", yr, " (ENTSO)HOURLY.csv")
                y2 = read.csv(fn2, stringsAsFactor=FALSE) 
                colnames(y2) = paste0(cty, colnames(y2))
                # browser()
            }
            # get rid of date column in yearly data "y", and then add "y" as column to "xYr"
            xYr = if(cty == first(ctyNames)) y else cbind(xYr, y[,-1]);
            if (load) { xYr = cbind(xYr, y2[,-1]) # y2[,-1] looses name because 1 col -> vec
                names(xYr)[length(names(xYr))] = paste0(cty,"load")
            }
        }
        x = if (yearRead[1]==yr) xYr else rbind(x, xYr); # append "xyR" as rows to "x"
    }
    # Do not use directly dates in file because DST not handled correctly. Use instead: seq()
    # seq() produces slightly offest point-in-times such that the doubled hours in autumn are different
    frmt = "%Y-%m-%d %H:%M"   #format of date
    idx = seq(
        from = as.POSIXct(first(x[,1]),format=frmt),
        to =   as.POSIXct(last(x[,1]), format=frmt), by="hour")
    x = xts(x[,-1], idx) # substitue new date column, convert df to xts

    # we only want solar, wind, and load columns
    if (!load) {
        col = c( paste0(rep(ctyNames, each=2), c(colSol,colWiOn)),
                colWiDEOff )
    } else {
        col = c( paste0(rep(ctyNames, each=3), c(colSol,colWiOn,colLoad)),
                colWiDEOff )
    }
    
    x = x[,col]
    if(sum(is.na(x))>0) {
        x[which(is.na(x))]
        stop("NA Values in time series x. Plese remove in loaded .csv-file.") 
    }

    #We have onshore and offshore wind: Make a column with total wind
    print("check colnames odering, because it will be altered! Ordering in x:")
    print(colnames(x))
    x = cbind(x, x[,"DEWind.Onshore."]+x[,"DEWind.Offshore."])
    colnames(x) =
        if (!load)
            c(
                "ATsolar", "ATwind",
                "CHsolar", "CHwind",
                "DEsolar", "DEwind",
                "FRsolar", "FRwind",
                "ITsolar", "ITwind",
                "DEOFFwind","DEONOFFwind")
        else
            c(
                "ATsolar", "ATwind", "ATload",
                "CHsolar", "CHwind", "CHload",
                "DEsolar", "DEwind", "DEload",
                "FRsolar", "FRwind", "FRload",
                "ITsolar", "ITwind", "ITload",
                "DEOFFwind","DEONOFFwind")

    return(x)
}


#write.csv(
#    rbind(colSums(x['2017']),colSums(x['2018']),colSums(x['2019'])),
#    "temp.csv")



#' Convert x into availability factors (divide by capacity):
#'
#' Generation is divided by end-of-year capacity. Load is dvidided by
#' a high quantile of each year. Currently only DE has also Offshore
#' wind
#'
#' @param x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.; as returned by \code{\link{x.read}}
#' @param fn File name of solar and wind capacities in csv format,
#'     with rows as years (e.g. 2017), and columns as country and tech
#'     (e.g. "ATsolar")
#'
x.normalize <- function(x
                       ,yearRead = c(2017,2018,2019)
                       ,fn = "SolarWindcapacity.csv"
                       ,load = FALSE
                       ,loadQuantile = 0.96
   ) {
    cap = read.csv(fn, row.names=1)
    for (cty in ctyNames) {
        for (yy in yearRead) {
            y = as.character(yy)
            for (t in c("solar", "wind")) {
                z = paste0(cty,t)
                x[,z][y] = x[,z][y] / cap[y,z]
                # DEwind is onshore, treat offshore and total:
                if (cty == "DE" & t == "wind") {
                    for(cty2 in c("DEOFF","DEONOFF")){
                        z = paste0(cty2,t)
                        x[,z][y] = x[,z][y] / cap[y,z]
                    }
                }
            }  
        }
    }

    if (load)
        for (cty in ctyNames) {
            for (yy in yearRead) {
                y = as.character(yy)
                for (t in c("load")) {
                    z = paste0(cty,t)
                    x[,z][y] = x[,z][y] /  quantile(x[,z][y],loadQuantile)
                } 
            }
        }
  return(x)
} # end of function abs2avl




#' Plot wind, solar (and load, if read-in) per country and per year
#'
#' @param  x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.; as returned by \code{\link{x.normalize}} or \code{\link{x.read}} 
#' @param cty Country
#' @param yr Year (or a range of years in xts-notation)
#'
#' @examples
#' \dontrun{
#' plotSolarWind(x,"CH","2017\2019")
#' plotSolarWind(x,"CH","2018")
#' }
#' 
x.plot <- function(x, cty, yr) {
    # Convert to zoo object to allow plotting of multivariate time-series
    plot.zoo(x[yr,grepl(cty, names(x))],
             main = cty,
             # ylab = c("solar (avl.)", "wind (avl.)"),
             xlab = paste(yr, collapse="/")
             )
}




#' Select for a country the wind: Onshore, Offshore, or On+Off?
#'
#' 
#' @param  x Time series (xts class) containing hourly data with column
#'     names \ldots, "DEONwind", "DEOFFwind", "DEONOFFwind"; as returned by \code{\link{x.normalize}} or \code{\link{x.read}} 
#' @param windType 1: onshore, 2: offhsore, 3: total
#' @param cty Country (default: "DE")
#' 
#' @return Time series with selected column (named "DEwind")
#' 
sel.wind <- function (x
                     ,windType
                     ,cty = "DE"
                      ,load = FALSE)
{
    offw = paste0(cty, "OFFwind")
    onw = paste0(cty, "wind")
    totw = paste0(cty, "ONOFFwind")
    cols =  !names(x) %in% switch(windType,
                                  c(offw, totw),
                                  c(onw, totw),
                                  c(onw, offw))
    
    x = x[, cols]
    print(names(x))
    # just name it wind
    regExp = paste0(cty, '[A-Z]*wind')
    colnames(x) = gsub(regExp, paste0(cty,'wind'), colnames(x))
    #colnames(x)[if(load) 8 else 6] = "DEwind" 
 return(x)   
}    




#' Correlation between solar, wind, and (optionally) load
#'
#' Prints correlation for different years. Currently, it must be three years.
#' If load is read-in: Print correlation-test of load/wind of first year
#' 
#' @param  x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.; as returned by \code{\link{sel.wind}} 
#' @param yr List of three years
#'
x.cor <- function(x
                 ,load = FALSE
                 ,yr = c('2017','2018','2019')
                  ){
    for (cty in ctyNames) {
        cat(cty, "cor solar wind",'\n')
        print( c(
            cor(x[yr[1]])[paste0(cty,"wind"),paste0(cty,"solar")],
            cor(x[yr[2]])[paste0(cty,"wind"),paste0(cty,"solar")],
            cor(x[yr[3]])[paste0(cty,"wind"),paste0(cty,"solar")]))
    }
    if (load)
        for (cty in ctyNames) {
            cat(cty, "cor load solar", '\n')
            print( c(
                cor(x[yr[1]])[paste0(cty,"load"),paste0(cty,"solar")],
                cor(x[yr[2]])[paste0(cty,"load"),paste0(cty,"solar")],
                cor(x[yr[3]])[paste0(cty,"load"),paste0(cty,"solar")]))
            cat(cty, "cor load wind", '\n')
            print( c (
                cor(x[yr[1]])[paste0(cty,"load"),paste0(cty,"wind")],
                cor(x[yr[2]])[paste0(cty,"load"),paste0(cty,"wind")],
                cor(x[yr[3]])[paste0(cty,"load"),paste0(cty,"wind")]))
        }
    if (load)
        for (cty in ctyNames) {
            cat(cty, "cor load wind",'\n')
            print(cor.test(x[yr[1],paste0(cty,"load")],
                           x[yr[1],paste0(cty,"wind")], method = c("pearson"))$conf)
        }

}





### Plot hourly patterns per season






#' Plot the seasonal (averaged) day-hour data
#'
#' The number of seasons is currently hardcoded: 4 seasons
#' 
#' @param  x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.; as returned \code{\link{sel.wind}} 
#' @param cols column to plot
#' @param yr list of years (default: \code{c("2017","2018","2019")})
#' @param xpd expand the dispay (required for plots of wind)
#' 
#' @examples
#' \dontrun{
#' windows(); plot24(x,"ATsolar"); dev.print(pdf, "ATsolar.pdf") 
#' windows(); plot24(x,"ATwind", TRUE)
#' }
plot24 <- function(x
                  ,cols
                  ,yr = c("2017","2018","2019") 
                  ,xpd = FALSE) {
    y = x24SeasonYears(x, cols, yr)
    if (xpd) par(xpd=TRUE, mar=par()$mar+c(0,0,0,7))
    plot(y, plot.type = "single",
         col = rep(1:4, each = 3), lty = 1:3, lwd = 2, 
         main = cols, ylab = "Availability", xlab = "hour of day")
    axis(1, at = time(y), labels = FALSE)
    if (xpd) {
        legend("right", inset=c(-0.4,0),
               colnames(y), lty=1:3, lwd=2, col=rep(1:4, each = 3), xpd=TRUE)
    } else
        legend("topleft", colnames(y), lty = 1:3, lwd=2, col = rep(1:4, each = 3))
    par(mar=c(5, 4, 4, 2) + 0.1, xpd=FALSE)  
}


#' Get hour and year of a date
#'
#' Helper function to plot seasonal day-hour value
#' 
#' @param t Date-time 
#' @return Hour and year (e.g. "02 2017")
#' 
hourOfDay <- function(t) format(t, "%H %Y")


#' Select columns from time series and simplify time index
#'
#' Select by year, by col-names, and reduce the time-index to the day-hour only 
#'
#' Helper function to plot seasonal day-hour value
#' 
#' @param z Time-series with time index like "23 2012"
#' @param cols Vector of column names to select
#' @param year Year (e.g. "2017")
#'   
#' @return Data in specified year and columns, where the time-index is the day hour
x24Seas <- function(z, cols, year) {
    z = z[grepl(year,time(z)), cols] # get data in year
    time(z) = substr(time(z),1,2) # take only the hour (2 digits) of time index
    return(z)
}


#' Calculate mean-hourly values per season over selected columns and years
#'
#' Helper function to plot seasonal day-hour value
#' 
#' @param  x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.; as returned \code{\link{sel.wind}} 
#' @param yr Years
#' @param cols Column-names of the time series to select (e.g. "ATsolar")
#'
#' @return List mean values over years and seasons. List elements have names like "SP2017", contain a time series, with columns "ATsolar", etc., and with rows the mean values over the day hours
#' 
x24SeasonYears <- function(x
                          ,cols
                          ,years = c("2017","2018", "2019"))
{
    # split into seasons
    # for each season: mean for each day-hour and year
    xSeasMeanDayhour = list() # list of seasonal time series
    for (s in seasNames) {
        z = x[(.indexmon(x)+1) %in% seasons[,s]] #.indexmonth starts at 0
        xSeasMeanDayhour[[s]] = aggregate(z, hourOfDay, mean)
    }
    # bind the years and season together over columns
    z = list()
    i=1
    for (s in seasNames)
        for (yr in years) {
            z[[paste0(s,yr)]] = x24Seas(xSeasMeanDayhour[[s]], cols, yr)
            i = i + 1
            }
    z = do.call(cbind, z)
    return(z)
}









### Prepare data for PCA (no. of columns = dimension of PCA)






#' Get data.frame of hourly wind+solar avl. per season in wide-format (hours are in columns)
#'
#' @param  x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.
#' @param cty Country 
#' @param s Seasons given as vector of month numbers (1-12)
#'
#' @return 
#' \itemize{
#'  \item 48-column data.frame of solar and wind for specified months
#'  \item if also load: 72-column df of solar, wind, and load)
#'  \item column names: "solarDE01", "solarDE02",...,"windDE24" ("loadDE24")
#'  \item Note: Length of columns can be different for different choices of seasons
#'}
#' 
#' @examples
#' \dontrun{
#' dayHourPerSeason("AT",3)
#' dayHourPerSeason("DE", 4:9) # summer halfyear
#' dayHourPerSeason("DE", c(1,2,3,10,11,12)) # winter halfyear
#' }
dayHourPerSeason <- function(x, cty, s, load = FALSE) {
    l = list()
    techList = if (load) c("solar", "wind", "load") else c("solar","wind"); 
    for (i in 0:23) {
        for (tech in techList) {
            h = paste0(tech, cty, sprintf("%02d", i+1))
            l[[h]] = as.numeric(x[hour(x)==i & (month(x) %in% s),
                                  paste0(cty,tech)]
                                )
            print(paste(length(l[[h]]),
                        "hours in data", h,
                        ",", tech, ", months=", paste(s, collapse=",")))
        }
    }
    l = as.data.frame(l)
    l = l[, order(names(l))]
}




#'Prepare data for PCA
#'
#' Assumption: The input time series is in CET, and we convert first to UTC.
#' The conversion is required to get rid of the missing DST hours in the night.
#' We do this before subsetting to a single year (else the last hour would be missed)
#'
#' @param  x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.
#' 
#' @return 
#' List over seasons. A list element contains data.frame of all countries in wide forma (the day hours are in different columns)
#' \itemize{
#'  \item data.frame has 240 columns (2*24 hours (solar+wind) * 5 countries) with names: "solarCH01", "solarCH02",\ldots,"windIT24"
#' }
#' 
x.24wide <- function(x
                    ,yr = '2017/2019')
{
    cat("time zone:"); print(tzone(x))
    x = x[yr]
    tzone(x) = 'UTC'
    time(x) = time(x)+3600 # shift a day back to start in year
    print("Data at first date:"); print(first(x))
    print("Data at last date:"); print(last(x)) # Checks

    xSWH = list() #  list over seasons
    for (s in seasNames) {
        print(s)
        print(seasons[,s])
        z = list() # list over countries
        for (cty in ctyNames)
            z[[cty]] = dayHourPerSeason(x, cty, seasons[,s])
        xSWH[[s]] = as.data.frame(do.call(cbind, z))
    }
    return(xSWH)
}




#' Data over seasons for all countries  
#'
#' Used to prepare data for external use, e.g. for external clustering. Convert time series first to UTC to avoid duplicated hours 
#' 
#' @param x Time series with columns "ATsolar","ATwind",\ldots,"ITwind" (xts-class)
#' 
#' @return
#' List of data.frames over seasons. Given a season, the columns of the data.frame are "ATsolar", \ldots, "ITwind"; the rows are the hourly data in the season
#'
l.stacked.All <- function(x)
{
    tzone(x) = 'UTC'
    time(x) = time(x)+3600 # shift a day back to start in year
    y = list()
    for (s in seasNames) {
        print(s)
        print(seasons[,s])
        z = list() # list over countries
        for (cty in ctyNames)
            z[[cty]] = daysPerSeason(x, cty, seasons[,s])
        y[[s]] = as.data.frame(do.call(cbind, z))
    }
    return(y)
}
    

#' @describeIn l.stacked.All
#'
#' @param  x Time series with columns "ATsolar","ATwind",\ldots,"ITwind" (xts-class)
#' @param cty,seas Country and season
#' 
daysPerSeason <- function(x, cty, seas)
{
    l = list()
    for (tech in c("solar","wind")) {
        h = paste0(tech, cty)
        l[[h]] = as.numeric(x[(month(x) %in% seas), paste0(cty,tech)])
        print(paste(length(l[[h]]),
                    "hours in data", h,
                    ",", tech, ", months=", paste(seas, collapse=",")))
    }
    l = as.data.frame(l)
    rownames(l) = time(x[(month(x) %in% seas), paste0(cty,tech)])
    l = l[, order(names(l))]
    return(l)
}




#' Contour plot of correlation matrix
#'
#' @param correl Correlation matrix 
#' @param title.cont Title of contour plot
#' @param cty Country or |-separated list of countries, used for labeling
#' @param s Season, used for labeling
#' @param triang Plot only upper triangle? (default: TRUE)
#' @param oldVer Plot old version with other color codes, whereas new version  asfixed levels? (default: TRUE)
#' 
#' @examples
#' \dontrun{
#' }
contour.plot <- function (correl
                         ,title.cont
                         ,cty
                         ,s
                         ,triang = TRUE
                         ,oldVer = FALSE
                          )
{
    m.r = 1
    if (triang) correl[lower.tri(correl)] = NA
    ctyVec = strsplit(cty, "\\|")[[1]]
    if (triang) par(mar=c(5,5,9,1)+0.1) else par(mar=c(5,5,2,1)+0.1)
    nColPos = if (oldVer) 25 else 10
    nColNeg = if (oldVer) 8 else 10
    xAxisPos = if (triang) 48 else 0
    xAxisSide = if (triang) 3 else 1
    titleLine = if (triang) 6 else 1
    keyLine = if (triang) 2 else 1
    nData = nrow(correl)
    labl = c(paste0("S h",1:24), paste0("W h",1:24))
    if (nData == 24) labl =  paste0("S h",1:24)
    labels.vecX = paste(rep(ctyVec,each=nData%/%length(ctyVec)),labl)
    #browser()
    if (nData == 24) labl =  paste0("W h",1:24)
    labels.vecY = paste(rep(ctyVec,each=nData%/%length(ctyVec)),labl)
    if (oldVer) {
        filled.contour(
            x=1:(nData*m.r),
            y=1:(nData*m.r),
            correl,
            plot.title = title(main = title.cont, line = titleLine),
            xlab="", ylab="", axes=FALSE,
            nlevels = 20,
            plot.axes = {
                axis(xAxisSide, at = 1:(nData*m.r), labels.vecX[1:nData],
                     las = 2, pos = xAxisPos)     
                axis(2, at = 1:(nData*m.r), labels.vecY[1:nData])
            },
            col = c(hcl.colors(nColNeg, "Oslo", rev=TRUE),
                    hcl.colors(nColPos, "Reds 3", rev=TRUE)),
            key.title = title(main = "Corr.:", line = keyLine),
            key.axes = axis(4, seq(-0.4, 1, by = 0.2)),
            asp = 1 # squared plot
        ) } else {
              filled.contour(
                  x=1:(nData*m.r),
                  y=1:(nData*m.r),
                  correl,
                  plot.title = title(main = title.cont, line = titleLine),
                  xlab="", ylab="", axes=FALSE,
                  levels =
                      c(-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,
                        0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                  plot.axes = {
                      axis(xAxisSide, at = 1:(nData*m.r),
                           labels.vecX[1:nData], las = 2, pos = xAxisPos)     
                      axis(2, at = 1:(nData*m.r),
                           labels.vecY[1:nData])
                  },
                  col = c(hcl.colors(nColNeg, "Oslo", rev=TRUE),
                          hcl.colors(nColPos, "Reds 3", rev=TRUE)),
                  key.title = title(main = "Corr.:", line = keyLine),
                  key.axes = axis(4, seq(-0.4, 1, by = 0.2)),
                  asp = 1 # squared plot
                  )
          }
    par(mar=c(5,4,4,2)+0.1)
    # divide solar and wind:
    if (nData > 24) {
        lines(c(-2, if (triang) 34.5/2-0.9 else 34.5),c(24,24), lty = 2, lwd=2) 
        lines(c(16,16),c(48,if (triang) 24 else 0), lty = 2, lwd=2)
    } 
}

#' Contour plot of correlations of data in wide format
#'
#' @param xSWH Data in wide format, as returned by \code{\link{x.24wide}}
#' @param cty Country, single or multiple countries (default: \code{"DE|IT"})
#' @param seas Season (single or multiple seasons (default: \code{c("SP","SU")})
#' @param yrRange Used in title of plot (default: \code{"2017-2019"})
#' @param subrange Plot only correlation between solar and wind (and not e.g. between hours of solar)? (default: FALSE)
#' @param upper.triangle Plot only upper triangle? (default: TRUE)
#' @param old.version Plot old version with other color codes, whereas new version  asfixed levels? (default: TRUE)
#' 
#' @examples
#' \dontrun{
#' contour.wide24(l.wide24, "DE", "SU", "2017-2019") 
#' }
#' 
contour.wide24 <- function (xSWH
                            ,cty = "DE|IT"
                            ,seas = c("SP","SU")
                            ,yrRange = "2017-2019"
                            ,upper.triangle = TRUE
                            ,subrange = FALSE
                            ,old.version = FALSE
                            ,plot2pdf = FALSE
                             ) {
    z = list()
    for (s in seas) 
        z[[s]] = xSWH[[s]][,grepl(cty,names(xSWH[[s]]))] 
    sw = do.call(rbind, z)
    par(pty="m") # type of plot region; "s" square, "m" maximal plotting region.
    contour.plot(correl = if (subrange) cor(sw)[1:24,25:48] else cor(sw)
                ,title.cont = paste(cty,s,yrRange,sep=", ")
                ,cty = cty
                ,s = seas
                ,triang = upper.triangle
                ,oldVer = old.version)
    if (plot2pdf) dev.print(pdf,
                            paste0("corHistContour_",
                                   cty, "_",
                                   paste(seas, collapse = ','), ".pdf")
                            )
}








### Scenarios by Clustering





#' Get the medoid of each cluster
#'
#' @param i Index
#' @param distmat Distance matrix between points
#' @param clusters Clusters
#' 
clust.medoid <- function(i, distmat, clusters) {
    ind = (clusters == i)
    #browser()
    names(which.min(rowSums( distmat[ind, ind, drop = FALSE] )))
}


#' Cluster by hierarchical method using medoids
#'
#' Single countries or multiple countries can be selected, whereas only a single season can be specified.
#' The agglomerative, hierarchical clustering uses a \code{method}, which is a
#' dissimilarity measure between clusters (so-called Linkage):
#' \describe{
#'    \item{\code{"ward.D2"} (default)}{Ward's method.
#'        Minimizes variance of (potentially) merged clusters
#'        (ward "likes" to produce clusters of equal sizes;
#'         "ward.D" was coded wrongly)}
#'  \item{\code{"average"}}{(= UPGMA):
#'             average distance between all points of two clusters}
#'  \item{\code{"complete"}}{Highest distance between point-pairs
#'            of two clusters}
#'  \item{\code{"single"}}{Smallest distance between
#'                 point-pairs of two clusters}
#' }
#' 
#' 
#' @param xSWH List over seasons of data in wide format (day hours are
#'     in columns). Given a season, a list element is a data.frame
#'     with columns "ATsolar01",
#'     "ATsolar02' etc., and rows = observations in the season; as returned by \code{\link{x.24wide}}
#' @param nclus Number of clusters
#' @param cty Country (default: \code{"DE"}; several countries \code{"AT|CH|DE|FR|IT"})
#' @param s Season (default: \code{"SU"})
#' @param interactive Plot detailed diagonstic of clustering (dendogramm etc.)? (default: FALSE)
#' @param crossCountry Cluster across countries (\code{cty} must have several regions)? (default: FALSE)
#' @param clustmethod Cluster method, that is, the linkage (default: Ward's method)
#' @param distmeasure Measure of distance between data points (default: \code{"euclidean"}).  Other possible distances are e.g. "maximum" and "manhattan"
#'
#' @return
#' For a single region, a data.frame (called \emph{scnBEM}) of scenarios
#' in stacked format is returned (i.e. the 24 day hours of scenario
#' values are in rows):
#' \tabular{cccrrr}{
#' \strong{scn} \tab \strong{country} \tab \strong{loadp} \tab
#'           \strong{prob} \tab \strong{windAvl} \tab \strong{solarAvl} \cr
#' "o1" \tab "AT" \tab "WI-D-01" \tab 0.01 \tab 0.3    \tab  0.0    \cr 
#' "o1" \tab "AT" \tab "WI-D-02" \tab 0.01 \tab 0.3    \tab  0.0    \cr
#' "o1" \tab "AT" \tab \ldots    \tab 0.01 \tab \ldots \tab  \ldots \cr
#' "o1" \tab "AT" \tab "WI-D-24" \tab 0.01 \tab 0.3    \tab  0.1    \cr
#' "o2" \tab "AT" \tab "WI-D-01" \tab 0.02 \tab 0.3    \tab  0.0
#' }
#' For multiple regions, a list over the countries is returned, where each list element is data.frame of format \emph{scnBEM}
#'
#' 
scn.clust <- function(xSWH
                     ,nclus = 20
                     ,cty = "DE"
                     ,s = "SU"
                     ,interactive = FALSE
                     ,crossCountry = FALSE
                     ,clustmethod = "ward.D2"
                     ,distmeasure = "euclidean")
{
    d = xSWH[[s]][,grepl(cty,names(xSWH[[s]]))] # select country and season 
    # distance measure between rows of d: "dist" returns lower-tringular matrix
    mydist = dist(d, method = distmeasure)
    hc = hclust(mydist, method=clustmethod) # hierarchical clustering
    # cutree returns vector of cluster membership in the order of original data rows
    mycl = cutree(hc, k=nclus) # split into k clusters
    print("Cluster sizes:"); print(table(mycl)) # count members in the clusters
    if (interactive) {
        print(mycl) # vector of cluster membership in order of original data rows
        plot(hc, cex=0.3) # plot the dendrogram
        print(mycl[hc$order]) # vector of cluster membership ordered as in dendogram
        ## grab cluster i:
        #i = 2  
        #cluster = d[mycl == i,]
        dID = cbind(d, clusterID = mycl) # add a cluster ID to original data 
        print(dID[hc$order,]) # data ordered as in dendogram
    }
    ## Find the medoids:
    distMatrix = as.matrix(mydist) # convert lower-triang to full matrix
    # Get index of medoids:
    idx.medoid = sapply(unique(mycl), clust.medoid, distMatrix, mycl)
    print("Index of medoids:"); print(idx.medoid)
    sw.medoid = d[idx.medoid,] # get the medoids
    #browser()
    # the probability of the medoid is proportional to cluster size:
    p.medoid = as.numeric(table(mycl)/nrow(d))
    print(sum(p.medoid)) # check: probabilities should sum to 1
    ## make a data-frame with the medoid as scenarios
    n = nclus*24
    scnBEM = data.frame(
        scn = character(n),
        country = character(n),
        loadp = character(n),
        prob = numeric(n),
        windAvl = numeric(n),
        solarAvl = numeric(n),
        stringsAsFactors = FALSE)

    if (crossCountry) {
        ctylist = list()
        for (ctry in ctyNames) {
            scnBEM$country = ctry
            scnBEM$scn =  paste0("o",rep(1:nclus, each=24))
            scnBEM$loadp = rep(load.24(s), nclus)
            for (i in 1:nclus) {
                scnBEM$prob[(24*(i-1)+1):(24*i)] = p.medoid[i]
                scnBEM$solarAvl[(24*(i-1)+1):(24*i)] =
                    as.numeric(sw.medoid[i,grepl(paste0("solar",ctry),names(sw.medoid))])
                scnBEM$windAvl[(24*(i-1)+1):(24*i)] =
                    as.numeric(sw.medoid[i,grepl(paste0("wind",ctry),names(sw.medoid))])   
            }
            ctylist[[ctry]] = scnBEM
        }
        return(do.call(rbind, ctylist))
    } else {
        scnBEM$country = cty
        scnBEM$scn =  paste0("o",rep(1:nclus, each=24))
        scnBEM$loadp = rep(load.24(s), nclus)
        for (i in 1:nclus) {
            scnBEM$prob[(24*(i-1)+1):(24*i)] = p.medoid[i]
            scnBEM$solarAvl[(24*(i-1)+1):(24*i)] = as.numeric(sw.medoid[i,1:24])
            scnBEM$windAvl[(24*(i-1)+1):(24*i)] = as.numeric(sw.medoid[i,25:48])   
        }
        return(scnBEM)
    }
}



#windows()
#contour.plot(cov.wt(sw.medoid, p.medoid, cor=T)$cor,
#             paste0(cty,", ",s,", 2017-2019, ",nclus," clusters"),cty,m.r=1,s)
#dev.print(pdf,"corScenContour_DE_SU(Hclust).pdf")



#' Clustering of cross-country scenarios over all seasons
#'
#' Stack the cross-country scenarios over all the seasons. Hence, increase the scenario index, because different seasons have independent scenarios in the electrictiy market model BEM
#' 
#' @param xSWH List over seasons of data in wide format (day hours are in columns).Given a season,  a list element is a data.frame with columns "ATsolar01", "ATsolar02' etc., and rows = observations in the season; as returned by \code{\link{x.24wide}}
#' @param nclust Number of clusters
#' 
#' @return scenarios over all countries over all seasons
#'
scn.clust.cross <- function(xSWH
                           ,nclust = 20)
{
    xClus = list() #cross-cluster, across countries
    scnNumOffset = 0
    for (s in seasNames) {
        print(s)
        xClus[[s]] = scn.clust(xSWH, nclus = nclust, cty = "AT|CH|DE|FR|IT", s=s, crossCountry = T)
        # scenarios in other seasons have increased scenario number
        xClus[[s]]$scn = paste0("o", as.numeric(substring(xClus[[s]]$scn, 2)) + scnNumOffset)
        scnNumOffset = scnNumOffset + nclust
    }
    clustAllCty = do.call(rbind, xClus)
    return(clustAllCty)
}



#' Average wind and solar availability to the 96 loadperiods of BEM per country
#'
#' A convenience function for BEM. Not required for PCA or clustering.
#'
#' @param xSWH Data in wide format, as returned by \code{\link{x.24wide}}
#' 
#' @return Data.frame of solar and wind availability for each load period and for each country
#'  \itemize{
#'     \item rows: 5 * 96 (countries x seasons x hours), e.g., "AT.WI-D-01"
#'     \item columns: "solar", "wind"
#' }
#' 
agggregate.avl.BEM <- function(xSWH)
{
    nrows = length(ctyNames)*length(seasNames)*24
    SWloadperiods = matrix(ncol = 2, nrow = nrows,
                           dimnames=list(seq(nrows),
                                         c("solar", "wind")))
    for (tech in c("solar", "wind")) {
        v = numeric(0)
        # stack-up the mean over the countries and seasons in vector "v"
        for (cty in ctyNames)
            for (s in seasNames) {
                avgEachHour = colMeans(
                    xSWH[[s]][,grepl(cty,names(xSWH[[s]])) &
                               grepl(tech,names(xSWH[[s]]))]) 
                names(avgEachHour) = paste0(cty, ".", load.24(s))
                v = c(v,avgEachHour)
            }
        SWloadperiods[,tech] = v
    }
    rownames(SWloadperiods) = names(v)
    return(SWloadperiods)
}






### PCA Single







#' Fit a Weibull distribution from data
#'
#' Before the fitting, the input data is shifted such that most of the values becomes positive, which may alos require an optional sign-reversion (if data was mainly nonnegative)
#' 
#' @param data Vector of data
#' 
#' @return Vector with components:
#' \enumerate{
#'  \item Shape and scale parameters of fitted Weibull distribution
#'  \item Sign reversion and shift (applied before fitting)
#'} 
fitWeibull <- function (data) {
    signF1 = if(sum(data)<0) -1 else 1
    fac1.mod = signF1*data
    # for Weibull, we have to shift to positive distribution, and cut>0
    fac1Shift = min(quantile(fac1.mod,0.01), 0)
    fac1.mod = fac1.mod - fac1Shift
    fac1.mod = fac1.mod[fac1.mod>0.000001]
    fit.dist1 = fitdist(fac1.mod, "weibull") #  lower = c(0, 0))
    # plot(fit.dist1)
    #shW=coef(fit.dist1)["shape"]
    #scW=coef(fit.dist1)["scale"]
    c(coef(fit.dist1), sign=signF1, shift=fac1Shift)
}


#' Get values and probabilities of a discretized Weibull distribution
#'
#' Intervals are equally spaced from quantiles: \code{0.1/i} up to \code{1-0.1/i}.
#' Values are the the mid-points of the intervals.
#' 
#' @param i Number of values 
#' @param  pWeib Parameters from fitting returned by \code{fitWeibull}
#'
#' @return: Matrix with two columns: probability, value
#' 
discrWeibull <- function (i,pWeib) {
    # equaly spaced intervals of values:
    v = seq(from=qweibull(0.1/i,pWeib[1],pWeib[2]),to=qweibull(1-0.1/i,pWeib[1],pWeib[2]), length.out=i+1)
    p = diff(pweibull(v,pWeib[1],pWeib[2]))
    v = (v[-1] + v[-length(v)])/2 # middle points of intervals
    # apply scaling
    v = pWeib[3]*(v + pWeib[4])
    cbind(p/sum(p),v)
}





#' PCA for a single season and country
#'
#' There can be several seasons (together) and countries (additional dimension)
#'
#' Log-normal factors gave usually bad results. Avoid!
#'
#' @param xSWH List over seasons of data in wide format (day hours are in columns).Given a season, a list element is a data.frame with columns "ATsolar01", "ATsolar02' etc., and rows = observations in the season; as returned by \code{\link{x.24wide}}
#' @param cty Country
#' @param season Season
#' @param m Number of loadings to plot
#' @param deMean Should data first be de-meaned? (default: FALSE)
#' @param plot2pdf Should plots be copied to pdf? (default: FALSE)
#' @param load If load is included make a different pdf file name (default: FALSE)
#'
#' @return list: 1. Factors implied by PCA over the input data points (rows), that is, a daily time series of factors (class: matrix); 2. pcaloads
#' 
#' @examples
#' \dontrun{
#' PCAsingle("DE|IT","SU")
#' PCAsingle("AT","SU")
#' PCAsingle("DE","WI|SP|SU|FA")
#' }
PCAfac <- function(xSWH
                  ,cty = "DE"
                  ,season = "SP|SU"
                  ,m = 4
                  ,deMean = FALSE
                  ,plot2pdf = FALSE
                  ,load = FALSE
                   )
{
    seasVec = strsplit(season, "\\|")[[1]]    
    z = list()
    for(s in seasVec)
        z[[s]] =  xSWH[[s]][,grepl(cty,names(xSWH[[s]]))]
    sw = do.call(rbind, z)
    ## Start PCA:
    ## de-mean time series: (USUALLY BAD RESULTS)
    m.sw = colMeans(sw)
    if (deMean) sw = scale(sw, scale = FALSE) else m.sw = 0

    ## log time series: (USUALLY BAD RESULTS)
    # prc = princomp(log(xSWH+0.00001))

    ## logit time series: (USUALLY BAD RESULTS)
    # k = 1/20000
    # x0 = 0
    # g.inv.logit = function(x,L,k,x0) L/(1+exp(-k*(x-x0)))  
    # g.logit     = function(x,L,k,x0) x0 + 1/k*log(x/(L-x))
    # logit.sw = g.logit(sw+1, max(sw)+2,k,0)
    # inv.sw   = g.inv.logit(logit.sw, max(sw)+2,k,0)
    # prc = princomp(logit.sw)

    prc = princomp(sw)
    print(summary(prc, loadings = TRUE, digits = 2))
    windows()
    screeplot(prc, main="Variance of Principal Components", col="blue")

    if(plot2pdf) {
        fn = paste0("screeplot_",
                    gsub("\\|","_",cty), "_",
                    gsub("\\|","_",season),
                    if (load) "_with_load", ".pdf")
        dev.print(pdf,fn)
    }

    pcaloads = loadings(prc)
    windows()
    layout(matrix(seq(1,m), nrow = 1, ncol = m, byrow = TRUE),
           widths = c(1.3,rep(1,m-1)))
    for (i in 1:m) {
        par(mar=c(5,ifelse(i==1,6,1),4,1)+.1)
        barplot(
            pcaloads[,i], 
            names = substr(rownames(pcaloads),1,10), 
            horiz = TRUE,
            col = "blue",
            yaxt = ifelse(i==1,"s","n"), # n = normal, s = supress
            xlab = c("First", "Second", "Third", "Fourth")[i],
            las = 2, # labels are perpendicular
            cex.names = 1.2,
            cex.lab = 1.2
        )
    }
    
    if(plot2pdf) {
        fn = paste0("PCA_",
                    gsub("\\|","_",cty), "_",
                    gsub("\\|","_",season),
                    if (load) "_with_load", ".pdf")
        dev.print(pdf,fn)
    }

    fac = as.matrix(sw) %*% pcaloads
    l = list()
    l[["fac"]] = fac
    l[["pcaloads"]] = pcaloads
    l[["m.sw"]] = m.sw
    return(l)
}




#' Factor-Analysis
#'
#' Factor-Analysis of first three PC
#'
#' @param fac Daily time series of factors as returned by \code{\link{PCAfac}} 
#' 
#' @return List \emph{fitParam} with elements:
#' \describe{
#'  \item{"mean"}{c(mF1, mF2, mF3)}
#'  \item{"std"}{c(sF1, sF2, sF3)}
#'  \item{"weibull"}{parWeibull (as returned by \code{\link{fitWeibull}}}
#' }
#' 
facAnalysis <- function(fac
                       ,logFac1 = FALSE
                       ,logFac2 = FALSE
                       ,logFac3 = FALSE
                       ,plot2pdf = FALSE)
{  
    fac1 = fac[,1]
    fac2 = fac[,2]
    fac3 = fac[,3]

    # if log-factors: Take log and shift to positive values:
    min1 = max(-fac1,0) + 0.0001
    min2 = max(-fac2,0) + 0.0001
    min3 = max(-fac3,0) + 0.0001
    if (logFac1) fac1 = log(fac1 + min1) 
    if (logFac2) fac2 = log(fac2 + min2) 
    if (logFac3) fac3 = log(fac3 + min3)

    # Mean and and sd of factors
    mF1 = mean(fac1)
    sF1 = sd(fac1)
    mF2 = mean(fac2)
    sF2 = sd(fac2)
    mF3 = mean(fac3)
    sF3 = sd(fac3)
    print(c(mF1, mF2, mF3, sF1, sF2, sF3))

    # Alternative is to fit the distribution directly (NOT USED):
    #windows()
    #descdist(fac1)
    #plotdist(fac1)
    #descdist(fac2)
    #plotdist(fac2)
    #descdist(fac3)
    #plotdist(fac3)

    ## Fit first factor with Weibull distribution
    
    parWeibull = fitWeibull(fac1) 

    fac1.mod.p = -fac1 - min(quantile(-fac1,0.01), 0)
    fac1.mod.p = fac1.mod.p[fac1.mod.p>0.000001]
    fit.dist1 = fitdist(fac1.mod.p, "weibull") #  lower = c(0, 0))
    #plot(fit.dist1)
    #fit.dist1 = fitdist(fac1, "norm")
    windows()
    qqcomp(fit.dist1,  main = "1stPC")
    if (plot2pdf) dev.print(pdf, "qqcomp_1stPC.pdf")
    windows()
    denscomp(fit.dist1, xlab="1st PC")
    if (plot2pdf) dev.print(pdf, "denscomp_1stPC.pdf")

    d = discrWeibull(10, parWeibull)
    windows()
    plot(d)
    windows()
    plot(d[,2],d[,1])

    ## Fit second and third factor with normal distribution
    
    fit.dist2 = fitdist(fac2, "norm")
    fit.dist3 = fitdist(fac3, "norm")
    #plot(fit.dist1)
    #plot(fit.dist2)
    #plot(fit.dist3)
     windows()
    qqcomp(fit.dist2, main = "2ndPC")
    if (plot2pdf) dev.print(pdf, "qqcomp_2ndPC.pdf")
    windows()
    denscomp(fit.dist2, xlab="2nd PC")
    if (plot2pdf) dev.print(pdf, "denscomp_2ndPC.pdf")
    windows()
    qqcomp(fit.dist3,  main = "3rdPC")
    if (plot2pdf) dev.print(pdf, "qqcomp_3rdPC.pdf")
     windows()
    denscomp(fit.dist3, xlab="3rd PC")
    if (plot2pdf) dev.print(pdf, "denscomp_3rdPC.pdf")

    fitParam = list()
    fitParam[["mean"]] = c(mF1, mF2, mF3)
    fitParam[["std"]] = c(sF1, sF2, sF3)
    fitParam[["weibull"]] = parWeibull
    
    return(fitParam)
}





### Scenario generation of PCA for selected country and season 




#' Discretized standard binomial distribution
#'
#' @param j Number of disretization (yields j+1 values) 
#' 
#' @return
#' \code{binomSd} returns values (j=0 yields 0);
#' \code{binomSdProb} returns probabilities
#'
#' @examples
#' bimonSd(2)
#' binomSdProb(2)
#' 
binomSd <- function(j) if(j>0) 2/sqrt(j)*(0:j - 0.5*j) else 0
#' @rdname binomSd
binomSdProb <- function(j) choose(j,0:j)*0.5^j



#' Generate scenarios given a PCA factor model
#'
#' @param J1 J1+1 = number of scenarios for factor 1
#' @param J2 J2+1 = number of scenarios for factor 2
#' @param J3 J3+1 = number of scenarios for factor 3
#' @param threshold Bound to identify a pathological, negative avl. scenario
#' @param dropNeg Should such negative scenario be dropped? (default: TRUE)
#' @param pcaloads Loadings as returned by \code{\link{PCAfac}}
#' @param m.sw Mean values as returned by \code{\link{PCAfac}}
#' 
#' @return data-frame "scenData"
#' \itemize{
#'   \item row 1: probability of scenario
#'   \item row 2,...,48+1: scneario realization of wind and solar
#'   \item columns: "scn1", "scn2", ... "scn(smax)"
#' }
scn.PCAfac <- function(fitParam
                      ,pcaloads
                      ,m.sw
                      ,J1 = 5  # J1+1 scenarios for factor 1
                      ,J2 = 1  # J2+1 scenarios for factor 2
                      ,J3 = 1  # J3+1 scenarios for factor 3
                      ,threshold = -0.01 #threshold = -10
                      ,dropNeg = TRUE
                      ,logFac1 = FALSE
                      ,logFac2 = FALSE
                      ,logFac3 = FALSE
                       )
{
    parWeibull = fitParam[["weibull"]]
    mF1 = fitParam[["mean"]][1]
    mF2 = fitParam[["mean"]][2]
    mF3 = fitParam[["mean"]][3]
    sF1 = fitParam[["std"]][1]
    sF2 = fitParam[["std"]][2]
    sF3 = fitParam[["std"]][3] 
    smax = (J1+1)*(J2+1)*(J3+1) 
    binom1 = binomSd(J1)
    binom2 = binomSd(J2)
    binom3 = binomSd(J3)
    # p1 = binomSdProb(J1)
    p1 = discrWeibull(J1+1,parWeibull)[,1]
    p2 = binomSdProb(J2)
    p3 = binomSdProb(J3)
    # plot(binom1, p1)
    plot(binom2, p2)
    plot(binom3, p3)

    # Discretized values factor distribution:
    # scnF1 = mF1 + sF1*binom1 # normal 
    scnF1 = discrWeibull(J1+1,parWeibull)[,2] # Weibull
    scnF2 = mF2 + sF2*binom2 # normal
    scnF3 = mF3 + sF3*binom3 # normal
    windows()
    plot(scnF1,p1)

    # Discretized values of shifted lognormal distribution:
    if(logFac1) scnF1 = exp(mF1 + sF1*binom1) - min1
    if(logFac2) scnF2 = exp(mF2 + sF2*binom2) - min2
    if(logFac3) scnF3 = exp(mF3 + sF3*binom3) - min3

    ## Init scnData (based on PCA) using binomial distribution
    n = 48;  # number of observations per single country (wind+solar)
    scnData = data.frame(matrix(ncol = smax, nrow = 1 + n ))
    colnames(scnData) = paste0("scn", 1:smax)
    scnCnt = 1
    for(i1 in 1:(J1+1)) {  # first factor
        for(i2 in 1:(J2+1)) { # second factor
            for(i3 in 1:(J3+1)) { # third factor 
                v = pcaloads[,1:3] %*% c(scnF1[i1], scnF2[i2], scnF3[i3])
                scnData[1,paste0("scn", scnCnt)] = p1[i1]*p2[i2]*p3[i3]
                scnData[2:(length(v)+1), paste0("scn", scnCnt)] =
                    v[1:length(v)]
                scnCnt = scnCnt + 1
            }
        }
    }

    # Shrink scenarios with values outside of [0,1]:
    scnData[-1,] = m.sw + scnData[-1,]
    windows()
    v = as.numeric(as.matrix(scnData[,-1]))
    hist(v, plot=TRUE)
    length(as.matrix(scnData[,-1][scnData>1]))/length(v)
    length(as.matrix(scnData[,-1][scnData< threshold]))/length(v)
    #scnData[scnData>1] = 1
    #scnData[scnData<0] = 0
    # Drop scenarios with negative values:
    if (dropNeg) {
        drop = names(scnData)[which(apply(scnData,2,min) < threshold)]
        scnData = scnData[!(colnames(scnData) %in% drop)]
        scnData[1,] = scnData[1,] / sum(scnData[1,]) # update probabilities to sum again to 1
    }
 return(scnData)
}



#' Plot scenarios generated from a PCA factor model
#'
#' The probability of a scenarios is proportional to the line width.
#' 
#' For the used rainbow palette you can also select start/end color
#' (red = 0, yellow = 1/6, green = 2/6, cyan = 3/6, blue= 4/6 and magenta = 5/6)
#' and saturation (s) and value (v):
#' rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
#'
#' @param scnData Scenarios as returnd by \code{\link{scn.PCAfac}}
#' 
plotScnData <- function(scnData, cty, season, plot2pdf = FALSE) {
    windows()
    n = 5
    layout(matrix(1:2, ncol = 2, byrow = TRUE),
           widths = c(1.2,1), heights = c(0.6,1))
    par(oma=c(0,0,3,0))
    plot.ts((scnData[-1,])[1:24,],
            plot.type = "single", 
            xlab="hour",
            ylab="Availability",
            col=rainbow(n, s = 1, v = 1, start = 0, end = 5/6),
            lty=1,
            lwd=scnData[1,]*70,
            main = "solar")
    par(mar=c(4,1,2,1)+0.1)
    plot.ts((scnData[-1,])[25:48,],
            plot.type = "single",
            xlab="hour",
            ylab="Availability",
            col=rainbow(n),
            lwd=scnData[1,]*70,
            main = "wind")
    title(paste(cty,season), outer=TRUE)
    par(mar=c(5,4,4,2)+0.1) # restore default

    if (plot2pdf) dev.print(pdf, "scenarios_632_DE_SU.pdf")
}




### Scenarios of PCA factor model with loop over all four seasons 


#' Create scenarios for all countries and all seasons
#'
#' Combination of previous functions: \code{\link{scn.PCAfac}} etc. Sensitivity analysis of a single country w.r.t different J is also possible
#' 
#' @param J.cases List of vectors of factor realizations,
#'     e.g. \code{list(c(5,2,1))}; usually the list has just a
#'     element. Sensitvity analysis is also possible with more list elements
#' @param comparePCA Sensitivity analysis? (default: FALSE)
#' @param compareCty Country of sensitivity analysis (default: "DE")
#' @param logFac Switches for log-normal factors (default: FALSE)
#' @param weibFac1 The first factor is fitted to weibull? (default: TRUE)
#' @param thresMin Bound to identify a pathological, negative avl. scenario
#' @param thresMinDeMean Different bound if data was de-meaned
#'
#' @return List of scenarios over countries of type \emph{scnBEM} (see help \code{\link{scn.clust}}). If sensitivity analysis: List of single-country scenarios over the different elements of \code{J.cases}
#' 
#' @examples
#' \dontrun{
#' scn.All(J.cases = list(c(5,2,1)))
#' scn.All(J.cases = list(c(5,0), c(5,2,0), c(5,2,1)))
#' }
#' 
scn.All <- function(xSWH,
                    J.cases =  list(c(5,2,1))
                   ,thresMin = -0.1
                   ,thresMinDeMean = -0.01
                   ,weibFac1 = TRUE 
                   ,comparePCA = FALSE
                   ,compareCty = "DE"
                   ,logFac = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE) 
                    )
{     
    loadh = 24 #  loadperiods (hours) per scenario (in fact, is always = 24)
    # collect scenarios of countries in a list
    scn.ws = list()
    # sensitivity analysis:
    # collect scenarios of a single country (compareCty)
    # with different J in a list, and count them
    scn.ws.comp = list(); icomp = 1

    ## Init parameters: 

    for (J in J.cases) {
        print(paste("case:",paste(J,collapse=",")))
        for (cty in ctyNames) {
            print(paste("country:",cty))    
            Jseq = seq_along(J) # sequence:1,2,3,...,number of factor
            deMean = switch (cty, AT = TRUE, FALSE) 
            if (deMean) thresMin = thresMinDeMean
            binom = lapply(J, binomSd ) # binomial distr values
            p = lapply(J, binomSdProb)# binomial distr probabilities
            smax = prod(J+1)  # number of scenarios per season
            n = 4*smax*loadh; # scenario values (if all would be positive)

            ## Fill data-frame scnBEM for all seasons and a specific country:
            scnBEM = data.frame(
                scn = character(n),
                country = character(n),
                loadp = character(n),
                prob = numeric(n),
                windAvl = numeric(n),
                solarAvl = numeric(n),
                stringsAsFactors = FALSE)           
            rowN   = 1 # row number in df scnBEM
            scnCnt = 1 # scenario counter
            for (seas in seasNames) {
                sw = xSWH[[seas]][,grepl(cty,names(xSWH[[seas]]))]
                m.sw = colMeans(sw)
                if (deMean) sw = scale(sw, scale = FALSE) else m.sw = 0
                prc = princomp(sw)
                pcaloads = loadings(prc)
                fac = as.matrix(sw) %*% pcaloads
                # fac-columns = factors 1,2,3,
                # Log-factors: shift to positive values:
                adjLog = as.numeric(length(J))
                for (i in Jseq) if (logFac[i]) {
                                    adjLog[i] = (-min(fac[,i],0) + 0.00001)
                                    fac[,i]=log(fac[,i]+adjLog[i])
                                }
                mF = colMeans(fac) # vector of means across factors
                sF = apply(fac,2,sd) # vector of sd across factors
                # first factor is fitted to Weibull 
                if (weibFac1) { 
                    parWeibull = fitWeibull(fac[,1])
                    p[[1]] = discrWeibull(J[1]+1,parWeibull)[,1]
                }
                # Disrecized values of normal distributions (list across factors)
                # (log-factors: values of shifted lognormal distribution):
                scnF = lapply(Jseq, function(i)
                    if(logFac[i]) {
                        exp(mF[i] + sF[i]*binom[[i]]) - adjLog[i]
                    } else { if (i==1 & weibFac1) discrWeibull(J[1]+1,parWeibull)[,2]
                             else  mF[i] + sF[i]*binom[[i]]
                    })
                df = expand.grid(lapply(J+1,seq)) # rows with all combination of factors
                # for each row of df (combination of factor discretiztaions) get wind+solar:
                # columns of v contain vector: c(prob,scenario values)
                v = apply(df, 1, function(x)
                    c( prod(sapply(Jseq, function(j) p[[j]][x[j]])),
                      pcaloads[,Jseq] %*% sapply(Jseq,
                                                 function(j) scnF[[j]][x[j]])
                      ))
                # write only a new scenario if above threshold
                ncolOld = ncol(v)
                print(paste("Season =", seas,
                            ": lowest value =", min(v[-1,]+m.sw)))
                v =v[,apply(v,2, function(x) min(x[-1]+m.sw) >= thresMin)]
                print(paste(ncolOld - ncol(v),
                            "scenarios removed (out of", ncolOld,")"))
                ## fill the next 24 rows in scnBEM:
                apply(v,2, function(x) {
                    rng = rowN:(rowN + loadh - 1)
                    scnBEM$scn[rng] <<- paste0("o", scnCnt)
                    scnBEM$country[rng] <<- cty
                    scnBEM$prob[rng] <<- x[1]
                    scnBEM$loadp[rng] <<- load.24(seas)
                    # make scn in [0,1], add mean, just savety measure.
                    x2 = pmax(x[-1]+m.sw, 0) # parallel maximum (>=0)
                    x2 = pmin(x2, 1)         # parallel minimum (<=1)
                    scnBEM$solarAvl[rng] <<- x2[1:loadh]
                    scnBEM$windAvl[rng]  <<- x2[(loadh+1):(2*loadh)]
                    rowN <<- rowN + loadh
                    scnCnt <<- scnCnt + 1 
                }) 
            }
            
            scnCnt = scnCnt-1
            scnCnt
            dim(scnBEM)
            # not all rows may be filled because of cutoff
            nrow(scnBEM[scnBEM$scn == "",])
            nrow(scnBEM[scnBEM$scn != "",])
            scnBEM = scnBEM[scnBEM$scn != "",]  
            # increase probabilities again to sum to 1:
            for (s in seasNames) {
                selScn = grepl(s,scnBEM$loadp)
                scnBEM$prob[selScn] = scnBEM$prob[selScn] /
                    (sum(scnBEM$prob[selScn]/24))
            }

            print(dim(scnBEM))
            scn.ws[[cty]] = scnBEM            
            if (comparePCA & cty==compareCty) {
                scn.ws.comp[[icomp]] = scnBEM
                icomp = icomp + 1
            }
            
        }
    }
    if(comparePCA) return(scn.ws.comp) else return(scn.ws)
}



    
#' Write scenarios over seasons and countries into files
#'
#'@param scn.ws List of scenarios over countries of type \emph{scnBEM} (see help \code{\link{scn.clust}}), as returned by \code{\link{scn.All}}
#' 
writeScnAll <- function(scn.ws) {
    lapply(scn.ws,
           function (x)
               write.csv(x,
                         paste0("scnWindSolar", x[1,"country"], ".csv"),
                         row.names = FALSE))
}




### Plot scenarios of PCA factor model



#' Plot scenarios of PCA factor model
#'
#' Plot scenarios for selected country and season.
#'
#'@param scnAll List of scenarios over countries of type \emph{scnBEM} (see help \code{\link{scn.clust}}), as returned by \code{\link{scn.All}}
#' 
#' @param cty Country
#' @param s Season
#' 
plotScn <- function(scnAll
                   ,cty = "DE"
                   ,s = "SU")
{
    scnBEM = scnAll[[cty]]
    y = scnBEM[scnBEM$country==cty & grepl(s,scnBEM$loadp),
               c("scn","solarAvl","windAvl")]
    # data is in stacked form by scn: each scenario should be in a column
    # (mean was added earlier)
    y = rbind(unstack(y, solarAvl ~ scn), unstack(y, windAvl ~ scn))
    # stopifnot(length(y)==prod(J+1))
    windows()
    layout(matrix(1:2, ncol = 2, byrow = TRUE),
           widths = c(1.3,1))
    par(oma=c(0,0,3,0))
    plot.ts(y[1:24,],
            plot.type = "single", 
            xlab="hour",
            ylab="Availability",
            col=rainbow(length(y)), 
            main = "solar")
    par(mar=c(4,1,2,1)+0.1)
    plot.ts(y[25:48,],
            plot.type = "single",
            xlab="hour",
            ylab="Availability",
            col=rainbow(length(y)),
            main = "wind")
    title(paste(cty,s), outer=TRUE)
}



### Correlation: Contourplot



#' Plot correlations of empirical data
#'
#' @param xSWH List of scenarios over countries of type \emph{scnBEM} (see help \code{\link{scn.clust}}), as returned by \code{\link{scn.All}}
#' @param cty Country
#' @param s Season
#' 
plotCor.Emp <- function(xSWH
                       ,cty = "DE"
                       ,s = "SU"
                       ,oldV = FALSE
                       ,plot2pdf = FALSE) {
    sw = xSWH[[s]][,grepl(cty,names(xSWH[[s]]))]
    windows()
    contour.plot(cor(sw),
                 title.cont = paste(cty,s,"2017-2019",sep=", "),
                 cty, s, triang = TRUE, oldVer=oldV)
    if (plot2pdf) dev.print(pdf,
                            paste0("corContourHist_",cty,"_", s,".pdf"))
    return(cov(sw))
}



#' Clustering: Generate scenarios and plot correlations
#'
#' @param  xSWH List over seasons of data in wide format (day hours are in columns).Given a season,  a list element is a data.frame with columns "ATsolar01", "ATsolar02' etc., and rows = observations in the season.
#' @param cty Country
#' @param s Season
#'
#' @return covariance-matrix
#'
plotCor.genClust <- function(xSWH
                            ,cty = "DE"
                            ,s = "SU"
                            ,n = 20
                            ,oldV = FALSE
                            ,plot2pdf = FALSE)
{
    scnBEM = scn.clust(xSWH, nclus=20, cty, s)
    y = scnBEM[scnBEM$country==cty & grepl(s,scnBEM$loadp),
               c("scn","solarAvl","windAvl", "prob")]
    y2 = rbind(unstack(y, solarAvl ~ scn), unstack(y, windAvl ~ scn))
    w2 = unstack(y, prob ~ scn)
    weighted.corr = cov.wt(t(y2), wt = as.numeric(w2[1,]), cor = TRUE)
    windows()
    contour.plot(weighted.corr$cor,
                 title.cont = paste0(cty, ", ", s, ", ", n, " clusters"),
                 cty, s, triang = TRUE, oldVer=oldV)
    if (plot2pdf) dev.print(pdf,
                            paste0("corContourClust_",cty,"_",s,".pdf"))
    return(weighted.corr$cov)
}


#' PCA factor model: Generate scenarios and plot correlations
#'
#' @param  xSWH List over seasons of data in wide format (day hours are in columns).Given a season,  a list element is a data.frame with columns "ATsolar01", "ATsolar02' etc., and rows = observations in the season.
#' @param cty Country
#' @param s Season
#' @param textArg Optional argument to insert in title of plot
#' @param J J+1 = Number of realizations of factors
#'
#' @return covariance-matrix
#' 
plotCor.genPCA <- function(xSWH
                         ,cty = "DE"
                         ,s = "SU"
                         ,J = c(5,2,1)
                         ,oldV = FALSE
                         ,plot2pdf = FALSE
                         ,textArg)
{ 
    scn.ws = scn.All(xSWH, J.cases = list(J))
    scnBEM = scn.ws[[cty]]
    y = scnBEM[scnBEM$country==cty & grepl(s,scnBEM$loadp),
               c("scn","solarAvl","windAvl", "prob")]
    # data is in stacked form by scn: each scenario should be in a column
    # (mean was added earlier)
    y2 = rbind(unstack(y, solarAvl ~ scn), unstack(y, windAvl ~ scn))
    w2 = unstack(y, prob ~ scn)
    weighted.corr = cov.wt(t(y2), wt = as.numeric(w2[1,]), cor = TRUE)
    windows()
    if (missing(textArg))
        textArg = paste0("PCA=(",paste(J+1,collapse=","),")")
    contour.plot(weighted.corr$cor,
                 title.cont =  paste0(cty, ", ", s, ", ", textArg),
                 cty, s, oldVer = oldV)
    if (plot2pdf)
        dev.print(pdf, paste0("corContourPCA_",cty,"_",s,".pdf"))
    return(weighted.corr$cov)
}







#' Check with empirical distribution function
#'
#' Plot empirical distribution function of wind and solar.
#' Select a season or entire year.
#'
#' @param scn.ws.or.comp List of scenarios over countries of type \emph{scnBEM} (see help \code{\link{scn.clust}}). If sensitivity analysis (comparison): List of single-country scenarios over the different elements of \code{J.cases}. As returned by \code{\link{scn.All}}
#' @param cty Country
#' @param s Season (default is fully year: s=seasNames)
#' @param yrs Used for plot title, range of years
#' @param comparePCA Is it a sensitivity analysis (default: FALSE)
#' @param J.cases Parameters of sensitivity analysis, used in title of plot
#' 
#' @examples
#' \dontrun{
#' plotEmpirical(scn.ws) # full year
#' plotEmpirical(scn.ws, s="SU") # SU
#' }
#' 
plotECDF <- function(scn.ws.or.comp
                    ,cty = "DE"
                    ,s = seasNames
                    ,yrs = "2017-2019"
                    ,comparePCA = FALSE
                    ,plot2pdf = FALSE
                    ,J.cases 
                     )
{
    if (comparePCA) {
        strCompPCA = paste("Discretization",
                           gsub("c","",
                                as.character(lapply(J.cases,function(x) x+1))))
    } else scnBEM = scn.ws.or.comp[[cty]]

    sBar = paste(s, collapse="|")
    for (tech in c("wind", "solar")) {
        windows()
        par(ps=9)
        v = as.numeric(x[(month(x) %in% seasons[,s]),paste0(cty,tech)])
        plot(ecdf(v),
             xlab=paste(tech, "capacity factor"),
             ylab="Probability",
             lwd=3,
             cex.axis = 2,
             cex.lab = 2,
             cex.main = 2,
             main =  paste(cty, sBar) 
             )
        
        if (comparePCA) {
            colCmp = 2 # colors: 1=black, 2=red, 3=green, 4=blue, 5=magenta
            for (z in scn.ws.or.comp) {
                v = z[c("prob",paste0(tech,"Avl"))][grepl(s,z$loadp),]
                lines(ewcdf(v[,2],weights = v$prob/96),
                      col=colCmp,
                      lwd=1
                      );
                colCmp = colCmp + 1
            }
        } else {
            v = scnBEM[c("prob",paste0(tech,"Avl"))][grepl(sBar,scnBEM$loadp),]
            lines(ewcdf(v[,2], weights = v$prob/96),
                  col="blue",
                  lwd=3
                  )
        }
        if (comparePCA) {
            legend("bottomright",
                   c(paste("Historical",yrs),
                     strCompPCA),
                   col=1:(colCmp-1),
                   lwd=3, inset=0.1, cex=2
                   )
            if (plot2pdf)
                dev.print(pdf,
                          paste0("compcdf_J_", tech, "_", cty,"_", s,".pdf"))
        }
        else {
            legend("bottomright",
                   c(paste("Historical",yrs), "Principal Comp."),
                   col=c("black", "blue"),
                   lwd=3, inset=0.1, cex=2
                   )
            if (plot2pdf)
                dev.print(pdf,
                          paste0("compcdf_", tech, "_", cty, "_", sBar, ".pdf"))
        }
    }   
}




### Comparison of availability factors




#' Seasonal availability factors implied by PCA
#'
#' @param scnBEM Data-frame  
#' @return avlPCA Array
#' 
avl.PCA <- function(scnBEM)
{
    avlPCA = array(dim = c(4,2,2),
                   dimnames = list(seasNames,
                                   c("solar", "wind"),
                                   c("mean", "sd"))
                   )
    for (seas in seasNames) {
        r = grepl(seas,scnBEM$loadp)
        p = scnBEM$prob[r]/24 
        avlPCA[seas, "solar", "mean"] = weighted.mean(scnBEM$solarAvl[r], p)
        avlPCA[seas, "wind", "mean"] = weighted.mean(scnBEM$windAvl[r], p)
        avlPCA[seas, "solar", "sd"] =
            sqrt(sum(p * (scnBEM$solarAvl[r] -avlPCA[seas,"solar","mean"])^2))
        avlPCA[seas, "wind", "sd"] =
            sqrt(sum(p * (scnBEM$windAvl[r] - avlPCA[seas,"wind","mean"])^2))
    }
    return(avlPCA)
}





#'  Seasonal availability factors from empirical data
#'
#' @param x time-series of data (xts)
#'
avl.emp <- function(x)
{
    avl = array(dim=c(5,4,2,2),
                dimnames=list(
                    ctyNames,
                    seasNames,
                    c("solar", "wind"),
                    c("mean", "sd")))
    for (ctyN in ctyNames) {
        for (s in seasNames) {
            for (tch in c("wind", "solar")) {
                avl[ctyN, s, tch, "mean"] =
                    mean(x[month(x) %in% seasons[,s], paste0(ctyN,tch)])
                avl[ctyN, s, tch, "sd"] =
                    sd(x[month(x) %in% seasons[,s], paste0(ctyN,tch)])
            }
        }
    }
    return(avl)
}




#' Difference in Availability: PCA and Empirical (Analysis only for DE)
#'
#' Print differences for DE and plot
#' 
#' @param avlPCA Array returned from \code{\link{avl.PCA}}
#' @param avl Array returned from \code{\link{avl.emp}}
#' 
diff.avl <- function(avlPCA
                    ,avl
                    ,plot2pdf = FALSE
                     )
{
    print("Sesonal availability (rel. difference to empirical):")
    print((avlPCA - avl["DE",,,])/avl["DE",,,]) 
    print("Yearly fraction (rel. difference to empirical):")
    print("Wind:")
    print(
    ( mean(avlPCA[,"wind","mean"]) - mean(avl["DE",,"wind","mean"]) )
    / mean(avlPCA[,"wind","mean"]))
    print("Solar:")
    print(
    ( mean(avlPCA[,"solar","mean"]) - mean(avl["DE",,"solar","mean"]) )
    / mean(avlPCA[,"solar","mean"])) 

    windows()
    y1 = stack(data.frame(avl["DE",,,"mean"]))[,1]
    sd1 = stack(data.frame(avl["DE",,,"sd"]))[,1]
    y2 = stack(data.frame(avlPCA[,,"mean"]))[,1]
    sd2 =  stack(data.frame(avlPCA[,,"sd"]))[,1]
    n = length(y1)
    x1 = 1:n
    x2 = x1 + 0.5
    par(mar=c(6.5,4,4,2)+0.1)
    plot(y1,
         ylim=c(0, 0.5),
         xlim=c(0.9, n+0.6),
         xaxt="n",
         pch=16,
         xlab="",
         cex=1.2,
         ylab="Seasonal availability"
         )
    lines(x2, y2, type="p", pch=16, cex=1.2, col="blue")
    techSeasName = paste(
        stack(data.frame(avl["DE",,,"mean"]))[,2],
        seasNames, sep=" ")
    axis(1,x1, rep("hist",n),las=2)
    axis(1,x2, rep("scn",n),las=2)
    axis(1,x1+0.25,techSeasName,las=2, line=2, tick=FALSE)
    arrows(x1,y1-sd1,x1,y1+sd1, code=3, length=0.02, angle = 90, lwd=2)
    arrows(x2,y2-sd2,x2,y2+sd2, code=3, length=0.02, angle = 90, col="blue", lwd=2)

    if (plot2pdf) dev.print(pdf,"avgAvlDiff_DE.pdf")
}
    






### Scenarios across countries for a season










#' Daily mean per season
#'
#' @param x xts-object 
#' @param seas Season. Season can be whole year (s = seasNames) 
#'
x.daily <- function(x,seas) {
    z = x[(month(x) %in% seasons[,seas]),]
    return(apply.daily(z, mean))
}




### Daily correlation matrix




#' Daily cor-matrix of wind a solar between countries
#'
#' @param  x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.
#' @param seas Season (seas = seasNames gives whole year)
#' 
#' @return correlation matrix
#' 
#' @examples
#' \dontrun{
#' daily.cor.cross(x.d, "SU")
#' daily.cor.cross(x.d, seasNames)
#' }
#'
daily.cor.cross <- function(x
                           ,seas = seasNames
                           ,plot2pdf = FALSE
                            )
{
    x.d = x.daily(x, seas)
    reOrder = c(paste0(ctyNames,"solar"), paste0(ctyNames,"wind"))
    c.d = cor(x.d)[reOrder,][,reOrder]
    windows()
    corrplot.mixed(c.d, tl.col=1,lower.col = "black", tl.cex = .8)
   
    #c.d2 = cor(x.d, method="kendall")[reOrder,][,reOrder] 
    #corrplot.mixed(c.d2, tl.col=1,lower.col = "black", tl.cex = .8)
    #round(c.d,2)
   
    if (plot2pdf) dev.print(pdf,"corDailyAllCty.pdf")
    return(c.d)
}



#' Reduce correlation matrix (remove cross-regional solar correlations)
#'
#' @param  x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.
#' @param seas Season (seas = seasNames gives whole year)
#' 
reduce.cor <- function(x, seas)
{
    x.d = x.daily(x, seas)
    reOrder = c(paste0(ctyNames,"solar"), paste0(ctyNames,"wind"))
    c.d = cor(x.d)[reOrder,][,reOrder] 
    # only correlation between solar and wind within a country:
    c.sw.d = c(
        c.d["ATwind","ATsolar"],
        c.d["CHwind","CHsolar"],
        c.d["DEwind","DEsolar"],
        c.d["FRwind","FRsolar"],
        c.d["ITwind","ITsolar"]
    )
    # cor-matrix only between winds:
    c.w.d = c.d[paste0(ctyNames, "wind"),
                paste0(ctyNames, "wind")]
    # build samller cor-matrix
    # where solar cross-country correlation are neglected:
    c.d.small = cbind(solar = c(1, c.sw.d), rbind(solar = c.sw.d, c.w.d) )

    return(c.d.small)
}



#' Simplified daily correlation matrix of wind and solar between countries
#'
#' Solar is no longer country-specific
#' 
#' @param  x Time series (xts class) containing hourly data with column
#'     names "ATsolar", etc.
#' @param seas Season (seas = seasNames gives whole year)
#'
#' @return correlation matrix
#'
daily.cor.cross.simple <- function(x
                                  ,seas = seasNames
                                  ,plot2pdf = FALSE
                                   )
{
    c.d.small = reduce.cor(x, seas)
    corrplot.mixed(c.d.small, tl.col=1,lower.col = "black", tl.cex = .8)
    if (plot2pdf) dev.print(pdf,paste0("corDailySimple_",seas,".pdf"))
    return(c.d.small)
}





### Lambda








#' Non-parametric fit of lambda
#'
#' Best fits are for single season, e.g. "SU"
#'
#' @param x An xts-time series (for the paper, this is a daily time-series, usually denoted by x.d)
#'
#' @examples
#' \dontrun{
#' lambda.non.param(x.d)
#' }
lambda.nonpar <- function(x) {
    flo = fitLambda(coredata(x), method="Schmid.Schmidt")
    fup = fitLambda(coredata(x), method="Schmid.Schmidt",
                     lower.tail = FALSE)
    corrplot.mixed(f.display(fup,names(x)))
    windows()
    corrplot.mixed(f.display(flo,names(x)))
}



#' Lambda from t-distribution
#'
#' Lambda from t-distribution is symmetric: counterfactual
#'
#' @param x An xts-time series (for the paper, this is a daily time-series, usually denoted by x.d)
#'
#' @examples
#' \dontrun{
#' lambda.t(x.d)
#' }
lambda.t <- function(x) {
    flo = fitLambda(coredata(x), method="t", verbose=TRUE)
    fup = fitLambda(coredata(x), method="t", verbose=TRUE,
                     lower.tail = FALSE)
    corrplot.mixed(f.display(fup$Lambda,names(x)))
    windows()
    corrplot.mixed(f.display(flo$Lambda,names(x))) # the same
    # dispersion matrix P = Sigma = neq correlation
    # f.display(fup$P,names(x.d)) # dispersion matrix of t-distribution
    # f.display(flo$P,names(x.d))
}

#' @describeIn lambda.t Round numeric argument and set row and column names
#'
#' @param y numeric argument
#' @return rounded argument with row and column names set
f.display = function (y, row.col.names) {
    z = round(y,2)
    rownames(z) = row.col.names
    colnames(z) = row.col.names
    return(z)
}





### Fit copulas



#'  Replace some data columns by a single, averaged column 
#'
#' @param x xts time-series
#' @param avgCols Data columns to be merged
#' @param newColname Name of merged column
#' 
avg.Cols <-  function(x
                     ,avgCols = c("ATsolar","DEsolar",
                                  "CHsolar","ITsolar","FRsolar")
                     ,newColname = "solar"
                      )
{
    z = cbind(x[, -which(names(x) %in% avgCols)],
              rowMeans(x[,avgCols]))
    names(z)[ncol(z)] = newColname
    return(z)
}


#' Fit t-copula to data, and sample
#'
#' Works only well with per season data. Attention: Symmetric copula.
#' 
#' @param x.d.small Time-series data (in paper: daily data)
#' @param sample.size Number of random samples 
#'
fit.tC.and.sample <- function (x.d.small
                              ,sample.size = 20
                               )
{
    tc.ml = fitCopula(tCopula(dim=6, dispstr="un")
                     ,as.matrix(x.d.small), method="ml"
                      #,start = c(0,0,0,46)
                      )
    # returns dispersion matrix P = Sigma = neq correlation 
    print(summary(tc.ml))
    #round(p2P(lambda(tc.ml@copula)),2) # wrong: fct. lambda only for bivariate

    ## sample t-copula from fitted copula:
    set.seed(5640)
    rtC = rCopula(sample.size, tc.ml@copula)
    colnames(rtC) = names(x.d.small)
    return(rtC)
}

 


#' Fit normal-copula from correlation matrix c.d.small, and sample
#'
#' @param c.d.small Time-series data (in paper: daily data)
#' @param sample.size Number of random samples 
#' 
fit.nC.and.sample <- function (c.d.small
                              ,sample.size = 20
                               )
{
    cop = normalCopula(P2p(c.d.small), dim = nrow(c.d.small), dispstr = "un")
    set.seed(5640)
    rnC = rCopula(sample.size, cop)
    colnames(rnC) = rownames(c.d.small)
    return(rnC)
}
    




#' Generate cross-country scenarios for seasonal data
#' 
#' @param scn.ws List of scenarios over countries of type \emph{scnBEM} (see help \code{\link{scn.clust}})
#' @param s Season
#' @param rC random samples as returned by \code{\link{fit.nC.and.sample}} or \code{\link{fit.tC.and.sample}} 
#' @param sample.size Sample size corresponding to rC (can in fact be determined by the dimensions of rC)
#' 
scn.cross <- function (scn.ws
                      ,s = "SU"
                      ,rC
                      ,sample.size)
{
    ## Get quantiles of scenarios per country per season
    
    # Bind the scenarios into a single data.frame:
    scn.ws.df = as.data.frame(do.call(rbind, scn.ws))
    cty.scn = names(scn.ws)
    # Split scenarios per country in into wind and solar, extract season,
    # and sum to daily data, and add rank column:
    scn.w.d = list()
    scn.s.d = list()
    for (ct in cty.scn) {
        # filter season and country:
        y = scn.ws.df[grepl(s, scn.ws.df$loadp) &
                      grepl(ct,scn.ws.df$country),]
        # split into wind and solar and sum over each day:
        yw = aggregate(cbind(prob,windAvl)~scn,
                       subset(y, select=-solarAvl), sum)
        ys = aggregate(cbind(prob,solarAvl)~scn,
                       subset(y, select=-windAvl), sum)
        yw$prob = yw$prob /24
        ys$prob = ys$prob /24
        # reorder by rank
        yw = yw[order(yw$windAvl),]
        ys = ys[order(ys$solarAvl),]
        # prepend the cumulative probability
        scn.w.d[[ct]] = as.data.frame(cbind(cumprob = cumsum(yw$prob), yw))
        scn.s.d[[ct]] = as.data.frame(cbind(cumprob = cumsum(ys$prob), ys))
    }
    sum(yw$prob)

    ## Match marginal distribution to sample copula values
    ## 1. Get scenario names per country matching the quantile (rank) of copual sample:
    ## 2. Collect them as a a list of countries in vectors "str.scn.w/s[[ct]]"

    str.scn.w = list()
    str.scn.s = list()
    for (ct in cty.scn) {
        p.w = scn.w.d[[ct]]$cumprob
        p.s = scn.s.d[[ct]]$cumprob 
        idx.w = sapply(rC[,paste0(ct, "wind")], get.nearest, y = p.w)
        idx.s = sapply(rC[,"solar"], get.nearest, y = p.s)
        str.scn.w[[ct]] = scn.w.d[[ct]] [idx.w, "scn"]
        str.scn.s[[ct]] = scn.s.d[[ct]] [idx.s, "scn"]
    }
    ## collect scenario values in scn.sw.new:
    scn.ws.new = list()
    j = 1
    for (i in 1:sample.size){
        for (ct in cty.scn) {
            # given season, country, get i'th 24h-scenario:
            # first filter season & country:
            y = scn.ws.df[grepl(s, scn.ws.df$loadp) & grepl(ct, scn.ws.df$country),]
            # get wind scenarios: 
            z1 =  y[str.scn.w[[ct]][i] == y$scn,
                    -which(names(y) %in% c("solarAvl")) ]
            # get solar scenarios: 
            z2 = y[str.scn.s[[ct]][i] == y$scn, c("scn","prob","solarAvl")]
            # bind columns of wind and solar 24h scenarios together
            scn.ws.new[[j]] = cbind(scn.new = paste0("o",rep(i,24)), z1, z2,
                                    stringsAsFactors=FALSE)
            j = j + 1
        }}
    j-1
    scn.ws.new = do.call(rbind, scn.ws.new)
    return(scn.ws.new)
}


#' describeIn scn.cross
get.nearest <- function(x,y) which.min(abs(x-y))


#' Plot cross-regional scenarios 
#'
#' @param scn.ws.new Cross-regional scenarios as returned by \code{\link{scn.cross}} 
#' 
scn.cross.plot <- function(scn.ws.new
                          ,s = "SU"
                          ,sample.size = 20
                          ,plot2pdf = FALSE)
{
    for (tech in c("solarAvl", "windAvl")) {
        # Reshape to wide format for wind and solar:
        dropTech = if (tech == "solarAvl") "windAvl" else "solarAvl" 
        scn.wide = reshape(scn.ws.new,
                           timevar="scn.new",
                           idvar = c("country","loadp"), direction="wide",
                           drop=c(dropTech, "prob", "scn"))
        windows()
        layout(matrix(seq_along(ctyNames),nrow = 1, ncol = 5),
               widths = c(1.4,1,1,1,1))
        par(oma=c(0,0,4,0))
        ymax = max(scn.wide[,c(-1,-2)])
        for (ct in ctyNames) {
            par(mar = if (ct=="AT") c(5, 4, 4, 1)+0.1 else c(5,1,4,1) + 0.1)
            plot.ts(scn.wide[scn.wide$country == ct,],
                    plot.type=c("single"),
                    lty=1:sample.size,
                    ylab= if (ct=="AT") tech else "",
                    xlab="",
                    ylim =c(0,ymax),
                    xaxt = 'n',
                    main = ct)
            axis(1,1:24,
                 label = paste0(" h", 1:24),
                 las = 2) # 2: perpendicular labels
            title(s, outer=TRUE)
        }
        if (plot2pdf) dev.copy(pdf,paste0("scn_", tech ,"_allCty_",s,".pdf"))
    }
}


#' Make cross-regional scenarios for all seasons
#'
#' @param scn.ws List of scenarios over countries of type \emph{scnBEM} (see help \code{\link{scn.clust}})
#' @param rC random samples from copula, as  rnC or rtC
#' @param sample.size Sample size corresponding to rC (CAN PERHAPS BE DETERMINED FROM rC)
#' @param tCop Take t-copula instead of guassian copula? (default: TRUE)
#'
#' @return List of cross-regional scenarios over the seasons
#' 
scn.cross.year <- function(x
                          ,scn.ws
                          ,sample.size
                          ,tCop = TRUE
                           )
{
    scn.ws.new.year = list() # list over seasons
    for (s in seasNames) {
        #browser()
        rC = if (tCop) {
                  x.d = x.daily(x, s)
                  x.d.small = avg.Cols(x.d)
                  fit.tC.and.sample(x.d.small, sample.size)
             } else {
                 c.d.small = reduce.cor(x, s)
                 fit.nC.and.sample(c.d.small, sample.size)
             }
        scn.ws.new = scn.cross(scn.ws, s, rC, sample.size) 
        scn.ws.new.year[[s]] = cbind(
            scn.ws.new[,c("scn.new","country","loadp")],   
            prob = rep(1/sample.size,nrow(scn.ws.new)),
            scn.ws.new[,c("windAvl","solarAvl")])
    }
    return(scn.ws.new.year)
}


#' Concatenate cross-regional scenarios over seasons and write to file
#'
#' @param scn.ws.new.year List of cross-regional scenarios over the seasons
#' @param sample.size Sample size: required to make correct scenario
#'     numbering: "o1", "o2", \ldots
#' @param tag String appended to file name
#'
#' @examples
#' \dontrun{
#' scn.cross.write(scn.ws.new.year, 20, "_tcopula")
#' }
scn.cross.write <- function(scn.ws.new.year
                            ,sample.size = 20
                            ,tag = "")
{
    z = as.data.frame(do.call(rbind, scn.ws.new.year))
    z$scn.new = paste0("o", rep(1:(sample.size*(length(seasNames))),
                                each=24*(length(ctyNames))))
    write.csv(z,
              paste0("scnWindSolarAllCty_", sample.size, tag, ".csv"),
              row.names = FALSE)
}








					
