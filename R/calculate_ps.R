#' Estimation of the Main Parameters of the Pollen Season
#'
#' Function to calculate the main parameters of the pollen season with regard to phenology and pollen intensity from a historical database of several pollen types. The function can use the most common methods in the definition of the pollen season.
#'
#' @param data A \code{data.frame} object including the general database where calculation of the pollen season must be applied. This \code{data.frame} must include a first column in \code{Date} format and the rest of columns in \code{numeric} format belonging to each pollen type by column.
#' @param method A \code{character} string specifying the method applied to calculate the pollen season and the main parameters. The implemented methods that can be used are: \code{"percentage"}, \code{"logistic"}, \code{"moving"}, \code{"clinical"} or \code{"grains"}. A more detailed information about the different methods for defining the pollen season may be consulted in \strong{Details}.
#' @param th.day A \code{numeric} value. The number of days whose pollen concentration is bigger than this threshold is calculated for each year and pollen type. This value will be obtained in the results of the function. The \code{th.day} argument will be \code{100} by default.
#' @param perc A \code{numeric} value ranging \code{0_100}. This argument is valid only for \code{method = "percentage"}. This value represents the percentage of the total annual pollen included in the pollen season, removing \code{(100_perc)/2\%} of the total pollen before and after of the pollen season. The \code{perc} argument will be \code{95} by default.
#' @param def.season A \code{character} string specifying the method for selecting the best annual period to calculate the pollen season. The pollen seasons may occur within the natural year between two years. The implemented options that can be used are: \code{"natural"}, \code{"interannual"} or \code{"peak"}. The \code{def.season} argument will be \code{"natural"} by default. A more detailed information about the different methods for selecting the best annual period to calculate the pollen season may be consulted in \strong{Details}.
#' @param reduction A \code{logical} value. This argument is valid only for the \code{"logistic"} method. If \code{FALSE} the reduction of the pollen data is not applicable. If \code{TRUE} a reduction of the peaks above a certain level (\code{red.level} argument) will be carried out before the definition of the pollen season. The \code{reduction} argument will be \code{FALSE} by default. A more detailed information about the reduction process may be consulted in \strong{Details}.
#' @param red.level A \code{numeric} value ranging \code{0_1} specifying the percentile used as level to reduce the peaks of the pollen series before the definition of the pollen season. This argument is valid only for the \code{"logistic"} method. The \code{red.level} argument will be \code{0.90} by default, specifying the percentile 90.
#' @param derivative A \code{numeric} (\code{integer}) value belonging to options of \code{4}, \code{5} or \code{6} specifying the derivative that will be applied to calculate the asymptotes which determines the pollen season using the \code{"logistic"} method. This argument is valid only for the \code{"logistic"} method. The \code{derivative} argument will be \code{5} by default. More information may be consulted in \strong{Details}.
#' @param man A \code{numeric} (\code{integer}) value specifying the order of the moving average applied to calculate the pollen season using the \code{"moving"} method. This argument is valid only for the \code{"moving"} method. The \code{man} argument will be \code{11} by default.
#' @param th.ma A \code{numeric} value specifying the threshold used for the \code{"moving"} method for defining the beginning and the end of the pollen season. This argument is valid only for the \code{"moving"} method. The \code{th.ma} argument will be \code{5} by default.
#' @param n.clinical A \code{numeric} (\code{integer}) value specifying the number of days which must exceed a given threshold (\code{th.pollen} argument) for defining the beginning and the end of the pollen season. This argument is valid only for the \code{"clinical"} method. The \code{n.clinical} argument will be \code{5} by default.
#' @param window.clinical A \code{numeric} (\code{integer}) value specifying the time window during which the conditions must be evaluated for defining the beginning and the end of the pollen season using the \code{"clinical"} method. This argument is valid only for the \code{"clinical"} method. The \code{window.clinical} argument will be \code{7} by default.
#' @param window.grains A \code{numeric} (\code{integer}) value specifying the time window during which the conditions must be evaluated for defining the beginning and the end of the pollen season using the \code{"grains"} method. This argument is valid only for the \code{"grains"} method. The \code{window.grains} argument will be \code{5} by default.
#' @param th.pollen A \code{numeric} value specifying the threshold that must be exceeded during a given number of days (\code{n.clinical} or \code{window.grains} arguments) for defining the beginning and the end of the pollen season using the \code{"clinical"} or \code{"grains"} methods. This argument is valid only for the \code{"clinical"} or \code{"grains"} methods. The \code{th.pollen} argument will be \code{10} by default.
#' @param th.sum A \code{numeric} value specifying the pollen threshold that must be exceeded by the sum of daily pollen during a given number of days (\code{n.clinical} argument) exceeding a given daily threshold (\code{th.pollen} argument) for defining the beginning and the end of the pollen season using the \code{"clinical"} method. This argument is valid only for the \code{"clinical"} method. The \code{th.sum} argument will be \code{100} by default.
#' @param type A \code{character} string specifying the parameters considered according to a specific pollen type for calculating the pollen season using the \code{"clinical"} method. The implemented pollen types that may be used are: \code{"birch"}, \code{"grasses"}, \code{"cypress"}, \code{"olive"} or \code{"ragweed"}. As result for selecting any of these pollen types the parameters \code{n.clinical}, \code{window.clinical}, \code{th.pollen} and \code{th.sum} will be automatically adjusted for the \code{"clinical"} method. If no pollen types are specified (\code{type = "none"}), these parameters will be considered by the user. This argument is valid only for the \code{"clinical"} method. The \code{type} argument will be \code{"none"} by default.
#' @param interpolation A \code{logical} value. If \code{FALSE} the interpolation of the pollen data is not applicable. If \code{TRUE} an interpolation of the pollen series will be applied to complete the gaps with no data before the calculation of the pollen season. The \code{interpolation} argument will be \code{TRUE} by default. A more detailed information about the interpolation method may be consulted in \strong{Details}.
#' @param int.method A \code{character} string specifying the method selected to apply the interpolation method in order to complete the pollen series. The implemented methods that may be used are: \code{"lineal"}, \code{"movingmean"}, \code{"spline"} or \code{"tseries"}. The \code{int.method} argument will be \code{"lineal"} by default.
#' @param maxdays A \code{numeric} (\code{integer} value) specifying the maximum number of consecutive days with missing data that the algorithm is going to interpolate. If the gap is bigger than the argument value, the gap will not be interpolated. Not valid with \code{int.method = "tseries"}. The \code{maxdays} argument will be \code{30} by default.
#' @param result A \code{character} string specifying the output for the function. The implemented outputs that may be obtained are: \code{"table"} and \code{"list"}. The argument \code{result} will be \code{"table"} by default.
#' @param plot A \code{logical} value specifying if a set of plots based on the definition of the pollen season will be shown in the plot history or not. If \code{FALSE} graphical results will not be shown. If \code{TRUE} a set of plots will be shown in the plot history. The \code{plot} argument will be \code{TRUE} by default.
#' @param export.plot A \code{logical} value specifying if a set of plots based on the definition of the pollen season will be saved in the workspace. If \code{TRUE} a \emph{pdf} file for each pollen type showing graphically the definition of the pollen season for each studied year will be saved within the \emph{plot_AeRobiology} directory created in the working directory. The \code{plot} argument will be \code{FALSE} by default.
#' @param export.result A \code{logical} value specifying if a excel file including all parameters for the definition of the pollen season saved in the working directoty will be required or not. If \code{FALSE} the results will not exported. If \code{TRUE} the results will be exported as \emph{xlsx} file including all parameters calculated from the definition of the pollen season within the \emph{table_AeRobiology} directory created in the working directory. The \code{export.result} argument will be \code{FALSE} by default.
#' @details This function allows to calculate the pollen season using five different methods which are described below. After calculating the start_date, end_date and peak_date for the pollen season all rest of parameters have been calculated as detailed in \strong{Value} section.
#' \itemize{
#'   \item \code{"percentage"} method. This is a commonly used method for defining the pollen season based on the elimination of a certain percentage in the beginning and the end of the pollen season \emph{(Nilsson and Persson, 1981; Andersen, 1991)}. For example if the pollen season is based on the 95\% of the total annual pollen (\code{"perc" = 95}), the start_date of the pollen season is marked as the day in which 2.5\% of the total pollen is registered and the end_date of the pollen season is marked as the day in which 97.5\% of the total pollen is registered.
#'   \item \code{"logistic"} method. This method was developed by \emph{Ribeiro et al. (2007)} and modified by \emph{Cunha et al. (2015)}. It is based on fitting annually a non_linear logistic regression model to the daily accumulated curve for each pollen type. This logistic function and the different derivatives were considered to calculate the start_date and end_date of the pollen season, based on the asymptotes when pollen amounts are stabilized on the beginning and the end of the accumulated curve. For more information about the method to see \emph{Ribeiro et al. (2007)} and \emph{Cunha et al. (2015)}. Three different derivatives may be used (\code{derivative} argument) \code{4}, \code{5} or \code{6} that represent from higher to lower restrictive criterion for defining the pollen season. This method may be complemented with an optional method for reduction the peaks values (\code{reduction = TRUE}), thus avoiding the effect of the great influence of extreme peaks. In this sense, peaks values will be cut below a certain level that the user may select based on a percentile analysis of peaks. For example, \code{red.level = 0.90} will cut all peaks above the percentile \code{90}.
#'   \item \code{"moving"} method. This method is proposed the first time by the authors of this package. The definition of the pollen season is based on the application of a moving average to the pollen series in order to obtain the general seasonality of the pollen curve avoiding the great variability of the daily fluctuations. Thus, the start_date and the end_date will be established when the curve of the moving average reaches a given pollen threshold (\code{th.ma} argument). Also the order of the moving average may be customized by the user (\code{man} argument). By default, \code{man} = 11 and \code{th.ma} = 5.
#'   \item \code{"clinical"} method. This method was proposed by \emph{Pfaar et al. (2017)}. It is based on the expert consensus in relation to pollen exposure and the relationship with allergic symptoms derived of the literature. Different periods may be defined by this method: the pollen season, the high pollen season and the high pollen days. The start_date and end_date of the pollen season were defined as a certain number of days (\code{n.clinical} argument) within a time window period (\code{window.clinical} argument) exceeding a certain pollen threshold (\code{th.pollen} argument) which summation is above a certain pollen sum (\code{th.sum} argument). All these parameters are established for each pollen type according to \emph{Pfaar et al. (2017)} and using the \code{type} argument these parameters may be automatically adjusted for the specific pollen types (\code{"birch"}, \code{"grasses"}, \code{"cypress"}, \code{"olive"} or \code{"ragweed"}). Furthermore, the user may change all parameters to do a customized definition of the pollen season. The start_date and end_date of the high pollen season were defined as three consecutive days exceeding a certain pollen threshold (\code{th.day} argument). The number of high pollen days will also be calculated exceeding this pollen threshold (\code{th.day}). For more information about the method to see \emph{Pfaar et al. (2017)}.
#' \item \code{"grains"} method. This method was proposed by \emph{Galan et al. (2001)} originally in olive pollen but after also applied in other pollen types. The start_date and end_date of the pollen season were defined as a certain number of days (\code{window.grains} argument) exceeding a certain pollen threshold (\code{th.pollen} argument). For more information about the method to see \emph{Galan et al. (2001)}.
#' }
#' The pollen season of the species may occur during the natural year (Calendar year: from 1. January to 31. December) or the start_date and the end_date of the pollen season may happen in two different natural years (or calendar years). This consideration has been taken into account and in this package different method for defining the period for calculating the pollen season have been implemented. In this sense, the \code{def.season} argument has been incorporated in three options:
#' \itemize{
#' \item \code{"natural"}: considering the pollination year as natural year from 1st January to 31th December for defining the start_dates and end_dates of the pollen season for each pollen types.
#' \item \code{"interannual"}: considering the pollination year from 1st June to 31th May for defining the start_dates and end_dates of the pollen season for each pollen types.
#' \item \code{"peak"}: considering a customized pollination year for each pollen types calculated as 6 previous months and 6 later months from the average peak_date.
#' }
#' Pollen time series frequently have different gaps with no data and this fact could be a problem for the calculation of specific methods for defining the pollen season even providing incorrect results. In this sense by default a linear interpolation will be carried out to complete these gaps before to define the pollen season (\code{interpolation = TRUE}). Additionally, the users may select other interpolation methods using the \code{int.method} argument, as \code{"lineal"}, \code{"movingmean"}, \code{"spline"} or \code{"tseries"}. For more information to see the \code{\link{interpollen}} function.
#' @return This function returns different results:\cr
#' \code{data.frame} when \code{result = "table"} including the main parameters of the pollen season with regard to phenology and pollen intensity as:
#' \itemize{
#'    \item \strong{type}: pollen type
#'    \item \strong{seasons}: year of the beginning of the season
#'    \item \strong{st.dt}: start_date (date)
#'    \item \strong{st.jd}: start_date (day of the year)
#'    \item \strong{en.dt}: end_date (date)
#'    \item \strong{en.jd}: end_date (day of the year)
#'    \item \strong{ln.ps}: length of the season
#'    \item \strong{sm.tt}: total sum
#'    \item \strong{sm.ps}: pollen integral
#'    \item \strong{pk.val}: peak value
#'    \item \strong{pk.dt}: peak_date (date)
#'    \item \strong{pk.jd}: peak_date (day of year)
#'    \item \strong{ln.prpk}: length of the pre_peak period
#'    \item \strong{sm.prpk}: pollen integral of the pre_peak period
#'    \item \strong{ln.pspk}: length of the post_peak period
#'    \item \strong{sm.pspk}: pollen integral of the post_peak period
#'    \item \strong{daysth}: number of days with more than 100 pollen grains
#'    \item \strong{st.dt.hs}: start_date of the High pollen season (date, only for clinical method)
#'    \item \strong{st.jd.hs}: start_date of the High pollen season (day of the year, only for clinical method)
#'    \item \strong{en.dt.hs}: end_date of the High pollen season (date, only for clinical method)
#'    \item \strong{en.jd.hs}: end_date of the High pollen season (day of the year, only for clinical method)
#' }
#' \code{list} when \code{result = "list"} including the main parameters of the pollen season, one pollen type by element
#' \code{plots} when \code{plot = TRUE} showing graphically the definition of the pollen season for each studied year in the plot history.
#' If \code{export.result = TRUE} this \code{data.frame} will also be exported as \emph{xlsx} file within the \emph{table_AeRobiology} directory created in the working directory. If \code{export.result = FALSE} the results will also be showed as list object named \code{list.ps}. \cr
#' If \code{export.plot = TRUE} a \emph{pdf} file for each pollen type showing graphically the definition of the pollen season for each studied year will be saved within the \emph{plot_AeRobiology} directory created in the working directory.
#' @references Andersen, T.B., 1991. A model to predict the beginning of the pollen season. \emph{Grana}, 30(1), pp.269_275.
#' @references Cunha, M., Ribeiro, H., Costa, P. and Abreu, I., 2015. A comparative study of vineyard phenology and pollen metrics extracted from airborne pollen time series. \emph{Aerobiologia}, 31(1), pp.45_56.
#' @references  Galan, C., Garcia_Mozo, H., Carinanos, P., Alcazar, P. and Dominguez_Vilches, E., 2001. The role of temperature in the onset of the \emph{Olea europaea} L. pollen season in southwestern Spain. \emph{International Journal of Biometeorology}, 45(1), pp.8_12.
#' @references Nilsson, S. and Persson, S., 1981. Tree pollen spectra in the Stockholm region (Sweden), 1973_1980. \emph{Grana}, 20(3), pp.179_182.
#' @references Pfaar, O., Bastl, K., Berger, U., Buters, J., Calderon, M.A., Clot, B., Darsow, U., Demoly, P., Durham, S.R., Galan, C., Gehrig, R., Gerth van Wijk, R., Jacobsen, L., Klimek, L., Sofiev, M., Thibaudon, M. and Bergmann, K.C., 2017. Defining pollen exposure times for clinical trials of allergen immunotherapy for pollen_induced rhinoconjunctivitis_an EAACI position paper. \emph{Allergy}, 72(5), pp.713_722.
#' @references Ribeiro, H., Cunha, M. and Abreu, I., 2007. Definition of main pollen season using logistic model. \emph{Annals of Agricultural and Environmental Medicine}, 14(2), pp.259_264.
#' @seealso \code{\link{interpollen}}, \code{\link{plot_ps}}
#' @examples data("munich_pollen")
#' @examples calculate_ps(munich_pollen, plot = TRUE, result = "table")
#' @importFrom circular circular mean.circular
#' @importFrom graphics abline lines  par plot
#' @importFrom grDevices dev.copy dev.off graphics.off recordPlot
#' @importFrom lubridate is.POSIXt yday
#' @importFrom stats aggregate coef complete.cases getInitial na.omit nls quantile SSlogis
#' @importFrom utils data head globalVariables
#' @importFrom writexl write_xlsx
#' @importFrom tidyr %>%
#' @export
calculate_ps <- function(data,
                          method = "percentage",
                          th.day = 100,
                          perc = 95,
                          def.season = "natural",
                          reduction = FALSE,
                          red.level = 0.90,
                          derivative = 5,
                          man = 11,
                          th.ma = 5,
                          n.clinical = 5,
                          window.clinical = 7,
                          window.grains = 5,
                          th.pollen = 10,
                          th.sum = 100,
                          type = "none",
                          interpolation = TRUE,
                          int.method = "lineal",
                          maxdays = 30,
                          result = "table",
                          plot = TRUE,
                          export.plot = FALSE,
                          export.result = FALSE) {

  ######################################################################################################

  if(class(export.plot) != "logical") stop ("Please include only logical values for export.plot argument")

  if(class(export.result) != "logical") stop ("Please include only logical values for export.result argument")

  if(export.plot == TRUE){ifelse(!dir.exists(file.path("plot_AeRobiology")), dir.create(file.path("plot_AeRobiology")), FALSE); plot = TRUE}

  if(export.result == TRUE){ifelse(!dir.exists(file.path("table_AeRobiology")), dir.create(file.path("table_AeRobiology")), FALSE)}

  data<-data.frame(data)
  if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

  if(class(th.day) != "numeric" | th.day < 0) stop ("Please include only numeric values >= 0 for th.day argument")

  if(class(perc) != "numeric" | perc < 0 | perc > 100) stop ("Please include only numeric values between 0-100 for perc argument")

  if(class(reduction) != "logical") stop ("Please include only logical values for reduction argument")

  if(class(red.level) != "numeric" | red.level < 0 | red.level > 1) stop ("Please include only numeric values between 0-1 for red.level argument")

  if(derivative != 4 & derivative != 5 & derivative != 6) stop ("Please derivative only accept values: 4, 5 or 6")

  if(class(man) != "numeric" | man < 0) stop ("Please include only numeric values > 0 for man argument")

  if(class(th.ma) != "numeric" | th.ma < 0) stop ("Please include only numeric values > 0 for th.ma argument")

  if(result != "table" & result != "list") stop ("Please result only accept values: 'table' or 'list'")

  if(class(plot) != "logical") stop ("Please include only logical values for plot argument")

  if(class(n.clinical) != "numeric" | n.clinical < 0) stop ("Please include only numeric values >= 0 for n.clinical argument")

  if(class(window.clinical) != "numeric" | window.clinical < 0) stop ("Please include only numeric values >= 0 for window.clinical argument")

  if(class(window.grains) != "numeric" | window.grains < 0) stop ("Please include only numeric values >= 0 for window.grains argument")

  if(class(th.pollen) != "numeric" | th.pollen < 0) stop ("Please include only numeric values >= 0 for th.pollen argument")

  if(class(th.sum) != "numeric" | th.sum < 0) stop ("Please include only numeric values >= 0 for th.sum argument")

  if(class(interpolation) != "logical") stop ("Please include only logical values for interpolation argument")

  if(int.method != "lineal" & int.method != "movingmean" & int.method != "spline" & int.method != "tseries") stop ("Please int.method only accept values: 'lineal', 'movingmean', 'spline' or 'tseries'")

  if(def.season != "natural" & def.season != "interannual" & def.season != "peak") stop ("Please def.season only accept values: 'natural', 'interannual' or 'peak'")

  if(method != "percentage" & method != "logistic" & method != "moving" & method != "clinical" & method != "grains") stop ("Please method only accept values: 'percentage', 'logistic', 'moving', 'clinical' or 'grains'")

  if(type != "none" & type != "birch" & type != "grasses" & type != "cypress" & type != "olive" & type != "ragweed") stop ("Please def.season only accept values: 'none', 'birch', 'grasses', 'cypress', 'olive' or 'ragweed'")

  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])

  ######################################################################################################
  SS<-c()
  a<-c()
  g<-c()
  b<-c()
  list.ps<-c()
  if(method == "clinical" & type == "birch") {n.clinical = 5; window.clinical = 7; th.pollen = 10; th.sum = 100; th.day = 100}
  if(method == "clinical" & type == "grasses") {n.clinical = 5; window.clinical = 7; th.pollen = 3; th.sum = 30; th.day = 50}
  if(method == "clinical" & type == "cypress") {n.clinical = 5; window.clinical = 7; th.pollen = 20; th.sum = 200; th.day = 100}
  if(method == "clinical" & type == "olive") {n.clinical = 5; window.clinical = 7; th.pollen = 20; th.sum = 200; th.day = 100}
  if(method == "clinical" & type == "ragweed") {n.clinical = 5; window.clinical = 7; th.pollen = 3; th.sum = 30; th.day = 50}

  th.rem = 100 - perc
  if(method == "clinical") {window.clinical <- window.clinical - 1}
  if(method == "grains") {window.grains <- window.grains - 1}

  if(interpolation == TRUE){data <- interpollen(data, method = int.method, maxdays = maxdays, plot = F)}

  type.name <- colnames(data)

  cols <- ncol(data)

  list.results <- list()

  for (t in 2: cols){

    ####################################################################################################
    ################################## Definition of the season ########################################

    if (def.season == "interannual"){

      n.year <- rep(NA, nrow(data))
      for(y in 1:length(n.year)){
        if (data[y, 1] <= as.Date(paste0(strftime(data[y, 1], "%Y"),"-05-31"))) {n.year[y] <- as.numeric(strftime(data[y, 1], "%Y"))-1} else {n.year[y] <- as.numeric(strftime(data[y, 1], "%Y"))}}
      data$n.year <- n.year

    } else if (def.season == "peak") {

      peaks <- (aggregate(data[ ,type.name[t]] ~ strftime(data[ ,1], "%Y"), FUN = max)) %>%
      .[which(.[ ,2] != 0), ]


      date.peak <- rep(NA, nrow(peaks))
      for (p in 1:nrow(peaks)){
        date.peak[p] <- data[which(strftime(data[ ,1], "%Y") == peaks[p,1] & data[ , type.name[t]] == peaks[p,2]), 1][1]}
      date.peak <- as.Date(date.peak, origin = "1970-01-01")
      date.peak <- as.numeric(strftime(date.peak, "%j")) * 360 / 365
      mean.cir <- round((as.numeric(mean.circular(circular(date.peak, units = "degrees")))) * 365 / 360)
      if(mean.cir <= 0) {mean.cir <- mean.cir + 365}

      date.year <- c(mean.cir - 182, mean.cir + 182); date.year <- date.year[date.year > 0 & date.year <= 365]
      n.year <- rep(NA, nrow(data))
      for(y in 1:length(n.year)){
        if (data[y, 1] <= as.Date(paste0(strftime(data[y, 1], "%Y"),substring(as.Date(strptime(date.year, "%j")), 5)))) {n.year[y] <- as.numeric(strftime(data[y, 1], "%Y")) - 1} else {n.year[y] <- as.numeric(strftime(data[y, 1], "%Y"))}}
      data$n.year <- n.year

    } else if (def.season == "natural") {
      n.year <- as.numeric(strftime(data[, 1], "%Y"))
      data$n.year <- n.year}

    ####################################################################################################

    seasons <- unique(data$n.year)

    st.ps <-NA; mx.ps <- NA; en.ps <- NA; st.hs <- NA; en.hs <- NA

    result.ps <- data.frame(seasons, st.dt = NA, st.jd = NA, en.dt = NA, en.jd = NA, ln.ps = NA, sm.tt = NA, sm.ps = NA, pk.val = NA, pk.dt = NA, pk.jd = NA, ln.prpk = NA, sm.prpk = NA, ln.pspk = NA, sm.pspk = NA, daysth = NA)

    if(method == "clinical") {result.ps <- data.frame(seasons, st.dt = NA, st.jd = NA, en.dt = NA, en.jd = NA, ln.ps = NA, sm.tt = NA, sm.ps = NA, pk.val = NA, pk.dt = NA, pk.jd = NA, ln.prpk = NA, sm.prpk = NA, ln.pspk = NA, sm.pspk = NA, daysth = NA, st.dt.hs = NA, st.jd.hs = NA, en.dt.hs = NA, en.jd.hs = NA)}

    #if(plot == TRUE){graphics.off()}
    #if(plot == TRUE){par(mfrow = c(ceiling(length(seasons)/4),4), mar = c(0.5,0.5,1,0.5))}
    par(mfrow = c(ceiling(length(seasons)/4),4), mar = c(0.5,0.5,1,0.5))

    for (j in 1:length(seasons)) {
      tryCatch({

        ye <- seasons[j]

        if(length(which(complete.cases(data[which(data$n.year==ye), c(1,t)]))) == 0 | sum(data[which(data$n.year==ye), type.name[t]], na.rm = T) == 0) {print(paste("Error: Year", ye, type.name[t], ". The entire aerobiological data are NA or 0")); next()}

        pollen.s <- na.omit(data[which(data$n.year==ye), c(1,t)]) # Pollen data without NAs
        #pollen.s <- data[which(data$n.year==ye), c(1,t)]
        pollen.s$jdays <- as.numeric(strftime(pollen.s[, 1], "%j"))
        
        pollen.s3 <- pollen.s

        if (method == "moving") { ### Method using the moving average limited by a threshold

          pollen.s1 <- data.frame(Date = seq(from = pollen.s[, 1][1], to = pollen.s[, 1][nrow(pollen.s)], by = "day"), pollen = 0, jdays = 0)
          colnames(pollen.s1) <- colnames(pollen.s)
          pollen.s1$jdays <- as.numeric(strftime(pollen.s1[, 1], "%j"))
          pollen.s1[which(pollen.s1[, 1] %in% pollen.s[, 1]), 2] <- pollen.s[which(pollen.s[, 1] %in% pollen.s1[, 1]), 2]
          pollen.s <- pollen.s1

              pollen.s$ma <- ma(pollen.s[, 2], man=man)
              pollen.s[which(is.nan(pollen.s[, "ma"])), "ma"] <- NA
              Peakmod <- pollen.s[which(pollen.s[, "ma"] == max(pollen.s[, "ma"], na.rm = T)), 1][1]

              if (is.na(Peakmod)|max(pollen.s[, "ma"], na.rm = T)<th.ma) {
                mx.ps <- NA
                st.ps <- NA
                en.ps <- NA

                #if(is.na(mx.ps) | is.na(st.ps) | is.na(en.ps)) {print(paste("Error: Year", ye, type.name[t], ". Try to check the aerobiological data for this year")); next()}

              } else{
                mx.ps <- pollen.s[which(pollen.s[, 2] == max(pollen.s[, 2], na.rm = T)), 1][1]%>%yday(.)
                #st.ps <- pollen.s[which((pollen.s[, "ma"] <= th.ma) & (yday(pollen.s[, 1]) < yday(Peakmod))), 1] %>% .[length(.)]%>%yday(.)
                st.ps <- pollen.s[which((pollen.s[, "ma"] < th.ma) & (pollen.s[, 1] < Peakmod)), 1] %>% .[length(.)]%>%yday(.)+1
                #en.ps <- pollen.s[which((pollen.s[, "ma"] <= th.ma) & (yday(pollen.s[, 1]) > yday(Peakmod))), 1] %>% .[1]%>%yday(.)
                en.ps <- pollen.s[which((pollen.s[, "ma"] < th.ma) & (pollen.s[, 1] > Peakmod)), 1] %>% .[1]%>%yday(.)-1

                if (is.numeric(st.ps) &
                    length(st.ps) >=1){} else{
                      st.ps <- NA }
                if (is.numeric(mx.ps) &
                    length(mx.ps) >=1) {} else{
                      mx.ps <- NA }
                if (is.numeric(en.ps) &
                    length(en.ps) >=1) {} else{
                      en.ps <- NA }

                #if(is.na(mx.ps) | is.na(st.ps) | is.na(en.ps)) {print(paste("Error: Year", ye, type.name[t], ". Try to check the aerobiological data for this year")); next()}
              }
  }

        if (method == "percentage") { ### Method using the percentage for removing data (Andersen et al. 1991)

          pollen.s$ndays <- 1:nrow(pollen.s)
          lim <- sum(pollen.s[type.name[t]])/100*(th.rem/2)

          pollen.s$acum1 <- NA
          for (i in 1:nrow(pollen.s)) {
            if (i == 1) {
              pollen.s$acum1[i] = pollen.s[i, type.name[t]]
            }  else {
              pollen.s$acum1[i] = pollen.s[i, type.name[t]] + pollen.s$acum1[i-1]
            }
          }

          pollen.s$acum2 <- NA
          for (i in nrow(pollen.s):1) {
            if (i == nrow(pollen.s)) {
              pollen.s$acum2[i] = pollen.s[i, type.name[t]]
            }  else {
              pollen.s$acum2[i] = pollen.s[i, type.name[t]] + pollen.s$acum2[i+1]
            }
          }

          st.ps <- pollen.s$jdays[which(pollen.s$acum1 >= lim)][1]
          en.ps <- pollen.s$jdays[which(pollen.s$acum2 < lim)][1] - 1
          mx.ps <- pollen.s$jdays[which(pollen.s[type.name[t]] == max(pollen.s[type.name[t]]))[1]]
        }


        if (method == "clinical"){ ### Method using the clinical method (Pfaar et al. 2017)

          #for (f in 1: (nrow(pollen.s) - window.clinical)){
          #  if(pollen.s[f, type.name[t]] >= th.pollen &
          #    length(which(pollen.s[f:(f + window.clinical), type.name[t]] >= th.pollen)) >= n.clinical &
          #    sum(head(sort(pollen.s[f:(f + window.clinical), type.name[t]][which(pollen.s[f:(f + window.clinical), type.name[t]] >= th.pollen)], decreasing = TRUE), n.clinical)) >= th.sum){
          #    st.ps <- pollen.s$jdays[f]; break()}
          #}
          
          for (f in 1: nrow(pollen.s)){
          if((nrow(pollen.s) - (f + window.clinical)) <= 0){
            if(pollen.s[f, type.name[t]] >= th.pollen &
               length(which(pollen.s[f:nrow(pollen.s), type.name[t]] >= th.pollen)) >= n.clinical &
               sum(head(sort(pollen.s[f:nrow(pollen.s), type.name[t]][which(pollen.s[f:nrow(pollen.s), type.name[t]] >= th.pollen)], decreasing = TRUE), n.clinical)) >= th.sum){
              st.ps <- pollen.s$jdays[f]; break()}
          } else {
            if(pollen.s[f, type.name[t]] >= th.pollen &
               length(which(pollen.s[f:(f + window.clinical), type.name[t]] >= th.pollen)) >= n.clinical &
               sum(head(sort(pollen.s[f:(f + window.clinical), type.name[t]][which(pollen.s[f:(f + window.clinical), type.name[t]] >= th.pollen)], decreasing = TRUE), n.clinical)) >= th.sum){
              st.ps <- pollen.s$jdays[f]; break()}
          }
          }
          
          #for (f in (1 + window.clinical) : nrow(pollen.s)){
          #    if(pollen.s[f, type.name[t]] >= th.pollen &
          #       length(which(pollen.s[f:(f - window.clinical), type.name[t]] >= th.pollen)) >= n.clinical &
          #       sum(head(sort(pollen.s[f:(f - window.clinical), type.name[t]][which(pollen.s[f:(f - window.clinical), type.name[t]] >= th.pollen)], decreasing = TRUE), n.clinical)) >= th.sum) {
          #      en.ps <- pollen.s$jdays[f]}
          #  }          
          
          for (f in 1 : nrow(pollen.s)){
            if((f - window.clinical) <= 0){
              
              if(pollen.s[f, type.name[t]] >= th.pollen &
                 length(which(pollen.s[f:1, type.name[t]] >= th.pollen)) >= n.clinical &
                 sum(head(sort(pollen.s[f:1, type.name[t]][which(pollen.s[f:1, type.name[t]] >= th.pollen)], decreasing = TRUE), n.clinical)) >= th.sum) {
                en.ps <- pollen.s$jdays[f]}
              
            } else {
              
              if(pollen.s[f, type.name[t]] >= th.pollen &
                 length(which(pollen.s[f:(f - window.clinical), type.name[t]] >= th.pollen)) >= n.clinical &
                 sum(head(sort(pollen.s[f:(f - window.clinical), type.name[t]][which(pollen.s[f:(f - window.clinical), type.name[t]] >= th.pollen)], decreasing = TRUE), n.clinical)) >= th.sum) {
                en.ps <- pollen.s$jdays[f]}
            }
          }

          
          for (h in 1: (nrow(pollen.s) - 2)){
            if(pollen.s[h, type.name[t]] >= th.day &
              length(which(pollen.s[h:(h + 2), type.name[t]] >= th.day)) == 3) {
              st.hs <- pollen.s$jdays[h]; break()}
          }

          for (h in (1 + 2) : nrow(pollen.s)){
            if(pollen.s[h, type.name[t]] >= th.day &
              length(which(pollen.s[h:(h - 2), type.name[t]] >= th.day)) == 3) {
              en.hs <- pollen.s$jdays[h]}

          }

          mx.ps <- pollen.s$jdays[which(pollen.s[type.name[t]] == max(pollen.s[type.name[t]]))[1]]
        }


        if (method == "grains"){ ### Method using the grains method (Galan et al. 1991)
          
          #for (f in 1: (nrow(pollen.s) - window.grains)){
          #  if(pollen.s[f, type.name[t]] >= th.pollen & length(which(pollen.s[f:(f + window.grains), type.name[t]] >= th.pollen)) >= (window.grains+1)) {st.ps <- pollen.s$jdays[f]; break()}
          #}
          
          for (f in 1: nrow(pollen.s)){
            if((nrow(pollen.s) - (f + window.grains)) <= 0){
              if(pollen.s[f, type.name[t]] >= th.pollen & length(which(pollen.s[f:nrow(pollen.s), type.name[t]] >= th.pollen)) >= (window.grains+1)) {st.ps <- pollen.s$jdays[f]; break()}
              
            } else {
              
              if(pollen.s[f, type.name[t]] >= th.pollen & length(which(pollen.s[f:(f + window.grains), type.name[t]] >= th.pollen)) >= (window.grains+1)) {st.ps <- pollen.s$jdays[f]; break()}
            }
              
          }
          
          #for (f in (1 + window.grains) : nrow(pollen.s)){
          #  if(pollen.s[f, type.name[t]] >= th.pollen & length(which(pollen.s[f:(f - window.grains), type.name[t]] >= th.pollen)) >= (window.grains+1)) {en.ps <- pollen.s$jdays[f]}
          #}
          
          for (f in 1 : nrow(pollen.s)){
            if((f - window.grains) <= 0){
              if(pollen.s[f, type.name[t]] >= th.pollen & length(which(pollen.s[f:1, type.name[t]] >= th.pollen)) >= (window.grains+1)) {en.ps <- pollen.s$jdays[f]}
            
            } else {
              
              if(pollen.s[f, type.name[t]] >= th.pollen & length(which(pollen.s[f:(f - window.grains), type.name[t]] >= th.pollen)) >= (window.grains+1)) {en.ps <- pollen.s$jdays[f]}
            }  
            
          }
          
          mx.ps <- pollen.s$jdays[which(pollen.s[type.name[t]] == max(pollen.s[type.name[t]]))[1]]
        }


        if (method == "logistic"){ ### Method using the accumulated curve (Ribeiro et al. 2004)

          pollen.s1 <- data.frame(Date = seq(from = pollen.s[, 1][1], to = pollen.s[, 1][nrow(pollen.s)], by = "day"), pollen = 0, jdays = 0)
          colnames(pollen.s1) <- colnames(pollen.s)
          pollen.s1$jdays <- as.numeric(strftime(pollen.s1[, 1], "%j"))
          pollen.s1[which(pollen.s1[, 1] %in% pollen.s[, 1]), 2] <- pollen.s[which(pollen.s[, 1] %in% pollen.s1[, 1]), 2]
          pollen.s <- pollen.s1

          pollen.s2 <- pollen.s

          if(reduction == TRUE){
            pollen.s[which( pollen.s[, type.name[t]] >= quantile(data[data[ ,type.name[t]] > 0, type.name[t]], red.level, na.rm = TRUE)), type.name[t]] <- quantile(data[data[ ,type.name[t]] > 0, type.name[t]], red.level, na.rm = TRUE)}

          #pollen.s[ , type.name[t]] <- round(ma(pollen.s[ , type.name[t]], order = 5))
          pollen.s <- na.omit(pollen.s)
          pollen.s$ndays <- 1:nrow(pollen.s)

          pollen.s$acum <- NA
          for (i in 1:nrow(pollen.s)) {
            if (i==1) {
              pollen.s$acum[i] = pollen.s[i, type.name[t]]
            }  else {
              pollen.s$acum[i] = pollen.s[i, type.name[t]] + pollen.s$acum[i-1]
            }
          }

          tryCatch({
            tryCatch({model.acum <- nls(acum~a*(1+exp(-(b+g*ndays)))^-1, start=list(a = max(pollen.s$acum),g = 1,b = -2), data=pollen.s)}, error = function(e){})
            tryCatch({
              SS <- getInitial(acum ~ SSlogis(ndays,alpha,xmid,scale), data = pollen.s)
              a <- SS["alpha"]
              g <- 1/SS["scale"]
              b <- SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
              model.acum <- nls(acum~a*(1+exp(-(b+g*ndays)))^-1, start=list(a = a, g = g, b = b), data=pollen.s)
            }, error = function(e){})
          }, error = function(e){print(paste("Error:", "Year:", ye, type.name[t]))})

          tryCatch({
            #summary(model.acum)
            c <- coef(model.acum) # a es c[1], g es c[2] y b es c[3]
            aj.acum <- function(t){c[1]*(1+exp(-(c[3]+c[2]*t)))^-1}
            pollen.s$model <- aj.acum(pollen.s$ndays)

            x2 <- -c[3]/c[2]; x2 <- round(x2)

            if(x2==0){
              SS <- getInitial(acum ~ SSlogis(ndays,alpha,xmid,scale), data = pollen.s)
              a <- SS["alpha"]
              g <- 1/SS["scale"]
              b <- SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
              model.acum <- nls(acum~a*(1+exp(-(b+g*ndays)))^-1, start=list(a = a, g = g, b = b), data=pollen.s)
              c <- coef(model.acum) # a es c[1], g es c[2] y b es c[3]
              aj.acum <- function(t){c[1]*(1+exp(-(c[3]+c[2]*t)))^-1}
              pollen.s$model <- aj.acum(pollen.s$ndays)
              x2 <- -c[3]/c[2]; x2 <- round(x2)
            }

            if (derivative == 4){
              # Resultados de la cuarta derivada
              x1 <- -(log(5+2*sqrt(6))+c[3])/c[2]; x1 <- floor(x1)
              x3 <- -(log(5-2*sqrt(6))+c[3])/c[2]; x3 <- ceiling(x3)
            } else if (derivative == 5) {
              # Resultados de la quinta derivada
              x1 <- -(c[3]+(log((sqrt(26*sqrt(105)+270)/2)+(sqrt(105)/2)+(13/2))))/c[2]; x1 <- floor(x1)
              x3 <- -(c[3]+(log(+(sqrt(105)/2)-(sqrt(26*sqrt(105)+270)/2)+(13/2))))/c[2]; x3 <- ceiling(x3)
            } else if (derivative == 6){
              # Resultados de la sexta derivada
              x1 <- -(c[3]+log(sqrt(336*sqrt(15)+1320)/2+3*sqrt(15)+14))/c[2]; x1 <- floor(x1)
              x3 <- -(c[3]+log((-sqrt(336*sqrt(15)+1320)/2)+3*sqrt(15)+14))/c[2]; x3 <- ceiling(x3)
            }

            st.ps <- pollen.s$jdays[which(pollen.s$ndays == x1)]
            en.ps <- pollen.s$jdays[which(pollen.s$ndays == x3)]
            mx.ps <- pollen.s2$jdays[which(pollen.s2[type.name[t]] == max(pollen.s2[type.name[t]]))[1]]

            if(reduction == TRUE){pollen.s <- pollen.s2}
            
          }, error = function(e){
            print(paste("Error: Year", ye, type.name[t], ". This method can not parametrize the model for this aerobiological data"))
          })
        }

        # Characteristics of the pollen season (PS)
        result.ps[j,"st.dt"] <- as.character(pollen.s[pollen.s$jdays == st.ps, 1]) # Start-date of PS
        result.ps[j,"st.jd"] <- pollen.s$jdays[pollen.s$jdays == st.ps] # Start-date of PS (JD)
        result.ps[j,"en.dt"] <- as.character(pollen.s[pollen.s$jdays == en.ps, 1]) # End-date of PS
        result.ps[j,"en.jd"] <- pollen.s$jdays[pollen.s$jdays == en.ps] # End-date of PS (JD)
        result.ps[j,"ln.ps"] <- round(as.numeric(pollen.s[pollen.s$jdays == en.ps, 1] - pollen.s[pollen.s$jdays == st.ps, 1])) + 1 # Length of PS
        result.ps[j,"sm.tt"] <- sum(pollen.s[type.name[t]]) # Total sum
        result.ps[j,"sm.ps"] <- sum(pollen.s[which(pollen.s$jdays == st.ps) : which(pollen.s$jdays == en.ps), type.name[t]]) # Pollen Index
        result.ps[j,"pk.val"] <- max(pollen.s[type.name[t]]) # Peak value
        result.ps[j,"pk.dt"] <- as.character(pollen.s[pollen.s$jdays == mx.ps, 1]) # Peak date
        result.ps[j,"pk.jd"] <- pollen.s$jdays[pollen.s$jdays == mx.ps] # Peak date (JD)
        result.ps[j,"ln.prpk"] <- round(as.numeric(pollen.s[pollen.s$jdays == mx.ps, 1] - pollen.s[pollen.s$jdays == st.ps, 1]))+1 # Lengh of pre-peak period
        result.ps[j,"sm.prpk"] <- sum(pollen.s[which(pollen.s$jdays == st.ps) : which(pollen.s$jdays == mx.ps), type.name[t]]) # Pre-peak sum
        result.ps[j,"ln.pspk"] <- round(as.numeric(pollen.s[pollen.s$jdays == en.ps, 1] - pollen.s[pollen.s$jdays == mx.ps, 1])) # Lengh of post-peak period
        result.ps[j,"sm.pspk"] <- sum(pollen.s[(which(pollen.s$jdays == mx.ps)+1) : which(pollen.s$jdays == en.ps), type.name[t]]) # Post-peak sum
        result.ps[j,"daysth"] <- length(which(pollen.s[type.name[t]] >= th.day))
        if (method == "clinical" & is.numeric(st.hs)){result.ps[j,"st.dt.hs"] <- as.character(pollen.s[pollen.s$jdays == st.hs, 1])} # Start-date of HPS
        if (method == "clinical" & is.numeric(st.hs)){result.ps[j,"st.jd.hs"] <- pollen.s$jdays[pollen.s$jdays == st.hs]} # Start-date of HPS (JD)
        if (method == "clinical" & is.numeric(en.hs)){result.ps[j,"en.dt.hs"] <- as.character(pollen.s[pollen.s$jdays == en.hs, 1])} # End-date of HPS
        if (method == "clinical" & is.numeric(en.hs)){result.ps[j,"en.jd.hs"] <- pollen.s$jdays[pollen.s$jdays == en.hs]} # End-date of HPS (JD)

        if (which(pollen.s$jdays == st.ps)-20 > 1) {pollen.s <- pollen.s[-(1:(which(pollen.s$jdays == st.ps)-20)), ]}
        if (which(pollen.s$jdays == en.ps)+20 < nrow(pollen.s)) {pollen.s <- pollen.s[-((which(pollen.s$jdays == en.ps)+20):nrow(pollen.s)), ]}

        if(plot == TRUE){

          if(method == "percentage" | method == "logistic" | method == "clinical" | method == "grains"){
          plot(y = pollen.s[ ,type.name[t]], x = 1:length(pollen.s$jdays), type = "l", main = paste(toupper(type.name[t]), ye), ylab = expression(paste("Granos de polen/m"^"3")), xlab = "Julian days", xaxt = "n", yaxt = "n", lwd = 1.5)
          #plot(y = pollen.s[ ,type.name[t]], x = pollen.s$jdays, type = "l", main = paste(toupper(type.name[t]), ye), ylab = expression(paste("Granos de polen/m"^"3")), xlab = "Julian days", xaxt = "n", yaxt = "n", lwd = 1.5)
          par(new=T)
          #plot(pollen.s$acum, x = pollen.s$ndays, type="l", col="red")
          #plot(pollen.s$acum, type="l", col="red")
          #lines(pollen.s$ndays, aj.acum(pollen.s$ndays), col="blue")
          ymax <- max(pollen.s[ ,type.name[t]])+100
          #lines(x = c(result.ps[j,"st.jd"], result.ps[j,"st.jd"]), y = c(0,ymax), col = 2, lty = 2, lwd = 2)
          lines(x = c(which(pollen.s$jdays == result.ps[j,"st.jd"]), which(pollen.s$jdays == result.ps[j,"st.jd"])), y = c(0,ymax), col = 2, lty = 2, lwd = 2)
          #lines(x=c(mx.ps, mx.ps), y=c(0,ymax), col=2, lty=2, lwd=2)
          #lines(x = c(result.ps[j,"en.jd"], result.ps[j,"en.jd"]), y = c(0,ymax), col = 2, lty = 2, lwd = 2)
          lines(x = c(which(pollen.s$jdays == result.ps[j,"en.jd"]), which(pollen.s$jdays == result.ps[j,"en.jd"])), y = c(0,ymax), col = 2, lty = 2, lwd = 2)
          #lines(x = c(result.ps[j,"pk.jd"], result.ps[j,"pk.jd"]), y = c(0,ymax), col = 3, lty = 2, lwd = 2)
          lines(x = c(which(pollen.s$jdays == result.ps[j,"pk.jd"]), which(pollen.s$jdays == result.ps[j,"pk.jd"])), y = c(0,ymax), col = 3, lty = 2, lwd = 2)
          } else {

            #plot(y = pollen.s[ ,2], x = yday(pollen.s[ ,1]), type = "l", main = paste(toupper(type.name[t])," ", ye), xaxt = "n", yaxt = "n", lwd = 1.5)
            plot(y = pollen.s[ ,2], x = 1:length(yday(pollen.s[ ,1])), type = "l", main = paste(toupper(type.name[t])," ", ye), xaxt = "n", yaxt = "n", lwd = 1.5)
            #lines(pollen.s[ ,"ma"], x = yday(pollen.s[ ,1]), col="red")
            lines(pollen.s[ ,"ma"], x = 1:length(yday(pollen.s[ ,1])), col="red")
            #abline(v = c(st.ps,en.ps), col = 4, lty = 2, lwd = 2)
            abline(v = c(which(yday(pollen.s[ ,1]) == st.ps), which(yday(pollen.s[ ,1]) == en.ps)), col = 4, lty = 2, lwd = 2)
            abline(h = th.ma, col = 3, lty = 2, lwd = 2)
          }

          #plot.season <- recordPlot()
        }

        print(paste(ye, type.name[t]))
      }, error = function(e){
        print(paste("Year", ye, type.name[t], ". Try to check the aerobiological data for this year"))
      })
        result.ps[j,"sm.tt"] <- sum(pollen.s3[type.name[t]])
        result.ps[j,"pk.val"] <- max(pollen.s3[type.name[t]])
        result.ps[j,"pk.dt"] <- as.character(pollen.s3[which(pollen.s3[type.name[t]] == max(pollen.s3[type.name[t]]))[1], 1]) # Peak date
        result.ps[j,"pk.jd"] <- pollen.s3$jdays[which(pollen.s3[type.name[t]] == max(pollen.s3[type.name[t]]))[1]] # Peak date (JD)


      st.ps <-NA; mx.ps <- NA; en.ps <- NA; st.hs <- NA; en.hs <- NA
    }
    result.ps$st.dt <- as.Date(strptime(result.ps$st.dt, "%Y-%m-%d"))
    result.ps$en.dt <- as.Date(strptime(result.ps$en.dt, "%Y-%m-%d"))
    result.ps$pk.dt <- as.Date(strptime(result.ps$pk.dt, "%Y-%m-%d"))

    if(export.plot == TRUE){plot.season <- recordPlot()}

    list.results[[type.name[t]]] <- result.ps

    par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))

    if(export.plot == TRUE){
      if(method == "percentage") {dev.copy(pdf, paste0("plot_AeRobiology/plot_",type.name[t],"_",method,perc,".pdf")); plot.season; dev.off()}
      if(method == "logistic") {dev.copy(pdf, paste0("plot_AeRobiology/plot_",type.name[t],"_",method,derivative,".pdf")); plot.season;  dev.off()}
      if(method == "moving") {dev.copy(pdf, paste0("plot_AeRobiology/plot_",type.name[t],"_",method,man,"_",th.ma,".pdf")); plot.season; dev.off()}
      if(method == "clinical") {dev.copy(pdf, paste0("plot_AeRobiology/plot_",type.name[t],"_",method,n.clinical,"_",window.clinical+1,"_",th.pollen,"_",th.sum,".pdf")); plot.season; dev.off()}
      if(method == "grains") {dev.copy(pdf, paste0("plot_AeRobiology/plot_",type.name[t],"_",method,window.grains+1,"_",th.pollen,".pdf")); plot.season; dev.off()}

      }
  }

  par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
  #if(plot == TRUE){graphics.off()}

  nam.list <- names(list.results)

  for(l in 1:length(nam.list)){
      if(l == 1) {df.results <- data.frame(type = nam.list[l], list.results[[l]])
    } else {
      df.results <- rbind(df.results, data.frame(type = nam.list[l], list.results[[l]]))
    }
  }

  list.ps <- list.results

  list.results [["Information"]] <- data.frame(
        Attributes = c("Pollen type", "seasons", "st.dt", "st.jd", "en.dt", "en.jd", "ln.ps", "sm.tt", "sm.ps", "pk.val", "pk.dt", "pk.jd", "ln.prpk", "sm.prpk", "ln.pspk", "sm.pspk", "daysth", "st.dt.hs", "st.jd.hs", "en.dt.hs", "en.jd.hs", "", "", "Method of the pollen season", "Interpolation", "Interpolation method", "Method for defining period", "", "Aditional arguments", "perc", "th.day", "reduction", "red.level", "derivative", "man", "th.ma", "n.clinical", "window.clinical", "window.grains", "th.pollen", "th.sum", "type", "", "", "Package", "Authors"),
        Description = c("Pollen type", "Year of the beginning of the season", "Start-date (date)", "Start-date (day of the year)", "End-date (date)", "End-date (day of the year)", "Length of the season", "Total sum", "Pollen integral", "Peak value", "Peak-date (date)", "Peak-date (day of year)", "Length of the pre-peak period", "Pollen integral of the pre-peak period", "Length of the post-peak period", "Pollen integral of the post-peak period", paste0("Number of days with more than ", th.day, " pollen grains"), "Start-date of the High pollen season (date, only for clinical method)", "Start-date of the High pollen season (day of the year, only for clinical method)", "End-date of the High pollen season (date, only for clinical method)", "End-date of the High pollen season (day of the year, only for clinical method)", "", "", method, interpolation, int.method, def.season, "", "", perc, th.day, reduction, red.level, derivative, man, th.ma, n.clinical, window.clinical, window.grains, th.pollen, th.sum, type,  "", "", "AeRobiology", "Jesus Rojo, Antonio Picornell & Jose Oteros"))

  if (export.result == TRUE) {
    if(method == "percentage") {write_xlsx(list.results, paste0("table_AeRobiology/table_ps_",method,perc,".xlsx"))}
    if(method == "logistic") {write_xlsx(list.results, paste0("table_AeRobiology/table_ps_",method,derivative,".xlsx"))}
    if(method == "moving") {write_xlsx(list.results, paste0("table_AeRobiology/table_ps_",method,man,"_",th.ma,".xlsx"))}
    if(method == "clinical") {write_xlsx(list.results, paste0("table_AeRobiology/table_ps_",method,n.clinical,"_",window.clinical+1,"_",th.pollen,"_",th.sum,".xlsx"))}
    if(method == "grains") {write_xlsx(list.results, paste0("table_AeRobiology/table_ps_",method,window.grains+1,"_",th.pollen,".xlsx"))}
    }

  par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))

  if (result == "table") {return(df.results)}
  if (result == "list") {return(list.ps)}
}
