#' Pollen Calendar by Different Methods from a Historical Pollen Database
#'
#' Function to calculate the pollen calendar from a historical database of several pollen types and using the most commonly used methods in the generation of the pollen calendars in the aerobiology field.
#'
#' @param data A \code{data.frame} object including the general database where calculation of the pollen calendar must be applied. This \code{data.frame} must include a first column in \code{Date} format and the rest of columns in \code{numeric} format belonging to each pollen type by column.
#' @param method A \code{character} string specifying the method applied to calculate and generate the pollen calendar. The implemented methods that can be used are: \code{"heatplot"}, \code{"violinplot"} or \code{"phenological"}. A more detailed information about the different methods for defining the pollen season may be consulted in \strong{Details}. The \code{method} argument will be \code{heatplot} by default.
#' @param n.types A \code{numeric} (\code{integer}) value specifying the number of the most abundant pollen types that must be represented in the pollen calendar. A more detailed information about the selection of the considered pollen types may be consulted in \strong{Details}. The \code{n.types} argument will be \code{15} types by default.
#' @param start.month A \code{numeric} (\code{integer}) value ranging \code{1_12} specifying the number of the month (January_December) when the beginning of the pollen calendar must be considered. This argument is only applicable for the \code{"heatplot"} method with \code{"daily"} period, for the \code{"phenological"} method with \code{"avg_before"} \code{average.method}, and for the \code{"violinplot"} method, and the rest of methods only may be generated from the January (\code{start.month = 1}). The \code{start.month} argument will be \code{1} (month of January) by default.
#' @param y.start,y.end A \code{numeric} (\code{integer}) value specifying the period selected to calculate the pollen calendar (start year _ end year). If \code{y.start} and \code{y.end} are not specified (\code{NULL}), the entire database will be used to generate the pollen calendar. The \code{y.start} and \code{y.end} arguments will be \code{NULL} by default.
#' @param perc1,perc2 A \code{numeric} value ranging \code{0_100}. These arguments are valid only for the \code{"phenological"} method. These values represent the percentage of the total annual pollen included in the pollen season, removing \code{(100_percentage)/2\%} of the total pollen before and after of the pollen season. Two percentages must be specified because of the definition of the "main pollination period" (\code{perc1}) and "early/late pollination" (\code{perc2}) based on the \code{"phenological"} method proposed by \emph{Werchan et al. (2018)}. The \code{perc1} argument will be \code{80} and \code{perc2} argument will be \code{99} by default. A more detailed information about the \code{phenological} method to generate the pollen calendar may be consulted in \strong{Details}.
#' @param th.pollen A \code{numeric} value specifying the minimum threshold of the average pollen concentration which will be used to generate the pollen calendar. Days below this threshold will not be considered. For the \code{"phenological"} method this value limits the "possible occurrence" period as proposed by \emph{Werchan et al. (2018)}. The \code{th.pollen} argument will be \code{1} by default. A more detailed information about the methods to generate the pollen calendar may be consulted in \emph{Details}.
#' @param average.method A \code{character} string specifying the moment of the application of the average. This argument is valid only for the \code{"phenological"} method. The implemented methods that can be used are: \code{"avg_before"} or \code{"avg_after"}. \code{"avg_before"} produces the average to the daily concentrations and then pollen season are calculated for all pollen types, this method is recommended as it is a more concordant method with respect to the rest of implemented methods. Otherwise, \code{"avg_after"} determines the pollen season for all years and all pollen types, and then a circular average is calculated for the start_dates and end_dates. The \code{average.method} argument will be \code{"avg_before"} by default.
#' @param period A \code{character} string specifying the interval time considered to generate the pollen calendar. This argument is valid only for the \code{"heatplot"} method. The implemented periods that can be used are: \code{"daily"} or \code{"weekly"}. \code{"daily"} selection produces a pollen calendar using daily averages during the year and \code{"weekly"} selection produces a pollen calendar using weekly averages during the year. The \code{period} argument will be \code{"daily"} by default.
#' @param method.classes A  \code{character} string specifying the method to define the classes used for classifying the average pollen concentrations to generate the pollen calendar. This argument is valid only for the  \code{"heatplot"} method. The implemented methods for defining classes are:  \code{"exponential"} and  \code{"custom"}. The  \code{method.classes} argument will be  \code{"exponential"} by default. A more detailed information about the methods to classify the average pollen concentrations to generate the pollen calendar may be consulted in \strong{Details}.
#' @param n.classes A \code{numeric} (\code{integer}) value specifying the number of classes that will be used for classifying the average pollen concentrations to generate the pollen calendar. This argument is valid only for the \code{"heatplot"} method and the classification by \code{method.classes = "custom"}. The \code{n.classes} argument will be \code{5} by default. A more detailed information about the methods to classify the average pollen concentrations to generate the pollen calendar may be consulted in \strong{Details}.
#' @param classes A \code{numeric} vector specifying the threshold established to define the different classes that will be used for classifying the average pollen concentrations to generate the pollen calendar. This argument is valid only for the \code{"heatplot"} method and the classification by \code{method.classes = "custom"}. The \code{classes} argument will be \code{c(25, 50, 100, 300)} by default. The number of specified classes must be equal to \code{n.classes _ 1} because of the maximum threshold will be automatically specified by the maximum value. A more detailed information about the methods to classify the average pollen concentrations to generate the pollen calendar may be consulted in \strong{Details}.
#' @param color A \code{character} string specifying the color to generate the graph showing the pollen calendar. This argument is valid only for the \code{"heatplot"} method. The implemented color palettes to generate the pollen calendar are: \code{"green"}, \code{"red"}, \code{"blue"}, \code{"purple"} or \code{"black"}. The \code{color} argument will be \code{"green"} by default.
#' @param interpolation A \code{logical} value. If \code{FALSE} the interpolation of the pollen data is not applicable. If \code{TRUE} an interpolation of the pollen series will be applied to complete the gaps before the calculation of the pollen calendar. The \code{interpolation} argument will be \code{TRUE} by default. A more detailed information about the interpolation method may be consulted in \strong{Details}.
#' @param int.method A \code{character} string specifying the method selected to apply the interpolation method in order to complete the pollen series. The implemented methods that may be used are: \code{"lineal"}, \code{"movingmean"}, \code{"spline"} or \code{"tseries"}. The \code{int.method} argument will be \code{"lineal"} by default.
#' @param na.remove A \code{logical} value specifying if \code{na} values must be remove for the pollen calendar or not. \code{n.remove = TRUE} by default.
#' @param legendname A \code{character} string specifying the title of the legend. By default is \code{"Pollen / m3"}.
#' @param result A \code{character} string specifying the output for the function. The implemented outputs that may be obtained are: \code{"plot"} and \code{"table"}. The argument \code{result} will be \code{"plot"} by default.
#' @param export.plot A \code{logical} value specifying if a plot with the pollen calendar saved in the working directory will be required or not. If \code{FALSE} graphical results will only be displayed in the active graphics window. If \code{TRUE} graphical results will be displayed in the active graphics window and also a \emph{pdf} file will be saved within the \emph{plot_AeRobiology} directory automatically created in the working directory. The \code{export.plot} will be \code{FALSE} by default.
#' @param export.format A \code{character} string specifying the format selected to save the pollen calendar plot. The implemented formats that may be used are: \code{"pdf"} and \code{"png"}. The \code{export.format} will be \code{"pdf"} by default.
#' @param ... Other additional arguments may be used to customize the exportation of the plots using \code{"pdf"} or \code{"png"} files and therefore arguments from \emph{pdf} and \emph{png} functions (\pkg{grDevices} package) may be implemented. For example, for \emph{pdf} files the user may custom the arguments: \code{width}, \code{height}, \code{family}, \code{title}, \code{fonts}, \code{paper}, \code{bg}, \code{fg}, \code{pointsize...}; and for \emph{png} files the user may custom the arguments: \code{width}, \code{height}, \code{units}, \code{pointsize}, \code{bg}, \code{res...}
#' @details This function allows to calculate and generate the pollen calendar using three different methods which are described below. The pollen calendar will be calculated and generated only for the period specified by the user using the \code{y.start} and \code{y.end} arguments, and for the specified number of the most abundant pollen types using the \code{n.types} argument by the user. The most abundant pollen types will be selected according to the highest average annual amounts of pollen registered by the pollen types during the considered period.
#' \itemize{
#' \item \code{"heatplot"} method. This pollen calendar is constructed based on the daily or weekly average of pollen concentrations, depending of the preferences of the user that may select \code{"daily"} or \code{"weekly"} as \code{period} argument. Then, these averages may be classified in different categories following different methods selected by the user according to the \code{method.classes} argument. If \code{method.classes = "exponential"} the user will apply the classification based on exponential classes proposed by \emph{Stix and Ferreti (1974)}, which has been commonly used in aerobiology for the generation of pollen calendars. The classification based on exponential method considers 11 classes (\code{1_2, 3_5, 6_11, 12_24, 25_49, 50_99, 100_199, 200_399, 400_799, 800_1600, >1600}). An example of this pollen calendar may be consulted in \emph{Rojo et al. (2016)}. This method to design pollen calendars is an adaptation from the pollen calendar proposed by \emph{Spieksma (1991)} who considered 10_day periods instead of daily or weekly periods. Otherwise, if \code{method.classes  = "custom"} the user may customize the classification according to the number of classes selected (\code{n.classes} argument) and the thresholds of the pollen concentrations used to define the classes (\code{classes} argument). Average values below the level of the \code{th.pollen} argument will be removed from the pollen calendar.
#' \item \code{"phenological"} method. This pollen calendar is based on phenological definition of the pollen season and adapted from the methodology proposed by \emph{Werchan et al. (2018)}. After to obtain daily average pollen concentrations for the most abundant pollen types different pollination periods were calculated using the daily averages. The main pollination period was calculated based on the percentage defined by \code{perc1} argument (selected by the user, 80\% by default) of the annual total pollen. For example, if \code{perc1 = 80} the beginning of the high season was marked when 10\% of the annual value was reached and the end was selected when 90\% was reached. In the case of the early/late pollination a total of the percentage defined by \code{perc2} argument (selected by the user, 99\% by default) of the annual total pollen will be registered during this period. For this kind of pollen calendar the \code{th.pollen} argument will define the "possible occurrence" period as adapted by \emph{Werchan et al. (2018)}, considering the entire period between the first and the last day when this pollen level is reached. In an alternative way the average may be carried out after to define the pollen seasons using \code{method_average = "avg_after"} (instead of \code{"avg_before"} by default). \code{"avg_after"} determines the pollen season for all years and all pollen types, and then an average for circular data is calculated from the start_dates and end_dates.
#' \item \code{"violinplot"} method. This pollen calendar is based on the pollen intensity and adapted from the pollen calendar published by \emph{ORourke (1990)}. In first time the daily averages of the pollen concentrations are calculated and then these averages are represented using the \emph{violin plot} graph. The shape of the \emph{violin plot} display the pollen intensity of the pollen types in a relative way i.e. the values will be calculated as relative measurements regarding to the most abundant pollen type in annual amounts. Therefore, this pollen calendar shows a relative comparison between the pollen intensity of the pollen types but without scales and units. Average values below the level of the \code{th.pollen} argument will be removed from the pollen calendar.
#' }
#' Pollen time series frequently have different gaps with no data and this fact could be a problem for the calculation of specific methods for defining the pollen season even providing incorrect results. In this sense by default a linear interpolation will be carried out to complete these gaps before to generate the pollen calendar. For more information to see the \code{\link{interpollen}} function.
#' @return This function returns different results:\cr
#' \code{plot} in the active graphics window displaying the pollen calendar generated by the user when \code{result = "plot"}. This plot may be included in an object by assignment operators.\cr
#' \code{data.frame} including the daily or weekly average pollen concentrations (according to the selection of the user) used to generate the pollen calendar. This \code{data.frame} will be returned when \code{result = "table"}.\cr
#' If \code{export.plot = TRUE} this plot displaying the pollen calendar will also be exported as file within the \emph{Plot_AeRobiology}" directory created in the working directory.\cr
#' If \code{export.plot = TRUE} and \code{export.format = pdf} a \emph{pdf} file of the pollen calendar will be saved within the \emph{plot_AeRobiology} directory created in the working directory. Additional characteristics may be incorporated to the exportation as \emph{pdf} file (see \pkg{grDevices} package)\cr
#' If \code{export.plot = TRUE} and \code{export.format = png} a \emph{png} file of the pollen calendar will be saved within the \emph{plot_AeRobiology} directory created in the working directory. Additional characteristics may be incorporated to the exportation as \emph{png} file (see \pkg{grDevices} package).
#' @references ORourke, M.K., 1990. Comparative pollen calendars from Tucson, Arizona: Durhamvs. Burkard samplers. \emph{Aerobiologia}, 6(2), p.136_140.
#' @references Rojo, J., Rapp, A., Lara, B., Sabariego, S., Fernandez_Gonzalez, F. and Perez_Badia, R., 2016. Characterisation of the airborne pollen spectrum in Guadalajara (central Spain) and estimation of the potential allergy risk. \emph{Environmental Monitoring and Assessment}, 188(3), p.130.
#' @references Spieksma, F.T.M., 1991. \emph{Regional European pollen calendars. Allergenic pollen and pollinosis in Europe}, pp.49_65.
#' @references Stix, E. and Ferretti, M.L., 1974. \emph{Pollen calendars of three locations in Western Germany. Atlas European des Pollens Allergisants}, pp.85_94.
#' @references Werchan, M., Werchan, B. and Bergmann, K.C., 2018. German pollen calendar 4.0_update based on 2011_2016 pollen data. \emph{Allergo Journal International}, 27, pp.69_71.
#' @seealso \code{\link{interpollen}}, \code{\link{calculate_ps}}
#' @examples data("munich_pollen")
#' @examples pollen_calendar(munich_pollen, method = "heatplot", interplation = FALSE)
#' @importFrom dplyr ungroup group_by mutate summarise_all
#' @importFrom circular circular mean.circular
#' @importFrom ggplot2 aes coord_flip element_text geom_tile geom_violin ggplot labs scale_fill_manual scale_x_continuous scale_x_date scale_y_date theme theme_bw theme_dark
#' @importFrom graphics abline barplot legend lines mtext par plot
#' @importFrom grDevices colorRampPalette dev.off pdf png recordPlot
#' @importFrom lubridate is.POSIXt
#' @importFrom scales date_format
#' @importFrom stats na.omit
#' @importFrom tidyr gather %>%
#' @export
pollen_calendar <- function (data,
          method = "heatplot", # "phenological" "violinplot" "heatplot"
          n.types = 15,
          start.month = 1,
          y.start = NULL,
          y.end = NULL,
          perc1 = 80,
          perc2 = 99,
          th.pollen = 1,
          average.method = "avg_before", # "avg_before" "avg_after"
          period = "daily", # "daily"  "weekly
          method.classes = "exponential", # "custom" "exponential"
          n.classes = 5,
          classes = c(25,50,100,300),
          color = "green", # "red" "green" "blue" "purple" "black"
          interpolation = TRUE,
          int.method = "lineal",
          na.remove = TRUE,
          result = "plot",
          export.plot = FALSE,
          export.format = "pdf",
          legendname = "Pollen grains / m3",...){

#############################################    CHECK THE ARGUMENTS       #############################

if(export.plot == TRUE){ifelse(!dir.exists(file.path("plot_AeRobiology")), dir.create(file.path("plot_AeRobiology")), FALSE)}
  data<-data.frame(data)
if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

if(method != "phenological" & method != "violinplot" & method != "heatplot") stop ("Please 'method' argument only accept values: 'phenological', 'violinplot' or 'heatplot'")

if(class(n.types) != "numeric") stop ("Please include only numeric values for 'n.types' argument indicating the number of pollen types which will be displayed")

if(class(start.month) != "numeric" | ((start.month %in% c(1,2,3,4,5,6,7,8,9,10,11,12)) == FALSE)) stop ("Please include only numeric integer values between 1-12 for 'start.month' argument indicating the start month for the pollen calendar")

if(average.method == "avg_after" & start.month != 1) stop ("'avg_after' only can be calculate when 'start.month' argument is equal to 1")

if(class(y.start) != "numeric" & !is.null(y.start)) stop ("Please include only numeric values for y.start argument indicating the start year considered")

if(class(y.end) != "numeric" & !is.null(y.end)) stop ("Please include only numeric values for 'y.end' argument indicating the end year considered")

if(class(perc1) != "numeric" | perc1 < 0 | perc1 > 100) stop ("Please include only numeric values between 0-100 for 'perc1' argument")

if(class(perc2) != "numeric" | perc2 < 0 | perc2 > 100) stop ("Please include only numeric values between 0-100 for 'perc2' argument")

if(perc1 > perc2) stop ("Please should be perc 1 < perc2")

if(class(th.pollen) != "numeric") stop ("Please include only numeric values for 'th.pollen' argument indicating the minimum averaged pollen concentration to be considered")

if(average.method != "avg_before" & average.method != "avg_after") stop ("Please average.method only accept values: 'avg_before' or 'avg_after'")

if(period != "daily" & period != "weekly") stop ("Please period only accept values: 'daily' or 'weekly'")

if(method.classes != "custom" & method.classes != "exponential") stop ("Please method.classes only accept values: 'custom' or 'exponential'")

if(class(n.classes) != "numeric") stop ("Please include only numeric values for n.classes argument indicating the number of classes which will be used for the plots")

if(class(classes) != "numeric") stop ("Please include only numeric values for classes argument indicating the thresholds used for classifying the average pollen concentration for the plots")

if ((length(classes) + 1) != n.classes) stop ("The number of specified classes must be equal to 'n.classes - 1' because of the maximum threshold will be automatically specified by the maximum value")

if(color != "red" & color != "blue" & color != "green" & color != "purple" & color != "black") stop ("Please 'color' argument only accept values: 'red', 'blue', 'green', 'purple' or 'black'")

if(class(interpolation) != "logical") stop ("Please include only logical values for interpolation argument")

if(int.method != "lineal" & int.method != "movingmean" & int.method != "spline" & int.method != "tseries") stop ("Please int.method only accept values: 'lineal', 'movingmean', 'spline' or 'tseries'")

if(result != "plot" & result != "table") stop ("Please result only accept values: 'plot' or 'table'")

if(class(export.plot) != "logical") stop ("Please include only logical values for export.plot argument")

  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])

if(ncol(data) - 1 == 1) stop ("The number of pollen types must be 2 at least")

if(ncol(data)-1 < n.types) {
  n.types = ncol(data)-1
  warning(paste("WARNING: the number of columns is smaller than 'n.types' argument. 'n.types' adjusted to", n.types))}

#############################################    MANAGEMENT OF THE DATABASE       #############################
  average_values<-data.frame()
perc1 <- 100 - perc1; perc2 <- 100 - perc2

if(interpolation == TRUE){data <- interpollen(data, method = int.method, plot = F)}

jd.month <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)

seasons <- unique(as.numeric(strftime(data[, 1], "%Y")))
if(is.null(y.start)){y.start <- min(seasons)}; if(is.null(y.end)){y.end <- max(seasons)}

data <- data[which(as.numeric(strftime(data[ ,1], "%Y")) >= y.start & as.numeric(strftime(data[ ,1], "%Y")) <= y.end), ]
data <- data.frame(date = data[ ,1], year = as.numeric(strftime(data[ ,1], "%Y")), jd = as.numeric(strftime(data[ ,1], "%j")), week = as.numeric(strftime(data[ ,1], "%W")), data[ ,-1])

#types <- ddply(data[ ,-c(1,3,4)], "year", function(x) colSums(x[-1], na.rm = T)) [-1] %>%
#  apply(., 2, function(x) mean(x, na.rm = T)) %>% .[order(., decreasing = TRUE)] %>%
#  names(.) %>%
#  .[1:n.types]
types <- data.frame(data[ ,-c(1,3,4)] %>%
  group_by(year) %>%
  summarise_all(sum, na.rm = TRUE))[-1] %>%
  apply(., 2, function(x) mean(x, na.rm = T)) %>% .[order(., decreasing = TRUE)] %>%
  names(.) %>%
  .[1:n.types]

#############################################    ORDER      ######################################################

  data.or <- data[ ,which(colnames(data) %in% c("jd", types))]

  #data.or <- ddply(data.or, "jd", function(x) colMeans(x[-1], na.rm = T)) %>%
  #.[ ,colnames(.) %in% c("jd", types)]
  data.or <- data.frame(data.or %>%
    group_by(jd) %>%
    summarise_all(mean, na.rm = T))

  data.or <- data.or[which(data.or$jd != 366), ]

  data.or[is.na(data.or)] <- NA #data <- data[-nrow(data), ]

  pc.df1 <- data.frame(type = NA, start = NA, end = NA)
  pc.df2 <- data.frame(type = NA, start = NA, end = NA)
  pc.df3 <- data.frame(type = NA, start = NA, end = NA)

  for (t in 1:length(types)){

    pollen.s <- na.omit(data.or[ ,which(colnames(data.or) %in% c("jd",types[t]))]) # Pollen data without NAs

    pollen.s$acum1 <- NA
    for (i in 1:nrow(pollen.s)) {
      if (i == 1) {
        pollen.s$acum1[i] = pollen.s[i, types[t]]
      }  else {
        pollen.s$acum1[i] = pollen.s[i, types[t]] + pollen.s$acum1[i-1]
      }
    }

    pollen.s$acum2 <- NA
    for (i in nrow(pollen.s):1) {
      if (i == nrow(pollen.s)) {
        pollen.s$acum2[i] = pollen.s[i, types[t]]
      }  else {
        pollen.s$acum2[i] = pollen.s[i, types[t]] + pollen.s$acum2[i+1]
      }
    }

    lim1 <- sum(pollen.s[types[t]])/100*(perc1/2)
    lim2 <- sum(pollen.s[types[t]])/100*(perc2/2)

    pc.df1[t,"type"] <- types[t]
    pc.df1[t,"start"] <- pollen.s$jd[which(pollen.s$acum1 >= lim1)][1]
    pc.df1[t,"end"] <- pollen.s$jd[which(pollen.s$acum2 < lim1)][1] - 1; pc.df1[t,"end"] <- pc.df1[t,"end"] - pc.df1[t,"start"] + 1

    pc.df2[t,"type"] <- types[t]
    pc.df2[t,"start"] <- pollen.s$jd[which(pollen.s$acum1 >= lim2)][1]
    pc.df2[t,"end"] <- pollen.s$jd[which(pollen.s$acum2 < lim2)][1] - 1; pc.df2[t,"end"] <- pc.df2[t,"end"] - pc.df2[t,"start"] + 1

    pc.df3[t,"type"] <- types[t]
    pc.df3[t,"start"] <- pollen.s$jd[which(pollen.s[ ,types[t]] >= th.pollen)][1]
    pc.df3[t,"end"] <- pollen.s$jd[which(pollen.s[ ,types[t]] >= th.pollen)][length(pollen.s$jd[which(pollen.s[ ,types[t]] >= th.pollen)])]; pc.df3[t,"end"] <- pc.df3[t,"end"] - pc.df3[t,"start"] + 1
  }

type.or <- pc.df2$type[order(pc.df2$start, decreasing = T)]

#############################################    METHODS OF CALENDAR      ###########################

if (method == "phenological" & average.method == "avg_before"){

  data <- data[ ,which(colnames(data) %in% c("jd", types))]

  #data <- ddply(data, "jd", function(x) colMeans(x[-1], na.rm = T)) %>%
  #  .[ ,colnames(.) %in% c("jd", types)]
  data <- data.frame(data %>%
        group_by(jd) %>%
          summarise_all(mean, na.rm = T))

  data <- data[which(data$jd != 366), ]

  if(start.month != 1){data <- rbind(data[which(data$jd == jd.month[start.month]):nrow(data), ], data[1:(which(data$jd == jd.month[start.month])-1), ])}

  data$jd1 <- 1:nrow(data)

  data[is.na(data)] <- NA #data <- data[-nrow(data), ]

  pc.df1 <- data.frame(type = NA, start = NA, end = NA)
  pc.df2 <- data.frame(type = NA, start = NA, end = NA)
  pc.df3 <- data.frame(type = NA, start = NA, end = NA)

  for (t in 1:length(types)){

  pollen.s <- na.omit(data[ ,which(colnames(data) %in% c("jd1",types[t]))]) # Pollen data without NAs

  pollen.s$acum1 <- NA
  for (i in 1:nrow(pollen.s)) {
    if (i == 1) {
      pollen.s$acum1[i] = pollen.s[i, types[t]]
    }  else {
      pollen.s$acum1[i] = pollen.s[i, types[t]] + pollen.s$acum1[i-1]
    }
  }

  pollen.s$acum2 <- NA
  for (i in nrow(pollen.s):1) {
    if (i == nrow(pollen.s)) {
      pollen.s$acum2[i] = pollen.s[i, types[t]]
    }  else {
      pollen.s$acum2[i] = pollen.s[i, types[t]] + pollen.s$acum2[i+1]
    }
  }

  lim1 <- sum(pollen.s[types[t]])/100*(perc1/2)
  lim2 <- sum(pollen.s[types[t]])/100*(perc2/2)


  pc.df1[t,"type"] <- types[t]
  pc.df1[t,"start"] <- pollen.s$jd1[which(pollen.s$acum1 >= lim1)][1]
  pc.df1[t,"end"] <- pollen.s$jd1[which(pollen.s$acum2 < lim1)][1] - 1; pc.df1[t,"end"] <- pc.df1[t,"end"] - pc.df1[t,"start"] + 1

  pc.df2[t,"type"] <- types[t]
  pc.df2[t,"start"] <- pollen.s$jd1[which(pollen.s$acum1 >= lim2)][1]
  pc.df2[t,"end"] <- pollen.s$jd1[which(pollen.s$acum2 < lim2)][1] - 1; pc.df2[t,"end"] <- pc.df2[t,"end"] - pc.df2[t,"start"] + 1

  pc.df3[t,"type"] <- types[t]
  pc.df3[t,"start"] <- pollen.s$jd1[which(pollen.s[ ,types[t]] >= th.pollen)][1]
  pc.df3[t,"end"] <- pollen.s$jd1[which(pollen.s[ ,types[t]] >= th.pollen)][length(pollen.s$jd[which(pollen.s[ ,types[t]] >= th.pollen)])]; pc.df3[t,"end"] <- pc.df3[t,"end"] - pc.df3[t,"start"] + 1
  }
}

######################################################################################################

if (method == "phenological" & average.method == "avg_after"){
  data <- data[ ,which(colnames(data) %in% c("date", types))]

  pc.df1 <- data.frame(type = NA, start = NA, end = NA)
  pc.df2 <- data.frame(type = NA, start = NA, end = NA)
  pc.df3 <- data.frame(type = NA, start = NA, end = NA)

  pollen.ps1 <- calculate_ps(data = data, method = "percentage", perc = (100 - perc1), plot = F, export.result = F, interpolation = F)
  pollen.ps2 <- calculate_ps(data = data, method = "percentage", perc = (100 - perc2), plot = F, export.result = F, interpolation = F)


  list.results <- list()

  for (t in 1:length(types)){

    seasons <- unique(as.numeric(strftime(data[, 1], "%Y")))

    st.ps <-NA; mx.ps <- NA; en.ps <- NA

    result.ps <- data.frame(seasons, st.jd = NA, en.jd = NA)

    for (j in 1:length(seasons)) {
      tryCatch({
      ye <- seasons[j]

      pollen.s <- na.omit(data[which(as.numeric(strftime(data[, 1], "%Y"))==ye), which(colnames(data) %in% c("date",types[t]))]) # Pollen data without NAs
      pollen.s$jdays <- as.numeric(strftime(pollen.s[, 1], "%j"))

      st.ps <- pollen.s$jdays[which(pollen.s[ ,types[t]] > 0)][1]
      en.ps <- pollen.s$jdays[which(pollen.s[ ,types[t]] > 0)][length(pollen.s$jdays[which(pollen.s[ ,types[t]] > 0)])]

      result.ps[j,"st.jd"] <- pollen.s$jdays[pollen.s$jdays == st.ps] # Start-date of PS (JD)
      result.ps[j,"en.jd"] <- pollen.s$jdays[pollen.s$jdays == en.ps] # End-date of PS (JD)

      print(paste(ye, types[t]))
      }, error = function(e){
        print(paste("Year", ye, types[t], ". Try to check the aerobiological data for this year"))})

    }
    list.results[[types[t]]] <- result.ps

  }

  nam.list <- names(list.results)
  for(l in 1:length(nam.list)){
    if(l == 1) {df.results <- data.frame(type = nam.list[l], list.results[[l]])
    } else {
      df.results <- rbind(df.results, data.frame(type = nam.list[l], list.results[[l]]))
    }
  }
  pollen.ps3 <- df.results

  for (t in 1:length(types)){
    tryCatch({
      pc.df1[t,"type"] <- types[t]
      st <- na.omit(pollen.ps1$st.jd[pollen.ps1$type == types[t]]) * 360 / 365
      pc.df1[t,"start"] <- round((as.numeric(mean.circular(circular(st, units = "degrees")))) * 365 / 360); if(pc.df1[t,"start"] <= 0) {pc.df1[t,"start"] <- pc.df1[t,"start"] + 365}
      en <- na.omit(pollen.ps1$en.jd[pollen.ps1$type == types[t]]) * 360 / 365
      pc.df1[t,"end"] <- round((as.numeric(mean.circular(circular(en, units = "degrees")))) * 365 / 360); if(pc.df1[t,"end"] <= 0) {pc.df1[t,"end"] <- pc.df1[t,"end"] + 365}; pc.df1[t,"end"] <- pc.df1[t,"end"] - pc.df1[t,"start"] + 1

      pc.df2[t,"type"] <- types[t]
      st <- na.omit(pollen.ps2$st.jd[pollen.ps2$type == types[t]]) * 360 / 365
      pc.df2[t,"start"] <- round((as.numeric(mean.circular(circular(st, units = "degrees")))) * 365 / 360); if(pc.df2[t,"start"] <= 0) {pc.df2[t,"start"] <- pc.df2[t,"start"] + 365}
      en <- na.omit(pollen.ps2$en.jd[pollen.ps2$type == types[t]]) * 360 / 365
      pc.df2[t,"end"] <- round((as.numeric(mean.circular(circular(en, units = "degrees")))) * 365 / 360); if(pc.df2[t,"end"] <= 0) {pc.df2[t,"end"] <- pc.df2[t,"end"] + 365}; pc.df2[t,"end"] <- pc.df2[t,"end"] - pc.df2[t,"start"] + 1

      pc.df3[t,"type"] <- types[t]
      st <- na.omit(pollen.ps3$st.jd[pollen.ps3$type == types[t]]) * 360 / 365
      pc.df3[t,"start"] <- round((as.numeric(mean.circular(circular(st, units = "degrees")))) * 365 / 360); if(pc.df3[t,"start"] <= 0) {pc.df3[t,"start"] <- pc.df3[t,"start"] + 365}
      en <- na.omit(pollen.ps3$en.jd[pollen.ps3$type == types[t]]) * 360 / 365
      pc.df3[t,"end"] <- round((as.numeric(mean.circular(circular(en, units = "degrees")))) * 365 / 360); if(pc.df3[t,"end"] <= 0) {pc.df3[t,"end"] <- pc.df3[t,"end"] + 365}; pc.df3[t,"end"] <- pc.df3[t,"end"] - pc.df3[t,"start"] + 1
    }, error = function(e){})
  }

}

######################################################################################################

if (method == "violinplot" | (method == "heatplot" & period == "daily")){

  data <- data[ ,which(colnames(data) %in% c("jd", types))]

  #data <- ddply(data, "jd", function(x) colMeans(x[-1], na.rm = T)) %>%
  #  .[ ,colnames(.) %in% c("jd", types)]
  data <- data.frame(data %>%
                       group_by(jd) %>%
                       summarise_all(mean, na.rm = T))

  data <- data[which(data$jd != 366), ]

  if(start.month != 1){data <- rbind(data[which(data$jd == jd.month[start.month]):nrow(data), ], data[1:(which(data$jd == jd.month[start.month])-1), ])}

  if(start.month != 1){data$jd <- seq(as.Date(strptime(paste0(as.character(data$jd),"-2017"), format = "%j-%Y"))[1], as.Date(strptime(paste0(as.character(data$jd),"-2018"), format = "%j-%Y"))[nrow(data)], by = "days")
  } else {
    data$jd <- seq(as.Date(strptime(paste0(as.character(data$jd),"-2017"), format = "%j-%Y"))[1], as.Date(strptime(paste0(as.character(data$jd),"-2017"), format = "%j-%Y"))[nrow(data)], by = "days")
  }

  #violin_data <- gather(data, key = variable, value = value, -jd, na.rm = TRUE)
  violin_data <- gather(data, key = variable, value = value, -jd)
  violin_data$value[violin_data$value < th.pollen] <- NA

  #if (method == "heatplot" & period == "daily"){violin_data <- na.omit(violin_data)}

  # create a weight variable for each variable.  dplyr will make this easy.
  violin_data <-
    violin_data %>%
    group_by(value) %>%
    mutate(wt = value / max(colSums(data[,-1], na.rm = TRUE))) %>%
    ungroup()

  data[ ,1] <- as.numeric(strftime(data[ ,1], "%j"))
}

######################################################################################################

if (method == "heatplot" & period == "weekly"){
  #data <- ddply(data[ ,-c(1:3)], "week", function(x) colMeans(x[-1], na.rm = TRUE)) %>%
  #  .[ ,colnames(.) %in% c("week", types)]
  data <- data.frame(data[ ,-c(1:3)] %>%
                          group_by(week) %>%
                          summarise_all(mean, na.rm = T)) %>%
    .[ ,colnames(.) %in% c("week", types)]

  heat.data <- gather(data, key = variable, value = value, -week)
  heat.data$value[heat.data$value < th.pollen] <- NA
  heat.data <- na.omit(heat.data)

  heat.data$variable <- factor(heat.data$variable, levels = as.character(pc.df2$type), ordered = TRUE)

}


#################################################################   PLOTS   ##########################


if (method == "phenological"){
  #ORDER

  pc.df3$type <- factor(pc.df3$type, levels = as.character(type.or), ordered = T); pc.df3 <- pc.df3[order(pc.df3$type), ]
  pc.df2$type <- factor(pc.df2$type, levels = as.character(type.or), ordered = T); pc.df2 <- pc.df2[order(pc.df2$type), ]
  pc.df1$type <- factor(pc.df1$type, levels = as.character(type.or), ordered = T); pc.df1 <- pc.df1[order(pc.df1$type), ]

  pc.df1 <- pc.df1[which(pc.df1$type %in% pc.df2$type & pc.df1$type %in% pc.df3$type), ]
  pc.df2 <- pc.df2[which(pc.df2$type %in% pc.df1$type & pc.df2$type %in% pc.df3$type), ]
  pc.df3 <- pc.df3[which(pc.df3$type %in% pc.df1$type & pc.df3$type %in% pc.df2$type), ]


  #PLOT

  par(mar = c(5,7,3,1), xpd = F)

  barplot(`colnames<-`(t(pc.df3[-1]), pc.df3[,1]), width = 1, space = 0, horiz = T, col=c("transparent","yellow"), border = NA, las = 2, xlim = c(1, 365), axes = F, main = paste0("Pollen calendar (Period ",y.start,"-",y.end, ")"), cex.main = 1.8, font.axis = 3)
  barplot(`colnames<-`(t(pc.df2[-1]), pc.df2[,1]), width = 1, space = 0, horiz = T, col=c("transparent","orange"), border = NA, add = T, las = 2, axes = F, font.axis = 3)
  barplot(`colnames<-`(t(pc.df1[-1]), pc.df1[,1]), width = 1, space = 0, horiz = T, col=c("transparent","red"), border = NA, add = T, las = 2, axes = F, font.axis = 3)

  abline (h=0:n.types, lwd = 1, col = "gray")

  lines(x = c(1,1) , y = c(0,nrow(pc.df1)), lwd = 3, col = 1); lines(x = c(365,365) , y = c(0,nrow(pc.df1)), lwd = 3, col = 1)
  lines(x = c(1,365) , y = c(0,0), lwd = 2, col = 1); lines(x = c(1,365) , y = c(nrow(pc.df1),nrow(pc.df1)), lwd = 2, col = 1)

  nam.mo <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  len.mo <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

  if(start.month == 1){
    for(m in c(start.month:12,1:(start.month-1))){
      if(m <= 11) {lines(x = c(sum(len.mo[1:m]),sum(len.mo[1:m])) , y = c(-1,n.types), lwd = 1, col = "gray", lty = 2)}
      if(m <= 11) {lines(x = c(sum(len.mo[1:m]),sum(len.mo[1:m])) , y = c(-1,0), lwd = 2, col = "black", lty = 1)}
      if(m == 1) {mtext (side = 1, adj = 0.028, text = nam.mo[m], cex = 1.2)}
      if(m > 1) {mtext (side = 1, adj = 0.028 + (0.086*(m-1)), text = nam.mo[m], cex = 1.2)}
    }
  } else {
    cont <- as.numeric()
    for(m in c(start.month:12,1:(start.month-1))){
      cont <-c(cont, m)
      if(length(cont) <= 11) {lines(x = c(sum(len.mo[cont]),sum(len.mo[cont])) , y = c(-1,n.types), lwd = 1, col = "gray", lty = 2)}
      if(length(cont) <= 11) {lines(x = c(sum(len.mo[cont]),sum(len.mo[cont])) , y = c(-1,0), lwd = 2, col = "black", lty = 1)}
      if(length(cont) == 1) {mtext (side = 1, adj = 0.028, text = nam.mo[m], cex = 1.2)}
      if(length(cont) > 1) {mtext (side = 1, adj = 0.028 + (0.086*(length(cont)-1)), text = nam.mo[m], cex = 1.2)}
    }
  }


  # legend("bottomleft", legend = c("Main pollination period", "Early/late pollination", "Possible occurrence"), col = c("red", "orange", "yellow"), pch = 15, bty = "n", cex = 1, x.intersp = 0.6, y.intersp = 0.5, inset = c(0.55,-0.2), xpd = TRUE)

plot.calendar <- recordPlot()

}

######################################################################################################

if (method == "violinplot"){

violin_data$variable <- factor(violin_data$variable, levels = as.character(type.or), ordered = T)

violin_data <- data.frame(violin_data)

if (na.remove == TRUE) {violin_data <- na.omit(violin_data)}

plot.calendar <- ggplot(violin_data) +
    aes(x = variable, y = jd, weight = wt) +
    geom_violin(fill = "yellow", colour = "yellow", size = 1)+
    coord_flip()+
    theme_dark()+
    labs(y = "", x = "", title = paste0("Pollen calendar (Period ",y.start,"-",y.end,")"))+
    scale_y_date(limits = c(min(violin_data[ ,1])[1], max(violin_data[ ,1])[1]), breaks = "1 month", labels = date_format("%b"))+
    theme(axis.text.y = element_text(size = 12, face = "bold.italic"), axis.text.x = element_text(size = 14, face = "bold"), title = element_text(size = 14, face = "bold"))
}

######################################################################################################

if (method == "heatplot" & period == "daily"){

violin_data$variable <- factor(violin_data$variable, levels = as.character(type.or), ordered = T)

if (color == "red") {colmin = "#fc9272"; colmax = "#67000d"}
if (color == "green") {colmin = "#c7e9c0"; colmax = "#00441b"}
if (color == "black") {colmin = "#d9d9d9"; colmax = "#000000"}
if (color == "blue") {colmin = "#d0d1e6"; colmax = "#023858"}
if (color == "purple") {colmin = "#dadaeb"; colmax = "#3f007d"}

if(method.classes == "custom"){
  classes = c(0, classes, max(violin_data$value, na.rm = TRUE)+100)
  lab.classes <- as.character()
  for (l in 2:length(classes)){
    if(l == 2) {lab.classes <- c(paste0("<",classes[l]))}
    if(l != 2 & l != length(classes)) {lab.classes <- c(lab.classes, paste0(classes[l-1],"-",classes[l]))}
    if(l == length(classes)) {lab.classes <- c(lab.classes, paste0(">",classes[l-1]))}
  }
}

if(method.classes == "exponential"){
  if (max(violin_data$value, na.rm = TRUE) > 1600 ) {n.classes = 11
  classes = c(0,2,5,11,24,49,99,199,399,799,1600,(max(violin_data$value, na.rm = TRUE)+100))
  lab.classes = c("1-2","3-5","6-11","12-24","25-49", "50-99", "100-199", "200-399", "400-799", "800-1600", ">1600")}

  if (max(violin_data$value, na.rm = TRUE) <= 1600 ) {n.classes = 10
  classes = c(0,2,5,11,24,49,99,199,399,799,1600)
  lab.classes = c("1-2","3-5","6-11","12-24","25-49", "50-99", "100-199", "200-399", "400-799", "800-1600")}

}

#violin_data <- data.frame(violin_data)
if (na.remove == TRUE) {violin_data <- na.omit(violin_data)}

plot.calendar <- ggplot(violin_data, aes(jd, variable)) +
    geom_tile(aes(fill = cut(value, breaks = classes, labels = lab.classes)))+
    scale_fill_manual(drop=FALSE, values=colorRampPalette(c(colmin, colmax))(n.classes), name = legendname, na.translate = F)+
    theme_bw()+
    labs(y = "", x = "", title = paste0("Pollen calendar (Period ",y.start,"-",y.end,")"))+
    #scale_x_date(limits = c(min(violin_data[ ,1], na.rm = T)[1], max(violin_data[ ,1], na.rm = T)[1]), breaks = "1 month", labels = date_format("%b"))+
    scale_x_date(breaks = "1 month", labels = date_format("%b"))+
    theme(axis.text.y = element_text(size = 12, face = "bold.italic"), axis.text.x = element_text(size = 14, face = "bold"), title = element_text(size = 14, face = "bold"))

}

######################################################################################################

if (method == "heatplot" & period == "weekly"){

heat.data$variable <- factor(heat.data$variable, levels = as.character(type.or), ordered = T)

if (color == "red") {colmin = "#fc9272"; colmax = "#67000d"}
if (color == "green") {colmin = "#c7e9c0"; colmax = "#00441b"}
if (color == "black") {colmin = "#d9d9d9"; colmax = "#000000"}
if (color == "blue") {colmin = "#d0d1e6"; colmax = "#023858"}
if (color == "purple") {colmin = "#dadaeb"; colmax = "#3f007d"}

if(method.classes == "custom"){
  classes = c(0, classes, max(heat.data$value, na.rm = TRUE)+100)
  lab.classes <- as.character()
  for (l in 2:length(classes)){
    if(l == 2) {lab.classes <- c(paste0("<",classes[l]))}
    if(l != 2 & l != length(classes)) {lab.classes <- c(lab.classes, paste0(classes[l-1],"-",classes[l]))}
    if(l == length(classes)) {lab.classes <- c(lab.classes, paste0(">",classes[l-1]))}
  }
}

if(method.classes == "exponential"){
  if (max(heat.data$value, na.rm = TRUE) > 1600 ) {n.classes = 11
  classes = c(0,2,5,11,24,49,99,199,399,799,1600,(max(heat.data$value, na.rm = TRUE)+100))
  lab.classes = c("1-2","3-5","6-11","12-24","25-49", "50-99", "100-199", "200-399", "400-799", "800-1600", ">1600")}

  if (max(heat.data$value, na.rm = TRUE) <= 1600 ) {n.classes = 10
  classes = c(0,2,5,11,24,49,99,199,399,799,1600)
  lab.classes = c("1-2","3-5","6-11","12-24","25-49", "50-99", "100-199", "200-399", "400-799", "800-1600")}

}

if (na.remove == TRUE) {heat.data <- na.omit(heat.data)}

plot.calendar <- ggplot(heat.data, aes(week, variable)) +
    geom_tile(aes(fill = cut(value, breaks = classes, labels = lab.classes)), colour = "white")+
    scale_fill_manual(drop=FALSE, values=colorRampPalette(c(colmin, colmax))(n.classes), name = legendname)+
    theme_bw()+
    labs(y = "", x = "Week of the year", title = paste0("Pollen calendar (Period ",y.start,"-",y.end,")"))+
    scale_x_continuous(breaks = seq(0,50,5), limits = c(0,53))+
    theme(axis.text.y = element_text(size = 12, face = "bold.italic"), axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16, face = "bold"), title = element_text(size = 14, face = "bold"))
}

#############################################    EXPORT RESULTS       #############################

if(export.plot == TRUE & export.format == "pdf") {
  pdf(paste0("plot_AeRobiology/pollen_calendar_",method,".pdf"), ...)
    if(method == "phenological"){print(plot.calendar)
      } else {
        plot(plot.calendar)
      }

  dev.off()}

if(export.plot == TRUE & export.format == "png") {
  png(paste0("plot_AeRobiology/pollen_calendar_",method,".png"), ...)
  if(method == "phenological"){print(plot.calendar)
  } else {
    plot(plot.calendar)
  }
  dev.off()
  png(paste0("plot_AeRobiology/credits.png"))

  dev.off()
  }

if (result == "plot") return(plot.calendar)
if (result == "table") {return(data)}

}
