#' Plot of the Relative Abundance of the Pollen Types
#'
#' Generates a barplot based on the relative abundance (as percentage) in the air of the pollen types with respect to the total amounts
#'
#' @param data A \code{data.frame} object including the general database where calculation of the pollen season must be applied. This \code{data.frame} must include a first column in \code{Date} format and the rest of columns in \code{numeric} format belonging to each pollen type by column.
#' @param n.types A \code{numeric} (\code{integer}) value specifying the number of the most abundant pollen types that must be represented in the plot of the relative abundance. A more detailed information about the selection of the considered pollen types may be consulted in \strong{Details}. The \code{n.types} argument will be \code{15} types by default.
#' @param y.start,y.end A \code{numeric} (\code{integer}) value specifying the period selected to calculate relative abundances of the pollen types (\code{start year _ end year}). If \code{y.start} and \code{y.end} are not specified (\code{NULL}), the entire database will be used to generate the pollen calendar. The \code{y.start} and \code{y.end} arguments will be \code{NULL} by default.
#' @param interpolation A \code{logical} value. If \code{FALSE} the interpolation of the pollen data is not applicable. If \code{TRUE} an interpolation of the pollen series will be applied to complete the gaps with no data before the calculation of the pollen season. The \code{interpolation} argument will be \code{TRUE} by default. A more detailed information about the interpolation method may be consulted in \strong{Details}.
#' @param int.method A \code{character} string specifying the method selected to apply the interpolation method in order to complete the pollen series. The implemented methods that may be used are: \code{"lineal"}, \code{"movingmean"}, \code{"spline"} or \code{"tseries"}. The \code{int.method} argument will be \code{"lineal"} by default.
#' @param col.bar A \code{character} string specifying the color of the bars to generate the graph showing the relative abundances of the pollen types. The \code{color} argument will be \code{#E69F00} by default, but any color may be selected.
#' @param type.plot A \code{character} string specifying the type of plot selected to show the plot showing the relative abundance of the pollen types. The implemented types that may be used are: \code{static} generates a static \strong{ggplot} object and \code{dynamic} generates a dynamic \strong{plotly} object.
#' @param result A \code{character} string specifying the output for the function. The implemented outputs that may be obtained are: \code{"plot"} and \code{"table"}. The argument \code{result} will be \code{"plot"} by default.
#' @param export.plot A \code{logical} value specifying if a plot saved in the working directory will be required or not. If \code{FALSE} graphical results will only be displayed in the active graphics window. If \code{TRUE} graphical results will be displayed in the active graphics window and also a \emph{pdf} or \emph{png} file (according to the \code{export.format} argument) will be saved within the \emph{plot_AeRobiology} directory automatically created in the working directory. This argument is applicable only for \code{"static"} plots. The \code{export.plot} will be \code{FALSE} by default.
#' @param export.format A \code{character} string specifying the format selected to save the plot showing the relative abundance of the pollen types. The implemented formats that may be used are: \code{"pdf"} and \code{"png"}. This argument is applicable only for \code{"static"} plots. The \code{export.format} will be \code{"pdf"} by default.
#' @param exclude A \code{character} string vector with the names of the pollen types to be excluded from the plot.
#' @param ... Other additional arguments may be used to customize the exportation of the plots using \code{"pdf"} or \code{"png"} files and therefore arguments from \code{\link[grDevices]{pdf}} and \code{\link[grDevices]{png}} functions (\pkg{grDevices} package) may be implemented. For example, for \emph{pdf} files the user may custom the arguments: \code{width}, \code{height}, \code{family}, \code{title}, \code{fonts}, \code{paper}, \code{bg}, \code{fg}, \code{pointsize...}; and for \emph{png} files the user may custom the arguments: \code{width}, \code{height}, \code{units}, \code{pointsize}, \code{bg}, \code{res...}
#' @details This function allows to calculate the relative abundance of the pollen types in the air from a database and to display a barplot with the percentage representation of the main pollen types as the graph reported by \emph{Rojo et al. (2016)}. This plot will be generated only for the specified number of the most abundant pollen types using the \code{n.types} argument by the user.\cr
#' \cr
#' Pollen time series frequently have different gaps with no data and this fact could be a problem for the calculation of specific methods for defining the pollen season even providing incorrect results. In this sense by default a linear interpolation will be carried out to complete these gaps before to define the pollen season (\code{interpolation = TRUE}). Additionally, the users may select other interpolation methods using the \code{int.method} argument, as \code{"lineal"}, \code{"movingmean"}, \code{"spline"} or \code{"tseries"}. For more information to see the \code{\link{interpollen}} function.
#' @return This function returns different results:\cr
#' \cr
#' \code{plot} in the active graphics window displaying the pollen calendar generated by the user when \code{result = "plot"}. This \code{plot} may be included in an object by assignment operators.\cr
#' \cr
#' \code{data.frame} including the yearly average pollen amounts for each pollen types used to generate the pollen of the relative abundance when \code{result = "table"}. This \code{data.frame} will be included in an object named \code{annual.sum.data}.\cr
#' \cr
#' If \code{export.plot = FALSE} graphical results will only be displayed in the active graphics window as \strong{ggplot} graph. Additional characteristics may be incorporated to the plot by \code{\link[ggplot2]{ggplot}} syntax (see \pkg{ggplot2} package).\cr
#' \cr
#' If \code{export.plot = TRUE} and \code{export.format = pdf} a \emph{pdf} file with the plot will be saved within the \emph{plot_AeRobiology} directory created in the working directory. This option is applicable only for \code{"static"} plots. Additional characteristics may be incorporated to the exportation as \emph{pdf} file (see \pkg{grDevices} package).\cr
#' \cr
#' If \code{export.plot = TRUE} and \code{export.format = png} a \emph{png} file with the plot will be saved within the \emph{plot_AeRobiology} directory created in the working directory. This option is applicable only for \code{"static"} plots. Additional characteristics may be incorporated to the exportation as \emph{png} file (see \pkg{grDevices} package).\cr
#' \cr
#' If \code{type.plot = dynamic} graphical results will be displayed in the active Viewer window as \strong{plotly} graph. Additional characteristics may be incorporated to the plot by \code{\link[plotly]{plotly}} syntax (see \pkg{plotly} package).
#' @references Rojo, J., Rapp, A., Lara, B., Sabariego, S., Fernandez_Gonzalez, F. and Perez_Badia, R., 2016. Characterisation of the airborne pollen spectrum in Guadalajara (central Spain) and estimation of the potential allergy risk. \emph{Environmental Monitoring and Assessment}, 188(3), p.130.
#' @seealso \code{\link{interpollen}}
#' @examples data("munich_pollen")
#' @examples iplot_abundance (munich_pollen, interpolation = FALSE, export.plot = FALSE)
#' @importFrom utils data
#' @importFrom lubridate is.POSIXt
#' @importFrom dplyr group_by summarise_all
#' @importFrom plotly ggplotly
#' @importFrom ggplot2 aes coord_flip element_blank element_text geom_bar geom_errorbar ggplot labs position_dodge theme theme_classic
#' @importFrom graphics plot
#' @importFrom grDevices dev.off pdf png
#' @importFrom stats sd
#' @importFrom tidyr %>%
#' @export
iplot_abundance <- function (data,
                             n.types = 15,
                             y.start = NULL,
                             y.end = NULL,
                             interpolation = TRUE,
                             int.method = "lineal",
                             col.bar = "#E69F00",
                             type.plot = "static",
                             result = "plot",
                             export.plot = FALSE,
                             export.format = "pdf",
                             exclude = NULL,
                             ...){

#############################################    CHECK THE ARGUMENTS       #############################

  if(export.plot == TRUE){ifelse(!dir.exists(file.path("plot_AeRobiology")), dir.create(file.path("plot_AeRobiology")), FALSE)}

  data<-data.frame(data)
  if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

  if(class(n.types) != "numeric") stop ("Please include only numeric values for 'n.types' argument indicating the number of pollen types which will be displayed")

  if(class(y.start) != "numeric" & !is.null(y.start)) stop ("Please include only numeric values for y.start argument indicating the start year considered")

  if(class(y.end) != "numeric" & !is.null(y.end)) stop ("Please include only numeric values for 'y.end' argument indicating the end year considered")

  if(class(interpolation) != "logical") stop ("Please include only logical values for interpolation argument")

  if(int.method != "lineal" & int.method != "movingmean" & int.method != "spline" & int.method != "tseries") stop ("Please int.method only accept values: 'lineal', 'movingmean', 'spline' or 'tseries'")

  if(class(col.bar) != "character") stop ("Please include only character values indicating the color selected in the bars for generating the plot")

  if(type.plot != "static" & type.plot != "dynamic") stop ("Please type.plot only accept values: 'static' or 'dynamic'")

  if(result != "plot" & result != "table") stop ("Please result only accept values: 'plot' or 'table'")

  if(class(export.plot) != "logical") stop ("Please include only logical values for export.plot argument")

  if(export.format != "pdf" & export.format != "png") stop ("Please export.format only accept values: 'pdf' or 'png'")

  if(class(exclude) != "character" & !is.null(exclude)) stop ("Please include only character values for exclude argument indicating the pollen type to be excluded")

  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])

if(interpolation == TRUE){data <- interpollen(data, method = int.method, plot = FALSE)}

#############################################    SELECT ABUNDANT TYPES       #############################
  #annual.sum.data<-data.frame()
data <- data.frame(date = data[ ,1], year = as.numeric(strftime(data[ ,1], "%Y")), data[ ,-1])

#types <- ddply(data[ ,-1], "year", function(x) colSums(x[-1], na.rm = T)) [-1] %>%
#  apply(., 2, function(x) mean(x, na.rm = T)) %>% .[order(., decreasing = TRUE)] %>%
#  names(.) %>%
#  .[1:n.types]
types <- data.frame(data[ ,-1] %>%
                      group_by(year) %>%
                      summarise_all(sum, na.rm = TRUE))[-1] %>%
  apply(., 2, function(x) mean(x, na.rm = T)) %>% .[order(., decreasing = TRUE)] %>%
  names(.) %>%
  .[1:n.types]

#data <- data[ ,which(colnames(data) %in% c("date", "year", types))]

#############################################    SELECT PERIOD       #############################

seasons <- unique(as.numeric(strftime(data[, 1], "%Y")))

if(is.null(y.start)){y.start <- min(seasons)}; if(is.null(y.end)){y.end <- max(seasons)}

data <- data[which(as.numeric(strftime(data[ ,1], "%Y")) >= y.start & as.numeric(strftime(data[ ,1], "%Y")) <= y.end), ]

##############################################################################

#sum.data <- ddply(data, c("year"), function(x) colSums(x[ ,-which(colnames(x) %in% c("date", "year"))], na.rm = T))
sum.data <- data.frame(data[ ,-1] %>%
                      group_by(year) %>%
                      summarise_all(sum, na.rm = TRUE))

sum.data$Total <- apply(sum.data[ ,-which(colnames(sum.data) %in% c("year"))], 1, sum, na.rm = T)

perc.df <- sum.data[ ,-which(colnames(sum.data) == "year")]
for (r in 1:nrow(perc.df)){
  for(c in c(1:ncol(perc.df))){
    perc.df[r,c] <- perc.df[r,c]*100/perc.df$Total[r]
  }
}

perc.df <- perc.df[which(colnames(perc.df) %in% types)]

mean.perc <- data.frame(types = colnames(perc.df),
        mean = apply(perc.df, 2, FUN = function(x) mean(x, na.rm=T)),
        sd = apply(perc.df, 2, FUN = function(x) sd(x, na.rm=T)))

mean.perc$types <- factor(mean.perc$types, levels = mean.perc$types[order(mean.perc$mean)])

if (!is.null(exclude)){
  mean.perc<-mean.perc[!(as.character(mean.perc$types)%in%exclude),]}

plot.abundance <- ggplot(mean.perc, aes(x = types, y = mean)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", fill = col.bar) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.2,position=position_dodge(.9))+
  coord_flip()+
  theme_classic()+
  labs(y="Relative abundance (%)", title = paste0("Relative Abundance in the Air (", y.start, "-", y.end, ")"), size = 14)+
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=10, face="bold"), axis.text.y = element_text(size=10, face="bold.italic"), axis.text = element_text(size=10), title=element_text(size=10, face="bold"), plot.title = element_text(hjust = 0.5, size = 16))

if(export.plot == TRUE & type.plot == "static" & export.format == "png") {
  png(paste0("plot_AeRobiology/abundance_plot", y.start, "-", y.end, ".png"), ...)
    plot(plot.abundance)
  dev.off()
  png(paste0("plot_AeRobiology/credits.png"))

  dev.off()
}

if(export.plot == TRUE & type.plot == "static" & export.format == "pdf") {
  pdf(paste0("plot_AeRobiology/abundance_plot", y.start, "-", y.end, ".pdf"), ...)
    plot(plot.abundance)

  dev.off()
}

if (result == "plot" & type.plot == "static") {return(plot.abundance)}
if (result == "plot" & type.plot == "dynamic") {return(ggplotly(plot.abundance))}

if (result == "table") {return(sum.data)}

}
