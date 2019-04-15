#' Phenological Plot
#'
#' Generates a boxplot based on phenological parameters (start_dates and end_dates) that are calculated by the estimation of the main parameters of the pollen season
#'
#' @param data A \code{data.frame} object including the general database where calculation of the pollen season must be applied in order to generate the phenological plot based on the start_dates and end_dates. This \code{data.frame} must include a first column in \code{Date} format and the rest of columns in \code{numeric} format belonging to each pollen type by column.
#' @param method A \code{character} string specifying the method applied to calculate the pollen season and the main parameters. The implemented methods that can be used are: \code{"percentage"}, \code{"logistic"}, \code{"moving"}, \code{"clinical"} or \code{"grains"}. A more detailed information about the different methods for defining the pollen season may be consulted in \code{\link{calculate_ps}} function.
#' @param n.types A \code{numeric} (\code{integer}) value specifying the number of the most abundant pollen types that must be represented in the pollen calendar. A more detailed information about the selection of the considered pollen types may be consulted in \strong{Details}. The \code{n.types} argument will be \code{15} types by default.
#' @param th.day A \code{numeric} value in order to calculate the number of days when this level is exceeded for each year and each pollen type. This value will be obtained in the results of the function. The \code{th.day} argument will be \code{100} by default.
#' @param perc A \code{numeric} value ranging \code{0_100}. This argument is valid only for \code{method = "percentage"}. This value represents the percentage of the total annual pollen included in the pollen season, removing \code{(100_perc)/2\%} of the total pollen before and after of the pollen season. The \code{perc} argument will be \code{95} by default.
#' @param def.season A \code{character} string specifying the method for selecting the best annual period to calculate the pollen season. The pollen seasons may occur within the natural year or otherwise may occur between two years which determines the best annual period considered. The implemented options that can be used are: \code{"natural"}, \code{"interannual"} or \code{"peak"}. The \code{def.season} argument will be \code{"natural"} by default. A more detailed information about the different methods for selecting the best annual period to calculate the pollen season may be consulted in \code{\link{calculate_ps}} function.
#' @param reduction A \code{logical} value. This argument is valid only for the \code{"logistic"} method. If \code{FALSE} the reduction of the pollen data is not applicable. If \code{TRUE} a reduction of the peaks above a certain level (\code{red.level} argument) will be carried out before the definition of the pollen season. The \code{reduction} argument will be \code{FALSE} by default. A more detailed information about the reduction process may be consulted in \code{\link{calculate_ps}} function.
#' @param red.level A \code{numeric} value ranging \code{0_1} specifying the percentile used as level to reduce the peaks of the pollen series before the definition of the pollen season. This argument is valid only for the \code{"logistic"} method. The \code{red.level} argument will be \code{0.90} by default, specifying the percentile 90.
#' @param derivative A \code{numeric} (\code{integer}) value belonging to options of \code{4}, \code{5} or \code{6} specifying the derivative that will be applied to calculate the asymptotes which determines the pollen season using the \code{"logistic"} method. This argument is valid only for the \code{"logistic"} method. The \code{derivative} argument will be \code{5} by default.
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
#' @param type.plot A \code{character} string specifying the type of plot selected to show the phenological plot. The implemented types that may be used are: \code{"static"} generates a static \strong{ggplot} object and \code{"dynamic"} generates a dynamic \strong{plotly} object.
#' @param export.plot A \code{logical} value specifying if a phenological plot saved in the working directory will be required or not. If \code{FALSE} graphical results will only be displayed in the active graphics window. If \code{TRUE} graphical results will be displayed in the active graphics window and also a \emph{pdf} or \emph{png} file (according to the \code{export.format} argument) will be saved within the \emph{plot_AeRobiology} directory automatically created in the working directory. This argument is applicable only for \code{"static"} plots. The \code{export.plot} will be \code{FALSE} by default.
#' @param export.format A \code{character} string specifying the format selected to save the phenological plot. The implemented formats that may be used are: \code{"pdf"} and \code{"png"}. This argument is applicable only for \code{"static"} plots. The \code{export.format} will be \code{"pdf"} by default.
#' @param ... Other additional arguments may be used to customize the exportation of the plots using \emph{pdf} or \emph{png} files and therefore arguments from \code{\link[grDevices]{pdf}} and \code{\link[grDevices]{png}} functions (\pkg{grDevices} package) may be implemented. For example, for \emph{pdf} files the user may custom the arguments: \code{width}, \code{height}, \code{family}, \code{title}, \code{fonts}, \code{paper}, \code{bg}, \code{fg}, \code{pointsize...}; and for \emph{png} files the user may custom the arguments: \code{width}, \code{height}, \code{units}, \code{pointsize}, \code{bg}, \code{res...}
#'
#' @details This function allows to calculate the pollen season using five different methods which are described in \code{\link{calculate_ps}} function. After calculating the start_date and end_date for each pollen type and each year a phenological plot will be generated using the boxplot approach where axis x represents the time (Day of the Year) and axis y includes the considered pollen types. The phenological plot will be generated only for the specified number of the most abundant pollen types using the \code{n.types} argument by the user. The implemented methods for defining the pollen season includes the most commonly used methodologies (\emph{Nilsson and Persson, 1981}; \emph{Andersen, 1991}; \emph{Galan et al., 2001}; \emph{Ribeiro et al., 2007}; \emph{Cunha et al., 2015}, \emph{Pfaar et al., 2017}) and a new implemented method (see \code{\link{calculate_ps}} function).\cr
#' Pollen time series frequently have different gaps with no data and this fact could be a problem for the calculation of specific methods for defining the pollen season even providing incorrect results. In this sense by default a linear interpolation will be carried out to complete these gaps before to define the pollen season (\code{interpolation = TRUE}). Additionally, the users may select other interpolation methods using the \code{int.method} argument, as \code{"lineal"}, \code{"movingmean"}, \code{"spline"} or \code{"tseries"}. For more information to see the \code{\link{interpollen}} function.
#'
#' @return This function returns different results:\cr
#' If \code{export.plot = FALSE} graphical results will only be displayed in the active graphics window as \strong{ggplot} graph. Additional characteristics may be incorporated to the plot by \code{\link[ggplot2]{ggplot}} syntax (see \pkg{ggplot2} package).\cr
#' If \code{export.plot = TRUE} and \code{export.format = pdf} a \emph{pdf} file of the phenological plot will be saved within the \emph{plot_AeRobiology} directory created in the working directory. This option is applicable only for \code{"static"} plots. Additional characteristics may be incorporated to the exportation as \emph{pdf} file (see \pkg{grDevices} package).\cr
#' If \code{export.plot = TRUE} and \code{export.format = png} a \emph{png} file of the phenological plot will be saved within the \emph{plot_AeRobiology} directory created in the working directory. This option is applicable only for \code{"static"} plots. Additional characteristics may be incorporated to the exportation \emph{png} file (see \pkg{grDevices} package).\cr
#' If \code{type.plot = dynamic} graphical results will be displayed in the active Viewer window as \strong{plotly} graph. Additional characteristics may be incorporated to the plot \code{\link[plotly]{plotly}} syntax (see \pkg{plotly} package).
#' @references Andersen, T.B., 1991. A model to predict the beginning of the pollen season. \emph{Grana}, 30(1), pp.269_275.
#' @references Cunha, M., Ribeiro, H., Costa, P. and Abreu, I., 2015. A comparative study of vineyard phenology and pollen metrics extracted from airborne pollen time series. \emph{Aerobiologia}, 31(1), pp.45_56.
#' @references  Galan, C., Garcia_Mozo, H., Carinanos, P., Alcazar, P. and Dominguez_Vilches, E., 2001. The role of temperature in the onset of the \emph{Olea europaea} L. pollen season in southwestern Spain. \emph{International Journal of Biometeorology}, 45(1), pp.8_12.
#' @references Nilsson, S. and Persson, S., 1981. Tree pollen spectra in the Stockholm region (Sweden), 1973_1980. \emph{Grana}, 20(3), pp.179_182.
#' @references Pfaar, O., Bastl, K., Berger, U., Buters, J., Calderon, M.A., Clot, B., Darsow, U., Demoly, P., Durham, S.R., Galan, C., Gehrig, R., Gerth van Wijk, R., Jacobsen, L., Klimek, L., Sofiev, M., Thibaudon, M. and Bergmann, K.C., 2017. Defining pollen exposure times for clinical trials of allergen immunotherapy for pollen_induced rhinoconjunctivitis_an EAACI position paper. \emph{Allergy}, 72(5), pp.713_722.
#' @references Ribeiro, H., Cunha, M. and Abreu, I., 2007. Definition of main pollen season using logistic model. \emph{Annals of Agricultural and Environmental Medicine}, 14(2), pp.259_264.
#' @seealso \code{\link{calculate_ps}}, \code{\link{interpollen}}
#' @examples data("munich_pollen")
#' @examples iplot_pheno (munich_pollen, interpolation = FALSE)
#' @importFrom plotly layout
#' @importFrom ggplot2 aes coord_flip element_text geom_boxplot ggplot labs position_dodge scale_fill_manual scale_y_continuous theme theme_bw
#' @importFrom graphics plot
#' @importFrom grDevices dev.off pdf png
#' @importFrom lubridate is.POSIXt
#' @importFrom plotly ggplotly
#' @importFrom dplyr group_by summarise_all
#' @importFrom stats na.omit
#' @importFrom utils data
#' @importFrom tidyr %>%
#' @export
iplot_pheno <- function(data,
                        method = "percentage",
                        n.types = 15,
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
                        type.plot = "static",
                        export.plot = FALSE,
                        export.format = "pdf",...){

#############################################    CHECK THE ARGUMENTS       #############################

  if(export.plot == TRUE){ifelse(!dir.exists(file.path("plot_AeRobiology")), dir.create(file.path("plot_AeRobiology")), FALSE)}

  data<-data.frame(data)
  if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

  if(method != "percentage" & method != "logistic" & method != "moving" & method != "clinical" & method != "grains") stop ("Please method only accept values: 'percentage', 'logistic', 'moving', 'clinical' or 'grains'")

  if(class(n.types) != "numeric") stop ("Please include only numeric values for 'n.types' argument indicating the number of pollen types which will be displayed")

  if(class(th.day) != "numeric" | th.day < 0) stop ("Please include only numeric values >= 0 for th.day argument")

  if(class(perc) != "numeric" | perc < 0 | perc > 100) stop ("Please include only numeric values between 0-100 for perc argument")

  if(def.season != "natural" & def.season != "interannual" & def.season != "peak") stop ("Please def.season only accept values: 'natural', 'interannual' or 'peak'")

  if(class(reduction) != "logical") stop ("Please include only logical values for reduction argument")

  if(class(red.level) != "numeric" | red.level < 0 | red.level > 1) stop ("Please include only numeric values between 0-1 for red.level argument")

  if(derivative != 4 & derivative != 5 & derivative != 6) stop ("Please derivative only accept values: 4, 5 or 6")

  if(class(man) != "numeric" | man < 0) stop ("Please include only numeric values > 0 for man argument")

  if(class(th.ma) != "numeric" | th.ma < 0) stop ("Please include only numeric values > 0 for th.ma argument")

  if(class(n.clinical) != "numeric" | n.clinical < 0) stop ("Please include only numeric values >= 0 for n.clinical argument")

  if(class(window.clinical) != "numeric" | window.clinical < 0) stop ("Please include only numeric values >= 0 for window.clinical argument")

  if(class(window.grains) != "numeric" | window.grains < 0) stop ("Please include only numeric values >= 0 for window.grains argument")

  if(class(th.pollen) != "numeric" | th.pollen < 0) stop ("Please include only numeric values >= 0 for th.pollen argument")

  if(class(th.sum) != "numeric" | th.sum < 0) stop ("Please include only numeric values >= 0 for th.sum argument")

  if(type != "none" & type != "birch" & type != "grasses" & type != "cypress" & type != "olive" & type != "ragweed") stop ("Please def.season only accept values: 'none', 'birch', 'grasses', 'cypress', 'olive' or 'ragweed'")

  if(class(interpolation) != "logical") stop ("Please include only logical values for interpolation argument")

  if(int.method != "lineal" & int.method != "movingmean" & int.method != "spline" & int.method != "tseries") stop ("Please int.method only accept values: 'lineal', 'movingmean', 'spline' or 'tseries'")

  if(type.plot != "static" & type.plot != "dynamic") stop ("Please type.plot only accept values: 'static' or 'dynamic'")

  if(class(export.plot) != "logical") stop ("Please include only logical values for export.plot argument")

  if(export.format != "pdf" & export.format != "png") stop ("Please export.format only accept values: 'pdf' or 'png'")

  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])

options(warn = -1)

#############################################    SELECT ABUNDANT TYPES       #############################

data <- data.frame(date = data[ ,1], year = as.numeric(strftime(data[ ,1], "%Y")), data[ ,-1])

#types <- ddply(data[ ,-1], "year", function(x) colSums(x[-1], na.rm = T)) [-1] %>%
  #apply(., 2, function(x) mean(x, na.rm = T)) %>% .[order(., decreasing = TRUE)] %>%
  #names(.) %>%
  #.[1:n.types]
types <- data.frame(data[ ,-1] %>%
                      group_by(year) %>%
                      summarise_all(sum, na.rm = TRUE))[-1] %>%
  apply(., 2, function(x) mean(x, na.rm = T)) %>% .[order(., decreasing = TRUE)] %>%
  names(.) %>%
  .[1:n.types]

data <- data[ ,which(colnames(data) %in% c("date", types))]

#############################################    CALCULATION OF THE POLLEN SEASON       #############################

list.ps <- calculate_ps(data = data, method = method, th.day = th.day, perc = perc, def.season = def.season, reduction = reduction, red.level = red.level, derivative = derivative, man = man, th.ma = th.ma, n.clinical = n.clinical, window.clinical = window.clinical, window.grains = window.grains, th.pollen = th.pollen, th.sum = th.sum, type = type, interpolation = interpolation, int.method = int.method, result = "list", plot = FALSE, export.plot = FALSE, export.result = FALSE)

season.data <- list.ps
type.name <- colnames(data)[-1]

dat.plot <- data.frame(); mean.st <- data.frame()

for(t in 1:(length(type.name))){
  dat.type <- season.data[[type.name[t]]]

  dat.plot <- rbind(dat.plot, data.frame(phen = as.numeric(c(rep(0, nrow(dat.type)), rep(1, nrow(dat.type)))), j.days = c(dat.type$st.jd, dat.type$en.jd), type = type.name[t])) # %>%  na.omit(.)

  mean.st <- rbind(mean.st, data.frame(type = type.name[t], mean.st = mean(dat.type$st.jd, na.rm = TRUE)))

}

  mean.st$mean.st[is.nan(mean.st$mean.st)] <- NA; mean.st <- na.omit(mean.st)#; mean.st <- mean.st[order(mean.st$mean.st), 1]

  dat.plot$type <- factor(dat.plot$type, levels = as.character(mean.st[order(mean.st$mean.st), 1]), ordered = T)

  dat.plot$phen[dat.plot$phen == 0] <- "Start-date"; dat.plot$phen[dat.plot$phen == 1] <- "End-date"

  dat.plot$phen <- factor(dat.plot$phen, levels = c("Start-date", "End-date"), ordered = T)

pheno.plot <- ggplot(data = na.omit(dat.plot), aes(y = j.days, x = type)) +
  #stat_boxplot(geom = "errorbar", width = 0.2)+
  geom_boxplot(aes(fill = phen), width = 1, varwidth = TRUE, position = position_dodge(width = 0))+
  labs(y = "Day of the year", x = "", title = "Phenological parameters", fill = "Phenophase")+
  scale_fill_manual(labels = c("Start-date", "End-date"), values = c("tan1", "lightskyblue")) +
  scale_y_continuous(breaks = seq(0,365,25))+
  coord_flip()+
  theme_bw()+
  theme(text = element_text(size = 14), legend.position = "bottom", axis.text.y = element_text(face="italic"))

if(export.plot == TRUE & type.plot == "static" & export.format == "png") {
  if(method == "percentage") {png(paste0("plot_AeRobiology/pheno_plot_",method,perc,".png"), ...); plot(pheno.plot); dev.off(); png(paste0("plot_AeRobiology/credits.png")); dev.off()}
  if(method == "logistic") {png(paste0("plot_AeRobiology/pheno_plot_",method,derivative,".png"), ...); plot(pheno.plot); dev.off(); png(paste0("plot_AeRobiology/credits.png")); dev.off()}
  if(method == "moving") {png(paste0("plot_AeRobiology/pheno_plot_",method,man,"_",th.ma,".png"), ...); plot(pheno.plot); dev.off(); png(paste0("plot_AeRobiology/credits.png")); dev.off()}
  if(method == "clinical") {png(paste0("plot_AeRobiology/pheno_plot_",method,n.clinical,"_",window.clinical+1,"_",th.pollen,"_",th.sum,".png"), ...); plot(pheno.plot); dev.off(); png(paste0("plot_AeRobiology/credits.png")); dev.off()}
  if(method == "grains") {png(paste0("plot_AeRobiology/pheno_plot_",method,window.grains+1,"_",th.pollen,".png"), ...); plot(pheno.plot); dev.off(); png(paste0("plot_AeRobiology/credits.png")); dev.off()}
}

if(export.plot == TRUE & type.plot == "static" & export.format == "pdf") {
  if(method == "percentage") {pdf(paste0("plot_AeRobiology/pheno_plot_",method,perc,".pdf"), ...); plot(pheno.plot); dev.off()}
  if(method == "logistic") {pdf(paste0("plot_AeRobiology/pheno_plot_",method,derivative,".pdf"), ...); plot(pheno.plot); dev.off()}
  if(method == "moving") {pdf(paste0("plot_AeRobiology/pheno_plot_",method,man,"_",th.ma,".pdf"), ...); plot(pheno.plot); dev.off()}
  if(method == "clinical") {pdf(paste0("plot_AeRobiology/pheno_plot_",method,n.clinical,"_",window.clinical+1,"_",th.pollen,"_",th.sum,".pdf"), ...); plot(pheno.plot); dev.off()}
  if(method == "grains") {pdf(paste0("plot_AeRobiology/pheno_plot_",method,window.grains+1,"_",th.pollen,".pdf"), ...); plot(pheno.plot); dev.off()}
}

if (type.plot == "static") {return(pheno.plot)}

options(warn = 0)

if (type.plot == "dynamic") {
  ggplotly(pheno.plot) %>%
  layout(legend = list(orientation = 'v', x = -0.2, y = -0.3))
  }
}
