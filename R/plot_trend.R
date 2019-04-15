#' Calculating and Plotting Trends of Pollen Data.
#'
#' Function to calculate the main seasonal indexes of the pollen season (\emph{Start Date}, \emph{Peak Date}, \emph{End Date} and \emph{Pollen Integral}). Trends analysis of the parameters over the seasons. Plots showing the distribution of the main seasonal indexes over the years.
#'
#' @param data A \code{data.frame} object. This \code{data.frame} should include a first column in format \code{Date} and the rest of columns in format \code{numeric} belonging to each pollen type by column.
#' @param interpolation A \code{logical} value specifying if the visualization shows the gaps in the inputs data (\code{interpolation = FALSE}) or if an interpolation method is used for filling the gaps (\code{interpolation = TRUE}). By default, \code{interpolation = TRUE}.
#' @param int.method A \code{character} string with the name of the interpolation method to be used. The implemented methods that may be used are: \code{"lineal"}, \code{"movingmean"}, \code{"tseries"} or \code{"spline"}. By default, \code{int.method = "lineal"}.
#' @param export.plot A \code{logical} value specifying if a plot will be exported or not. If \code{FALSE} graphical results will only be displayed in the active graphics window. If \code{TRUE} graphical results will be displayed in the active graphics window and also one pdf/png file will be saved within the \emph{plot_AeRobiology} directory automatically created in the working directory. By default, \code{export.plot = TRUE}.
#' @param export.format A \code{character} string specifying the format selected to save the plot. The implemented formats that may be used are: \code{"pdf"} or \code{"png"}. By default, \code{export.format = "pdf"}.
#' @param export.result A \code{logical} value. If \code{export.result = TRUE}, a table is exported with the extension \emph{.xlsx}, in the directory \emph{table_AeRobiology}. This table has the information about the \code{slope} \emph{"beta coefficient of a lineal model using as predictor the year and as dependent variable one of the main pollen season indexes"}. The information is referred to the main pollen season indexes: \emph{Start Date}, \emph{Peak Date}, \emph{End Date} and \emph{Pollen Integral}.
#' @param method A \code{character} string specifying the method applied to calculate the pollen season and the main seasonal parameters. The implemented methods that can be used are: \code{"percentage"}, \code{"logistic"}, \code{"moving"}, \code{"clinical"} or \code{"grains"}. By default, \code{method = "percentage"} (\code{perc = 95}\%). A more detailed information about the different methods for defining the pollen season may be consulted in the function \code{\link{calculate_ps}}.
#' @param ... Additional arguments for the function \code{\link{calculate_ps}} are also accepted.
#' @return This function returns several plots in the directory \emph{plot_AeRobiology/trend_plots} with the extension \emph{.pdf} or \emph{.png}.Also produces an object of the class \code{data.frame} and export a table with the extension \emph{.xlsx}, in the directory \emph{table_AeRobiology}.\cr
#'These tables have the information about the \code{slope} \emph{(beta coefficient of a lineal model using as predictor the year and as dependent variable one of the main pollen season indexes)}. The information is referred to the main pollen season indexes: \emph{Start Date}, \emph{Peak Date}, \emph{End Date} and \emph{Pollen Integral}.
#' @seealso \code{\link{calculate_ps}}; \code{\link{analyse_trend}}
#' @examples  data("munich_pollen")
#' @examples  plot_trend(munich_pollen, interpolation = FALSE, export.plot = FALSE, export.result = TRUE)
#' @importFrom graphics plot
#' @importFrom utils data
#' @importFrom ggplot2 aes element_text geom_point geom_smooth ggplot labs theme theme_classic theme_set scale_x_continuous
#' @importFrom grDevices dev.off pdf png
#' @importFrom grid grid.layout pushViewport viewport
#' @importFrom lubridate is.POSIXt
#' @importFrom stats as.formula complete.cases lm pf
#' @importFrom scales pretty_breaks
#' @importFrom writexl write_xlsx
#' @importFrom tidyr %>%
#' @export
plot_trend  <-   function   (data,
                            interpolation = TRUE,
                            int.method = "lineal",
                            export.plot = TRUE,
                            export.format = "pdf",
                            export.result=TRUE,
                            method="percentage",
                            ...){

  # Sys.setlocale(category = "LC_ALL", locale="english")

  #############################################    CHECK THE ARGUMENTS       #############################

  if(export.plot == TRUE){ifelse(!dir.exists(file.path("plot_AeRobiology")), dir.create(file.path("plot_AeRobiology")), FALSE)}
  if(export.plot == TRUE){ifelse(!dir.exists(file.path("plot_AeRobiology/trend_plots")), dir.create(file.path("plot_AeRobiology/trend_plots")), FALSE)}
  if(export.result == TRUE){ifelse(!dir.exists(file.path("table_AeRobiology")), dir.create(file.path("table_AeRobiology")), FALSE)}
  data<-data.frame(data)
  if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

  if(class(export.plot) != "logical") stop ("Please include only logical values for export.plot argument")

  if(export.format != "pdf" & export.format != "png") stop ("Please export.format only accept values: 'pdf' or 'png'")

  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])

  if(class(interpolation) != "logical") stop ("Please include only logical values for interpolation argument")

  # if(int.method != "lineal" & int.method != "movingmean" & int.method != "spline" & int.method != "tseries") stop ("Please int.method only accept values: 'lineal', 'movingmean', 'tseries' or 'spline'")

  # if(interpolation == TRUE){data <- interpollen(data, method = int.method)}


  colnames(data)[1]<-"Date"

  ## Function for p value
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }

  datafram<-calculate_ps(data, method=method, interpolation = interpolation, int.method=int.method,plot=FALSE,...)

  variables<-c("st.jd","pk.jd","en.jd","sm.ps")

  trendtime<-data.frame()

  data_summary<-
    for (t in 1:length(unique(datafram$type))){

      type<-unique(as.character(datafram$type))[t]



          for (v in 1:length(variables)){
            tryCatch({
            variable<-variables[v]

            temp<-datafram[which(datafram$type==type),c(1:2,which( colnames(datafram)==variable))]

            lm <- lm (as.formula(paste(variable,"~ seasons")), data= temp, x = TRUE, y = TRUE)

            tempframe<-data.frame(type=type,variable=variable,coef=summary(lm)$coefficients[2,1],p=lmp(lm))

            trendtime<-rbind(trendtime,tempframe)
            }, error=function(e){
              print(paste(type, variable, ": Error, linear model not calculated. Probably due to insufficient amount of years"))
            })
          }
      }




datafram<-datafram[complete.cases(datafram),]


### Now start the rock and roll

  for (p in 1:length(unique(as.character(datafram$type)))){
      pollen<-unique(as.character(datafram$type))[p]
      dataframtemp<-datafram[which(datafram$type==pollen),]

      dataframtemp$seasons<-as.integer(dataframtemp$seasons)

  slope<-trendtime[which(trendtime$type==pollen & trendtime$variable=="st.jd"),"coef"]
  slope<-round(slope, 1)
  slope<-paste0("slope: ",slope)
  p<-trendtime[which(trendtime$type==pollen & trendtime$variable=="st.jd"),"p"]
  p<-round(p, 3)
  if (p < 0.001) {
    pvalue <- "p<0.001"
  } else if (p < 0.01) {
    pvalue <- "p<0.01"
  } else if (p < 0.05) {
    pvalue <- "p<0.05"
  } else {
    pvalue <- "p>0.05"
  }
  p<-pvalue
  comb<-paste0(slope,", ",p)

  p1 <- ggplot(dataframtemp, aes(x=seasons, y=st.jd)) +
    theme_set(theme_classic()) +
    geom_smooth(method="loess", size=1, colour="blue", fill="light blue") +
    geom_smooth(method="lm", size=1, se=FALSE, colour="red", linetype="dashed")  +
    geom_point(shape=21, stroke=1.5,size = 3.5, fill="darkgrey", colour="black") +
    labs(x="", y="Day of the Year (DOY)", title= paste("Start Date", pollen), subtitle=comb)+
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10)) +
    theme(plot.title=element_text(size=13, face = "bold.italic"))+
    scale_x_continuous(breaks=pretty_breaks())

  p1 <- ggplot(dataframtemp, aes(x=seasons, y=st.jd)) +
    theme_set(theme_classic()) +
    geom_smooth(method="loess", size=1, colour="blue", fill="light blue") +
    geom_smooth(method="lm", size=1, se=FALSE, colour="red", linetype="dashed")  +
    geom_point(shape=21, stroke=1.5,size = 3.5, fill="darkgrey", colour="black") +
    labs(x="", y="Day of the Year (DOY)", title= paste("Start Date", pollen), subtitle=comb)+
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10)) +
    theme(plot.title=element_text(size=13, face = "bold.italic"))+
    scale_x_continuous(breaks=pretty_breaks())


  slope<-trendtime[which(trendtime$type==pollen & trendtime$variable=="pk.jd"),"coef"]
  slope<-round(slope, 1)
  slope<-paste0("slope: ",slope)
  p<-trendtime[which(trendtime$type==pollen & trendtime$variable=="pk.jd"),"p"]
  p<-round(p, 3)
  if (p < 0.001) {
    pvalue <- "p<0.001"
  } else if (p < 0.01) {
    pvalue <- "p<0.01"
  } else if (p < 0.05) {
    pvalue <- "p<0.05"
  } else {
    pvalue <- "p>0.05"
  }
  p<-pvalue
  comb<-paste0(slope,", ",p)
  p2 <- ggplot(dataframtemp, aes(x=seasons, y=pk.jd)) +
    theme_set(theme_classic()) +
    geom_smooth(method="loess",size=1, colour="blue", fill="light blue") +
    geom_smooth(method="lm", size=1, se=FALSE, colour="red", linetype="dashed")  +
    geom_point(shape=21, stroke=1.5,size = 3.5,  fill="darkgrey", colour="black") +
    labs(x="", y="Day of the Year (DOY)", title=paste("Peak Date", pollen), subtitle=comb)+
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10)) +
    theme(plot.title=element_text(size=13, face = "bold.italic"))+
    scale_x_continuous(breaks=pretty_breaks())

  slope<-trendtime[which(trendtime$type==pollen & trendtime$variable=="en.jd"),"coef"]
  slope<-round(slope, 1)
  slope<-paste0("slope: ",slope)
  p<-trendtime[which(trendtime$type==pollen & trendtime$variable=="en.jd"),"p"]
  p<-round(p, 3)
  if (p < 0.001) {
    pvalue <- "p<0.001"
  } else if (p < 0.01) {
    pvalue <- "p<0.01"
  } else if (p < 0.05) {
    pvalue <- "p<0.05"
  } else {
    pvalue <- "p>0.05"
  }
  p<-pvalue
  comb<-paste0(slope,", ",p)
  p3 <- ggplot(dataframtemp, aes(x=seasons, y=en.jd)) +
    theme_set(theme_classic()) +
    geom_smooth(method="loess",size=1, colour="blue", fill="light blue") +
    geom_smooth(method="lm", size=1, se=FALSE, colour="red", linetype="dashed")  +
    geom_point(shape=21, stroke=1.5,size = 3.5,  fill="darkgrey", colour="black") +
    labs(x="", y="Day of the Year (DOY)", title= paste("End Date", pollen), subtitle=comb)+
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10)) +
    theme(plot.title=element_text(size=13, face = "bold.italic"))+
    scale_x_continuous(breaks=pretty_breaks())

  slope<-trendtime[which(trendtime$type==pollen & trendtime$variable=="sm.ps"),"coef"]
  slope<-round(slope, 1)
  slope<-paste0("slope: ",slope)
  p<-trendtime[which(trendtime$type==pollen & trendtime$variable=="sm.ps"),"p"]
  p<-round(p, 3)
  if (p < 0.001) {
    pvalue <- "p<0.001"
  } else if (p < 0.01) {
    pvalue <- "p<0.01"
  } else if (p < 0.05) {
    pvalue <- "p<0.05"
  } else {
    pvalue <- "p>0.05"
  }
  p<-pvalue
  comb<-paste0(slope,", ",p)
  p4 <- ggplot(dataframtemp, aes(x=seasons, y=sm.ps)) +
    theme_set(theme_classic()) +
    geom_smooth(method="loess", size=1, colour="blue", fill="light blue") +
    geom_smooth(method="lm", size=1, se=FALSE, colour="red", linetype="dashed")  +
    geom_point(shape=21, stroke=1.5,size = 3.5,  fill="darkgrey", colour="black") +
    labs(x="", y="Pollen grains", title=paste("Total Pollen", pollen),subtitle=comb)+
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10)) +
    theme(plot.title=element_text(size=13, face = "bold.italic"))+
    scale_x_continuous(breaks=pretty_breaks())


 if(export.plot == TRUE  & export.format == "png") {
   png(paste0("plot_AeRobiology/trend_plots/plot_trend_",pollen,".png"))
   pushViewport(viewport(layout=grid.layout(2,2)))
   vplayout<-function(x,y) viewport(layout.pos.row = x, layout.pos.col=y)
   print(p1, vp = vplayout(1,1))
   print(p2, vp = vplayout(1,2))
   print(p3, vp = vplayout(2,1))
   print(p4, vp = vplayout(2,2))
   dev.off()

 }

 if(export.plot == TRUE  & export.format == "pdf") {
   pdf(paste0("plot_AeRobiology/trend_plots/plot_trend_",pollen, ".pdf"))
   pushViewport(viewport(layout=grid.layout(2,2)))
   vplayout<-function(x,y) viewport(layout.pos.row = x, layout.pos.col=y)
   print(p1, vp = vplayout(1,1))
   print(p2, vp = vplayout(1,2))
   print(p3, vp = vplayout(2,1))
   print(p4, vp = vplayout(2,2))

   dev.off()
 }

}

lista<-list()
lista[["plot_trend"]]<-trendtime


lista [["Information"]] <- data.frame(
  Attributes = c("st.jd", "pk.jd", "en.jd", "sm.ps", "coef", "p", "", "", "Package", "Authors"),
  Description = c("Start-date (day of the year)","Peak-date (day of year)",  "End-date (day of the year)", "Pollen integral", "Slope of the linear trend", "Significance level of the linear trend",  "", "", "AeRobiology", "Jesus Rojo, Antonio Picornell & Jose Oteros"))





if (export.result == TRUE) {
  write_xlsx(lista, "table_AeRobiology/summary_of_plot_trend.xlsx")
}
return(trendtime)
}
