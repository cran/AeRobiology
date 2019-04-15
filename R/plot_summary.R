#' Plotting Several Pollen Seasons.
#'
#' Function to plot the pollen data during several seasons. Also plots the averaged pollen season over the study period. It is possible to plot the relative abundance per day and smoothing the pollen season by calculating a moving average.
#'
#' @param data A \code{data.frame} object. This \code{data.frame} should include a first column in format \code{Date} and the rest of columns in format \code{numeric} belonging to each pollen type by column.
#' @param pollen A \code{character} string with the name of the particle to show. This \code{character} must match with the name of a column in the input database. This is a mandatory argument.
#' @param mave An \code{integer} value specifying the order of the moving average applied to the data. By default, \code{mave = 1}.
#' @param normalized A \code{logical} value specifying if the visualization shows real pollen data (\code{normalized = FALSE}) or the percentage of every day over the whole pollen season (\code{normalized = TRUE}). By default, \code{normalized = FALSE}.
#' @param interpolation A \code{logical} value specifying if the visualization shows the gaps in the inputs data (\code{interpolation = FALSE}) or if an interpolation method is used for filling the gaps (\code{interpolation = TRUE}). By default, \code{interpolation = TRUE}.
#' @param int.method A \code{character} string with the name of the interpolation method to be used. The implemented methods that may be used are: \code{"lineal"}, \code{"movingmean"}, \code{"tseries"} or \code{"spline"}. By default, \code{int.method = "lineal"}.
#' @param export.plot A \code{logical} value specifying if a plot will be exported or not. If \code{FALSE} graphical results will only be displayed in the active graphics window. If \code{TRUE} graphical results will be displayed in the active graphics window and also one pdf/png file will be saved within the \emph{plot_AeRobiology} directory automatically created in the working directory. By default, \code{export.plot = FALSE}.
#' @param export.format A \code{character} string specifying the format selected to save the plot. The implemented formats that may be used are: \code{"pdf"} or \code{"png"}. By default, \code{export.format = "pdf"}.
#' @param axisname A \code{character} string specifying the title of the y axis. By default, \code{axisname =  "Pollen grains / m3"}.
#' @param ... Other additional arguments may be used to customize the exportation of the plots using \code{"pdf"} or \code{"png"} files and therefore arguments from functions \code{\link[grDevices]{pdf}} and \code{\link[grDevices]{png}} may be implemented. For example, for pdf files the user may custom the arguments: width, height, family, title, fonts, paper, bg, fg, pointsize...; and for png files the user may custom the arguments: width, height, units, pointsize, bg, res...
#'
#' @details This function allows to summarize the pollen season by a simple plot. Even though the package was originally designed to treat aeropalynological data, it can be used to study many other atmospheric components (e.g., bacteria in the air, fungi, insects ...) \emph{(Buters et al., 2018; Oteros et al., 2019)}.
#' @return This function returns plot of class \pkg{ggplot2}. User are able to customize the output as a \pkg{ggplot2} object.
#' @references Buters, J. T. M., Antunes, C., Galveias, A., Bergmann, K. C., Thibaudon, M., Galan, C. & Oteros, J. (2018). Pollen and spore monitoring in the world. \emph{Clinical and translational allergy}, 8(1), 9.
#' @references Oteros, J., Bartusel, E., Alessandrini, F., Nunez, A., Moreno, D. A., Behrendt, H., ... & Buters, J. (2019). Artemisia pollen is the main vector for airborne endotoxin. \emph{Journal of Allergy and Clinical Immunology}.
#' @seealso \code{\link{calculate_ps}}; \code{\link{plot_normsummary}}
#' @examples  data("munich_pollen")
#' @examples  plot_summary(munich_pollen, pollen = "Betula", export.plot = FALSE, interpolation = FALSE)
#' @importFrom graphics plot
#' @importFrom utils data
#' @importFrom ggplot2 aes element_text geom_area geom_line ggplot ggtitle labs theme theme_bw theme_classic theme_set
#' @importFrom grDevices dev.off pdf png
#' @importFrom lubridate is.POSIXt yday year
#' @importFrom stats aggregate
#' @importFrom tidyr %>%
#' @export
plot_summary<-function (data,
                      pollen,
                      mave=1,
                      normalized=FALSE,
                      interpolation = TRUE,
                      int.method = "lineal",
                      export.plot = FALSE,
                      export.format = "pdf",
                      axisname="Pollen grains / m3",
                      ...){

  # Sys.setlocale(category = "LC_ALL", locale="english")

  #############################################    CHECK THE ARGUMENTS       #############################

  if(export.plot == TRUE){ifelse(!dir.exists(file.path("plot_AeRobiology")), dir.create(file.path("plot_AeRobiology")), FALSE)}
  data<-data.frame(data)
  if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

  if(class(pollen) != "character") stop ("Please include only charactr values for 'pollen'")

  if(class(export.plot) != "logical") stop ("Please include only logical values for export.plot argument")

  if(export.format != "pdf" & export.format != "png") stop ("Please export.format only accept values: 'pdf' or 'png'")

  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])


  if(class(axisname) != "character") stop ("Please include only character values for 'axisname'")

  # if(class(interpolation) != "logical") stop ("Please include only logical values for interpolation argument")
  #
  # if(int.method != "lineal" & int.method != "movingmean" & int.method != "tseries" &  int.method != "spline") stop ("Please int.method only accept values: 'lineal', 'movingmean', 'tseries' or 'spline'")
  #


  ncolpollen<-which( colnames(data)==pollen)
  data<-data[,c(1,ncolpollen)]


if(interpolation == TRUE){data <- interpollen(data, method = int.method)}

  data[, 1] <- as.Date(data[, 1])
  colnames(data)[1]<-"Date"
  data$DOY<-yday(data[, 1] )
  data$Year<-year(data[, 1] )
  data[,2]<-ma(data[,2],man=mave)

  if (normalized==T){
    sum<-aggregate(data[,pollen], by=list(data[,"Year"]), FUN=sum, na.rm=T)
    colnames(sum)[1]<-"Year"
    data<-merge(data,sum,by="Year")
    colnames(data)[length(data)]<-"sum"

    data$percent<-data[,pollen]*100/data[,"sum"]
    data<-data[which(data$sum!=0),]
    meant<-aggregate(data[,"percent"], by=list(data[,"DOY"]), FUN=mean, na.rm=T)
    colnames(meant)[1]<-"DOY"
    data_pollen<-merge(data, meant, by ="DOY")
    data_pollen$Year<-as.character(year(data_pollen[,"Date"]))
    # data_pollen<-data_pollen[,c(colnames(data_pollen)[1:2],pollen,"Year","x")]


    plot.summary<-ggplot(data_pollen, aes(x=DOY, y=data_pollen[,"percent"],fill=Year)) +
      theme_set(theme_bw()) +
      geom_area(position = "identity", alpha=0.3)+
      geom_line(data=data_pollen, aes(y=x), size=1)+
      labs(x="Day of the year", y="percentage (%)") +
      theme(axis.text=element_text(size=10), axis.title=element_text(size=10)) +
      ggtitle(pollen) +
      theme_classic()+
      theme(plot.title=element_text( face="bold.italic", size=13))


  }else{

    sum<-aggregate(data[,pollen], by=list(data[,"Year"]), FUN=sum, na.rm=T)
    colnames(sum)[1]<-"Year"
    data<-merge(data,sum,by="Year")
    colnames(data)[length(data)]<-"sum"

    data$percent<-data[,pollen]*100/data[,"sum"]
    data<-data[which(data$sum!=0),]
    meant<-aggregate(data[,pollen], by=list(data[,"DOY"]), FUN=mean, na.rm=T)
    colnames(meant)[1]<-"DOY"
    data_pollen<-merge(data, meant, by ="DOY")
    data_pollen$Year<-as.character(year(data_pollen[,"Date"]))
    data_pollen<-data_pollen[,c(colnames(data_pollen)[1:2],pollen,"Year","x")]


    plot.summary<-ggplot(data_pollen, aes(x=DOY, y=data_pollen[,pollen],fill=Year)) +
      theme_set(theme_bw()) +
      geom_area(position = "identity", alpha=0.3)+
      geom_line(data=data_pollen, aes(y=x), size=1)+
      labs(x="Day of the year", y=axisname) +
      theme(axis.text=element_text(size=10), axis.title=element_text(size=10)) +
      ggtitle(pollen) +
      theme_classic()+
      theme(plot.title=element_text(face="bold.italic", size=13))

  }




if(export.plot == TRUE  & export.format == "png") {
  png(paste0("plot_AeRobiology/plot_summary_", pollen,".png"), ...)
  plot(plot.summary)
  dev.off()
  png(paste0("plot_AeRobiology/credits.png"))

  dev.off()
}

if(export.plot == TRUE  & export.format == "pdf") {
  pdf(paste0("plot_AeRobiology/plot_summary_",  pollen, ".pdf"), ...)
  plot(plot.summary)

  dev.off()
}

return(plot.summary)

}
