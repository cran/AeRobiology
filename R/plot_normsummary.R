#' Plotting the Amplitude of Several Pollen Seasons.
#'
#' Function to plot the pollen data amplitude during several seasons: daily average pollen concentration over the study period, maximum pollen concentration of each day over the study period and minimum pollen concentration of each day value over the study period. It is possible to plot the relative abundance per day and smoothing the pollen season by calculating a moving average.
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
#' @param color.plot A \code{character} string. The argument defines the color to fill the plot. Will be \code{"orange2"} by default.
#' @param ... Other additional arguments may be used to customize the exportation of the plots using \code{"pdf"} or \code{"png"} files and therefore arguments from functions \code{\link[grDevices]{pdf}} and  \code{\link[grDevices]{png}} may be implemented. For example, for pdf files the user may custom the arguments: width, height, family, title, fonts, paper, bg, fg, pointsize...; and for png files the user may custom the arguments: width, height, units, pointsize, bg, res...
#' @return This function returns plot of class \pkg{ggplot2}. User are able to customize the output as a \pkg{ggplot2} object.
#' @seealso \code{\link{calculate_ps}}; \code{\link{plot_summary}}
#' @examples data("munich_pollen")
#' @examples plot_normsummary(munich_pollen, pollen = "Betula", interpolation = FALSE, export.plot = FALSE)
#' @importFrom graphics plot
#' @importFrom utils data
#' @importFrom ggplot2 aes element_text geom_line geom_ribbon ggplot labs scale_colour_manual scale_x_date theme theme_classic
#' @importFrom grDevices dev.off pdf png
#' @importFrom lubridate is.POSIXt yday year
#' @importFrom scales date_format
#' @importFrom stats aggregate
#' @importFrom tidyr %>%
#' @export
plot_normsummary<-function (data,
                        pollen,
                        mave=1,
                        normalized=FALSE,
                        interpolation = TRUE,
                        int.method = "lineal",
                        export.plot = FALSE,
                        export.format = "pdf",
                        color.plot="orange2",
                        axisname="Pollen grains / m3",
                        ...){

  # Sys.setlocale(category = "LC_ALL", locale="english")

  #############################################    CHECK THE ARGUMENTS       #############################

  if(export.plot == TRUE){ifelse(!dir.exists(file.path("plot_AeRobiology")), dir.create(file.path("plot_AeRobiology")), FALSE)}
  data<-data.frame(data)
  if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

  if(class(pollen) != "character") stop ("Please include only character values for 'pollen'")

  if(class(export.plot) != "logical") stop ("Please include only logical values for export.plot argument")

  if(export.format != "pdf" & export.format != "png") stop ("Please export.format only accept values: 'pdf' or 'png'")

  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])

  if(class(interpolation) != "logical") stop ("Please include only logical values for interpolation argument")

  if(class(axisname) != "character") stop ("Please include only character values for 'axisname'")

  if(class(color.plot) != "character") stop ("Please include only character values indicating the name of a colour for 'color.plot'")

  # if(int.method != "lineal" & int.method != "movingmean" & int.method != "spline") stop ("Please int.method only accept values: 'lineal', 'movingmean' or 'spline'")

  ncolpollen<-which( colnames(data)==pollen)
  data<-data[,c(1,ncolpollen)]


  if(interpolation == TRUE){data <- interpollen(data, method = int.method)}
  data[, 1] <- as.Date(data[, 1])
  colnames(data)[1]<-"Date"
  data$DOY<-yday(data[, "Date"])
  data$Year<-year(data[, "Date"])

  sum<-aggregate(data[,pollen], by=list(data[,"Year"]), FUN=sum, na.rm=T)
  colnames(sum)[1]<-"Year"
  data<-merge(data,sum,by="Year")
  colnames(data)[length(data)]<-"sum"

  max<-aggregate(data[,pollen], by=list(data[,"Year"]), FUN=max, na.rm=T)
  colnames(max)[1]<-"Year"
  data<-merge(data,max,by="Year")
  colnames(data)[length(data)]<-"max"

  min<-aggregate(data[,pollen], by=list(data[,"Year"]), FUN=min, na.rm=T)
  colnames(min)[1]<-"Year"
  data<-merge(data,min,by="Year")
  colnames(data)[length(data)]<-"min"

  data$Norm<-data[,pollen]*100/data$sum
  data$noNorm<-data[,pollen]


  if (normalized == TRUE){
    meant<-aggregate(data[,"Norm"], by=list(data[,"DOY"]), FUN=mean, na.rm=T)

    mint<-aggregate(data[,"Norm"], by=list(data[,"DOY"]), FUN=min, na.rm=T)
    maxt<-aggregate(data[,"Norm"], by=list(data[,"DOY"]), FUN=max, na.rm=T)

    frame<-data.frame(DOY=meant$Group.1,Mean=meant$x,Min=mint$x,Max=maxt$x, date=as.Date(meant$Group.1, origin = "2000-01-01"))

    frame<-frame[!is.infinite(frame[,3]),]
    frame$Mean<-ma(frame$Mean,man=mave)
    frame$Min<-ma(frame$Min,man=mave)
    frame$Max<-ma(frame$Max,man=mave)

    plot.summary.norm<-ggplot(data = frame, aes(x = date, y = Mean)) +
      labs(x="", y="percentage (%)", title=pollen) +
      geom_ribbon(aes(ymin=Min,ymax=Max),alpha=0.3,fill=color.plot)+
      geom_line(size = 0.8)+
      scale_colour_manual("Parameters",values=c("black"))+
      theme_classic()+
      theme(text = element_text(size = 13), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))+
      scale_x_date(labels=date_format("%b"), date_breaks = "1 month")+
      theme(plot.title=element_text(face="bold.italic", size=13))


  } else {
    meant<-aggregate(data[,"noNorm"], by=list(data[,"DOY"]), FUN=mean, na.rm=T)
    mint<-aggregate(data[,"noNorm"], by=list(data[,"DOY"]), FUN=min, na.rm=T)
    maxt<-aggregate(data[,"noNorm"], by=list(data[,"DOY"]), FUN=max, na.rm=T)

    frame<-data.frame(DOY=meant$Group.1,Mean=meant$x,Min=mint$x,Max=maxt$x, date=as.Date(meant$Group.1, origin = "2000-01-01"))

    frame<-frame[!is.infinite(frame[,3]),]
    frame$Mean<-ma(frame$Mean,man=mave)
    frame$Min<-ma(frame$Min,man=mave)
    frame$Max<-ma(frame$Max,man=mave)

    plot.summary.norm<-ggplot(data = frame, aes(x = date, y = Mean)) +
      labs(x="", y=axisname, title=pollen) +
      geom_ribbon(aes(ymin=Min,ymax=Max),alpha=0.3,fill=color.plot)+
      geom_line(size = 0.8)+
      scale_colour_manual("Parameters",values=c("black"))+
      theme_classic()+
      theme(text = element_text(size = 13), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))+
      scale_x_date(labels=date_format("%b"), date_breaks = "1 month")+
      theme(plot.title=element_text(face="bold.italic", size=13))


  }




  if(export.plot == TRUE  & export.format == "png") {
    png(paste0("plot_AeRobiology/plot_normsummary_", pollen,".png"), ...)
    plot(plot.summary.norm)
    dev.off()
    png(paste0("plot_AeRobiology/credits.png"))

    dev.off()
  }

  if(export.plot == TRUE  & export.format == "pdf") {
    pdf(paste0("plot_AeRobiology/plot_normsummary_", pollen,".pdf"), ...)
    plot(plot.summary.norm)

    dev.off()
  }

  return(plot.summary.norm)

}
