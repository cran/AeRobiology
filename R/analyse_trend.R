#' Calculating and Plotting Trends of Pollen Data (summary plot).
#'
#' Function to calculate the main seasonal indexes of the pollen season (\emph{Start Date}, \emph{Peak Date}, \emph{End Date} and \emph{Pollen Integral}). Trends analysis of the parameters over the seasons. Summary dot plot showing the distribution of the main seasonal indexes over the years.
#'
#' @param data A \code{data.frame} object. This \code{data.frame} should include a first column in format \code{Date} and the rest of columns in format \code{numeric} belonging to each pollen type by column.
#' @param interpolation A \code{logical} value specifying if the visualization shows the gaps in the inputs data (\code{interpolation = FALSE}) or if an interpolation method is used for filling the gaps (\code{interpolation = TRUE}). By default, \code{interpolation = TRUE}.
#' @param int.method A \code{character} string with the name of the interpolation method to be used. The implemented methods that may be used are: \code{"lineal"}, \code{"movingmean"}, \code{"tseries"} or \code{"spline"}. By default, \code{int.method = "lineal"}.
#' @param export.plot A \code{logical} value specifying if a plot will be exported or not. If \code{FALSE} graphical results will only be displayed in the active graphics window. If \code{TRUE} graphical results will be displayed in the active graphics window and also one pdf/png file will be saved within the \emph{plot_AeRobiology} directory automatically created in the working directory. By default, \code{export.plot = TRUE}.
#' @param export.format A \code{character} string specifying the format selected to save the plot. The implemented formats that may be used are: \code{"pdf"} or \code{"png"}. By default, \code{export.format = "pdf"}.
#' @param export.result A \code{logical} value. If \code{export.result = TRUE}, a table is exported with the extension \emph{.xlsx}, in the directory \emph{table_AeRobiology}. This table has the information about the \code{slope} \emph{"beta coefficient of a lineal model using as predictor the year and as dependent variable one of the main pollen season indexes"}. The information is referred to the main pollen season indexes: \emph{Start Date}, \emph{Peak Date}, \emph{End Date} and \emph{Pollen Integral}.
#' @param method A \code{character} string specifying the method applied to calculate the pollen season and the main seasonal parameters. The implemented methods that can be used are: \code{"percentage"}, \code{"logistic"}, \code{"moving"}, \code{"clinical"} or \code{"grains"}. By default, \code{method = "percentage"} (\code{perc = 95}\%). A more detailed information about the different methods for defining the pollen season may be consulted in the function \code{\link{calculate_ps}}.
#' @param quantil A \code{numeric} value (between 0 and 1) indicating the quantile of data to be displayed in the graphical output of the function. \code{quantil = 1} would show all the values, however a lower quantile will exclude the most extreme values of the sample. To split the parameters using a different sampling units (e.g. dates and pollen concentrations) can be used low vs high values of \code{quantil} argument (e.g. 0.5 vs 1). Also can be used an extra argument: \code{split}. By default, \code{quantil = 0.75}. \code{quantil} argument can only be applyed when \code{split = FALSE}.
#' @param significant A \code{numeric} value indicating the significant level to be considered in the linear trends analysis. This \emph{p} level is displayed in the graphical output of the function. By default, \code{significant = 0.05}.
#' @param split A \code{logical} argument. If \code{split = TRUE}, the plot is separated in two according to the nature of the variables (i.e. dates or pollen concentrations). By default, \code{split = TRUE}.
#'@param result A \code{character} object with the definition of the object to be produced by the function. If \code{result == "plot"}, the function returns a list of objects of class \pkg{ggplot2}; if \code{result == "table"}, the function returns a \pkg{data.frame}. By default, \code{result = "table"}.
#' @param ... Additional arguments for the function \code{\link{calculate_ps}} are also accepted.
#' @details This function allows to study time series trends of the pollen season. Even though the package was originally designed to treat aeropalynological data, it can be used to study many other atmospheric components (e.g., bacteria in the air, fungi, insects ...) \emph{(Buters et al., 2018; Oteros et al., 2019)}. The study of trends in pollen time series is a common approach to study the impact of climate change or other environmental factors on vegetation (Galan et al., 2016; Garcia_Mozo et al., 2016;  Recio et al., 2018). This tool can also be useful for studying trends in other fields (Oteros et al., 2015).
#' @return If \code{result == "plot"}, the function returns a list of objects of class \pkg{ggplot2}; if \code{result == "table"}, the function returns a \pkg{data.frame} with the hourly patterns.
#' The plot is of the class \pkg{ggplot2} or a list of plots of the class \pkg{ggplot2} (depending on the argument \code{split}). This is a combined dot plot showing the trends (\emph{slope} and \emph{p} value) of the main seasonal features.\cr
#' The object of the class \code{data.frame} has the information about the \code{slope} \emph{(beta coefficient of a lineal model using as predictor the year and as dependent variable one of the main pollen season indexes)}. The information is referred to the main pollen season indexes: \emph{Start Date}, \emph{Peak Date}, \emph{End Date} and \emph{Pollen Integral}.
#'
#' @references Buters, J. T. M., Antunes, C., Galveias, A., Bergmann, K. C., Thibaudon, M., Galan, C., ... & Oteros, J. (2018). Pollen and spore monitoring in the world. \emph{Clinical and translational allergy}, 8(1), 9.
#' @references Galan, C., Alcazar, P., Oteros, J., Garcia_Mozo, H., Aira, M. J., Belmonte, J., ... & Perez_Badia, R. (2016). Airborne pollen trends in the Iberian Peninsula. \emph{Science of the Total Environment}, 550, 53_59.
#' @references Garcia_Mozo, H., Oteros, J. A., & Galan, C. (2016). Impact of land cover changes and climate on the main airborne pollen types in Southern Spain. \emph{Science of the Total Environment}, 548, 221_228.
#' @references Oteros, J., Garcia_Mozo, H., Botey, R., Mestre, A., & Galan, C. (2015). Variations in cereal crop phenology in Spain over the last twenty_six years (1986_2012). \emph{Climatic Change}, 130(4), 545_558.
#' @references Oteros, J., Bartusel, E., Alessandrini, F., Nunez, A., Moreno, D. A., Behrendt, H., ... & Buters, J. (2019). Artemisia pollen is the main vector for airborne endotoxin. \emph{Journal of Allergy and Clinical Immunology}.
#' @references Recio, M., Picornell, A., Trigo, M. M., Gharbi, D., Garcia_Sanchez, J., & Cabezudo, B. (2018). Intensity and temporality of airborne Quercus pollen in the southwest Mediterranean area: Correlation with meteorological and phenoclimatic variables, trends and possible adaptation to climate change. \emph{Agricultural and Forest Meteorology}, 250, 308_318.
#' @seealso \code{\link{calculate_ps}}; \code{\link{plot_trend}}
#' @examples data("munich_pollen")
#' @examples analyse_trend(munich_pollen, interpolation = FALSE, export.result = FALSE, export.plot = FALSE)
#' @importFrom graphics plot
#' @importFrom utils data
#' @importFrom ggplot2 aes element_blank element_text geom_point geom_vline ggplot guide_legend guides labs scale_colour_manual scale_fill_discrete scale_fill_gradientn scale_x_continuous theme theme_bw
#' @importFrom grDevices dev.off pdf png rainbow
#' @importFrom grid grid.layout pushViewport viewport
#' @importFrom lubridate is.POSIXt
#' @importFrom stats as.formula complete.cases lm pf quantile
#' @importFrom writexl write_xlsx
#' @importFrom tidyr %>%
#' @export
analyse_trend  <-   function   (data,
                            interpolation = TRUE,
                            int.method = "lineal",
                            export.plot = TRUE,
                            export.format = "pdf",# pdf, png
                            export.result=TRUE,
                            method="percentage",
                            quantil=0.75,
                            significant=0.05,
                            split=TRUE, result="table",
                            ...){

  # Sys.setlocale(category = "LC_ALL", locale="english")

  #############################################    CHECK THE ARGUMENTS       #############################

  if(export.plot == TRUE){ifelse(!dir.exists(file.path("plot_AeRobiology")), dir.create(file.path("plot_AeRobiology")), FALSE)}

  if(export.result == TRUE){ifelse(!dir.exists(file.path("table_AeRobiology")), dir.create(file.path("table_AeRobiology")), FALSE)}

   data<-data.frame(data)
  if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

  if(class(export.plot) != "logical") stop ("Please include only logical values for export.plot argument")

  if(export.format != "pdf" & export.format != "png") stop ("Please export.format only accept values: 'pdf' or 'png'")

  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])

  if(class(interpolation) != "logical") stop ("Please include only logical values for interpolation argument")

  if(class(split) != "logical") stop ("Please include only logical values for split argument")

  if(class(quantil) != "numeric") stop ("Please include only logical values for quantil argument")

  if(class(significant) != "numeric") stop ("Please include only logical values for significant argument")

  if(method != "percentage" & method != "logistic" & method != "moving" & method != "clinical" & method != "grains") stop ("Please method only accept values: 'percentage', 'logistic', 'moving', 'clinical' or 'grains'")

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

  datafram<-calculate_ps(data,method=method, interpolation =  interpolation, int.method=int.method,plot=FALSE,...)

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
      print(type)
      }



trendtimeresult<-trendtime


### Now start the rock and roll

trendtime$variable<-as.factor(trendtime$variable)

  trendtime2<-trendtime
  trendtime2<-trendtime2[complete.cases(trendtime2),]
  sig<-paste("Significant (p<",significant,")", sep="")
  trendtime2$significant<-sig
  #trendtime2[which(trendtime2$p<=significant),"significant"]<-sig
  trendtime2[which(trendtime2$p>significant),"significant"]<-"Non Significant"
  trendtime2[which(trendtime2$p>significant),"p"]<-0.99
  trendtime2[which(trendtime2$p!=0.99),"p"]<-0.01

trendtime2$variable<-as.character(trendtime2$variable)

trendtime2[which(trendtime2$variable=="st.jd"),"variable"]  <- "Start Date"
trendtime2[which(trendtime2$variable=="pk.jd"),"variable"]  <- "Peak Date"
trendtime2[which(trendtime2$variable=="en.jd"),"variable"]  <- "End Date"
trendtime2[which(trendtime2$variable=="sm.ps"),"variable"]  <- "Total Pollen"

trendtime2$variable<-as.factor(as.character(trendtime2$variable))

trend_p1 <- ggplot(trendtime2, aes(x=coef, y=type))+
   theme_bw()+
   geom_vline(aes(xintercept=0), colour="red", linetype = "dashed", size=1)+
   geom_point(aes(fill=as.numeric(variable), colour=trendtime2$significant, stroke=2-trendtime2$p*2),size=6, alpha=0.9,shape=21)+
   scale_fill_gradientn(colours = rainbow(6), labels = levels(trendtime2$variable))+
   scale_colour_manual(values = c("gray80","black"))+
   scale_x_continuous(limits=c(-quantile(abs(trendtime2$coef),quantil),quantile(abs(trendtime2$coef),quantil)))+
   labs(x= "slope", y="")+
   guides(fill=guide_legend(title=NULL))+
   theme(legend.position="top", legend.title = element_blank(), axis.text.y = element_text(face="italic"))


 filtertrend<-trendtime2[which(trendtime2$variable!="Total Pollen"),]
 filtertrend$variable<-as.factor(as.character(filtertrend$variable))
 filtertrend$type<-as.factor(as.character(filtertrend$type))

 trend_p2 <- ggplot(filtertrend, aes(x=coef, y=type))+
   theme_bw()+
   geom_vline(aes(xintercept=0), colour="red", linetype = "dashed", size=1)+
   geom_point(aes(fill=filtertrend$variable, colour=filtertrend$significant, stroke=2-filtertrend$p*2),size=6, alpha=0.9,shape=21)+
   scale_colour_manual(values = c("gray80","black"))+
   scale_x_continuous(limits=c(-max(abs(filtertrend$coef), na.rm=T),max(abs(filtertrend$coef), na.rm=T)))+
   labs(x= "slope", y="")+
   guides(fill=guide_legend(title=NULL))+
   theme(legend.position="top", legend.title = element_blank(), axis.text.y = element_text(face="italic"))


 filtertrend2<-trendtime2[which(trendtime2$variable=="Total Pollen"),]
 filtertrend2$variable<-as.factor(as.character(filtertrend2$variable))
 filtertrend2$type<-as.factor(as.character(filtertrend2$type))

 trend_p3 <- ggplot(filtertrend2, aes(x=coef, y=type))+
   theme_bw()+
   geom_vline(aes(xintercept=0), colour="red", linetype = "dashed", size=1)+
   geom_point(aes( colour=filtertrend2$significant, stroke=2-filtertrend2$p*2),size=6, alpha=0.9,shape=21,fill="lightpink")+
   scale_colour_manual(values = c("gray80","black"), name = "Total Pollen: ")+
   scale_x_continuous(limits=c(-max(abs(filtertrend2$coef), na.rm=T),max(abs(filtertrend2$coef), na.rm=T)))+
   labs(x= "slope", y="")+
   theme(legend.position="top", axis.text.y = element_text(face="italic"))+
   guides(fill=guide_legend(title=""))


 if (split == F){
   analyse_trend<-trend_p1
   if(export.plot == TRUE  & export.format == "png") {
     png(paste0("plot_AeRobiology/analyse_trend",".png"), ...)
     plot(analyse_trend)
     dev.off()

   }

   if(export.plot == TRUE  & export.format == "pdf") {
     pdf(paste0("plot_AeRobiology/analyse_trend", ".pdf"))
     plot(analyse_trend)

     dev.off()
   }

 }else{
   analyse_trend<-list()
   analyse_trend[[1]]<-trend_p2
   analyse_trend[[2]]<-trend_p3


   if(export.plot == TRUE  & export.format == "png") {
     png("plot_AeRobiology/analyse_trend_split.png")
     pushViewport(viewport(layout=grid.layout(2,1)))
     vplayout<-function(x,y) viewport(layout.pos.row = x, layout.pos.col=y)
     print(trend_p2, vp = vplayout(1,1))
     print(trend_p3, vp = vplayout(2,1))
     dev.off()

   }

   if(export.plot == TRUE  & export.format == "pdf") {
     pdf("plot_AeRobiology/analyse_trend_split.pdf", ...)
     pushViewport(viewport(layout=grid.layout(2,1)))
     vplayout<-function(x,y) viewport(layout.pos.row = x, layout.pos.col=y)
     print(trend_p2, vp = vplayout(1,1))
     print(trend_p3, vp = vplayout(2,1))

     dev.off()
   }
   }

 lista<-list()
 lista[["trends"]]<-trendtime


 lista [["Information"]] <- data.frame(
   Attributes = c("st.jd", "pk.jd", "en.jd", "sm.ps", "coef", "p", "", "", "Package", "Authors"),
   Description = c("Start-date (day of the year)","Peak-date (day of year)",  "End-date (day of the year)", "Pollen integral", "Slope of the linear trend", "Significance level of the linear trend",  "", "", "AeRobiology", "Jesus Rojo, Antonio Picornell & Jose Oteros"))

 if (export.result == TRUE) {
   write_xlsx(lista, "table_AeRobiology/summary_of_trends.xlsx")
 }

 if (result=="table"){
   return(trendtimeresult)
 }else if (result=="plot"){
   return(analyse_trend)
 }
}
