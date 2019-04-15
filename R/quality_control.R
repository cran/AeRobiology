#' Quality Control of a Pollen Database
#'
#' Function to check the quality of an historical database of several pollen types.
#'@param data A \code{data.frame} object including the general database where quality must be checked. This data.frame must include a first column in \code{Date} format and the rest of columns in \code{numeric} format belonging to each pollen type by column. It is not necessary to insert the missing gaps; the function will automatically detect them.
#'@param int.window A  \code{numeric (interger)}value bigger or equal to \code{1}. The argument specifies the number of days of each side of the start, peak or end date of the main pollen season which will be checked during the quality control. If any of these days has been interpolated, the current season will not pass the quality control. The \code{int.window} argument will be \code{2} by default.
#'@param perc.miss A \code{numeric (interger)} value between \code{0} and \code{100}. The argument specifies the maximal percentage of interpolated days which is allowed inside the main pollen season to pass the quality control. The \code{perc.miss} argument will be \code{20} by default.
#'@param ps.method A \code{character} string specifying the method applied to calculate the pollen season and the main parameters. The implemented methods that can be used are: \code{"percentage"}, \code{"logistic"}, \code{"moving"}, \code{"clinical"} or \code{"grains"}. A more detailed information about the different methods for defining the pollen season may be consulted in \code{\link{calculate_ps}} function. The \code{ps.method} argument will be \code{"percentage"} by default.
#'@param result A \code{character} string specifying the format of the results. Only \code{"plot"} or \code{"table"} available. If \code{"plot"}, graphical resume of the quality control will be plotted. If \code{"table"}, a \code{data.frame} will be created indicating the filters passed by each pollen type and season. Consult 'Return' for more information .The \code{result} argument will be \code{"plot"} by default.
#'@param th.day See \code{\link{calculate_ps}} for more details.
#'@param perc See \code{\link{calculate_ps}} for more details.
#'@param def.season See \code{\link{calculate_ps}} for more details.
#'@param reduction See \code{\link{calculate_ps}} for more details.
#'@param red.level See \code{\link{calculate_ps}} for more details.
#'@param derivative See \code{\link{calculate_ps}} for more details.
#'@param man See \code{\link{calculate_ps}} for more details.
#'@param th.ma See \code{\link{calculate_ps}} for more details.
#'@param n.clinical See \code{\link{calculate_ps}} for more details.
#'@param window.clinical See \code{\link{calculate_ps}} for more details.
#'@param window.grains See \code{\link{calculate_ps}} for more details.
#'@param th.pollen See \code{\link{calculate_ps}} for more details.
#'@param th.sum See \code{\link{calculate_ps}} for more details.
#'@param type See \code{\link{calculate_ps}} for more details.
#'@param int.method See \code{\link{calculate_ps}} for more details.
#'@param ... Other arguments passed on to the pollen season calculation as specified in \code{\link{calculate_ps}} function.
#'@details Quality control is a relevant topic for aerobiology (Oteros et al., 2013). This function is another approach to improve the quality control management in the field. \cr\code{quality_control} function checks the quality of the pollen data of each pollen type and season. The filters applied by the function are: \cr
#'\itemize{
#'\item If the main pollen season (Galan et al., 2017) cannot be calculated according to \code{\link{calculate_ps}} function minimal requirements (lack of data for these pollen type and year). Filter named \code{"Complete"} in the \code{"quality_control"} \code{data.frame}.
#'\item If the start, end or peak date of the main pollen season has been interpolated or a day near to it (number of days specified by \code{int.window} argument). If a day near to these dates is missing, the selected date could not be the right one. Filters named \code{"Start"}, \code{"Peak"} and \code{"End"} in the \code{"quality_control"} \code{data.frame}.
#'\item The percentage of missing data inside the main pollen season. It calculates the number of days which have been interpolated by the algorithm and their percentage inside the main pollen season. If a high percentage of the main pollen season has been interpolated, the information of these season could not be reliable. Filter named \code{"Comp.MPS"} in the \code{"quality_control"} \code{data.frame}.
#'}
#'@return This function can return different results: \cr
#'\itemize{
#'\item If \code{result = "plot"}: Graphical resume of the Quality Control results showing the seasons of each pollen type and their quality (the risk assumed if they are included in further studies). The legend indicates the number of filter that have been unsuccessfully passed for each case. Object of class \code{\link[ggplot2]{ggplot}}. For graphical customization, see \code{\link[ggplot2]{ggplot}} function.
#'\item If \code{result = "table"}: \code{data.frame} with \code{logical} values for each pollen type and season. If \code{TRUE}, the filter has been successfully passed for this case. If FALSE, this case does not fit the minimal requirements of this filter.
#'}
#'@references Galan, C., Ariatti, A., Bonini, M., Clot, B., Crouzy, B., Dahl, A., Fernandez_Gonzalez, D., Frenguelli, G., Gehrig, R., Isard, S., Levetin, E., Li, D.W., Mandrioli, P., Rogers, C.A., Thibaudon, M., Sauliene, I., Skjoth, C., Smith, M., Sofiev, M., 2017. Recommended terminology for aerobiological studies. Aerobiologia (Bologna). 293_295.
#'@references Oteros, J., Galan, C., Alcazar, P., & Dominguez_Vilches, E. (2013). Quality control in bio_monitoring networks, Spanish Aerobiology Network. Science of the Total Environment, 443, 559_565.
#'@seealso \code{\link{calculate_ps}}, \code{\link{interpollen}}, \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{ggsave}}
#' @examples data("munich_pollen")
#' @examples quality_control(munich_pollen[,c(1:4)])
#' @importFrom utils data
#' @importFrom lubridate is.POSIXt
#' @importFrom ggplot2 aes element_blank element_text geom_tile ggplot ggsave labs scale_fill_gradient scale_x_discrete theme theme_classic
#' @importFrom graphics plot
#' @importFrom grDevices dev.off png
#' @importFrom tidyr %>%
#' @export
quality_control<-function(data,
                         int.window=2,
                         perc.miss=20,
                         ps.method="percentage",
                         result = "plot",
                         th.day=100,
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
                         int.method = "lineal",
                         ...){
  data<-data.frame(data)
  if (class(data) != "data.frame"& !is.null(data)){
    stop ("Please include a data.frame: first column with date, and the rest with pollen types")}
  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])
  if(class(int.window)!="numeric" | int.window %% 1 != 0){stop("int.window: Please, insert only an entire number bigger than 1")}
  if(int.window<1){stop("int.window: Please, insert only an entire number bigger or equal to 1")}
  if(class(perc.miss)!="numeric"){stop("perc.miss: Please, insert only a number between 0 and 100")}
  if(perc.miss<0 | perc.miss > 100){stop("perc.miss: Please, insert only a number between 0 and 100")}
  if(result!="plot" & result!="table"){stop("result: Please, insert only 'plot' or 'table'.")}

  Pollen<-calculate_ps(data=data, method=ps.method, th.day=th.day, perc=perc, def.season=def.season, reduction=reduction, red.level=red.level, derivative=derivative, man=man, th.ma=th.ma, n.clinical=n.clinical, window.clinical=window.clinical, window.grains=window.grains, th.pollen=th.pollen, th.sum=th.sum, type=type, interpolation=TRUE, int.method=int.method, plot=FALSE, maxdays = 300 )
  Interpolated<-interpollen(data=data, method=int.method, maxdays = 300, plot=FALSE, result="long")
  Qualitycontrol<-data.frame()
  Dataframe<-Pollen[,c(1,2)]
  Dataframe$Complete<-TRUE
  Dataframe$Start<-TRUE
  Dataframe$Peak<-TRUE
  Dataframe$End<-TRUE
  Dataframe$Comp.MPS<-TRUE
  nrow<-nrow(Dataframe)
  for(a in 1:nrow){
    if(any(is.na(Pollen[a,]))){
      Dataframe[a,-c(1,2)]<-FALSE
    }else{
      #Start
      Interpolated$Type<-as.character(Interpolated$Type)
      Interwindow<-Interpolated[which(Interpolated$Type==as.character(Pollen[a,1]) & Interpolated$Date>=(Pollen[a,3]-int.window) & Interpolated$Date<=(Pollen[a,3]+int.window) ),]
      if(sum(Interwindow[,4],na.rm = T)!=0){
        Dataframe[a,4]<-FALSE
      }
      #Peak
      Interwindow2<-Interpolated[which(Interpolated$Type==as.character(Pollen[a,1]) & Interpolated$Date>=(Pollen[a,11]-int.window) & Interpolated$Date<=(Pollen[a,11]+int.window) ),]
      if(sum(Interwindow2[,4],na.rm = T)!=0){
        Dataframe[a,5]<-FALSE
      }
      #End
      Interwindow3<-Interpolated[which(Interpolated$Type==as.character(Pollen[a,1]) & Interpolated$Date>=(Pollen[a,5]-int.window) & Interpolated$Date<=(Pollen[a,5]+int.window) ),]
      if(sum(Interwindow3[,4],na.rm = T)!=0){
        Dataframe[a,6]<-FALSE
      }
      MPS<-Interpolated[which(Interpolated$Type==as.character(Pollen[a,1]) & Interpolated$Date>=Pollen[a,3] & Interpolated$Date<=Pollen[a,5]),]
      Day.dif<-as.numeric(Pollen[a,5]-Pollen[a,3])
      if(sum(MPS[,4],na.rm=T)>round(Day.dif*(perc.miss/100))){
        Dataframe[a,7]<-FALSE
      }
    }
  }
  Dataframe$Risk<-apply(Dataframe[,-c(1,2)], 1, function(x) length(x[x==FALSE]) )
  graph<-ggplot(Dataframe, aes(seasons, type))+
    geom_tile(aes(fill=Risk), colour="grey")+
    scale_fill_gradient(low="white", high="#cb181d", limits=c(0,5))+
    scale_x_discrete(limits = unique(Dataframe$seasons))+
    theme_classic()+
    labs(title="Quality Control")+
    theme(plot.title = element_text(hjust=0.5), title = element_text(size = 16, face="bold"),axis.title=element_blank(), axis.text = element_text(size = 12, face="bold", colour= "black"), axis.text.y = element_text(face="bold.italic"), legend.title = element_text(size = 14), legend.text = element_text(size = 12))
  if(result=="plot"){
    return(graph)
  }
  if(result=="table"){
    return(Dataframe)
  }
}
