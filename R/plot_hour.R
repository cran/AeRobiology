#'Plotting hourly patterns
#'
#'Function to plot pollen data expressed in concentrations with time resolution higher than 1 day (e.g. hourly, bi-hourly concentrations).
#'
#'@param data A \code{data.frame} object with the structure \code{long}. Where the first two columns are \code{factors} indicating the \code{pollen} and the \code{location}. The 3 and 4 columns are \code{POSIXct}, showing the hour. Where the third column is the beguinning of the concentration \code{from} and the fourth column is the end time of the concentrtion data \code{to}. The fifth column shows the concentrations of the different pollen types as \code{numeric}. Please see the example 3-hourly data from the automatic pollen monitor BAA500 from Munich and Viechtach in Bavaria (Germany) \code{data("POMO_pollen")}, supplied by ePIN Network supported by the Bavarian Government.
#'@param result A \code{character} object with the definition of the object to be produced by the function. If \code{result == "plot"}, the function returns a list of objects of class \pkg{ggplot2}; if \code{result == "table"}, the function returns a \pkg{data.frame} with the hourly patterns. By default, \code{result = "plot"}.
#'@param locations A \code{logical} object with the specification if the different locations will be displayed in the plot. Argument only used when \code{result == "plot"}. By default, \code{locations = FALSE}.
#'@return If \code{result == "plot"}, the function returns a list of objects of class \pkg{ggplot2}; if \code{result == "table"}, the function returns a \pkg{data.frame} with the hourly patterns.
#'@references Oteros, J., Pusch, G., Weichenmeier, I., Heimann, U., Mueller, R., Roeseler, S., ... & Buters, J. T. (2015). Automatic and online pollen monitoring. \emph{International archives of allergy and immunology}, 167(3), 158-166.
#'@examples data("POMO_pollen")
#'@examples plot_hour(POMO_pollen, result="plot", locations = FALSE)
#'@examples plot_hour(POMO_pollen, result="plot", locations = TRUE)
#'@importFrom stats aggregate
#'@importFrom data.table data.table
#'@importFrom lubridate hour yday
#'@importFrom dplyr group_by summarise
#'@importFrom tidyr spread
#'@importFrom ggplot2 ggplot scale_fill_brewer
#'@export
plot_hour <-
  function (data, result="plot", locations = FALSE) {
    data<-data.frame(data)
    if(class(data) != "data.frame") stop ("Please include a data.frame: first column with factor indicating the pollen, second column with factor indicating the locaiton, third column with POSIXct indicating the (from), fourth column with POSIXct indicating the (to) and fifth column with numbers indicating the concentration")
    if(result != "plot" & result != "table") stop ("Please result only accept values: 'table' or 'plot'")
    if(class(locations) != "logical") stop ("Please include only logical values for locations argument")

    colnames(data)<-c("pollen","location","from","to","value")
# HOUR = percent = location = NULL
    data$date<-as.Date(as.character(data$from),"%Y-%m-%d")
    data$hour<-hour(strptime(as.character(data$from), format = '%Y-%m-%d %H:%M', 'GMT'))
    data$hour2<-hour(strptime(as.character(data$to), format = '%Y-%m-%d %H:%M', 'GMT'))
    data$Hour<-paste0(data$hour,"-",data$hour2)
    # data$Hour<-factor(data$Hour,levels = c("0-3","3-6","6-9","9-12","12-15","15-18","18-21","21-0"))

    ### Data to 24h concentrations AND POLLENS selection ###

    dayconc<-data.table(data)
    colnames(dayconc)[1]<-"pollen"

    table<-aggregate( value ~ date+pollen+location, data=dayconc, FUN=mean)
    table2<-aggregate( value ~ date+pollen+location, data=dayconc, FUN=max)
    table$max<-table2$value

    mainpollen<-data%>%group_by(pollen)%>%summarise(total=sum(value, na.rm=TRUE))
    mainpollen<-as.data.frame(mainpollen)

    pollens <-mainpollen[,"pollen"]


    pollens<-as.character(pollens)

    table<-table[which(as.character(table$pollen) %in% pollens),]
    table$location<-as.character(table$location)
    table$location<-as.factor(table$location)
    ##### Calculate hourly frames
    frame3<-data.frame()
    data_end<-data.frame()
    if (locations == TRUE){
    for (loc in levels(table$location)){

      tabletemp<-dayconc[which(dayconc$location==loc),]
      tabletemp$location<-as.character(tabletemp$location)
      data_wide <- spread(tabletemp, pollen, value)
      data_wide$location<-loc
      data_end<-rbind(data_end,data_wide)
      tabletemp<-as.data.frame(tabletemp)
      for ( a in 1:length(pollens)){
        tempo<-tabletemp[which(tabletemp$pollen==pollens[a]),]

        if(nrow(tempo)==0){next}
        for ( b in 1:nrow(tempo)){

          tempo[b,"percent"]<-tempo[b,"value"]*100/table[which(table$pollen==pollens[a] & table$date==tempo[b,"date"] & table$location==tempo[b,"location"]),"value"]
        }
        tempo$percent<-tempo$percent/length(unique(data$Hour))

        frame3<-rbind(frame3,tempo)
      }

      frame3<-frame3[complete.cases(frame3),]
    }
   }else {

      tabletemp<-dayconc

      data_wide <- spread(tabletemp, pollen, value)
      data_wide$location<-"All"
      data_end<-data_wide
      tabletemp<-as.data.frame(tabletemp)
      for ( a in 1:length(pollens)){
        tempo<-tabletemp[which(tabletemp$pollen==pollens[a]),]

        if(nrow(tempo)==0){next}
        for ( b in 1:nrow(tempo)){

          tempo[b,"percent"]<-tempo[b,"value"]*100/table[which(table$pollen==pollens[a] & table$date==tempo[b,"date"] & table$location==tempo[b,"location"]),"value"]
        }
        tempo$percent<-tempo$percent/length(unique(data$Hour))

        frame3<-rbind(frame3,tempo)
      }

      frame3<-frame3[complete.cases(frame3),]

  }

    # frame3<-frame3[which(yday(frame3$date)%in%summerseries),]

    frame3 = frame3[order(frame3$hour), ]

    frame3$Hour = factor(frame3$Hour, levels = unique(frame3$Hour))

if (result == "plot" & locations==FALSE) {
  plottotal<-ggplot(frame3,aes(x=as.factor(Hour),y=percent))+
    geom_bar(stat = "summary", fun.y = "mean", color="black", fill="gray")+
    theme_bw()+
    labs(x="hour", y="daily pollen (%)", title="Total pollen")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pollenlist<-list()
  pollenlist[[1]]<-plottotal
  for ( a in 1:length(pollens)){
    tempo<-frame3[which(frame3$pollen==pollens[a]),]
    pollenlist[[a+1]]<- ggplot(tempo,aes(x=as.factor(Hour),y=percent))+
      geom_bar(stat = "summary", fun.y = "mean", color="black", fill="gray")+
      theme_bw()+
      labs(x="hour", y="daily pollen (%)", title=pollens[a])+
      theme(axis.text.x = element_text(angle = 45, hjust  = 1))


  }
  pollenlist
}else if (result == "plot" & locations==TRUE){
  plottotal<-ggplot(frame3,aes(x=as.factor(Hour),y=percent,fill=location))+
    geom_bar(stat = "summary", fun.y = "mean", color="black",position='dodge')+
    theme_bw()+
    labs(x="hour", y="daily pollen (%)", title="Total pollen")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_brewer(palette="Set1")

  pollenlist<-list()
  pollenlist[[1]]<-plottotal
  for ( a in 1:length(pollens)){
    tempo<-frame3[which(frame3$pollen==pollens[a]),]
    pollenlist[[a+1]]<- ggplot(tempo,aes(x=as.factor(Hour),y=percent, fill=location))+
      geom_bar(stat = "summary", fun.y = "mean", color="black",position='dodge')+
      theme_bw()+
      labs(x="hour", y="daily pollen (%)", title=pollens[a])+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_fill_brewer(palette="Set1")


  }
  pollenlist
}

    if (result=="plot"){
      return(pollenlist)
    } else if (result=="table"){
      return(frame3)
    }


  }
