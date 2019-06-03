#' Interpolation of Missing Data in a Pollen Database by Different Methods
#'
#' Function to simultaneously replace all missing data of an historical database of several pollen types by using different methods of interpolation.
#'
#' @param data A \code{data.frame} object including the general database where interpollation must be performed. This \code{data.frame} must include a first column in \code{Date} format and the rest of columns in \code{numeric} format. Each column must contain information of one pollen type. It is not necessary to insert missing gaps; the function will automatically detect them.
#' @param method A \code{character} string specifying the method applied to calculate and generate the pollen missing data. The implemented methods that can be used are: \code{"lineal"}, \code{"movingmean"}, \code{"spline"}, \code{"tseries"} or \code{"neighbour"}. A more detailed information about the different methods may be consulted in Details. The \code{method} argument will be \code{"lineal"} by default.
#' @param maxdays A \code{numeric (interger)} value specifying the maximum number of consecutive days with missing data that the algorithm is going to interpolate. If the gap is bigger than the argument value, the gap will not be interpolated. Not valid with \code{"tseries"} method. The \code{maxdays} argument will be \code{30} by default.
#' @param plot A \code{logical} argument. If \code{TRUE}, graphical previews of the input database will be plot at the end of the interpolation process. All the interpolated gaps will be marked in red. The \code{plot} argument will be \code{TRUE} by default.
#' @param factor A \code{numeric (interger)} value bigger than \code{1}. Only valid if the \code{"movingmean"} method is chosen. The argument specifies the factor which will multiply the gap size to stablish the range of the moving mean that will fulfill the gap. A more detailed information about the selection of the factor may be consulted in Details. The argument \code{factor} will be \code{1} by default.
#' @param ndays A \code{numeric (interger)} value bigger than \code{1}. Only valid if the \code{"spline"} method is chosen. Specifies the number of days beyond each side of the gap which are used to perform the spline regression. The argument \code{ndays} will be \code{3} by default.
#' @param spar A \code{numeric (double)} value ranging \code{0_1} specifying the degree of smoothness of the spline regression adjustment. As smooth as the adjustment is, more data are considered as outliers for the spline regression. Only valid if the \code{"spline"} method is chosen. The argument \code{"spar"} will be \code{0.5} by default.
#' @param data2,data3,data4,data5 A \code{data.frame} object (each one) including database of a neighbour pollen station which will be used to interpolate missing data in the target station. Only valid if the "neighbour" method is chosen. This \code{data.frame} must include a first column in \code{Date} format and the rest of columns in \code{numeric} format belonging to each pollen type by column. It is not necessary to insert the missing gaps; the function will automatically detect them. The arguments will be \code{NULL} by default.
#' @param mincorr A \code{numeric (double)} value ranging \code{0_1}. It specifies the minimal correlation coefficient (Spearman correlations) that neighbour stations must have with the target station to be taken into account for the interpolation. Only valid if the \code{"neighbour"} method is chosen. The argument \code{"mincorr"} will be \code{0.6} by default.
#' @param result A \code{character} string specifying the format of the resulting \code{data.frame}. Only \code{"wide"} or \code{"long"}. The \code{result} argument will be \code{"wide"} by default.
#' @details This function allows to interpolate missing data in a pollen database using 4 different methods which are described below. Interpolation for each pollen type will be automatically done for gaps smaller than the \code{"maxdays"} argument. \cr
#' \itemize{
#' \item \code{"lineal"} method. The interpolation will be carried out by tracing a straight line between the gap extremes.
#' \item \code{"movingmean"} method. It calculates the moving mean of the pollen daily concentrations with a window size of the gap size multiplicated by the \code{factor} argument and replace the missing data with the moving mean for these days. It is a dynamic function and for each gap of the database, the window size of the moving mean changes depending of each gap size.
#' \item \code{"spline"} method. The interpolation will be carried out by performing a spline regression with the previous and following days to the gap. The number of days of each side of the gap that will be taken into account for calculating the spline regression are specified by \code{ndays} argument. The smoothness of the adjustment of the spline regression can be specified by the \code{spar} argument.
#' \item \code{"tseries"} method. The interpolation will be carried out by analysing the time series of pollen database. It performs a seasonal_trend decomposition based on LOESS (Cleveland et al., 1990). The seasonality of the historical database is extracted and used to predict the missing data by performing a linear regression with the target year.
#' \item \code{"neighbour"} method. Other near stations provided by the user are used to interpolate the missing data of the target station. First of all, a Spearman correlation is performed between the target station and the neighbour stations to discard the neighbour stations with a correlation coefficient smaller than \code{mincorr} value. For each gap, a linear regression is performed between the neighbour stations and the target stations to determine the equation which converts the pollen concentrations of the neighbour stations into the pollen concentration of the target station. Only neighbour stations without any missing data during the gap period are taken into account for each gap.
#' }
#' @return This function returns different results:
#' \itemize{
#' \item If \code{result = "wide"}, returns a \code{data.frame} including the original data and completed with the interpolated data.
#' \item If \code{result = "long"}, returns a \code{data.frame} containing your data in long format (the first column for date, the second for pollen type, the third for concentration and an additional fourth column with \code{1} if this data has been interpolated or \code{0} if not).
#' \item If \code{plot = TRUE}, plots for each year and pollen type with daily values are represented in the active graphic window. Interpolated values are marked in red. If \code{method} argument is \code{"tseries"}, the seasonality is also represented in grey.
#' }
#' @references Cleveland RB, Cleveland WS, McRae JE, Terpenning I (1990) STL: a seasonal_trend decomposition procedure based on loess. J Off Stat 6(1):3_33.
#' @seealso \code{\link{ma}}
#' @examples data("munich_pollen")
#' @examples interpollen(munich_pollen, method = "lineal", plot = FALSE)
#' @importFrom lubridate is.POSIXt
#' @importFrom graphics legend lines par plot points
#' @importFrom stats aggregate cor.test lm na.omit predict predict.lm smooth.spline stl ts window
#' @importFrom zoo na.approx
#' @importFrom tidyr %>% gather spread
#' @export
interpollen <-
  function(data,
           method = "lineal" ,
           maxdays = 30,
           plot = TRUE,
           factor = 2,
           ndays = 3,
           spar = 0.5,
           data2=NULL,
           data3=NULL,
           data4=NULL,
           data5=NULL,
           mincorr=0.6,
           result="wide") {
    if(class(maxdays)!="numeric"){stop("maxday: Please, insert an entire number bigger than 1")}
    if(class(maxdays)!="numeric"| maxdays%%1!=0| maxdays<=0){stop("maxday: Please insert an entire number bigger than 1")}
    if(class(method)!="character"){stop("method: Please, only type 'lineal' , 'movingmean' or 'spline'")}
    if (method != "lineal" &
        method != "movingmean" &
        method != "spline" & method != "neighbour" & method != "tseries"){
      stop("method: Please, only type 'lineal' , 'movingmean', 'spline', 'tseries' or 'neighbour'")}
    if (result != "wide" &
        result != "long"){
      stop("result: Please, only type 'wide' or 'long'")}
    if (plot != TRUE &
        plot != FALSE){
      stop("plot: Please, only type 'TRUE' or 'FALSE'")}
    if (class(factor) != "numeric"){
      stop("factor: Please, insert only a number bigger than 1")}
    if (factor <= 1){
      stop("factor: Please, insert only a number bigger than 1")}
    if (class(ndays)!="numeric" | ndays %% 1 != 0){
      stop("ndays:Please, insert only an entire number bigger than 1")}
    if (class(spar)!="numeric" | spar <= 0 |
        spar > 1){
      stop("spar: Please, insert only a number between 0 (not included) and 1")}
    if (ndays <= 0){
      stop("ndays: Please, insert only an entire number bigger than 1")}
    if (ndays >= 15){
      warning("WARNING: Results of spline may not be reliable with more than 15 days")}
    if(mincorr<0 | mincorr>1){stop("mincorr: Please insert only a number between 0 and 1 ")}

    data<-data.frame(data)
    colnames(data)[1]<-"date"
    if (class(data) != "data.frame"& !is.null(data)){
      stop ("Please, include a data.frame: first column with date, and the rest with pollen types")}
    if (class(data2) != "data.frame" & !is.null(data2)){
      stop ("data2: Please, include a data.frame: first column with date, and the rest with pollen types")}
    if (class(data3) != "data.frame"& !is.null(data3)){
      stop ("data3: Please, include a data.frame: first column with date, and the rest with pollen types")}
    if (class(data4) != "data.frame"& !is.null(data4)){
      stop ("data4: Please, include a data.frame: first column with date, and the rest with pollen types")}
    if (class(data5) != "data.frame" & !is.null(data5)){
      stop ("data5: Please, include a data.frame: first column with date, and the rest with pollen types")}
    if (method=="neighbour" & is.null(data2)){stop("You need at least another station loaded in 'data2' to use this method")}
    if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
    data[,1]<-as.Date(data[,1])

    if (method == "tseries"){
      warning("WARNING: 'maxdays' argument doesn't work with this method. All gaps will be interpolated.")

      date.st <- min(data[ ,1], na.rm = TRUE)
      date.en <- max(data[ ,1], na.rm = TRUE)
      date.seq <- seq(from = date.st, to = date.en, "days")

      data1 <- data.frame(matrix(data = NA, nrow = length(date.seq), ncol = length(colnames(data))))
      colnames(data1) <- colnames(data)
      data1[ ,1] <- as.Date(date.seq)
      data1[data1[ ,1] %in% data[ ,1], -1] <- data[data[ ,1] %in% data1[ ,1], -1]
      data <- data1

      data.full <- data

      type.name <- colnames(data)


      for (cols in 2: length(type.name)){

        ma.data <- data.frame(date = data[ ,1],
                              year = strftime(data[ ,1], "%Y"),
                              ma = ma(data[ ,type.name[cols]], man = 15))

        ag.data <- aggregate(ma ~ year, data = ma.data, max)

        ag.data$peak <- as.Date("2017-01-01")

        for (r in 1:nrow(ag.data)){
          ag.data$peak[r] <- ma.data$date[which(ma.data$year == ag.data$year[r] & ma.data$ma == ag.data$ma[r])][1]
        }

        ag.data$st <- ag.data$peak - 182
        ag.data$en <- ag.data$peak + 182

        for (d in 1:nrow(ag.data)){
          if (d == 1) {
            date.ps <- seq(from = ag.data$st[d], to = ag.data$en[d], "days")
          } else {
            date.ps <- c(date.ps, seq(from = ag.data$st[d], to = ag.data$en[d], "days"))
          }

        }

        ps.df <- data.frame(date = date.ps, pollen = NA)
        for (p in 1:nrow(ps.df)){
          if (ps.df$date[p] %in% data[ ,1] == TRUE) {ps.df$pollen[p] <- data[data[ ,1] == ps.df$date[p], type.name[cols]]}
        }

        stl.df <- (stl(ts(ps.df$pollen, frequency = 365), s.window = "periodic", na.action=na.approx)$time.series[,1])

        if(min(as.numeric(window(stl.df, start=c(2,1), end=c(2,365)))) < 0){
          pollen.year <- floor(as.numeric(window(stl.df, start=c(2,1), end=c(2,365))) + abs(min(as.numeric(window(stl.df, start=c(2,1), end=c(2,365))))))
        } else {
          pollen.year <- floor(as.numeric(window(stl.df, start=c(2,1), end=c(2,365))))
        }

        par(mfrow = c(ceiling(length(unique(ag.data$year))/4),4), mar = c(0.5,0.5,1,0.5))
        for (y in 1:nrow(ag.data)){
          date.pred <- seq(from = ag.data$st[y], to = ag.data$en[y], "days")
          pollen.sel <- data.frame(date = date.pred, pollen = NA)
          pollen.sel$pollen[pollen.sel$date %in% data[, 1]] <- data[data[ ,1] %in% pollen.sel$date, type.name[cols]]
          pollen.pred <- pollen.year * summary(lm(pollen.sel$pollen ~ pollen.year + 0))$coefficients[1]

          data.full[which(data.full$date %in% date.pred & is.na(data.full[, type.name[cols]])), type.name[cols]] <- pollen.pred[which(date.pred %in% data.full$date[is.na(data.full[, type.name[cols]])])]

          pollen.full <- pollen.sel
          pollen.full$pollen <- NA

          pollen.full$pollen[pollen.full$date %in% data.full[ ,1]] <- data.full[data.full[ ,1] %in% pollen.full$date, type.name[cols]]
          if(plot==TRUE){
          plot(pollen.year, col = "gray", lwd = 2, type = "l", xaxt = "n", yaxt = "n", main = paste(toupper(type.name[cols]), ag.data$year[y]))
          lines(pollen.full$pollen, col = 2)
          lines(pollen.sel$pollen, type = "l", col = 1, lwd = 2)
          legend("topleft", legend = c("original data", "interpolation", "seasonality"), col = c("black", "red", "gray"), pch = 19, cex = 0.8, bty = "n")
        }}
      }

      Interpolated <- gather(data.full, "type", "pollen", colnames(data.full)[-1], factor_key = TRUE)
      Non.interpolated <- gather (data, "type", "pollen", colnames(data.full)[-1], factor_key = TRUE)


      Interpolated$Interpolated <- 1
      Interpolated$Interpolated[which(Interpolated$pollen == Non.interpolated$pollen)] <- 0
      colnames(Interpolated) <- c("Date", "Type", "Pollen", "Interpolated")
      if(result=="long"){
        return(Interpolated)
        print(paste(
          "Process completed: If the data has been interpolated it is marked with 1"
        ))
      }else{return(data.full)}

    }else{


    ######CONVERSION DE DATOS WIDE-LONG
    ConvertData1 <- function(Dataf) {
      cols <- ncol(Dataf)
      Nombres <- colnames(Dataf)
      idvariables <- Nombres[1]
      Tipos <- Nombres[2:cols]
      Result <- gather(Dataf, "Type", "Pollen", Tipos, factor_key = TRUE)
      return(Result)
    }

    Long <- ConvertData1(data)
    Warning <- FALSE
    RawData <- na.omit(Long)
    Interpolated <- rep(0, nrow(RawData))
    RawData <- data.frame(RawData, Interpolated)
    RawData[, 1] <- as.Date(RawData[, 1], origin = "1970-01-01")
    First <- min(RawData[, 1])
    Last <- max(RawData[, 1])
    PollenTypes <- unique(RawData[, 2])
    nPollenTypes <- length(PollenTypes)
    for (a in 1:nPollenTypes) {
      Active_Pollen <- as.character(PollenTypes[a])
      Datatemp <- RawData[RawData[, 2] == Active_Pollen, ]
      Days <- nrow(Datatemp)
      Interpolate <- data.frame()
      for (c in 1:(Days - 1)) {
        if ((Datatemp[c, 1]) != (Datatemp[c + 1, 1] - 1)) {
          Interpolate <-
            rbind(Interpolate, Datatemp[c, 1:4], Datatemp[c + 1, 1:4])
        }
      }
      nInterpolate <- nrow(Interpolate)
      if (nInterpolate >= 2) {
        pares <- seq(1, nInterpolate, by = 2)
        npares <- length(pares)
        #Bucle que va cogiendo todos los pares de fechas que definen el hueco que hay. Va dando saltos de 2 en 2.

        for (d in 1:npares) {
          Active_InterpolStart <-
            pares[d] #Posicion en la que empieza el salto
          Active_InterpolEnd <-
            pares[d] + 1#Posicion en la que acaba el salto
          Datatemp2 <-
            Interpolate[Active_InterpolStart:Active_InterpolEnd, 1:3]#Dataframe que contiene las dos lineas del salto.
          colnames(Datatemp2) <- c("Date", "Type", "Pollen")
          ncasos <-
            as.numeric(Datatemp2[nrow(Datatemp2), 1] - Datatemp2[1, 1])#Numero de dias que faltan
          ncasos <- ncasos - 1 #Numero de datos a introducir.
          Date <- c()
          if (ncasos <= maxdays) {
            #No interpola si los datos que faltan son mas de 30.
            Type <- as.character(c())
            prediccion <- data.frame()
            for (l in 1:ncasos) {
              ##Va creando los datos vacios para rellenar con el metodo
              Active_Date <-
                as.Date((as.numeric(Datatemp2[1, 1]) + l), origin = "1970-01-01")
              Date <-
                as.Date(c(Date, Active_Date), origin = "1970-01-01")
              Type <- as.character(c(Type, Active_Pollen))
            }
            Pollen <- rep(NA, ncasos)
            prediccion <- data.frame(Date, Type, Pollen)
            #######METER AQUi TODOS LOS MeTODOS
            #LINEAL

            if (method == "lineal") {
              modelo <- lm(Pollen ~ Date, data = Datatemp2)#Modelo lineal
              prediccion$Pollen <- predict.lm(modelo, prediccion)

            }

            #MEDIA MoVIL

            if (method == "movingmean") {
              ma1 <- function(arr, n = 10) {
                if (n %% 2 == 0) {
                  n = n + 1
                  try(warning (paste("WARNING! moving average is calculated for n:", n)))
                }
                temp  <-  arr
                for (i in 1:length(arr)) {
                  if (i  <=  (n  -  1)  /  2) {
                    init  <-  1
                  } else{
                    init <- i - (n - 1) / 2
                  }
                  if (i > length(arr) - (n - 1) / 2) {
                    end <- length(arr)
                  } else{
                    end <- i + (n - 1) / 2
                  }

                  temp[i] = mean(arr[init:end], na.rm = T)
                }
                return(temp)
              }

              ftr <- (ncasos + 1) * factor
              Dt <-
                Datatemp[Datatemp[, 1] >= as.Date(Datatemp2[1, 1] - ftr) &
                           Datatemp[, 1] <= as.Date(Datatemp2[2, 1] + ftr), 1:3]
              colnames(Dt) <- colnames(prediccion)
              prediccion[, 3] <- as.numeric(prediccion[, 3])
              Dt <- rbind.data.frame(Dt, prediccion)
              Dt <- Dt[order(Dt[, 1]), ]
              Dt$Interpolated <- ifelse(is.na(Dt[, 3]), 1, 0)
              Dt[, 3] <- ma1(Dt[, 3], ftr)
              solution <- Dt[which(Dt[, 4] == 1), 3]
              prediccion$Pollen <- solution

            }

            #SPLINE
            if (method == "spline") {
              Df <-
                Datatemp[which(Datatemp[, 1] >= (Datatemp2[1, 1] - ndays) &
                                 Datatemp[, 1] <= (Datatemp2[2, 1] + ndays)), ]
              modelo <-
                smooth.spline(x = as.numeric(Df[, 1]),
                              y = Df[, 3],
                              spar = spar)
              prediccion$Pollen <-
                predict(modelo, x = as.numeric(prediccion[, 1]))$y
              prediccion[which(prediccion$Pollen < 0), 3] <- 0
            }

            ########Near Station
            if(method=="neighbour"){
              tryCatch({
                mother<-data[,c(1,which(colnames(data)==Active_Pollen))]
                mother[,1]<-as.Date(mother[,1])
                nmother<-nrow(mother)
                emptiness<-data.frame(as.Date(seq(mother[1,1],mother[nmother,1], by="days")))
                colnames(emptiness)<-colnames(mother)[1]
                motherframe<-merge(emptiness, mother, by.x=c(colnames(emptiness)[1]), all.x = TRUE)
                data2[,1]<-as.Date(data2[,1])
                colnames(data2)[1]<-colnames(data)[1]
                nameneigh<-colnames(data)[1]
                neighbours<-merge(motherframe, data2[,c(nameneigh,Active_Pollen)], by= colnames(Datatemp)[1], all.x=T)
                colnames(neighbours)<-c(colnames(neighbours)[1], "station1", "station2")
                if(!is.null(data3)){
                  data3[,1]<-as.Date(data3[,1])
                  colnames(data3)[1]<-colnames(data)[1]
                  neighbours<-merge(neighbours, data3[,c(nameneigh, Active_Pollen)], by= colnames(Datatemp)[1], all.x=T)
                  colnames(neighbours)<-c(colnames(neighbours)[1], "station1", "station2", "station3")
                }
                if(!is.null(data4)){
                  data4[,1]<-as.Date(data4[,1])
                  colnames(data4)[1]<-colnames(data)[1]
                  neighbours<-merge(neighbours, data4[,c(nameneigh, Active_Pollen)], by= colnames(Datatemp)[1], all.x=T)
                  colnames(neighbours)<-c(colnames(neighbours)[1], "station1", "station2", "station3", "station4")
                }
                if(!is.null(data5)){
                  data5[,1]<-as.Date(data5[,1])
                  colnames(data5)[1]<-colnames(data)[1]
                  neighbours<-merge(neighbours, data5[,c(nameneigh,Active_Pollen)], by= colnames(Datatemp)[1], all.x=T)
                  colnames(neighbours)<-c(colnames(neighbours)[1], "station1", "station2", "station3", "station4", "station5")
                }
                neighbourscomplete<-neighbours#Guardo una copia de los datos sin filtrar
                neighbours[neighbours<=5]<-NA#Quito todos los datos que tengan menos de 5 granos. No son fiables para hacer regresiones
                ncolsneig<-ncol(neighbours)
                for(q in 3:ncolsneig){
                  Correlval<-cor.test(neighbours[,2], neighbours[,q], method = "spearman")#Hago correlaciones para ver si las estaciones se relacionan lo suficiente como para aplicar el metodo
                  if(as.numeric(Correlval$estimate)<mincorr){neighbourscomplete[,-q]}#Quito las correlaciones que no llegan al minimo
                }
                #Si no hay estaciones que tengan suficiente indice de correlacion te imprime un error
                if(ncol(neighbourscomplete)==2){stop(paste("ERROR: None of the neighbour stations have enough correlation with your station to use this method for a mincorr=", mincorr))}
                regresion<-neighbours[,-1]
                ###Hago las regresiones
                regresion<-na.omit(regresion)#Datos para entrenar que sean >5 y sin NA
                npred<-nrow(prediccion)#Numero de datos a predecir
                #Datos originales de todas las estaciones
                gooddata<-neighbourscomplete[which(neighbourscomplete[,1]>=prediccion[1,1] & neighbourscomplete[,1]<=prediccion[npred,1]),]
                columnas<-colnames(regresion)
                ncolumn<-length(columnas)
                for(x in 2:ncolumn){
                  Active_column<-columnas[x]
                  if(anyNA(gooddata[,Active_column])){regresion<-regresion[,-which(colnames(regresion)==Active_column)]}
                }
                if(class(regresion)=="data.frame"){
                  modelo<-lm(station1~.+0, data = regresion)
                  prediccion$Pollen<- predict.lm(modelo, gooddata[,-1])
                }
              }, error = function(e){
                print("Warning: Some of the gaps may not be interpolated due to lack of data")

              })}







            ########FIN MeTODO
            Interpolated <- rep(1, ncasos)
            prediccion <- data.frame(prediccion, Interpolated)
            prediccion<-na.omit(prediccion)
            colnames(prediccion) <- colnames(RawData)
            RawData <- rbind(RawData, prediccion)
          } else{
            Warning <- TRUE
          }
        }
      }

      if (plot == TRUE) {
        Years <- as.numeric(unique(strftime(Datatemp[, 1], "%Y")))
        nyears <- length(Years)
        par(mfrow = c(ceiling(nyears / 4), 4), mar = c(1, 1, 1, 0))
        for (f in 1:nyears) {
          Active_Year <- Years[f]
          Dataplot <-
            RawData[which(as.numeric(strftime(RawData[, 1], "%Y")) == Active_Year &
                            RawData[, 2] == Active_Pollen), ]
          Dataplot <- Dataplot[order(Dataplot[, 1]), ]
          ndataplot<-nrow(Dataplot)
          date<-as.Date(as.Date(Dataplot[1,1]):as.Date(Dataplot[ndataplot,1]), origin = "1970-01-01")
          Emptydata<-data.frame(date, Type=Active_Pollen)
          colnames(Emptydata)<-colnames(Dataplot)[c(1:2)]
          Dataplot1<-merge(Emptydata,
                           Dataplot,
                           by.x = c(colnames(Emptydata)[1], colnames(Emptydata)[2]),
                           all.x = TRUE)
          Dataplot1<-Dataplot1[,1:3]
          Predicted <- Dataplot[which(Dataplot[, 4] == 1), 1:3]
          Prediplot <-
            merge(
              Emptydata,
              Predicted,
              by.x = c(colnames(Predicted)[1], colnames(Predicted)[2]),
              all.x = TRUE
            )

          plot(y=Dataplot1[,3], x=1:nrow(Dataplot1), type = "l", col="black", lwd=1, main = paste(toupper(Active_Pollen) , Active_Year),
               xaxt = "n",
               yaxt = "n")
          lines(
            x = 1:nrow(Dataplot1),
            y = Prediplot[, 3],
            col = "red",
            lwd = 2,
            type = "l"
          )
          points(
            x = 1:nrow(Dataplot1),
            y = Prediplot[, 3],
            col = "red",
            cex = 0.7,
            pch = 19
          )
        }

      }
    }
    colnames(RawData) <- c("Date", "Type", "Pollen", "Interpolated")
    RawData <- RawData[order(RawData$Type, RawData$Date), ]
    Date <- seq.Date(First, Last, by = "day")
    n <- length(Date)
    #Creo una matriz de datos vacia sobre la que meter los otros datos para que no elimine los NA interestacionales
    CompleteData <- data.frame()
    for (a in 1:nPollenTypes) {
      Active_Type <- as.character(PollenTypes[a])
      Type <- rep(Active_Type, n)
      Complet <- cbind.data.frame(Date, Type)
      CompleteData <- rbind.data.frame(CompleteData, Complet)
    }
    colnames(CompleteData) <- colnames(RawData)[c(1, 2)]
    #####OUTPUT
    DataINTER <-
      merge(CompleteData,
            RawData,
            by.x = c(colnames(RawData)[1], colnames(RawData)[2]),
            all = TRUE)
    if (Warning == TRUE){
      warning (paste("WARNING: Gaps with more than", maxdays, "missing data have not been interpolated"))}
    if(result=="long"){
    return(DataINTER)
      print(paste(
        "Process completed: If the data has been interpolated it is marked with 1"
      ))
    }else{
      DataINTER2<-spread(DataINTER[,1:3], "Type", "Pollen")
      return(DataINTER2)}
    }}

