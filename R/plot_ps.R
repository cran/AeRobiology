#'Pollen Season Plot
#'
#'Function to plot the main pollen season of a single pollen type.
#'
#'@param data A \code{data.frame} object including the general database where interpollation must be performed. This \code{data.frame} must include a first column in \code{Date} format and the rest of columns in \code{numeric} format. Each column must contain information of one pollen type. It is not necessary to insert missing gaps; the function will automatically detect them.
#'@param pollen.type A \code{character} string specifying the name of the pollen type which will be plotted. The name must be exactly the same that appears in the column name. Mandatory argument with no default.
#'@param year A \code{numeric (interger)} value specifying the season to be plotted. The season does not necessary fit a natural year. See \code{\link{calculate_ps}} for more details. Mandatory argument with no default.
#'@param days A \code{numeric (interger)} specifying the number of days beyond each side of the main pollen season that will be represented. The \code{days} argument will be \code{30} by default.
#'@param fill.col A \code{character} string specifying the name of the color to fill the main pollen season (Galan et al., 2017) in the plot. See \code{\link[ggplot2]{ggplot}} function for more details. The \code{fill.col} argument will be \code{"turquoise4"} by default.
#'@param int.method A \code{character} string specifying the method selected to apply the interpolation method in order to complete the pollen series. The implemented methods that may be used are: \code{"lineal"},  \code{"movingmean"}, \code{"spline"} or  \code{"tseries"}. See \code{\link{interpollen}} function for more details. The \code{int.method} argument will be \code{"lineal"} by default.
#'@param axisname A \code{character} string or an expression specifying the y axis title of the plot. The \code{axisname} argument will be \code{expression(paste("Pollen grains / m" ^ "3"))} by default.
#'@param ... Other arguments passed on to the pollen season calculation as specified in \code{\link{calculate_ps}} function.
#'@details \code{plot_ps} function is designed to easily plot the main pollen season (Galan et al., 2017). The pre_peak period and the post_peak period are marked with different color intensity in the graph. The user must choose a single pollen type and season to plot.
#'@return The function returns an object of class \code{\link[ggplot2]{ggplot}} with a graphical representation of the main pollen season of the selected pollen type. The pre_peak and post_peak periods are marked with different color intensity.
#'@references Galan, C., Ariatti, A., Bonini, M., Clot, B., Crouzy, B., Dahl, A., Fernandez_Gonzalez, D., Frenguelli, G., Gehrig, R., Isard, S., Levetin, E., Li, D.W., Mandrioli, P., Rogers, C.A., Thibaudon, M., Sauliene, I., Skjoth, C., Smith, M., Sofiev, M., 2017. Recommended terminology for aerobiological studies. Aerobiologia (Bologna). 293_295.
#'@seealso \code{\link{calculate_ps}}, \code{\link{interpollen}}, \code{\link[ggplot2]{ggplot}}.
#'@examples data("munich_pollen")
#'@examples plot_ps(munich_pollen, year = 2013, pollen.type = "Betula")
#'@importFrom utils data
#'@importFrom scales date_format
#'@importFrom lubridate is.POSIXt year
#'@importFrom ggplot2 aes element_blank element_text geom_area geom_line ggplot labs scale_x_date theme theme_classic ylab
#'@importFrom graphics plot
#'@importFrom grDevices dev.off png
#'@importFrom tidyr %>%
#'@export
plot_ps <- function(data,
           pollen.type = NULL ,
           year = NULL ,
           days = 30,
           fill.col = "turquoise4",
           int.method="lineal",
           axisname= expression(paste("Pollen grains / m" ^ "3")),
           ...) {
    if(class(axisname)!="expression" & class(axisname)!="character"){stop("axisname: Please, insert only a character string or an expression")}
    if (class(fill.col) != "character") {
      stop("fill.col: Please, insert only a character string defining an existing color")
    }
  data<-data.frame(data)
    if (class(data) != "data.frame" &
        !is.null(data)) {
      stop ("Please include a data.frame: first column with date, and the rest with pollen types")
    }
    if (is.null(year)) {
      stop("Please, select a year to plot")
    }
    if (is.null(pollen.type)) {
      stop("Please, select a pollen.type to plot")
    }
    if (class(pollen.type) != "character") {
      stop("pollen.type: Please, insert only a character string")
    }
    if (class(year) != "numeric") {
      stop("year: Please, insert only a number")
    }
    if (class(days) != "numeric") {
      stop("days: Please, insert only a number bigger than 1")
    }
    if (days < 1 | days %% 1 != 0 ) {
      stop("days: Please, insert only an entire number bigger than 1")
    }
    if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
    data[,1]<-as.Date(data[,1])
    namecolumn <- colnames(data)[-1]
    years <- unique(year(data[, 1]))
    if (!(pollen.type %in% namecolumn)) {
      stop(
        "pollen.type: Please, insert only a pollen type which is in your database (type the exact name)"
      )
    }
    if (!(year %in% years)) {
      stop("pollen.type: Please, insert only a year which is in your database")
    }
    datename <- colnames(data)[1]
    dataframe <- data[, c(datename, pollen.type)]
    referencetable <-
      calculate_ps(data = dataframe,
                   plot = FALSE,
                   export.result = FALSE,
                   int.method = int.method,
                   interpolation = TRUE, ...)
    dataframe <- interpollen(data = dataframe, method=int.method, plot = FALSE)

    Start <-
      referencetable[which(referencetable$seasons == year), 3] - (days)
    End <-
      referencetable[which(referencetable$seasons == year), 5] + (days)

    StartMPS <-
      referencetable[which(referencetable$seasons == year), 3]
    EndMPS <-
      referencetable[which(referencetable$seasons == year), 5]
    Peak <-
      referencetable[which(referencetable$seasons == year), 11]

    dataplot <-
      dataframe[which(dataframe[, 1] >= Start &
                        dataframe[, 1] <= End),]
    colnames(dataplot) <- c("date", "pollen")
    Preframe <- dataplot
    Postframe <- dataplot
    Preframe[which(Preframe$date < StartMPS |
                     Preframe$date > Peak), 2] <- NA
    Postframe[which(Postframe$date < Peak |
                      Preframe$date > EndMPS), 2] <- NA
    graph <- ggplot() +
      geom_area(data = dataplot,
                aes(date, pollen),
                color = "grey90",
                fill = "grey90") +
      geom_area(data = Preframe, aes(date, pollen), fill = fill.col) +
      geom_area(data = Postframe,
                aes(date, pollen),
                fill = fill.col,
                alpha = 0.5) +
      geom_line(data = dataplot,
                aes(date, pollen),
                color = "grey10",
                size = 0.3) +
      scale_x_date(labels = date_format("%d-%b"), breaks = '7 days') +
      theme_classic() +
      labs(title = paste(pollen.type, year)) +
      ylab(axisname) +
      theme(
        plot.title = element_text(face = "bold.italic", size = 16),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    return(graph)
  }
