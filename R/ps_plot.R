#'Pollen Season Plot
#'
#'Function to plot the main pollen season of a single pollen type.
#'
#'@param data A \code{data.frame} object including the general database where interpollation must be performed. This \code{data.frame} must include a first column in \code{Date} format and the rest of columns in \code{numeric} format. Each column must contain information of one pollen type. It is not necessary to insert missing gaps; the function will automatically detect them.
#'@param pollen.type A \code{character} string specifying the name of the pollen type which will be plotted. The name must be exactly the same that appears in the column name. Mandatory argument with no default.
#'@param year A \code{numeric (interger)} value specifying the season to be plotted. The season does not necessary fit a natural year. See \code{\link{calculate_ps}} for more details. Mandatory argument with no default.
#'@param export.plot  A \code{logical} argument. If \code{TRUE}, the graph created will be exported in a new folder (\emph{"plot_AeRobiology"}) in your working directory. The \code{export.plot} argument will be \code{TRUE} by default.
#'@param export.format A \code{character} string specifying the extension of the file in which the graph will be saved. The implemented extensions that can be used are: \code{"pdf", "png", "bmp", "jpeg"} and \code{"tiff"}. See \code{\link[ggplot2]{ggsave}} function for more details. Only valid if \code{export.plot = TRUE}. The \code{export.format} argument will be \code{"pdf"} by default.
#'@param days A \code{numeric (interger)} specifying the number of days beyond each side of the main pollen season that will be represented. The \code{days} argument will be \code{30} by default.
#'@param fill.col A \code{character} string specifying the name of the color to fill the main pollen season (Galan et al., 2017) in the plot. See \code{\link[ggplot2]{ggplot}} function for more details. The \code{fill.col} argument will be \code{"turquoise4"} by default.
#'@param export.width A \code{numeric (double)} value specifying the width of the graph in inches. Only valid if \code{export.plot = TRUE}. See \code{\link[ggplot2]{ggsave}} function for more details. The \code{export.width} argument will be \code{6} by default.
#'@param export.height A \code{numeric (double)} value specifying the height of the graph in inches. Only valid if \code{export.plot = TRUE}. See \code{\link[ggplot2]{ggsave}} function for more details. The \code{export.width} argument will be \code{6} by default.
#'@param int.method A \code{character} string specifying the method selected to apply the interpolation method in order to complete the pollen series. The implemented methods that may be used are: \code{"lineal"},  \code{"movingmean"}, \code{"spline"} or  \code{"tseries"}. See \code{\link{interpollen}} function for more details. The \code{int.method} argument will be \code{"lineal"} by default.
#'@param axisname A \code{character} string or an expression specifying the y axis title of the plot. The \code{axisname} argument will be \code{expression(paste("Pollen grains / m" ^ "3"))} by default.
#'@param ... Other arguments passed on to the pollen season calculation as specified in \code{\link{calculate_ps}} function.
#'@details \code{ps_plot} function is designed to easily plot the main pollen season (Galan et al., 2017). The pre_peak period and the post_peak period are marked with different color intensity in the graph. The user must choose a single pollen type and season to plot.
#'@return The function returns an object of class \code{\link[ggplot2]{ggplot}} with a graphical representation of the main pollen season of the selected pollen type. The pre_peak and post_peak periods are marked with different color intensity. If \code{export.plot = TRUE}, the graph will be exported into \emph{"plot_AeRobiology"} folder (created in the working directory). The extension of the file will be the specified by \code{export.format} argument.
#'@references Galan, C., Ariatti, A., Bonini, M., Clot, B., Crouzy, B., Dahl, A., Fernandez_Gonzalez, D., Frenguelli, G., Gehrig, R., Isard, S., Levetin, E., Li, D.W., Mandrioli, P., Rogers, C.A., Thibaudon, M., Sauliene, I., Skjoth, C., Smith, M., Sofiev, M., 2017. Recommended terminology for aerobiological studies. Aerobiologia (Bologna). 293_295.
#'@seealso \code{\link{calculate_ps}}, \code{\link{interpollen}}, \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{ggsave}}
#'@examples data("munich")
#'@examples ps_plot(munich, year = 2013, pollen.type = "Betula", export.plot = FALSE)
#'@importFrom utils data
#'@importFrom scales date_format
#'@importFrom lubridate is.POSIXt year
#'@importFrom ggplot2 aes element_blank element_text geom_area geom_line ggplot ggsave labs scale_x_date theme theme_classic ylab
#'@importFrom graphics plot
#'@importFrom grDevices dev.off png
#' @importFrom tidyr %>%
#'@export
ps_plot <- function(data,
           pollen.type = NULL ,
           year = NULL ,
           export.plot = TRUE,
           export.format = "pdf",
           days = 30,
           fill.col = "turquoise4",
           export.width = 6,
           export.height = 6,
           int.method="lineal",
           axisname= expression(paste("Pollen grains / m" ^ "3")),
           ...) {
    if(class(axisname)!="expression" & class(axisname)!="character"){stop("axisname: Please, insert only a character string or an expression")}
    if (class(fill.col) != "character") {
      stop("fill.col: Please, insert only a character string defining an existing color")
    }
    if (class(data) != "data.frame" &
        !is.null(data)) {
      stop ("Please include a data.frame: first column with date, and the rest with pollen types")
    }
    if (is.null(year)) {
      stop("Please, select a year to plot")
    }
    if (class(export.format) != "character") {
      stop(
        "export.format: Please, insert a character string with the file format (pdf, png, bmp, jpeg, ...)"
      )
    }
    if (is.null(pollen.type)) {
      stop("Please, select a pollen.type to plot")
    }
    if (class(export.width) != "numeric") {
      stop("export.width: Please, insert only a number")
    }
    if (class(export.height) != "numeric") {
      stop("export.height: Please, insert only a number")
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
    if (export.plot == TRUE) {
      ifelse(!dir.exists(file.path("plot_AeRobiology")), dir.create(file.path("plot_AeRobiology")), FALSE)

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
    if (export.plot == TRUE) {
      graph
      ggsave(
        paste("plot_AeRobiology/ps_plot_",pollen.type, year,".",export.format, sep = ""),
        plot = graph,
        width = export.width,
        height = export.height
      )
      png(paste0("plot_AeRobiology/credits.png"))

      dev.off()
    }
    return(graph)
  }
