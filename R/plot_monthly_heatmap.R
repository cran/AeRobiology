#'Monthly Pollen Heatmap
#'
#'Function to plot monthly sums of daily pollen concentrations as a heatmap. Months
#'are displayed in the x axis and years in the y axis for one selected pollen type.
#'
#'@param data A \code{data.frame} object. The first column must be in \code{Date}
#'format and pollen types must be provided as numeric columns. If the data include
#'several locations, one column can identify the location.
#'@param pollen.type A \code{character} string with the name of the pollen type to
#'be plotted. This name must match one numeric column in \code{data}.
#'@param location A \code{character} string with the location to be selected. If
#'\code{NULL} and a location column is detected, the first location in the data is
#'used. If the data do not include a location column, this argument is ignored.
#'@param location.col A \code{character} string with the name of the column
#'containing the locations. If \code{NULL}, common names such as \code{"location"},
#'\code{"city"} or \code{"station"} are automatically detected.
#'@param start.month A \code{numeric} value ranging \code{1_12} specifying the
#'first month shown in the x axis. The \code{start.month} argument will be
#'\code{1} by default.
#'@param color.A A \code{character} object with the specification of the color of
#'the lowest value for the scale. By default, \code{color.A = "yellow"}.
#'@param color.B A \code{character} object with the specification of the color of
#'the highest value for the scale. By default, \code{color.B = "red"}.
#'@param result A \code{character} object with the definition of the object to be
#'produced by the function. If \code{result == "plot"}, the function returns an
#'object of class \pkg{ggplot2}; if \code{result == "table"}, the function returns
#'a \pkg{data.frame} with the monthly sums. By default, \code{result = "plot"}.
#'@return If \code{result == "plot"}, the function returns an object of class
#'\pkg{ggplot2}; if \code{result == "table"}, the function returns a
#'\pkg{data.frame} with the monthly sums.
#'@examples data("munich_pollen")
#'@examples plot_monthly_heatmap(munich_pollen, pollen.type = "Betula")
#'@examples plot_monthly_heatmap(munich_pollen, pollen.type = "Poaceae", start.month = 6)
#'@importFrom ggplot2 aes element_text geom_tile ggplot labs scale_fill_gradient scale_x_discrete scale_y_continuous theme theme_bw
#'@importFrom lubridate is.POSIXt
#'@export
plot_monthly_heatmap <- function(data,
                                 pollen.type,
                                 location = NULL,
                                 location.col = NULL,
                                 start.month = 1,
                                 color.A = "yellow",
                                 color.B = "red",
                                 result = "plot") {

  data <- data.frame(data)

  if(!is.data.frame(data)) stop ("Please include a data.frame: first column with date, and the rest with pollen types")
  if(missing(pollen.type)) stop ("Please include one pollen type as character value for pollen.type argument")
  if(!is.character(pollen.type) | length(pollen.type) != 1) stop ("Please include one pollen type as character value for pollen.type argument")
  if(!is.null(location) & (!is.character(location) | length(location) != 1)) stop ("Please include one location as character value for location argument")
  if(!is.null(location.col) & (!is.character(location.col) | length(location.col) != 1)) stop ("Please include one column name as character value for location.col argument")
  if(!is.numeric(start.month) | length(start.month) != 1 | start.month < 1 | start.month > 12) stop ("Please include only numeric values between 1-12 for start.month argument")
  if(!is.character(color.A) | length(color.A) != 1) stop ("Please include only character values for color.A argument, introduce valid color names.")
  if(!is.character(color.B) | length(color.B) != 1) stop ("Please include only character values for color.B argument, introduce valid color names.")
  if(result != "plot" & result != "table") stop ("Please result only accept values: 'table' or 'plot'")
  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  if(!pollen.type %in% colnames(data)) stop ("Please pollen.type must match one column name in data")
  if(!is.numeric(data[, pollen.type])) stop ("Please pollen.type must match a numeric column in data")

  tryCatch(grDevices::col2rgb(color.A), error = function(e) stop("Please introduce a valid color name for color.A argument"))
  tryCatch(grDevices::col2rgb(color.B), error = function(e) stop("Please introduce a valid color name for color.B argument"))

  data[,1] <- as.Date(data[,1])

  if(is.null(location.col)){
    location.candidates <- c("location", "Location", "city", "City", "station", "Station", "site", "Site")
    location.col <- location.candidates[location.candidates %in% colnames(data)][1]
  } else if(!location.col %in% colnames(data)){
    stop("Please location.col must match one column name in data")
  }

  selected.location <- NA
  if(!is.na(location.col)){
    data[, location.col] <- as.character(data[, location.col])
    if(is.null(location)){
      selected.location <- unique(data[, location.col])[1]
    } else {
      selected.location <- location
    }
    if(!selected.location %in% data[, location.col]) stop ("Please location must match one value in the location column")
    data <- data[which(data[, location.col] == selected.location), ]
  }

  if(nrow(data) == 0) stop("No data available for the selected location")

  month.order <- c(start.month:12, 1:(start.month - 1))
  if(start.month == 1) month.order <- 1:12

  heat.data <- data.frame(
    date = data[,1],
    pollen = data[, pollen.type],
    year = as.numeric(strftime(data[,1], "%Y")),
    month = as.numeric(strftime(data[,1], "%m"))
  )

  sum_fun <- function(x) {
    if(all(is.na(x))) return(NA)
    sum(x, na.rm = TRUE)
  }

  monthly <- aggregate(pollen ~ year + month, data = heat.data, FUN = sum_fun, na.action = stats::na.pass)
  colnames(monthly) <- c("year", "month", "sum.pollen")
  monthly$month.label <- factor(month.abb[monthly$month], levels = month.abb[month.order])
  monthly$year <- as.numeric(monthly$year)
  monthly <- monthly[order(monthly$year, match(monthly$month, month.order)), ]
  rownames(monthly) <- NULL

  if(!is.na(selected.location)){
    monthly$location <- selected.location
  }

  if(result == "table"){
    return(monthly)
  }

  plot.title <- pollen.type
  if(!is.na(selected.location)){
    plot.title <- paste(pollen.type, selected.location, sep = " - ")
  }

  ggplot(monthly, aes(x = month.label, y = year, fill = sum.pollen)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = color.A, high = color.B, na.value = "grey90") +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = sort(unique(monthly$year))) +
    labs(x = "", y = "", fill = "Monthly sum", title = plot.title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
