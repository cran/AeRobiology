#' Interactive Plotting Pollen Data (one season).
#'
#' Function to plot the pollen data during one season. The plots are fully interactive.
#'
#' @param data A \code{data.frame} object. This \code{data.frame} should include a first column in format \code{Date} and the rest of columns in format \code{numeric} belonging to each pollen type by column.
#' @param year An \code{integer} value specifying the year to display. This is a mandatory argument.
#' @return An interactive plot of the class \pkg{ggvis}.
#' @seealso \code{\link{iplot_years}}
#' @examples data("munich_pollen")
#' @examples iplot_pollen(data = munich_pollen, year = 2012)
#' @importFrom dplyr filter
#' @importFrom ggvis add_axis axis_props ggvis input_checkboxgroup input_slider layer_lines layer_points scale_numeric %>%
#' @importFrom lubridate is.POSIXt yday year
#' @importFrom tidyr gather
#' @importFrom stats complete.cases
#' @export
iplot_pollen<-function(data, year){

  data<-data.frame(data)
  if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

  if(class(year) != "numeric") stop ("Please include only numeric values for 'year' (including only years in your database)")
  if(class(data[,1])!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  colnames(data)[1]<-"date"
  data[,1]<-as.Date(data[,1])

  datalong <-
    gather(
      data,
      colnames(data[2:length(data)]),
      key = "variable",
      value = "value"
    )

## Compare_pollen
datalong_pollen <- datalong[which(year(datalong[, 1]) == year), ]
datalong_pollen$date <- yday(datalong_pollen$date)
datalong_pollen <-
  datalong_pollen[which(complete.cases(datalong_pollen)), ]

iplot_pollen <- datalong_pollen %>%
  ggvis(x =  ~ date,
        y =  ~ value,
        stroke = ~ variable) %>%
  filter(variable %in% eval(input_checkboxgroup(unique(
    datalong_pollen$variable
  )))) %>%
  layer_lines(opacity := 0.6, strokeWidth := 1) %>%
  layer_points(size := 20, opacity := 0.8, strokeWidth := 5) %>%
  add_axis("y", title = "Pollen concentration", title_offset = 50) %>%
  scale_numeric("x", domain = input_slider(1, 365, c(1, 365)), clamp = T) %>%
  add_axis("x", title = "Day of the year") %>%
  add_axis("x", orient = "top", ticks = 0, title = paste("Year",year),
           properties = axis_props(
             axis = list(stroke = "white"),
             labels = list(fontSize = 0)))

return(iplot_pollen)

}
