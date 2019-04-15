#' Interactive Plotting Pollen Data (one pollen type).
#'
#' Function to plot the pollen data of a pollen type during several seasons. The plots are fully interactive.
#'
#' @param data A \code{data.frame} object. This \code{data.frame} should include a first column in format \code{Date} and the rest of columns in format \code{numeric} belonging to each pollen type by column.
#' @param pollen A \code{character} string with the name of the particle to show. This \code{character} must match with the name of a column in the input database. This is a mandatory argument.
#' @return An interactive plot of the class \pkg{ggvis}.
#' @seealso \code{\link{iplot_pollen}}
#' @examples data("munich_pollen")
#' @examples iplot_years(data = munich_pollen, pollen = "Betula")
#' @importFrom dplyr filter
#' @importFrom ggvis add_axis axis_props ggvis input_checkboxgroup input_slider layer_lines layer_points scale_numeric %>%
#' @importFrom lubridate is.POSIXt yday year
#' @importFrom tidyr gather
#' @importFrom stats complete.cases
#' @export
iplot_years<-function(data, pollen){

  data<-data.frame(data)
  if(class(data) != "data.frame") stop ("Please include a data.frame: first column with date, and the rest with pollen types")

  if(class(pollen) != "character") stop ("Please include only character values for 'pollen' (including only pollen types in your database)")
  if(class(data[,1])[1]!="Date" & !is.POSIXt(data[,1])) {stop("Please the first column of your data must be the date in 'Date' format")}
  data[,1]<-as.Date(data[,1])
  colnames(data)[1]<-"date"

datalong <-
  gather(
    data,
    colnames(data[2:length(data)]),
    key = "variable",
    value = "value"
  )



## Compare_year
datalong_year <- datalong[which(datalong[, 2] == pollen), ]
datalong_year$year <- as.character(year(datalong_year$date))
datalong_year$date <- yday(datalong_year$date)
datalong_year <-
  datalong_year[which(complete.cases(datalong_year)), ]

iplot_year <- datalong_year %>%
  ggvis(x =  ~ date,
        y =  ~ value,
        stroke = ~ year) %>%
  filter(year %in% eval(input_checkboxgroup(unique(datalong_year$year)))) %>%
  layer_lines(opacity := 0.6, strokeWidth := 1) %>%
  layer_points(size := 20, opacity := 0.8, strokeWidth := 5) %>%
  add_axis("y", title = "Pollen concentration", title_offset = 50) %>%
  scale_numeric("x", domain = input_slider(1, 365, c(1, 365)), clamp = T) %>%
  add_axis("x", title = "Day of the year") %>%
  add_axis("x", orient = "top", ticks = 0, title = paste("Pollen",pollen),
           properties = axis_props(
             axis = list(stroke = "white"),
             labels = list(fontSize = 0)))

return(iplot_year)
}
