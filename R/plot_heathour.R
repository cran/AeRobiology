#'Plotting hourly patterns with heatplot
#'
#'Function to plot pollen data expressed in concentrations with time resolution higher than 1 day (e.g. hourly, bi-hourly concentrations). As heatplot.
#'
#'@param data A \code{data.frame} object with the structure \code{long}. Where the first two columns are \code{factors} indicating the \code{pollen} and the \code{location}. The 3 and 4 columns are \code{POSIXct}, showing the hour. Where the third column is the beguinning of the concentration \code{from} and the fourth column is the end time of the concentrtion data \code{to}. The fifth column shows the concentrations of the different pollen types as \code{numeric}. Please see the example 3-hourly data from the automatic pollen monitor BAA500 from Munich and Viechtach in Bavaria (Germany) \code{data("POMO_pollen")}, supplied by ePIN Network supported by the Bavarian Government.
#'@param locations A \code{logical} object with the specification if the different locations will be displayed in the plot. Argument only used when \code{result == "plot"}. By default, \code{locations = FALSE}.
#'@param low.col A \code{character} object with the specification of the color of the lowest value for the scale. By default, \code{low.col = "blue"}.
#'@param mid.col A \code{character} object with the specification of the color of the medium value for the scale. By default, \code{mid.col = "white"}.
#'@param high.col A \code{character} object with the specification object with the specification of the color of the highest value for the scale. By default, \code{high.col = "red"}.
#'@return The function returns an object or a list of objects of class \pkg{ggplot2}.
#'@references Oteros, J., Pusch, G., Weichenmeier, I., Heimann, U., Mueller, R., Roeseler, S., ... & Buters, J. T. (2015). Automatic and online pollen monitoring. \emph{International archives of allergy and immunology}, 167(3), 158-166.
#'@examples data("POMO_pollen")
#'@examples plot_heathour(POMO_pollen)
#'@importFrom ggplot2 ggplot scale_fill_brewer geom_tile scale_fill_gradientn
#'@importFrom dplyr %>%
#'@export
plot_heathour <-
  function (data,
            locations = FALSE,
            low.col = "blue",
            mid.col = "white",
            high.col = "red") {
    data <- data.frame(data)
    if (class(data) != "data.frame"){
      stop (
        "Please include a data.frame: first column with factor indicating the pollen, second column with factor indicating the locaiton, third column with POSIXct indicating the (from), fourth column with POSIXct indicating the (to) and fifth column with numbers indicating the concentration"
      )}

    if (class(locations) != "logical"){
      stop ("Please include only logical values for locations argument")}
    if (class(low.col) != "character"){
      stop ("Please include only character values for low.col argument, introduce valid color names.")}
    if (class(mid.col) != "character"){
      stop ("Please include only character values for mid.col argument, introduce valid color names.")}
    if (class(high.col) != "character"){
      stop ("Please include only character values for high.col argument, introduce valid color names.")}

    frame3 <- plot_hour(data, result = "table", locations = locations)

    frame3$location <- as.character(frame3$location)

    summaryhour <-
      frame3 %>% group_by(pollen, location, Hour) %>% summarise(percent = mean(percent, na.rm = TRUE))


    if (locations == TRUE) {
      pollenlist2 <- list()
      for (loca in 1:length(unique(summaryhour$location))) {
        locat <- unique(summaryhour$location)[loca]
        temp <- summaryhour[which(summaryhour$location == locat), ]
        plo <- ggplot(temp, aes(x = Hour, y = pollen, fill = percent)) +
          geom_tile() +
          scale_fill_gradientn(colours = c(low.col, mid.col, high.col),
                               limits = c(0, NA)) +
          labs(x = "", y = "", title = locat) +
          theme(axis.text.y = element_text(face = "italic")) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        pollenlist2[[loca]] <- plo
      }
      return(pollenlist2)
    } else{
      temp<-summaryhour
      plott <- ggplot(temp, aes(x = Hour, y = pollen, fill = percent)) +
        geom_tile() +
        scale_fill_gradientn(colours = c(low.col, mid.col, high.col),
                             limits = c(0, NA)) +
        labs(x = "", y = "", title = "") +
        theme(axis.text.y = element_text(face = "italic")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      return(plott)
    }
  }
