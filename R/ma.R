#' Moving Average Calculator
#'
#' Function to calculate the moving average of a given numeric vector. The order of the moving average may be customized by the user (\code{man} argument).
#'
#' @param data A \code{numeric} vector (e.g. a \code{numeric} column of a \code{dataframe}).
#' @param man An \code{integer} value specifying the order of the moving average applied to the data. By default, \code{man = 10}.
#' @param warnings A \code{logical} value specifying the show of warning messages. By default, \code{warnings = FALSE}.
#' @return This function returns a vector with the moving average of the input data.
#' @examples data("munich_pollen")
#' @examples ma(data = munich_pollen$Betula, man = 10, warnings = FALSE)
#' @export
ma <- function(data, man = 10, warnings = FALSE) {
  if (man %% 2 == 0) {
    man = man + 1
    try(if (warnings == TRUE)
      warning (paste("WARNING! moving average is calculated for man:", man)))
  }
  temp  <-  data
  for (i in 1:length(data)) {
    if (i  <=  (man  -  1)  /  2) {
      init  <-  1
    } else{
      init <- i - (man - 1) / 2
    }
    if (i > length(data) - (man - 1) / 2) {
      end <- length(data)
    } else{
      end <- i + (man - 1) / 2
    }

    temp[i] = mean(data[init:end], na.rm = T)
  }
  return(temp)
}
