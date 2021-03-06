% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_hour.R
\name{plot_hour}
\alias{plot_hour}
\title{Plotting hourly patterns}
\usage{
plot_hour(data, result = "plot", locations = FALSE)
}
\arguments{
\item{data}{A \code{data.frame} object with the structure \code{long}. Where the first two columns are \code{factors} indicating the \code{pollen} and the \code{location}. The 3 and 4 columns are \code{POSIXct}, showing the hour. Where the third column is the beguinning of the concentration \code{from} and the fourth column is the end time of the concentrtion data \code{to}. The fifth column shows the concentrations of the different pollen types as \code{numeric}. Please see the example 3-hourly data from the automatic pollen monitor BAA500 from Munich and Viechtach in Bavaria (Germany) \code{data("POMO_pollen")}, supplied by ePIN Network supported by the Bavarian Government.}

\item{result}{A \code{character} object with the definition of the object to be produced by the function. If \code{result == "plot"}, the function returns a list of objects of class \pkg{ggplot2}; if \code{result == "table"}, the function returns a \pkg{data.frame} with the hourly patterns. By default, \code{result = "plot"}.}

\item{locations}{A \code{logical} object with the specification if the different locations will be displayed in the plot. Argument only used when \code{result == "plot"}. By default, \code{locations = FALSE}.}
}
\value{
If \code{result == "plot"}, the function returns a list of objects of class \pkg{ggplot2}; if \code{result == "table"}, the function returns a \pkg{data.frame} with the hourly patterns.
}
\description{
Function to plot pollen data expressed in concentrations with time resolution higher than 1 day (e.g. hourly, bi-hourly concentrations).
}
\examples{
data("POMO_pollen")
plot_hour(POMO_pollen, result="plot", locations = FALSE)
plot_hour(POMO_pollen, result="plot", locations = TRUE)
}
\references{
Oteros, J., Pusch, G., Weichenmeier, I., Heimann, U., Mueller, R., Roeseler, S., ... & Buters, J. T. (2015). Automatic and online pollen monitoring. \emph{International archives of allergy and immunology}, 167(3), 158-166.
}
