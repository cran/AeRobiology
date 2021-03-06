% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heathour.R
\name{plot_heathour}
\alias{plot_heathour}
\title{Plotting hourly patterns with heatplot}
\usage{
plot_heathour(data, locations = FALSE, low.col = "blue",
  mid.col = "white", high.col = "red")
}
\arguments{
\item{data}{A \code{data.frame} object with the structure \code{long}. Where the first two columns are \code{factors} indicating the \code{pollen} and the \code{location}. The 3 and 4 columns are \code{POSIXct}, showing the hour. Where the third column is the beguinning of the concentration \code{from} and the fourth column is the end time of the concentrtion data \code{to}. The fifth column shows the concentrations of the different pollen types as \code{numeric}. Please see the example 3-hourly data from the automatic pollen monitor BAA500 from Munich and Viechtach in Bavaria (Germany) \code{data("POMO_pollen")}, supplied by ePIN Network supported by the Bavarian Government.}

\item{locations}{A \code{logical} object with the specification if the different locations will be displayed in the plot. Argument only used when \code{result == "plot"}. By default, \code{locations = FALSE}.}

\item{low.col}{A \code{character} object with the specification of the color of the lowest value for the scale. By default, \code{low.col = "blue"}.}

\item{mid.col}{A \code{character} object with the specification of the color of the medium value for the scale. By default, \code{mid.col = "white"}.}

\item{high.col}{A \code{character} object with the specification object with the specification of the color of the highest value for the scale. By default, \code{high.col = "red"}.}
}
\value{
The function returns an object or a list of objects of class \pkg{ggplot2}.
}
\description{
Function to plot pollen data expressed in concentrations with time resolution higher than 1 day (e.g. hourly, bi-hourly concentrations). As heatplot.
}
\examples{
data("POMO_pollen")
plot_heathour(POMO_pollen)
}
\references{
Oteros, J., Pusch, G., Weichenmeier, I., Heimann, U., Mueller, R., Roeseler, S., ... & Buters, J. T. (2015). Automatic and online pollen monitoring. \emph{International archives of allergy and immunology}, 167(3), 158-166.
}
