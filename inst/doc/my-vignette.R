## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, fig.width = 7)
Sys.setlocale(locale="english")

## ----echo=FALSE----------------------------------------------------------
library ("knitr")

## ------------------------------------------------------------------------
library ("AeRobiology")

## ----eval=FALSE, echo = TRUE---------------------------------------------
#  install.packages("AeRobiology")
#  library (AeRobiology)

## ----echo = TRUE---------------------------------------------------------
data("munich_pollen")

## ----eval=FALSE, echo = TRUE---------------------------------------------
#  install.packages("readxl")
#  library (readxl)

## ----eval=FALSE, echo = TRUE---------------------------------------------
#  Mydata<-read_xlsx("C:/Users/Antonio/Desktop/Prueba Markdown/mydata.xlsx")
#  

## ----eval=FALSE, echo = TRUE---------------------------------------------
#  Mydata<-read_xlsx("C:/Users/Antonio/Desktop/Prueba Markdown/mydata.xlsx", sheet=2)
#  

## ----echo=TRUE, results='hold'-------------------------------------------
str(munich_pollen)

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
QualityControl<-quality_control(munich_pollen, result = "table")

## ----echo=TRUE, results='hold'-------------------------------------------

head(QualityControl)

## ----echo=TRUE,  fig.keep='first', results='hide'------------------------
quality_control(munich_pollen, result = "plot")

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
quality_control(munich_pollen, int.window = 4, perc.miss = 50, ps.method = "percentage", perc = 80)

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
Interpolated<-interpollen(munich_pollen[,c(1,6)], method="lineal", plot = TRUE, result = "long")


## ----echo=TRUE, results='hold'-------------------------------------------
head(Interpolated)

## ----echo=TRUE, results='hide'-------------------------------------------
CompleteData<-interpollen(munich_pollen[,c(1,6)], method="lineal", plot = FALSE, result = "wide")

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  calculate_ps(munich_pollen, method="percentage", interpolation=TRUE, int.method = "lineal", plot = F)

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
i<-interpollen(munich_pollen[,c(1,6)], method="movingmean", factor = 2, plot = TRUE)

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
i2<-interpollen(munich_pollen[,c(1,6)], method="spline", ndays=3, spar=0.7, plot = FALSE)

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
i3<-interpollen(munich_pollen[,c(1,6)], method="spline", ndays=5, spar=0.2, plot = TRUE)

## ----echo=TRUE, results='hide', fig.keep='last'--------------------------
i4<-interpollen(munich_pollen, method="tseries", plot = TRUE)

## ----echo = TRUE, results='hide', fig.keep='last', warning=FALSE---------
pollen_season <- calculate_ps(munich_pollen)

## ----echo = FALSE, fig.keep='all', warning=FALSE-------------------------
knitr::kable(pollen_season[24:31, ] , format = "html", booktabs = TRUE)

## ----echo = TRUE, fig.keep='last', warning=FALSE-------------------------
calculate_ps(munich_pollen[,c(1,6)], plot = TRUE)

## ----echo=TRUE, results='hide', fig.keep='first', eval=FALSE-------------
#  calculate_ps(munich_pollen[,c(1,6)], method = "percentage", perc = 90, export.result = TRUE)

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
calculate_ps(munich_pollen[,c(1,6)], method = "percentage", perc = 90, export.result = FALSE, int.method = "spline")

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
calculate_ps(munich_pollen[,c(1,6)], method = "percentage", perc = 75, export.result = FALSE, interpolation = FALSE)

## ----echo=TRUE, results='hide', fig.keep='first', warning=FALSE----------
pollen_season<-calculate_ps(munich_pollen[,c(1,6)], method = "logistic", derivative = 5, reduction=FALSE)

## ----echo = FALSE, fig.keep='all', warning=FALSE-------------------------
knitr::kable(pollen_season, format = "html", booktabs = TRUE )

## ----echo=TRUE, results='hide', fig.keep='first', warning=FALSE----------
pollen_season<-calculate_ps(munich_pollen[,c(1,6)], method = "logistic", derivative = 6, reduction=FALSE, red.level = 0.8)

## ----echo = FALSE, fig.keep='all', warning=FALSE-------------------------
knitr::kable(pollen_season, format = "html", booktabs = TRUE)

## ----echo=TRUE, results='hide', fig.keep='first', warning=FALSE----------
pollen_season<-calculate_ps(munich_pollen[,c(1,6)], method = "clinical", type = "birch")

## ----echo = FALSE, fig.keep='all', warning=FALSE-------------------------
knitr::kable(pollen_season, format = "html", booktabs = TRUE)

## ----echo=TRUE, results='hide', fig.keep='first', warning=FALSE----------
pollen_season<-calculate_ps(munich_pollen[,c(1,6)], method = "clinical", n.clinical = 5, window.clinical = 7, th.pollen = 10, th.sum = 100, th.day = 100)

## ----echo = FALSE, fig.keep='all', warning=FALSE-------------------------
knitr::kable(pollen_season, format = "html", booktabs = TRUE)

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
calculate_ps(munich_pollen[,c(1,6)], method = "grains", window.grains = 3, th.pollen = 2 )

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
calculate_ps(munich_pollen[,c(1,6)], method = "moving", man = 7, th.ma = 4)

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
calculate_ps(munich_pollen[,c(1,3)], method = "moving", man = 7, th.ma = 4, def.season = "interannual")

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
calculate_ps(munich_pollen[,c(1,6)], method = "percentage", perc=95, def.season = "peak")

## ----echo=TRUE, results='hide', fig.keep='first'-------------------------
CompleteData<-interpollen(munich_pollen, method="spline", ndays=3, spar=0.7, plot = TRUE, maxdays = 3, result = "wide")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_ps(CompleteData, pollen.type="Alnus", year=2013, fill.col = "orange", axisname = "AeRobiology custom units")

## ----echo = TRUE, results='hide',fig.keep='all', warning=FALSE-----------
pollen_calendar(munich_pollen, method = "heatplot", period = "daily", color = "green", interpolation = FALSE)

## ----echo = TRUE, fig.keep='all', warning=FALSE--------------------------
average_values<-pollen_calendar(munich_pollen, method = "heatplot", period = "daily", color = "green", interpolation = FALSE, result = "table")
knitr::kable(average_values[82:90, ], format = "html", booktabs = TRUE)

## ----echo = TRUE, results='hide',fig.keep='all', warning=FALSE-----------
pollen_calendar(data = munich_pollen, method = "heatplot", period = "daily", color = "red", method.classes = "custom", n.classes = 5, classes = c(5, 25, 50, 200), interpolation = FALSE)

## ----echo = TRUE, results='hide',fig.keep='all', warning=FALSE-----------
pollen_calendar(data = munich_pollen, method = "heatplot", period = "daily", color = "purple", method.classes = "custom", n.classes = 5, classes = c(5, 25, 50, 200), start.month = 11, na.remove = FALSE, interpolation = FALSE)

## ----echo = TRUE, results='hide',fig.keep='all', warning=FALSE-----------
pollen_calendar(data = munich_pollen, method = "heatplot", period = "weekly", color = "blue", method.classes = "exponential", n.types = 4, y.start = 2011, y.end = 2014, interpolation = FALSE)

## ----echo = TRUE, results='hide',fig.keep='all', warning=FALSE-----------
pollen_calendar(data = munich_pollen, method = "phenological", n.types = 5, y.start = 2011, y.end = 2014, interpolation = FALSE)

## ----echo = TRUE, results='hide',fig.keep='all', warning=FALSE-----------
pollen_calendar(data = munich_pollen,  method = "phenological", perc1 = 90, perc2 = 95, th.pollen = 5, interpolation = FALSE)

## ----echo = TRUE, results='hide', fig.keep='all', warning=FALSE----------
pollen_calendar(data = munich_pollen, method = "violinplot", y.start = 2012, y.end = 2015, interpolation = FALSE)

## ----echo = TRUE, results='hide', fig.keep='all', warning=FALSE----------
pollen_calendar(data = munich_pollen, method = "violinplot", th.pollen = 10, na.rm = FALSE, interpolation = FALSE)

## ----echo = TRUE, results='hide',fig.keep='all', eval=FALSE--------------
#  iplot_pollen(munich_pollen, year = 2012)
#  iplot_years(munich_pollen, pollen = "Betula")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_summary(munich_pollen, pollen = "Betula")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_summary(munich_pollen, pollen = "Betula", mave = 5)

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_summary(munich_pollen, pollen = "Betula", mave = 5, normalized = TRUE)

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_normsummary(munich_pollen, pollen = "Betula", color.plot = "red")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_normsummary(munich_pollen, pollen = "Betula", color.plot = "green", mave = 5, normalized = TRUE)

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
analyse_trend(munich_pollen, interpolation = FALSE, export.result = FALSE, export.plot = FALSE, result="plot")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
analyse_trend(munich_pollen, interpolation = FALSE, export.result = FALSE, export.plot = FALSE, split = FALSE, quantil = 1, result="plot")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
analyse_trend(munich_pollen, interpolation = FALSE, export.result = FALSE, export.plot = FALSE, split=FALSE, quantil = 0.5, result="plot")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
analyse_trend(munich_pollen, interpolation = FALSE, export.result = FALSE, export.plot = FALSE, significant = 1, result = "plot")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
analyse_trend(munich_pollen, interpolation = FALSE, export.result = FALSE, export.plot = FALSE, significant = 1, result="table")

## ----echo = TRUE, results='hide',fig.keep='all', eval=FALSE--------------
#  plot_trend(munich_pollen, interpolation = FALSE, export.plot = TRUE, export.result = TRUE)

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
iplot_abundance(munich_pollen, interpolation = FALSE, export.plot = FALSE, export.result = FALSE)

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
iplot_abundance(munich_pollen, interpolation = FALSE, export.plot = FALSE, export.result = FALSE, n.types = 3)

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
iplot_abundance(munich_pollen, interpolation = FALSE, export.plot = FALSE, export.result = FALSE, n.types = 3, y.start = 2011, y.end = 2011, col.bar = "#d63131")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
iplot_abundance(munich_pollen, interpolation = FALSE, export.plot = FALSE, export.result = FALSE, n.types = 3, y.start = 2011, y.end = 2011, col.bar = "#d63131", type.plot = "dynamic")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
iplot_pheno(munich_pollen, method= "percentage", perc=80, int.method="spline", n.types = 8)

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
iplot_pheno(munich_pollen, method= "clinical", n.clinical = 3, int.method="spline", n.types = 4)

## ----echo = TRUE, results='hide',fig.keep='all', eval=FALSE--------------
#  iplot_pheno(munich_pollen, method= "clinical", n.clinical = 3, int.method="spline", n.types = 4, type.plot = "dynamic")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_ps(munich_pollen, pollen.type="Alnus", year=2011)

## ----echo = TRUE, results='hold', error=TRUE-----------------------------
plot_ps(munich_pollen, pollen.type="Alnuscdscscr", year=2011)

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_ps(munich_pollen, pollen.type="Alnus", year=2013, method= "percentage", perc=95 ,int.method = "lineal")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_ps(munich_pollen, pollen.type="Alnus", year=2013, method= "percentage", perc=95 ,int.method = "movingmean")

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_ps(munich_pollen, pollen.type="Alnus", year=2013, days = 90)

## ----echo = TRUE, results='hide',fig.keep='all'--------------------------
plot_ps(munich_pollen, pollen.type="Alnus", year=2013, fill.col = "orange", axisname = "AeRobiology custom units")

## ----echo = FALSE, warning=FALSE, message=FALSE--------------------------
library(ggplot2)
library(dplyr)

## ----echo = TRUE---------------------------------------------------------
data("POMO_pollen")

## ----echo = TRUE, message=FALSE------------------------------------------
plot_hour(POMO_pollen)

## ----message=FALSE-------------------------------------------------------
TO<-plot_hour(POMO_pollen, result ="table")
knitr::kable(TO[1:10,], caption = "3-Hourly patterns", row.names = FALSE, digits = 1, format = "html", booktabs = TRUE)

## ---- message=FALSE, echo=TRUE-------------------------------------------
plot_hour(POMO_pollen, locations = TRUE)

## ---- message=FALSE, echo=TRUE-------------------------------------------
plot_heathour(POMO_pollen)

## ---- message=FALSE, echo=TRUE-------------------------------------------
plot_heathour(POMO_pollen, low.col = "darkgreen", mid.col = "moccasin", high.col = "brown")

## ---- message=FALSE, echo=TRUE-------------------------------------------
plot_heathour(POMO_pollen, locations = TRUE)

