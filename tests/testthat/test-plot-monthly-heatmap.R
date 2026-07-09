context("plot_monthly_heatmap")

test_that("plot_monthly_heatmap aggregates monthly sums for selected location", {
  dates <- seq(as.Date("2025-01-01"), as.Date("2025-02-03"), by = "day")
  heat_data <- rbind(
    data.frame(date = dates, location = "A", Betula = 1),
    data.frame(date = dates, location = "B", Betula = 2)
  )

  default_table <- plot_monthly_heatmap(
    heat_data,
    pollen.type = "Betula",
    result = "table"
  )

  selected_table <- plot_monthly_heatmap(
    heat_data,
    pollen.type = "Betula",
    location = "B",
    result = "table"
  )

  expect_true(all(default_table$location == "A"))
  expect_true(all(selected_table$location == "B"))
  expect_equal(default_table$sum.pollen[default_table$month == 1], 31)
  expect_equal(default_table$sum.pollen[default_table$month == 2], 3)
  expect_equal(selected_table$sum.pollen[selected_table$month == 1], 62)
  expect_equal(selected_table$sum.pollen[selected_table$month == 2], 6)
})

test_that("plot_monthly_heatmap uses the selected first month", {
  dates <- seq(as.Date("2025-01-01"), as.Date("2025-12-31"), by = "day")
  heat_data <- data.frame(date = dates, Betula = 1)

  result <- plot_monthly_heatmap(
    heat_data,
    pollen.type = "Betula",
    start.month = 6,
    result = "table"
  )

  expect_equal(as.character(result$month.label[1]), "Jun")
  expect_equal(result$year[result$month == 1], 2025)
})
