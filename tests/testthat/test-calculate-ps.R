context("calculate_ps")

test_that("calculate_ps returns all pollen types", {
  dates <- seq(as.Date("2025-01-01"), as.Date("2025-01-20"), by = "day")
  pollen_data <- data.frame(
    date = dates,
    pollen_a = c(rep(0, 4), 1, 3, 8, 10, 8, 3, 1, rep(0, 9)),
    pollen_b = c(rep(0, 9), 2, 4, 9, 12, 9, 4, 2, rep(0, 4))
  )

  result <- calculate_ps(
    pollen_data,
    method = "percentage",
    perc = 90,
    interpolation = FALSE,
    plot = FALSE,
    export.result = FALSE
  )

  expect_equal(unique(result$type), c("pollen_a", "pollen_b"))
  expect_equal(nrow(result), 2)
})
