test_that("test_DTU() works faultlessly.", {
  data("Kang_subset", package = "distinct")

  set.seed(61217)
  res = distinct_test(
    x = Kang_subset, 
    name_assays_expression = "counts",
    name_cluster = "cell",
    name_group = "stim",
    name_sample = "sample_id",
    P = 10^2, 
    min_non_zero_cells = 20)
  
  expect_is(res, "data.frame")
})
