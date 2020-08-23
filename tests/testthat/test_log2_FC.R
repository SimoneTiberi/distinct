test_that("log2_FC() works faultlessly.", {
  # load pre-computed results (obtaines via `distinct_test`)
  data("res", package = "distinct")
  # load the input data:
  data("Kang_subset", package = "distinct")
  
  # We can optionally add the fold change (FC) and log2-FC between groups:
  res = log2_FC(res = res,
                x = Kang_subset, 
                name_assays_expression = "cpm",
                name_group = "stim",
                name_cluster = "cell")
  
  expect_is(res, "data.frame")
})
