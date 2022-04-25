test_that("top_results() works faultlessly.", {
  # load pre-computed results (obtaines via `distinct_test`)
  data("res", package = "distinct")
  
  # Visualize significant results:
  top_1 = top_results(res)

  # Visualize significant results:
  top_2 = top_results(res, significance = 0.1)

  # Visualize significant results:
  top_3 = top_results(res, global = FALSE)
    
  # Visualize significant results from a specified cluster of cells:
  top_4 =  top_results(res, cluster = "Dendritic cells")
  
  expect_is(top_1, "data.frame")
  expect_is(top_2, "data.frame")
  expect_is(top_3, "data.frame")
  expect_is(top_4, "data.frame")
})
