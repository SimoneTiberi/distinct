test_that("plot_densities() works faultlessly.", {
  data("Kang_subset", package = "distinct")
  Kang_subset
  
  res = plot_densities(x = Kang_subset,
            gene = "ISG15",
            cluster = "Dendritic cells",
            name_assays_expression = "logcounts",
            name_cluster = "cell",
            name_sample = "sample_id",
            name_group = "stim")
  
  expect_is(res, "ggplot")
})
