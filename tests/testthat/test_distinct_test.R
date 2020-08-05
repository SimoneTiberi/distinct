test_that("distinct_test() works faultlessly.", {
  data("Kang_subset", package = "distinct")

  # create the design of the study:
  samples = metadata(Kang_subset)$experiment_info$sample_id
  group = metadata(Kang_subset)$experiment_info$stim
  design = model.matrix(~group)
  # rownames of the design must indicate sample ids:
  rownames(design) = samples

  set.seed(61217)
  res = distinct_test(
    x = Kang_subset, 
    name_assays_expression = "logcounts",
    name_cluster = "cell",
    design = design,
    column_to_test = 2,
    min_non_zero_cells = 20,
    n_cores = 2)
  
  expect_is(res, "data.frame")
})
