test_that("test_DTU() works faultlessly.", {
  set.seed(61217)
  
  main = HDCytoData::Weber_BCR_XL_sim_main_SE()
  main
  
  # normalize data:
  main <- diffcyt::transformData(main, cofactor = 5)
  
  # select all State markers:
  # sel_cols = colData(main)$marker_class == "state"
  # select 1 marker only:
  sel_cols = colData(main)$marker_name == "pNFkB"
  
  # create a SummarizedExperiment object with seleceted columns, and inverting row-column structure:
  sce <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = t(assays(main)$exprs[, sel_cols]) ), 
    colData = rowData(main), 
    rowData = colnames(main)[sel_cols], 
    metadata = list(experiment_info =  metadata(main)$experiment_info ))
  
  # Perform differential analyses, within each cluster of cells, between conditions:
  res = distinct_test(
    sce, 
    name_assays_expression = "exprs",
    name_cluster = "population_id",
    name_group = "group_id",
    name_sample = "sample_id",
    P = 10^2, 
    min_non_zero_cells_per_group = 0)
  
  expect_is(res, "data.frame")
})
