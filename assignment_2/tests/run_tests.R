# Run with Rscript tests/run_tests.R from the assignment root.
for (file in c("model.R", "theory.R", "simulation.R", "importance_sampling.R", "experiments.R", "stopped_survivors.R")) {
  source(file.path("project_A", file))
}
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)
