context("dodgr_dists_proportional")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

if (!test_all) {
  RcppParallel::setThreadOptions(numThreads = 2)
}

test_that("proportional dists", {

  expect_silent(graph <- weight_streetnet(hampi))

  nf <- 100
  nt <- 50
  set.seed (1)
  from <- sample(graph$from_id, size = nf)
  to <- sample(graph$to_id, size = nt)

  expect_error (d <- dodgr_dists_proportional (graph,
                                               from,
                                               to),
                "graph must have a column named 'edge_type'")

  graph$edge_type <- "A"
  expect_error (d <- dodgr_dists_proportional (graph,
                                               from,
                                               to),
                "'edge_type' column must contain integer values only")

  graph$edge_type <- -1L
  expect_error (d <- dodgr_dists_proportional (graph,
                                               from,
                                               to),
                "'edge_type' must contain non-negative values only")

  types <- sort (table (graph$highway), decreasing = TRUE)
  graph$edge_type <- 0L
  graph$edge_type [graph$highway != "path"] <- 1L

  expect_silent(d <- dodgr_dists_proportional(graph, from = from, to = to))
  expect_equal(nrow(d), nf)
  expect_equal(ncol(d), nt * 2) # twice as many rows as dodgr_dists!!
  expect_message(
    d2 <- dodgr_dists_proportional(graph, from = from, to = to, quiet = FALSE),
    "Calculating shortest paths ..."
  )
  expect_identical(d, d2)

  d0 <- d [, seq (nt)]
  d1 <- d [, nt + seq (nt)]
  index <- which (is.na (d1))
  d0 <- d0 [-index]
  d1 <- d1 [-index]
  expect_true (all (d1 <= d0))
  expect_true (mean (d0 - d1) > 0) # proportional dists are less
})
