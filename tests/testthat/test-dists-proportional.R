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
  graph$edge_type <- match (graph$highway, names (types))
  graph$edge_type [graph$edge_type == 2] <- 22L
  expect_error (d <- dodgr_dists_proportional (graph,
                                               from,
                                               to),
                "'edge_type' values must be sequential integers")

  graph$edge_type <- 0L
  graph$edge_type [graph$highway != "path"] <- 1L

  expect_silent(d <- dodgr_dists_proportional(graph, from = from, to = to))
  expect_type (d, "list")
  expect_length (d, 2L) # distance + 1L; 0 is ignored

  # All dimensions should be equal:
  dims <- vapply (d, dim, integer (2))
  expect_identical (sum (diff (dims [1, ])), 0L)
  expect_identical (sum (diff (dims [2, ])), 0L)

  expect_message(
    d2 <- dodgr_dists_proportional(graph, from = from, to = to, quiet = FALSE),
    "Calculating shortest paths ..."
  )
  expect_identical(d, d2)
})
