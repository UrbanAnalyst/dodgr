test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") ||
    identical (Sys.getenv ("GITHUB_JOB"), "test-coverage"))

testthat::skip_on_cran ()

test_that ("optalloc input validation", {

    from <- c ("a", "b", "c")
    to <- c ("x", "y")

    expect_error (
        check_optalloc_inputs (from, to, c (1, 2), c (1, 1)),
        "'source_densities' must have the same length as 'from'"
    )
    expect_error (
        check_optalloc_inputs (from, to, c (1, 1, 1), c (1)),
        "'target_capacities' must have the same length as 'to'"
    )
    expect_error (
        check_optalloc_inputs (from, to, c (1, -1, 1), c (1, 1)),
        "'source_densities' must be finite and non-negative"
    )
    expect_error (
        check_optalloc_inputs (from, to, c (1, 1, 1), c (1, NA)),
        "'target_capacities' must be finite and non-negative"
    )
    expect_error (
        check_optalloc_inputs (from, to, c (2, 2, 2), c (1, 1)),
        "sum\\(source_densities\\) = 6 exceeds sum\\(target_capacities\\) = 2"
    )
    expect_silent (
        check_optalloc_inputs (from, to, c (1, 1, 1), c (2, 2))
    )
})

test_that ("optalloc control validation", {

    expect_error (
        check_optalloc_control ("sinkhorn"),
        "'control' must be a list"
    )
    expect_error (
        check_optalloc_control (list (lambda = 1)),
        "'control' must include an 'algorithm' element"
    )
    expect_error (
        check_optalloc_control (list (algorithm = "nope")),
        "control\\$algorithm must be one of 'sinkhorn', 'lp'"
    )
    expect_error (
        check_optalloc_control (list (algorithm = "sinkhorn", lambda = -1)),
        "control\\$lambda must be a single numeric value > 0"
    )
    expect_error (
        check_optalloc_control (list (algorithm = "sinkhorn", tol = c (1, 2))),
        "control\\$tol must be a single numeric value > 0"
    )
    expect_error (
        check_optalloc_control (list (algorithm = "sinkhorn", maxiter = "a")),
        "control\\$maxiter must be a single numeric value > 0"
    )
    expect_warning (
        check_optalloc_control (list (algorithm = "sinkhorn", lamdba = 0.1)),
        "Unrecognized 'control' names ignored: lamdba"
    )

    control <- check_optalloc_control (list (algorithm = "sinkhorn"))
    expect_identical (control$lambda, 1)
    expect_identical (control$tol, 1e-8)
    expect_identical (control$maxiter, 1000)

    control <- check_optalloc_control (list (algorithm = "lp"))
    expect_identical (control, list (algorithm = "lp"))
})

test_that ("sinkhorn_optalloc marginals", {

    set.seed (2)
    dist <- matrix (runif (12, 1, 10), nrow = 3, ncol = 4)
    source_densities <- c (2, 3, 1)
    target_capacities <- c (2, 2, 2, 3) # sum 9 >= sum 6

    control <- check_optalloc_control (list (algorithm = "sinkhorn"))
    p <- sinkhorn_optalloc (dist, source_densities, target_capacities, control)

    expect_equal (dim (p), c (3, 4))
    expect_equal (rowSums (p), source_densities, tolerance = 1e-6)
    expect_true (all (colSums (p) <= target_capacities + 1e-6))
})

test_that ("lp_optalloc marginals", {

    testthat::skip_if_not_installed ("lpSolve")

    set.seed (2)
    dist <- matrix (runif (12, 1, 10), nrow = 3, ncol = 4)
    source_densities <- c (2, 3, 1)
    target_capacities <- c (2, 2, 2, 3)

    p <- lp_optalloc (dist, source_densities, target_capacities)

    expect_equal (dim (p), c (3, 4))
    expect_equal (rowSums (p), source_densities, tolerance = 1e-6)
    expect_true (all (colSums (p) <= target_capacities + 1e-6))
})

test_that ("sinkhorn approximates lp optimum", {

    testthat::skip_if_not_installed ("lpSolve")

    set.seed (3)
    dist <- matrix (runif (20, 1, 10), nrow = 4, ncol = 5)
    source_densities <- c (2, 1, 3, 1)
    target_capacities <- c (2, 2, 1, 1, 2)

    control <- check_optalloc_control (
        list (algorithm = "sinkhorn", lambda = 0.05, maxiter = 5000)
    )
    p_sinkhorn <- sinkhorn_optalloc (
        dist, source_densities, target_capacities, control
    )
    p_lp <- lp_optalloc (dist, source_densities, target_capacities)

    cost_sinkhorn <- sum (p_sinkhorn * dist)
    cost_lp <- sum (p_lp * dist)

    expect_true (cost_sinkhorn >= cost_lp)
    expect_lt ((cost_sinkhorn - cost_lp) / cost_lp, 0.2)
})

test_that ("dodgr_flows_optalloc end to end", {

    graph <- weight_streetnet (hampi)
    graphc <- dodgr_contract_graph (graph)
    set.seed (1)
    from <- sample (graphc$from_id, size = 10)
    to <- sample (graphc$to_id, size = 5)
    to <- to [!to %in% from]

    source_densities <- runif (length (from))
    target_capacities <- runif (length (to))
    target_capacities <- target_capacities *
        1.5 * sum (source_densities) / sum (target_capacities)

    expect_error (
        dodgr_flows_optalloc (
            graph,
            from = from,
            to = to,
            source_densities = source_densities * 1e6,
            target_capacities = target_capacities
        ),
        "exceeds sum\\(target_capacities\\)"
    )

    graph_sinkhorn <- dodgr_flows_optalloc (
        graph,
        from = from,
        to = to,
        source_densities = source_densities,
        target_capacities = target_capacities
    )
    expect_identical (ncol (graph_sinkhorn) - ncol (graph), 1L)
    expect_true ("flow" %in% names (graph_sinkhorn))
    if (test_all) { # fails on CRAN
        expect_gt (sum (graph_sinkhorn$flow, na.rm = TRUE), 0)
    }

    testthat::skip_if_not_installed ("lpSolve")

    graph_lp <- dodgr_flows_optalloc (
        graph,
        from = from,
        to = to,
        source_densities = source_densities,
        target_capacities = target_capacities,
        control = list (algorithm = "lp")
    )
    expect_identical (ncol (graph_lp) - ncol (graph), 1L)
    expect_true ("flow" %in% names (graph_lp))
    if (test_all) { # fails on CRAN
        expect_gt (sum (graph_lp$flow, na.rm = TRUE), 0)
    }
})
