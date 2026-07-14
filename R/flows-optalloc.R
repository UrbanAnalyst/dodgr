#' @title Optimal allocation of flows from sources to capacity-limited targets
#'
#' @description Solve the transportation problem of allocating flows from `N`
#' source ('from') points with associated densities to `M` target ('to')
#' points with finite capacities, such that total flow from each source is
#' fully allocated, no target receives more than its capacity, and total
#' allocation cost (source-to-target network distance) is minimized. The
#' resultant optimal allocation is then aggregated on to the network with
#' \link{dodgr_flows_aggregate}, so this function returns a graph with an
#' additional `flow` column, just like \link{dodgr_flows_aggregate},
#' \link{dodgr_flows_disperse}, and \link{dodgr_flows_si}.
#'
#' @details This function performs an initial call to \link{dodgr_dists} to
#' obtain the `N x M` matrix of shortest-path distances between all `from` and
#' `to` points. The optimal allocation is then obtained by numerical
#' optimization over that matrix alone (see `control`, below), with no further
#' path-finding required, before finally calling \link{dodgr_flows_aggregate}
#' with the resultant allocation matrix as its `flows` argument.
#'
#' Because targets may collectively have more capacity than sources have
#' density, `sum(source_densities) <= sum(target_capacities)` must hold. This
#' is checked before any allocation is attempted.
#'
#' @inheritParams dodgr_flows_aggregate
#' @param source_densities Numeric vector of densities at each `from` point,
#' with `length(source_densities) == length(from)`.
#' @param target_capacities Numeric vector of capacities at each `to` point,
#' with `length(target_capacities) == length(to)`.
#' @param control A named list controlling the allocation algorithm. Must
#' include an `algorithm` entry, either `"sinkhorn"` (default) or `"lp"`:
#' \itemize{
#' \item `"sinkhorn"` solves an entropic-regularized approximation to the
#' optimal allocation via iterative matrix scaling. This is generally much
#' faster for large numbers of source/target points, at the cost of only
#' approximating the true optimum. Recognizes the following additional
#' `control` entries, all optional:
#' \itemize{
#' \item `lambda` Entropic regularization strength (default `1`). Smaller
#' values approach the exact optimum more closely, at the cost of slower and
#' less numerically stable convergence.
#' \item `tol` Convergence tolerance on row/column marginal errors (default
#' `1e-8`).
#' \item `maxiter` Maximum number of scaling iterations (default `1000`).
#' }
#' \item `"lp"` solves the exact transportation linear program via
#' \pkg{lpSolve}, which must be installed (it is only a "Suggested", not
#' "Imported" dependency). No additional `control` entries are used.
#' }
#'
#' @return The input `graph`, with an additional `flow` column added, similar to
#' behaviour of \link{dodgr_flows_aggregate}.
#'
#' @note The `"sinkhorn"` algorithm is generally the faster choice, especially
#' for large numbers of source/target points, but only approximates the true
#' optimal allocation, with accuracy controlled by `control$lambda`. Use
#' `control = list(algorithm = "lp")` for the exact optimum, at the cost of
#' both requiring the \pkg{lpSolve} package and scaling less well to large
#' numbers of points.
#'
#' @family flows
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' graphc <- dodgr_contract_graph (graph)
#' set.seed (1)
#' from <- sample (graphc$from_id, size = 10)
#' to <- sample (graphc$to_id, size = 5)
#' to <- to [!to %in% from]
#' source_densities <- runif (length (from))
#' target_capacities <- runif (length (to))
#' # scale target_capacities to ensure sum(source_densities) <=
#' # sum(target_capacities):
#' target_capacities <- target_capacities *
#'     1.5 * sum (source_densities) / sum (target_capacities)
#' graph <- dodgr_flows_optalloc (
#'     graph,
#'     from = from,
#'     to = to,
#'     source_densities = source_densities,
#'     target_capacities = target_capacities
#' )
#' # graph then has an additional 'flow' column, exactly as for
#' # 'dodgr_flows_aggregate'
#'
#' \dontrun{
#' # The exact optimum can be obtained instead with the 'lpSolve' package:
#' graph <- dodgr_flows_optalloc (
#'     graph,
#'     from = from,
#'     to = to,
#'     source_densities = source_densities,
#'     target_capacities = target_capacities,
#'     control = list (algorithm = "lp")
#' )
#' }
dodgr_flows_optalloc <- function (graph,
                                  from,
                                  to,
                                  source_densities,
                                  target_capacities,
                                  control = list (algorithm = "sinkhorn"),
                                  contract = TRUE,
                                  heap = "BHeap",
                                  norm_sums = TRUE,
                                  quiet = TRUE) {

    check_optalloc_inputs (from, to, source_densities, target_capacities)
    control <- check_optalloc_control (control)

    dist <- dodgr_dists (
        graph,
        from = from,
        to = to,
        heap = heap,
        quiet = quiet
    )
    dist <- as.matrix (dist)

    if (control$algorithm == "sinkhorn") {
        flows <- sinkhorn_optalloc (
            dist,
            source_densities,
            target_capacities,
            control
        )
    } else {
        flows <- lp_optalloc (dist, source_densities, target_capacities)
    }

    dodgr_flows_aggregate (
        graph,
        from = from,
        to = to,
        flows = flows,
        contract = contract,
        heap = heap,
        norm_sums = norm_sums,
        quiet = quiet
    )
}

#' Validate and default the `control` list of `dodgr_flows_optalloc`
#'
#' @param control The `control` argument as passed by the user
#' @return Validated version of `control`, with defaults filled in for any
#' missing `"sinkhorn"`-specific entries
#' @noRd
check_optalloc_control <- function (control) {

    if (!is.list (control)) {
        stop ("'control' must be a list")
    }
    if (!"algorithm" %in% names (control)) {
        stop ("'control' must include an 'algorithm' element")
    }

    algorithms <- c ("sinkhorn", "lp")
    if (!control$algorithm %in% algorithms) {
        stop (
            "control$algorithm must be one of '",
            paste (algorithms, collapse = "', '"), "'"
        )
    }

    defaults <- list (lambda = 1, tol = 1e-8, maxiter = 1000)

    unrecognized <- setdiff (names (control), c ("algorithm", names (defaults)))
    if (length (unrecognized) > 0) {
        warning (
            "Unrecognized 'control' names ignored: ",
            paste (unrecognized, collapse = ", ")
        )
    }

    if (control$algorithm == "sinkhorn") {
        for (n in names (defaults)) {
            if (n %in% names (control)) {
                val <- control [[n]]
                if (!is.numeric (val) || length (val) != 1 || val <= 0) {
                    stop ("control$", n, " must be a single numeric value > 0")
                }
            } else {
                control [[n]] <- defaults [[n]]
            }
        }
    }

    return (control)
}

#' Validate the non-`control` inputs of `dodgr_flows_optalloc`
#'
#' @inheritParams dodgr_flows_optalloc
#' @return Nothing; called for its side-effect of raising errors on invalid
#' input
#' @noRd
check_optalloc_inputs <- function (from,
                                   to,
                                   source_densities,
                                   target_capacities) {

    if (length (source_densities) != length (from)) {
        stop ("'source_densities' must have the same length as 'from'")
    }
    if (length (target_capacities) != length (to)) {
        stop ("'target_capacities' must have the same length as 'to'")
    }
    if (any (!is.finite (source_densities)) || any (source_densities < 0)) {
        stop ("'source_densities' must be finite and non-negative")
    }
    if (any (!is.finite (target_capacities)) || any (target_capacities < 0)) {
        stop ("'target_capacities' must be finite and non-negative")
    }

    sum_dens <- sum (source_densities)
    sum_cap <- sum (target_capacities)
    if (sum_dens > sum_cap) {
        stop (
            "sum(source_densities) = ", sum_dens,
            " exceeds sum(target_capacities) = ", sum_cap
        )
    }
}

#' Approximate optimal allocation via entropic-regularized optimal transport
#'
#' @param dist `N x M` matrix of source-target costs (distances)
#' @param source_densities Numeric vector of source densities, length `N`
#' @param target_capacities Numeric vector of target capacities, length `M`
#' @param control Validated `control` list, with `lambda`, `tol`, `maxiter`
#' all present
#' @return `N x M` allocation matrix
#' @noRd
sinkhorn_optalloc <- function (dist,
                               source_densities,
                               target_capacities,
                               control) {

    n <- nrow (dist)
    finite_dist <- dist [is.finite (dist)]
    penalty <- if (length (finite_dist) > 0) max (finite_dist) * 1e3 else 1
    dist [!is.finite (dist)] <- penalty

    slack <- sum (target_capacities) - sum (source_densities)
    dist_aug <- rbind (dist, rep (0, ncol (dist)))
    a <- c (source_densities, slack)
    b <- target_capacities

    dmax <- max (dist_aug)
    if (dmax <= 0) {
        dmax <- 1
    }
    kernel <- exp (-dist_aug / dmax / control$lambda)

    u <- rep (1, nrow (kernel))
    v <- rep (1, ncol (kernel))
    converged <- FALSE

    for (iter in seq_len (control$maxiter)) {

        u <- a / as.vector (kernel %*% v)
        v <- b / as.vector (t (kernel) %*% u)

        err_row <- max (abs (u * as.vector (kernel %*% v) - a))
        err_col <- max (abs (v * as.vector (t (kernel) %*% u) - b))

        if (max (err_row, err_col) < control$tol) {
            converged <- TRUE
            break
        }
    }

    if (!converged) {
        warning (
            "Sinkhorn algorithm did not converge within ",
            control$maxiter, " iterations"
        )
    }

    p_full <- kernel * outer (u, v)
    p_full [seq_len (n), , drop = FALSE]
}

#' Exact optimal allocation via the transportation linear program
#'
#' @inheritParams sinkhorn_optalloc
#' @return `N x M` allocation matrix
#' @noRd
lp_optalloc <- function (dist, source_densities, target_capacities) {

    if (!requireNamespace ("lpSolve", quietly = TRUE)) {
        stop (
            "control$algorithm = 'lp' requires the 'lpSolve' package; ",
            "please install it with install.packages('lpSolve')"
        )
    }

    finite_dist <- dist [is.finite (dist)]
    penalty <- if (length (finite_dist) > 0) max (finite_dist) * 1e3 else 1
    dist [!is.finite (dist)] <- penalty

    n <- nrow (dist)
    m <- ncol (dist)

    res <- lpSolve::lp.transport (
        cost.mat = dist,
        direction = "min",
        row.signs = rep ("=", n),
        row.rhs = source_densities,
        col.signs = rep ("<=", m),
        col.rhs = target_capacities,
        integers = NULL
    )

    if (res$status != 0) {
        stop ("lpSolve::lp.transport failed to find an optimal solution")
    }

    res$solution
}
