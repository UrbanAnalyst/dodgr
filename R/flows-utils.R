check_for_flow_col <- function (graph) {
    if ("flow" %in% names (graph)) {
        warning (
            "graph already has a 'flow' column; ",
            "this will be overwritten"
        )
    }
}

check_k <- function (k, from) {

    nk <- 1

    if (is.data.frame (k)) {
        k <- as.matrix (k)
    }

    if (is.matrix (k)) {
        if (nrow (k) != length (from)) {
            stop ("nrow(k) must equal length of 'from' points")
        }
        nk <- ncol (k)
    } else if (is.numeric (k)) {
        if (length (k) == 1) {
            k <- rep (k, length (from))
        } else if (length (k) != length (from)) {
            # convert to matrix
            nk <- length (k)
            k <- array (rep (k, each = length (from)),
                dim = c (length (from), nk)
            )
        }
    } else {
        stop ("'k' must be either a single value, a vector, or a matrix")
    }

    list (k = k, nk = nk)
}

check_pairwise_from_to <- function (from, to, flows) {

    nf <- ifelse (is.vector (from), length (from), nrow (from))
    nt <- ifelse (is.vector (to), length (to), nrow (to))

    if (nf != nt) {
        stop (
            "'from' and 'to' must be the same length or dimensions",
            "when using 'pairwise = TRUE'.",
            call. = FALSE
        )
    }
    if (nrow (flows) != nf) {
        stop (
            "'flows' must be the same length or dimensions as 'from' and 'to'",
            call. = FALSE
        )
    }
}

get_random_prefix <- function (prefix = "flow",
                               n = 5) {

    charvec <- c (letters, LETTERS, 0:9)
    rand_prefix <- paste (sample (charvec, n, replace = TRUE), collapse = "")
    fs::path (fs::path_temp (), paste0 (prefix, "_", rand_prefix))
}

nodes_arg_to_pts <- function (nodes,
                              graph) {

    if (!is.matrix (nodes) && is.numeric (nodes)) {
        nodes <- matrix (nodes, ncol = 2)
    }
    if (is.vector (nodes)) {
        # non-numeric, so must be vector of node IDs
        return (nodes)
    }
    if (ncol (nodes) == 2) {
        verts <- dodgr_vertices (graph)
        nodes <- verts$id [match_pts_to_verts (verts, nodes)]
    }
    return (nodes)
}


# keep from and to routing points in contracted graph
contract_graph_with_pts <- function (graph,
                                     from,
                                     to) {

    pts <- NULL
    if (!missing (from)) {
        pts <- c (pts, from)
    }
    if (!missing (to)) {
        pts <- c (pts, to)
    }
    dodgr_contract_graph (graph, unique (pts))
}
