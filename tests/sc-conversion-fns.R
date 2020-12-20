genhash <- function (len = 10) {
    paste0 (sample (c (letters, LETTERS, 0:9), size = len), collapse = "")
}

sf_to_sc <- function (x) {

    pts <- do.call (rbind, x$geometry)
    pts <- data.frame ("x_" = pts [, 1],
                       "y_" = pts [, 2],
                       "vertex_" = rownames (pts),
                       stringsAsFactors = FALSE)

    edge <- lapply (x$geometry, function (i)
                     cbind (rownames (i) [1:(nrow (i) - 1)],
                            rownames (i) [2:nrow (i)]))
    for (e in seq (edge))
        edge [[e]] <- cbind (edge [[e]], names (x$geometry) [e])
    edge <- data.frame (do.call (rbind, edge), stringsAsFactors = FALSE)
    edge$edge_ <- vapply (seq (nrow (edge)), function (i) genhash (10),
                           character (1))

    object_link_edge <- data.frame (edge_ = edge$edge_,
                                    object_ = edge$X3,
                                    native_ = TRUE,
                                    stringsAsFactors = FALSE)
    edge <- data.frame (".vx0" = edge$X1,
                        ".vx1" = edge$X2,
                        "edge_" = edge$edge_,
                        stringsAsFactors = FALSE)

    x_no_g <- x
    x_no_g$geometry <- NULL
    osm_id <- as.character (x_no_g$osm_id)
    x_no_g$osm_id <- NULL
    for (i in names (x_no_g))
        x_no_g [[i]] <- as.character (x_no_g [[i]])
    x_no_g <- as.list (x_no_g)
    x_no_g <- lapply (x_no_g, function (i) {
                          res <- cbind (osm_id, i)
                          res [which (!is.na (res [, 2])), , drop = FALSE]
                                    })
    for (i in seq (x_no_g))
        x_no_g [[i]] <- cbind (x_no_g [[i]], names (x_no_g) [i])
    x_no_g <- data.frame (do.call (rbind, x_no_g),
                          stringsAsFactors = FALSE)
    object <- data.frame ("object_" = x_no_g$osm_id,
                          key = x_no_g$V3,
                          value = x_no_g$i,
                          stringsAsFactors = FALSE)
    object <- object [order (object$object_), ]

    res <- list (nodes = NULL,
                 object = object,
                 object_link_edge = object_link_edge,
                 edge = edge,
                 vertex = pts)
    class (res) <- c ("SC", "sc", "osmdata_sc")
    return (res)
}
