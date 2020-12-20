vg_check <- function () {

    vg <- system2 (command = 'R',
                   args = c ('-d "valgrind --tool=memcheck --leak-check=full"',
                             '-f valgrind-script.R'),
                   stdout = TRUE, stderr = TRUE)

    lost <- NULL
    types <- c ("definitely lost", "indirectly lost", "possibly lost")
    for (ty in types) {
        lost_type <- which (grepl (ty, vg))
        n <- regmatches (vg [lost_type],
                         gregexpr("[[:digit:]]+", vg [lost_type]))
        lost <- c (lost, as.numeric (n [[1]] [2:3]))
    }
    if (any (lost > 0))
        stop ("valgrind memory leaks detected!")

    if (attr (vg, "status") != 0)
        stop ("valgrind error")
}

# The valgrind test doesn't work with current (3.11) version because of this:
# https://www.mail-archive.com/kde-bugs-dist@kde.org/msg115932.html
# -> that's actually okay now, but valgrind detects memory leaks because of
# https://github.com/RcppCore/RcppParallel/issues/81
if (identical (Sys.getenv ("TRAVIS"), "true")) {
    #vv <- system2 (command = "valgrind", args = "--version", stdout = TRUE)
    #if (strsplit (vv, "valgrind-") [[1]] [2] >= "3.12.0")
    #    vg_check ()
}
