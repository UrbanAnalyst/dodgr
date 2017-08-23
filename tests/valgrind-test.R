vg_check <- function ()
{
    vg <- system2 (command = 'R',
                   args = c ('-d "valgrind --tool=memcheck --leak-check=full"',
                             '-f valgrind-script.R'),
                   stdout = TRUE, stderr = TRUE)

    lost <- NULL
    types <- c ("definitely lost", "indirectly lost", "possibly lost")
    for (ty in types)
    {
        lost_type <- which (grepl (ty, vg))
        n <- regmatches(vg [lost_type], gregexpr("[[:digit:]]+", vg [lost_type]))
        lost <- c (lost, as.numeric (n [[1]] [2:3]))
    }
    if (any (lost > 0))
        stop ("valgrind memory leaks detected!")

    if (attr (vg, "status") != 0)
        stop ("valgrind error")
}

# TODO: Re-instate this test once valgrind = 3.12, not current 3.11; see
# https://www.mail-archive.com/kde-bugs-dist@kde.org/msg115932.html
if (identical (Sys.getenv ("TRAVIS"), "true"))
{
    #vg_check ()
}
