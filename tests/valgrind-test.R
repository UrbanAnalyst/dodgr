vg_check <- function ()
{
    vg <- system2 (command = 'R',
                   args = c ('-d "valgrind --tool=memcheck --leak-check=full"',
                             '-f test.R'),
                   stdout = TRUE, stderr = TRUE)
    if (attr (vg, "status") != 0)
        stop ("valgrind error")
}
if (identical (Sys.getenv ("TRAVIS"), "true"))
{
    vg_check ()
}
