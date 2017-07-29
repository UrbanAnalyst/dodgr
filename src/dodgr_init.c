#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _dodgr_rcpp_get_sp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dodgr_rcpp_lines_as_network(SEXP, SEXP);
extern SEXP _dodgr_rcpp_make_compact_graph(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_dodgr_rcpp_get_sp",             (DL_FUNC) &_dodgr_rcpp_get_sp,             5},
    {"_dodgr_rcpp_lines_as_network",   (DL_FUNC) &_dodgr_rcpp_lines_as_network,   2},
    {"_dodgr_rcpp_make_compact_graph", (DL_FUNC) &_dodgr_rcpp_make_compact_graph, 2},
    {NULL, NULL, 0}
};


void R_init_dodgr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
