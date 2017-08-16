#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _dodgr_rcpp_get_components(SEXP);
extern SEXP _dodgr_rcpp_get_sp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dodgr_rcpp_insert_vertices(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dodgr_rcpp_make_compact_graph(SEXP, SEXP, SEXP);
extern SEXP _dodgr_rcpp_points_index(SEXP, SEXP);
extern SEXP _dodgr_rcpp_sample_graph(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dodgr_rcpp_sf_as_network(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_dodgr_rcpp_get_components",     (DL_FUNC) &_dodgr_rcpp_get_components,     1},
    {"_dodgr_rcpp_get_sp",             (DL_FUNC) &_dodgr_rcpp_get_sp,             5},
    {"_dodgr_rcpp_insert_vertices",    (DL_FUNC) &_dodgr_rcpp_insert_vertices,    4},
    {"_dodgr_rcpp_make_compact_graph", (DL_FUNC) &_dodgr_rcpp_make_compact_graph, 3},
    {"_dodgr_rcpp_points_index",       (DL_FUNC) &_dodgr_rcpp_points_index,       2},
    {"_dodgr_rcpp_sample_graph",       (DL_FUNC) &_dodgr_rcpp_sample_graph,       4},
    {"_dodgr_rcpp_sf_as_network",      (DL_FUNC) &_dodgr_rcpp_sf_as_network,      2},
    {NULL, NULL, 0}
};

void R_init_dodgr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
