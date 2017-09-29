#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _dodgr_rcpp_contract_graph(SEXP, SEXP);
extern SEXP _dodgr_rcpp_get_component_vector(SEXP);
extern SEXP _dodgr_rcpp_get_sp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dodgr_rcpp_points_index(SEXP, SEXP);
extern SEXP _dodgr_rcpp_sample_graph(SEXP, SEXP);
extern SEXP _dodgr_rcpp_sf_as_network(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_dodgr_rcpp_contract_graph",       (DL_FUNC) &_dodgr_rcpp_contract_graph,       2},
    {"_dodgr_rcpp_get_component_vector", (DL_FUNC) &_dodgr_rcpp_get_component_vector, 1},
    {"_dodgr_rcpp_get_sp",               (DL_FUNC) &_dodgr_rcpp_get_sp,               5},
    {"_dodgr_rcpp_points_index",         (DL_FUNC) &_dodgr_rcpp_points_index,         2},
    {"_dodgr_rcpp_sample_graph",         (DL_FUNC) &_dodgr_rcpp_sample_graph,         2},
    {"_dodgr_rcpp_sf_as_network",        (DL_FUNC) &_dodgr_rcpp_sf_as_network,        2},
    {NULL, NULL, 0}
};

void R_init_dodgr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
