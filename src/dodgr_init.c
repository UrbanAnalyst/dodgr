#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _dodgr_rcpp_aggregate_flows(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dodgr_rcpp_contract_graph(SEXP, SEXP);
extern SEXP _dodgr_rcpp_get_component_vector(SEXP);
extern SEXP _dodgr_rcpp_get_paths(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dodgr_rcpp_get_sp_dists(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dodgr_rcpp_merge_flows(SEXP);
extern SEXP _dodgr_rcpp_points_index(SEXP, SEXP);
extern SEXP _dodgr_rcpp_sample_graph(SEXP, SEXP);
extern SEXP _dodgr_rcpp_sf_as_network(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_dodgr_rcpp_aggregate_flows",      (DL_FUNC) &_dodgr_rcpp_aggregate_flows,      6},
    {"_dodgr_rcpp_contract_graph",       (DL_FUNC) &_dodgr_rcpp_contract_graph,       2},
    {"_dodgr_rcpp_get_component_vector", (DL_FUNC) &_dodgr_rcpp_get_component_vector, 1},
    {"_dodgr_rcpp_get_paths",            (DL_FUNC) &_dodgr_rcpp_get_paths,            5},
    {"_dodgr_rcpp_get_sp_dists",         (DL_FUNC) &_dodgr_rcpp_get_sp_dists,         5},
    {"_dodgr_rcpp_merge_flows",          (DL_FUNC) &_dodgr_rcpp_merge_flows,          1},
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
