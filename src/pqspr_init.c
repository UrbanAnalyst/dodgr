#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _pqspr_rcpp_get_sp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
        {"_pqspr_rcpp_get_sp",       (DL_FUNC) &_pqspr_rcpp_get_sp,       2},
        {NULL, NULL, 0}
};


void R_init_pqspr(DllInfo *dll)
{
        R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
            R_useDynamicSymbols(dll, FALSE);
}
