#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _GroupAssignment_optimalAssignment(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_GroupAssignment_optimalAssignment", (DL_FUNC) &_GroupAssignment_optimalAssignment, 4},
  {NULL, NULL, 0}
};

void R_init_GroupAssignment(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}