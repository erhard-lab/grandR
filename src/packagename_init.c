#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP fastsparsematcompmult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fastsparsematcompmult1m(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fastsparsematdiv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"fastsparsematcompmult",   (DL_FUNC) &fastsparsematcompmult,   6},
  {"fastsparsematcompmult1m", (DL_FUNC) &fastsparsematcompmult1m, 6},
  {"fastsparsematdiv",        (DL_FUNC) &fastsparsematdiv,        7},
  {NULL, NULL, 0}
};

void R_init_grandR(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
