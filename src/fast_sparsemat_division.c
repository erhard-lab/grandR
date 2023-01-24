// In C ----------------------------------------
#include <R.h>
#include <Rinternals.h>

SEXP fastsparsematdiv(SEXP Xi, SEXP Xj, SEXP Xx, SEXP Yi, SEXP Yj, SEXP Yx, SEXP dr) {
  int n = length(Xi);
  int *pXi, *pXj, *pYi, *pYj;
  double *pXx, *pYx, *pout;

  SEXP out = PROTECT(allocVector(REALSXP, n));

  pXi = INTEGER(Xi);
  pXj = INTEGER(Xj);
  pXx = REAL(Xx);
  pYi = INTEGER(Yi);
  pYj = INTEGER(Yj);
  pYx = REAL(Yx);
  pout = REAL(out);
  double pdr = REAL(dr)[0];

  int j = 0;
  for (int i = 0; i < n; i++) {
    for (; pXi[i]!=pYi[j] || pXj[i]!=pYj[j]; ++j);
    pout[i] = (pXx[i]/pdr)/pYx[i];
    if (pout[i]>1) pout[i]=1;
  }
  UNPROTECT(1);

  return out;
}
