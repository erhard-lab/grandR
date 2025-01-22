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
    pout[i] = (pXx[i]/pdr)/pYx[j];
    if (pout[i]>1) pout[i]=1;
  }
  UNPROTECT(1);

  return out;
}

SEXP fastsparsematsum(SEXP Xi, SEXP Xj, SEXP Xx, SEXP Yi, SEXP Yj, SEXP Yx) {
  int n = length(Xi);
  int m = length(Yi);
  int cmp, k, x, y;
  int *pXi, *pXj, *pYi, *pYj, *pRi, *pRj;
  double *pXx, *pYx, *pRx;

  int len = n<m?n:m;
  SEXP Ri = PROTECT(allocVector(INTSXP, len));
  SEXP Rj = PROTECT(allocVector(INTSXP, len));
  SEXP Rx = PROTECT(allocVector(REALSXP, len));

  pXi = INTEGER(Xi);
  pXj = INTEGER(Xj);
  pXx = REAL(Xx);
  pYi = INTEGER(Yi);
  pYj = INTEGER(Yj);
  pYx = REAL(Yx);
  pRi = INTEGER(Ri);
  pRj = INTEGER(Rj);
  pRx = REAL(Rx);

  x = 0;
  y = 0;
  k = 0;
  while (x<n && y<m) {
    cmp = pXj[x]-pYj[y];
    if (cmp==0) cmp = pXi[x]-pYi[y];
    if (cmp==0) {
      pRi[k] = pXi[x];
      pRj[k] = pXj[x];
      pRx[k] = pXx[x]+pYx[y];
      ++k; ++x; ++y;
    }
    else if (cmp<0) {
      ++x;
    }
    else {
      ++y;
    }
  }
//  SETLENGTH(Ri, k);
//  SETLENGTH(Rj, k);
//  SETLENGTH(Rx, k);
  Ri = Rf_xlengthgets(Ri, k);
  Rj = Rf_xlengthgets(Rj, k);
  Rx = Rf_xlengthgets(Rx, k);

  SEXP re = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(re, 0, Ri);
  SET_VECTOR_ELT(re, 1, Rj);
  SET_VECTOR_ELT(re, 2, Rx);
  UNPROTECT(4);

  return re;
}


// X is count (or norm or similar), Y is ntr
SEXP fastsparsematcompmult(SEXP Xi, SEXP Xj, SEXP Xx, SEXP Yi, SEXP Yj, SEXP Yx) {
  int n = length(Xi);
  int m = length(Yi);
  int cmp, k, x, y;
  int *pXi, *pXj, *pYi, *pYj, *pRi, *pRj;
  double *pXx, *pYx, *pRx;

  int len = n<m?n:m;
  SEXP Ri = PROTECT(allocVector(INTSXP, len));
  SEXP Rj = PROTECT(allocVector(INTSXP, len));
  SEXP Rx = PROTECT(allocVector(REALSXP, len));

  pXi = INTEGER(Xi);
  pXj = INTEGER(Xj);
  pXx = REAL(Xx);
  pYi = INTEGER(Yi);
  pYj = INTEGER(Yj);
  pYx = REAL(Yx);
  pRi = INTEGER(Ri);
  pRj = INTEGER(Rj);
  pRx = REAL(Rx);

  x = 0;
  y = 0;
  k = 0;
  while (x<n && y<m) {
    cmp = pXj[x]-pYj[y];
    if (cmp==0) cmp = pXi[x]-pYi[y];
    if (cmp==0) {
      pRi[k] = pXi[x];
      pRj[k] = pXj[x];
      pRx[k] = pXx[x]*pYx[y];
      ++k; ++x; ++y;
    }
    else if (cmp<0) {
      ++x;
    }
    else {
      ++y;
    }
  }
//SETLENGTH(Ri, k);
//SETLENGTH(Rj, k);
//SETLENGTH(Rx, k);
  Ri = Rf_xlengthgets(Ri, k);
  Rj = Rf_xlengthgets(Rj, k);
  Rx = Rf_xlengthgets(Rx, k);

  SEXP re = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(re, 0, Ri);
  SET_VECTOR_ELT(re, 1, Rj);
  SET_VECTOR_ELT(re, 2, Rx);
  UNPROTECT(4);

  return re;
}



// X is count (or norm or similar), Y is ntr, compute count*(1-ntr)
SEXP fastsparsematcompmult1m(SEXP Xi, SEXP Xj, SEXP Xx, SEXP Yi, SEXP Yj, SEXP Yx) {
  int n = length(Xi);
  int m = length(Yi);
  int cmp, k, x, y;
  int *pXi, *pXj, *pYi, *pYj, *pRi, *pRj;
  double *pXx, *pYx, *pRx;

  int len = n;
  SEXP Ri = PROTECT(allocVector(INTSXP, len));
  SEXP Rj = PROTECT(allocVector(INTSXP, len));
  SEXP Rx = PROTECT(allocVector(REALSXP, len));

  pXi = INTEGER(Xi);
  pXj = INTEGER(Xj);
  pXx = REAL(Xx);
  pYi = INTEGER(Yi);
  pYj = INTEGER(Yj);
  pYx = REAL(Yx);
  pRi = INTEGER(Ri);
  pRj = INTEGER(Rj);
  pRx = REAL(Rx);

  x = 0;
  y = 0;
  k = 0;
  while (x<n && y<m) {
    cmp = pXj[x]-pYj[y];
    if (cmp==0) cmp = pXi[x]-pYi[y];
    if (cmp==0) {
      if (pYx[y]!=1) {
        pRi[k] = pXi[x];
        pRj[k] = pXj[x];
        pRx[k] = pXx[x]*(1-pYx[y]);
        ++k;
      }
      ++x; ++y;
    }
    else if (cmp<0) {
      pRi[k] = pXi[x];
      pRj[k] = pXj[x];
      pRx[k] = pXx[x];
      ++k;
      ++x;
    }
    else {
      ++y;
    }
  }
//  SETLENGTH(Ri, k);
//  SETLENGTH(Rj, k);
//  SETLENGTH(Rx, k);
  Ri = Rf_xlengthgets(Ri, k);
  Rj = Rf_xlengthgets(Rj, k);
  Rx = Rf_xlengthgets(Rx, k);

  SEXP re = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(re, 0, Ri);
  SET_VECTOR_ELT(re, 1, Rj);
  SET_VECTOR_ELT(re, 2, Rx);
  UNPROTECT(4);

  return re;
}

