#include <R.h>
#include <Rinternals.h>

SEXP sparse2dense(SEXP sparse_mat, SEXP fill_val) {
  if (!isS4(sparse_mat)) error("Input must be an S4 object");

  SEXP i = GET_SLOT(sparse_mat, install("i"));
  SEXP j = GET_SLOT(sparse_mat, install("j"));
  SEXP x = GET_SLOT(sparse_mat, install("x"));
  SEXP dims = GET_SLOT(sparse_mat, install("Dim"));
  
  int nrow = INTEGER(dims)[0];
  int ncol = INTEGER(dims)[1];
  int nnz = LENGTH(x);

  SEXP out = PROTECT(allocMatrix(REALSXP, nrow, ncol));
  double *out_ptr = REAL(out);
  double fill = REAL(fill_val)[0];

  for (int k = 0; k < nrow * ncol; k++) {
    out_ptr[k] = fill;
  }

  int *i_ptr = INTEGER(i);
  int *j_ptr = INTEGER(j);
  double *x_ptr = REAL(x);

  for (int k = 0; k < nnz; k++) {
    out_ptr[i_ptr[k] + j_ptr[k] * nrow] = x_ptr[k];
  }

  SEXP dimnames = GET_SLOT(sparse_mat, install("Dimnames"));
  if (!isNull(dimnames)) {
    setAttrib(out, R_DimNamesSymbol, dimnames);
  }

  UNPROTECT(1);
  return out;
}
