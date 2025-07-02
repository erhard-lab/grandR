#include <R.h>
#include <Rinternals.h>
#include <string.h>

SEXP sparse2dense(SEXP sparse_mat, SEXP fill_val) {
  if (!isS4(sparse_mat)) error("Input must be an S4 object");
  
  // Get class name to determine matrix type
  SEXP class_name = R_do_slot(sparse_mat, install("class"));
  const char* class_str = CHAR(STRING_ELT(class_name, 0));
  
  // Get dimensions
  SEXP dims = R_do_slot(sparse_mat, install("Dim"));
  int nrow = INTEGER(dims)[0];
  int ncol = INTEGER(dims)[1];
  
  // Create output dense matrix filled with the fill value
  SEXP out = PROTECT(allocMatrix(REALSXP, nrow, ncol));
  double *out_ptr = REAL(out);
  double fill = REAL(fill_val)[0];
  
  // Initialize all elements with fill value
  for (int k = 0; k < nrow * ncol; k++) {
    out_ptr[k] = fill;
  }
  
  // Process based on matrix type
  if (strcmp(class_str, "dgCMatrix") == 0) {
    // Compressed Sparse Column format
    SEXP i = R_do_slot(sparse_mat, install("i"));       // Row indices
    SEXP p = R_do_slot(sparse_mat, install("p"));       // Column pointers
    SEXP x = R_do_slot(sparse_mat, install("x"));       // Non-zero values
    
    int *i_ptr = INTEGER(i);
    int *p_ptr = INTEGER(p);
    double *x_ptr = REAL(x);
    
    // Fill in the non-zero elements
    for (int j = 0; j < ncol; j++) {
      for (int k = p_ptr[j]; k < p_ptr[j+1]; k++) {
        int row = i_ptr[k];
        out_ptr[row + j * nrow] = x_ptr[k];
      }
    }
  } 
  else if (strcmp(class_str, "dgTMatrix") == 0) {
    // Triplet/COO format
    SEXP i = R_do_slot(sparse_mat, install("i"));       // Row indices
    SEXP j = R_do_slot(sparse_mat, install("j"));       // Column indices
    SEXP x = R_do_slot(sparse_mat, install("x"));       // Non-zero values
    int nnz = LENGTH(x);
    
    int *i_ptr = INTEGER(i);
    int *j_ptr = INTEGER(j);
    double *x_ptr = REAL(x);
    
    // Fill in the non-zero elements directly using row and column indices
    for (int k = 0; k < nnz; k++) {
      int row = i_ptr[k];
      int col = j_ptr[k];
      out_ptr[row + col * nrow] = x_ptr[k];
    }
  }
  else {
    UNPROTECT(1);
    error("Unsupported sparse matrix format: %s", class_str);
  }
  
  // Copy dimnames if they exist
  SEXP dimnames = R_do_slot(sparse_mat, install("Dimnames"));
  if (!isNull(dimnames)) {
    setAttrib(out, R_DimNamesSymbol, dimnames);
  }
  
  UNPROTECT(1);
  return out;
}
