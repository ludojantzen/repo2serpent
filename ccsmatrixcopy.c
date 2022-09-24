/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ccsmatrixcopy.c                                */
/*                                                                           */
/* Created:       2011/05/02 (MPu)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Kopioi matriisin sisältöineen                                */
/*                                                                           */
/* Comments: - Matriiseilla pitää olla samat speksit                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ccsMatrixCopy:"

/*****************************************************************************/

void ccsMatrixCopy(struct ccsMatrix *cmat, struct ccsMatrix *cmat_copy)
{
  long n, nnz; 

  if( (cmat->n != cmat_copy->n) ||  (cmat->m != cmat_copy->m) ||
      (cmat->nnz != cmat_copy->nnz))
    Die(FUNCTION_NAME, "Matrix parameters not consistent\n"); 

  n   = cmat->n; 
  nnz = cmat->nnz; 

  memcpy(cmat_copy->colptr, cmat->colptr, (n+1)*sizeof(long));
  memcpy(cmat_copy->rowind, cmat->rowind,  nnz *sizeof(long));
  memcpy(cmat_copy->values, cmat->values,  nnz *sizeof(complex));

}

/*****************************************************************************/
