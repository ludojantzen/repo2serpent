/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ccsmatrixprint.c                               */
/*                                                                           */
/* Created:       2011/05/02 (MPu)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Printtaa matriisin                                           */
/*                                                                           */
/* Comments: - Käyttö debuggaukseen                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ccsMatrixPrint:"

/*****************************************************************************/

void ccsMatrixPrint(struct ccsMatrix *cmat)
{
  long k, j;

  fprintf(outp, "ccs matrix\n\n");

  for (k = 0; k < cmat->m; k++)
    {
      fprintf(outp, "column = %ld\n", k + 1); 
    
      for (j = cmat->colptr[k]; j < cmat->colptr[k + 1]; j++)
        fprintf(outp, "(%+-12.5e, %+-12.5e), %ld\n", cmat->values[j].re, 
               cmat->values[j].im, cmat->rowind[j] + 1);
    }
}

/*****************************************************************************/
