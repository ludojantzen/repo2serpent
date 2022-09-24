/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : matProduct.c                                   */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2014/06/07 (MPu)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Product of two complex matrices                              */
/*                                                                           */
/* Comments: A = n x m, B = m x p => AB = n x p                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "matProduct:"

/*****************************************************************************/

void matProduct(long n, long m, long p, complex **A, complex **B, complex **AB)
{
  long i, j, k; 
  complex z; 

  /* Clear product matrix */

  for (i = 0; i < n; i++)
    memset(AB[i], 0, p*sizeof(complex)); 

  /* Loop over values */

  for(i = 0; i < n; i++)
    for (j = 0; j < p; j++)
      for (k = 0; k < m; k++)
        {
          /* Check values */
          
          CheckValue(FUNCTION_NAME, "A[i][k].re","", A[i][k].re,-INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "A[i][k].im","", A[i][k].im,-INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "B[k][j].re","", B[k][j].re,-INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "B[k][j].im","", B[k][j].im,-INFTY, INFTY);

          /* Calculate product */
          
          z = c_mul(A[i][k], B[k][j]); 

          /* Add to result */

          AB[i][j] = c_add(AB[i][j], z);  
        }      
}

/*****************************************************************************/



