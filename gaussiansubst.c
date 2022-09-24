/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : gaussiansubst.c                                */
/*                                                                           */
/* Created:       2011/05/04 (MPu)                                           */
/* Last modified: 2014/08/22 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Gaussin eliminaatio takaisinsijoitusvaihe                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "GaussianSubst:"

/*****************************************************************************/

long GaussianSubst(struct ccsMatrix *lu, complex *b, complex *x, complex *diag)
{
  /* gaussin eliminaation sijoitusvaihe */

  long i,j,k,n;
  long *row, *col; 
  complex z, *val, *sum; 

  n       = lu->n;
  row     = lu->rowind; 
  col     = lu->colptr; 
  val     = lu->values;

  /* Alusta */
  
  for(i = 0 ;i < n; i++)
    {
      x[i].re = 0.0; 
      x[i].im = 0.0; 
    }
  
  sum = (complex *)Mem(MEM_ALLOC, n, sizeof(complex));

  /* kaikki on 0-pohjaista! */

  for (j=n-1 ; j >= 0 ; j--)
    {
      z = c_sub(b[j], sum[j]);

      if (c_norm(diag[j]) > 1.0e-20)
        x[j] = c_div(z, diag[j]);
      else
        {
          /* Print warning */

          Warn(FUNCTION_NAME, "matrix singular\n"); 

          /* Return error */

          return -1;
        }

      for (k = col[j]; k < col[j+1]; k++)
        { 
          i = row[k];
          if (i < 0)
            break;
      
          if (j > i)
            { 
              z = c_mul(val[k], x[j]);
              sum[i] = c_add(z,sum[i]);
            }
        }
    }  
  
  Mem(MEM_FREE, sum); 

  /* Return OK */

  return 0;
}

/*****************************************************************************/
