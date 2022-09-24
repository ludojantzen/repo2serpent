/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : vectorNorm.c                                   */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2016/09/30 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Compute 2-norm of complex-valued vector                      */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "vectorNorm:"

/*****************************************************************************/

double vectorNorm(long n, complex *v)
{
  long i; 
  double nrm;  

  /* Reset sum */

  nrm = 0.0; 
  
  /* Loop over values */

  for (i = 0; i < n; i++)
    {
      /* Check real and imaginary parts */
  
      CheckValue(FUNCTION_NAME, "v[i].re","", v[i].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "v[i].re","", v[i].im, -INFTY, INFTY);
  
      /* Add to sum */
      
      nrm = nrm + c_norm(v[i])*c_norm(v[i]); 
    }
  
  /* Square root */

  nrm = sqrt(nrm); 
  CheckValue(FUNCTION_NAME, "nrm","", nrm, 0.0, INFTY);

  /* Return value */

  return nrm; 
}

/*****************************************************************************/
