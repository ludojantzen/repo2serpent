/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : trapz.c                                        */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2014/06/07 (MPu)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Trapezoidal integration of complex-valued function           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "trapz:"

/*****************************************************************************/

complex trapz(long n, double *t, complex *F)
{
  long i; 
  double c; 
  complex I, z;

  /* Reset values */
  
  I.re = 0.0; 
  I.im = 0.0; 

  /* Loop over points */
  
  for (i = 0; i < n - 1; i++)
    {
      /* Calculate sum of values and width */

      z = c_add(F[i], F[i + 1]); 
      c = 0.5*(t[i + 1] - t[i]);  
      
      z.re = c*z.re; 
      z.im = c*z.im; 

      /* Check values */

      CheckValue(FUNCTION_NAME, "z.re","", z.re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "z,im","", z.im, -INFTY, INFTY);
      
      /* Add to result */

      I = c_add(I, z); 
    }    

  /* Return value of integral */

  return I;   
}

/*****************************************************************************/
