/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ttachain.c                                     */
/*                                                                           */
/* Created:       2011/06/10 (AIs)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Perusratkaisu TTA-ketjuille, joissa ei ole toistoja.         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "TTAChain:"

/*****************************************************************************/

double TTAChain(long n, double t, double *l, double *Pout)
{
  long i, j;
  double a, b, sum;
  double sumP;

  /* Lasketaan alphan osoittaja */
  
  b = 1.0;
  
  for (i = 0; i < n-1; i++)
    b = b*l[i];

  /* Lasketaan summaus */  

  sum  = 0.0;
  sumP = 0.0;

  for (i = 0; i < n; i++)
    {
      /* Lasketaan alphan nimittäjä (alpha=a*b) */
      
      a = 1.0;
      
      for (j = 0; j < n; j++)
        if (i != j)
          a = a/(l[j] - l[i]);

      sum  = sum  + a*exp(-l[i]*t);
      sumP = sumP + a*exp(-l[i]*t)/l[i];
    }

  /* Palautetaan arvot */
  
  if (l[n-1]==0)
    *Pout = 0.0;
  else
    *Pout = 1 - l[n-1]*b*sumP;

  return b*sum;
}

/*****************************************************************************/
