/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : truncate.c                                     */
/*                                                                           */
/* Created:       2012/01/05 (JLe)                                           */
/* Last modified: 2012/01/05 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Truncates a double to requested precision                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "Truncate:"

/*****************************************************************************/

double Truncate(double x0, long prec)
{
  long n;
  double x1, nex, f;

  /* Check zero */

  if (x0 == 0.0)
    return x0;

  /* Reset variables */

  nex = 0.0;
  x1 = x0;

  /* Check range */

  if (fabs(x0) < 1.0)
    {
      /* Adjust */

      do 
	{
	  x1 = x1*10.0;
	  nex = nex - 1.0;
	}
      while (fabs(x1) < 1.0);
    }
  else if (fabs(x0) > 10.0)
    {
      /* Adjust */

      do 
	{
	  x1 = x1/10.0;
	  nex = nex + 1.0;
	}
      while (fabs(x1) >= 10.0);
    }

  /* Check range */

  CheckValue(FUNCTION_NAME, "x1", "", fabs(x1), 1.0, 10.0);

  /* Calculate multiplier */

  f = pow(10.0, (double)prec);

  /* Truncate decimal part */

  n = round(x1*f);
  x1 = ((double)n)/f;

  /* Put exponent */

  x1 = x1*pow(10.0, nex);

  /* Return truncated value */

  return x1;
}

/*****************************************************************************/
