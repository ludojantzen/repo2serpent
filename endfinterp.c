/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : endfinterp.c                                   */
/*                                                                           */
/* Created:       2011/08/13 (JLe)                                           */
/* Last modified: 2018/11/14 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Interpolates between tabulated values according to ENDF      */
/*              interpolation schemes.                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ENDFInterp:"

/*****************************************************************************/

double ENDFInterp(long type, double x, double x0, double x1, double y0, 
                  double y1)
{
  double y;

  /* Check input values */

  CheckValue(FUNCTION_NAME, "x", "", x, 0.0, INFTY);
  CheckValue(FUNCTION_NAME, "x0", "", x0, 0.0, INFTY);
  CheckValue(FUNCTION_NAME, "x1", "", x1, 0.0, INFTY);
  CheckValue(FUNCTION_NAME, "y0", "", y0, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y1", "", y1, -INFTY, INFTY);

  /* Avoid compiler warning */

  y = 0.0;

  /* Check type */
  
  switch (type)
    {
    case 1:
      {
        /* Histogram */

        y = y0;

        /* Break loop */

        break;
      }
    case 0:
    case 2:
      {
        /* Lin-lin */

        if (x1 != x0)
          y = (y1 - y0)*(x - x0)/(x1 - x0) + y0;
        else if ((y1 == y0) || (x == x0))
          y = y0;
        else
          Die(FUNCTION_NAME, "Division by zero (type %ld)", type);

        /* Break loop */

        break;
      }
    case 3:
      {
        /* Lin-log */

        y = (y1 - y0)*log(x/x0)/log(x1/x0) + y0;
        
        /* Break loop */

        break;
      }
    case 4:
      {
        /* Log-lin */

        if ((y0 != 0.0) && (x1 != x0))
          y = pow(y1/y0, (x - x0)/(x1 - x0))*y0;
        else
          Die(FUNCTION_NAME, "Division by zero (type %ld)", type);

        /* Break loop */

        break;
      }
    case 5:
      {
        /* Log-log */

        if (y0 != 0.0)
          y = pow(y1/y0, log(x/x0)/log(x1/x0))*y0;
        else
          Die(FUNCTION_NAME, "Division by zero (type %ld)", type);

        /* Break loop */

        break;
      }
    default:
      Die(FUNCTION_NAME, "Invalid interpolation mode %ld");
    }

  /* Check output value */

  CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);

  /* Return value */

  return y;
}

/*****************************************************************************/
