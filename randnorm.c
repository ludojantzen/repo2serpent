/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : randnorm.c                                     */
/*                                                                           */
/* Created:       2019/05/03 (JLe)                                           */
/* Last modified: 2019/09/16 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Samples a random number from normal distribution.            */
/*                                                                           */
/* Comments: - Based on the rejection sampling method in Lux & Koblinger     */
/*             Sec. 2.I.F.                                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "RandNorm:"

/*****************************************************************************/

double RandNorm(double m, double v, long id)
{
  double x, y, mu, s;

  if (1 == 2)
    {
      /***********************************************************************/

      /***** Normal distribution *********************************************/

      /* Rejection loop */

      do
        {
          /* Sample two values */

          x = -log(RandF(id));
          y = -log(RandF(id));
        }
      while (0.5*(x - 1.0)*(x - 1.0) > y);

      /* Sample sign */

      if (RandF(id) < 0.5)
        x = -x;

      /* Scale */

      x = x*v + m;
      CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Log-normal distribution *****************************************/

      /* Rejection loop */

      do
        {
          /* Sample two values */

          x = -log(RandF(id));
          y = -log(RandF(id));
        }
      while (0.5*(x - 1.0)*(x - 1.0) > y);

      /* Sample sign */

      if (RandF(id) < 0.5)
        x = -x;

      /* Calculate mu and s */

      mu = log(m/sqrt(1 + v*v/m/m));
      s = sqrt(log(1 + v*v/m/m));

      /* Scale */

      x = exp(mu + s*x);
      CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);

      /***********************************************************************/
    }

  /* Return value */

  return x;
}

/*****************************************************************************/
