/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sampletabular.c                                */
/*                                                                           */
/* Created:       2017/03/06 (JLe)                                           */
/* Last modified: 2020/06/24 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Samples random variable from a 1D tabular distribution.      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleTabular:"

/*****************************************************************************/

double SampleTabular(const double *x, const double *pdf, const double *cdf,
                     long ne, long i, long id)
{
  long n, k;
  double val, x0, x1, p0, p1, c0, c1, rnd, a, b;

  /* Avoid compiler warning */

  val = 0.0;

  /* Re-sampling loop */

  for (k = 0; k < 10; k++)
    {
      /* Sample random variable */

      rnd = RandF(id);

      /* search interval */

      if ((n = SearchArray(cdf, rnd, ne)) < 0)
        Die(FUNCTION_NAME, "Search failed");
      else if (n > ne - 2)
        Die(FUNCTION_NAME, "WTF?");

      /* Get values */

      x0 = x[n];
      x1 = x[n + 1];

      c0 = cdf[n];
      c1 = cdf[n + 1];

      p0 = pdf[n];
      p1 = pdf[n + 1];

      /* Check */

      CheckValue(FUNCTION_NAME, "rnd", "", rnd, c0, c1);

      /* Check distribution type */

      if (i == 0)
        {
          /* Line spectrum */

          val = x1;

          /* Break loop */

          break;
        }
      else if (i == 1)
        {
          /* Histogram */

          if (p0 > 0.0)
            val = x0 + (rnd - c0)/p0;
          else
            {
              /* Resample */

              continue;
            }

          /* Break loop */

          break;
        }
      else if (i == 2)
        {
          /* Linear-linear */

          a = (p1 - p0)/(x1 - x0);
          b = p0*p0 + 2.0*a*(rnd - c0);

          /* Check */

          if ((a > 0.0) && (b >= 0.0))
            {
              /* Interpolate */

              val = x0 + (sqrt(b) - p0)/a;
            }
          else
            {
              /* Resample */

              continue;
            }

          /* Break loop */

          break;
        }
      else if (i == 4)
        {
          /* Log-linear */

          a = log(p1/p0);
          b = x1 - x0;

          /* Check values */

          CheckValue(FUNCTION_NAME, "a", "", a, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "b", "", b, ZERO, INFTY);

          /* Interpolate */

          val = x0 + (log((a*(rnd - c0))/(p0*b) + 1.0)*b)/a;

          /* Break loop */

          break;
        }
      else
        Die(FUNCTION_NAME, "Unsupported interpolation type %ld", i);
    }

  /* Check number of resamples */

  if (k == 10)
    {
      /* Failed */

      Warn(FUNCTION_NAME, "Sampling failed");

      /* Return error */

      return -1.0;
    }
  else if (k > 0)
    Warn(FUNCTION_NAME, "Re-sampled %ld times", k);

  /* Check value */

  CheckValue(FUNCTION_NAME, "val", "", val, ZERO, INFTY);

  /* Return value */

  return val;
}

/*****************************************************************************/
