/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : interpolatemupdf.c                             */
/*                                                                           */
/* Created:       2017/04/05 (VVa)                                           */
/* Last modified: 2017/04/05 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Interpolates the PDF value for a certain scattering cosine   */
/*              based on linear-linear tabular format                        */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InterpolateMuPDF:"

/*****************************************************************************/

double InterpolateMuPDF(long ptr, double mu)
{

  double d1, d2;
  double p1, p2;
  long np, k;
  double m, l, P;

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Initialize the d1 and d2 values */

  d1 = 0.0;
  d2 = 1.0;

  /* Get number of tabulated points */

  np = (long)RDB[ptr];

  /* Move pointer to the first tabulated cosine value */

  ptr++;

  /* Find mu-interval from the tabulated cosine values*/

  for (k = 0; k < (np - 1); k++)
    {
      /* Get limiting cosines */

      d1 = RDB[ptr + k];
      d2 = RDB[ptr + k + 1];

      /* Check if in the correct interval */

      if (mu < d2)
        break;
    }

  /* Check whether mu is smaller than the lower limit of the interval */

  if (mu < d1)
    Die(FUNCTION_NAME, "(mu < d1) np=%ld, k=%ld, mu=%e, d1=%e, d2=%e",
        np, k, mu, d1, d2);

  /* Limiting PDF values for this interval */

  p1 = RDB[ptr + np + k];
  p2 = RDB[ptr + np + k + 1];

  /* Calculate dx for interpolation */

  l = (mu - d1);

  /* Calculate slope for interpolation */

  m = (p2 - p1)/(d2 - d1);

  /* Interpolate PDF value linearly */

  P = p1 + m*l;

  return P;
}
