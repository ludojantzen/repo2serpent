/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : lagrangeinterp.c                               */
/*                                                                           */
/* Created:       2017/10/01 (TKa)                                           */
/* Last modified: 2018/01/30 (TKa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Lagrange interpolation                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/


#include "header.h"

double LagrangeInterpCubic(const double *x, const double *y, long nx,
                           double xx, long type) {
  /* Performs cubic (four-point) Lagrangian interpolation.
   *
   * Two points below and two above xx are used for interpolation. If xx is
   * in the first or last interval, first four or last four points are used,
   * respectively. If xx is not found in x, Die() is called.
   * */
  static char * const FUNCTION_NAME = "ProcessElectrons:";
  long idx, i, i0, i1, j;
  double yy, yyi;
  double xp[4], yp[4];

  /* Avoid compiler warning */
  i0 = -10;

  if (nx < 4)
    Die(FUNCTION_NAME, "Number of points %ld less than 4", nx);

  /* Find the value from the array */
  idx = SearchArray(x, xx, nx);

  /* Check the value and set starting index */
  if (idx < 0)
    Die(FUNCTION_NAME, "Value %E not found in the array", xx);
  else if (xx == x[0])
    return y[0];
  else if (xx == x[nx-1])
    return y[nx-1];
  else if (idx == 0)
    i0 = 0;
  else if (idx == nx - 2)
    i0 = idx - 2;
  else
    i0 = idx - 1;

  i1 = i0 + 4;

  /* Check values */
  CheckValue(FUNCTION_NAME, "i0", "", i0, 0, nx - 4);


  /* Check interpolation mode */
  switch (type) {
    case 1:
      {
        Die(FUNCTION_NAME, "Invalid interpolation mode %ld", type);
      }
    case 0:
    case 2:
      {
        /* Lin-lin */

        xp[0] = x[i0];
        xp[1] = x[i0 + 1];
        xp[2] = x[i0 + 2];
        xp[3] = x[i0 + 3];
        yp[0] = y[i0];
        yp[1] = y[i0 + 1];
        yp[2] = y[i0 + 2];
        yp[3] = y[i0 + 3];

        break;
      }
    case 3:
      {
        /* Lin-log */

        if ((x[i0] <= 0) || (x[i0 + 1] <= 0) || (x[i0 + 2] <= 0) ||
            (x[i0 + 3] <= 0))
          Die(FUNCTION_NAME, "Non-positive x value (lin-log)");

        xp[0] = log(x[i0]);
        xp[1] = log(x[i0 + 1]);
        xp[2] = log(x[i0 + 2]);
        xp[3] = log(x[i0 + 3]);
        yp[0] = y[i0];
        yp[1] = y[i0 + 1];
        yp[2] = y[i0 + 2];
        yp[3] = y[i0 + 3];

        xx = log(xx);

        break;
      }
    case 4:
      {
        /* Log-lin */

        if ((y[i0] <= 0) || (y[i0 + 1] <= 0) || (y[i0 + 2] <= 0) ||
            (y[i0 + 3] <= 0))
          Die(FUNCTION_NAME, "Non-positive y value (log-lin)");

        xp[0] = x[i0];
        xp[1] = x[i0 + 1];
        xp[2] = x[i0 + 2];
        xp[3] = x[i0 + 3];
        yp[0] = log(y[i0]);
        yp[1] = log(y[i0 + 1]);
        yp[2] = log(y[i0 + 2]);
        yp[3] = log(y[i0 + 3]);

        break;
      }
    case 5:
      {
        /* Log-log */

        if ((x[i0] <= 0) || (x[i0 + 1] <= 0) || (x[i0 + 2] <= 0) ||
            (x[i0 + 3] <= 0))
          Die(FUNCTION_NAME, "Non-positive x value (log-log)");
        if ((y[i0] <= 0) || (y[i0 + 1] <= 0) || (y[i0 + 2] <= 0) ||
            (y[i0 + 3] <= 0))
          Die(FUNCTION_NAME, "Non-positive y value (log-log)");

        xp[0] = log(x[i0]);
        xp[1] = log(x[i0 + 1]);
        xp[2] = log(x[i0 + 2]);
        xp[3] = log(x[i0 + 3]);
        yp[0] = log(y[i0]);
        yp[1] = log(y[i0 + 1]);
        yp[2] = log(y[i0 + 2]);
        yp[3] = log(y[i0 + 3]);

        xx = log(xx);

        break;
      }
    default:
      Die(FUNCTION_NAME, "Invalid interpolation mode %ld", type);
    }

  yy = 0.0;

  /* Loop over Lagrange basis polynomials */
  for (i = 0; i < 4; i++) {
    yyi = yp[i];

    /* Calculate the Lagrange polynomial */
    for (j = 0; j < 4; j++) {
      if (j == i)
        continue;
      yyi *= (xx - xp[j])/(xp[i] - xp[j]);

#ifdef DEBUG
      /* Check j */
      if ((j < 0) || (j > 4))
        Die(FUNCTION_NAME, "j out of limits [%ld, %ld]", 0, 4);
#endif
    }

#ifdef DEBUG
      /* Check i */
      if ((i < 0) || (i > 4))
        Die(FUNCTION_NAME, "i out of limits [%ld, %ld]", 0, 4);
#endif

    yy += yyi;
  }

  if ((type == 4) || (type == 5))
    yy = exp(yy);

  return yy;
}
