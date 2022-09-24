/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : feifc.c                                        */
/*                                                                           */
/* Created:       2018/01/29 (BWe)                                           */
/* Last modified: 2018/02/23 (BWe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: FE-defined material factors for multi-physics interface      */
/*                                                                           */
/* Comments: - This subroutine is intended to be used for user-defined       */
/*             density and temperature distributions with the multi-physics  */
/*             interface (types 31 & 32).                                    */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ExpandFE:"

/*****************************************************************************/

double ExpandFE(const double *const params, const double *const coefs,
                double x, double y, double z, long threadID)
{
  long a, ncoef, deg0, deg1, deg2, idx0, idx1, idx2, vector0, vector1, vector2;
  double axial, center0, center1, min0, max0, min1, max1, min2, max2,
       norm0, norm1, norm2, pos0, pos1, radius, theta, part0, part1, value;

  value = 0.0;

  switch ((long)params[FET_PARAM_TYPE])
  {
    case FET_TYPE_CARTESIAN:
    {
      /* Get X-parameters and check location */

      min0 = params[FET_CART_MIN_X];
      max0 = params[FET_CART_MAX_X];
      if ((x < min0) || (max0 < x))
        return INFTY;
      norm0 = 2.0 * (x - min0) / (max0 - min0) - 1;
      deg0 = params[FET_CART_ORDER_X];

      /* Get Y-parameters and check location */

      min1 = params[FET_CART_MIN_Y];
      max1 = params[FET_CART_MAX_Y];
      if ((y < min1) || (max1 < y))
        return INFTY;
      norm1 = 2.0 * (y - min1) / (max1 - min1) - 1;
      deg1 = params[FET_CART_ORDER_Y];

      /* Get Z-parameters and check location */

      min2 = params[FET_CART_MIN_Z];
      max2 = params[FET_CART_MAX_Z];
      if ((z < min2) || (max2 < z))
        return INFTY;
      norm2 = 2.0 * (z - min2) / (max2 - min2) - 1;
      deg2 = params[FET_CART_ORDER_Z];

      /* Get the locations of the coefficient arrays */

      vector0 = (long)params[FET_CART_CALC_ARRAY_PTR_X];
      vector1 = (long)params[FET_CART_CALC_ARRAY_PTR_Y];
      vector2 = (long)params[FET_CART_CALC_ARRAY_PTR_Z];

      /* Calculate the functional expansions */

      PolynomialLegendre(deg0, norm0, vector0, threadID, FET_CALCULATION_STANDARD);
      PolynomialLegendre(deg1, norm1, vector1, threadID, FET_CALCULATION_STANDARD);
      PolynomialLegendre(deg2, norm2, vector2, threadID, FET_CALCULATION_STANDARD);

      /* Loop over all the coefficients to expand the series at this location */

      for (idx0 = 0, a = 0; idx0 <= deg0; idx0++)
      {
        part0 = GetPrivateData(vector0 + idx0, threadID);

        for (idx1 = 0; idx1 <= deg1; idx1++)
        {
          part1 = part0 * GetPrivateData(vector1 + idx1, threadID);

          for (idx2 = 0; idx2 <= deg2; idx2++, a++)\
            value += coefs[a] * part1 * GetPrivateData(vector2 + idx2, threadID);
        }
      }

      break;
    }

    case FET_TYPE_CYLINDRICAL:
    {
      /* Get the coordinates of the axial center */

      center0 = params[FET_CYL_CENTER_A0];
      center1 = params[FET_CYL_CENTER_A1];
      pos0 = 0;
      pos1 = 0;
      axial = 0;

      switch ((long)params[FET_CYL_ORIENTATION_A])
      {
        /* Axial in x direction */

        case FET_ORIENTATION_X:
        {
          axial = x;
          pos0 =  y - center0;
          pos1 =  z - center1;

          break;
        }

        /* Axial in y direction */

        case FET_ORIENTATION_Y:
        {
          pos0 =  x - center0;
          axial = y;
          pos1 =  z - center1;

          break;
        }

        /* Axial in z direction */

        case FET_ORIENTATION_Z:
        {
          pos0 =  x - center0;
          pos1 =  y - center1;
          axial = z;

          break;
        }

#ifdef DEBUG
        default:
        {
          Die(FUNCTION_NAME, "Unknown orientation %ld", (long)params[FET_CYL_ORIENTATION_A]);
        }
#endif /* DEBUG */
      }

      /* Get axial position and check location */

      min1 = params[FET_CYL_MIN_A];
      max1 = params[FET_CYL_MAX_A];
      if (axial < min1 || max1 < axial)
        return INFTY;
      norm1 = 2.0 * (axial - min1) / (max1 - min1) - 1;
      deg1 = params[FET_CYL_ORDER_A];

      /* Get radius and check location */

      radius = params[FET_CYL_MAX_R];
      norm0 = sqrt(pos0 * pos0 + pos1 * pos1) / radius;
      if (1 <= norm0)
        return INFTY;
      theta = PolarAngle(pos0, pos1);
      deg0 = params[FET_CYL_ORDER_R];
      ncoef = params[FET_CYL_NCOEF_R];

      /* Get the locations of the coefficient arrays */

      vector0 = (long)params[FET_CYL_CALC_ARRAY_PTR_R];
      vector1 = (long)params[FET_CYL_CALC_ARRAY_PTR_A];

      /* Calculate the functional expansions. */

      PolynomialZernike(deg0, norm0, theta, vector0, threadID, FET_CALCULATION_STANDARD);
      PolynomialLegendre(deg1, norm1, vector1, threadID, FET_CALCULATION_STANDARD);

      /* Loop over all the coefficients to expand the series at this location */

      for (idx0 = 0, a = 0; idx0 < ncoef; idx0++) /* Zernike */
      {
        part0 = GetPrivateData(vector0 + idx0, threadID);

        for (idx1 = 0; idx1 <= deg1; a++, idx1++) /* Legendre */
          value += coefs[a] * part0 * GetPrivateData(vector1 + idx1, threadID);
      }

      break;
    }

#ifdef DEBUG
    default:
    {
      Die(FUNCTION_NAME, "Unsupported FET type %ld", (long)params[FET_PARAM_TYPE]);
    }
#endif /* DEBUG */
  }

  return value;
}
