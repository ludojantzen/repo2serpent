/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : allocfetcache.c                                */
/*                                                                           */
/* Created:       2018/02/08 (BWe)                                           */
/* Last modified: 2018/02/23 (BWe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Allocates the memory for FET calculation caches.             */
/*                                                                           */
/* Comments: - ei mitään vielä                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AllocFETCache:"

void AllocFETCache(long params, long usage, long detector, long update)
{
  long bufferSize, localCoef, totalCoef, type;
  long ptr;

  switch (type = (long)RDB[params + FET_PARAM_TYPE])
  {
    case FET_TYPE_CARTESIAN:
    {
      localCoef = (long)RDB[params + FET_CART_ORDER_X];
      WDB[params + FET_CART_NCOEF_X] =
        (localCoef == -1 ? FET_DEFAULT_ORDER + 1 : localCoef + 1);

      localCoef = (long)RDB[params + FET_CART_ORDER_Y];
      WDB[params + FET_CART_NCOEF_Y] =
        (localCoef == -1 ? FET_DEFAULT_ORDER + 1 : localCoef + 1);

      localCoef = (long)RDB[params + FET_CART_ORDER_Z];
      WDB[params + FET_CART_NCOEF_Z] =
        (localCoef == -1 ? FET_DEFAULT_ORDER + 1 : localCoef + 1);

      /* Calculate the detector volumes */

      WDB[params + FET_PARAM_FET_VOLUME] = 2.0*2.0*2.0;
      WDB[params + FET_PARAM_PHYS_VOLUME] =
                (RDB[params + FET_CART_MAX_X] - RDB[params + FET_CART_MIN_X])
                * (RDB[params + FET_CART_MAX_Y] - RDB[params + FET_CART_MIN_Y])
                * (RDB[params + FET_CART_MAX_Z] - RDB[params + FET_CART_MIN_Z]);

      break;
    }

    case FET_TYPE_CYLINDRICAL:
    {
      /* Zernike polynomial orders follow the triangle number relationship. */
      /* Add one to the terms for 'n' in the formula to account for the     */
      /* zeroth-order polynomial in the number of coefficients.             */

      localCoef = (long)RDB[params + FET_CYL_ORDER_R];
      localCoef = (localCoef == -1 ? FET_DEFAULT_ORDER + 1 : localCoef + 1);
      WDB[params + FET_CYL_NCOEF_R] = (localCoef * (localCoef + 1)) / 2;

      /* Legendre polynomials in the axial direction */

      localCoef = (long)RDB[params + FET_CYL_ORDER_A];
      WDB[params + FET_CYL_NCOEF_A] = (localCoef == -1 ? FET_DEFAULT_ORDER + 1 : localCoef + 1);

      /* Calculate the detector volumes */

      WDB[params + FET_PARAM_FET_VOLUME] = M_PI*2.0;
      WDB[params + FET_PARAM_PHYS_VOLUME] =
              M_PI*(RDB[params + FET_CYL_MAX_R]*RDB[params + FET_CYL_MAX_R])
              * (RDB[params + FET_CYL_MAX_A] - RDB[params + FET_CYL_MIN_A]);

      break;
    }

#ifdef DEBUG
    default:
    {
      Die(FUNCTION_NAME, "Unsupported FET type %ld", type);
    }
#endif /* DEBUG */
  }

  /* Load the detector volume into the detector parameters, if needed */

  if (detector > 0)
    WDB[detector + DET_VOL] = RDB[params + FET_PARAM_PHYS_VOLUME];

  /* Allocate memory for temporary arrays and set the total number of bins.   */

  totalCoef = 1;

  /* Check coefficient 0 */

  localCoef = (long)RDB[params + FET_PARAM_NCOEF0];
  if (localCoef > 0)
  {
    totalCoef *= localCoef;

    if (localCoef > (long)RDB[params + FET_PARAM_CALC_ARRAY_SIZE0])
    {
#ifdef DEBUG
      if ((update == YES) && ((long)RDB[params + FET_PARAM_CALC_ARRAY_SIZE0] > 0))
        Warn(FUNCTION_NAME, "Resizing the FET calculation array for series 0; the previous array is being discarded but no deallocated");
#endif /* DEBUG */
      bufferSize = ceil(localCoef * 1.5);
      ptr = AllocPrivateData(bufferSize, PRIVA_ARRAY);
      WDB[params + FET_PARAM_CALC_ARRAY_SIZE0] = (double)bufferSize;
      WDB[params + FET_PARAM_CALC_ARRAY_PTR0] = (double)ptr;
    }
  }

  /* Check coefficient 1 */

  localCoef = (long)RDB[params + FET_PARAM_NCOEF1];
  if (localCoef > 0)
  {
    totalCoef *= localCoef;

    if (localCoef > (long)RDB[params + FET_PARAM_CALC_ARRAY_SIZE1])
    {
#ifdef DEBUG
      if ((update == YES) && ((long)RDB[params + FET_PARAM_CALC_ARRAY_SIZE1] > 0))
        Warn(FUNCTION_NAME, "Resizing the FET calculation array for series 1; the previous array is being discarded but no deallocated");
#endif /* DEBUG */
      bufferSize = ceil(localCoef * 1.5);
      ptr = AllocPrivateData(bufferSize, PRIVA_ARRAY);
      WDB[params + FET_PARAM_CALC_ARRAY_SIZE1] = (double)bufferSize;
      WDB[params + FET_PARAM_CALC_ARRAY_PTR1] = (double)ptr;
    }
  }

  /* Check coefficient 2 */

  localCoef = (long)RDB[params + FET_PARAM_NCOEF2];
  if (localCoef > 0)
  {
    totalCoef *= localCoef;

    if (localCoef > (long)RDB[params + FET_PARAM_CALC_ARRAY_SIZE2])
    {
#ifdef DEBUG
      if ((update == YES) && ((long)RDB[params + FET_PARAM_CALC_ARRAY_SIZE2] > 0))
        Warn(FUNCTION_NAME, "Resizing the FET calculation array for series 2; the previous array is being discarded but no deallocated");
#endif /* DEBUG */
      bufferSize = ceil(localCoef * 1.5);
      ptr = AllocPrivateData(bufferSize, PRIVA_ARRAY);
      WDB[params + FET_PARAM_CALC_ARRAY_SIZE2] = (double)bufferSize;
      WDB[params + FET_PARAM_CALC_ARRAY_PTR2] = (double)ptr;
    }
  }

  WDB[params + FET_PARAM_NCOEF_TOTAL] = totalCoef;

  if (usage == FET_USE_COEF)
    return; /* Our work here is done, no need to generate uncertainties. */
  else if (totalCoef + 1 > (long)RDB[params + FET_PARAM_UNC_ARRAY_SIZE])
  {
    /* Allocate the space for processing FET uncertainties. Unfortunately the */
    /* current statistics processing scheme does not provide the required     */
    /* per-particle granularity, so this alternate fix has been implemented.  */
    /*                                                                        */
    /* Specifically, we will re-purpose of each buffer value. One is to track */
    /* the contributions of each event, the second to track the contributions */
    /* from the current particle, and the third is to track the sum of the    */
    /* squares of the per-particle contributions.                             */

#ifdef DEBUG
    if ((update == YES) && ((long)RDB[params + FET_PARAM_UNC_ARRAY_SIZE] > 0))
      Warn(FUNCTION_NAME, "Resizing the FET uncertainty array; the previous array is being discarded but no deallocated");
#endif /* DEBUG */

    if ((long)RDB[DATA_OMP_MAX_THREADS] > 1
        && (long)RDB[DATA_OPTI_SHARED_BUF] == YES)
      Die("ProcessDetectors", "OMP shared buffers cannot be used with FETs. "
                              "Change the input option \"set shbuf\" to \"0\"");

    /* One extra detector bin is used to count the number of particles        */
    /* contributing to the FET. One extra private data slot is used to track  */
    /* and the current particle ID.                                           */

    bufferSize = ceil(totalCoef * 1.5) + 1;

    ptr = AllocPrivateData(bufferSize, PRIVA_ARRAY);
    WDB[params + FET_PARAM_UNC_ARRAY_SIZE] = (double)bufferSize;
    WDB[params + FET_PARAM_UNC_ARRAY_PTR] = (double)ptr;
  }
  else
  {
    ptr = RDB[params + FET_PARAM_UNC_ARRAY_PTR];
    memset(&WDB[ptr],
           0.0,
           (long)RDB[params + FET_PARAM_UNC_ARRAY_SIZE]*sizeof(double));
  }
}
