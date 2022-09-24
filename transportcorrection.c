/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : transportcorrection.c                          */
/*                                                                           */
/* Created:       2016/06/11 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Returns material-wise transport correction factor            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TransportCorrection:"

/*****************************************************************************/

double TransportCorrection(long mat, double E, long id)
{
  long loc0, ng, nt, ptr, n, i;
  double f, T, d0, d1;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "mat", DATA_ARRAY, mat);

  /* Check if transport correction is used */

  if ((loc0 = (long)RDB[mat + MATERIAL_PTR_TRANSP_CORR]) < VALID_PTR)
    return 1.0;

  /* Get number of energy groups and temperatures */

  ng = (long)RDB[loc0 + TRANSP_CORR_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 1000000);

  nt = (long)RDB[loc0 + TRANSP_CORR_NT];
  CheckValue(FUNCTION_NAME, "nt", "", nt, 0, 1000);

  /* Pointer to energy array */

  ptr = (long)RDB[loc0 + TRANSP_CORR_PTR_ENE];
  CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

  /* Get energy group */

  if ((n = SearchArray(&RDB[ptr], E, ng + 1)) < 0)
    return 1.0;

  /* Check index */

  CheckValue(FUNCTION_NAME, "n", "", n, 0, ng);
  
  /* Check if multiple temperatures */

  if (nt > 1)
    {
      /* Pointer to temperature array */

      ptr = (long)RDB[loc0 + TRANSP_CORR_PTR_TEMP];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Get temperature (use first provided if temperature is not defined) */
             
      if ((T = GetTemp(mat, id)) < ZERO)
        if ((T = RDB[mat + MATERIAL_DOPPLER_TEMP]) < ZERO)
          T = RDB[ptr];

      /* Check limits */

      if (T == RDB[ptr])
        {
          i = 0;
          f = 0.0;
        }
      else if (T == RDB[ptr + nt - 1])
        {
          i = nt - 2;
          f = 1.0;
        }
      else 
        {
          /* Search interval */

          if ((i = SearchArray(&RDB[ptr], T, nt)) < 0)
            return 1.0;
      
          /* Check index and temperature */
          
          CheckValue(FUNCTION_NAME, "i", "", i, 0, nt - 2);
          CheckValue(FUNCTION_NAME, "T", "", T, RDB[ptr + i], RDB[ptr + i + 1]);

          /* Calculate factor */

          f = (T - RDB[ptr + i])/(RDB[ptr + i + 1] - RDB[ptr + i]);
          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
        }

      /* Pointer to data */

      ptr = (long)RDB[loc0 + TRANSP_CORR_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Get values */

      d0 = RDB[ptr + i*ng + n];
      d1 = RDB[ptr + (i + 1)*ng + n];

      /* Interpolate */

      f = f*(d1 - d0) + d0;
      CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 100.0);
    }
  else
    {
      /* Pointer to data */

      ptr = (long)RDB[loc0 + TRANSP_CORR_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Get value */

      f = RDB[ptr + n];
      CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 100.0);
    }

  /* Return factor */

  return f;
}

/*****************************************************************************/
