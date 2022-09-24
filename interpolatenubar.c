/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : interpolatenubar.c                             */
/*                                                                           */
/* Created:       2014/04/02 (JLe)                                           */
/* Last modified: 2018/11/13 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Creates temporary nubar array on unionized energy grid       */
/*                                                                           */
/* Comments: - Used for calculating pre-generated macroscopic nsf            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InterpolateNubar:"

/*****************************************************************************/

void InterpolateNubar(double *nubar, long rea)
{
  long loc0, erg, ptr, ne, n, np, m;
  double x, *dnubar;
  const double *E;

  /* Check mode and pointers */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == NO)
    Die(FUNCTION_NAME, "Macroscopic cross sections not reconstructed");
  else if (nubar == NULL)
    Die(FUNCTION_NAME, "Null pointer");

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get pointer to total nubar data */

  if ((loc0 = (long)RDB[rea + REACTION_PTR_TNUBAR]) < VALID_PTR)
    Die(FUNCTION_NAME, "No nubar data");

  /* Pointer to energy grid */

  ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of energy points */

  ne = (long)RDB[ptr + ENERGY_GRID_NE];

  /* Get pointer to energy array */

  ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  E = &RDB[ptr];

  /***************************************************************************/

  /***** Construct total nubar data ******************************************/

   /* Check type */

  if ((long)RDB[loc0 + NUBAR_DATA_TYPE] == 1)
    {
      /***********************************************************************/

      /***** Polynomial data *************************************************/

      /* Pointer to nubar data */

      ptr = (long)RDB[loc0 + NUBAR_PTR_POLY_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Number of coefficients */

      np = (long)RDB[ptr++];

      /* Loop over energy array */

      for (n = 0; n < ne; n++)
        {
          /* Initial value */

          nubar[n] = RDB[ptr];
          x = E[n];

          /* Loop over coefficients */

          for (m = 1; m < np; m++)
            {
              nubar[n] = nubar[n] + RDB[ptr + m]*x;
              x = x*E[n];
            }

          /* Check value */

          CheckValue(FUNCTION_NAME, "nubar", "", nubar[n], 1.5, 60.0);
        }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Tabular data ****************************************************/

      /* Pointer to energy grid */

      erg = (long)RDB[loc0 + NUBAR_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get number of points */

      np = (long)RDB[erg + ENERGY_GRID_NE];

      /* Pointer to energy array */

      erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Pointer to data */

      ptr = (long)RDB[loc0 + NUBAR_PTR_PTS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Interpolate data */

      InterpolateData(E, nubar, ne, &RDB[erg], &RDB[ptr], np, 0, NULL, NULL, NO);

      /* Fill leading zeros */

      /* Check values (HUOM: arvo voi olla nolla fission energiakynnysen */
      /* alapuolella. Huomaa myös kaarisulkeet) */

#ifdef DEBUG

      for (n = 0; n < ne; n++)
        if (E[n] >= RDB[erg])
          {
            /* Check value */

            CheckValue(FUNCTION_NAME, "nubar", "", nubar[n], 0.0, 60.0);
          }
#endif

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Subtract delayed nubar data *****************************************/

  /* Check flag */

  if ((long)RDB[DATA_USE_DELNU] != NO)
    return;

  /* Get pointer to delayed nubar data */

  if ((loc0 = (long)RDB[rea + REACTION_PTR_DNUBAR]) < VALID_PTR)
    return;

   /* Check type */

  if ((long)RDB[loc0 + NUBAR_DATA_TYPE] == 1)
    {
      /***********************************************************************/

      /***** Polynomial data *************************************************/

      /* Pointer to nubar data */

      ptr = (long)RDB[loc0 + NUBAR_PTR_POLY_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Number of coefficients */

      np = (long)RDB[ptr++];

      /* Loop over energy array */

      for (n = 0; n < ne; n++)
        {
          /* Initial value */

          nubar[n] = nubar[n] - RDB[ptr];
          x = E[n];

          /* Loop over coefficients */

          for (m = 1; m < np; m++)
            {
              nubar[n] = nubar[n] - RDB[ptr + m]*x;
              x = x*E[n];
            }

          /* Check value */

          CheckValue(FUNCTION_NAME, "nubar", "", nubar[n], 1.5, 60.0);
        }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Tabular data ****************************************************/

      /* Allocate memory for temporary array */

      dnubar = (double *)Mem(MEM_ALLOC, ne, sizeof(double));

      /* Pointer to energy grid */

      erg = (long)RDB[loc0 + NUBAR_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get number of points */

      np = (long)RDB[erg + ENERGY_GRID_NE];

      /* Pointer to energy array */

      erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Pointer to data */

      ptr = (long)RDB[loc0 + NUBAR_PTR_PTS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Interpolate data */

      InterpolateData(E, dnubar, ne, &RDB[erg], &RDB[ptr], np, 0, NULL, NULL, NO);

      /* Subtract data */

      for (n = 0; n < ne; n++)
        {
          /* Check with minimum and maximum energy */

          if ((E[n] < RDB[rea + REACTION_EMIN]) ||
              (E[n] > RDB[rea + REACTION_EMAX]))
            nubar[n] = 0.0;
          else
            {
              nubar[n] = nubar[n] - dnubar[n];

              /* Check value */

              CheckValue(FUNCTION_NAME, "nubar", "", nubar[n], 1.5, 60.0);
            }
        }

      /* Free temporary array */

      Mem(MEM_FREE, dnubar);

      /***********************************************************************/
    }

  /***************************************************************************/
}

/*****************************************************************************/
