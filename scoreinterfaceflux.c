/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoreinterfaceflux.c                           */
/*                                                                           */
/* Created:       2013/06/26 (VVa)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Scores fast flux for multi-physics interface                 */
/*                                                                           */
/* Comments: - Polttoaineinterfacen aksiaalijako lisätty 3.4.2013            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreInterfaceFlux:"

/*****************************************************************************/

void ScoreInterfaceFlux(double wgt, double flx, double E, double x, double y,
                        double z, double t, long id)
{
  long loc0, loc1, ptr, nz, na, nr, i, j, k, uni;
  double r2, phi;

  /* Check that interfaces are defined */

  if ((loc0 = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  /***************************************************************************/

  /***** Interface to fuel performance codes *********************************/

  /* Get collision universe */

  ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
  uni = (long)GetPrivateData(ptr, id);

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Check pointer to interface */

  if ((loc1 = (long)RDB[uni + UNIVERSE_PTR_IFC_FUEP]) > VALID_PTR)
    {
      /* Get coordinates */

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_X];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      x = GetPrivateData(ptr, id);

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Y];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      y = GetPrivateData(ptr, id);

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Z];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      z = GetPrivateData(ptr, id);

      ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_T];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      t = GetPrivateData(ptr, id);

      /* Get the size of the statistics */

      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FLIM];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      nz = (long)RDB[ptr + FUEP_NZ];
      na = (long)RDB[ptr + FUEP_NA];
      nr = (long)RDB[ptr + FUEP_NR];

      /* Check energy limits for scoring */

      if((E < RDB[ptr + FUEP_EMIN]) || (E > RDB[ptr + FUEP_EMAX]))
          return;

      /* Coordinate transformation to cold system */

      CoordExpans(loc1, &x, &y, &z, t, 1);

      /* Get pointer to axial output zones */

      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FZ];
      CheckPointer(FUNCTION_NAME, "(Zptr)", DATA_ARRAY, ptr);

      /* Find interval */

      i = SearchArray(&RDB[ptr], z, nz + 1);

      /* Check */

      if (i < 0)
        {

          return;
        }

      /* Calculate angle */

      phi = PolarAngle(x, y);

      /* Check phi */

      CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 2.0*PI);

      /* Get pointer to angular output zones */

      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FPHI];
      CheckPointer(FUNCTION_NAME, "(Aptr)", DATA_ARRAY, ptr);

      /* Rotate if needed */

      if (phi > 2.0*PI+RDB[ptr])
        phi -= 2.0*PI;

      /* Find interval */

      j = SearchArray(&RDB[ptr], phi, na + 1);

      /* Check */

      if (j < 0)
        {

          return;
        }

      /* Calculate square radius */

      r2 = x*x + y*y;

      /* Get pointer to radial output zones */

      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FR2];
      CheckPointer(FUNCTION_NAME, "(Rptr)", DATA_ARRAY, ptr);

      /* Find interval */

      k = SearchArray(&RDB[ptr], r2, nr + 1);

      /* Check */

      if (k < 0)
        {

          return;
        }

      /* Score */

      ptr = (long)RDB[loc1 + IFC_FUEP_PTR_FLUX];
      CheckPointer(FUNCTION_NAME, "(Fptr)", DATA_ARRAY, ptr);
      AddBuf(flx, wgt, ptr, id, -1, i, j, k);

    }

  /***************************************************************************/
}

/*****************************************************************************/
