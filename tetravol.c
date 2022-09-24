/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : tetravol.c                                     */
/*                                                                           */
/* Created:       2015/02/19 (VVa)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Calculates volume of a tetrahedral unstructured mesh cell    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TetraVol:"

/*****************************************************************************/

double TetraVol(long tet)
{
  long p0, p1, p2, p3, prnt;
  double vol, a1, a2, a3, b1, b2, b3;
  double c1, c2, c3;

  /* Check tet pointer */

  CheckPointer(FUNCTION_NAME, "(tet)", DATA_ARRAY, tet);

  /* Get pointers to points */

  p0 = (long)RDB[tet + TET_POINTS + 0];
  p1 = (long)RDB[tet + TET_POINTS + 1];
  p2 = (long)RDB[tet + TET_POINTS + 2];
  p3 = (long)RDB[tet + TET_POINTS + 3];

  /* Calculate vector components */

  a1 = RDB[p1+0] - RDB[p0+0];
  a2 = RDB[p1+1] - RDB[p0+1];
  a3 = RDB[p1+2] - RDB[p0+2];

  b1 = RDB[p2+0] - RDB[p0+0];
  b2 = RDB[p2+1] - RDB[p0+1];
  b3 = RDB[p2+2] - RDB[p0+2];

  c1 = RDB[p3+0] - RDB[p0+0];
  c2 = RDB[p3+1] - RDB[p0+1];
  c3 = RDB[p3+2] - RDB[p0+2];

   /* Calculate vector triple product and volume */

  vol = a1*(b2*c3 - b3*c2) - a2*(b1*c3 - b3*c1) + a3*(b1*c2 - b2*c1);
  vol = fabs(vol)/6.0;

  /* Get pointer to parent (if exists) or self */

  prnt = (long)RDB[tet + TET_PTR_PARENT];

  /* Check */

  if (vol == 0.0)
    Warn(FUNCTION_NAME, "Zero volume in cell %ld",
        (long)RDB[prnt + IFC_TET_PRNT_IDX]);
  else if (vol > 1E+9)
    Warn(FUNCTION_NAME, "Suspiciously large volume %1.5E in cell %ld", vol,
        (long)RDB[prnt + IFC_TET_PRNT_IDX]);

  /* Return volume */

  return vol;

}

/*****************************************************************************/
