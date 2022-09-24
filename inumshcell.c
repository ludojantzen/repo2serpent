/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : inumshcell.c                                   */
/*                                                                           */
/* Created:       2018/01/26 (VVa)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Checks if point is inside an umsh cell                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InUMSHCell:"

/*****************************************************************************/

long InUMSHCell(long ifc, long cgns, double x, double y, double z, long on)
{
  long n, ptr, surf, side, facelist, surflist, sidelist, retval;
  long pt0, pt1, pt2, nf, np;
  long i;
  double x0, x1, x2, y0, y1, y2, z0, z1, z2, d;
  double dotp, len1sq, len2sq, cossq, toler;

  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Get pointer to interface surfaces */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Get pointer to cell's face list */

  facelist = (long)RDB[cgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(facelist)", DATA_ARRAY, facelist);

  /* Get pointer to cell's side list */

  sidelist = (long)RDB[cgns + IFC_TET_MSH_PTR_SIDES];
  CheckPointer(FUNCTION_NAME, "(sidelist)", DATA_ARRAY, sidelist);

  /* Loop over cell faces */

  nf = (long)RDB[cgns + IFC_TET_MSH_NF];
  /*
  CheckValue(FUNCTION_NAME, "nf", "", nf, 4, 4);
  */
  for (i = 0; i < nf; i++)
    {

      /* Get index of face */

      n = (long)RDB[facelist + i];

      /* Get side for this face */

      side = (long)RDB[sidelist + i];

      /* Get pointer to face surface */

      surf = (long)RDB[surflist + n];

      /*      printf("face %ld, side %ld\n", n, side);*/

      /* Get number of points on the face */

      np = (long)RDB[surf + UMSH_SURF_N_POINTS];

      /* Get pointer to surface points */

      ptr = (long)RDB[surf + UMSH_SURF_PTR_POINTS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointers to points */

      pt0 = (long)RDB[ptr + 0];
      pt1 = (long)RDB[ptr + 1];
      pt2 = (long)RDB[ptr + 2];

      /* Three points, calculate vectors */

      x0 = -RDB[pt1 + 0] + x;
      y0 = -RDB[pt1 + 1] + y;
      z0 = -RDB[pt1 + 2] + z;

      x1 = -RDB[pt0 + 0] + RDB[pt1 + 0];
      y1 = -RDB[pt0 + 1] + RDB[pt1 + 1];
      z1 = -RDB[pt0 + 2] + RDB[pt1 + 2];

      x2 = -RDB[pt1 + 0] + RDB[pt2 + 0];
      y2 = -RDB[pt1 + 1] + RDB[pt2 + 1];
      z2 = -RDB[pt1 + 2] + RDB[pt2 + 2];

      /* Check angle between perimeter vectors */
      /* (too small angle can yield erroneous results for the check) */

      dotp = x1*x2 + y1*y2 + z1*z2;

      if (dotp > 0)
        {
          /* For negative dot-products the angle is above 90 degrees */
          /* Positive ones may be small */

          /* Get square of perimeter vector lengths */

          len1sq = x1*x1 + y1*y1 + z1*z1;
          len2sq = x2*x2 + y2*y2 + z2*z2;

          /* Based on definition of dot product calculate square of */
          /* cosine */

          cossq = dotp*dotp/(len1sq*len2sq);

          /* Test square of cosine */

          /* Tolerance of 1e-6 corresponds to angles smaller than 0.057 degrees */

          toler = 1e-6;

          if (cossq > 1 - toler)
            {
              /* Get next two points from perimeter to calculate the vectors from */

              /* Get pointers to points */

              pt0 = (long)RDB[ptr + np-1];
              pt1 = (long)RDB[ptr + 0];
              pt2 = (long)RDB[ptr + 1];

              /* Three points, calculate vectors */

              x0 = -RDB[pt1 + 0] + x;
              y0 = -RDB[pt1 + 1] + y;
              z0 = -RDB[pt1 + 2] + z;

              x1 = -RDB[pt0 + 0] + RDB[pt1 + 0];
              y1 = -RDB[pt0 + 1] + RDB[pt1 + 1];
              z1 = -RDB[pt0 + 2] + RDB[pt1 + 2];

              x2 = -RDB[pt1 + 0] + RDB[pt2 + 0];
              y2 = -RDB[pt1 + 1] + RDB[pt2 + 1];
              z2 = -RDB[pt1 + 2] + RDB[pt2 + 2];

              /* NB: These vectors are not checked for a small angle! */
              /* TODO: Recursive checking until two acceptable vectors are */
              /* found? */
            }
        }

      /* Scalar triple product */

      d = x0*(y1*z2 - y2*z1) - y0*(x1*z2 - x2*z1) + z0*(x1*y2 - x2*y1);

      /* Check if surface itself is included as inside (NOTE: this is */
      /* used with tet-mesh structures to avoid undefined regions on  */
      /* cell boundaries, otherwise the boundary is included in the   */
      /* outside of the surface) */

      if (on == YES)
        {
          /* Check */

          if (d <= 0.0)
            retval = YES;
          else
            retval = NO;
        }
      else
        {
          /* Check */

          if (d < 0.0)
            retval = YES;
          else
            retval = NO;
        }

      /* Test surface */

      if (retval == YES)
        side = -side;

      /* Check result */

      if (side < 0)
        return NO;
    }

  /* Point is inside all surfaces */

  return YES;

}

/*****************************************************************************/
