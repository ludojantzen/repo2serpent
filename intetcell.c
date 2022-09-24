/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : intetcell.c                                    */
/*                                                                           */
/* Created:       2010/10/12 (JLe)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Checks if point is inside a tet cell                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InTetCell:"

/*****************************************************************************/

long InTetCell(long tet, double x, double y, double z, long on)
{
  long pt0, pt1, pt2;
  long i;
  long *face[4];
  long face1[3] = {0,2,1};
  long face2[3] = {1,2,3};
  long face3[3] = {0,3,2};
  long face4[3] = {0,1,3};
  double x0, x1, x2, y0, y1, y2, z0, z1, z2, d;

  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(tet)", DATA_ARRAY, tet);

  /***********************************************************/
  /* Numbering of tet faces:                                 */
  /* First face  (0,2,1) is out of cell (out of face).       */
  /* Second face (1,2,3) is forward on old face perimeter    */
  /* Third face  (0,3,2) is backward on old face perimeter   */
  /* Fourth face (0,1,3) is inside cell but not on this face */
  /***********************************************************/

  /* Store point list for easy access */

  face[0] = face1;
  face[1] = face2;
  face[2] = face3;
  face[3] = face4;

  for (i = 0; i < 4; i++)
    {

      /* Get pointers to points */

      pt0 = (long)RDB[tet + TET_POINTS + face[i][0]];
      pt1 = (long)RDB[tet + TET_POINTS + face[i][1]];
      pt2 = (long)RDB[tet + TET_POINTS + face[i][2]];

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

      /* Scalar triple product */

      d = x0*(y1*z2 - y2*z1) - y0*(x1*z2 - x2*z1) + z0*(x1*y2 - x2*y1);

      /* Check */

      if (d > 0.0)
        return NO;

    }

  /* Point is inside all surfaces */

  return YES;

}

/*****************************************************************************/
