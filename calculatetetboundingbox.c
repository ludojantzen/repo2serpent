/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculatetetboundingbox.c                      */
/*                                                                           */
/* Created:       2018/01/26 (VVa)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Calculates bounding box of tetrahedron based on its points   */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateTetBoundingBox:"

/*****************************************************************************/

void CalculateTetBoundingBox(long tet, double bb[6])
{
  long pt0, pt1, pt2, pt3;

  /* Reset bounding box for cell */

  bb[0] =  INFTY;
  bb[1] = -INFTY;
  bb[2] =  INFTY;
  bb[3] = -INFTY;
  bb[4] =  INFTY;
  bb[5] = -INFTY;

  /* Get pointers to points */

  pt0 = (long)RDB[tet + TET_POINTS + 0];
  pt1 = (long)RDB[tet + TET_POINTS + 1];
  pt2 = (long)RDB[tet + TET_POINTS + 2];
  pt3 = (long)RDB[tet + TET_POINTS + 3];

  /* Put bounding box */

  if (bb[0] > RDB[pt0 + 0])
    bb[0]   = RDB[pt0 + 0];
  if (bb[1] < RDB[pt0 + 0])
    bb[1]   = RDB[pt0 + 0];

  if (bb[2] > RDB[pt0 + 1])
    bb[2]   = RDB[pt0 + 1];
  if (bb[3] < RDB[pt0 + 1])
    bb[3]   = RDB[pt0 + 1];

  if (bb[4] > RDB[pt0 + 2])
    bb[4]   = RDB[pt0 + 2];
  if (bb[5] < RDB[pt0 + 2])
    bb[5]   = RDB[pt0 + 2];

  /***/

  if (bb[0] > RDB[pt1 + 0])
    bb[0]   = RDB[pt1 + 0];
  if (bb[1] < RDB[pt1 + 0])
    bb[1]   = RDB[pt1 + 0];

  if (bb[2] > RDB[pt1 + 1])
    bb[2]   = RDB[pt1 + 1];
  if (bb[3] < RDB[pt1 + 1])
    bb[3]   = RDB[pt1 + 1];

  if (bb[4] > RDB[pt1 + 2])
    bb[4]   = RDB[pt1 + 2];
  if (bb[5] < RDB[pt1 + 2])
    bb[5]   = RDB[pt1 + 2];

  /***/

  if (bb[0] > RDB[pt2 + 0])
    bb[0]   = RDB[pt2 + 0];
  if (bb[1] < RDB[pt2 + 0])
    bb[1]   = RDB[pt2 + 0];

  if (bb[2] > RDB[pt2 + 1])
    bb[2]   = RDB[pt2 + 1];
  if (bb[3] < RDB[pt2 + 1])
    bb[3]   = RDB[pt2 + 1];

  if (bb[4] > RDB[pt2 + 2])
    bb[4]   = RDB[pt2 + 2];
  if (bb[5] < RDB[pt2 + 2])
    bb[5]   = RDB[pt2 + 2];

  /***/

  if (bb[0] > RDB[pt3 + 0])
    bb[0]   = RDB[pt3 + 0];
  if (bb[1] < RDB[pt3 + 0])
    bb[1]   = RDB[pt3 + 0];

  if (bb[2] > RDB[pt3 + 1])
    bb[2]   = RDB[pt3 + 1];
  if (bb[3] < RDB[pt3 + 1])
    bb[3]   = RDB[pt3 + 1];

  if (bb[4] > RDB[pt3 + 2])
    bb[4]   = RDB[pt3 + 2];
  if (bb[5] < RDB[pt3 + 2])
    bb[5]   = RDB[pt3 + 2];
}
