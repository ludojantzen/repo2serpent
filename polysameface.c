/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : polysameface.c                                 */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Checks if two lists of points contain the same points        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PolySameFace:"

/*****************************************************************************/

long PolySameFace(long *face0, long *face1, long np)
{
  long i;

  /* If one of the points in face 0 is not in face 1 */
  /* return NO */

  for (i = 0; i < np; i++)
    if (!PolyPInF(face0[i], face1, np))
      return NO;

  /* Otherwise return YES */

  return YES;

}
