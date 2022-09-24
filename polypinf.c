/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : polypinf.c                                     */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Check if point index is in polyhdron face                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PolyPInF:"

/*****************************************************************************/

long PolyPInF(long idx, long *face, long np)
{
  long i;

  /* Loop over face points and check if one of them is equal to idx */

  for (i = 0; i < np; i++)
    if (idx == face[i])
      {
        /* Found a match */
        return YES;
      }

  /* Did not find a match*/

  return NO;
  
}
