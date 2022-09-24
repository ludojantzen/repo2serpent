/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hexrotateface.c                                */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2016/02/17 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Rotates a face consisting of four points                     */
/*                                                                           */
/* Comments:   -Rotation of face 0-1-2-3 by +1 gives 3-0-1-2                 */
/*             -Rotation of face 0-1-2-3 by -1 gives 1-2-3-0                 */
/*                                                                           */
/*             7-----6                                                       */
/*            /|    /|                                                       */
/*           4-----5 |                                                       */
/*           | 3---|-2                                                       */
/*           |/    |/                                                        */
/*           0-----1                                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FixHexMesh:"

/*****************************************************************************/

void HexRotateFace(long *face, long steps, long np)
{
  long temp;
  long i, j;

  if (steps > 0)
    {

      /* Positive rotation */
      /* Rotate for requested number of steps */

      for (i = 0; i < steps; i++)
        {
          temp = face[np - 1];
          for (j = np - 1; j > 0; j--)
              face[j] = face[j-1];
          /*
              face[3] = face[2];
              face[2] = face[1];
              face[1] = face[0];
              face[0] = temp;
          */
          face[0] = temp;
        }
    }
  else
    {
      /* Negative rotation */

      steps = -steps;

      /* Rotate for requested number of steps */

      for (i = 0; i < steps; i++)
        {
          temp = face[0];

          for (j = 0; j < np - 1; j++)
            face[j] = face[j+1];

          face[np-1] = temp;
        }
    }
}
