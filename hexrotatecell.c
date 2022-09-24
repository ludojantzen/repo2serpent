/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hexrotatecell.c                                */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: "Rotates" hex around 0<->6 diagonal so that the single face  */
/*               with/without diagonal connected to 6 is on the right        */
/*                                                                           */
/* Comments:                                                                 */
/*             -Based on "J. Dompierre, P. Labbe, M. Vallet and R. Camarero, */
/*              How to Subdivide Pyramids, Prisms, and Hexahedra             */ 
/*              into Tetrahedra, Proceedings,                                */
/*              8th International Meshing Roundtable,                        */
/*              South Lake Tahoe, CA, U.S.A., pp.195-204, October 1999       */
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

#define FUNCTION_NAME "HexRotateCell:"

/*****************************************************************************/

void HexRotateCell(long *hex, long (*diags)[8])
{
  long nd, i, rot, temp;

  /* Check number of diagonal connections to position 6       */
  /* Possible connections with 1 (right), 3 (back) or 4 (top) */

  nd = 0;
  
  for (i = 0; i < 8; i++)
    {
      /* If diagonals matrix is set to 1 for [6][i] */
      /* there is a connection between the two      */

      if (diags[6][i] == 1)
        nd++;
    }

  /* nd now contains the number of face diagonals connected to 6 */
  /* (0 <= nd <= 3) */

  /* Get number of rotations needed based on nd */
  /* and the position of face diagonals         */

  rot = 0;

  if ((nd == 0) || (nd == 3))
    {
      /* No need to rotate if 0 or 3 diagonals */
      
      return;

    }
  else if (nd == 1)
    {
      /* Only one diagonal */

      /* Get number of needed rotations */

      if (diags[6][1] == 1)
        return;  /* diag already on the right ->no rotation needed */
      else if (diags[6][3] == 1)
        rot = 2; /* diag on the back -> two rotations of 120 deg  */
      else if (diags[6][4] == 1)
        rot = 1; /* diag on the top -> one rotation of 120 deg   */
      else
        Die(FUNCTION_NAME, "Strange connection 1 diag");
                
    }
  else if (nd == 2)
    {
      /* Only one side without diagonal */

      if (diags[6][1] == 0)
        return;  /* no diag on the right -> no rotation needed */
      else if (diags[6][3] == 0)
        rot = 2; /* no diag on the back -> two rotations of 120 deg */
      else if (diags[6][4] == 0)
        rot = 1; /* no diag on the top  one rotation of 120 deg */
      else
        Die(FUNCTION_NAME, "Strange connection 2 diag");

    }
  else
    Die(FUNCTION_NAME, "Strange number of diagonals %ld", nd);

  /* Do rotations */

  for (i = 0; i < rot; i++)
    {

      temp   = hex[1];
      hex[1] = hex[4];
      hex[4] = hex[3];
      hex[3] = temp;

      temp   = hex[5];
      hex[5] = hex[7];
      hex[7] = hex[2];
      hex[2] = temp;

    }

  return;

}

