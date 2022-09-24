/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshcellindirection.c                          */
/*                                                                           */
/* Created:       2015/08/31 (VVa)                                           */
/* Last modified: 2015/08/31 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: "Rotates" a hex or a prism so that the corner point with     */
/*              smallest index is in bottom left front/back corner           */
/*                                                                           */
/* Comments:   -Used for fixhexmesh.c                                        */
/*             -Based on "J. Dompierre, P. Labbe, M. Vallet and R. Camarero, */
/*              How to Subdivide Pyramids, Prisms, and Hexahedra             */ 
/*              into Tetrahedra, Proceedings,                                */
/*              8th International Meshing Roundtable,                        */
/*              South Lake Tahoe, CA, U.S.A., pp.195-204, October 1999       */
/*                                                                           */
/*                                                                           */
/*                            3---5                                          */
/*             7-----6        |\ /|                                          */
/*            /|    /|        | 4 |           4                              */
/*           4-----5 |        | | |          /|\                             */
/*           | 3---|-2        0-|-2        3-----2                           */
/*           |/    |/          \|/        /     /                            */
/*           0-----1            1        0-----1                             */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshCellIndirection:"

/*****************************************************************************/

void MeshCellIndirection(long *hex, long celltype )
{
  long temp[8], i, minidx, minpos;

  /* Set minimum point index to first hex point */

  minidx = hex[0];
  minpos = 0;

  /* Copy initial positions to temp */

  for (i = 0; i < 8; i++)
    {
      temp[i] = hex[i];
    }

  switch (celltype)
    {
    case MESH_CELL_TYPE_TET:
      return;

      break;

    case MESH_CELL_TYPE_PYRAMID:
      return;

      break;

    case MESH_CELL_TYPE_PRISM:

      /* Check if one of the other points has a smaller index */

      for (i = 0; i < 6; i++)
        if (hex[i] < minidx)
          {
            minidx = hex[i];
            minpos = i;
          }

      /* Rotate hex based on indirection table and minpos */

      switch (minpos)
        {
        case 0:
          /* Minimum index already at position 0 */
          /* We do not have to rotate            */

          break;

        case 1:

          /* Rotate */

          hex[0] = temp[1];
          hex[1] = temp[2];
          hex[2] = temp[0];
          hex[3] = temp[4];
          hex[4] = temp[5];
          hex[5] = temp[3];

          break;

        case 2:

          /* Rotate */

          hex[0] = temp[2];
          hex[1] = temp[0];
          hex[2] = temp[1];
          hex[3] = temp[5];
          hex[4] = temp[3];
          hex[5] = temp[4];

          break;

        case 3:

          /* Rotate */

          hex[0] = temp[3];
          hex[1] = temp[5];
          hex[2] = temp[4];
          hex[3] = temp[0];
          hex[4] = temp[2];
          hex[5] = temp[1];

          break;

        case 4:

          /* Rotate */

          hex[0] = temp[4];
          hex[1] = temp[3];
          hex[2] = temp[5];
          hex[3] = temp[1];
          hex[4] = temp[0];
          hex[5] = temp[2];

          break;

        case 5:

          /* Rotate */

          hex[0] = temp[5];
          hex[1] = temp[4];
          hex[2] = temp[3];
          hex[3] = temp[2];
          hex[4] = temp[1];
          hex[5] = temp[0];

          break;

        default:
          Die(FUNCTION_NAME, "bad smallest index");

          break;

        }

      break;

    case MESH_CELL_TYPE_HEX:

      /* Check if one of the other points has a smaller index */

      for (i = 0; i < 8; i++)
        if (hex[i] < minidx)
          {
            minidx = hex[i];
            minpos = i;
          }

      /* Rotate hex based on indirection table and minpos */

      switch (minpos)
        {
        case 0:
          /* Minimum index already at position 0 */
          /* We do not have to rotate            */

          break;

        case 1:

          /* Rotate */

          hex[0] = temp[1];
          hex[1] = temp[0];
          hex[2] = temp[4];
          hex[3] = temp[5];
          hex[4] = temp[2];
          hex[5] = temp[3];
          hex[6] = temp[7];
          hex[7] = temp[6];       

          break;

        case 2:

          /* Rotate */

          hex[0] = temp[2];
          hex[1] = temp[1];
          hex[2] = temp[5];
          hex[3] = temp[6];
          hex[4] = temp[3];
          hex[5] = temp[0];
          hex[6] = temp[4];
          hex[7] = temp[7];       

          break;

        case 3:

          /* Rotate */

          hex[0] = temp[3];
          hex[1] = temp[0];
          hex[2] = temp[1];
          hex[3] = temp[2];
          hex[4] = temp[7];
          hex[5] = temp[4];
          hex[6] = temp[5];
          hex[7] = temp[6];       

          break;

        case 4:

          /* Rotate */

          hex[0] = temp[4];
          hex[1] = temp[0];
          hex[2] = temp[3];
          hex[3] = temp[7];
          hex[4] = temp[5];
          hex[5] = temp[1];
          hex[6] = temp[2];
          hex[7] = temp[6];       

          break;

        case 5:

          /* Rotate */

          hex[0] = temp[5];
          hex[1] = temp[1];
          hex[2] = temp[0];
          hex[3] = temp[4];
          hex[4] = temp[6];
          hex[5] = temp[2];
          hex[6] = temp[3];
          hex[7] = temp[7];       

          break;

        case 6:

          /* Rotate */

          hex[0] = temp[6];
          hex[1] = temp[2];
          hex[2] = temp[1];
          hex[3] = temp[5];
          hex[4] = temp[7];
          hex[5] = temp[3];
          hex[6] = temp[0];
          hex[7] = temp[4];       

          break;

        case 7:

          /* Rotate */

          hex[0] = temp[7];
          hex[1] = temp[3];
          hex[2] = temp[2];
          hex[3] = temp[6];
          hex[4] = temp[4];
          hex[5] = temp[0];
          hex[6] = temp[1];
          hex[7] = temp[5];       

          break;

        default:
          Die(FUNCTION_NAME, "bad smallest index");

          break;

        }

      break;

    default:
      Die(FUNCTION_NAME, "Non-supported cell type %ld", celltype);

    }
  return;

}
