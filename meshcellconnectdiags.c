/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshcellconnectdiags.c                         */
/*                                                                           */
/* Created:       2015/08/31 (VVa)                                           */
/* Last modified: 2015/08/31 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Connects a diagonal to faces based on point-indices          */
/*              of face points                                               */
/*                                                                           */
/* Comments:    Diagonal on e.g. 0 1 4 3 face can be 0-4 or 1-3              */
/*              the diagonal that contains the smallest point index is       */
/*              chosen, this way the same diagonal is always chosen          */
/*             -Based on "J. Dompierre, P. Labbe, M. Vallet and R. Camarero, */
/*              How to Subdivide Pyramids, Prisms, and Hexahedra             */ 
/*              into Tetrahedra, Proceedings,                                */
/*              8th International Meshing Roundtable,                        */
/*              South Lake Tahoe, CA, U.S.A., pp.195-204, October 1999       */
/*                                                                           */
/*                              |                                            */
/*                              v                                            */
/*                            3---5                                          */
/*             7-----6        |\ /|                                          */
/*            /|    /|        | 4 |           4                              */
/*           4-----5 |        | | |          /|\                             */
/*           | 3---|-2        0-|-2        3-----2                           */
/*           |/    |/          \|/        /     /                            */
/*           0-----1            1        0-----1                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshCellConnectDiags:"

/*****************************************************************************/

void MeshCellConnectDiags(long hex[8], long (*diags)[8], long celltype)
{
  long pos[4], i00, i01, i10, i11, i;

  if ((celltype == MESH_CELL_TYPE_TET) || 
      (celltype == MESH_CELL_TYPE_PYRAMID))
    return;

  for (i = 0; i < 6; i++)
    {
      /* Prisms do not have a frontface, and bot and top don't have */
      /* a diagonal, skip those faces */

      if (celltype == MESH_CELL_TYPE_PRISM)
        if ((i == MESH_CELL_FACE_FRONT) || (i == MESH_CELL_FACE_BOTTOM) || (i == MESH_CELL_FACE_TOP))
          continue;

      /* Get node positions for face i */
          
      MeshCellGetFacePos(pos, i, celltype);

      /* Get first possible diagonal */
      i00 = pos[0];
      i01 = pos[2];

      /* Get second possible diagonal */
      i10 = pos[1];
      i11 = pos[3];

      /* If smaller index of first possible diagonal      */
      /* is smaller than that of second possible diagonal */

      if ((hex[i00] < hex[i01] ? hex[i00] : hex[i01]) 
          < (hex[i10] < hex[i11] ? hex[i10] : hex[i11]))
        {
          /* Make connection to first diagonal */

          diags[i00][i01] = 1;
          diags[i01][i00] = 1;

        }
      else if ((hex[i00] < hex[i01] ? hex[i00] : hex[i01]) 
               > (hex[i10] < hex[i11] ? hex[i10] : hex[i11]))
        {
          /* Otherwise make connection to second diagonal */

          diags[i10][i11] = 1;
          diags[i11][i10] = 1;

        }
      else
        Die("FUNCTION_NAME", "Two same indices in hex?!");

    }
}
