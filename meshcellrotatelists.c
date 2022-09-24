/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshcellrotatelists.c                          */
/*                                                                           */
/* Created:       2015/08/31 (VVa)                                           */
/* Last modified: 2017/11/29 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Rotates the auxilliary directional lists for mesh cells      */
/*              after the cell itself has been rotated to final position     */
/*                                                                           */
/* Comments:   -Needed for fixhexmesh.c                                      */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*                            3---5                                          */
/*             7-----6        |\ /|                                          */
/*            /|    /|        | 4 |           4            3                 */
/*           4-----5 |        | | |          /|\          /|\                */
/*           | 3---|-2        0-|-2        3-----2       0-|-2               */
/*           |/    |/          \|/        /     /         \|/                */
/*           0-----1            1        0-----1           1                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshCellRotateLists:"

/*****************************************************************************/

void MeshCellRotateLists(long hex[8], long (*initFaces)[4], long *hexfaces,
                    long celltype)
{
  long i, j, np, maxf;
  long tmpfaces[6], face0[4], *face1;

  /* Store all the lists to temporary lists */

  for (i = 0; i < 6; i++)
    {
      tmpfaces[i] = hexfaces[i];
    }

  /* Loop over face directions */

  for (i = 0; i < 6; i++)
    {
      /* Tets do not have a top face or a front face */

      if (celltype == MESH_CELL_TYPE_TET)
        {
          /* Tets do not have a top or front face */

          if (( i == MESH_CELL_FACE_TOP ) || ( i == MESH_CELL_FACE_FRONT ))
            {
              hexfaces[i] = -1;

              continue;
            }
        }
      else if (celltype == MESH_CELL_TYPE_PYRAMID)
        {
          /* Pyramids do not have a top face */

          if ( i == MESH_CELL_FACE_TOP )
            {
              hexfaces[i] = -1;

              continue;
            }

        }
      else if (celltype == MESH_CELL_TYPE_PRISM)
        {
          /* Prisms do not have a front face */

          if ( i == MESH_CELL_FACE_FRONT )
            {
              hexfaces[i] = -1;

              continue;
            }

        }

      /* Get next face from hex */

      MeshCellGetFace(hex, face0, i, celltype);

      /* Get number of initial faces */
      /* And number of points in this face */

      switch (celltype)
        {
        case MESH_CELL_TYPE_TET:

          maxf = 4;
          np = 3;

          break;
        case MESH_CELL_TYPE_PYRAMID:

          maxf = 5;

          if (i == MESH_CELL_FACE_BOTTOM)
            np = 4;
          else
            np = 3;


          break;
        case MESH_CELL_TYPE_PRISM:

          maxf = 5;

          if ((i == MESH_CELL_FACE_TOP) || (i == MESH_CELL_FACE_BOTTOM))
            np = 3;
          else
            np = 4;

          break;
        case MESH_CELL_TYPE_HEX:

          maxf = 6;
          np = 4;

          break;
        default:

          maxf = -1;
          np = -1;

          Die(FUNCTION_NAME, "Invalid cell type %ld", celltype);

        }

      /* Loop over initial faces and find match */

      for (j = 0; j < maxf; j++)
        {

          face1 = initFaces[j];

          /* Check if faces match */

          if (PolySameFace(face0,face1, np))
            {
              /* If faces matched store list data for this face */

              hexfaces[i] = tmpfaces[j];

              break;
            }

        }

      /* If this rotated face did not match any of the initial faces */
      /* something is wrong in fixhexmesh.c                          */

      if (j == maxf)
        Die(FUNCTION_NAME, "Could not match faces, cell type %ld", celltype);
    }

}
