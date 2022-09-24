/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshcellgetfacepos.c                           */
/*                                                                           */
/* Created:       2015/08/31 (VVa)                                           */
/* Last modified: 2015/08/31 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Gets face positions of a cell based on face direction        */
/*                                                                           */
/* Comments:   -Positions for hex face on right are 1 2 6 5 etc.             */
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

#define FUNCTION_NAME "MeshCellGetFacePos:"

/*****************************************************************************/

void MeshCellGetFacePos(long *pos, long facedirection, long celltype)
{
  switch (celltype)
    {
    case (MESH_CELL_TYPE_TET):

      switch (facedirection)
        {

        case (MESH_CELL_FACE_BOTTOM):

          pos[0] = 0;
          pos[1] = 1;
          pos[2] = 2;
          pos[3] = -1;

          break;

        case (MESH_CELL_FACE_LEFT):

          pos[0] = 0;
          pos[1] = 3;
          pos[2] = 1;
          pos[3] = -1;

          break;

        case (MESH_CELL_FACE_RIGHT):

          pos[0] = 1;
          pos[1] = 2;
          pos[2] = 3;
          pos[3] = -1;

          break;

        case (MESH_CELL_FACE_BACK):

          pos[0] = 0;
          pos[1] = 3;
          pos[2] = 2;
          pos[3] = -1;

          break;

        default:
          Die(FUNCTION_NAME, "Invalid face index for tet %ld\n", facedirection);      

        }

      break;

    case MESH_CELL_TYPE_PYRAMID:

      switch (facedirection)
        {

        case (MESH_CELL_FACE_BOTTOM):

          pos[0] = 0;
          pos[1] = 1;
          pos[2] = 2;
          pos[3] = 3;

          break;

        case (MESH_CELL_FACE_LEFT):

          pos[0] = 0;
          pos[1] = 3;
          pos[2] = 4;
          pos[3] = -1;

          break;

        case (MESH_CELL_FACE_RIGHT):

          pos[0] = 1;
          pos[1] = 2;
          pos[2] = 4;
          pos[3] = -1;

          break;

        case (MESH_CELL_FACE_FRONT):

          pos[0] = 0;
          pos[1] = 4;
          pos[2] = 1;
          pos[3] = -1;

          break;

        case (MESH_CELL_FACE_BACK):

          pos[0] = 2;
          pos[1] = 3;
          pos[2] = 4;
          pos[3] = -1;

          break;

        default:
          Die(FUNCTION_NAME, "Invalid face index for pyram %ld\n", facedirection);      

        }

      break;

    case MESH_CELL_TYPE_PRISM:

      switch (facedirection)
        {

        case (MESH_CELL_FACE_BOTTOM):

          pos[0] = 0;
          pos[1] = 1;
          pos[2] = 2;
          pos[3] = -1;

          break;

        case (MESH_CELL_FACE_TOP):

          pos[0] = 3;
          pos[1] = 4;
          pos[2] = 5;
          pos[3] = -1;

          break;

        case (MESH_CELL_FACE_LEFT):

          pos[0] = 0;
          pos[1] = 3;
          pos[2] = 4;
          pos[3] = 1;

          break;

        case (MESH_CELL_FACE_RIGHT):

          pos[0] = 1;
          pos[1] = 2;
          pos[2] = 5;
          pos[3] = 4;

          break;

        case (MESH_CELL_FACE_BACK):

          pos[0] = 0;
          pos[1] = 3;
          pos[2] = 5;
          pos[3] = 2;

          break;

        default:
          Die(FUNCTION_NAME, "Invalid face index for prism %ld\n", facedirection);      

        }

      break;

    case MESH_CELL_TYPE_HEX:

      switch (facedirection)
        {

        case (MESH_CELL_FACE_BOTTOM):

          pos[0] = 0;
          pos[1] = 1;
          pos[2] = 2;
          pos[3] = 3;

          break;

        case (MESH_CELL_FACE_TOP):

          pos[0] = 4;
          pos[1] = 5;
          pos[2] = 6;
          pos[3] = 7;

          break;

        case (MESH_CELL_FACE_LEFT):

          pos[0] = 0;
          pos[1] = 3;
          pos[2] = 7;
          pos[3] = 4;

          break;

        case (MESH_CELL_FACE_RIGHT):

          pos[0] = 1;
          pos[1] = 2;
          pos[2] = 6;
          pos[3] = 5;

          break;

        case (MESH_CELL_FACE_FRONT):

          pos[0] = 0;
          pos[1] = 4;
          pos[2] = 5;
          pos[3] = 1;

          break;

        case (MESH_CELL_FACE_BACK):

          pos[0] = 2;
          pos[1] = 3;
          pos[2] = 7;
          pos[3] = 6;

          break;

        default:
          Die(FUNCTION_NAME, "Invalid face index %ld\n", facedirection);      

        }

      break;

    default:
      Die(FUNCTION_NAME, "Invalid cell type %ld\n", celltype);
      
    }
}
