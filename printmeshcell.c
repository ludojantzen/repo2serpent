/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printmeshcell.c                                */
/*                                                                           */
/* Created:       2015/08/31 (VVa)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Prints a mesh cell as shown below, based on the point list   */
/*                                                                           */
/*                                                                           */
/* Comments:    -For debugging of fixhexmesh.c                               */
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

#define FUNCTION_NAME "MeshCellFromCGNS:"

/*****************************************************************************/

void  PrintMeshCell(long hex[8], long celltype) 
{
  
  switch (celltype)
    {
    case MESH_CELL_TYPE_TET:
      fprintf(outp,"::\n");
      fprintf(outp,"     %2ld \n",hex[3]);
      fprintf(outp,"     /|\\     \n");
      fprintf(outp,"    / | \\   \n");
      fprintf(outp,"   /  |  \\   \n");
      fprintf(outp,"  /   |   \\   \n");
      fprintf(outp,"%2ld ---|-- %2ld\n",hex[0],hex[2]);
      fprintf(outp,"  \\   |   /\n");
      fprintf(outp,"   \\  |  /\n");
      fprintf(outp,"    \\ | /\n");
      fprintf(outp,"     \\|/\n");
      fprintf(outp,"     %2ld \n",hex[1]);

      break;

    case MESH_CELL_TYPE_PYRAMID:
      fprintf(outp,"::\n");
      fprintf(outp,"         %2ld \n",hex[4]);
      fprintf(outp,"         /|\\     \n\n");      
      fprintf(outp,"    %2ld ------ %2ld\n",hex[3],hex[2]);
      fprintf(outp,"    /         /\n");
      fprintf(outp,"   /         / \n");
      fprintf(outp,"  /         /  \n");
      fprintf(outp,"%2ld ------ %2ld\n",hex[0],hex[1]);

      break;

    case MESH_CELL_TYPE_PRISM:
      fprintf(outp,"    %2ld ------ %2ld\n",hex[3],hex[5]);
      fprintf(outp,"     |\\       /|\n");
      fprintf(outp,"     | \\     / |\n");
      fprintf(outp,"     |  \\   /  |\n");
      fprintf(outp,"     |   %2ld    |\n",hex[4]);
      fprintf(outp,"     |    |    |\n");      
      fprintf(outp,"     |    |    |\n");
      fprintf(outp,"    %2ld----|---%2ld\n",hex[0],hex[2]);
      fprintf(outp,"      \\   |   /\n");
      fprintf(outp,"       \\  |  / \n");
      fprintf(outp,"        \\ | /  \n");
      fprintf(outp,"         \\|/   \n");
      fprintf(outp,"          %2ld  \n",hex[1]);

      break;

    case MESH_CELL_TYPE_HEX:
      fprintf(outp,"    %2ld ------ %2ld\n",hex[7],hex[6]);
      fprintf(outp,"    /|        /|\n");
      fprintf(outp,"   / |       / |\n");
      fprintf(outp,"  /  |      /  |\n");
      fprintf(outp,"%2ld ------ %2ld   |\n",hex[4],hex[5]);
      fprintf(outp," |   |     |   |\n");      
      fprintf(outp," |  %2ld ----|- %2ld\n",hex[3],hex[2]);
      fprintf(outp," |  /      |  /\n");
      fprintf(outp," | /       | / \n");
      fprintf(outp," |/        |/  \n");
      fprintf(outp,"%2ld ------ %2ld\n",hex[0],hex[1]);

      break;
    default:

      Die(FUNCTION_NAME, "Invalid cell type %ld", celltype);
    }
}

/*****************************************************************************/
