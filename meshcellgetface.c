/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshcellgetface.c                              */
/*                                                                           */
/* Created:       2015/08/31 (VVa)                                           */
/* Last modified: 2015/08/31 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Gets face points from a mesh cell based on face direction    */
/*              and cell type                                                */
/*                                                                           */
/* Comments:   -Points for hex face on right are hex[1] hex[2] hex[6] hex[5] */
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

#define FUNCTION_NAME "MeshCellGetFace:"

/*****************************************************************************/

void MeshCellGetFace(long hex[8], long *face, long facedirection, long celltype)
{
  long pos[4], i;

  /* Check direction value */

  CheckValue(FUNCTION_NAME, "face direction", "", facedirection, 0, 5);

  /* Get face positions */

  MeshCellGetFacePos(pos, facedirection, celltype);

  /* Store corresponding points from hex to face */

  for (i = 0; i < 4; i++)
    face[i] = hex[pos[i]];

}
