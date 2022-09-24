/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshcelltype.c                                 */
/*                                                                           */
/* Created:       2015/07/27 (VVa)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Finds out type of mesh cell (hexahedral, tet, smth else)     */
/*                                                                           */
/* Comments:                                                                 */
/*             -Based on "J. Dompierre, P. Labbe, M. Vallet and R. Camarero, */
/*              How to Subdivide Pyramids, Prisms, and Hexahedra             */
/*              into Tetrahedra, Proceedings,                                */
/*              8th International Meshing Roundtable,                        */
/*              South Lake Tahoe, CA, U.S.A., pp.195-204, October 1999       */
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

#define FUNCTION_NAME "MeshCellType:"

/*****************************************************************************/

long MeshCellType(long cgns, long ifc)
{
  long nd, ptr, surf, nfacepts[3];
  long surflist;
  long np, n, i;

  /* Get number of faces */

  nd = (long)RDB[cgns + IFC_TET_MSH_NF];

  /* Reset number of different face types */
  /* 3 point faces */
  nfacepts[0] = 0;
  /* 4 point faces */
  nfacepts[1] = 0;
  /* other faces */
  nfacepts[2] = 0;

  /* Get pointer to parents surfacelist */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Get pointer to face list */

  ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];

  /* Loop over faces */

  for (i = 0; i < nd ; i++)
    {
      /* Get surface index */

      n = (long)RDB[ptr + i];

      /* Get pointer to surface */

      surf = (long)RDB[surflist + n];
      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

      /* Get number of points on the face */

      np = (long)RDB[surf + UMSH_SURF_N_POINTS];

      /* Store number of points */

      if (np == 3)
        nfacepts[0]++;
      else if (np == 4)
        nfacepts[1]++;
      else
        nfacepts[2]++;

    }

  /* Check if polyhedral cell */

  /* Less than four or more than six faces */
  if ((nd < 4) || (nd > 6))
    return MESH_CELL_TYPE_POLY;
  /* Contains faces that have less than three or more than */
  /* four points */
  else if (nfacepts[2] > 0)
    return MESH_CELL_TYPE_POLY;
  /* Four sides with three points in each */
  else if ((nd == 4) && (nfacepts[0] == 4))
    return MESH_CELL_TYPE_TET;
  /* Five sides with three points in 4, four points in 1 */
  else if ((nd == 5) && ((nfacepts[0] == 4) && (nfacepts[1] == 1)))
    return MESH_CELL_TYPE_PYRAMID;
  /* Five sides with three points in 2, four points in 3 */
  else if ((nd == 5) && ((nfacepts[0] == 2) && (nfacepts[1] == 3)))
    return MESH_CELL_TYPE_PRISM;
  /* Six sides with four points in each */
  else if ((nd == 6) && (nfacepts[1] == 6))
    return MESH_CELL_TYPE_HEX;
  else
    return MESH_CELL_TYPE_POLY;

}
