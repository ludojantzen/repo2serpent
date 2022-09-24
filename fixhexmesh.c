/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fixhexmesh.c                                   */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2018/01/12 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Splits hexahedral meshes to tetrahedrons                     */
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

#define FUNCTION_NAME "FixHexMesh:"

/*****************************************************************************/

void FixHexMesh(long ifc, long nfaces, long ncells)
{
  long cgns, cdone, myk;
  long tetlist, tet, type, fidx;
  long ownrlist, nbrlist, *newnbrlist[4], *newnbrfacelist[4];
  long hexfaces[6];
  long jo, otet, ptr;
  long nf, nc, i, j, k, P, Po;
  long V[8], initFaces[6][4];
  long diags[8][8];
  long *face[4];
  long face1[4] = {0,2,1,3};
  long face2[4] = {1,2,3,0};
  long face3[4] = {0,1,3,2};
  long face4[4] = {0,3,2,1};

  /* Store point list for easy access */

  face[0] = face1;
  face[1] = face2;
  face[2] = face3;
  face[3] = face4;

  /* Get number of parent cells */

  nc = (long)RDB[ifc + IFC_NC_PARENTS];

  /* Get number of parent faces */

  nf = (long)RDB[ifc + IFC_NF_PARENTS];

  /* Get pointer to owner list */

  ownrlist = (long)RDB[ifc + IFC_PTR_OWNR_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(ownrlist)", DATA_ARRAY, ownrlist);

  /* Get pointer to neighbour list */

  nbrlist  = (long)RDB[ifc + IFC_PTR_NBR_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(nbrlist)", DATA_ARRAY, nbrlist);

  /* This list contains the tet cells that have neighbours through surface X */
  /* Since each face is divided at most into two triangles there are at most */
  /* 4 tets that share some parts of surface X (2 on both sides)             */

  for (i = 0; i < 4; i++)
    {
      newnbrlist[i] = (long *)Mem(MEM_ALLOC, nf, sizeof(long));
      newnbrfacelist[i] = (long *)Mem(MEM_ALLOC, nf, sizeof(long));
    }

  /* Get pointer to tet list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Allocate array to store the new tet-cells to */

  ptr = ReallocMem(DATA_ARRAY, ncells);
  WDB[ifc + IFC_PTR_TET_LIST] = (double)ptr;

  /* Reset number of created tets */

  cdone = 0;

  /****** Divide cells *******/

  /* Get pointer to tet list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Loop over cells to prepare them for division */

  for (i = 0; i < nc; i++)
    {

      /* Get pointer to cell */

      cgns = ListPtr(cgns, i);
      CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

      /* Get mesh cell type */

      type = MeshCellType(cgns, ifc);

      /* Reset V */

      for (j = 0; j < 8; j++)
        V[j] = -1;

#ifdef SPLITPRINT
      printf("Cell type %ld\n", type);
#endif

      /* Create the point list V and auxiliary lists */

      MeshCellFromCGNS(ifc, cgns, V, hexfaces, initFaces,
                       type);

      /* Rotate the hexahedron so that the smallest index is bottom left, front */

      MeshCellIndirection(V, type);

#ifdef SPLITPRINT
      PrintMeshCell(V, type);
#endif

      /* Reset diagonals */

      for (j = 0; j < 8; j++)
        for (k = 0; k < 8; k++)
          diags[j][k] = 0;

      /* Get diagonals */

      MeshCellConnectDiags(V, diags, type);

      /* Print connectivity */
#ifdef SPLITPRINT

      for (j = 0; j < 8; j++)
        {
          for (k = 0; k < 8; k++)
            printf("%ld ", diags[j][k]);

          printf("\n");
        }

#endif

      /* Rotate along 0-6 diagonal if needed */

      if (type == MESH_CELL_TYPE_HEX)
        HexRotateCell(V, diags);


#ifdef SPLITPRINT
      PrintMeshCell(V, type);
#endif
      /* Now we don't have to rotate the cell anymore       */
      /* We can map neighbours to the current top, bot etc. */

      MeshCellRotateLists(V, initFaces, hexfaces, type);

#ifdef SPLITPRINT

      printf("Neighbours:  ");

      for (j = 0; j < 6; j++)
        {
          printf("%ld ", hexnbrs[0][j]);
        }

      printf("\n");

      printf("Neighbours2: ");

      for (j = 0; j < 6; j++)
        {
          printf("%ld ", hexnbrs[1][j]);
        }

      printf("\n");

      printf("Faces: ");

      for (j = 0; j < 6; j++)
        {
          printf("%ld ", hexfaces[j]);
        }
      printf("\n");

      printf("Points: ");

      for (j = 0; j < 8; j++)
        {
          printf("%ld ", V[j]);
        }
      printf("\n");
#endif

      /* Create new cells and faces      */

      DivideMeshCell(ifc, cgns, V, diags, hexfaces,
                     newnbrlist, newnbrfacelist, &cdone, type);

    }

  tetlist = (long)RDB[ifc + IFC_PTR_TET_LIST];

  /* Connect out of cell faces */

  for (i = 0; i < cdone; i++)
    {
      /* Get pointer to tet */

      tet = (long)RDB[tetlist + i];

#ifdef SPLITPRINT
      fprintf(outp, "Tet %ld:\n", i);


      fprintf(outp, "Points: ");

      for (j = 0; j < 4; j++)
        fprintf(outp, "%ld ", (long)RDB[tet + TET_POINTS + j]);

      fprintf(outp, "\nNeighbours: ");

      for (j = 0; j < 4; j++)
        fprintf(outp, "%ld ", (long)RDB[tet + TET_NEIGHBOURS + j]);

      fprintf(outp, "\n");
#endif

      /* Loop over faces */

      for (j = 0; j < 4; j++)
        {
          /* Negative values indicate non-linked face */

          if ((fidx = (long)RDB[tet + TET_NEIGHBOURS + j]) >= 0)
            continue;

          /* Get real face idx (we subtracted one in DivideMeshCell */
          /* to keep even face index 0 negative) */

          fidx = -(fidx + 1);

          /* Check that this is one of the owners/neighbours of that face */

          myk = -1;

          for (k = 0; k < 4; k++)
            {
              /* Tet numbers in newnbrlist have been incremented by one to keep them */
              /* nonzero */

              if (newnbrlist[k][fidx] == tet)
                {
                  myk = k;
                }
              else if (newnbrlist[k][fidx] == 0)
                {
                  /* Zero value indicates nothing stored */

                  break;
                }
            }

          /* Reverse k by one since we incremented it to a bad value */

          k--;

          /* Check that this tet was found */

          if (myk < 0)
            Die(FUNCTION_NAME, "Could not find tet %ld from owners/neighbours of face %ld\n", tet, fidx);

          /* Loop over stored tets that share this initial surface */

          for (; k >= 0; k--)
            {
              /* Skip this tet */

              if (k == myk)
                continue;

              /* Get pointer to other tet */

              otet = newnbrlist[k][fidx];

              /* Get face of other tet */

              jo = newnbrfacelist[k][fidx];

              /* Check that my points are the same as points of that */
              /* We don't know the order of points for sure so we'll brute force */

              for (P = 0; P < 3; P++)
                {
                  /* Break the inner loop if the points match */

                  for (Po = 0; Po < 3; Po++)
                    if ((long)RDB[tet  + TET_POINTS + face[j ][P ]] ==
                        (long)RDB[otet + TET_POINTS + face[jo][Po]])
                      break;

                  /* Break if current P was not found */

                  if (Po == 3)
                    break;
                }

              if (P == 3)
                {
                  /* Link these and break if all points were found */

                  WDB[tet + TET_NEIGHBOURS + j] = (double)otet;
                  WDB[otet + TET_NEIGHBOURS + jo] = (double)tet;

                  break;
                }

            }

          /* If k went negative we did not find the other side of this surface */
          /* i.e. the tet should be neighbourless through this surface */

          WDB[tet + TET_NEIGHBOURS + j] = (double)(-1);
        }
    }

  /* This list contains the tet cells that have neighbours through surface X */

  for (i = 0; i < 4; i++)
    {
      Mem(MEM_FREE, newnbrlist[i]);
      Mem(MEM_FREE, newnbrfacelist[i]);
    }

  /* Store number of child cells */

  WDB[ifc + IFC_NC] = (double)cdone;

  return;
}

/*****************************************************************************/
