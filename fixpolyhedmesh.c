/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fixpolyhedmesh.c                               */
/*                                                                           */
/* Created:       2014/02/24 (VVa)                                           */
/* Last modified: 2018/01/12 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Splits polyhedral meshes to tetrahedrons                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FixPolyhedMesh:"

/*****************************************************************************/

void FixPolyhedMesh(long ifc)
{
  long cgns, sdone;
  long ptr, ptr1, surf;
  long sslist, cellpts, facepts, **childlists;
  long nf, nc, np, i, j, k, n, newcells, mychildren;
  long tet, tet1, nbr, p0a, p0b, p1a, p1b, p2a, p2b, prevprint;
  double newfaces;

  /*******************************************************/
  /* Calculate centerpoints of parent cells and their    */
  /* faces. They will form some of the vertices of the   */
  /* child tets so they need to be stored.               */
  /*******************************************************/

  /* Get number of parent cells */

  nc = (long)RDB[ifc + IFC_NC_PARENTS];

  /* Allocate memory for cell centerpoints */

  cellpts = ReallocMem(DATA_ARRAY, nc*3);

  /* Store cell centerpoint list */

  WDB[ifc + IFC_PTR_CELL_CP_LIST_PARENTS] = (double)cellpts;

  /* Get number of parent faces */

  nf = (long)RDB[ifc + IFC_NF_PARENTS];

  /* Allocate memory for face centerpoints */

  facepts = ReallocMem(DATA_ARRAY, nf*3);

  /* Store face centerpoint list */

  WDB[ifc + IFC_PTR_FACE_CP_LIST_PARENTS] = (double)facepts;

  /* Get pointer to surfacelist */

  sslist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(sslist)", DATA_ARRAY, sslist);

  /* Get pointer to tet list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];
  CheckPointer(FUNCTION_NAME, "(cgns00)", DATA_ARRAY, cgns);

  /* Calculate cell and face centerpoints */

  CalculateTetCenter(cgns, sslist, cellpts, facepts);

  /******* Calculate number of new faces and cells **********/

  /* Allocate temporary memory for a child cell list ********/

  childlists = (long **)Mem(MEM_ALLOC, nc, sizeof(long *));

  /* Get pointer to tet list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];
  CheckPointer(FUNCTION_NAME, "(cgns0)", DATA_ARRAY, cgns);

  newfaces = 0;
  newcells = 0;

  /* Loop over cells to calculate number of new faces */

  for (i = 0; i < nc; i++)
    {
      /* Get pointer to cell */

      cgns = ListPtr(cgns, i);
      CheckPointer(FUNCTION_NAME, "(cgns1)", DATA_ARRAY, cgns);

      /* Get number of faces */

      nf = (long)RDB[cgns + IFC_TET_MSH_NF];

      /* Get pointer to face list */

      ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, ptr);

      /* Get pointer to side list */

      ptr1 = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
      CheckPointer(FUNCTION_NAME, "(PTR1)", DATA_ARRAY, ptr1);

      /* Reset number of child cells from current cell */

      mychildren = 0;

      /* Loop over faces */

      for (j = 0; j < nf; j++)
        {
          /* Get index of face */

          n = (long)RDB[ptr + j];

          /* Get pointer to surface */

          surf = (long)RDB[sslist + n];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get number of points on the face */

          np = (long)RDB[surf + UMSH_SURF_N_POINTS];

          /* Add to number of new cells */

          newcells += np;
          mychildren += np;

          /* Add to number of external faces */

          if ((long)RDB[ptr1 + j] == -1)
            {

              /* This cell owns the child faces */
              /* Sides out of the cell */

              newfaces += (double)np;
            }

          /* The sides between the cells on this face*/

          newfaces += (double)np;

          /* The sides between this faces cells and other faces in this cell */

          newfaces += (double)np/2.0;
        }

      /* Allocate memory for a zero-terminated list of child-cells */

      childlists[i] = (long *)Mem(MEM_ALLOC, mychildren+1, sizeof(long));
    }

  /* Check that there are no half faces */

  if(newfaces - (double)((long)newfaces) != 0)
    Die(FUNCTION_NAME, "Something wrong with number of new faces");

  /* Allocate memory for new owners and neighbors list    */
  /* The number of child faces has to be calculated first */

  fprintf(outp, "New TETs should take %f GB (%ld tets)\n", newcells*TET_BLOCK_SIZE*8.0/GIGA, newcells);
  fprintf(outp, "Temporary neighbour lists should take %f GB.\n", 2*(newcells+nc*2)*8.0/GIGA);

  /* Allocate array to store the new tet-cells to */

  ptr = ReallocMem(DATA_ARRAY, newcells);
  WDB[ifc + IFC_PTR_TET_LIST] = (double)ptr;

  /* Loop over parents and divide the cells */
  /* Print */

  sdone = 0;

  fprintf(outp, "   0%% complete\n");
  prevprint = 0;

  for (i = 0; i < nc; i++)
    {
      if ((100*i)/nc > prevprint + 10)
        {
          prevprint += 10;

          fprintf(outp, " %3ld%% complete\n", prevprint);
        }
      DividePolyhedCell(ifc, i, childlists[i], &sdone);
    }

  fprintf(outp, " 100%% complete\n");

  /* Link first faces (out-of-cell faces) */
  fprintf(outp, "\nLinking out-of-cell faces of tetrahedrons.\n");
  fprintf(outp, "   0%% complete\n");
  prevprint = 0;
  for (i = 0; i < nc; i++)
    {
      if ((100*i)/nc > prevprint + 10)
        {
          prevprint += 10;

          fprintf(outp, " %3ld%% complete\n", prevprint);
        }

      /* Loop over children of i'th cell */

      j = 0;
      while ((tet = childlists[i][j]) > VALID_PTR)
        {
          /* Get pre-linked first neighbour */

          nbr = (long)RDB[tet + TET_NEIGHBOURS + 0];

          /* may be parent index (0 -- nc), linked tet (first tet ptr -- last tet ptr) or outside (-1) */

          if ((nbr >= 0) && (nbr < nc))
            {
              /* parent index of neighbour */

              /* Check that children of that neighbour haven't been processed yet */

              if (nbr < i)
                Die(FUNCTION_NAME, "Should link to previously processed tet, but did not link then? Now %ld nbr %ld.", i, nbr);

              /* Check that the neighbour is not the parent of this */

              if (nbr == i)
                Die(FUNCTION_NAME, "Should link to child of same parent?");

              /* Get first three points of this tet (outwards facing face) */

              p0a = (long)RDB[tet + TET_POINTS + 0];
              p1a = (long)RDB[tet + TET_POINTS + 1];
              p2a = (long)RDB[tet + TET_POINTS + 2];

              /* Loop over children of that neighbour */

              k = 0;
              while ((tet1 = childlists[nbr][k]) > VALID_PTR)
                {
                  /* Get third point of this tet (old face centerpoint) and compare */
                  /* to the third point of the other tet */

                  p2b = (long)RDB[tet1 + TET_POINTS + 2];

                  if (p2a == p2b)
                    {
                      /* The points match -> these tets are on the same old face   */
                      /* But might not be opposite to each other so test the other */
                      /* two points                                                */

                      p0b = (long)RDB[tet1 + TET_POINTS + 0];
                      p1b = (long)RDB[tet1 + TET_POINTS + 1];

                      /* Test match (face is reversed) */

                      if ((p1b == p0a) && (p0b == p1a))
                        break;
                    }

                  k++;
                }

              if (tet1 < VALID_PTR)
                {
                  fprintf(outp, "Face to be linked (tet %ld parent %ld):\n", tet, i);
                  fprintf(outp, "%ld %ld %ld\n", p0a, p1a, p2a);

                  fprintf(outp, "\nFaces in neighbour %ld:\n", nbr);

                  k = 0;
                  while ((tet1 = childlists[nbr][k]) > VALID_PTR)
                    {
                      p0b = (long)RDB[tet1 + TET_POINTS + 0];
                      p1b = (long)RDB[tet1 + TET_POINTS + 1];
                      p2b = (long)RDB[tet1 + TET_POINTS + 2];

                      fprintf(outp, "%ld %ld %ld\n", p0b, p1b, p2b);
                      k++;
                    }
                  fprintf(outp, "\n\n");

                Die(FUNCTION_NAME, "Could not find outward facing neighbour from children of %ld\n", nbr);
                }

              /* Link these tets as neighbours */

              WDB[tet  + TET_NEIGHBOURS + 0] = (double)tet1;
              WDB[tet1 + TET_NEIGHBOURS + 0] = (double)tet;
            }

          j++;
        }
    }
  fprintf(outp, " 100%% complete\n");

  /* Free temporary child-list */

  for (i = 0; i < nc; i++)
    Mem(MEM_FREE, childlists[i]);

  Mem(MEM_FREE, childlists);

  /* Put new number of cells  */

  WDB[ifc + IFC_NC] = (double)newcells;

}

/*****************************************************************************/
