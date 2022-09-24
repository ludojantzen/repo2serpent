/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dividepolyhedcell.c                            */
/*                                                                           */
/* Created:       2014/02/24 (VVa)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Splits a polyhedral cell into tetrahedrons by adding points  */
/*              to face centerpoints as well as to the centerpoint of        */
/*              face centerpoints ("cell centerpoint")                       */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DividePolyhedCell:"

/*****************************************************************************/

void DividePolyhedCell(long ifc, long idx, long *mychildren, long *subidx)
{
  long surf, np, i, j, n, cgns;
  long facelist, surflist, sidelist, side;
  long ptr, *newcells, tetlist;
  long ndone, nsides, nc, tet, tet1, p0a, p1a, p0b, p1b, p3a;

  /* Initialize number of created subcells */

  ndone = 0;

  /* Get pointer to interface surface list */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  nc = 0;

  /* Get pointer to parent list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PARENTS];
  CheckPointer(FUNCTION_NAME, "(cgnslist)", DATA_ARRAY, cgns);

  /* Get pointer to cell to be divided */

  cgns = ListPtr(cgns, idx);
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Get number of faces */

  nsides = (long)RDB[cgns + IFC_TET_MSH_NF];

  /* Get pointer to face list */

  ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, ptr);

  /* Loop over faces to get number of sides and points */

  for (j = 0; j < nsides; j++)
    {
      /* Get index of face */

      n = (long)RDB[ptr + j];

      /* Get pointer to face */

      surf = (long)RDB[surflist + n];

      /* Get number of points for face */

      np = (long)RDB[surf + UMSH_SURF_N_POINTS];

      /* Add to number of cells to be created */

      nc += np;
    }

  /*printf("Dividing cell %ld to %ld cells (%ld faces)\n", idx, nc, nsides);*/

  /* Get pointer to face list */

  facelist = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, facelist);

  /* Get pointer to side list */

  sidelist = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
  CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, sidelist);

  /* Get pointer to tet list */

  tetlist = (long)RDB[ifc + IFC_PTR_TET_LIST];
  CheckPointer(FUNCTION_NAME, "(tetlist)", DATA_ARRAY, tetlist);

  /* Loop over faces to divide each face separately */

  for (j = 0; j < nsides; j++)
    {
      /* Get index of face */

      n = (long)RDB[facelist + j];

      /* Get pointer to face */

      surf = (long)RDB[surflist + n];

      /* Get side for this cell */

      side = (long)RDB[sidelist + j];

      /* Get points on the face */

      np = (long)RDB[surf + UMSH_SURF_N_POINTS];

      /* Allocate memory for a list of new face cells */
      /* DividePolyhedFace will use it to return a list of created cgns-cells*/

      newcells = (long *)Mem(MEM_ALLOC, np, sizeof(long));

      /* Divide face into np cells */

      DividePolyhedFace(ifc, idx, n, side, np, newcells);

      /* Copy pointers of new face cells to list of children for this cell */

      for (i = 0; i < np; i++)
        {
          mychildren[ndone + i] = newcells[i];
          WDB[tetlist + *subidx] = (double)newcells[i];
          *subidx = *subidx + 1;
        }

      ndone += np;
      /* Free memory */

      Mem(MEM_FREE,newcells);
    }

  /* Link fourth faces (inside cells but not on same face) */

  for (i = 0; i < nc; i++)
    {
      /* Get tet to link */

      tet = mychildren[i];

      /* Skip if already linked (if the other tet was processed first) */

      if ((long)RDB[tet + TET_NEIGHBOURS + 3] > VALID_PTR)
        continue;

      /* Not already linked, take two first points of this tet */

      p0a = (long)RDB[tet + TET_POINTS + 0];
      p1a = (long)RDB[tet + TET_POINTS + 1];

      /* Take apex also, just to be sure */

      p3a = (long)RDB[tet + TET_POINTS + 3];

      /* Loop over other tets to find matching face */

      for (j = 0; j < nc; j++)
        {
          /* SKip this tet */

          if (j == i)
            continue;

          /* Get other tet */

          tet1 = mychildren[j];

          /* The tets share the two first points (but in reverse order) */
          /* take two first points of other tet */

          p1b = (long)RDB[tet1 + TET_POINTS + 0];
          p0b = (long)RDB[tet1 + TET_POINTS + 1];

          /* Check the two points */

          if ((p0a == p0b) && (p1a == p1b))
            {
              /* These tets should be neighbours */

              /* Test apex, just to be sure */

              if ((long)RDB[tet1 + TET_POINTS + 3] != p3a)
                Die(FUNCTION_NAME, "A major failure in the mesh splitting algorithm.");

              /* Check that other tet does not already have a neighbour */

              if ((long)RDB[tet1 + TET_NEIGHBOURS + 3] > VALID_PTR)
                Die(FUNCTION_NAME, "Other tet already has a neighbour.");

              /* Link these tets as neighbours */

              WDB[tet  + TET_NEIGHBOURS + 3] = (double)tet1;
              WDB[tet1 + TET_NEIGHBOURS + 3] = (double)tet;
            }
        }
    }

  /* Put a NULLPTR to terminate the child cell list */

  mychildren[nc] = NULLPTR;
}

/*****************************************************************************/
