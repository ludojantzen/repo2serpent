/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : createumshcells.c                              */
/*                                                                           */
/* Created:       2017/08/03 (VVa)                                           */
/* Last modified: 2017/11/29 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Creates UMSH geometry/interface cells based on points, faces */
/*              and owner/neighbour information.                             */
/*                                                                           */
/* Comments:   -                                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CreateUMSHCells:"

/*****************************************************************************/

void CreateUMSHCells(long ifc)
{
  long n, fi, solist, snlist, i, nf, nc, np, nd, sslist, **surfaces, **sides, *nsurf;
  long cgns, nbr, ptr, ptr1, surf, j;
  double x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;

  /* Get number of cells */

  nc = (long)RDB[ifc + IFC_NC];

  /* Get number of faces */

  nf = (long)RDB[ifc + IFC_NF_PARENTS];

  /* Check number of cells */

  CheckValue(FUNCTION_NAME, "nc", "", nc, 1, 100000000000);

  /* Get owner and neighbor lists */

  solist = (long)RDB[ifc + IFC_PTR_OWNR_LIST_PARENTS];
  snlist = (long)RDB[ifc + IFC_PTR_NBR_LIST_PARENTS];

  /* Preallocate memory for data */

  PreallocMem(nc*(IFC_TET_MSH_LIST_BLOCK_SIZE), DATA_ARRAY);

  /* Avoid compiler warning */

  cgns = -1;

  /* Create cells */

  for (n = 0; n < nc; n++)
    {
      /* Allocate memory for tet cell */

      cgns = NewItem(ifc + IFC_PTR_TET_MSH, IFC_TET_MSH_LIST_BLOCK_SIZE);

      /* Put index */

      WDB[cgns + IFC_TET_PRNT_IDX] = (double)n;

      /* Reset cell boundaries */

      WDB[cgns + IFC_TET_PRNT_XMIN] = INFTY;
      WDB[cgns + IFC_TET_PRNT_XMAX] = -INFTY;
      WDB[cgns + IFC_TET_PRNT_YMIN] = INFTY;
      WDB[cgns + IFC_TET_PRNT_YMAX] = -INFTY;
      WDB[cgns + IFC_TET_PRNT_ZMIN] = INFTY;
      WDB[cgns + IFC_TET_PRNT_ZMAX] = -INFTY;

    }

  /* Close list */

  CloseList(cgns);

  /* Copy list pointer to parent list */

  WDB[ifc + IFC_PTR_TET_MSH_PARENTS] = RDB[ifc + IFC_PTR_TET_MSH];

  /***********************************************************************/

  /***** Put geometry data to cells **************************************/

  /* Get pointer to surface list */

  sslist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PARENTS];

  /* Check number of faces */

  CheckValue(FUNCTION_NAME, "nf", "", nf, 4, 100000000000);

  /* Put number of faces */

  WDB[ifc + IFC_NF_PARENTS] = (double)nf;

  /* Put number of cells */

  WDB[ifc + IFC_NC] = (double)nc;
  WDB[ifc + IFC_NC_PARENTS] = (double)nc;

  /* Allocate memory for face lists for cells */

  /* List of cells */

  surfaces = (long **)Mem(MEM_ALLOC, nc, sizeof(long*));
  sides = (long **)Mem(MEM_ALLOC, nc, sizeof(long*));

  /* Allocate memory for four faces initially */

  for (n = 0; n < nc; n++)
    {
      surfaces[n] = (long *)Mem(MEM_ALLOC, 4, sizeof(long));
      sides[n] = (long *)Mem(MEM_ALLOC, 4, sizeof(long));
      /* Reset faces */

      for (i = 0; i < 4; i++)
        {
          surfaces[n][i] = -1;
          sides[n][i] = -1;
        }
    }

  /* Number of surfaces per cell */

  nsurf = (long *)Mem(MEM_ALLOC, nc, sizeof(long*));

  /* Loop over faces */

  for (n = 0; n < nf; n++)
    {

      /* Repeat for owner and neighbour files */

      for (fi = 1; fi < 3; fi++)
        {
          /* Get cell index */
          if(fi == 1)
            i = (long)RDB[solist + n];
          else
            i = (long)RDB[snlist + n];

          /* Check */

          if (i < 0)
            continue;
          else if (i > nc - 1)
            Error(ifc, "Cell index exceeds maximum");

          /* Get pointer to tet cell */

          cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];
          CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

          cgns = ListPtr(cgns, i);
          CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

          /* Store straight pointer of cell to owner/neighbor list */

          if(fi == 1)
            WDB[solist + n] = (double)cgns;
          else
            WDB[snlist + n] = (double)cgns;

          /* Check index */

          if ((long)RDB[cgns + IFC_TET_PRNT_IDX] != i)
            Die(FUNCTION_NAME, "Indexing error");

          /* Get neighbour index */

          if(fi == 1)
            {
              /* Neighbour list contains index (not pointer) */

              i = (long)RDB[snlist + n];

              if (i > -1)
                {
                  if (i > nc - 1)
                    Error(ifc, "Cell index exceeds maximum");

                  nbr = ListPtr(cgns, i);
                  CheckPointer(FUNCTION_NAME, "(nbr)", DATA_ARRAY, nbr);

                  /* Check index */

                  if ((long)RDB[nbr + IFC_TET_PRNT_IDX] != i)
                    Die(FUNCTION_NAME, "Indexing error (%ld %ld)",
                        (long)RDB[nbr + IFC_TET_PRNT_IDX], i);

                }
              else
                nbr = -1;
            }
          else
            {
              /* Owner list already contains direct pointer */

              nbr = (long)RDB[solist + n];
            }

          /* Get owner / neighbour */

          if(fi == 1)
            cgns = (long)RDB[solist + n];
          else
            cgns = (long)RDB[snlist + n];

          /* Get index of owner */

          i = (long)RDB[cgns + IFC_TET_PRNT_IDX];

          /* Store index of face to cell's face list */
          /* Check that cells surface list is not full */

          if (nsurf[i] >= 4)
            {

              /* Allocate memory for one more surface */

              surfaces[i] = (long*)Mem(MEM_REALLOC, surfaces[i], (nsurf[i]+1)*sizeof(long));
              sides[i] = (long*)Mem(MEM_REALLOC, sides[i], (nsurf[i]+1)*sizeof(long));

              /* Store face index */

              surfaces[i][nsurf[i]] = n;

              if(fi == 1)
                sides[i][nsurf[i]] = -1;
              else
                sides[i][nsurf[i]] = 1;

              /* Increase number of stored surfaces */

              nsurf[i]++;

            }
          else
            {

              /* Store face pointer */

              surfaces[i][nsurf[i]] = n;

              if(fi == 1)
                sides[i][nsurf[i]] = -1;
              else
                sides[i][nsurf[i]] = 1;

              /* Increase number of stored surfaces */

              nsurf[i]++;

            }

          /* Get cell bounding box */

          xmin = RDB[cgns + IFC_TET_PRNT_XMIN];
          xmax = RDB[cgns + IFC_TET_PRNT_XMAX];
          ymin = RDB[cgns + IFC_TET_PRNT_YMIN];
          ymax = RDB[cgns + IFC_TET_PRNT_YMAX];
          zmin = RDB[cgns + IFC_TET_PRNT_ZMIN];
          zmax = RDB[cgns + IFC_TET_PRNT_ZMAX];

          /* Get pointer to surface */

          surf = (long)RDB[sslist + n];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get pointer to points */

          ptr = (long)RDB[surf + UMSH_SURF_PTR_POINTS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get number of points */

          np = (long)RDB[surf + UMSH_SURF_N_POINTS];
          CheckValue(FUNCTION_NAME, "np", "", np, 3, 900);

          /* Loop over values */

          for (j = 0; j < np; j++)
            {

              /* Get point */
              ptr1 = (long)RDB[ptr + j];

              /* Get coordinates */

              x = RDB[ptr1++];
              y = RDB[ptr1++];
              z = RDB[ptr1++];

              /* Compare to bounding box */

              if (x < xmin)
                xmin = x;
              if (x > xmax)
                xmax = x;

              if (y < ymin)
                ymin = y;
              if (y > ymax)
                ymax = y;

              if (z < zmin)
                zmin = z;
              if (z > zmax)
                zmax = z;
            }

          /* Put cell boundaries */

          WDB[cgns + IFC_TET_PRNT_XMIN] = xmin;
          WDB[cgns + IFC_TET_PRNT_XMAX] = xmax;
          WDB[cgns + IFC_TET_PRNT_YMIN] = ymin;
          WDB[cgns + IFC_TET_PRNT_YMAX] = ymax;
          WDB[cgns + IFC_TET_PRNT_ZMIN] = zmin;
          WDB[cgns + IFC_TET_PRNT_ZMAX] = zmax;
        }
    }

  /* Loop over cell list */

  for (n = 0; n < nc; n++)
    {
      /* Get pointer to tet cell */

      cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];
      CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

      cgns = ListPtr(cgns, n);
      CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

      /* Get number of faces */

      nd = nsurf[n];

      /* Put number of faces */

      WDB[cgns + IFC_TET_MSH_NF] = (double)nd;

      /* Allocate memory for face list */

      ptr = ReallocMem(DATA_ARRAY, nd);

      /* Put pointer to face list */

      WDB[cgns + IFC_TET_MSH_PTR_FACES] = (double)ptr;

      /* Allocate memory for side list */

      ptr1 = ReallocMem(DATA_ARRAY, nd);

      /* Put pointer to side list */

      WDB[cgns + IFC_TET_MSH_PTR_SIDES] = (double)ptr1;

      /* Loop over faces to store them */

      for (i = 0; i < nd ; i++)
        {
          WDB[ptr + i] = (double)surfaces[n][i];
          WDB[ptr1 + i] = (double)sides[n][i];
        }

    }

  /* Free temporary lists */

  for (n = 0; n < nc; n++)
    Mem(MEM_FREE, surfaces[n]);

  Mem(MEM_FREE, surfaces);

  Mem(MEM_FREE, nsurf);
}
