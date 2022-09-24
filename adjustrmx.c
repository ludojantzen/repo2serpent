/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : vrcycle.c                                      */
/*                                                                           */
/* Created:       2018/10/09 (JLe)                                           */
/* Last modified: 2019/04/11 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Adjusts RMX mesh in self-adaptive mode.                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

long TestNeighbour(double *lims0, double *lims1);

void TestRMX(long rmx);

#define FUNCTION_NAME "AdjustRMX:"

/*****************************************************************************/

long AdjustRMX(long rmx, long crit)
{
  long msh0, msh1, itp, loc0, loc1, loc2, loc3, ptr, ng, n, idx, i, j, k, l;
  long nx, ny, nz, i0, j0, k0, nmax1, nmax2, det, nsplit, splitmax;
  double lims0[6], lims1[6], f, *tot, *flx, dmin, xmin, ymin, zmin;

  /* Allow memory allocation */

  Mem(MEM_ALLOW);

  /* Re-open cell list */

  ptr = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);
  ReopenList(ptr);

  /* Maximum number of splits */

  splitmax = (long)RDB[DATA_MAX_RMX_SPLIT_FLAGS];
  CheckValue(FUNCTION_NAME, "splitmax", "", splitmax, 1, 100);

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /* Get pointer to iteration data */

  itp = (long)RDB[rmx + RMX_PTR_ITER];
  CheckPointer(FUNCTION_NAME, "(itp)", DATA_ARRAY, itp);

  /***************************************************************************/

  /***** Select mesh cells for splitting *************************************/

  /* Loop over mesh and reset split flags */

  msh0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh0 > VALID_PTR)
    {
      /* Reset common split flag if not density mode */

      if (crit > 0)
        WDB[msh0 + RMX_CELL_SPLIT_FLAG] = (double)NO;

      /* Next */

      msh0 = NextItem(msh0);
    }

  /* Put splits by contribution */

  if (crit == 4)
    ContribSplit(rmx);

  /* Loop over mesh */

  msh0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (msh0 > VALID_PTR)
    {
      /* Get pointer to mesh */

      ptr = (long)RDB[msh0 + RMX_CELL_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);

      /* Get index */

      i = (long)RDB[msh0 + RMX_CELL_MESH_IDX];
      CheckValue(FUNCTION_NAME, "i", "", i, 0, 1000000000);

      /* Get boundaries */

      MeshCellBounds(ptr, i, &lims0[0], &lims0[1], &lims0[2], &lims0[3],
                     &lims0[4], &lims0[5]);

      /* Calculate minimum dimension */

      dmin = INFTY;

      if (lims0[1] - lims0[0] < dmin)
        dmin = lims0[1] - lims0[0];

      if (lims0[3] - lims0[2] < dmin)
        dmin = lims0[3] - lims0[2];

      if (lims0[5] - lims0[4] < dmin)
        dmin = lims0[5] - lims0[4];

      /* Check criterion */

      if (crit == 1)
        {
          /*******************************************************************/

          /***** Split by mfp ************************************************/

          /* NOTE: Tässä on vain yksi energiaryhmä */

          /* Get split criterion */

          f = RDB[itp + WWD_ITER_SPLIT_MFP_CRIT];
          CheckValue(FUNCTION_NAME, "f", "", f, 1.0, 1E18);

          /* Get total reaction rate and flux for thick materials */

          ptr = (long)RDB[msh0 + RMX_CELL_MC_MFP_FLX0];
          CheckPointer(FUNCTION_NAME, "(ptr3)", RES2_ARRAY, ptr);
          flx = &RES2[ptr];

          ptr = (long)RDB[msh0 + RMX_CELL_MC_MFP_TOTRR0];
          CheckPointer(FUNCTION_NAME, "(ptr4)", RES2_ARRAY, ptr);
          tot = &RES2[ptr];

          /* Compare mfp and dimension */

          if ((tot[0] > 0.0) && (flx[0] > 0.0))
            if (flx[0]/tot[0] < f*dmin)
              {
                /* Put split flag (use common) */

                WDB[msh0 + RMX_CELL_SPLIT_FLAG] = (double)YES;
              }

          /* Reset data */

          tot[0] = 0.0;
          flx[0] = 0.0;

          /* Get total reaction rate and flux for thin materials */

          ptr = (long)RDB[msh0 + RMX_CELL_MC_MFP_FLX1];
          CheckPointer(FUNCTION_NAME, "(ptr5)", RES2_ARRAY, ptr);
          flx = &RES2[ptr];

          ptr = (long)RDB[msh0 + RMX_CELL_MC_MFP_TOTRR1];
          CheckPointer(FUNCTION_NAME, "(ptr6)", RES2_ARRAY, ptr);
          tot = &RES2[ptr];

          /* Compare mfp and dimension */

          if ((tot[0] > 0.0) && (flx[0] > 0.0))
            if (flx[0]/tot[0] < f*dmin)
              {
                /* Put split flag (use common) */

                WDB[msh0 + RMX_CELL_SPLIT_FLAG] = (double)YES;
              }

          /* Reset data */

          tot[0] = 0.0;
          flx[0] = 0.0;

          /*******************************************************************/
        }
      else if (crit == 2)
        {
          /*******************************************************************/

          /***** Split by importance *****************************************/

          /* Get split criterion */

          f = RDB[itp + WWD_ITER_SPLIT_IMP_CRIT];
          CheckValue(FUNCTION_NAME, "f", "", f, 1.0, 1E18);

          /* Pointer to importances */

          loc0 = (long)RDB[msh0 + RMX_CELL_IMP_CURR];
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Loop over neighbours */

          loc2 = (long)RDB[msh0 + RMX_CELL_PTR_BOUNDS];
          while (loc2 > VALID_PTR)
            {
              /* Pointer to mesh cell */

              msh1 = (long)RDB[loc2 + RMX_CELL_BOUND_PTR_CELL];
              CheckPointer(FUNCTION_NAME, "msh1", DATA_ARRAY, msh1);

              /* Pointer to importances */

              loc1 = (long)RDB[msh1 + RMX_CELL_IMP_CURR];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

              /* Loop over energy groups and set split flag */

              for (n = 0; n < ng; n++)
                if ((RDB[loc0 + n] > 0.0) && (RDB[loc1 + n] > 0.0))
                  if (((RDB[loc0 + n] > RDB[loc1 + n]) &&
                       (RDB[loc0 + n]/RDB[loc1 + n] > f)) ||
                      ((RDB[loc1 + n] > RDB[loc0 + n]) &&
                       (RDB[loc1 + n]/RDB[loc0 + n] > f)))
                    {
                      /* Put split flag (use common) */

                      WDB[msh0 + RMX_CELL_SPLIT_FLAG] = (double)YES;
                    }

              /* Next neighbour */

              loc2 = NextItem(loc2);
            }

          /*******************************************************************/
        }
      else if (crit == 3)
        {
          /*******************************************************************/

          /***** Split by number of neighbours *******************************/

          /* Split if too many neighbours */

          if ((long)RDB[msh0 + RMX_CELL_MTX_SIZE] >
              (long)RDB[itp + WWD_ITER_SPLIT_MAX_NEIGH])
            {
              /* Put split flag (use common) */

              WDB[msh0 + RMX_CELL_SPLIT_FLAG] = (double)YES;
            }

          /*******************************************************************/
        }

      /* Next cell */

      msh0 = NextItem(msh0);
    }

  /***************************************************************************/

  /***** Handle splits *******************************************************/

  /* Reset count */

  nsplit = 0;

  /* Pointer to mesh data */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get next index */

  idx = 0;

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Compare index */

      if (idx < (long)RDB[loc0 + RMX_CELL_IDX])
        idx = (long)RDB[loc0 + RMX_CELL_IDX];

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Loop over data */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Reset processed flag (this is used to prevent multiple */
      /* memory allocations later on) */

      WDB[loc0 + RMX_CELL_MTX_ADA_PROC] = (double)NO;

      /* Get split size */

      nx = (long)RDB[itp + WWD_ITER_SPLIT_NX];
      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 100);

      ny = (long)RDB[itp + WWD_ITER_SPLIT_NY];
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 100);

      nz = (long)RDB[itp + WWD_ITER_SPLIT_NZ];
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 100);

      /* Get pointer to mesh */

      ptr = (long)RDB[loc0 + RMX_CELL_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(ptr7)", DATA_ARRAY, ptr);

      /* Get index */

      i = (long)RDB[loc0 + RMX_CELL_MESH_IDX];
      CheckValue(FUNCTION_NAME, "i", "", i, 0, 1000000000);

      /* Get boundaries */

      MeshCellBounds(ptr, i, &lims0[0], &lims0[1], &lims0[2], &lims0[3],
                     &lims0[4], &lims0[5]);

      /***********************************************************************/

      /***** Split by density and check dimensions ***************************/

      /* NOTE: noi cellikohtaiset vektorit on luotu probevrmesh.c:ssä */
      /* tässä loopataan niiden yli ja etsitään ensimmäinen matchi,   */
      /* jota verrataan solun dimensioihin. Muilla crit-arvoilla se   */
      /* flagi asetetaan jo tuolla edellisessä silmukassa. */

      if (crit == 0)
        {
          /* Pointer to split density vector */

          loc1 = (long)RDB[itp + WWD_ITER_SPLIT_PTR_DENS];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Loop over split flag vector */

          ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over flags and find first matching criterion */

          for (n = 0; n < splitmax; n++)
            {
              /* Check flag */

              if ((long)RDB[ptr] == YES)
                break;

              /* Update pointer */

              ptr++;

              /* Next criteria */

              if ((loc1 = NextItem(loc1)) < VALID_PTR)
                break;
            }

          /* Check pointer */

          if (loc1 > VALID_PTR)
            {
              /* Compare to minimum dimensions */

              if ((lims0[1] - lims0[0])/((double)nx)
                  < RDB[loc1 + WWD_ITER_SPLIT_DENS_XMIN])
                nx = 1;
              if ((lims0[3] - lims0[2])/((double)ny)
                  < RDB[loc1 + WWD_ITER_SPLIT_DENS_YMIN])
                ny = 1;
              if ((lims0[5] - lims0[4])/((double)nz)
                  < RDB[loc1 + WWD_ITER_SPLIT_DENS_ZMIN])
                nz = 1;

              /* Set flag */

              if (nx*ny*nz > 1)
                WDB[loc0 + RMX_CELL_SPLIT_FLAG] = (double)YES;
            }
        }
      else if ((long)RDB[loc0 + RMX_CELL_SPLIT_FLAG] == YES)
        {
          /* Get minimum dimensions for other criteria */

          xmin = INFTY;
          ymin = INFTY;
          zmin = INFTY;

          /* Loop over vector */

          loc1 = (long)RDB[itp + WWD_ITER_SPLIT_PTR_DENS];
          while (loc1 > VALID_PTR)
            {
              /* Compare */

              if (RDB[loc1 + WWD_ITER_SPLIT_DENS_XMIN] < xmin)
                xmin = RDB[loc1 + WWD_ITER_SPLIT_DENS_XMIN];

              if (RDB[loc1 + WWD_ITER_SPLIT_DENS_YMIN] < ymin)
                ymin = RDB[loc1 + WWD_ITER_SPLIT_DENS_YMIN];

              if (RDB[loc1 + WWD_ITER_SPLIT_DENS_ZMIN] < zmin)
                zmin = RDB[loc1 + WWD_ITER_SPLIT_DENS_ZMIN];

              /* Next */

              loc1 = NextItem(loc1);
            }

          /* Compare to minimum dimensions */

          if ((lims0[1] - lims0[0])/((double)nx) < xmin)
            nx = 1;
          if ((lims0[3] - lims0[2])/((double)ny) < ymin)
            ny = 1;
          if ((lims0[5] - lims0[4])/((double)nz) < zmin)
            nz = 1;

          /* Set flag */

          if (nx*ny*nz > 1)
            WDB[loc0 + RMX_CELL_SPLIT_FLAG] = (double)YES;
        }

      /* Check split flag */

      if ((long)RDB[loc0 + RMX_CELL_SPLIT_FLAG] == NO)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /***********************************************************************/

      /***** Create a new mesh cell ******************************************/

      /* Get pointer to mesh */

      ptr = (long)RDB[loc0 + RMX_CELL_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(ptr8)", DATA_ARRAY, ptr);

      /* Get index */

      i = (long)RDB[loc0 + RMX_CELL_MESH_IDX];
      CheckValue(FUNCTION_NAME, "i", "", i, 0, 1000000000);

      /* Get boundaries */

      MeshCellBounds(ptr, i, &lims0[0], &lims0[1], &lims0[2], &lims0[3],
                     &lims0[4], &lims0[5]);

      /* Create new mesh */

      msh0 = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_PTR, -1, nx, ny, nz,
                        lims0, -1);
      /* Add to count */

      nsplit++;

      /* Change type to adaptive */

      WDB[msh0 + MESH_TYPE] = (double)MESH_TYPE_ADAPTIVE;

      /* Get pointer to pointer */

      ptr = (long)RDB[ptr + MESH_PTR_PTR] + i;
      CheckPointer(FUNCTION_NAME, "(ptr9)", DATA_ARRAY, ptr);

      /* Put pointer */

      WDB[ptr] = -msh0;

      /* Reset index */

      i = 0;

      /* Loop over mesh */

      for (k0 = 0; k0 < nz; k0++)
      for (j0 = 0; j0 < ny; j0++)
      for (i0 = 0; i0 < nx; i0++)
        {
          /* Create structure */

          loc1 = NewItem(rmx + RMX_PTR_MESH_DATA, RMX_CELL_BLOCK_SIZE);

          /* Put pointer */

          ptr = ReadMeshPtr(msh0, i0, j0, k0);
          CheckPointer(FUNCTION_NAME, "(ptr10)", DATA_ARRAY, ptr);
          WDB[ptr] = (double)loc1;

          /* Put index and pointer */

          WDB[loc1 + RMX_CELL_IDX] = (double)(++idx);
          WDB[loc1 + RMX_CELL_MESH_IDX] = (double)(i++);
          WDB[loc1 + RMX_CELL_PTR_MESH] = (double)msh0;

          /* Calculate volume */

          WDB[loc1 + RMX_CELL_RVOL] = RDB[loc0 + RMX_CELL_RVOL]
            /((double)(nx*ny*nz));

          /* Copy source importances */

          loc2 = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          ptr = ReallocMem(DATA_ARRAY, ng);
          for (n = 0; n < ng; n++)
            WDB[ptr + n] = RDB[loc2 + n];

          WDB[loc1 + RMX_CELL_IMP_SRC] = (double)ptr;

          loc2 = (long)RDB[loc0 + RMX_CELL_IMP_SRC_KEEP];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          ptr = ReallocMem(DATA_ARRAY, ng);
          for (n = 0; n < ng; n++)
            WDB[ptr + n] = RDB[loc2 + n];

          WDB[loc1 + RMX_CELL_IMP_SRC_KEEP] = (double)ptr;

          /* Copy current importances */

          loc2 = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          ptr = ReallocMem(DATA_ARRAY, ng);
          for (n = 0; n < ng; n++)
            WDB[ptr + n] = RDB[loc2 + n];

          WDB[loc1 + RMX_CELL_IMP_CURR] = (double)ptr;

          loc2 = (long)RDB[loc0 + RMX_CELL_IMP_CURR_KEEP];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          ptr = ReallocMem(DATA_ARRAY, ng);
          for (n = 0; n < ng; n++)
            WDB[ptr + n] = RDB[loc2 + n];

          WDB[loc1 + RMX_CELL_IMP_CURR_KEEP] = (double)ptr;

          /* Loop over detectors */

          det = (long)RDB[rmx + RMX_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Create structure */

              ptr = NewItem(loc1 + RMX_CELL_PTR_DET, RMX_DET_BLOCK_SIZE);

              /* Copy data */

              for (j = LIST_DATA_SIZE; j < RMX_DET_BLOCK_SIZE - 1; j++)
                WDB[ptr + j] = RDB[det + j];

              /* Put pointer to parent */

              WDB[ptr + RMX_DET_PTR_DET0] = (double)det;

              /* Next detector */

              det = NextItem(det);
            }
        }

      /***********************************************************************/

      /***** Find neighbours *************************************************/

      /* Loop over mesh */

      for (k0 = 0; k0 < nz; k0++)
      for (j0 = 0; j0 < ny; j0++)
      for (i0 = 0; i0 < nx; i0++)
        {
          /* Pointer to data */

          loc1 = ReadMeshPtr(msh0, i0, j0, k0);
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          loc1 = (long)RDB[loc1];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Mesh index */

          l = (long)RDB[loc1 + RMX_CELL_MESH_IDX];

          /* Get boundaries */

          MeshCellBounds(msh0, l, &lims0[0], &lims0[1], &lims0[2], &lims0[3],
                         &lims0[4], &lims0[5]);

          /* Reset maximum number of neighbours */

          nmax1 = 0;

          /* Add new cells as neighbours for the created cell */

          for (i = i0 - 1; i < i0 + 2; i++)
          for (j = j0 - 1; j < j0 + 2; j++)
          for (k = k0 - 1; k < k0 + 2; k++)
            if ((i != i0) || (j != j0) || (k != k0))
              {
                /* Check boundaries */

                if (((i > -1) && (i < nx) && (j > -1) && (j < ny) &&
                     (k > -1) && (k < nz)) &&
                    (((i != i0) && (j == j0) && (k == k0)) ||
                     ((i == i0) && (j != j0) && (k == k0)) ||
                     ((i == i0) && (j == j0) && (k != k0))))
                  {
                    /* Pointer to data */

                    loc2 = ReadMeshPtr(msh0, i, j, k);
                    CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

                    loc2 = (long)RDB[loc2];
                    CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

                    /* Add to bounds vector */

                    ptr = NewItem(loc1 + RMX_CELL_PTR_BOUNDS,
                                  RMX_CELL_BOUND_BLOCK_SIZE);

                    /* Put pointer and reset indexes */

                    WDB[ptr + RMX_CELL_BOUND_PTR_CELL] = (double)loc2;
                    WDB[ptr + RMX_CELL_BOUND_FWD_IDX] = -1.0;
                    WDB[ptr + RMX_CELL_BOUND_ADJ_IDX] = -1.0;

                    /* Add to maximum number of neighbours */

                    nmax1++;
                  }
              }

          /* Add old cells as neighbours for the created cell */

          loc2 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
          while (loc2 > VALID_PTR)
            {
              /* Pointer to cell */

              ptr = (long)RDB[loc2 + RMX_CELL_BOUND_PTR_CELL];
              CheckPointer(FUNCTION_NAME, "(ptr11)", DATA_ARRAY, ptr);

              /* Get pointer to mesh */

              msh1 = (long)RDB[ptr + RMX_CELL_PTR_MESH];
              CheckPointer(FUNCTION_NAME, "(msh1)", DATA_ARRAY, msh1);

              /* Mesh index */

              l = (long)RDB[ptr + RMX_CELL_MESH_IDX];

              /* Get boundaries */

              MeshCellBounds(msh1, l, &lims1[0], &lims1[1], &lims1[2],
                             &lims1[3], &lims1[4], &lims1[5]);

              /* Check joining surface */

              if (TestNeighbour(lims0, lims1))
                {
                  /* Create new neighour */

                  ptr = NewItem(loc1 + RMX_CELL_PTR_BOUNDS,
                                RMX_CELL_BOUND_BLOCK_SIZE);

                  /* Copy pointer and reset indexes */

                  WDB[ptr + RMX_CELL_BOUND_PTR_CELL] =
                    RDB[loc2 + RMX_CELL_BOUND_PTR_CELL];
                  WDB[ptr + RMX_CELL_BOUND_FWD_IDX] = -1.0;
                  WDB[ptr + RMX_CELL_BOUND_ADJ_IDX] = -1.0;

                  /* Add to maximum number of neighbours */

                  nmax1++;
                }

              /* Next */

              loc2 = NextItem(loc2);
            }

          /* Put matrix size */

          if (nmax1 == 0)
            Die(FUNCTION_NAME, "No neighbours");
          else
            WDB[loc1 + RMX_CELL_MTX_SIZE] = (double)nmax1;
        }

      /* Avoid compiler warning */

      loc2 = -1;

      /* Add new cells as neighbours for the neighbours of the old cell */

      loc1 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
      while (loc1 > VALID_PTR)
        {
          /* Pointer to cell */

          loc2 = (long)RDB[loc1 + RMX_CELL_BOUND_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Get pointer to mesh */

          msh1 = (long)RDB[loc2 + RMX_CELL_PTR_MESH];
          CheckPointer(FUNCTION_NAME, "(msh1)", DATA_ARRAY, msh1);

          /* Mesh index */

          l = (long)RDB[loc2 + RMX_CELL_MESH_IDX];

          /* Get boundaries */

          MeshCellBounds(msh1, l, &lims1[0], &lims1[1], &lims1[2],
                         &lims1[3], &lims1[4], &lims1[5]);

          /* Loop over neighbours */

          loc3 = (long)RDB[loc2 + RMX_CELL_PTR_BOUNDS];
          while (loc3 > VALID_PTR)
            {
              /* Compare pointers */

              if ((long)RDB[loc3 + RMX_CELL_BOUND_PTR_CELL] == loc0)
                break;

              /* Next */

              loc3 = NextItem(loc3);
            }

          /* Check */

          if (loc3 < VALID_PTR)
            Warn(FUNCTION_NAME, "Old neighbour not found");
          else
            {
              /* Remove cell from list */

              RemoveItem(loc3);

              /* Adjust size */

              WDB[loc2 + RMX_CELL_MTX_SIZE] =
                RDB[loc2 + RMX_CELL_MTX_SIZE] - 1.0;
            }

          /* Loop over mesh */

          for (k0 = 0; k0 < nz; k0++)
          for (j0 = 0; j0 < ny; j0++)
          for (i0 = 0; i0 < nx; i0++)
            {
              /* Pointer to data */

              loc3 = ReadMeshPtr(msh0, i0, j0, k0);
              CheckPointer(FUNCTION_NAME, "(loc3)", DATA_ARRAY, loc3);

              loc3 = (long)RDB[loc3];
              CheckPointer(FUNCTION_NAME, "(loc3)", DATA_ARRAY, loc3);

              /* Mesh index */

              l = (long)RDB[loc3 + RMX_CELL_MESH_IDX];

              /* Get boundaries */

              MeshCellBounds(msh0, l, &lims0[0], &lims0[1], &lims0[2],
                             &lims0[3], &lims0[4], &lims0[5]);

              /* Check joining surface */

              if (TestNeighbour(lims0, lims1))
                {
                  /* Create current */

                  ptr = NewItem(loc2 + RMX_CELL_PTR_BOUNDS,
                                RMX_CELL_BOUND_BLOCK_SIZE);

                  /* Put pointer and reset indexes */

                  WDB[ptr + RMX_CELL_BOUND_PTR_CELL] = (double)loc3;
                  WDB[ptr + RMX_CELL_BOUND_FWD_IDX] = -1.0;
                  WDB[ptr + RMX_CELL_BOUND_ADJ_IDX] = -1.0;

                  /* Adjust size */

                  WDB[loc2 + RMX_CELL_MTX_SIZE] =
                    RDB[loc2 + RMX_CELL_MTX_SIZE] + 1.0;
                }
            }

          /* Next */

          loc1 = NextItem(loc1);
        }

      /***********************************************************************/

      /***** Allocate memory for data ****************************************/

      /* Loop over mesh */

      for (k0 = 0; k0 < nz; k0++)
      for (j0 = 0; j0 < ny; j0++)
      for (i0 = 0; i0 < nx; i0++)
        {
          /* Pointer to data */

          loc1 = ReadMeshPtr(msh0, i0, j0, k0);
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          loc1 = (long)RDB[loc1];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Number of cells */

          nmax1 = (long)RDB[loc1 + RMX_CELL_MTX_SIZE];

          /* Check size */

          ptr = (long)RDB[loc1 + RMX_CELL_PTR_BOUNDS];
          CheckPointer(FUNCTION_NAME, "(ptr12)", DATA_ARRAY, ptr);

          if (nmax1 != ListSize(ptr))
            Die(FUNCTION_NAME, "Mismatch in size");

          /* Allocate memory */

          RecycleRMXData(rmx, loc1, ng, nmax1, NO);
        }

      /* Next cell */

      loc0 = NextItem(loc0);

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Reallocate memory for neighbour cells *******************************/

  /* NOTE: Tässä hukataan nyt varattua muistia, mutta en onnistunut */
  /* toteuttamaan tota siten, että nuo vanhat structured dumpataan. */

  /* Loop over data */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Check split */

      if ((long)RDB[loc0 + RMX_CELL_SPLIT_FLAG] == NO)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Loop over neighbours */

      loc1 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
      while (loc1 > VALID_PTR)
        {
          /* Pointer to cell */

          loc2 = (long)RDB[loc1 + RMX_CELL_BOUND_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Check if already processed */

          if ((long)RDB[loc2 + RMX_CELL_MTX_ADA_PROC] == NO)
            {
              /* Number of neighbous */

              nmax1 = (long)RDB[loc2 + RMX_CELL_MTX_SIZE];

              /* Check size */

              ptr = (long)RDB[loc2 + RMX_CELL_PTR_BOUNDS];
              CheckPointer(FUNCTION_NAME, "(ptr13)", DATA_ARRAY, ptr);

              if (nmax1 != ListSize(ptr))
                Die(FUNCTION_NAME, "Mismatch in size");

              /* Allocate memory */

              RecycleRMXData(rmx, loc2, ng, nmax1, NO);

              /* Set flag */

              WDB[loc2 + RMX_CELL_MTX_ADA_PROC] = (double)YES;
            }

          /* Next */

          loc1 = NextItem(loc1);
        }

      /***********************************************************************/

      /* Copy pointer */

      loc1 = loc0;

      /* Next cell */

      loc0 = NextItem(loc0);

      /* Remove */

      RemoveItem(loc1);

      /* Move to dump */

      AddItem(DATA_RMX_ADA_PTR_CELL_DUMP, loc1);
    }

  /***************************************************************************/

  /***** Put neighbour indexes ***********************************************/

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Check split flag */

      if ((long)RDB[loc0 + RMX_CELL_SPLIT_FLAG] == YES)
        Die(FUNCTION_NAME, "Splitted cell");

      /* Maximum number of neighbours */

      nmax1 = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax1", "", nmax1, 1, 10000);

      /* Pointer to bounds */

      loc1 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Check size */

      if (nmax1 != ListSize(loc1))
        Die(FUNCTION_NAME, "Mismatch in size");

      /* Loop over boundaries */

      for (j = 0; j < nmax1; j++)
        {
          /* Pointer to cell */

          ptr = (long)RDB[loc1 + RMX_CELL_BOUND_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(ptr14)", DATA_ARRAY, ptr);

          /* Maximum number of neighbours */

          nmax2 = (long)RDB[ptr + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax2", "", nmax2, 1, 10000);

          /* Pointer to bounds */

          loc2 = (long)RDB[ptr + RMX_CELL_PTR_BOUNDS];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Check size */

          if (nmax2 != ListSize(loc2))
            Die(FUNCTION_NAME, "Mismatch in size");

          /* Loop over boundaries */

          for (i = 0; i < nmax2; i++)
            {
              /* Compare pointers */

              if ((long)RDB[loc2 + RMX_CELL_BOUND_PTR_CELL] == loc0)
                break;

              /* Next boundary */

              loc2 = NextItem(loc2);
            }

          /* Check */

          if (i == nmax2)
            Die(FUNCTION_NAME, "Neighbour not found %ld %ld", loc0, i);

          /* Put indexes */

          WDB[loc1 + RMX_CELL_BOUND_FWD_IDX] = (double)i;
          WDB[loc2 + RMX_CELL_BOUND_ADJ_IDX] = (double)j;

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Check connectivity **************************************************/

  /* Loop over mesh data */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Loop over bounds */
      break;
      loc1 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
      while (loc1 > VALID_PTR)
        {
          /* Pointer to cell */

          ptr = (long)RDB[loc1 + RMX_CELL_BOUND_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(ptr15)", DATA_ARRAY, ptr);

          /* Check connection to self */

          if (ptr == loc0)
            Die(FUNCTION_NAME, "Connection to self");

          /* Loop over bounds */

          loc2 = (long)RDB[ptr + RMX_CELL_PTR_BOUNDS];
          while (loc2 > VALID_PTR)
            {
              /* Compare */

              if ((long)RDB[loc2 + RMX_CELL_BOUND_PTR_CELL] == loc0)
                break;

              /* Next */

              loc2 = NextItem(loc2);
            }

          /* Check pointer */

          if (loc2 < VALID_PTR)
            Die(FUNCTION_NAME, "Error in connectivity");

          /* Check that connected cell is part of list */

          ptr = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          while (ptr > VALID_PTR)
            {
              /* Compare */

              if ((long)RDB[loc2 + RMX_CELL_BOUND_PTR_CELL] == ptr)
                break;

              /* Next */

              ptr = NextItem(ptr);
            }

          /* Check pointer */

          if (ptr < VALID_PTR)
            Die(FUNCTION_NAME, "Error in connectivity");

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Check duplicate cell indexes ****************************************/

#ifdef DEBUG2

  /* Loop over data */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Loop over remaining */

      loc1 = NextItem(loc0);
      while (loc1 > VALID_PTR)
        {
          /* Check index */

          if (RDB[loc0 + RMX_CELL_IDX] == RDB[loc1 + RMX_CELL_IDX])
            Die(FUNCTION_NAME, "Duplicate index %ld",
                (long)RDB[loc0 + RMX_CELL_IDX]);

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

#endif

  /***************************************************************************/

  /* Close cell list */

  ptr = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr16)", DATA_ARRAY, ptr);
  CloseList(ptr);

  /* Deny memory allocation */

  Mem(MEM_DENY);

  /* Expand private arrays */

  ExpandPrivateArrays();

  /* Test */

  TestRMX(rmx);

  /* Return number of splits */

  return nsplit;

  /***************************************************************************/
}

/*****************************************************************************/


/***** Function to test if cells are neighbours ******************************/

long TestNeighbour(double *lims0, double *lims1)
{
  /* Check if common surface */

  if ((fabs(lims0[0] - lims1[1]) < 1E-6) ||
      (fabs(lims0[1] - lims1[0]) < 1E-6))
    {
      /* Check other common surfaces */

       if ((fabs(lims0[2] - lims1[3]) < 1E-6) ||
           (fabs(lims0[3] - lims1[2]) < 1E-6))
         return NO;
       else if ((fabs(lims0[4] - lims1[5]) < 1E-6) ||
                (fabs(lims0[5] - lims1[4]) < 1E-6))
          return NO;

      /* Check that surfaces overlap */

       if ((lims0[3] < lims1[2]) || (lims1[3] < lims0[2]) ||
           (lims0[5] < lims1[4]) || (lims1[5] < lims0[4]))
         return NO;

       return YES;
    }
  else if ((fabs(lims0[2] - lims1[3]) < 1E-6) ||
           (fabs(lims0[3] - lims1[2]) < 1E-6))
    {
      /* Check other common surfaces */

       if ((fabs(lims0[0] - lims1[1]) < 1E-6) ||
           (fabs(lims0[1] - lims1[0]) < 1E-6))
         return NO;
       else if ((fabs(lims0[4] - lims1[5]) < 1E-6) ||
                (fabs(lims0[5] - lims1[4]) < 1E-6))
          return NO;

      /* Check that surfaces overlap */

       if ((lims0[1] < lims1[0]) || (lims1[1] < lims0[0]) ||
           (lims0[5] < lims1[4]) || (lims1[5] < lims0[4]))
         return NO;

       return YES;
    }
  else if ((fabs(lims0[4] - lims1[5]) < 1E-6) ||
           (fabs(lims0[5] - lims1[4]) < 1E-6))
    {
      /* Check other common surfaces */

       if ((fabs(lims0[0] - lims1[1]) < 1E-6) ||
           (fabs(lims0[1] - lims1[0]) < 1E-6))
         return NO;
       else if ((fabs(lims0[2] - lims1[3]) < 1E-6) ||
                (fabs(lims0[3] - lims1[2]) < 1E-6))
          return NO;

      /* Check that surfaces overlap */

       if ((lims0[3] < lims1[2]) || (lims1[3] < lims0[2]) ||
           (lims0[1] < lims1[0]) || (lims1[1] < lims0[0]))
         return NO;

       return YES;
    }

  /* Not a neighbour */

  return NO;
}

/*****************************************************************************/

/***** Test routine for adjusted mesh ****************************************/

void TestRMX(long rmx)
{
  long loc0, loc1, ptr, msh, msh0, i, dim, k;
  double lims0[6], lims1[6], x, y, z, dx, dy, dz;
  unsigned long seed;
  return;
  fprintf(outp, "Testing...\n");

  /* Dimensions */

  dim = 2;

  /* Initialize rng */

  seed = ReInitRNG(423235);
  *SEED = seed;

  /* Pointer to mesh */

  if ((msh = (long)RDB[rmx + RMX_PTR_MESH]) < VALID_PTR)
    Die(FUNCTION_NAME, "No mesh");

  /* Loop over mesh data */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Get pointer to mesh */

      msh0 = (long)RDB[loc0 + RMX_CELL_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh0)", DATA_ARRAY, msh0);

      /* Get index */

      i = (long)RDB[loc0 + RMX_CELL_MESH_IDX];
      CheckValue(FUNCTION_NAME, "i", "", i, 0, 1000000000);

      /* Get boundaries */

      MeshCellBounds(msh0, i, &lims0[0], &lims0[1], &lims0[2], &lims0[3],
                     &lims0[4], &lims0[5]);

      /* Sample random points in cell */

      for (i = 0; i < 1000; i++)
        {
          /* Sample point */

          x = RandF(0)*(lims0[1] - lims0[0]) + lims0[0];
          y = RandF(0)*(lims0[3] - lims0[2]) + lims0[2];

          if (dim == 3)
            z = RandF(0)*(lims0[5] - lims0[4]) + lims0[4];
          else
            z = 0.0;

          /* Find mesh cell */

          if ((loc1 = MeshPtr(msh, x, y, z)) < VALID_PTR)
            Die(FUNCTION_NAME, "Mesh error");

          /* Pointer to structure */

          loc1 = (long)RDB[loc1];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Compare */

          if (loc0 != loc1)
            Die(FUNCTION_NAME, "Mismatch in cell 1");
        }

      /* Calculate dimensions */

      dx = lims0[1] - lims0[0];
      dy = lims0[3] - lims0[2];
      dz = lims0[5] - lims0[4];

      /* Sample random points outside cell */

      for (i = 0; i < 1000; i++)
        {
          /* Sample point */

          x = RandF(0)*(lims0[1] - lims0[0]) + lims0[0];
          y = RandF(0)*(lims0[3] - lims0[2]) + lims0[2];

          if (dim == 3)
            z = RandF(0)*(lims0[5] - lims0[4]) + lims0[4];
          else
            z = 0.0;

          k = (long)(dim*RandF(0));

          if (k == 0)
            {
              if (RandF(0) < 0.5)
                x = 0.5*(lims0[0] + lims0[1]) + 0.5001*dx;
              else
                x = 0.5*(lims0[0] + lims0[1]) - 0.5001*dx;
            }
          else if (k == 1)
            {
              if (RandF(0) < 0.5)
                y = 0.5*(lims0[2] + lims0[3]) + 0.5001*dy;
              else
                y = 0.5*(lims0[2] + lims0[3]) - 0.5001*dy;
            }
          else
            {
              if (RandF(0) < 0.5)
                z = 0.5*(lims0[4] + lims0[5]) + 0.5001*dz;
              else
                z = 0.5*(lims0[4] + lims0[5]) - 0.5001*dz;
            }

          /* Find mesh cell */

          if ((loc1 = MeshPtr(msh, x, y, z)) < VALID_PTR)
            continue;

          /* Pointer to structure */

          loc1 = (long)RDB[loc1];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Loop over neighbours */

          ptr = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
          while (ptr > VALID_PTR)
            {
              /* Compare pointers */

              if ((long)RDB[ptr + RMX_CELL_BOUND_PTR_CELL] == loc1)
                break;

              /* Next */

              ptr = NextItem(ptr);
            }

          if (ptr < VALID_PTR)
            Die(FUNCTION_NAME, "Mismatch in cell 2");
        }

      /* Loop over neighbours */

      ptr = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
      while (ptr > VALID_PTR)
        {
          /* Sample random points outside cell */

          for (i = 0; i < 100000; i++)
            {
              /* Sample point */

              x = RandF(0)*(lims0[1] - lims0[0]) + lims0[0];
              y = RandF(0)*(lims0[3] - lims0[2]) + lims0[2];

              if (dim == 3)
                z = RandF(0)*(lims0[5] - lims0[4]) + lims0[4];
              else
                z = 0.0;

              k = (long)(dim*RandF(0));

              if (k == 0)
                {
                  if (RandF(0) < 0.5)
                    x = 0.5*(lims0[0] + lims0[1]) + 0.5001*dx;
                  else
                    x = 0.5*(lims0[0] + lims0[1]) - 0.5001*dx;
                }
              else if (k == 1)
                {
                  if (RandF(0) < 0.5)
                    y = 0.5*(lims0[2] + lims0[3]) + 0.5001*dy;
                  else
                    y = 0.5*(lims0[2] + lims0[3]) - 0.5001*dy;
                }
              else
                {
                  if (RandF(0) < 0.5)
                    z = 0.5*(lims0[4] + lims0[5]) + 0.5001*dz;
                  else
                    z = 0.5*(lims0[4] + lims0[5]) - 0.5001*dz;
                }

              /* Find mesh cell */

              if ((loc1 = MeshPtr(msh, x, y, z)) < VALID_PTR)
                continue;

              /* Pointer to structure */

              loc1 = (long)RDB[loc1];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

              /* Compare pointers */

              if ((long)RDB[ptr + RMX_CELL_BOUND_PTR_CELL] == loc1)
                break;
            }

          /* Check count */

          if (i == 100000)
            {
              ptr = (long)RDB[ptr + RMX_CELL_BOUND_PTR_CELL];

              /* Get pointer to mesh */
              printf("%ld %ld %ld\n", msh,
                     msh0, (long)RDB[ptr + RMX_CELL_PTR_MESH]);
              msh0 = (long)RDB[ptr + RMX_CELL_PTR_MESH];
              CheckPointer(FUNCTION_NAME, "(msh0)", DATA_ARRAY, msh0);

              /* Get index */

              i = (long)RDB[ptr + RMX_CELL_MESH_IDX];
              CheckValue(FUNCTION_NAME, "i", "", i, 0, 1000000000);

              /* Get boundaries */

              MeshCellBounds(msh0, i, &lims1[0], &lims1[1], &lims1[2],
                             &lims1[3], &lims1[4], &lims1[5]);

              printf("%12.5E %12.5E : %12.5E %12.5E : %12.5E %12.5E\n",
                     lims0[0], lims0[1], lims0[2],
                     lims0[3], lims0[4], lims0[5]);

              printf("%12.5E %12.5E : %12.5E %12.5E : %12.5E %12.5E\n",
                     lims1[0], lims1[1], lims1[2],
                     lims1[3], lims1[4], lims1[5]);

              printf("%E %E %E\n",
                     lims0[1] - lims1[0],
                     lims0[3] - lims1[2],
                     lims0[5] - lims1[4]);

              printf("%E %E %E\n",
                     lims1[1] - lims0[0],
                     lims1[3] - lims0[2],
                     lims1[5] - lims0[4]);

              Die(FUNCTION_NAME, "Mismatch in cell 3");
            }

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
