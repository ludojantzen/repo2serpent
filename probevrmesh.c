/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : probevrmesh.c                                  */
/*                                                                           */
/* Created:       2018/11/06 (JLe)                                           */
/* Last modified: 2019/10/08 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Samples random points in geometry to adapt VR mesh.          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProbeVRMesh:"

/*****************************************************************************/

void ProbeVRMesh(long rmx)
{
  long msh, loc0, loc1, itp, idx, id, cell, mat, n, nt, sz, i, j, nmax, ptr;
  long splitmax, nx, ny, nz;
  unsigned long seed;
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, u, v, w, rho, d;
  double *split, *split0, rhmax;

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */
  /* (this is needed to enable call to cell test routine) */

  ExpandPrivateArrays();

  /* Get pointer to iteration data */

  itp = (long)RDB[rmx + RMX_PTR_ITER];
  CheckPointer(FUNCTION_NAME, "(itp)", DATA_ARRAY, itp);

  /* Get number of tracks */

  nt =  RDB[itp + WWD_ITER_SPLIT_TRACKS];
  CheckValue(FUNCTION_NAME, "nt", "", nt, 1, INFTY);

  /* Maximum number of splits */

  splitmax = (long)RDB[DATA_MAX_RMX_SPLIT_FLAGS];
  CheckValue(FUNCTION_NAME, "splitmax", "", splitmax, 1, 100);

  /* Get mesh size */

  idx = (long)RDB[DATA_VR_PROBE_IDX];
  WDB[DATA_VR_PROBE_IDX] = RDB[DATA_VR_PROBE_IDX] + 1.0;

  /***************************************************************************/

  /***** Prepare mesh ********************************************************/

  /* Pointer to mesh data */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get number of cells */

  nmax = ListSize(loc0);
  CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 1000000000);

  /* Calculate number of points per cell */

  nt = (long)((double)nt/((double)nmax));

  /* Put minimum (average 5 per direction) */

  if (nt < 15)
    nt = 15;

  /* Loop over data and reset flags */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Reset common flag */

      WDB[loc0 + RMX_CELL_SPLIT_FLAG] = (double)NO;

      /* Loop over split flags */

      ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      for (i = 0; i < splitmax; i++)
        WDB[ptr++] = (double)NO;

      /* Get index */

      i = (long)RDB[loc0 + RMX_CELL_MESH_IDX];
      CheckValue(FUNCTION_NAME, "i", "", i, 0, nmax - 1);

      /* Pointer to mesh */

      msh = (long)RDB[loc0 + RMX_CELL_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get boundaries */

      MeshCellBounds(msh, i, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

      /* Get split size */

      nx = (long)RDB[itp + WWD_ITER_SPLIT_NX];
      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 100);

      ny = (long)RDB[itp + WWD_ITER_SPLIT_NY];
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 100);

      nz = (long)RDB[itp + WWD_ITER_SPLIT_NZ];
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 100);

      /* Loop over split criteria */

      loc1 = (long)RDB[itp + WWD_ITER_SPLIT_PTR_DENS];
      while (loc1 > VALID_PTR)
        {
          /* Check */

          if ((xmax - xmin)/((double)nx)
              < RDB[loc1 + WWD_ITER_SPLIT_DENS_XMIN])
            nx = 1;

          if ((ymax - ymin)/((double)ny)
              < RDB[loc1 + WWD_ITER_SPLIT_DENS_YMIN])
            ny = 1;

          if ((zmax - zmin)/((double)nz)
              < RDB[loc1 + WWD_ITER_SPLIT_DENS_ZMIN])
            nz = 1;

          /* Check if split */

          if (nx*ny*nz == 1)
            {
              /* Set flag to -1 to indicate skip */

              WDB[loc0 + RMX_CELL_SPLIT_FLAG] = -1.0;

              /* Break loop */

              break;
            }

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Sample tracks inside cells ******************************************/

  /* Loop over cells */

#ifdef OPEN_MP
#pragma omp parallel for private(n, id, seed, x, y, z, u, v, w, cell, mat, loc0, rho, xmin, xmax, ymin, ymax, zmin, zmax, i, j, msh, ptr ,d, loc1)
#endif
  for (j = 0; j < nmax; j++)
    {
      /* Get Open MP thread id */

      id = OMP_THREAD_NUM;

      /* Init random number sequence */

      seed = ReInitRNG(idx*nmax + j);
      SEED[id*RNG_SZ] = seed;

      /* Get pointer to data */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      loc0 = ListPtr(loc0, j);
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Check common flag */

      if ((long)RDB[loc0 + RMX_CELL_SPLIT_FLAG] == -1)
        {
          /* Reset flag */

          WDB[loc0 + RMX_CELL_SPLIT_FLAG] = (double)NO;

          /* Cycle loop */

          continue;
        }

      /* Loop over split flags */

      ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      for (i = 0; i < splitmax; i++)
        if ((long)RDB[ptr++] == YES)
          break;

      /* Check if already split (previous outer cycles) */

      if (i < splitmax)
        continue;

      /* Get index */

      i = (long)RDB[loc0 + RMX_CELL_MESH_IDX];
      CheckValue(FUNCTION_NAME, "i", "", i, 0, nmax - 1);

      /* Pointer to mesh */

      msh = (long)RDB[loc0 + RMX_CELL_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get boundaries */

      MeshCellBounds(msh, i, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

      /* Loop over tracks */

      for (n = 0; n < nt; n++)
        {
          /* Sample starting point */

          if ((i = (long)(3.0*RandF(id))) == 0)
            {
              /* X-axis */

              u = 1.0;
              v = 0.0;
              w = 0.0;

              x = xmin + EXTRAP_L;
              y = RandF(id)*(ymax - ymin) + ymin;
              z = RandF(id)*(zmax - zmin) + zmin;
            }
          else if (i == 1)
            {
              /* Y-Axis */

              u = 0.0;
              v = 1.0;
              w = 0.0;

              x = RandF(id)*(xmax - xmin) + xmin;
              y = ymin + EXTRAP_L;
              z = RandF(id)*(zmax - zmin) + zmin;
            }
          else if (i == 2)
            {
              /* Z-axis */

              u = 0.0;
              v = 0.0;
              w = 1.0;

              x = RandF(id)*(xmax - xmin) + xmin;
              y = RandF(id)*(ymax - ymin) + ymin;
              z = zmin + EXTRAP_L;
            }
          else
            {
              /* Avoid compiler warning */

              x = 0.0;
              y = 0.0;
              z = 0.0;

              u = 0.0;
              v = 0.0;
              w = 0.0;

              /* Shouldn't be here */

              Die(FUNCTION_NAME, "WTF?");
            }

#ifdef DEBUG

          /* Find mesh cell (needed for debuggin only) */

          if ((ptr = MeshPtr(msh, x, y, z)) < VALID_PTR)
            Error(rmx, "Point (%E, %E, %E) is outside mesh", x, y, z);

          /* Check pointer */

          if ((long)RDB[ptr] != loc0)
            Die(FUNCTION_NAME, "Pointer error");

#endif
          /* Loop until outside */

          do
            {
              /* Check true symmetry */

              if (TestUniSym(x, y, z) == YES)
                {
                  /* Get pointer to root universe */

                  ptr =  (long)RDB[DATA_PTR_ROOT_UNIVERSE];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Get pointer to symmetry */

                  ptr = (long)RDB[ptr + UNIVERSE_PTR_SYM];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Calculate distance to boundary */

                  if ((d = SymmetryBoundary(ptr, x, y, z, u, v, w)) < 0.0)
                    break;

                  /* Move to new position */

                  x = x + (d + EXTRAP_L)*u;
                  y = y + (d + EXTRAP_L)*v;
                  z = z + (d + EXTRAP_L)*w;

                  /* Check */

                  if ((x > xmax) || (y > ymax) || (z > zmax))
                    break;
                }

              /* Get cell */

              if ((cell = WhereAmI(x, y, z, u, v, w, id)) < 0)
                Error(0, "Geometry error at %E %E %E", x, y, z);

              /* Check if cell has material */

              if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                {
                  /* Get material pointer */

                  mat = MatPtr(mat, id);

                  /* Pointer to split vector */

                  ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Loop over split criteria */

                  loc1 = (long)RDB[itp + WWD_ITER_SPLIT_PTR_DENS];
                  while (loc1 > VALID_PTR)
                    {
                      /* Get maximum density */

                      rhmax = RDB[loc1 + WWD_ITER_SPLIT_DENS_RHO];

                      /* Get density */

                      if (rhmax < 0.0)
                        rho = RDB[mat + MATERIAL_MDENS];
                      else
                        rho = RDB[mat + MATERIAL_ADENS];

                      /* Check material density */

                      if (rho >= fabs(rhmax))
                        {
                          /* Set split flag */

                          WDB[ptr] = (double)YES;

                          /* Break loop */

                          break;
                        }

                      /* Update pointer */

                      ptr++;

                      /* Next */

                      loc1 = NextItem(loc1);
                    }
                }

              /* Check if already split (tässä pitää katsoa kaikki */
              /* kriteerit, jotta flagi menee varmasti myös sinne  */
              /* tiheimpään). */

              ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              for (i = 0; i < splitmax; i++)
                if ((long)RDB[ptr++] == NO)
                  break;

              /* Check if already split */

              if (i == splitmax)
                break;

              /* Calculate distance to boundary */

              d = NearestBoundary(id);

              /* Move to new position */

              x = x + (d + EXTRAP_L)*u;
              y = y + (d + EXTRAP_L)*v;
              z = z + (d + EXTRAP_L)*w;
            }
          while ((x < xmax) && (y < ymax) && (z < zmax));

          /* Check if already split */

          ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          for (i = 0; i < splitmax; i++)
            if ((long)RDB[ptr++] == NO)
              break;

          /* Check if already split */

          if (i == splitmax)
            break;
        }
    }

  /***************************************************************************/

  /***** Distribute split flags to MPI tasks *********************************/

  /* Avoid compiler warning */

  split0 = NULL;
  split = NULL;
  sz = 0;

#ifdef MPI

  /* Pointer to mesh data */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get size */

  sz = ListSize(loc0);
  CheckValue(FUNCTION_NAME, "sz", "", sz, 1, 100000000);

  /* Allocate memory for data */

  split0 = (double *)Mem(MEM_ALLOC, splitmax*sz, sizeof(double));
  split = (double *)Mem(MEM_ALLOC, splitmax*sz, sizeof(double));

  /* Read flags */

  n = 0;
  while (loc0 > VALID_PTR)
    {
      /* Pointer to splits */

      ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Put flags */

      for (i = 0; i < splitmax; i++)
        split[n++] = RDB[ptr + i];

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Reduce data */

  MPITransfer(split, split0, sz, 0, MPI_METH_RED);

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Broadcast data to other tasks */

  MPITransfer(split0, NULL, sz, 0, MPI_METH_BC);

  /* Put flags */

  n = 0;

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to splits */

      ptr = (long)RDB[loc0 + RMX_CELL_PTR_SPLIT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over flags */

      for (i = 0; i < splitmax; i++)
        if (split0[n++] > 0.0)
          WDB[ptr + i] = (double)YES;

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Free temporary arrays */

  Mem(MEM_FREE, split);
  Mem(MEM_FREE, split0);

  /* Synchronise */

  MPI_Barrier(my_comm);

#endif

  /***************************************************************************/
}

/*****************************************************************************/
