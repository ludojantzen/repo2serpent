/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processrmx.c                                   */
/*                                                                           */
/* Created:       2016/04/25 (JLe)                                           */
/* Last modified: 2019/10/03 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Processes response matrix stuff                              */
/*                                                                           */
/* Comments: - Ton muistinvarauksen voisi hoitaa samalla aliohjelmalla       */
/*             kuin adjustrmx.c:ssä.                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessRMX:"

/*****************************************************************************/

void ProcessRMX()
{
  long msh, loc0, loc1, loc2, ptr, nx, ny, nz, i0, j0, k0, i, j, k, rmx, tgt;
  long det, nmax1, nmax2, type, erg, ng;
  double min0, max0, min1, max1, min2, max2, p, r0, r1, lims[6];
  double a0, a1, x0, x1, y0, y1, z0, z1;

  /* Check if response matrix calculation is on */

  if ((long)RDB[DATA_RMTX_CALC] == NO)
    return;

  fprintf(outp, "Processing response matrix structures...\n");

  /* Memory count */

  WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] + (double)MemCount();

  /* Create cell dump list for memory management in self-adaptive mode */

  loc0 = NewItem(DATA_RMX_ADA_PTR_CELL_DUMP, RMX_CELL_BLOCK_SIZE);
  WDB[loc0 + RMX_CELL_MTX_SIZE] = -1.0;

  /***************************************************************************/

  /***** Link detectors ******************************************************/

  /* Loop over response matrixes */

  rmx = (long)RDB[DATA_PTR_RMX0];
  while (rmx > VALID_PTR)
    {
      /* Reset TLE flag */

      WDB[rmx + RMX_USE_TLE] = (double)NO;

      /* Check convergence acceleration */

      if ((long)RDB[DATA_RMX_CONVG_ACC] == YES)
        {
          /* Detector shouldn't be defined */

          if ((long)RDB[rmx + RMX_PTR_DET] > VALID_PTR)
            Die(FUNCTION_NAME, "Detector structure already exists");

          /* Create detector structure */

          det = NewItem(rmx + RMX_PTR_DET, RMX_DET_BLOCK_SIZE);

          /* Reset name and weight */

          WDB[det + RMX_DET_PTR_DET] = -1;
          WDB[det + RMX_DET_WGT] = 1.0;

          /* Set particle type */

          WDB[rmx + RMX_PARTICLE_TYPE] = (double)PARTICLE_TYPE_NEUTRON;
        }
      else
        {
          /* Pointer to detector array */

          if ((loc0 = (long)RDB[rmx + RMX_PTR_DET]) < VALID_PTR)
            {
              /* Check type */

              if ((long)RDB[rmx + RMX_MODE] != RMX_MODE_GVR)
                Die(FUNCTION_NAME, "No responses");

              /* Create detector structure for GVR */

              loc0 = NewItem(rmx + RMX_PTR_DET, RMX_DET_BLOCK_SIZE);

              /* Reset name and weight */

              WDB[loc0 + RMX_DET_PTR_DET] = -1;
              WDB[loc0 + RMX_DET_WGT] = 1.0;

              /* Set TLE flag */

              WDB[rmx + RMX_USE_TLE] = (double)YES;

              /* Set rmx particle type */

              if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
                {
                  /* Check if photon transport mode is on */

                  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
                    Note(rmx, "Weight window mesh generated for neutrons");

                  /* Set type */

                  WDB[rmx + RMX_PARTICLE_TYPE] = (long)PARTICLE_TYPE_NEUTRON;
                }
              else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
                WDB[rmx + RMX_PARTICLE_TYPE] = (long)PARTICLE_TYPE_GAMMA;
              else
                Die(FUNCTION_NAME, "Transport mode not set");
            }
          else
            {
              /* Check for multiple responses in GVR mode */

              if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_GVR)
                if (ListSize(loc0) > 1)
                  Error(rmx, "Multiple responses in GVR mode not allowed");

              /* Check for multiple responses in multi-response mode */

              if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_MULTI)
                if (ListSize(loc0) == 1)
                  Error(rmx, "Multi-response mode with single response");

              /* Loop over array */

              while (loc0 > VALID_PTR)
                {
                  /* Find detector */

                  det = (long)RDB[DATA_PTR_DET0];
                  while (det > VALID_PTR)
                    {
                      /* Compare names */

                      if (CompareStr(det + DET_PTR_NAME,
                                     loc0 + RMX_DET_PTR_DET))
                        break;

                      /* Next */

                      det = NextItem(det);
                    }

                  /* Check pointer */

                  if (det < VALID_PTR)
                    Error(rmx, "Detector %s not defined",
                          GetText(loc0 + RMX_DET_PTR_DET));

                  /* Put pointer */

                  WDB[loc0 + RMX_DET_PTR_DET] = (double)det;

                  /* Set rmx particle type */

                  if ((long)RDB[rmx + RMX_PARTICLE_TYPE] == 0)
                    WDB[rmx + RMX_PARTICLE_TYPE] = RDB[det + DET_PARTICLE];
                  else if (RDB[rmx + RMX_PARTICLE_TYPE] !=
                           RDB[det + DET_PARTICLE])
                    Error(rmx, "Multiple particle types in responses");

                  /* Next */

                  loc0 = NextItem(loc0);
                }
            }
        }

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Loop over detectors */

          det = (long)RDB[rmx + RMX_PTR_DET];
          while (det > VALID_PTR)
            {
              /* Create structure */

              ptr = NewItem(loc0 + RMX_CELL_PTR_DET, RMX_DET_BLOCK_SIZE);

              /* Copy data */

              for (i = LIST_DATA_SIZE; i < RMX_DET_BLOCK_SIZE - 1; i++)
                WDB[ptr + i] = RDB[det + i];

              /* Link parent */

              WDB[ptr + RMX_DET_PTR_DET0] = (double)det;

              /* Next detector */

              det = NextItem(det);
            }

          /* Next cell */

          loc0 = NextItem(loc0);
        }

      /* Next */

      rmx = NextItem(rmx);
    }

  /****************************************************************************/

  /***** Link energy grid *****************************************************/

  /* Loop over response matrixes */

  rmx = (long)RDB[DATA_PTR_RMX0];
  while (rmx > VALID_PTR)
    {
      /* Check if data is from file */

      if ((long)RDB[rmx + RMX_FROM_FILE] == YES)
        {
          /* Pointer to next */

          rmx = NextItem(rmx);

          /* Cycle loop */

          continue;
        }

      /* Check if defined */

      if ((long)RDB[rmx + RMX_PTR_EGRID] < VALID_PTR)
        {
          /* Set number of groups */

          WDB[rmx + RMX_NG] = 1.0;

          /* Pointer to next */

          rmx = NextItem(rmx);

          /* Cycle loop */

          continue;
        }

      /* Loop over energy grids */

      erg = (long)RDB[DATA_PTR_ENE0];
      while (erg > VALID_PTR)
        {
          /* Compare */

          if (CompareStr(rmx + RMX_PTR_EGRID, erg + ENE_PTR_NAME))
            break;

          /* Next */

          erg = NextItem(erg);
        }

      /* Check if found */

      if (erg < VALID_PTR)
        Error(rmx, "Energy grid %s not defined", GetText(rmx + RMX_PTR_EGRID));

      /* Set number of groups */

      if ((ng = (long)RDB[erg + ENE_NB]) < 1)
        Error(erg, "No energy bins in grid %s\n", GetText(erg + ENE_PTR_NAME));
      else
        WDB[rmx + RMX_NG] = (double)ng;

      /* Pointer to grid */

      erg = (long)RDB[erg + ENE_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Put pointer */

      WDB[rmx + RMX_PTR_EGRID] = (double)erg;

      /* Next */

      rmx = NextItem(rmx);
    }

  /***************************************************************************/

  /***** Process mesh ********************************************************/

  /* Loop over response matrixes */

  rmx = (long)RDB[DATA_PTR_RMX0];
  while (rmx > VALID_PTR)
    {
      /* Check if data is from file */

      if ((long)RDB[rmx + RMX_FROM_FILE] == YES)
        {
          /* Pointer to next */

          rmx = NextItem(rmx);

          /* Cycle loop */

          continue;
        }

      /* Get pointer to mesh */

      if ((msh = (long)RDB[rmx + RMX_PTR_MESH]) < VALID_PTR)
        Error(rmx, "Mesh parameters not defined");

      /* Get type */

      type = (long)RDB[msh + MESH_TYPE];
      CheckValue(FUNCTION_NAME, "type", "", type, 1, 8);

      /* Get size */

      nx = (long)RDB[msh + MESH_N0];
      ny = (long)RDB[msh + MESH_N1];
      nz = (long)RDB[msh + MESH_N2];

      /* Get limits */

      min0 = RDB[msh + MESH_MIN0];
      max0 = RDB[msh + MESH_MAX0];
      min1 = RDB[msh + MESH_MIN1];
      max1 = RDB[msh + MESH_MAX1];
      min2 = RDB[msh + MESH_MIN2];
      max2 = RDB[msh + MESH_MAX2];

      /* Get number of energy groups */

      ng = (long)RDB[rmx + RMX_NG];
      CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

      /* Loop over mesh */

      for (i0 = 0; i0 < nx; i0++)
      for (j0 = 0; j0 < ny; j0++)
      for (k0 = 0; k0 < nz; k0++)
        {
          /* Pointer to data */

          loc0 = ReadMeshPtr(msh, i0, j0, k0);
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          loc0 = (long)RDB[loc0];
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Reset maximum number of neighbours */

          nmax1 = 0;

          /* Loop over neighbour cells */

          for (i = i0 - 1; i < i0 + 2; i++)
          for (j = j0 - 1; j < j0 + 2; j++)
          for (k = k0 - 1; k < k0 + 2; k++)
            if ((i != i0) || (j != j0) || (k != k0))
              {
                /* Check type */

                if ((type == MESH_TYPE_CARTESIAN) ||
                    (type == MESH_TYPE_ORTHOGONAL))
                  {
                    /* Square mesh, check boundaries */

                    if (((i > -1) && (i < nx) && (j > -1) && (j < ny) &&
                         (k > -1) && (k < nz)) &&
                        (((i != i0) && (j == j0) && (k == k0)) ||
                         ((i == i0) && (j != j0) && (k == k0)) ||
                         ((i == i0) && (j == j0) && (k != k0))))
                      {
                        /* Pointer to data */

                        loc1 = ReadMeshPtr(msh, i, j, k);
                        CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                        loc1 = (long)RDB[loc1];
                        CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                        /* Create current */

                        ptr = NewItem(loc0 + RMX_CELL_PTR_BOUNDS,
                                      RMX_CELL_BOUND_BLOCK_SIZE);

                        /* Put pointer and reset indexes */

                        WDB[ptr + RMX_CELL_BOUND_PTR_CELL] = (double)loc1;
                        WDB[ptr + RMX_CELL_BOUND_FWD_IDX] = -1.0;
                        WDB[ptr + RMX_CELL_BOUND_ADJ_IDX] = -1.0;

                        /* Add to maximum number of neighbours */

                        nmax1++;
                      }
                  }
                else if ((type == MESH_TYPE_HEXX) ||
                         (type == MESH_TYPE_HEXY))
                  {
                    /* Hexagonal mesh, check boundaries */

                    if ((i > -1) && (i < nx) && (j > -1) && (j < ny) &&
                        (k > -1) && (k < nz))
                      {
                        /* Pointer to data */

                        loc1 = ReadMeshPtr(msh, i, j, k);
                        CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                        loc1 = (long)RDB[loc1];
                        CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                        /* Put face (W-SW-SE-E-NE-NW-B-T) */

                        if ((i - i0 == -1) && (j - j0 == 0) && (k - k0 == 0))
                          tgt = 0;
                        else if ((i - i0 == 0) && (j - j0 == -1) &&
                                 (k - k0 == 0))
                          tgt = 1;
                        else if ((i - i0 == 1) && (j - j0 == -1) &&
                                 (k - k0 == 0))
                          tgt = 2;
                        else if ((i - i0 == 1) && (j - j0 == 0) &&
                                 (k - k0 == 0))
                          tgt = 3;
                        else if ((i - i0 == 0) && (j - j0 == 1) &&
                                 (k - k0 == 0))
                          tgt = 4;
                        else if ((i - i0 == -1) && (j - j0 == 1) &&
                                 (k - k0 == 0))
                          tgt = 5;
                        else if ((i - i0 == 0) && (j - j0 == 0) &&
                                 (k - k0 == -1))
                          tgt = 6;
                        else if ((i - i0 == 0) && (j - j0 == 0) &&
                                 (k - k0 == 1))
                          tgt = 7;
                        else
                          tgt = -1;

                        /* Check */

                        if (tgt > -1)
                          {
                            /* Create current */

                            ptr = NewItem(loc0 + RMX_CELL_PTR_BOUNDS,
                                          RMX_CELL_BOUND_BLOCK_SIZE);

                            /* Put pointer and reset indexes */

                            WDB[ptr + RMX_CELL_BOUND_PTR_CELL] = (double)loc1;
                            WDB[ptr + RMX_CELL_BOUND_FWD_IDX] = -1.0;
                            WDB[ptr + RMX_CELL_BOUND_ADJ_IDX] = -1.0;

                            /* Add to maximum number of neighbours */

                            nmax1++;
                          }
                      }
                  }
                else if ((type == MESH_TYPE_CYLINDRICAL) ||
                         (type == MESH_TYPE_ICYL))
                  {
                    /* Cylindrical mesh, check boundaries */

                    if (((i > -1) && (i < nx) && (k > -1) && (k < nz)) &&
                        (((i != i0) && (j == j0) && (k == k0)) ||
                         ((i == i0) && (j != j0) && (k == k0)) ||
                         ((i == i0) && (j == j0) && (k != k0))))
                      {
                        /* Check angular bin and get pointer to data */
                        /* NOTE: uloimmat binit kytketään toisiinsa  */
                        /* vain jos binejä on enemmän kuin yksi. */

                        if ((ny == 1) && (i == i0) && (k == k0))
                          continue;
                        else if (j == -1)
                          loc1 = ReadMeshPtr(msh, i, ny - 1, k);
                        else if (j == ny)
                          loc1 = ReadMeshPtr(msh, i, 0, k);
                        else
                          loc1 = ReadMeshPtr(msh, i, j, k);

                        CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                        loc1 = (long)RDB[loc1];
                        CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                        /* Create current */

                        ptr = NewItem(loc0 + RMX_CELL_PTR_BOUNDS,
                                      RMX_CELL_BOUND_BLOCK_SIZE);

                        /* Put pointer and reset indexes */

                        WDB[ptr + RMX_CELL_BOUND_PTR_CELL] = (double)loc1;
                        WDB[ptr + RMX_CELL_BOUND_FWD_IDX] = -1.0;
                        WDB[ptr + RMX_CELL_BOUND_ADJ_IDX] = -1.0;

                        /* Add to maximum number of neighbours */

                        nmax1++;
                      }
                  }
                else
                  Die(FUNCTION_NAME, "Unsupported mesh type %ld", type);
              }

          /* Put matrix size */

          if (nmax1 == 0)
            Die(FUNCTION_NAME, "No neighbours");
          else if (nmax1 > 26)
            Die(FUNCTION_NAME, "WTF? %ld", nmax1);
          else
            WDB[loc0 + RMX_CELL_MTX_SIZE] = (double)nmax1;

          /* Calculate relative volumes */

          if ((type == MESH_TYPE_CARTESIAN) ||
              (type == MESH_TYPE_HEXX) ||
              (type == MESH_TYPE_HEXY))
            {
              /* Calculate volume */

              WDB[loc0 + RMX_CELL_RVOL] = 1.0/((double)(nx*ny*nz));
            }
          else if (type == MESH_TYPE_ICYL)
            {
              /* Pointer to radial boundaries */

              ptr = (long)RDB[msh + MESH_ORTHO_PTR_XLIM];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get bounds */

              r0 = RDB[ptr + i0];
              r1 = RDB[ptr + i0 + 1];

              CheckValue(FUNCTION_NAME, "r0", "", r0, 0.0, INFTY);
              CheckValue(FUNCTION_NAME, "r1", "", r1, r1, INFTY);

              /* Pointer to sector boundaries */

              ptr = (long)RDB[msh + MESH_ORTHO_PTR_YLIM];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get bounds */

              a0 = RDB[ptr + j0];
              a1 = RDB[ptr + j0 + 1];

              CheckValue(FUNCTION_NAME, "a0", "", a0, -2.0*PI, 2.0*PI);
              CheckValue(FUNCTION_NAME, "a1", "", a1, a0, 2.0*PI);

              /* Pointer to axial boundaries */

              ptr = (long)RDB[msh + MESH_ORTHO_PTR_ZLIM];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get bounds */

              z0 = RDB[ptr + k0];
              z1 = RDB[ptr + k0 + 1];

              CheckValue(FUNCTION_NAME, "z0", "", z0, -INFTY, INFTY);
              CheckValue(FUNCTION_NAME, "z1", "", z1, z1, INFTY);

              /* Calculate volume */

              WDB[loc0 + RMX_CELL_RVOL] =
                ((r1*r1 - r0*r0)*(z1 - z0)*(a1 - a0)/(2.0*PI))/
                (max0*max0*(max2 - min2));
            }
          else if (type == MESH_TYPE_ORTHOGONAL)
            {
              /* Pointer to x-boundaries */

              ptr = (long)RDB[msh + MESH_ORTHO_PTR_XLIM];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get bounds */

              x0 = RDB[ptr + i0];
              x1 = RDB[ptr + i0 + 1];

              CheckValue(FUNCTION_NAME, "x0", "", x0, -INFTY, INFTY);
              CheckValue(FUNCTION_NAME, "x1", "", x1, x1, INFTY);

              /* Pointer to y-boundaries */

              ptr = (long)RDB[msh + MESH_ORTHO_PTR_YLIM];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get bounds */

              y0 = RDB[ptr + j0];
              y1 = RDB[ptr + j0 + 1];

              CheckValue(FUNCTION_NAME, "y0", "", y0, -INFTY, INFTY);
              CheckValue(FUNCTION_NAME, "y1", "", y1, y1, INFTY);

              /* Pointer to z-boundaries */

              ptr = (long)RDB[msh + MESH_ORTHO_PTR_ZLIM];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get bounds */

              z0 = RDB[ptr + k0];
              z1 = RDB[ptr + k0 + 1];

              CheckValue(FUNCTION_NAME, "z0", "", z0, -INFTY, INFTY);
              CheckValue(FUNCTION_NAME, "z1", "", z1, z1, INFTY);

              /* Calculate volume */

              WDB[loc0 + RMX_CELL_RVOL] = (x1 - x0)*(y1 - y0)*(z1 - z0)/
                (max0 - min0)/(max1 - min1)/(max2 - min2);
            }
          else if (type == MESH_TYPE_CYLINDRICAL)
            {
              /* Calculate radii */

              r0 = ((double)i0)/((double)nx)*(max0 - min0) + min0;
              r1 = ((double)(i0 + 1))/((double)nx)*(max0 - min0) + min0;

              /* Calculate volume */

              WDB[loc0 + RMX_CELL_RVOL] =
                (r1*r1 - r0*r0)/((double)(ny*nz))/(max0*max0);
            }
          else
            Die(FUNCTION_NAME, "Unsupported mesh type %ld", type);
        }

      /* Put neighbour indexes */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Maximum number of neighbours */

          nmax1 = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax1", "", nmax1, 1, 10000);

          /* Loop over boundaries */

          loc1 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
          for (j = 0; j < nmax1; j++)
            {
              /* Pointer to cell */

              ptr = (long)RDB[loc1 + RMX_CELL_BOUND_PTR_CELL];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Maximum number of neighbours */

              nmax2 = (long)RDB[ptr + RMX_CELL_MTX_SIZE];
              CheckValue(FUNCTION_NAME, "nmax2", "", nmax2, 1, 10000);

              /* Loop over boundaries */

              loc2 = (long)RDB[ptr + RMX_CELL_PTR_BOUNDS];
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
                Die(FUNCTION_NAME, "Neighbour not found");

              /* Put indexes */

              WDB[loc1 + RMX_CELL_BOUND_FWD_IDX] = (double)i;
              WDB[loc2 + RMX_CELL_BOUND_ADJ_IDX] = (double)j;

              /* Next */

              loc1 = NextItem(loc1);
            }

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Close list */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      CloseList(loc0);

      /* Next */

      rmx = NextItem(rmx);
    }

  /***************************************************************************/

  /***** Additional stuff ****************************************************/

  /* Loop over response matrixes */

  rmx = (long)RDB[DATA_PTR_RMX0];
  while (rmx > VALID_PTR)
    {
      /* Check wwg-mode */

      if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_WWG)
        {
          /* Allocate memory for event buffer */

          ptr = AllocPrivateData(2*MAX_RMX_BUFF, PRIVA_ARRAY);
          WDB[rmx + RMX_PTR_EVENT_BUFF] = (double)ptr;

          /* Check that OMP load balancing is not in use */

          if ((long)RDB[DATA_COMMON_QUE_LIM] != 0)
            Error(rmx, "OMP load balancing not allowed in mode %ld",
                  RMX_MODE_WWG);
        }

      /* Get number of energy groups */

      ng = (long)RDB[rmx + RMX_NG];
      CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

      /* Loop over cells and allocate memory for data */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get size */

          nmax1 = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];

          /* Allocate memory for results and coefficients */

          RecycleRMXData(rmx, loc0, ng, nmax1, (long)RDB[rmx + RMX_FROM_FILE]);

          /* Allocate memory for importances (NOTE: noita ei voi */
          /* johtuen siitä miten importancen kopiointi jaetulle  */
          /* mesh cellille on toteutettu adjustrmx.c:ssä. */

          if ((long)RDB[rmx + RMX_FROM_FILE] == NO)
            {
              ptr = ReallocMem(DATA_ARRAY, ng);
              WDB[loc0 + RMX_CELL_IMP_SRC] = (double)ptr;

              ptr = ReallocMem(DATA_ARRAY, ng);
              WDB[loc0 + RMX_CELL_IMP_SRC_KEEP] = (double)ptr;

              ptr = ReallocMem(DATA_ARRAY, ng);
              WDB[loc0 + RMX_CELL_IMP_CURR] = (double)ptr;

              /* NOTE: tätä ei käytetä, mutta arvot kirjoitetaan */
              /*       fileen niin roikkukoon vielä mukana.      */

              ptr = ReallocMem(DATA_ARRAY, ng);
              WDB[loc0 + RMX_CELL_IMP_CURR_KEEP] = (double)ptr;
            }

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Check that type is set */

      if ((long)RDB[rmx + RMX_PARTICLE_TYPE] == 0)
        Die(FUNCTION_NAME, "Particle type not set");

      /* Next */

      rmx = NextItem(rmx);
    }

  /***************************************************************************/

  /***** Change mesh type for adaptive iteration *****************************/

  /* Loop over weight window structures */

  loc0 = (long)RDB[DATA_PTR_WWD0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over iteratinons */

      loc1 = (long)RDB[loc0 + WWD_PTR_ITER];
      while (loc1 > VALID_PTR)
        {
          /* Get pointer to response matrix */

          rmx = (long)RDB[loc1 + WWD_ITER_PTR_RMX];
          CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

          /* Put particle type */

          if ((long)RDB[loc0 + WWD_PARTICLE_TYPE] == 0)
            WDB[loc0 + WWD_PARTICLE_TYPE] = RDB[rmx + RMX_PARTICLE_TYPE];
          else if (RDB[loc0 + WWD_PARTICLE_TYPE] !=
                   RDB[rmx + RMX_PARTICLE_TYPE])
            Error(loc0, "Multiple particle types not allowed in iteration");

          /* Get pointer to mesh */

          msh = (long)RDB[rmx + RMX_PTR_MESH];
          CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

          /* Check adaptive iterations */

          if (((long)RDB[loc1 + WWD_ITER_TYPE] == WWIN_ITER_GEO) ||
              ((long)RDB[loc1 + WWD_ITER_TYPE] == WWIN_ITER_TRK) ||
              ((long)RDB[loc1 + WWD_ITER_TYPE] == WWIN_ITER_CON))
            {
              /* Change type */

              if (((long)RDB[msh + MESH_TYPE] != MESH_TYPE_CARTESIAN) &&
                  ((long)RDB[msh + MESH_TYPE] != MESH_TYPE_ADAPTIVE))
                Error(loc0, "Invalid mesh type for adaptive iterations");
              else
                WDB[msh + MESH_TYPE] = (double)MESH_TYPE_ADAPTIVE;
            }

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Additional stuff for convergence acceleration ***********************/

  if ((long)RDB[DATA_RMX_CONVG_ACC] == YES)
    {
      /* Loop over response matrixes */

      rmx = (long)RDB[DATA_PTR_RMX0];
      while (rmx > VALID_PTR)
        {
          /* Get sub-mesh size */

          nx = (long)RDB[rmx + RMX_CONVG_SUBMESH_NX];
          ny = (long)RDB[rmx + RMX_CONVG_SUBMESH_NY];
          nz = (long)RDB[rmx + RMX_CONVG_SUBMESH_NZ];

          /* Pointer to global mesh */

          msh = (long)RDB[rmx + RMX_PTR_MESH];
          CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

          /* Check if 2D geometry */

          if ((long)RDB[DATA_GEOM_DIM] == 2)
            {
              /* Check mesh and submesh */

              if ((long)RDB[msh + MESH_N2] > 1)
                Error(rmx, "Axial mesh size must be 1 in 2D geometry");
              else if ((nz > 1))
                Error(rmx, "Axial submesh size must be 1 in 2D geometry");
            }

          /* Loop over data */

          loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          while (loc0 > VALID_PTR)
            {
              /* Number of neighbours */

              nmax1 = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];

              /* Allocate memory for form factors */

              ptr = AllocPrivateData(nmax1*nx*ny*nz, RES2_ARRAY);
              WDB[loc0 + RMX_CELL_MC_SURF_FF] = (double)ptr;

              ptr = AllocPrivateData(nx*ny*nz, RES2_ARRAY);
              WDB[loc0 + RMX_CELL_MC_SRC_FF] = (double)ptr;

              /* Allocate memory for calculated coefficients */

              ptr = ReallocMem(DATA_ARRAY, nmax1*nx*ny*nz);
              WDB[loc0 + RMX_CELL_COEF_SURF_FF] = (double)ptr;

              ptr = ReallocMem(DATA_ARRAY, nmax1*nx*ny*nz);
              WDB[loc0 + RMX_CELL_COEF_SRC_FF] = (double)ptr;

              /* Get cell index */

              i = (long)RDB[loc0 + RMX_CELL_IDX];
              CheckValue(FUNCTION_NAME, "i", "", i, 0, 1E+9);

              /* Get boundaries */

              MeshCellBounds(msh, i, &lims[0], &lims[1], &lims[2], &lims[3],
                             &lims[4], &lims[5]);

              /* Check type */

              if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXX)
                {
                  /* Check number of cells */

                  if (nx != ny)
                    Die(FUNCTION_NAME, "Error in submesh definition");

                  /* Calculate origin */

                  x0 = 0.5*(lims[0] + lims[1]);
                  y0 = 0.5*(lims[2] + lims[3]);

                  /* Get pitch */

                  p = RDB[rmx + RMX_CONVG_SUBMESH_P];
                  CheckValue(FUNCTION_NAME, "p", "", p, ZERO, INFTY);

                  /* Put values */

                  lims[0] = x0;
                  lims[2] = y0;
                  lims[1] = p;
                  lims[3] = 0.0;

                  /* Create sub-mesh */

                  if (nx > 1)
                    ptr = CreateMesh(MESH_TYPE_HEXY, MESH_CONTENT_DAT, -1,
                                     nx, ny, nz, lims, -1);
                  else
                    ptr = CreateMesh(MESH_TYPE_HEXX, MESH_CONTENT_DAT, -1,
                                     nx, ny, nz, lims, -1);


                  WDB[loc0 + RMX_CELL_PTR_SUBMESH] = (double)ptr;
                }
              else if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXY)
                {
                  /* Check number of cells */

                  if (nx != ny)
                    Die(FUNCTION_NAME, "Error in submesh definition");

                  /* Calculate origin */

                  x0 = 0.5*(lims[0] + lims[1]);
                  y0 = 0.5*(lims[2] + lims[3]);

                  /* Calculate pitch */

                  p = RDB[rmx + RMX_CONVG_SUBMESH_P];

                  /* Put values */

                  lims[0] = x0;
                  lims[2] = y0;
                  lims[1] = p;
                  lims[3] = 0.0;

                  /* Create sub-mesh */

                  if (nx > 1)
                    ptr = CreateMesh(MESH_TYPE_HEXX, MESH_CONTENT_DAT, -1,
                                     nx, ny, nz, lims, -1);
                  else
                    ptr = CreateMesh(MESH_TYPE_HEXY, MESH_CONTENT_DAT, -1,
                                     nx, ny, nz, lims, -1);

                  WDB[loc0 + RMX_CELL_PTR_SUBMESH] = (double)ptr;
                }
              else
                {
                  /* Adjust */

                  lims[0] = lims[0] - 1E-6;
                  lims[1] = lims[1] + 1E-6;
                  lims[2] = lims[2] - 1E-6;
                  lims[3] = lims[3] + 1E-6;
                  lims[4] = lims[4] - 1E-6;
                  lims[5] = lims[5] + 1E-6;

                  /* Create sub-mesh */

                  ptr = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_DAT, -1,
                                   nx, ny, nz, lims, -1);
                  WDB[loc0 + RMX_CELL_PTR_SUBMESH] = (double)ptr;
                }

              /* Create sampling structure */

              for (i = 0; i < nx*ny*nz; i++)
                NewItem(loc0 + RMX_CELL_PTR_SUBMESH_SAMPLE,
                        RMX_SUBMESH_SAMPLE_BLOCK_SIZE);

              /* Next */

              loc0 = NextItem(loc0);
            }

          /* Next */

          rmx = NextItem(rmx);
        }
    }

  /* Memory count */

  WDB[DATA_TOT_RMX_BYTES] = RDB[DATA_TOT_RMX_BYTES] + (double)MemCount();

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
