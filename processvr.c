/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processvr.c                                    */
/*                                                                           */
/* Created:       2011/05/15 (JLe)                                           */
/* Last modified: 2019/10/08 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Processes some variance reduction stuff                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessVR:"

/*****************************************************************************/

void ProcessVR()
{
  long wwd, msh, src, loc0, loc1, ptr, nx, ny, nz, lat, ne, n, type;
  long rmx;
  double lims[6], p, q, c, x, y, z, E, min, max, sum, mean;

  /***************************************************************************/

  /****** Uniform fission source method **************************************/

  if (((long)RDB[DATA_UFS_MODE] != UFS_MODE_NONE) &&
      ((long)RDB[DATA_UFS_MODE] != UFS_MODE_RMX))
    {
      /* Check pointer to lattice */

      if (((long)RDB[DATA_UFS_PTR_LAT]) > VALID_PTR)
        {
          /* Find lattice */

          lat = (long)RDB[DATA_PTR_L0];
          while (lat > VALID_PTR)
            {
              /* Compare names */

              if (CompareStr(lat + LAT_PTR_NAME, DATA_UFS_PTR_LAT))
                break;

              /* Next lattice */

              lat = NextItem(lat);
            }

          /* Check */

          if (lat < VALID_PTR)
            Error(0, "Lattice \"%s\" used for UFS not found in geometry",
                  GetText(DATA_UFS_PTR_LAT));
          else
            WDB[DATA_UFS_PTR_LAT] = (double)lat;

          /* Check type */

          if (((long)RDB[lat + LAT_TYPE] != LAT_TYPE_S) &&
              ((long)RDB[lat + LAT_TYPE] != LAT_TYPE_HX) &&
              ((long)RDB[lat + LAT_TYPE] != LAT_TYPE_HY))
            Error(0, "Lattice \"%s\" is wrong type for UFS",
                  GetText(lat + LAT_PTR_NAME));

          /* Put x- and y-bins */

          WDB[DATA_UFS_NX] = RDB[lat + LAT_NX];
          WDB[DATA_UFS_NY] = RDB[lat + LAT_NY];
        }

      /* Set boundaries if not set in input */

      if (RDB[DATA_UFS_XMIN] == -INFTY)
        WDB[DATA_UFS_XMIN] = RDB[DATA_GEOM_MINX];

      if (RDB[DATA_UFS_XMAX] == INFTY)
        WDB[DATA_UFS_XMAX] = RDB[DATA_GEOM_MAXX];

      if (RDB[DATA_UFS_YMIN] == -INFTY)
        WDB[DATA_UFS_YMIN] = RDB[DATA_GEOM_MINY];

      if (RDB[DATA_UFS_YMAX] == INFTY)
        WDB[DATA_UFS_YMAX] = RDB[DATA_GEOM_MAXY];

      if (RDB[DATA_UFS_ZMIN] == -INFTY)
        WDB[DATA_UFS_ZMIN] = RDB[DATA_GEOM_MINZ];

      if (RDB[DATA_UFS_ZMAX] == INFTY)
        WDB[DATA_UFS_ZMAX] = RDB[DATA_GEOM_MAXZ];

      /* Get mesh size */

      nx = (long)RDB[DATA_UFS_NX];
      ny = (long)RDB[DATA_UFS_NY];
      nz = (long)RDB[DATA_UFS_NZ];

      /* Check values */

      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 10000);
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 10000);
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 10000);

      /* Put mesh variables */

      lims[0] = RDB[DATA_UFS_XMIN];
      lims[1] = RDB[DATA_UFS_XMAX];
      lims[2] = RDB[DATA_UFS_YMIN];
      lims[3] = RDB[DATA_UFS_YMAX];
      lims[4] = RDB[DATA_UFS_ZMIN];
      lims[5] = RDB[DATA_UFS_ZMAX];

      /* Create mesh */

      ptr = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_RES, -1,
                       nx, ny, nz, lims, -1);

      /* Put pointer */

      WDB[DATA_UFS_PTR_SRC_MESH] = (double)ptr;

      /* Allocate memory for factors */

      ptr = ReallocMem(DATA_ARRAY, nx*ny*nz);
      WDB[DATA_UFS_PTR_FACTORS] = (double)ptr;
    }

  /***************************************************************************/

  /***** Geometry importances ************************************************/

  if ((long)RDB[DATA_USE_GEOM_IMP] == YES)
    {
      /* Check delta-tracking */

      if ((((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES) &&
           (RDB[DATA_DT_PTHRESH] < 1.0)) ||
          (((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES) &&
           (RDB[DATA_DT_NTHRESH] < 1.0)))
        if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO)
          Error(0, "Importances not allowed with delta-tracking");

      /* Put unit cell importances */

      ptr = (long)RDB[DATA_PTR_C0];
      while (ptr > VALID_PTR)
        {
          /* Reset */

          WDB[ptr + CELL_IMP] = 1.0;

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Loop over data */

      loc0 = (long)RDB[DATA_PTR_IMP0];
      while (loc0 > VALID_PTR)
        {
          /* Check cell pointer */

          if ((long)RDB[loc0 + IMP_PTR_CELL] > VALID_PTR)
            {
              /* Find cell */

              ptr = (long)RDB[DATA_PTR_C0];
              if ((ptr = SeekListStr(ptr, CELL_PTR_NAME,
                                     GetText(loc0 + IMP_PTR_CELL))) < VALID_PTR)
                Error(loc0, "Cell %s not found",
                      GetText(loc0 + IMP_PTR_CELL));
              else
                WDB[ptr + CELL_IMP] = RDB[loc0 + IMP_I];
            }

          /* Next */

          loc0 = NextItem(loc0);
        }
    }

  /***************************************************************************/

  /***** Process weight windows **********************************************/

  /* Set default weight window bounds (set to -1 in initdata.c) */

  if (RDB[DATA_WWD_LOWER_BOUND] < 0.0)
    {
      WDB[DATA_WWD_LOWER_BOUND] = 0.5;
      WDB[DATA_WWD_UPPER_BOUND] = 2.0;
    }

  /* Loop over weight windows */

  wwd = (long)RDB[DATA_PTR_WWD0];
  while (wwd > VALID_PTR)
    {
      /* Get type */

      type = (long)RDB[wwd + WWD_TYPE];

      /* Check for iterations */

      if ((loc0 = (long)RDB[wwd + WWD_PTR_ITER]) > VALID_PTR)
        {
          /* Check number of tracks in adaptive mode (NOTE: Tää ei */
          /* välttämättä ole edes tarpeellinen) */

          if ((long)RDB[loc0 + WWD_ITER_TYPE] == WWIN_ITER_TRK)
            if (RDB[DATA_SRC_POP] < RDB[loc0 + WWD_ITER_SPLIT_TRACKS])
              Error(wwd, "Maximum number of tracks is %ld\n",
                    (long)RDB[DATA_SRC_POP]);

          /* Put total number of cycles */

          WDB[DATA_TOT_VR_ITER] = (double)ListSize(loc0);

          /* Loop over iterations */

          while (loc0 > VALID_PTR)
            {
              /* Find response matrix definition */

              rmx = (long)RDB[DATA_PTR_RMX0];
              while (rmx > VALID_PTR)
                {
                  /* Compare */

                  if (CompareStr(loc0 + WWD_ITER_PTR_RMX, rmx + RMX_PTR_NAME))
                    {
                      /* Put pointer */

                      WDB[loc0 + WWD_ITER_PTR_RMX] = (double)rmx;

                      /* Break loop */

                      break;
                    }

                  /* Next */

                  rmx = NextItem(rmx);
                }

              /* Check if found */

              if (rmx < VALID_PTR)
                Error(wwd, "RMX %s not defined",
                      GetText(loc0 + WWD_ITER_PTR_RMX));
              /* Next */

              loc0 = NextItem(loc0);
            }

          /* Check if first mesh if read from file */

          if (RDB[wwd + WWD_TYPE] == 0.0)
            {
              /* Put type */

              WDB[wwd + WWD_TYPE] = (double)WWD_MESH_TYPE_ITER;

              /* Put scaling and normalization factors */

              WDB[wwd + WWD_MULT] = -1.0;
              WDB[wwd + WWD_POW] = -1.0;
              WDB[wwd + WWD_NORM_FACT] = 1.0;
            }

          /* Check that no more structures are defined */

          if ((ptr = NextItem(wwd)) > VALID_PTR)
            Error(ptr, "Multiple weight window definitions in iteration mode");

          /* Break loop */

          if ((long)RDB[wwd + WWD_TYPE] == WWD_MESH_TYPE_ITER)
            break;
        }

      /* Check file */

      if ((long)RDB[wwd + WWD_PTR_FNAME] > VALID_PTR)
        {
          /* Read file */

          ReadWWDMesh(wwd);

          /* Apply smoothing */

          AvgRMX(wwd);

          /* Number of energy group is zero if data was skipped */

          if (((long)RDB[wwd + WWD_NE] == 0) && (type == WWD_MESH_TYPE_MCNP))
            {
              /* Copy pointer */

              ptr = wwd;

              /* Pointer to next */

              wwd = NextItem(wwd);

              /* Remove */

              RemoveItem(ptr);

              /* Cycle loop */

              continue;
            }
        }
      else if (type != WWD_MESH_TYPE_MCNP)
        Error(wwd, "Weight window mesh file or iteration must be defined");

      /***********************************************************************/

      /***** Source biasing **************************************************/

      /* Check source biasing */

      if (((long)RDB[wwd + WWD_SRC_BIAS] == YES) &&
          (type != WWD_MESH_TYPE_MCNP))
        {
          Die(FUNCTION_NAME, "Tähän tarvii sampling listan");

          /* Check if separate bounds */

          if ((long)RDB[wwd + WWD_BOUNDS_TYPE] == WW_MESH_BOUNDS_SEP)
            Error(wwd, "Segment-importances not allowed with source biasing");

          /* Loop over sources and check illegal types */

          src = (long)RDB[DATA_PTR_SRC0];
          while (src > VALID_PTR)
            {
              /* Check point source */

              if ((RDB[src + SRC_X0] > -INFTY) &&
                  (RDB[src + SRC_Y0] > -INFTY) &&
                  (RDB[src + SRC_Z0] > -INFTY))
                Error(wwd, "Source biasing cannot be used with point sources");

              /* Check surface source */

              if ((long)RDB[src + SRC_PTR_SURF] > VALID_PTR)
                Error(wwd,
                      "Source biasing cannot be used with surface sources");

              /* Check plasma source */

              if ((long)RDB[src + SRC_READ_FILE_TYPE] ==
                  SRC_FILE_TYPE_FUSION_PLASMA)
                Error(wwd,
                      "Source biasing cannot be used with plasma sources");

              /* Check source files */

              if ((long)RDB[src + SRC_READ_PTR_FILE] > VALID_PTR)
                Error(wwd, "Source biasing cannot be used with source files");

              /* Next */

              src = NextItem(src);
            }

          /* Get number of energy groups */

          ne = (long)RDB[wwd + WWD_NE];
          CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000);

          /* Allocate memory for maximum importances */

          loc1 = ReallocMem(DATA_ARRAY, ne);
          WDB[wwd + WWD_PTR_SUM_SRC_IMP] = (double)loc1;

          /* Loop over energy groups */

          for (n = 0; n < ne; n++)
            {
              /* Reset sum and mean */

              sum = 0.0;
              mean = 0.0;

              /* Loop over data */

              loc0 = (long)RDB[wwd + WWD_PTR_MESH_DATA];
              while (loc0 > VALID_PTR)
                {
                  /* Pointer to current importances */

                  ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
                  CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

                  /* Add ot sum and mean */

                  if (RDB[ptr + n] > 0.0)
                    {
                      sum = sum + RDB[ptr + n];
                      mean = mean + 1.0;
                    }

                  /* Next */

                  loc0 = NextItem(loc0);
                }

              /* Check count */

              if (mean > 0.0)
                {
                  /* Calculate mean */

                  mean = sum/mean;

                  /* Renormalize importances */

                  loc0 = (long)RDB[wwd + WWD_PTR_MESH_DATA];
                  while (loc0 > VALID_PTR)
                    {
                      /* Source importance */

                      ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      WDB[ptr + n] = RDB[ptr + n]/mean;

                      /* Current importance */

                      ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      WDB[ptr + n] = RDB[ptr + n]/mean;

                      /* Next */

                      loc0 = NextItem(loc0);
                    }

                  /* Put total */

                  WDB[loc1 + n] = sum/mean;
                }
            }
        }

      /***********************************************************************/

      /***** Normalization ***************************************************/

      if (RDB[wwd + WWD_NORM_FACT] > 0.0)
        {
          /* Get point */

          x = RDB[wwd + WWD_NORM_X];
          y = RDB[wwd + WWD_NORM_Y];
          z = RDB[wwd + WWD_NORM_Z];
          E = RDB[wwd + WWD_NORM_E];

          /* Avoid compiler warning */

          n = -1;

          /* Pointer to energy array */

          if ((ptr = (long)RDB[wwd + WWD_PTR_ERG]) < VALID_PTR)
            n = 0;
          else
            {
              /* Get number of energy groups */

              ne = (long)RDB[wwd + WWD_NE];
              CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000);

              /* Get energy group (maksimi ei oo mukana --> ongelmia) */

              if ((n = SearchArray(&RDB[ptr], E, ne + 1)) < 0)
                Error(wwd, "Normalization energy outside group boundaries");
            }

          /* Check type */

          if (type == WWD_MESH_TYPE_MCNP)
            {
              /* Pointer to mesh */

              msh = (long)RDB[wwd + WWD_PTR_MESH];
              CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

              /* Pointer to data */

              if ((ptr = MeshPtr(msh, x, y, z)) < VALID_PTR)
                Error(wwd, "Normalization point is outside mesh");

              ptr = (long)RDB[ptr];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Pointer to importance vector */

              ptr = (long)RDB[ptr + WWD_MESH_PTR_IMP];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get value */

              p = RDB[ptr + n];
            }
          else
            {
              /* Pointer to mesh */

              msh = (long)RDB[wwd + WWD_PTR_MESH];
              CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

              /* Pointer to data */

              if ((ptr = MeshPtr(msh, x, y, z)) < VALID_PTR)
                Error(wwd, "Normalization point is outside mesh");

              ptr = (long)RDB[ptr];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Pointer to importance vector */

              ptr = (long)RDB[ptr + RMX_CELL_IMP_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get value */

              p = RDB[ptr + n];
            }

          /* Get coefficient and exponential */

          if (((c = RDB[wwd + WWD_MULT]) > 0.0) &&
              ((q = RDB[wwd + WWD_POW]) > 0.0))
            p = c*pow(p, q);

          /* Calculate normalization factor */

          WDB[wwd + WWD_NORM_FACT] = RDB[wwd + WWD_NORM_FACT]/p;
        }
      else
        {
          /* Relative factor given (negative entry) */

          WDB[wwd + WWD_NORM_FACT] = -RDB[wwd + WWD_NORM_FACT];
        }

      /* Reset minimum and maximum */

      min = INFTY;
      max = -INFTY;

      /* Check type */

      if (type == WWD_MESH_TYPE_MCNP)
        {
          /* Loop over mesh */

          loc1 = (long)RDB[wwd + WWD_PTR_MESH_DATA];
          while (loc1 > VALID_PTR)
            {
              /* Pointer to values */

              ptr = (long)RDB[loc1 + WWD_MESH_PTR_IMP];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get number of energy groups */

              ne = (long)RDB[wwd + WWD_NE];
              CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000);

              /* Loop over data */

              for (n = 0; n < ne; n++)
                {
                  /* Compare to limits */

                  if ((RDB[wwd + WWD_LIM_MIN] > 0.0) &&
                      (RDB[ptr] > 0.0) &&
                      (RDB[ptr] < RDB[wwd + WWD_LIM_MIN]))
                    WDB[ptr] = RDB[wwd + WWD_LIM_MIN];
                  else if ((RDB[wwd + WWD_LIM_MAX] > 0.0) &&
                           (RDB[ptr] > 0.0) &&
                           (RDB[ptr] > RDB[wwd + WWD_LIM_MAX]))
                    WDB[ptr] = RDB[wwd + WWD_LIM_MAX];

                  /* Get value */

                  p = RDB[ptr++];

                  /* Get coefficient and exponential */

                  if (((c = RDB[wwd + WWD_MULT]) > 0.0) &&
                      ((q = RDB[wwd + WWD_POW]) > 0.0))
                    p = c*pow(p, q);

                  /* Multiply by normalization factor */

                   if (RDB[wwd + WWD_NORM_FACT] > 0.0)
                    p = p*RDB[wwd + WWD_NORM_FACT];

                  /* Compare to minimum and maximum */

                  if (p > max)
                    max = p;
                  if ((p > 0.0) && (p < min))
                    min = p;
                }

              /* Next */

              loc1 = NextItem(loc1);
            }
              }
      else if (type == WWD_MESH_TYPE_SSS)
        {
          /* Get number of energy groups */

          ne = (long)RDB[wwd + WWD_NE];
          CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000);

          /* Loop over data */

          loc0 = (long)RDB[wwd + WWD_PTR_MESH_DATA];
          while (loc0 > VALID_PTR)
            {
              /* Pointer to current importances */

              ptr = (long)RDB[loc0 + RMX_CELL_IMP_CURR];
              CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

              /* Loop over energy groups */

              for (n = 0; n < ne; n++)
                {
                  /* Compare to limits */

                  if ((RDB[wwd + WWD_LIM_MIN] > 0.0) &&
                      (RDB[ptr] > 0.0) &&
                      (RDB[ptr] < RDB[wwd + WWD_LIM_MIN]))
                    WDB[ptr] = RDB[wwd + WWD_LIM_MIN];
                  else if ((RDB[wwd + WWD_LIM_MAX] > 0.0) &&
                           (RDB[ptr] > 0.0) &&
                           (RDB[ptr] > RDB[wwd + WWD_LIM_MAX]))
                    WDB[ptr] = RDB[wwd + WWD_LIM_MAX];

                  /* Get value */

                  p = RDB[ptr++];

                  /* Get coefficient and exponential */

                  if (((c = RDB[wwd + WWD_MULT]) > 0.0) &&
                      ((q = RDB[wwd + WWD_POW]) > 0.0))
                    p = c*pow(p, q);

                  /* Multiply by normalization factor */

                   if (RDB[wwd + WWD_NORM_FACT] > 0.0)
                    p = p*RDB[wwd + WWD_NORM_FACT];

                  /* Compare to minimum and maximum */

                  if (p > max)
                    max = p;
                  if ((p > 0.0) && (p < min))
                    min = p;
                }

              /* Next */

              loc0 = NextItem(loc0);
            }
        }
      else
        Die(FUNCTION_NAME, "Unsupported mesh type");

      /* Check zero result */

      if (min == INFTY)
        Error(wwd, "No nonzero importances");

      /***********************************************************************/

      /***** Criticality source calculation **********************************/

#ifdef mmmmmmmmmmmmmmmm
      /* Check criticality source simulation */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        {
          /* Reset count and sum */

          n = 0;
          sum = 0.0;

          /* Loop over energy groups */

          loc0 = (long)RDB[wwd + WWD_PTR_EDATA];
          while (loc0 > VALID_PTR)
            {
              /* Loop over mesh data */

              msh = (long)RDB[loc0 + WWD_EDATA_PTR_MESH_DATA];
              while (msh > VALID_PTR)
                {
                  /* Add to sum */

                  sum = sum + RDB[msh + WWD_MESH_SRC_IMP];

                  /* Add to count */

                  if (RDB[msh + WWD_MESH_SRC_IMP] > 0.0)
                    n++;

                  /* Next */

                  msh = NextItem(msh);
                }

              /* Next */

              loc0 = NextItem(loc0);
            }

          /* Put normalization factor */

          if (n > 0)
            WDB[wwd + WWD_NORM_FACT] = ((double)n)/sum;
          else
            Error(wwd, "No nonzero source importances");
        }

#endif
      /***********************************************************************/

      /* Print */

      fprintf(outp, "Importance mesh:\n\n");

      if (type == WWD_MESH_TYPE_MCNP)
        fprintf(outp, " - MCNP weight window mesh");
      else if (type == WWD_MESH_TYPE_SSS)
        fprintf(outp, " - Serpent weight window mesh");
      else
        Die(FUNCTION_NAME, "Type not supported");

      /* Check if particle type is given */

      if ((long)RDB[wwd + WWD_PARTICLE_TYPE] == PARTICLE_TYPE_NEUTRON)
        fprintf(outp, " for neutrons\n");
      else if ((long)RDB[wwd + WWD_PARTICLE_TYPE] == PARTICLE_TYPE_GAMMA)
        fprintf(outp, " for photons\n");
      else
        Die(FUNCTION_NAME, "Particle type is not defined");

      if ((long)RDB[wwd + WWD_SRC_BIAS] == YES)
        fprintf(outp, " - Source biasing in use\n");

      /* Pointer to energies */

      if ((ptr = (long)RDB[wwd + WWD_PTR_ERG]) > VALID_PTR)
        {
          /* Get number of energy groups */

          ne = (long)RDB[wwd + WWD_NE];
          CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000);

          /* Print */

          if (type == WWD_MESH_TYPE_SSS)
            fprintf(outp,
                    " - %ld energy groups between %1.5E and %1.5E MeV\n", ne,
                    RDB[ptr + ENERGY_GRID_EMIN], RDB[ptr + ENERGY_GRID_EMAX]);
          else
            fprintf(outp,
                    " - %ld energy groups between %1.5E and %1.5E MeV\n", ne,
                    RDB[ptr], RDB[ptr + ne]);
         }
      else
        fprintf(outp, " - No energy binning\n");

      /* Pointer to mesh */

      msh = (long)RDB[wwd + WWD_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Check type */

      if (((long)RDB[msh + MESH_TYPE] == MESH_TYPE_CARTESIAN) ||
          ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_ORTHOGONAL))
        {
          /* Print */

          fprintf(outp, " - %ld X-bins between %1.5E and %1.5E cm\n",
                  (long)RDB[msh + MESH_N0], RDB[msh + MESH_MIN0],
                  RDB[msh + MESH_MAX0]);

          fprintf(outp, " - %ld Y-bins between %1.5E and %1.5E cm\n",
                  (long)RDB[msh + MESH_N1], RDB[msh + MESH_MIN1],
                  RDB[msh + MESH_MAX1]);

          fprintf(outp, " - %ld Z-bins between %1.5E and %1.5E cm\n",
                  (long)RDB[msh + MESH_N2], RDB[msh + MESH_MIN2],
                  RDB[msh + MESH_MAX2]);
        }
      else if (((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXX) ||
               ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXY))
        {
          /* Print */

          if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXX)
            fprintf(outp, " - X-type hexagonal mesh\n");
          else
            fprintf(outp, " - Y-type hexagonal mesh\n");

          fprintf(outp, " - Centered at x = %1.5E cm, y = %1.5E cm\n",
                  RDB[msh + MESH_MIN0], RDB[msh + MESH_MIN1]);

          fprintf(outp, " - %ld x %ld cells, %1.5E cm pitch\n",
                  (long)RDB[msh + MESH_N0], (long)RDB[msh + MESH_N1],
                  RDB[msh + MESH_MAX0]);

          fprintf(outp, " - %ld Z-bins between %1.5E and %1.5E cm\n",
                  (long)RDB[msh + MESH_N2], RDB[msh + MESH_MIN2],
                  RDB[msh + MESH_MAX2]);
        }
      else if (((long)RDB[msh + MESH_TYPE] == MESH_TYPE_CYLINDRICAL) ||
               ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_ICYL))
        {
          /* Print */

          fprintf(outp, " - %ld Radial bins between %1.5E and %1.5E cm\n",
                  (long)RDB[msh + MESH_N0], RDB[msh + MESH_MIN0],
                  RDB[msh + MESH_MAX0]);

          fprintf(outp, " - %ld Angular bins between %1.1f and %1.1f deg\n",
                  (long)RDB[msh + MESH_N1],
                  RDB[msh + MESH_MIN1]*180.0/PI,
                  RDB[msh + MESH_MAX1]*180.0/PI);

          fprintf(outp, " - %ld Axial bins between %1.5E and %1.5E cm\n",
                  (long)RDB[msh + MESH_N2], RDB[msh + MESH_MIN2],
                  RDB[msh + MESH_MAX2]);
        }
      else if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_ADAPTIVE)
        {
          ptr = (long)RDB[wwd + WWD_PTR_MESH_DATA];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          fprintf(outp, " - Adaptive mesh with %ld cells\n", ListSize(ptr));

          fprintf(outp, " - X-bins between %1.5E and %1.5E cm\n",
                  RDB[msh + MESH_MIN0], RDB[msh + MESH_MAX0]);

          fprintf(outp, " - Y-bins between %1.5E and %1.5E cm\n",
                  RDB[msh + MESH_MIN1], RDB[msh + MESH_MAX1]);

          fprintf(outp, " - Z-bins between %1.5E and %1.5E cm\n",
                  RDB[msh + MESH_MIN2], RDB[msh + MESH_MAX2]);
        }
      else
        Die(FUNCTION_NAME, "TODO: %ld", (long)RDB[msh + MESH_TYPE]);

      if (RDB[wwd + WWD_NORM_FACT] != 1.0)
        fprintf(outp, " - Renormalization factor %1.5E\n",
                RDB[wwd + WWD_NORM_FACT]);

      /* Check cut-off limits */

      if (min == RDB[wwd + WWD_LIM_MIN])
        fprintf(outp, " - Minimum cell importance: %1.5E (lim)\n", min);
      else
        fprintf(outp, " - Minimum cell importance: %1.5E\n", min);

      if (max == RDB[wwd + WWD_LIM_MAX])
        fprintf(outp, " - Maximum cell importance: %1.5E (lim)\n", max);
      else
        fprintf(outp, " - Maximum cell importance: %1.5E\n", max);

      /* Put minimum and maximum values */

      WDB[wwd + WWD_MESH_MIN] = min;
      WDB[wwd + WWD_MESH_MAX] = max;

      /* Check source biasing */

      if ((long)RDB[wwd + WWD_SRC_BIAS] == YES)
        {
          /* Check renormalization */

          if (RDB[wwd + WWD_NORM_FACT] != 1.0)
            Error(wwd, "Renormalization not allowed with source biasing");

          /* Check adjustment */

          if (((c = RDB[wwd + WWD_MULT]) > 0.0) ||
              ((q = RDB[wwd + WWD_POW]) > 0.0))
            Error(wwd, "Adjustment not allowed with source biasing");

          /* Check type */

          if (type != WWD_MESH_TYPE_SSS)
            Die(FUNCTION_NAME, "Invalid mesh type");
        }

      /* Find minimum and maximum source importance */

      min = INFTY;
      max = -INFTY;

      if (type == WWD_MESH_TYPE_SSS)
        {
          /* Get number of energy groups */

          ne = (long)RDB[wwd + WWD_NE];
          CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000);

          /* Loop over data */

          loc0 = (long)RDB[wwd + WWD_PTR_MESH_DATA];
          while (loc0 > VALID_PTR)
            {
              /* Pointer to current importances */

              ptr = (long)RDB[loc0 + RMX_CELL_IMP_SRC];
              CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

              /* Loop over energy groups */

              for (n = 0; n < ne; n++)
                {
                  /* Compare to limits */

                  if ((RDB[wwd + WWD_LIM_MIN] > 0.0) && (RDB[ptr + n] > 0.0) &&
                      (RDB[ptr + n] < RDB[wwd + WWD_LIM_MIN]))
                    WDB[ptr + n] = RDB[wwd + WWD_LIM_MIN];
                  else if ((RDB[wwd + WWD_LIM_MAX] > 0.0) &&
                           (RDB[ptr + n] > 0.0) &&
                           (RDB[ptr + n] > RDB[wwd + WWD_LIM_MAX]))
                    WDB[ptr + n] = RDB[wwd + WWD_LIM_MAX];


                 /* Get value */

                  p = RDB[ptr++];

                  /* Get coefficient and exponential */

                  if (((c = RDB[wwd + WWD_MULT]) > 0.0) &&
                      ((q = RDB[wwd + WWD_POW]) > 0.0))
                    p = c*pow(p, q);

                  /* Multiply by normalization factor */

                   if (RDB[wwd + WWD_NORM_FACT] > 0.0)
                    p = p*RDB[wwd + WWD_NORM_FACT];

                  /* Compare to minimum and maximum */

                  if (p > max)
                    max = p;
                  if ((p > 0.0) && (p < min))
                    min = p;
                }

              /* Next */

              loc0 = NextItem(loc0);
            }
        }

      /* Print */

      if (min < max)
        {
          /* Check cut-off limits */

          if (min == RDB[wwd + WWD_LIM_MIN])
            fprintf(outp, " - Minimum source importance: %1.5E (lim)\n", min);
          else
            fprintf(outp, " - Minimum source importance: %1.5E\n", min);

          if (min == RDB[wwd + WWD_LIM_MAX])
            fprintf(outp, " - Maximum source importance: %1.5E (lim)\n", max);
          else
            fprintf(outp, " - Maximum source importance: %1.5E\n", max);
        }

      /* Newline */

      fprintf(outp, "\n");

      /* Allocate memory for statistics in source importance calculation */

      if ((long)RDB[DATA_SRC_IMP_CALC] == YES)
        {
          /* Get number of energy groups */

          if ((long)RDB[wwd + WWD_PTR_ERG] > VALID_PTR)
            ne = (long)RDB[wwd + WWD_NE];
          else
            ne = 1;

          /* Get number of spatial bins */

          n = (long)(RDB[msh + MESH_N0]*RDB[msh + MESH_N1]*RDB[msh + MESH_N2]);

          /* Allocate memory */

          ptr = NewStat("WWD_SRC_IMP", 2, ne, n);
          WDB[wwd + WWD_PTR_SRC_IMP_DIS] = (double)ptr;
        }

      /* Next */

      wwd = NextItem(wwd);
    }

  /***************************************************************************/
}

/*****************************************************************************/
