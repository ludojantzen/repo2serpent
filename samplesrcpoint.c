/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : samplesrcpoint.c                               */
/*                                                                           */
/* Created:       2011/03/02 (JLe)                                           */
/* Last modified: 2020/06/29 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Samples neutron source point                                 */
/*                                                                           */
/* Comments: - Toi cell-haku ei toimi ihan oikein sillä whereami palauttaa   */
/*             alimman cellin. Universumihakua ei oo edes tehty vielä.       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleSrcPoint:"
#define MAX_SRC_RESAMPLE 1000000

/*****************************************************************************/

long SampleSrcPoint(long id, long np, long idx)
{
  long n, cell, mat, src, surf, lst, loc0, loc1, rea, ptr, cell0, mat0, ncol;
  long loc2, type, stp, part, npmin, npmax, idx0, i, Ecutflag, m, ns, in;
  unsigned long seed;
  double rnd, mu, x, y, z, u, v, w, E, wgt, xmin, xmax, ymin, ymax, zmin, zmax;
  double rmin, rmax, xt, yt, zt, t, dummy, wgt0, wgt1, meanw, meanb, val, r;
  double theta;

  /***************************************************************************/

  /***** Divide source to MPI tasks ******************************************/

#ifdef MPI_MODE1

  /* Check number of tasks */
  /* Division isn't done for transient modes (indicated with np < 0) */

  if ((mpitasks > 1) && (np >= 0))
    {
      /* Calculate number of particles per task */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        npmax = (long)(RDB[DATA_CRIT_POP]/((double)mpitasks));
      else if ((ptr = (long)RDB[DATA_PTR_PREC_DET]) > VALID_PTR)
        npmax = (long)RDB[ptr + PRECDET_N_LIVE];
      else
        npmax = (long)(RDB[DATA_SRC_POP]/((double)mpitasks));

      /* Calculate minimum and maximum particle index for this task */

      npmin = mpiid*npmax;

      if (mpiid < mpitasks - 1)
        npmax = (mpiid + 1)*npmax;
      else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        npmax = (long)RDB[DATA_CRIT_POP];
      else if ((ptr = (long)RDB[DATA_PTR_PREC_DET]) > VALID_PTR)
        npmax = (long)RDB[ptr + PRECDET_N_LIVE];
      else
        npmax = (long)RDB[DATA_SRC_POP];

      /* Check with task number */

      if ((np < npmin) || (np > npmax - 1))
        return -1;
    }

#endif

  /***************************************************************************/

  /***** Source sampling starts here *****************************************/

  /* Reset material pointer and source point coordinates passed to */
  /* BoundaryConditions() */

  mat = -1;
  xt = 0.0;
  yt = 0.0;
  zt = 0.0;

  /* Init random number sequence */

  seed = ReInitRNG(idx);
  SEED[id*RNG_SZ] = seed;

  /* Restore initial seed for debugging */

  SEED0[id*RNG_SZ] = seed;

  /* Store history index for debugging */

  ptr = (long)RDB[DATA_PTR_PRIVA_HIS_IDX];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  PutPrivateData(ptr, idx, id);

  /* Set cell search list option */

  ptr = (long)RDB[DATA_CELL_SEARCH_LIST];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  PutPrivateData(ptr, (double)CELL_SEARCH_LIST_SRC, id);

  /* Init direction cosines for call to WhereAmI() */

  u = 0.0;
  v = 0.0;
  w = 1.0;

  /* Reset time */

  t = 0.0;

  /* Check type */

  if ((long)RDB[DATA_RMX_CONVG_ITER_IDX] > 0)
    {
      /***********************************************************************/

      /***** Distribution from response matrix solver ************************/

      /* NOTE: DD-moodissa oikea domaini varmistetaan nyt sillä, että    */
      /*       MatPtr():stä antaman materiaalin fissile flag on päällä.  */
      /*       Tuossa jälkimmäisessä silmukassa tarkistetaan myös parent */
      /*       materiaalin flagi jotta ei arvottaisi turhaan pisteitä    */
      /*       muualle kuin polttoaineeseen. Lopussa tehdään vielä yksi  */
      /*       tarkastus CheckDDDomai():lla. */

      /* Check iteration index */

      if (((long)RDB[DATA_RMX_CONVG_ITER_IDX] < (long)RDB[DATA_RMX_CONVG_OUT_N])
          || ((long)RDB[DATA_UFS_MODE] == UFS_MODE_RMX))
        {
          /* Pointer to mesh */

          loc0 = (long)RDB[DATA_PTR_RMX0];
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          loc1 = (long)RDB[loc0 + RMX_PTR_MESH];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Sampling loop */

          for (n = 0; n < MAX_SRC_RESAMPLE; n++)
            {
              /* Sample coordinates */

              x = RandF(id)*(RDB[DATA_GEOM_MAXX] - RDB[DATA_GEOM_MINX])
                + RDB[DATA_GEOM_MINX];
              y = RandF(id)*(RDB[DATA_GEOM_MAXY] - RDB[DATA_GEOM_MINY])
                + RDB[DATA_GEOM_MINY];

              if ((long)RDB[DATA_GEOM_DIM] == 3)
                z = RandF(id)*(RDB[DATA_GEOM_MAXZ] - RDB[DATA_GEOM_MINZ])
                  + RDB[DATA_GEOM_MINZ];
              else
                z = 0.0;

              /* Check true symmetry */

              if (TestUniSym(x, y, z) == YES)
                continue;

              /* Find location */

              if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
                {
                  /* Apply boundary conditions */

                  BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w,
                                     &xt, &yt, &zt, &dummy, id);

                  /* Check cell pointer */

                  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

                  /* Get material pointer */

                  if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                    mat = MatPtr(mat, id);

                  /* Check if cell contains fissile material */

                  if (mat > VALID_PTR)
                    if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT))
                      {
                        /* Pointer to global mesh */

                        ptr = MeshPtr(loc1, x, y, z);
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                        ptr = (long)RDB[ptr];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                        /* Pointer to local mesh */

                        ptr = (long)RDB[ptr + RMX_CELL_PTR_SUBMESH];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                        /* Set weight */

                        if ((wgt = MeshVal(ptr, x, y, z)) > 0.0)
                          break;
                      }
                }
            }
        }
      else
        {
          /* Loop over resamples */

          for (n = 0; n < MAX_SRC_RESAMPLE; n++)
            {
              /* Sample random number */

              rnd = RandF(id);

              /* Pointer to response matrix */

              loc0 = (long)RDB[DATA_PTR_RMX0];
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

              /* Loop over data */

              loc1 = (long)RDB[loc0 + RMX_PTR_MESH_DATA];
              while (loc1 > VALID_PTR)
                {
                  /* Compare probability */

                  if ((rnd = rnd - RDB[loc1 + RMX_CELL_SRC_PROB]) < 0.0)
                    break;

                  /* Next */

                  loc1 = NextItem(loc1);
                }

              /* Check */

              if (loc1 < VALID_PTR)
                Die(FUNCTION_NAME, "Source sampling failed");

              /* Sample random number */

              rnd = RandF(id);

              /* Loop over sub-mesh data */

              loc2 = (long)RDB[loc1 + RMX_CELL_PTR_SUBMESH_SAMPLE];
              while (loc2 > VALID_PTR)
                {
                  /* Compare probability */

                  if ((rnd = rnd - RDB[loc2 + RMX_SUBMESH_SAMPLE_P]) < 0.0)
                    break;

                  /* Next */

                  loc2 = NextItem(loc2);
                }

              /* Check */

              if (loc2 < VALID_PTR)
                Die(FUNCTION_NAME, "Source sampling failed");

              /* Get index */

              i = (long)RDB[loc2 + RMX_SUBMESH_SAMPLE_IDX];
              CheckValue(FUNCTION_NAME, "i", "", i, 0, 1E+9);

              /* Pointer to mesh */

              loc2 = (long)RDB[loc1 + RMX_CELL_PTR_SUBMESH];
              CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

              /* Get boundaries */

              MeshCellBounds(loc2, i, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

              /* Sample point in fissile material */

              for (m = 0; m < MAX_SRC_RESAMPLE; m++)
                {
                  /* Sample coordinates */

                  x = RandF(id)*(xmax - xmin) + xmin;
                  y = RandF(id)*(ymax - ymin) + ymin;

                  if ((long)RDB[DATA_GEOM_DIM] == 3)
                    z = RandF(id)*(zmax - zmin) + zmin;
                  else
                    z = 0.0;

                  /* Check true symmetry */

                  if (TestUniSym(x, y, z) == YES)
                    continue;

                  /* Check index and reject */

                  if (i != MeshIndex(loc2, x, y, z, 0.0))
                    continue;

                  /* Find location */

                  if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
                    {
                      /* Apply boundary conditions */

                      BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w,
                                     &xt, &yt, &zt, &dummy, id);

                      /* Check cell pointer */

                      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

                      /* Check fissile */

                      if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                        if (((long)RDB[mat + MATERIAL_OPTIONS] &
                             OPT_FISSILE_MAT))
                          break;
                    }
                }

              /* Check error */

              if (m == MAX_SRC_RESAMPLE)
                Error(0, "Unable to sample fission source");

              /* Get pointer */

              mat = MatPtr(mat, id);
              CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

              /* Check if fissile */

              if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT))
                break;
            }

          /* Set weight to unity */

          wgt = 1.0;
        }

      /* Check error */

      if (n == MAX_SRC_RESAMPLE)
        Error(0, "Unable to sample fission source");

      /* Sample isotropic direction */

      IsotropicDirection(&u, &v, &w, id);

      /* Energy from a maxwellian distribution */

      E = MaxwellEnergy(1.2895, id);

      /* Set weight */

      wgt1 = wgt;

      /* Set type to neutron */

      type = PARTICLE_TYPE_NEUTRON;

      /* Reset source pointer */

      src = -1;

      /***********************************************************************/
    }
  else if ((src = (long)RDB[DATA_PTR_SRC0]) < VALID_PTR)
    {
      /***********************************************************************/

      /***** Default uniform fission source **********************************/

      /* Check mode */

      if (((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT) &&
          ((long)RDB[DATA_RMX_CONVG_ACC] == NO))
        Die(FUNCTION_NAME, "No source definition");

      /* Loop until neutron is in fissile material */

      for (n = 0; n < MAX_SRC_RESAMPLE; n++)
        {
          /* Sample coordinates */

          x = RandF(id)*(RDB[DATA_GEOM_MAXX] - RDB[DATA_GEOM_MINX])
            + RDB[DATA_GEOM_MINX];
          y = RandF(id)*(RDB[DATA_GEOM_MAXY] - RDB[DATA_GEOM_MINY])
            + RDB[DATA_GEOM_MINY];

          if ((long)RDB[DATA_GEOM_DIM] == 3)
            z = RandF(id)*(RDB[DATA_GEOM_MAXZ] - RDB[DATA_GEOM_MINZ])
              + RDB[DATA_GEOM_MINZ];
          else
            z = 0.0;

          /* Check true symmetry */

          if (TestUniSym(x, y, z) == YES)
            continue;

          /* Find location */

          if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
            {
              /* Apply boundary conditions */

              BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w, &xt, &yt, &zt,
                                 &dummy, id);

              /* Check cell pointer */

              CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

              /* Get material pointer */

              if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                mat = MatPtr(mat, id);

              /* Check if cell contains fissile material */

              if (mat > VALID_PTR)
                if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
                    ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) ||
                    ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN) ||
                    ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN))
                  if (((long)RDB[DATA_DD_DECOMPOSE] == NO) ||
                      ((long)RDB[mat + MATERIAL_MPI_ID] < 0) ||
                      ((long)RDB[mat + MATERIAL_MPI_ID] == mpiid))
                    {
                      /* Sample isotropic direction */

                      IsotropicDirection(&u, &v, &w, id);

                      /* Energy from a maxwellian distribution */

                      E = MaxwellEnergy(1.2895, id);

                      /* Set weight to unity */

                      wgt = 1.0;
                      wgt1 = 1.0;

                      /* Break loop */

                      break;
                    }
            }
        }

      /* Set type to neutron */

      type = PARTICLE_TYPE_NEUTRON;

      /* Check error */

      if (n == MAX_SRC_RESAMPLE)
        Error(0, "Unable to sample fission source - try explicit definition");

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Source from user-defined distribution ***************************/

      /* Sample source definition */

      rnd = RandF(id);

      /* loop over sources */

      src = (long)RDB[DATA_PTR_SRC0];
      while (src > VALID_PTR)
        {
          /* Compare to weight */

          if ((rnd = rnd - RDB[src + SRC_WGT]) < 0.0)
            break;

          /* Next source definition */

          src = NextItem(src);
        }

      /* Check pointer */

      if (src < VALID_PTR)
        Die(FUNCTION_NAME, "Unable to sample source");

      /* Set type */

      type = (long)RDB[src + SRC_TYPE];

      /* Check */

      if ((type != PARTICLE_TYPE_NEUTRON) && (type != PARTICLE_TYPE_GAMMA))
        Die(FUNCTION_NAME, "Invalid particle type");

      /* Reset mean weight and boundary */

      meanw = 0.0;
      meanb = 0.0;

      /* Outer rejection loop */

      for (m = 0; m < MAX_SRC_RESAMPLE; m++)
        {
          /*******************************************************************/

          /***** Time ********************************************************/

          /* Check if interval is given */

          if (RDB[src + SRC_TMAX] > RDB[src + SRC_TMIN])
            {
              /* Sample */

              t = RandF(id)*(RDB[src + SRC_TMAX] - RDB[src + SRC_TMIN]) +
                RDB[src + SRC_TMIN];

              /* Check */

              if ((t < RDB[DATA_TIME_CUT_TMIN]) ||
                  (t > RDB[DATA_TIME_CUT_TMAX]))
                Error(src, "Emission time is outside cut-off boundaries");
            }

          /*******************************************************************/

          /***** Direction ***************************************************/

          /* Check if isotropic (this is set to INFTY in processsources.c) */

          if (RDB[src + SRC_U0] > 1.0)
            IsotropicDirection(&u, &v, &w, id);
          else
            {
              /* Set values */

              u = RDB[src + SRC_U0];
              v = RDB[src + SRC_V0];
              w = RDB[src + SRC_W0];
            }

          /*******************************************************************/

          /***** Energy ******************************************************/

          /* Check source energy */

          if ((E = RDB[src + SRC_E]) > ZERO)
            {
              /* Check particle type */

              if (type == PARTICLE_TYPE_NEUTRON)
                {
                  /* Compare to upper boundary */

                  if (E > RDB[DATA_NEUTRON_EMAX])
                    Error(src, "Source energy above maximum (%1.2f MeV)",
                          RDB[DATA_NEUTRON_EMAX]);

                  /* Compare to lower boundary */

                  if (E < RDB[DATA_NEUTRON_EMIN])
                    Error(src, "Source energy below minimum (%1.2E MeV)",
                          RDB[DATA_NEUTRON_EMIN]);
                }
              else
                {
                  /* Compare to upper boundary */

                  if (E > RDB[DATA_PHOTON_EMAX])
                    Error(src, "Source energy above maximum (%1.2f MeV)",
                          RDB[DATA_PHOTON_EMAX]);

                  /* Compare to lower boundary */

                  if (E < RDB[DATA_PHOTON_EMIN])
                    Error(src, "Source energy below minimum (%1.2E MeV)",
                          RDB[DATA_PHOTON_EMIN]);
                }
            }

          /* Initialize weight */

          wgt0 = 1.0;

          /* Check type */

          if ((lst = (long)RDB[src + SRC_PTR_EBINS]) > VALID_PTR)
            {
              /* Get pointers to vectors */

              ptr = (long)RDB[lst + SRC_EBIN_PTR_E];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              loc0 = (long)RDB[lst + SRC_EBIN_PTR_PDF];
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

              loc1 = (long)RDB[lst + SRC_EBIN_PTR_CDF];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

              /* Get number of points */

              n = (long)RDB[lst + SRC_EBIN_NE];
              CheckValue(FUNCTION_NAME, "n", "", n, 2, 100000000);

              /* Get interpolation */

              i = (long)RDB[lst + SRC_EBIN_INTT];
              CheckValue(FUNCTION_NAME, "i", "", i, 0, 4);

              /* Sample from table */

              E = SampleTabular(&RDB[ptr], &RDB[loc0], &RDB[loc1], n, i, id);

              /* Check */

              if (E < 0.0)
                continue;
            }
          else if ((rea = (long)RDB[src + SRC_PTR_REA]) > VALID_PTR)
            {
              /* Put incident energy */

              if ((E = RDB[src + SRC_E]) < ZERO)
                E = 1.000001*RDB[rea + REACTION_EMIN];

              /* Init cosine (-1.0, 0.0 and 1.0 case dubious error) */

              mu = 0.5;

              /* Sample energy and direction */

              SampleENDFLaw(rea, -1, E, &E, &mu, id);

              /* Check if cosine is changed */

              if (mu != 0.5)
                {
                  /* Sanity check for mu and direction vectors */
                  /* (for NAN's etc.) */

                  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.01, 1.01);
                  CheckValue(FUNCTION_NAME, "u", "", u, -1.01, 1.01);
                  CheckValue(FUNCTION_NAME, "v", "", v, -1.01, 1.01);
                  CheckValue(FUNCTION_NAME, "w", "", w, -1.01, 1.01);

                  /* Rotate */

                  AziRot(mu, &u, &v, &w, id);
                }
            }
          else if ((E = RDB[src + SRC_E]) < ZERO)
            E = 1.0;

          /*******************************************************************/

          /***** Position ****************************************************/

          /* Get boundaries */

          if ((xmin = RDB[src + SRC_XMIN]) == -INFTY)
            xmin = RDB[DATA_GEOM_MINX];

          if ((xmax = RDB[src + SRC_XMAX]) == INFTY)
            xmax = RDB[DATA_GEOM_MAXX];

          if ((ymin = RDB[src + SRC_YMIN]) == -INFTY)
            ymin = RDB[DATA_GEOM_MINY];

          if ((ymax = RDB[src + SRC_YMAX]) == INFTY)
            ymax = RDB[DATA_GEOM_MAXY];

          if ((zmin = RDB[src + SRC_ZMIN]) == -INFTY)
            zmin = RDB[DATA_GEOM_MINZ];

          if ((zmax = RDB[src + SRC_ZMAX]) == INFTY)
            zmax = RDB[DATA_GEOM_MAXZ];

          rmin = RDB[src + SRC_RMIN];
          rmax = RDB[src + SRC_RMAX];

          /* Pointer to cell and material */

          cell0 = (long)RDB[src + SRC_PTR_CELL];
          mat0 = (long)RDB[src + SRC_PTR_MAT];

          /* Position re-sampling loop */

          for (n = 0; n < MAX_SRC_RESAMPLE; n++)
            {
              /* Reset weight */

              wgt = wgt0;

              /* Add to track counter */

              ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
              AddPrivateData(ptr, 1.0, id);

              /* Check type */

              if ((long)RDB[src + SRC_READ_FILE_TYPE] ==
                  SRC_FILE_TYPE_FUSION_PLASMA)
                {
                  /* Fusion plasma source */

                  SamplePlasmaSrc(src, &x, &y, &z, &u, &v, &w, &E, &wgt,
                                  &t, id);
                }
              else if ((long)RDB[src + SRC_READ_PTR_FILE] > VALID_PTR)
                {
                  /* Read parameters from file (all variables are overriden) */

                  ReadSourceFile(src, &x, &y, &z, &u, &v, &w, &E, &wgt, &t);

                  /* Transfer point (10.3.2016 / 2.1.26) */

                  if ((RDB[src + SRC_X0] > -INFTY) &&
                      (RDB[src + SRC_Y0] > -INFTY) &&
                      (RDB[src + SRC_Z0] > -INFTY))
                    {
                      x = x + RDB[src + SRC_X0];
                      y = y + RDB[src + SRC_Y0];
                      z = z + RDB[src + SRC_Z0];
                    }

                  /* Check true symmetry */

                  if (TestUniSym(x, y, z) == YES)
                    continue;

                  /* Check that point is inside the geometry (x and y)*/

                  if (((x <= RDB[DATA_GEOM_MINX]) ||
                       (x >= RDB[DATA_GEOM_MAXX])) ||
                      ((y <= RDB[DATA_GEOM_MINY]) ||
                       (y >= RDB[DATA_GEOM_MAXY])))
                    continue;

                  /* Check z-coordinate only in 3D geometries */

                  if ((RDB[DATA_GEOM_DIM] == 3) &&
                      ((z <= RDB[DATA_GEOM_MINZ]) ||
                       (z >= RDB[DATA_GEOM_MAXZ])))
                    continue;
                }
              else if ((RDB[src + SRC_X0] > -INFTY) &&
                       (RDB[src + SRC_Y0] > -INFTY) &&
                       (RDB[src + SRC_Z0] > -INFTY) &&
                       (rmax < rmin))
                {
                  /* Check material and cell pointers */

                  if (mat0 > VALID_PTR)
                    Error(src, "Point source not allowed with material source");
                  else if (cell0 > VALID_PTR)
                    Error(src, "Point source not allowed with cell source");

                  /* Point source, set coordinates */

                  x = RDB[src + SRC_X0];
                  y = RDB[src + SRC_Y0];
                  z = RDB[src + SRC_Z0];
                }
              else if ((surf = (long)RDB[src + SRC_PTR_SURF]) > VALID_PTR)
                {
                  /* Sample point on surface */

                  SurfaceSrc(src, surf, &x, &y, &z, &u, &v, &w, id);
                }
              else
                {
                  /* Sample coordinates */

                  if (rmax < rmin)
                    {
                      /* Cartesian */

                      x = RandF(id)*(xmax - xmin) + xmin;
                      y = RandF(id)*(ymax - ymin) + ymin;
                    }
                  else
                    {
                      /* Cylindrical */

                      r = sqrt(RandF(id)*(rmax*rmax - rmin*rmin) + rmin*rmin);
                      theta = RandF(id)*2.0*PI;

                      x = r*cos(theta);
                      y = r*sin(theta);

                      /* Check if origin is shifted */

                      if (RDB[src + SRC_X0] > -INFTY)
                        x = x + RDB[src + SRC_X0];

                      if (RDB[src + SRC_Y0] > -INFTY)
                        y = y + RDB[src + SRC_Y0];
                    }

                  /* Sample axial coordinate */

                  if ((long)RDB[DATA_GEOM_DIM] == 3)
                    z = RandF(id)*(zmax - zmin) + zmin;
                  else
                    z = 0.0;
                }

              /* Sample using source biasing from weight window mesh */

              if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == YES)
                {
                  /* Sample coordinates */

                  WWinSrc(src, type, &x, &y, &z, E, &wgt, id);

                  /* Check dimension and override axial */

                  if ((long)RDB[DATA_GEOM_DIM] == 2)
                    z = 0.0;
                }

              /* Override everything with user-defined subroutine */

              if ((long)RDB[src + SRC_PTR_USR] > VALID_PTR)
                UserSrc(src, &x, &y, &z, &u, &v, &w, &E, &wgt, &t, id);

              /* Check boundaries */

              if ((x > RDB[src + SRC_XMAX]) || (x < RDB[src + SRC_XMIN]) ||
                  (y > RDB[src + SRC_YMAX]) || (y < RDB[src + SRC_YMIN]) ||
                  (z > RDB[src + SRC_ZMAX]) || (z < RDB[src + SRC_ZMIN]))
                {
                  /* Not within boundaries */

                  Die(FUNCTION_NAME, "Source point outside boundaries");
                }

              /* Check true symmetry */

              if (TestUniSym(x, y, z) == YES)
                continue;

              /* Check time cut-off (tää ei ehkä säilytä */
              /* kokonaislähdenopeutta?) */

              if ((t >= RDB[DATA_TIME_CUT_TMIN]) &&
                  (t < RDB[DATA_TIME_CUT_TMAX]))
                if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
                  {
                    /* Check if only surface or point source is used */

                    if ((cell0 < VALID_PTR) && (mat0 < VALID_PTR) &&
                        ((long)RDB[DATA_USE_DECAY_SRC] == NO) &&
                        ((long)RDB[cell + CELL_TYPE] != CELL_TYPE_OUTSIDE))
                      break;

                    /* Check cell */

                    if (cell0 > VALID_PTR)
                      {
                        /* Get collision number */

                        ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
                        ncol = (long)GetPrivateData(ptr, id);

                        /* Check collision */

                        if (TestValuePair(cell0 + CELL_COL_COUNT, (double)ncol,
                                          id) > 0.0)
                          break;
                      }

                    /* Get material pointer */

                    if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                      {
                        mat = MatPtr(mat, id);
                        CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
                      }

                    /* Radioactive decay source */

                    if ((mat > VALID_PTR) &&
                        ((long)RDB[DATA_USE_DECAY_SRC] == YES))
                      {
                        /* Check material and cell pointers */

                        if (mat0 > VALID_PTR)
                          Error(src,
                                "Decay source not allowed with material source");
                        else if (cell0 > VALID_PTR)
                          Error(src, "Decay source not allowed with cell source");

                        /* Call source routine */

                        if (RadGammaSrc(src, mat, &E, &wgt, &in, id) >
                            VALID_PTR)
                          {
                            /* Score source weight */

                            ptr = (long)RDB[mat + MATERIAL_SAMPLED_DECAY_SRC];
                            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,
                                         ptr);
                            AddBuf1D(wgt, 1.0, ptr, id, 2 - type);

                            /* Break loop */

                            break;
                          }
                      }

                    /* Check material */

                    if ((mat0 > VALID_PTR) && (mat > VALID_PTR))
                      {
                        /* Check match */

                        if (mat == mat0)
                          break;

                        /* Check match with parent */

                        if ((mat = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) >
                            VALID_PTR)
                          if (mat == mat0)
                            break;
                      }
                  }
            }

          /* Check failure */

          if (n == MAX_SRC_RESAMPLE)
            Error(0, "Source sampling failed because of low efficiency");

          /* Check cell pointer */

          if (cell < VALID_PTR)
            Die(FUNCTION_NAME, "Cell pointer is NULL");

          /* Check if weight windows are used */

          if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO)
            {
              /* Set weight */

              wgt1 = wgt;

              /* Score source splitting */

              ptr = (long)RDB[RES_SRC_WW_SPLIT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(1.0, 1.0, ptr, id, 2 - type);

              /* Break loop */

              break;
            }
          else
            {
              /* Score source efficiency */

              ptr = (long)RDB[RES_SRC_WW_EFF];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(1.0, 1.0, ptr, id, -1, 2 - type, 1);
            }

          /* Reject by importance */

          if ((wgt1 = SrcRoulette(type, x, y, z, E, wgt, &val, id)) < 0.0)
            {
              /* Rejected, check particle type */

              if (type == PARTICLE_TYPE_GAMMA)
                {
                  /* Score source rate */

                  stp = (long)RDB[RES_TOT_PHOTON_SRCRATE];
                  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
                  AddBuf1D(1.0, wgt, stp, id, 0);
                }
              else
                {
                  /* Score source rate */

                  stp = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
                  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
                  AddBuf1D(1.0, wgt, stp, id, 0);

                  /* Score source rate in fissile and non-fissile materials */

                  if (mat > VALID_PTR)
                    {
                      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
                        AddBuf1D(1.0, wgt, stp, id, 1);
                      else
                        AddBuf1D(1.0, wgt, stp, id, 2);
                    }
                }

              /* Score source efficiency */

              ptr = (long)RDB[RES_SRC_WW_EFF];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(1.0, 1.0, ptr, id, -1, 2 - type, 0);

              /* Add to mean weight and boundary */

              meanw = meanw + wgt;
              meanb = meanb + RDB[DATA_WWD_LOWER_BOUND]/val;
            }
          else
            {
              /* Calculate source splitting (this is applied later) */

              ns = (long)(wgt*val/RDB[DATA_WWD_UPPER_BOUND]) + 1;

              /* Score source splitting */

              ptr = (long)RDB[RES_SRC_WW_SPLIT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D((double)ns, 1.0, ptr, id, 2 - type);

              /* Check excessive splitting */

              if (ns > 1000)
                Note(0, "Excessive source splitting (%E)",
                     wgt*val/RDB[DATA_WWD_UPPER_BOUND]);

              /* Break loop */

              break;
            }
        }

      /* Check failure */

      if (m == MAX_SRC_RESAMPLE)
        Error(0, "Source sampling in ww mesh failed (%E %E)\n",
              meanw/((double)m), meanb/((double)m));

      /* Get material pointer */

      if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
        {
          mat = MatPtr(mat, id);
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
        }

      /* Score sampling efficiency */

      ptr = (long)RDB[RES_SRC_SAMPLING_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, 1.0, ptr, id, 2 - type);
      AddBuf1D((double)((n + 1)*(m + 1)), 1.0, ptr, id, 4 - type);

      /* Score mean weight */

      ptr = (long)RDB[RES_SRC_MEAN_WGT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(wgt1, 1.0, ptr, id, 2 - type);

      /***********************************************************************/
    }

  /* Check that particle is in correct domain */

  if (CheckDDDomain(mat) == NO)
    Die(FUNCTION_NAME, "Mismatch in domain");

  /* Reset energy cutoff flag */

  Ecutflag = NO;

  /* Adjust minimum and maximum energy */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      if (E < 1.0000001*RDB[DATA_NEUTRON_EMIN])
        E = 1.000001*RDB[DATA_NEUTRON_EMIN];
      else if (E > 0.999999*RDB[DATA_NEUTRON_EMAX])
        E = 0.999999*RDB[DATA_NEUTRON_EMAX];
    }
  else
    {
      /* Apply energy cutoff */

      if ((E < 1.0000001*RDB[DATA_PHOTON_EMIN]) ||
          ((mat > VALID_PTR) && (E < RDB[mat + MATERIAL_PHOTON_ECUT])))
        {
          /* Set cutoff flag */

          Ecutflag = YES;
        }
      else if ((E > 0.999999*RDB[DATA_PHOTON_EMAX]) &&
               (E <= RDB[DATA_PHOTON_EMAX]))
        {
          /* Adjust the energy */

          E = 0.999999*RDB[DATA_PHOTON_EMAX];
        }
      else if (E > RDB[DATA_PHOTON_EMAX])
        Die(FUNCTION_NAME, "Photon energy %E MeV above maximum %E MeV",
            E, RDB[DATA_PHOTON_EMAX]);
    }

  /* Check time (noi karsitaan jo tuolla ylempänä) */

  if (src > VALID_PTR)
    {
      /* Lower limit */

      if (t < RDB[DATA_TIME_CUT_TMIN])
        Error(src, "Source time %1.2E is below time cutoff %1.2E", t,
              RDB[DATA_TIME_CUT_TMIN]);

      /* Upper limit */

      if (t >= RDB[DATA_TIME_CUT_TMAX])
        Error(src, "Source time %1.2E is above time cutoff %1.2E", t,
              RDB[DATA_TIME_CUT_TMAX]);
    }

  /* Get particle from stack */

  part = FromStack(type, id);

  /* Put values */

  WDB[part + PARTICLE_X] = x;
  WDB[part + PARTICLE_Y] = y;
  WDB[part + PARTICLE_Z] = z;

  WDB[part + PARTICLE_U] = u;
  WDB[part + PARTICLE_V] = v;
  WDB[part + PARTICLE_W] = w;

  WDB[part + PARTICLE_E] = E;
  WDB[part + PARTICLE_WGT] = wgt1;
  WDB[part + PARTICLE_T0] = t;
  WDB[part + PARTICLE_T] = t;
  WDB[part + PARTICLE_TD] = 0.0;
  WDB[part + PARTICLE_TT] = 0.0;
  WDB[part + PARTICLE_COL_IDX] = 0.0;

  WDB[part + PARTICLE_PTR_MAT] = (double)mat;
  WDB[part + PARTICLE_RNG_IDX] = (double)idx;
  WDB[part + PARTICLE_HISTORY_IDX] = (double)idx;
  WDB[part + PARTICLE_PREV_I] = RDB[cell + CELL_IMP];

  /* Set photon emission type */

  if (type == PARTICLE_TYPE_GAMMA)
    WDB[part + PARTICLE_PHOTON_TYPE] = PHOTON_TYPE_SRC;

  /* Reset detector flags */

  WDB[part + PARTICLE_DET_FLAGS] = 0.0;

  /* Get fission matrix index */

  idx0 = FissMtxIndex(mat, id);

  /* Put fission matrix index */

  WDB[part + PARTICLE_FMTX_IDX] = (double)idx0;

  /* Score mesh plotter */

  ScoreMesh(part, mat, 0.0, -1.0, x, y, z, E, t, wgt, -1.0, id);

  /* Score source detector */

  SrcDet(part, mat, x, y, z, u, v, w, E, t, wgt, id);

  /* Kill photon if energy cutoff flag was set (11.3.2017 / 2.1.29 / TKa) */
  /* TODO: Onko tämä oikeassa kohdassa? */

  if ((Ecutflag == YES) && (type == PARTICLE_TYPE_GAMMA))
    {
      /* Score source rate */

      stp = (long)RDB[RES_TOT_PHOTON_SRCRATE];
      CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
      AddBuf1D(1.0, wgt, stp, id, 0);

      /* Score photon heat detector */

      ScorePhotonHeat(part, mat, E, x, y, z, E, t, wgt, id);

      /* Score pulse-height detector */

      PulseDet(part, mat, E, x, y, z, wgt, id);

      /* Put particle back in stack */

      ToStack(part, id);

      /* Exit subroutine */

      return -1;
    }

  /* Score response matrix source term */

  ScoreRMXSrc(part, x, y, z, E, wgt, id);

  /* Get maximum source point importance */

  MaxSrcImp(type, x, y, z, E, id);

  /* Check simulation mode */

  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) ||
      ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN) ||
      ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN))
    {
      /* Score source rate for fission matrix */

      if (idx0 > -1)
        {
          /* Pointer to matrix */

          ptr = (long)RDB[DATA_PTR_FMTX];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          ptr = (long)RDB[ptr + FMTX_PTR_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Score total */

          AddBuf(wgt, 1.0, ptr, id, -1, 0, idx0);

          /* Score prompt (no delayed in source) */

          AddBuf(wgt, 1.0, ptr, id, -1, 1, idx0);
        }

      /* Check particle type */

      if (type == PARTICLE_TYPE_GAMMA)
        {
          /* Score source rate */

          stp = (long)RDB[RES_TOT_PHOTON_SRCRATE];
          CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
          AddBuf1D(1.0, wgt, stp, id, 0);

          /* Score particle balance */

          ptr = (long)RDB[RES_G_BALA_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_SRC_EXT, 0);
          AddBuf(wgt1, 1.0, ptr, id, -1, BALA_G_SRC_EXT, 1);
          AddBuf(wgt1*E, 1.0, ptr, id, -1, BALA_G_SRC_EXT, 2);
        }
      else
        {
          /* Score source rate */

          stp = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
          CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
          AddBuf1D(1.0, wgt, stp, id, 0);

          /* Score initial source weight */

          ptr = (long)RDB[RES_INI_SRC_WGT];
          CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
          AddBuf1D(1.0, wgt, ptr, id, 0);

          /* Score source rate in fissile and non-fissile materials */

          if (mat > VALID_PTR)
            {
              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
                AddBuf1D(1.0, wgt, stp, id, 1);
              else
                AddBuf1D(1.0, wgt, stp, id, 2);
            }

          /* Score particle balance */

          ptr = (long)RDB[RES_N_BALA_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(1.0, 1.0, ptr, id, -1, BALA_N_SRC_EXT, 0);
          AddBuf(wgt1, 1.0, ptr, id, -1, BALA_N_SRC_EXT, 1);
        }

      /* Score time-dependent source rate detectors */

      ScoreTimeSource(part, mat, MT_PRIMARY_LIVE_SOURCE,
                      x, y, z, E, t+ZERO, 1, wgt, id);

      /* Put MPI id */

      WDB[part + PARTICLE_MPI_ID] = (double)mpiid;

      /* Reset generation index */

      WDB[part + PARTICLE_GEN_IDX] = 0.0;

      /* Apply weight window */

      if (WeightWindow(-1, part, type, x, y, z, u, v, w, E, &wgt1, t,
                       WWMESH_SRC, id) == TRACK_END_WCUT)
        {
          /* This should not happen */

          Die(FUNCTION_NAME, "Weight cut-off in source routine");
        }
      else
        {
          /* Store new weight */

          WDB[part + PARTICLE_WGT] = wgt1;

          /* Put particle in que */

          ToQue(part, id);
        }
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Add to simulated batch size (this is different from DATA_NBATCH in */
      /* MPI mode, and may be different for each parallel task) NOTE: tän   */
      /* voisi tehdä jotenkin fiksumminkin jos jossain vaiheessa katsotaan  */
      /* että miten noi historiat jaetaan eri taskeille. */

#ifdef OPEN_MP
#pragma omp atomic
#endif
      WDB[DATA_SIMUL_BATCH_SIZE]++;

      /* Number of independent populations */

      n = (long)RDB[DATA_N_POP_EIG];
      CheckValue(FUNCTION_NAME, "n", "", n, 1, 99);

      /* Eigenvalue idex */

      WDB[part + PARTICLE_EIG_IDX] = (double)(idx%n);

      /* Put particle to bank */

      ToBank(part, id);
    }
  else
    Die(FUNCTION_NAME, "Invalid simulation mode");

  /* Return pointer */

  return part;
}

/*****************************************************************************/
