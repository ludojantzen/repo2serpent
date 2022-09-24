/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : volumesmc.c                                    */
/*                                                                           */
/* Created:       2010/11/10 (JLe)                                           */
/* Last modified: 2020/06/23 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates volumes by Monte Carlo simulation                 */
/*                                                                           */
/* Comments: - Tän rutiinin ajaminen vaikuttaa cell search listien           */
/*             statistiikkaan.                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "VolumesMC"

/*****************************************************************************/

void VolumesMC()
{
  long cell, mat, mat0, ptr, loc0, nt, idx, m, n, nmax, nb, id, dim, mode;
  long algo, ptl, umsh, ifc;
  unsigned long seed;
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, u, v, w, vol, f, T, d;
  double max, err, tmax, emax, t, val, est, diff, df, lmax, div;
  char tmpstr[MAX_STR], fname[MAX_STR], fmt[MAX_STR];
  FILE *fp;

  /* Check mode */

  if ((long)RDB[DATA_VOLUME_MC] == NO)
    return;

  /* Mode (1 = Material volumes / 2 = cell volumes) */

  mode = 1;

  /* Cut-offs */

  nmax = (long)RDB[DATA_VOLUME_MC_NMAX];
  tmax = RDB[DATA_VOLUME_MC_TMAX];
  emax = RDB[DATA_VOLUME_MC_EMAX];
  algo = -1;

  /* Algorithm (1 = points / 2 = tracks) */

  if (nmax > 0)
    algo = 1;
  else if (nmax < 0)
    {
      algo = 2;
      nmax = -nmax;

      /* Loop over divisors */

      ptr = (long)RDB[DATA_PTR_DIV0];
      while (ptr > VALID_PTR)
        {
          /* Check for unsupported sub-division */

          if ((long)RDB[ptr + DIV_NSEG] > 1)
            Error(ptr,
                  "Volume calculation by tracks does not work with sub-division");

          /* next */

          ptr = NextItem(ptr);
        }
    }
  else
    Error(0, "Zero points / tracks in volume calculation routine");

  /* Check if calculation is required */

  if ((nmax == 0) && (tmax < 0.0) && (emax < 0.0))
    return;

  if (algo == 1)
    fprintf(outp, "Calculating material volumes by sampling points...\n");
  else
    fprintf(outp, "Calculating material volumes by sampling tracks...\n");

  /***************************************************************************/

  /***** Allocate memory for statistics **************************************/

  /* Reset stat pointer */

  loc0 = -1;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Allocate memory for volume and density */

      ptr = NewStat("VOLUME", 1, 1);
      WDB[mat + MATERIAL_PTR_MC_VOLUME] = (double)ptr;

      /* Remember pointer to first */

      if (loc0 < 0)
        loc0 = ptr;

      ptr = NewStat("DENSITY", 1, 1);
      WDB[mat + MATERIAL_PTR_MC_DENSITY] = (double)ptr;

      /* Next */

      mat = NextItem(mat);
    }

  /* Loop over cells */

  cell = (long)RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /* Allocate memory for volume and density */

      ptr = NewStat("VOLUME", 1, 1);
      WDB[cell + CELL_PTR_MC_VOLUME] = (double)ptr;

      /* Remember pointer to first */

      if (loc0 < 0)
        loc0 = ptr;

      ptr = NewStat("DENSITY", 1, 1);
      WDB[cell + CELL_PTR_MC_DENSITY] = (double)ptr;

      /* Next */

      cell = NextItem(cell);
    }

  /* Loop over unstructured mesh-based geometries (näiden cell-rakenteet) */
  /* eivät löydy tulta DATA_PTR_C0-alkavasta listasta). */

  umsh = (long)RDB[DATA_PTR_UMSH0];
  while (umsh > VALID_PTR)
    {
      /* Pointer to interface structure */

      ifc = (long)RDB[umsh + UMSH_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

      /* Loop over cells */

      cell = (long)RDB[ifc + IFC_PTR_GCELL_LIST];
      while (cell > VALID_PTR)
        {
          /* Allocate memory for volume and density */

          ptr = NewStat("VOLUME", 1, 1);
          WDB[cell + CELL_PTR_MC_VOLUME] = (double)ptr;

          /* Remember pointer to first */

          if (loc0 < 0)
            loc0 = ptr;

          ptr = NewStat("DENSITY", 1, 1);
          WDB[cell + CELL_PTR_MC_DENSITY] = (double)ptr;

          /* Next cell */

          cell = NextItem(cell);
        }

      /* Next */

      umsh = NextItem(umsh);
    }

  /* Total length */

  ptl = NewStat("LTOT", 1, 1);

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

  ExpandPrivateArrays();

  /***************************************************************************/

  /***** Monte carlo volume calculation **************************************/

  /* Get geometry boundaries */

  xmin = RDB[DATA_GEOM_MINX];
  xmax = RDB[DATA_GEOM_MAXX];
  ymin = RDB[DATA_GEOM_MINY];
  ymax = RDB[DATA_GEOM_MAXY];
  zmin = RDB[DATA_GEOM_MINZ];
  zmax = RDB[DATA_GEOM_MAXZ];

  /* Calculate total volume */

  vol = (xmax - xmin)*(ymax - ymin);

  /* Check dimension */

  if ((dim = (long)RDB[DATA_GEOM_DIM]) == 3)
    {
      /* Add axial dimension to volume */

      vol = vol*(zmax - zmin);
    }

  /* Print dimensions */

  fprintf(outp, "\nx-dimensions : %1.5E %1.5E cm\n", xmin, xmax);
  fprintf(outp, "y-dimensions : %1.5E %1.5E cm\n", ymin, ymax);

  if (dim == 3)
    fprintf(outp, "z-dimensions : %1.5E %1.5E cm\n", zmin, zmax);

  /* Start timer */

  StartTimer(TIMER_VOLUME_CALC);

  /* Avoid compiler warning */

  nt = -1;

  /* Number of samples per batch */

  if (nmax <= 0)
    nt = 100000;
  else if (nmax > 99)
    nt = (long)((double)nmax/100.0);
  else
    Error(0, "Minimum number of points for volume calculation is 100");

  /* Reset counters */

  nb = 0;

  /* Loop over batches */

  while (1 != 2)
    {
      /* Start parallel timer */

      StartTimer(TIMER_OMP_PARA);

      /* Check algorithm */

      if (algo == 1)
        {
          /*******************************************************************/

          /***** Volume calculation by sampling points ***********************/

#ifdef OPEN_MP
#pragma omp parallel private (m, idx, seed, ptr, x, y, z, u, v, w, cell, mat, id, f, T)
#endif
          {
            /* Loop over points */

#ifdef OPEN_MP
#pragma omp for
#endif
            for (m = 0; m < nt; m++)
              {
                /* Get OpenMP thread num */

                id = OMP_THREAD_NUM;

                /* Calculate index */

                idx = nb*nt + m + 1;

                /* Init random number sequence */

                seed = ReInitRNG(idx);
                SEED[id*RNG_SZ] = seed;

                /* Sample point */

                x = RandF(id)*(xmax - xmin) + xmin;
                y = RandF(id)*(ymax - ymin) + ymin;

                if (dim == 3)
                  z = RandF(id)*(zmax - zmin) + zmin;
                else
                  z = 0.0;

                /* Sample direction (this is necessary for STL geometries) */

                IsotropicDirection(&u, &v, &w, id);

                /* Find position */

                if ((cell = WhereAmI(x, y, z, u, v, w, id)) < 0)
                  Error(0, "Geometry error at %E %E %E", x, y, z);

                /* Score point estimator of cell volume */

                ptr = (long)RDB[cell + CELL_PTR_MC_VOLUME];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                AddBuf1D(vol/((double)nt), 1.0, ptr, id, 0);

                /* Check if cell has material */

                if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                  {
                    /* Get material pointer */

                    mat = MatPtr(mat, id);

                    /* Reset density and temperature */

                    f = 1.0;
                    T = 0.0;

                    /* Get point from interface */

                    IFCPoint(mat, &f, &T, -1.0, id);

                    /* Check for undefined density */

                    if (f < 0.0)
                      continue;

                    /* Score point estimators */

                    ptr = (long)RDB[mat + MATERIAL_PTR_MC_VOLUME];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    AddBuf1D(vol/((double)nt), 1.0, ptr, id, 0);

                    ptr = (long)RDB[mat + MATERIAL_PTR_MC_DENSITY];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    AddBuf1D(f*vol/((double)nt), 1.0, ptr, id, 0);

                    ptr = (long)RDB[cell + CELL_PTR_MC_DENSITY];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                    AddBuf1D(f*vol/((double)nt), 1.0, ptr, id, 0);

                    /* Check if material is divided */

                    if ((mat = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT])
                        > VALID_PTR)
                      {
                        /* Add to stat */

                        ptr = (long)RDB[mat + MATERIAL_PTR_MC_VOLUME];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                        AddBuf1D(vol/((double)nt), 1.0, ptr, id, 0);

                        ptr = (long)RDB[mat + MATERIAL_PTR_MC_DENSITY];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                        AddBuf1D(f*vol/((double)nt), 1.0, ptr, id, 0);
                      }
                  }
              }
          }

          /*******************************************************************/
        }
      else if (algo == 2)
        {
          /*******************************************************************/

          /***** Volume calculation by sampling tracks ***********************/

          /* Adjust axial boundaries in 2D geometries */

          if (dim == 2)
            {
              zmin = -INFTY;
              zmax = INFTY;
            }

#ifdef OPEN_MP
#pragma omp parallel private (m, idx, seed, ptr, x, y, z, u, v, w, cell, mat, id, f, T, d, lmax)
#endif
          {
            /* Loop over points */

#ifdef OPEN_MP
#pragma omp for
#endif
            for (m = 0; m < nt; m++)
              {
                /* Get OpenMP thread num */

                id = OMP_THREAD_NUM;

                /* Calculate index */

                idx = nb*nt + m + 1;

                /* Init random number sequence */

                seed = ReInitRNG(idx);
                SEED[id*RNG_SZ] = seed;

                /* Sample face */

                if (dim == 2)
                  n = (long)(RandF(id)*4.0);
                else
                  n = (long)(RandF(id)*6.0);

                /* Avoid compiler warning */

                x = 0.0;
                y = 0.0;
                z = 0.0;
                u = 0.0;
                v = 0.0;
                w = 0.0;
                lmax = 0.0;

                /* Check */

                if (n == 0)
                  {
                    /* West face, sample position */

                    x = xmin + EXTRAP_L;
                    y = RandF(id)*(ymax - ymin) + ymin;

                    if (dim == 3)
                      z = RandF(id)*(zmax - zmin) + zmin;
                    else
                      z = 0.0;

                    /* Set direction */

                    u = 1.0;
                    v = 0.0;
                    w = 0.0;

                    /* Set ray length */

                    lmax = xmax - xmin;
                  }
                else if (n == 1)
                  {
                    /* South face, sample position */

                    x = RandF(id)*(xmax - xmin) + xmin;
                    y = ymin + EXTRAP_L;

                    if (dim == 3)
                      z = RandF(id)*(zmax - zmin) + zmin;
                    else
                      z = 0.0;

                    /* Set direction */

                    u = 0.0;
                    v = 1.0;
                    w = 0.0;

                    /* Set ray length */

                    lmax = ymax - ymin;
                  }
                else if (n == 2)
                  {
                    /* East face, sample position */

                    x = xmax - EXTRAP_L;
                    y = RandF(id)*(ymax - ymin) + ymin;

                    if (dim == 3)
                      z = RandF(id)*(zmax - zmin) + zmin;
                    else
                      z = 0.0;

                    /* Set direction */

                    u = -1.0;
                    v = 0.0;
                    w = 0.0;

                    /* Set ray length */

                    lmax = xmax - xmin;
                  }
                else if (n == 3)
                  {
                    /* North face, sample position */

                    x = RandF(id)*(xmax - xmin) + xmin;
                    y = ymax - EXTRAP_L;

                    if (dim == 3)
                      z = RandF(id)*(zmax - zmin) + zmin;
                    else
                      z = 0.0;

                    /* Set direction */

                    u = 0.0;
                    v = -1.0;
                    w = 0.0;

                    /* Set ray length */

                    lmax = ymax - ymin;
                  }
                else if (n == 4)
                  {
                    /* Bottom face, sample position */

                    x = RandF(id)*(xmax - xmin) + xmin;
                    y = RandF(id)*(ymax - ymin) + ymin;
                    z = zmin + EXTRAP_L;

                    /* Set direction */

                    u = 0.0;
                    v = 0.0;
                    w = 1.0;

                    /* Set ray length */

                    lmax = zmax - zmin;
                  }
                else if (n == 5)
                  {
                    /* Top face, sample position */

                    x = RandF(id)*(xmax - xmin) + xmin;
                    y = RandF(id)*(ymax - ymin) + ymin;
                    z = zmax - EXTRAP_L;

                    /* Set direction */

                    u = 0.0;
                    v = 0.0;
                    w = -1.0;

                    /* Set ray length */

                    lmax = zmax - zmin;
                  }
                else
                  Die(FUNCTION_NAME, "Error");

                /* Score total length */

                AddBuf1D(lmax, 1.0, ptl, id, 0);

                /* Loop until outside */

                do
                  {
                    /* Get cell */

                    if ((cell = WhereAmI(x, y, z, u, v, w, id)) < 0)
                      Error(0, "Geometry error at %E %E %E", x, y, z);

                    /* Calculate distance to boundary */

                    d = NearestBoundary(id);

                    /* Score cell volume */

                    if ((long)RDB[cell + CELL_TYPE] != CELL_TYPE_OUTSIDE)
                      {
                        ptr = (long)RDB[cell + CELL_PTR_MC_VOLUME];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                        AddBuf1D(vol*d, 1.0, ptr, id, 0);
                      }

                    /* Check if cell has material */

                    if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                      {
                        /* Get material pointer */

                        mat = MatPtr(mat, id);

                        /* Score volume */

                        ptr = (long)RDB[mat + MATERIAL_PTR_MC_VOLUME];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                        AddBuf1D(vol*d, 1.0, ptr, id, 0);

                        /* Check if material is divided */

                        if ((mat = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT])
                            > VALID_PTR)
                          {
                            /* Add to stat */

                            ptr = (long)RDB[mat + MATERIAL_PTR_MC_VOLUME];
                            CheckPointer(FUNCTION_NAME, "(ptr)",
                                         DATA_ARRAY, ptr);
                            AddBuf1D(vol*d, 1.0, ptr, id, 0);
                          }
                      }

                    /* Move to new position */

                    x = x + (d + EXTRAP_L)*u;
                    y = y + (d + EXTRAP_L)*v;
                    z = z + (d + EXTRAP_L)*w;
                  }
                while ((x > xmin) && (x < xmax) && (y > ymin) && (y < ymax) &&
                       (z > zmin) && (z < zmax));


              }
          }

          /*******************************************************************/
        }
      else
        Die(FUNCTION_NAME, "Error in algorithm");

      /* Stop parallel timer */

      StopTimer(TIMER_OMP_PARA);

      /* Update batch number */

      nb++;

      /* Reduce scoring buffer */

      ReduceBuffer();

      /* Maximum error */

      max = 0.0;

      /* Loop over statistics */

      ptr = loc0;
      while (ptr > VALID_PTR)
        {
          /* Collect data */

          val = BufVal(ptr, 0);
          div = BufVal(ptl, 0);

          if (algo == 1)
            AddStat(val, ptr, 0);
          else
            AddStat(val/div, ptr, 0);

          /* Get relative error */

          err = RelErr(ptr, 0);

          /* Compare maximum error */

          if (err > max)
            max = err;

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Clear buffer */

      ClearBuf();

      /* Get time */

      t = TimerVal(TIMER_VOLUME_CALC);

      if (nb == 1)
        fprintf(outp, "\nEstimated calculation time: %s\n",
                TimeStr((long)((double)nmax*t/((double)nt))));

      /* Check batch cut-off */

      if ((nmax > 0) && (nb*nt >= nmax))
        {
          fprintf(outp, "Realized calculation time:  %s\n\n",
                  TimeStr((long)t));
          break;
        }

      /* Check time cut-off */

      if ((tmax > 0.0) && (t > tmax))
        break;

      /* Check error cut-off */

      if ((emax > 0.0) && (max > 0.0) && (max < emax))
        break;
    }

  /* Stop timer */

  StopTimer(TIMER_VOLUME_CALC);

  /***************************************************************************/

  /***** Print volumes *******************************************************/

  /* Check dimensions */

  if ((long)RDB[DATA_VOLUME_CALCULATION_MODE] == YES)
    {
      if (dim == 3)
        fprintf(outp, "Volumes (in cm3) :\n\n");
      else
        fprintf(outp, "Volumes (2D problem, the values are in cm2) :\n\n");
    }

  /* Check mode */

  if (((long)RDB[DATA_VOLUME_CALCULATION_MODE] == YES) && (mode == 1))
    {
      /* Set file name */

      sprintf(fname, "%s.mvol", GetText(DATA_PTR_INPUT_FNAME));

      /* Open file */

      if ((fp = fopen(fname, "w")) == NULL)
        Die(FUNCTION_NAME, "Unable to open file for writing");

      /* Print comments and card name */

      if (mode == 1)
        fprintf(fp, "%% --- Material volumes:\n\n");
      else
        fprintf(fp, "%% --- Cell volumes:\n\n");

      fprintf(fp, "%% Produced %s by MC volume calculation routine by\n",
              GetText(DATA_PTR_DATE));

      if (algo == 1)
        fprintf(fp, "%% sampling %ld random points in the geometry.\n\n",
                nb*nt);
      else
        fprintf(fp, "%% sampling %ld random tracks in the geometry.\n\n",
                nb*nt);
      fprintf(fp, "set mvol\n\n");
    }
  else
    fp = NULL;

  /* Check mode */

  if (mode == 1)
    {
      /***********************************************************************/

      /***** Material volumes ************************************************/

      /* Get maximum string length */

      nmax = -1;

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Compare string length */

          if ((n = strlen(GetText(mat + MATERIAL_PTR_NAME))) > nmax)
            nmax = n;

          /* Next material */

          mat = NextItem(mat);
        }

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Get density */

          ptr = (long)RDB[mat + MATERIAL_PTR_MC_DENSITY];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          df = Mean(ptr, 0);

          /* Pointer to volume */

          ptr = (long)RDB[mat + MATERIAL_PTR_MC_VOLUME];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Given or calculated volume */

          vol = RDB[mat + MATERIAL_VOLUME];

          /* Estimated volume, error and relative difference */

          est = Mean(ptr, 0);
          err = RelErr(ptr, 0);
          diff = est/vol - 1.0;

          /* Pointer to parent */

          mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT];

          /* Get number of zones */

          m = 0;
          if (mat0 > VALID_PTR)
            m = (long)RDB[mat0 + MATERIAL_DIV_N_TOT_ZONES];

          /* Check mode */

          if ((long)RDB[DATA_VOLUME_CALCULATION_MODE] == YES)
            {
              /* Print if no or more than 1 zones */

              if (m != 1)
                {
                  /* Get printed name */

                  sprintf(fmt, "Material %%-%lds", nmax);

                  sprintf(tmpstr, fmt,
                          GetText(mat + MATERIAL_PTR_NAME));

                  /* Print */

                  if ((vol > ZERO) && (vol < INFTY))
                    {
                      if (diff > 1.0)
                        fprintf(outp,
                                "%s : %1.4E %1.4E (%7.5f) :    > 1.0 ",
                                tmpstr, vol, est, err);
                      else if (diff < -1.0)
                        fprintf(outp, "%s : %1.4E %1.4E (%7.5f) :  < -1.0 ",
                                tmpstr, vol, est, err);
                      else
                        fprintf(outp, "%s : %1.4E %1.4E (%7.5f) : %8.5f ",
                                tmpstr, vol, est, err, diff);
                    }
                  else if (vol == 0.0)
                    fprintf(outp, "%s : %1.4E %1.4E (%7.5f) :      N/A ",
                            tmpstr, vol, est, err);
                  else
                    fprintf(outp,
                            "%s :        N/A %1.4E (%7.5f) :      N/A ",
                            tmpstr, est, err);

                  if (fabs(diff) > 1.96*err)
                    fprintf(outp, "* ");
                  else
                    fprintf(outp, "  ");

                  if ((algo == 1) && (est > 0.0) && ((long)RDB[DATA_PTR_IFC0] > VALID_PTR))
                    fprintf(outp, "(%5.1f %% den.)\n", 100*df/est);
                  else
                    fprintf(outp, "\n");
                }

              /* Print volume to file */

              if ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
                {
                  if ((mat0 =
                       (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
                fprintf(fp, "%-10s ", GetText(mat0 + MATERIAL_PTR_NAME));
                  else
                    fprintf(fp, "%-10s ", GetText(mat + MATERIAL_PTR_NAME));

                  fprintf(fp, "%6ld %1.5E %% (%5.3f)\n",
                          (long)RDB[mat + MATERIAL_DIV_ZONE_IDX], Mean(ptr, 0),
                          RelErr(ptr, 0));
                }
            }
          else
            {
              /* Put volume if not given */

              if (RDB[mat + MATERIAL_VOLUME_GIVEN] < 0.0)
                WDB[mat + MATERIAL_VOLUME] = Mean(ptr, 0);
            }

          /* Next */

          mat = NextItem(mat);
        }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Cell volumes ****************************************************/

      /* Get maximum string length */

      nmax = -1;

      cell = (long)RDB[DATA_PTR_C0];
      while (cell > VALID_PTR)
        {
          /* Compare string length */

          if ((n = strlen(GetText(cell + CELL_PTR_NAME))) > nmax)
            nmax = n;

          /* Next cell */

          cell = NextItem(cell);
        }

      /* Close file if open */

      if (fp != NULL)
        fclose(fp);

      /* Write in matlab format file */

      sprintf(fname, "%s_vol.m", GetText(DATA_PTR_INPUT_FNAME));

      /* Open file */

      if ((fp = fopen(fname, "w")) == NULL)
        Die(FUNCTION_NAME, "Unable to open file for writing");

      fprintf(fp, "vol = [\n");

      /* Loop over cells */

      cell = (long)RDB[DATA_PTR_C0];
      while (cell > VALID_PTR)
        {
          /* Get density */

          ptr = (long)RDB[cell + CELL_PTR_MC_DENSITY];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          df = Mean(ptr, 0);

          /* Pointer to volume */

          ptr = (long)RDB[cell + CELL_PTR_MC_VOLUME];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Given or calculated volume */

          vol = RDB[cell + CELL_VOLUME];

          /* Estimated volume, error and relative difference */

          est = Mean(ptr, 0);
          err = RelErr(ptr, 0);
          diff = est/vol - 1.0;

          /* Get printed name */

          sprintf(fmt, "Cell %%-%lds", nmax);
          sprintf(tmpstr, fmt,
                  GetText(cell + CELL_PTR_NAME));

          /* Print */

          if ((vol > ZERO) && (vol < INFTY))
            {
              if (diff > 1.0)
                fprintf(outp,
                        "%s : %1.4E %1.4E (%7.5f) :    > 1.0 ",
                        tmpstr, vol, est, err);
              else if (diff < -1.0)
                fprintf(outp, "%s : %1.4E %1.4E (%7.5f) :  < -1.0 ",
                        tmpstr, vol, est, err);
              else
                fprintf(outp, "%s : %1.4E %1.4E (%7.5f) : %8.5f ",
                        tmpstr, vol, est, err, diff);
            }
          else if (vol == 0.0)
            fprintf(outp, "%s : %1.4E %1.4E (%7.5f) :      N/A ",
                    tmpstr, vol, est, err);
          else
            fprintf(outp,
                    "%s :        N/A %1.4E (%7.5f) :      N/A ",
                    tmpstr, est, err);

          if (fabs(diff) > 1.96*err)
            fprintf(outp, "* ");
          else
            fprintf(outp, "  ");

          if ((algo == 1) && (est > 0.0))
            fprintf(outp, "(%5.1f %% den.)\n", 100*df/est);
          else
            fprintf(outp, "\n");

          /* Print to file */

          fprintf(fp, "%s %1.5E %7.5f\n", GetText(cell + CELL_PTR_NAME),
                  est, err);

          /* Next */

          cell = NextItem(cell);
        }

      fprintf(fp, "];\n\n");

      /***********************************************************************/
    }

  /* Close file */

  if (fp != NULL)
    fclose(fp);

  /***************************************************************************/

  if (((long)RDB[DATA_VOLUME_CALCULATION_MODE] == YES) && (mode == 1))
    fprintf(outp, "\nVolumes written in file \"%s\"\n", fname);

  /* Exit subroutine */

  fprintf(outp, "\nOK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
