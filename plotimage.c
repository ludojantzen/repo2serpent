/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : plotimage.c                                    */
/*                                                                           */
/* Created:       2020/06/17 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Produces a matrix of color codes used with geometry plotter. */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PlotImage:"

/*****************************************************************************/

long PlotImage(long gpl, long *mtx1, long xp, long yp, double xmin,
               double xmax, double ymin, double ymax, double zmin,
               double zmax, double pos, long bou, long par, long ax,
               long scale, double E, double fmin, double fmax, long qp)
{
#ifndef NO_GFX_MODE

  long old, cell, n, m, mat, ptr, id, nerr, wwp, i, *mtx2, nifc, nc0, ncol;
  long lin;
  double fmean, T, d, x, y, z, u, v, w, *spt, f, xt, yt, zt, dummy, pixw;

  /* Plot weight window boundaries */

  if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == YES)
    wwp = YES;
  else
    wwp = NO;

  /* Override for now */

  wwp = NO;

  /* Allocate memory for pointer matrix */

  mtx2 = (long *)Mem(MEM_ALLOC, xp*yp, sizeof(long));

  /* Allocate memory for source points */

  if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
    spt = (double *)Mem(MEM_ALLOC, xp*yp, sizeof(double));
  else
    spt = NULL;

  /* Read some values (nää kirjoitetaan putplotcolors.c:ssä)  */

  nifc = (long)RDB[DATA_PLOTTER_NIFC];
  ncol = (long)RDB[DATA_PLOTTER_NCOL];
  nc0 = (long)RDB[DATA_PLOTTER_NC0];

  /* Put quick plotter mode */

  WDB[DATA_QUICK_PLOT_MODE] = (double)qp;

  /* Set importance plot type */

  if ((scale == 1) || (scale == 3))
    lin = YES;
  else
    lin = NO;

  if (scale > 2)
    scale = 2;
  else if (scale > 0)
    scale = 1;

  /***************************************************************************/

  /***** Fill plot matrix with data ******************************************/

  /* Start parallel timer */

  StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private(n, m, x, y, z, u, v, w, cell, mat, ptr, f, T, id)
#endif
  {
    /* Get Open MP thread id */

    id = OMP_THREAD_NUM;

    /* Avoid compiler warnings */

    x = 0;
    y = 0;
    z = 0;
    u = 0;
    v = 0;
    w = 0;

    /* Loop over geometry */

#ifdef OPEN_MP
#pragma omp for
#endif
    for (n = 0; n < xp; n++)
      {
        for (m = 0; m < yp; m++)
          {
            /* Reset pixel value */

            mtx1[n*yp + m] = 0;
            mtx2[n*yp + m] = 0;

            /*****************************************************************/

            /***** Calculate Coordinates **************************************/

            switch (ax)
              {
              case PLOT_MODE_YZ:
                {
                  /* yz-plot */

                  if (pos > -INFTY)
                    x = pos;
                  else
                    x = (xmax - xmin)/2.0 + xmin;

                  y = (n/(xp - 1.0))*(ymax - ymin) + ymin;
                  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;

                  u = 0.0;
                  v = 0.0;
                  w = 1.0;

                  break;
                }
              case PLOT_MODE_XZ:
                {
                  /* xz-plot */

                  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;

                  if (pos > -INFTY)
                    y = pos;
                  else
                    y = (ymax - ymin)/2.0 + ymin;

                  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;

                  u = 0.0;
                  v = 0.0;
                  w = 1.0;

                  break;
                }
              case PLOT_MODE_XY:
                {

                  /* xy-plot */

                  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;
                  y = (m/(yp - 1.0))*(ymax - ymin) + ymin;

                  if (pos > -INFTY)
                    z = pos;
                  else
                    z = (zmax - zmin)/2.0 + zmin;

                  u = 0.0;
                  v = 1.0;
                  w = 0.0;

                  break;
                }
              default:
                Die(FUNCTION_NAME, "Invalid plot mode");
              }

            /* Check true symmetry */

            if (TestUniSym(x, y, z) == YES)
              continue;

            /*****************************************************************/

            /***** Handle importances ****************************************/

            /* Check if importances are plotted */

            if (scale > 0)
              {
                /* Check scale */

                if (scale == 1)
                  {
                    /* Get cell importance */

                    if (par > 0)
                      f = WWImportance(par, x, y, z,
                                       u, v, w, E, WWMESH_COLL);
                    else if (E > 0.0)
                      f = WWImportance(PARTICLE_TYPE_NEUTRON, x, y, z,
                                       u, v, w, E, WWMESH_COLL);
                    else
                      f = WWImportance(PARTICLE_TYPE_GAMMA, x, y, z,
                                       u, v, w, -E, WWMESH_COLL);
                  }
                else
                  {
                    /* Get source importance */

                    if (par > 0)
                      f = WWImportance(par, x, y, z,
                                       u, v, w, E, WWMESH_SRC);

                    else if (E > 0.0)
                      f = WWImportance(PARTICLE_TYPE_NEUTRON, x, y, z,
                                       u, v, w, E, WWMESH_SRC);
                    else
                      f = WWImportance(PARTICLE_TYPE_GAMMA, x, y, z,
                                       u, v, w, -E, WWMESH_SRC);
                  }

                /* Check if cell importances are given */

                if (f <= 0.0)
                  {
                    if ((cell = WhereAmI(x, y, z, u, v, w, id)) >
                        VALID_PTR)
                      BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w,
                                         &xt, &yt, &zt, &dummy, id);
                    if (cell > VALID_PTR)
                      f = RDB[cell + CELL_IMP];
                  }

                /* If no value provided, set to 1.0 */

                if (f <= 0.0)
                  f = 1.0;

                /* Interpolate */

                if (lin == YES)
                  f = (f - fmin)/(fmax - fmin);
                else
                  f = (log(f) - log(fmin))/(log(fmax) - log(fmin));

                /* Check */

                if (f < 0.0)
                  f = 0.0;
                else if (f > 1.0)
                  f = 1.0;

                /* Put index */

                mtx1[n*yp + m] = (long)(f*(ncol - nc0 - 1)) + nc0;
              }

            /*****************************************************************/

            /***** Standard geometry plot ************************************/

            /* Find location */

            if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
              BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w,
                                 &xt, &yt, &zt, &dummy, id);

            /* Check return value */

            if (cell < 0)
              {
                /* Put error condition */

                nerr = cell;

                /* Set colour */

                mtx1[n*yp + m] = -cell;
                mtx2[n*yp + m] = -cell;
              }
            else if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
              {
                /* Get pointer */

                mat = MatPtr(mat, id);
                CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

                /* Set colour */

                if (mtx1[n*yp + m] == 0)
                  mtx1[n*yp + m] = (long)RDB[mat + MATERIAL_COLOUR_IDX];

                mtx2[n*yp + m] = (long)mat;

                /* Reset density and temperature */

                f = 1.0;
                T = 0.0;

                /* Get point from interface */

                IFCPoint(mat, &f, &T, -1.0, id);

                /* Check for undefined region */

                if (f < 0.0)
                  T = -1.0;
                else if ((ptr = (long)RDB[mat + MATERIAL_PTR_IFC]) >
                         VALID_PTR)
                  {
                    /* Scaling */

                    /* JLe 26.8.2015: tää feilaa kun MATERIAL_ADENS */
                    /* ei enää rajoiteta interfacen maksimiin       */
                    /* version 2.1.24 jälkeen. */

                    /*
                      if (RDB[ptr + IFC_MAX_DENSITY] !=
                      RDB[ptr + IFC_MIN_DENSITY])
                      f = (f*RDB[mat + MATERIAL_ADENS] -
                      RDB[ptr + IFC_MIN_DENSITY])/
                      (RDB[ptr + IFC_MAX_DENSITY] -
                      RDB[ptr + IFC_MIN_DENSITY]);
                      else
                      f = 1.0;
                    */

                    if (RDB[ptr + IFC_MAX_DENSITY] ==
                        RDB[ptr + IFC_MIN_DENSITY])
                      f = 1.0;
                  }

                /* Check temperature (overrides density) */

                if (T > 0.0)
                  {
                    /* Calculate factor */

                    if (RDB[mat + MATERIAL_TMS_TMAX] -
                        RDB[mat + MATERIAL_TMS_TMIN] > 0.0)
                      f = 0.7*(T - RDB[mat + MATERIAL_TMS_TMIN])/
                        (RDB[mat + MATERIAL_TMS_TMAX] -
                         RDB[mat + MATERIAL_TMS_TMIN]) + 0.3;
                    else
                      f = 1.0;
                  }

                /* Check */

                if ((f < 0.0) || (f > 1.0))
                  mtx1[n*yp + m] = 4;
                else
                  mtx1[n*yp + m] = mtx1[n*yp + m] +
                    (long)((1.0 - f)*(nifc - 1));
              }

            /*****************************************************************/

            /***** Plot super-imposed source / detector region ***************/

#ifdef mmmmmmmm

            ptr = (long)RDB[DATA_PTR_SRC0];
            while (ptr > VALID_PTR)
              {
                if (InSuperCell((long)RDB[ptr + SRC_PTR_UNIV],
                                (long)RDB[ptr + SRC_PTR_CELL], x, y, z)
                    == YES)
                  mtx1[n*yp + m] = 3;

                /* Next */

                ptr = NextItem(ptr);
              }
#endif
            /*****************************************************************/
          }
      }
  }

  /* Stop parallel timer */

  StopTimer(TIMER_OMP_PARA);

  /***************************************************************************/

  /***** Add boundaries  *****************************************************/

  /* Start parallel timer */

  StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private(n, m, x, y, z, d, cell, mat, pixw, old, ptr, id)
#endif
  {
    /* Get Open MP thread id */

    id = OMP_THREAD_NUM;

    /* Avoid compiler warnings */

    pixw = 0;
    old = -1;

    /*************************************************************************/

    /***** Sweep 1 ***********************************************************/

#ifdef OPEN_MP
#pragma omp for
#endif

    for (n = 0; n < xp; n++)
      {
        /* Reset minimum distance */

        d = -1.0;

        for (m = 0; m < yp; m++)
          {
            /* Calculate Co-ordinates */

            switch (ax)
              {
              case PLOT_MODE_YZ:
                {
                  /* yz-plot */

                  if (pos > -INFTY)
                    x = pos;
                  else
                    x = (xmax - xmin)/2.0 + xmin;

                  y = (n/(xp - 1.0))*(ymax - ymin) + ymin;
                  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;

                  pixw = (zmax - zmin)/((double)yp);

                  u = 0.0;
                  v = 0.0;
                  w = 1.0;

                  break;
                }
              case PLOT_MODE_XZ:
                {
                  /* xz-plot */

                  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;

                  if (pos > -INFTY)
                    y = pos;
                  else
                    y = (ymax - ymin)/2.0 + ymin;

                  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;

                  pixw = (zmax - zmin)/((double)yp);

                  u = 0.0;
                  v = 0.0;
                  w = 1.0;

                  break;
                }
              case PLOT_MODE_XY:
                {
                  /* xy-plot */

                  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;
                  y = (m/(yp - 1.0))*(ymax - ymin) + ymin;

                  if (pos > -INFTY)
                    z = pos;
                  else
                    z = (zmax - zmin)/2.0 + zmin;

                  pixw = (ymax - ymin)/((double)yp);

                  u = 0.0;
                  v = 1.0;
                  w = 0.0;

                  break;
                }
              default:
                Die(FUNCTION_NAME, "Invalid plot mode");
              }

            /* Check true symmetry */

            if (TestUniSym(x, y, z) == YES)
              continue;

            /* Weight window boundary */

            if (wwp == YES)
              if (WWDis(par, x, y, z, u, v, w) < pixw)
                mtx1[n*yp + m] = 0;

            /* Check mode */

            if ((bou == 1) || (bou == 3))
              {
                /* Find location */

                if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
                  BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w,
                                     &xt, &yt, &zt, &dummy, id);

                /* Check pointer */

                if (cell > VALID_PTR)
                  {
                    /* Calculate distance to boundary */

                    d = NearestBoundary(id);

                    /* Compare minimum distance to pixel width */

                    if (d < pixw)
                      {
                        /* Set colour */

                        mtx1[n*yp + m] = 0;
                      }
                  }
              }

            if ((bou == 2) || (bou == 3))
              {
                /* Material index */

                mat = mtx2[n*yp + m];

                /* Compare to previous */

                if (old != mat)
                  {
                    /* Put colour */
                    /*
                      if ((old < nc0) || (mat < nc0) ||
                      ((old > nc0 - 1) && (mat > nc0 - 1) &&
                      (long)((double)(old - nc0)/((double)nifc)) !=
                      (long)((double)(mat - nc0)/((double)nifc))))
                    */
                    if (old != 0)
                      mtx1[n*yp + m] = 0;

                    /* Set previous pointer */

                    old = mat;
                  }
              }

            /* Draw boundary */

            if ((m == 0) || (m == yp - 1))
              mtx1[n*yp + m] = 0;
          }
      }

    /*************************************************************************/

    /***** Sweep 2 ***********************************************************/

#ifdef OPEN_MP
#pragma omp for
#endif
    for (m = 0; m < yp; m++)
      {
        /* Reset minimum distance */

        d = -1.0;

        for (n = 0; n < xp; n++)
          {
            /* Calculate Co-ordinates */

            switch (ax)
              {
              case PLOT_MODE_YZ:
                {
                  /* yz-plot */

                  if (pos > -INFTY)
                    x = pos;
                  else
                    x = (xmax - xmin)/2.0 + xmin;

                  y = (n/(xp - 1.0))*(ymax - ymin) + ymin;
                  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;

                  pixw = (ymax - ymin)/((double)yp);

                  u = 0.0;
                  v = 1.0;
                  w = 0.0;

                  break;
                }
              case PLOT_MODE_XZ:
                {
                  /* xz-plot */

                  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;

                  if (pos > -INFTY)
                    y = pos;
                  else
                    y = (ymax - ymin)/2.0 + ymin;

                  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;

                  pixw = (xmax - xmin)/((double)xp);

                  u = 1.0;
                  v = 0.0;
                  w = 0.0;

                  break;
                }
              case PLOT_MODE_XY:
                {

                  /* xy-plot */

                  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;
                  y = (m/(yp - 1.0))*(ymax - ymin) + ymin;

                  if (pos > -INFTY)
                    z = pos;
                  else
                    z = (zmax - zmin)/2.0 + zmin;

                  pixw = (xmax - xmin)/((double)xp);

                  u = 1.0;
                  v = 0.0;
                  w = 0.0;

                  break;
                }
              default:
                Die(FUNCTION_NAME, "Invalid plot mode");
              }

            /* Weight window boundary */

            if (wwp == YES)
              if (WWDis(par, x, y, z, u, v, w) < pixw)
                mtx1[n*yp + m] = 0;

            /* Check mode */

            if ((bou == 1) || (bou == 3))
              {
                /* Find location */

                if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
                  BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w,
                                     &xt, &yt, &zt, &dummy, id);

                /* Check pointer */

                if (cell > VALID_PTR)
                  {
                    /* Calculate distance to boundary */

                    d = NearestBoundary(id);

                    /* Compare minimum distance to pixel width */

                    if (d < pixw)
                      {
                        /* Set colour */

                        mtx1[n*yp + m] = 0;
                      }
                  }
              }

            if ((bou == 2) || (bou == 3))
              {
                /* Material index */

                mat = mtx2[n*yp + m];

                /* Compare to previous */

                if (old != mat)
                  {
                    /* Avoid double lines */
                    /*
                      if ((old < nc0) || (mat < nc0) ||
                      ((old > nc0 - 1) && (mat > nc0 - 1) &&
                      (long)((double)(old - nc0)/((double)nifc)) !=
                      (long)((double)(mat - nc0)/((double)nifc))))
                    */

                    if ((old != 0) && (n > 0) &&
                        (mtx1[(n - 1)*yp + m] != 0))
                      mtx1[n*yp + m] = 0;

                    /* Set previous pointer */

                    old = mat;
                  }
              }

            /* Draw boundary */

            if ((n == 0) || (n == xp - 1))
              mtx1[n*yp + m] = 0;
          }
      }
  }

  /* Stop parallel timer */

  StopTimer(TIMER_OMP_PARA);

  /***************************************************************************/

  /***** Plot source point distribution **************************************/

  /* Check pointer to source distribution */

  if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
    {
      /* Reset distribution */

      for (n = 0; n < xp; n++)
        for (m = 0; m < yp; m++)
          spt[n*yp + m] = 0.0;

      /* Get pointer */

      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to first after dummy */

      ptr = NextItem(ptr);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over source */

      while (ptr > VALID_PTR)
        {
          /* Get coordinates */

          x = RDB[ptr + PARTICLE_X];
          y = RDB[ptr + PARTICLE_Y];
          z = RDB[ptr + PARTICLE_Z];

          /* Avoid compiler warning */

          n = 0;
          m = 0;

          /* Check plot type */

          switch (ax)
            {
            case PLOT_MODE_YZ:
              {
                /* yz-plot */

                n = xp*(y - ymin)/(ymax - ymin);
                m = yp*(z - zmin)/(zmax - zmin);

                break;
              }
            case PLOT_MODE_XZ:
              {
                /* xz-plot */

                n = xp*(x - xmin)/(xmax - xmin);
                m = yp*(z - zmin)/(zmax - zmin);

                break;
              }
            case PLOT_MODE_XY:
              {
                /* yz-plot */

                n = xp*(x - xmin)/(xmax - xmin);
                m = yp*(y - ymin)/(ymax - ymin);

                break;
              }
            default:
              Die(FUNCTION_NAME, "Invalid plot mode");
            }

          /* Put point */

          if ((n > -1) && (n < xp) && (m > -1) && (m < yp))
            spt[n*yp + m] = spt[n*yp + m] + 1.0;

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Calculate mean */

      fmean = 0.0;
      i = 0;

      for (n = 0; n < xp; n++)
        for (m = 0; m < yp; m++)
          if (spt[n*yp + m] > 0.0)
            {
              fmean = fmean + spt[n*yp + m];
              i++;
            }

      if (i > 0)
        fmean = fmean/((double)i);

      /* distribution minimum and maximum */

      fmin = 0.5*fmean*(1.0 - RDB[DATA_SOURCE_PT_ANIM_F]);
      fmax = 0.5*fmean*(1.0 + RDB[DATA_SOURCE_PT_ANIM_F]);

      /* Put colours */

      if (fmean > 0)
        for (n = 0; n < xp; n++)
          for (m = 0; m < yp; m++)
            if (spt[n*yp + m] > 0.0)
              {
                /* Calculate factor */

                f = (spt[n*yp + m] - fmin)/(fmax - fmin);

                /* Adjust */

                if (f < 0.0)
                  f = 0.0;
                else if (f > 1.0)
                  f = 1.0;

                /* Put color */

                mtx1[n*yp + m] = (long)(f*((double)ncol - 1.0)) + ncol;
              }
    }

  /***************************************************************************/

  /* Free source point distribution */

  if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
    Mem(MEM_FREE, spt);

  /* Free memory */

  Mem(MEM_FREE, mtx2);

  /* Return error code */

  return nerr;

#endif

  /***************************************************************************/
}

/*****************************************************************************/
