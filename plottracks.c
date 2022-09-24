/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : plottracks.c                                   */
/*                                                                           */
/* Created:       2014/10/08 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Plots particle tracks                                        */
/*                                                                           */
/* Comments: - T‰t‰ sallitaan vissiin k‰ytrtt‰v‰n OpenMP-luupissa, mutta noi */
/*             piirto-operaatiot ei oo thread-safeja.                        */
/*                                                                           */
/*           - Toroidimuunnos puuttuu viel‰                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PlotTracks:"

#ifndef NO_GFX_MODE

/* Local subroutines */

void CoordPt(long, long *, long *, long, long, double, double, double, double,
             double, double, double, double, double);

void DrawTail(long, gdImagePtr, long, double **, long, long, long, double,
              double, double, double, double, double);

/*****************************************************************************/

#define MAX_IDX 10000

long PlotTracks(long gpl, gdImagePtr im, long nf, double xmin, double xmax,
                double ymin, double ymax, double zmin, double zmax, long xp,
                long yp)
{
  long nftot, id, evn, n, m, n0, m0, part, type, type0, i, idx, fiss, cyl;
  long pop, ptype, ncol;
  double t0, t, x0, y0, z0, x1, y1, z1, lmax, t1 ,f, r, E, spd;
  double **R;
  /*
  FILE *fp;
  */

  /* Read number of colors (t‰‰ kirjoitetaan putplotcolors.c:ss‰)  */

  ncol = (long)RDB[DATA_PLOTTER_NCOL];

  /* Cylindrical symmetry */

  cyl = NO;

  /* Get number of frames */

  if ((nftot = (long)RDB[DATA_TRACK_PLOT_FRAMES]) < 2)
    {
      /***********************************************************************/

      /***** No frames, plot entire track ************************************/

      /* Loop over threads and bank (NOTE: type argument is not used in  */
      /* FromTrkBank() at the moment, so this works for photons as well) */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        while ((part = FromTrkBank(id)) > VALID_PTR)
          {
            /* Put particle back to bank */

            ToBank(part, id);

            /* Reset indexes */

            n  = -1;
            m  = -1;
            type0  = -1;

            /* Loop over events */

            evn = (long)RDB[part + PARTICLE_PTR_EVENTS];
            while (evn > VALID_PTR)
              {
                /* Get type */

                type = (long)RDB[evn + EVENT_TYPE];

                /* Skip all but collisions, start, leak, bc and cut-offs */

                if ((type != TRACK_END_COLL) &&
                    (type != TRACK_END_STRT) &&
                    (type != TRACK_END_LEAK) &&
                    (type != TRACK_END_TCUT) &&
                    (type != TRACK_END_ECUT) &&
                    (type != TRACK_END_WCUT) &&
                    (type != TRACK_END_SURF) &&
                    (type != TRACK_END_WWIN) &&
                    (type != TRACK_END_BC))
                  {
                    /* Pointer to next */

                    evn = NextItem(evn);

                    /* Cycle loop */

                    continue;
                  }

                /* Get coordinates */

                x0 = RDB[evn + EVENT_X];
                CheckValue(FUNCTION_NAME, "x0", "", x0, -INFTY, INFTY);

                y0 = RDB[evn + EVENT_Y];
                CheckValue(FUNCTION_NAME, "y0", "", y0, -INFTY, INFTY);

                z0 = RDB[evn + EVENT_Z];
                CheckValue(FUNCTION_NAME, "z0", "", z0, -INFTY, INFTY);

                /* Conversion to cylindrical coordinates */

                if (((long)RDB[gpl + GPL_GEOM] != PLOT_MODE_XY) &&
                    (cyl == YES))
                  {
                    r = sqrt(x0*x0 + y0*y0);
                    x0 = r;
                    y0 = r;
                  }

                /* Get position */

                n0 = n;
                m0 = m;

                CoordPt(gpl, &n, &m, xp, yp, x0, y0, z0, xmin, xmax, ymin,
                        ymax, zmin, zmax);

                /* Check periodic */

                if ((type0 != TRACK_END_BC) || (type != TRACK_END_SURF))
                  {
                    /* Draw line (HUOM: t‰ss‰ tulee ongelma OMP:n kanssa */
                    /* kun useamman threadin piirt‰m‰t viivat riste‰‰)   */

                    if (type0 > -1)
                      gdImageLine(im, n0, m0, n, m, 0);
                  }

                if (type == TRACK_END_WWIN)
                  {
                    gdImageSetPixel(im, n-1, m, 2);
                    gdImageSetPixel(im, n+1, m, 2);
                    gdImageSetPixel(im, n, m-1, 2);
                    gdImageSetPixel(im, n, m+1, 2);
                  }

                /* Store previous type */

                type0 = type;

                /* Next event */

                evn = NextItem(evn);
              }
          }

      /* Flush bank to track plot bank */

      FlushBank();

      /* Exit */

      return NO;

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Track plot animation ************************************************/

  /* Allocate memory for point array */

  R = (double **)Mem(MEM_ALLOC, 4, sizeof(double *));
  for(n = 0; n < 5; n++)
    R[n] = (double *)Mem(MEM_ALLOC, MAX_IDX + 1, sizeof(double));

  /* Get number of frames and history lenght */

  nftot = (long)RDB[DATA_TRACK_PLOT_FRAMES];
  lmax = RDB[DATA_TRACK_PLOT_HIS_LENGTH];

  /* Reset population count */

  pop = 0;

  /* Calculate end-of-interval time */

  t = (RDB[DATA_TRACK_PLOT_TMAX] - RDB[DATA_TRACK_PLOT_TMIN])*
    ((double)(nftot - nf - 1)/((double)nftot - 1.0)) +
    RDB[DATA_TRACK_PLOT_TMIN];

  /* Reset fission flag */

  fiss = NO;

  /* Loop over threads and bank (NOTE: type argument is not used in  */
  /* FromTrkBank() at the moment, so this works for photons as well) */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    while ((part = FromTrkBank(id)) > VALID_PTR)
      {
        /* Put particle back to bank */

        ToBank(part, id);

        /* Reset point index */

        idx = 0;

        /* Reset coordinates and time */

        x0 = 0.0;
        y0 = 0.0;
        z0 = 0.0;
        t0 = 0.0;

        x1 = 0.0;
        y1 = 0.0;
        z1 = 0.0;
        t1 = 0.0;

        /* Avoid compiler warning */

        type = -1;
        E = -1.0;
        spd = -1.0;

        /* Get particle type */

        ptype = (long)RDB[part + PARTICLE_TYPE];

        /* Loop over events */

        evn = (long)RDB[part + PARTICLE_PTR_EVENTS];
        while (evn > VALID_PTR)
          {
            /* Get type */

            type = (long)RDB[evn + EVENT_TYPE];

            /* Skip all but collisions, start, leak, bc and cut-offs */

            if ((type != TRACK_END_COLL) && (type != TRACK_END_STRT) &&
                (type != TRACK_END_LEAK) && (type != TRACK_END_TCUT) &&
                (type != TRACK_END_ECUT) && (type != TRACK_END_WCUT) &&
                (type != TRACK_END_SURF) && (type != TRACK_END_BC))
              {
                /* Pointer to next */

                evn = NextItem(evn);

                /* Cycle loop */

                continue;
              }

            /* Update coordinates and time */

            x1 = x0;
            y1 = y0;
            z1 = z0;
            t1 = t0;

            x0 = RDB[evn + EVENT_X];
            CheckValue(FUNCTION_NAME, "x0", "", x0, -INFTY, INFTY);

            y0 = RDB[evn + EVENT_Y];
            CheckValue(FUNCTION_NAME, "y0", "", y0, -INFTY, INFTY);

            z0 = RDB[evn + EVENT_Z];
            CheckValue(FUNCTION_NAME, "z0", "", z0, -INFTY, INFTY);

            t0 = RDB[evn + EVENT_T];
            CheckValue(FUNCTION_NAME, "t0", "", t0, 0.0, INFTY);

            /* Get energy */

            E = RDB[evn + EVENT_E];

            /* Check interval */

            if (t0 <= t)
              break;

            /* Next event */

            evn = NextItem(evn);
          }

        if ((t0 == 0.0) && (t1 == 0.0))
          continue;

        /* Calculate interpolation factor */

        f = (t - t0)/(t1 - t0);
        CheckValue(FUNCTION_NAME, "f", "", f, -INFTY, INFTY);

        /* Check if particle is killed */

        if (f < 0.0)
          {
            /* Override speed */

            if (ptype == PARTICLE_TYPE_NEUTRON)
              spd = RDB[DATA_NEUTRON_SPD];
            else if (ptype == PARTICLE_TYPE_GAMMA)
              spd = RDB[DATA_PHOTON_SPD];
            else
              Die(FUNCTION_NAME, "Invalid particle type");

            /* Check and calculate */

            if (spd < 0.0)
              spd = Speed(ptype, E);

            /* Move to new position */

            x0 = x0 + spd*(t - t0);
          }
        else
          {
            /* Interpolate last point */

            x0 = x0 + f*(x1 - x0);
            y0 = y0 + f*(y1 - y0);
            z0 = z0 + f*(z1 - z0);
          }

        /* Store point in array */

        R[0][idx] = x0;
        R[1][idx] = y0;
        R[2][idx] = z0;

        /* Mark for drawing */

        if (t < t1)
          R[3][idx] = 1;
        else
          R[3][idx] = 0;

        /* Update population count */

        if ((RDB[part + PARTICLE_T0] <= t) && (t <= t1))
          pop++;

        /* Store type */

        R[4][idx] = (double)type;

        /* Update index */

        idx++;

        /* Loop over remaining events */

        while (evn > VALID_PTR)
          {
            /* Get type */

            type = (long)RDB[evn + EVENT_TYPE];

            /* Skip all but collisions, start, leak, bc and cut-offs */

            if ((type != TRACK_END_COLL) && (type != TRACK_END_STRT) &&
                (type != TRACK_END_LEAK) && (type != TRACK_END_TCUT) &&
                (type != TRACK_END_ECUT) && (type != TRACK_END_WCUT) &&
                (type != TRACK_END_SURF) && (type != TRACK_END_BC))
              {
                /* Pointer to next */

                evn = NextItem(evn);

                /* Cycle loop */

                continue;
              }

            /* Get coordinates */

            x0 = RDB[evn + EVENT_X];
            y0 = RDB[evn + EVENT_Y];
            z0 = RDB[evn + EVENT_Z];

            /* Store point in array */

            R[0][idx] = x0;
            R[1][idx] = y0;
            R[2][idx] = z0;
            R[3][idx] = 1;

            /* Store type */

            R[4][idx] = (double)type;

            /* Update index */

            idx++;

            /* Next event */

            evn = NextItem(evn);
          }

        /* Loop over points and draw track */

        for (i = 0; i < idx; i++)
          {
            /* Get coordinates */

            x0 = R[0][i];
            y0 = R[1][i];
            z0 = R[2][i];
          }

        /* Draw tail */

        DrawTail(gpl, im, ncol + (ptype - 1)*TRACK_PLOT_NCOL, R, idx, xp, yp,
                 xmin, xmax, ymin, ymax, zmin, zmax);
      }

  /* Free allocated memory */

  for(n = 0; n < 5; n++)
    Mem(MEM_FREE, R[n]);

  Mem(MEM_FREE, R);

  /* Flush bank to track plot bank */

  FlushBank();

  /* Print time and population */
  /*
  fp = fopen("burst.m", "a");
  fprintf(fp, "pop(%ld,i) = %ld;\n", nftot - nf, pop);
  fclose(fp);
  */
  /* Exit subroutine */

  return fiss;

  /***************************************************************************/
}

/*****************************************************************************/

/***** Coordinates to points *************************************************/

void CoordPt(long gpl, long *n, long *m, long xp, long yp, double x, double y,
             double z, double xmin, double xmax, double ymin, double ymax,
             double zmin, double zmax)
{
  /* Avoid compiler warning */

  *n = -1;
  *m = -1;

  /* Calculate indexes */

  switch ((long)RDB[gpl + GPL_GEOM])
    {
    case PLOT_MODE_YZ:
      {
        /* yz-plot */

        *n = (long)((y - ymin)/(ymax - ymin)*(xp - 1.0));
        *m = (long)((z - zmin)/(zmax - zmin)*(yp - 1.0));
        *m = yp - 1 - *m;

        break;
      }
    case PLOT_MODE_XZ:
      {
        /* xz-plot */

        *n = (long)((x - xmin)/(xmax - xmin)*(xp - 1.0));
        *m = (long)((z - zmin)/(zmax - zmin)*(yp - 1.0));
        *m = yp - 1 - *m;

        break;
      }
    case PLOT_MODE_XY:
      {
        /* xy-plot */

        *n = (long)((x - xmin)/(xmax - xmin)*(xp - 1.0));
        *m = (long)((y - ymin)/(ymax - ymin)*(yp - 1.0));
        *m = yp - 1 - *m;

        break;
      }
    default:
      Die(FUNCTION_NAME, "Invalid plot mode");
    }
}

/*****************************************************************************/

/***** Draw particle tail ****************************************************/

void DrawTail(long gpl, gdImagePtr im, long ncol, double **R, long idx,
              long xp, long yp, double xmin, double xmax, double ymin,
              double ymax, double zmin, double zmax)
{
  long n, m, i, j, prev[2][1000], k;
  double l, d, x0, y0, z0, x1, y1, z1, u, v, w, r, x, y, z, dl, lmax;

  /* Draw without tail for debugging */

  if (1 == 2)
    {
      /* Get point */

      x = R[0][0];
      y = R[1][0];
      z = R[2][0];

      /* Get pixels */

      CoordPt(gpl, &n, &m, xp, yp, x, y, z, xmin, xmax, ymin,
              ymax, zmin, zmax);

      /* Draw (HUOM: t‰ss‰ tulee ongelma OMP:n kanssa kun useamman */
      /* threadin piirt‰m‰t pisteet menee p‰‰llekk‰in)   */

      if (R[3][0] == 1)
        {
          gdImageSetPixel(im, n, m, 0);
          gdImageSetPixel(im, n - 1, m, 0);
          gdImageSetPixel(im, n + 1, m, 0);
          gdImageSetPixel(im, n, m - 1, 0);
          gdImageSetPixel(im, n, m + 1, 0);
        }
      else
        {
          gdImageSetPixel(im, n, m, 2);
          gdImageSetPixel(im, n - 1, m, 2);
          gdImageSetPixel(im, n + 1, m, 2);
          gdImageSetPixel(im, n, m - 1, 2);
          gdImageSetPixel(im, n, m + 1, 2);
        }

      /* Exit */

      return;
    }

  /* Set maximum length and distance between drawn points */

  lmax = RDB[DATA_TRACK_PLOT_HIS_LENGTH];
  dl = 0.001*lmax;

  /* Reset total length */

  l = 0.0;

  /* Reset index to previous pixel positions */

  k = 0;

  /* Loop over points */

  for (i = 1; i < idx; i++)
    {
      /* Get points */

      x0 = R[0][i - 1];
      CheckValue(FUNCTION_NAME, "x0", "", x0, -INFTY, INFTY);

      y0 = R[1][i - 1];
      CheckValue(FUNCTION_NAME, "y0", "", y0, -INFTY, INFTY);

      z0 = R[2][i - 1];
      CheckValue(FUNCTION_NAME, "z0", "", z0, -INFTY, INFTY);

      x1 = R[0][i];
      CheckValue(FUNCTION_NAME, "x1", "", x1, -INFTY, INFTY);

      y1 = R[1][i];
      CheckValue(FUNCTION_NAME, "y1", "", y1, -INFTY, INFTY);

      z1 = R[2][i];
      CheckValue(FUNCTION_NAME, "z1", "", z1, -INFTY, INFTY);

      /* Check periodic */

      if ((R[4][i - 1] == TRACK_END_BC) && (R[4][i] == TRACK_END_SURF))
        continue;

      /* Calculate distance between points */

      r = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) +
               (z0 - z1)*(z0 - z1));

      /* Check */

      if (r > 0.0)
        {
          /* Calculate direction cosines */

          u = (x1 - x0)/r;
          v = (y1 - y0)/r;
          w = (z1 - z0)/r;

          /* Put initial position */

          x = x0;
          y = y0;
          z = z0;

          /* Reset length between points */

          d = 0;

          while (1 != 2)
            {
              /* Move forward */

              x = x + u*dl;
              y = y + v*dl;
              z = z + w*dl;

              /* Check length between points */

              if ((d = d + dl) > r)
                break;

              /* Check total length */

              if ((l = l + dl) > lmax)
                return;

              /* Get pixels */

              CoordPt(gpl, &n, &m, xp, yp, x, y, z, xmin, xmax, ymin,
                      ymax, zmin, zmax);

              /* Compare to previous */

              for (j = 0; j < k; j++)
                if ((n == prev[0][j]) && (m == prev[1][j]))
                  break;

              /* Check */

              if (j == k)
                {
                  /* Store to vector */

                  prev[0][k] = n;
                  prev[1][k] = m;

                  /* Update count */

                  if (k++ == 1000)
                    Die(FUNCTION_NAME, "Need a bigger array");

                  /* Check if drawn (HUOM: t‰ss‰ tulee ongelma OMP:n */
                  /* kanssa kun useamman threadin piirt‰m‰t pisteet  */
                  /* menee p‰‰llekk‰in) */

                  if (R[3][i - 1] == 1)
                    {
                      /* Check mode */

                      if (1 != 2)
                        {
                          /* Calculate color index based on length */

                          j = (long)((double)TRACK_PLOT_NCOL*l/lmax) + ncol;

                          /* Draw pixel */

                          gdImageSetPixel(im, n, m, j);
                        }
                      else
                        {
                          /* Get color index from overlaps (t‰ss‰ on joku */
                          /* bugi, piirt‰‰ vaaleita viivoja) */

                          if ((j = gdImageGetPixel(im, n, m)) < ncol)
                            gdImageSetPixel(im, n, m, ncol);
                          else if (j == ncol + TRACK_PLOT_NCOL - 1)
                            gdImageSetPixel(im, n, m, j);
                          else
                            gdImageSetPixel(im, n, m, j + 1);
                        }
                    }
                }
               }
        }
    }
}

#endif


/*****************************************************************************/
