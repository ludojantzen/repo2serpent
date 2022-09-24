/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : geometryplotter.c                              */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Plots geometry                                               */
/*                                                                           */
/* Comments: - Lisää optio jolla saa valita rajapintojen vs. materiaali-     */
/*             rajojen korostuksen.                                          */
/*                                                                           */
/*           - Tätä kutsutaan monesta eri paikkaa ja eri tapausten käsittely */
/*             on kamala sekamelska. Jaa koko roska jossain vaiheessa perus- */
/*             osaan joka plottaa pelkän geometrian ja erillisiin aliohjel-  */
/*             miin jotka hoitaa muut hommat.                                */
/*                                                                           */
/*           - Siirrettiin iso pätkä paletteja valmistelevaa koodia          */
/*             silmukoiden sisään 9.10.2015. Kamala sekasotku, ihme jos      */
/*             vielä toimii.                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifndef NO_GFX_MODE
#include <gd.h>
#endif

#define FUNCTION_NAME "GeometryPlotter:"

/* Local subroutines */

#ifndef NO_GFX_MODE

void FlashPalette(gdImagePtr im, long *r, long *g, long *b, long type);

#endif

/*****************************************************************************/

void GeometryPlotter(long ini)
{
#ifndef NO_GFX_MODE

  long n, m, gpl, xp, yp, ptr, bor, par, ax, nf, nftot, qp;
  long R[MAX_PLOT_COLORS], G[MAX_PLOT_COLORS], B[MAX_PLOT_COLORS];
  long count, nplot, scale, *mtx;
  double xmin, xmax, ymin, ymax, zmin, zmax, tmp, pos, E, fmin, fmax;
  gdImagePtr im;
  FILE *fp;
  long palette[MAX_PLOT_COLORS];
  char fname[MAX_STR];

  /* Check if plotter is ignored */

  if ((long)RDB[DATA_IGNORE_GEOM_PLOTS] == YES)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check if tracks are printed in file */

  if ((long)RDB[DATA_TRACK_PLOTTER_FILE] == YES)
    {
     /* Check if done */

      if ((long)RDB[DATA_SIMULATION_COMPLETED] == NO)
        return;

      /* Print tracks into file */

      TrackFile();

      /* Return */

      return;
    }

  /* Pointer to geometry plot */

  if ((gpl = (long)RDB[DATA_PTR_GPL0]) < 1)
    {
      /* Check stop mode */

      if ((long)RDB[DATA_STOP_AFTER_PLOT] == STOP_AFTER_PLOT_GEOM)
        exit(-1);

      /* Return */

      return;
    }

  /* Check track plotter mode */

  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
    {
      /* Check if done */

      if ((long)RDB[DATA_SIMULATION_COMPLETED] == NO)
        return;
      else if ((long)RDB[DATA_TRACK_PLOT_ANIM] == YES)
        fprintf(outp, "Plotting particle tracks in geometry:\n\n");
    }
  else if (ini == YES)
    fprintf(outp, "Plotting geometry:\n\n");
  else if (((long)RDB[DATA_SOURCE_PT_ANIM] == NO) &&
           ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_DELDYN) &&
           ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_DYN) &&
           ((long)RDB[DATA_RUN_VR_ITER] == NO))
    return;

  /* Set plotter flag on */

  WDB[DATA_PLOTTER_MODE] = (double)YES;

  /***************************************************************************/

  /***** Loop over track plotter frames **************************************/

  /* Nuber of frames */

  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
    nftot = (long)RDB[DATA_TRACK_PLOT_FRAMES];
  else
    nftot = 1;

  /* Loop over frames */

  for (nf = 0; nf < nftot; nf++)
    {
      /***********************************************************************/

      /***** Print progress of frames ****************************************/

      /* Print progress */

      if (nftot > 200)
        {
          if (!(nf % 20))
            fprintf(outp, " %3.0f%% complete\n",
                    100.0*((double)nf)/((double)nftot));
        }
      else if (nftot > 100)
        {
          if (!(nf % 10))
            fprintf(outp, " %3.0f%% complete\n",
                    100.0*((double)nf)/((double)nftot));
        }
      else if (nftot > 50)
        {
          if (!(nf % 5))
            fprintf(outp, " %3.0f%% complete\n",
                    100.0*((double)nf)/((double)nftot));
        }
      else if (nftot > 1)
        fprintf(outp, " %3.0f%% complete\n",
                100.0*((double)nf)/((double)nftot));

      /***********************************************************************/

      /***** Loop over plot definitions **************************************/

      /* Get pointer to geometry plots */

      gpl = (long)RDB[DATA_PTR_GPL0];
      CheckPointer(FUNCTION_NAME, "(gpl)", DATA_ARRAY, gpl);

      /* Reset counter and get number of plots */

      count = 0;
      nplot = ListSize(gpl);

      /* Loop over plots */

      gpl = FirstItem(gpl);
      while (gpl > VALID_PTR)
        {
          /* Skip importance mesh plots if no weight windows */

          if ((long)RDB[gpl + GPL_IMP_SCALE] > 0)
            if ((long)RDB[DATA_PTR_WWD0] < VALID_PTR)
              {
                /* Pointer to next */

                gpl = NextItem(gpl);

                /* Cycle loop */

                continue;
              }

          /* Put RGB values for plot */

          if ((long)RDB[gpl + GPL_IMP_SCALE] > 0)
            PutPlotColors(YES, R, G, B);
          else
            PutPlotColors(NO, R, G, B);

          /* Print progress */

          if ((nftot == 1) && ((long)RDB[DATA_PART_PTR_SOURCE] < VALID_PTR))
            fprintf(outp, " %3.0f%% complete\n",
                    100.0*(count++)/((double)nplot));

          /*******************************************************************/

          /***** Set file name ***********************************************/

          /* Check VR iteration */

          if ((long)RDB[DATA_RUN_VR_ITER] == YES)
            {
              sprintf(fname, "%s_geom%ld_vr%ld.png",
                      GetText(DATA_PTR_INPUT_FNAME),
                      (long)RDB[gpl + GPL_IDX],
                      (long)RDB[DATA_VR_ITER_IDX]);

            }
          else if ((long)RDB[DATA_BURN_STEP] > 0)
            {
              if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
                {
                  /* Check number of frames */

                  if (nftot < 2)
                    sprintf(fname, "%s_trck%ld_bu%ld.png",
                            GetText(DATA_PTR_INPUT_FNAME),
                            (long)RDB[gpl + GPL_IDX],
                            (long)RDB[DATA_BURN_STEP]);
                  else
                    sprintf(fname, "%s_trck%ld_bu%ld_frame%s.png",
                            GetText(DATA_PTR_INPUT_FNAME),
                            (long)RDB[gpl + GPL_IDX], (long)RDB[DATA_BURN_STEP],
                            IdxStr(nftot - nf, nftot));
                }
              else
                sprintf(fname, "%s_geom%ld_bu%ld.png",
                        GetText(DATA_PTR_INPUT_FNAME),
                        (long)RDB[gpl + GPL_IDX], (long)RDB[DATA_BURN_STEP]);
            }
          else
            {
              if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
                {
                  /* Check number of frames */

                  if (nftot < 2)
                    sprintf(fname, "%s_trck%ld.png",
                            GetText(DATA_PTR_INPUT_FNAME),
                            (long)RDB[gpl + GPL_IDX]);
                  else
                    sprintf(fname, "%s_trck%ld_frame%s.png",
                            GetText(DATA_PTR_INPUT_FNAME),
                            (long)RDB[gpl + GPL_IDX],
                            IdxStr(nftot - nf, nftot));
                }
              else if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
                {
                  sprintf(fname, "%s_spt%ld_frame%s.png",
                          GetText(DATA_PTR_INPUT_FNAME),
                          (long)RDB[gpl + GPL_IDX],
                          IdxStr((long)RDB[DATA_CYCLE_IDX],
                        (long)RDB[DATA_CRIT_CYCLES] +
                                 (long)RDB[DATA_CRIT_SKIP]));
                }
              else if (((long)RDB[DATA_SIMULATION_MODE] ==
                        SIMULATION_MODE_DELDYN) ||
                       ((long)RDB[DATA_SIMULATION_MODE] ==
                        SIMULATION_MODE_DYN))
                {
                  sprintf(fname, "%s_geom%ld_t%ld.png",
                          GetText(DATA_PTR_INPUT_FNAME),
                          (long)RDB[gpl + GPL_IDX], (long)RDB[DATA_DYN_TB]);

                }
              else
                sprintf(fname, "%s_geom%ld.png", GetText(DATA_PTR_INPUT_FNAME),
                        (long)RDB[gpl + GPL_IDX]);
            }

          /********************************************************************/

          /***** Initialize plot **********************************************/

          /* Get boundaries of geometry */

          xmin = RDB[DATA_GEOM_MINX];
          xmax = RDB[DATA_GEOM_MAXX];
          ymin = RDB[DATA_GEOM_MINY];
          ymax = RDB[DATA_GEOM_MAXY];
          zmin = RDB[DATA_GEOM_MINZ];
          zmax = RDB[DATA_GEOM_MAXZ];

          /* Check if boundaries are set by user */

          if (RDB[gpl + GPL_XMIN] != -INFTY)
            xmin = RDB[gpl + GPL_XMIN];

          if (RDB[gpl + GPL_XMAX] !=  INFTY)
            xmax = RDB[gpl + GPL_XMAX];

          if (RDB[gpl + GPL_YMIN] != -INFTY)
            ymin = RDB[gpl + GPL_YMIN];

          if (RDB[gpl + GPL_YMAX] !=  INFTY)
            ymax = RDB[gpl + GPL_YMAX];

          if (RDB[gpl + GPL_ZMIN] != -INFTY)
            zmin = RDB[gpl + GPL_ZMIN];

          if (RDB[gpl + GPL_ZMAX] !=  INFTY)
            zmax = RDB[gpl + GPL_ZMAX];

          /* Check boundaries and swap */

          if (xmin > xmax)
            {
              tmp = xmax;
              xmax = xmin;
              xmin = tmp;
            }

          if (ymin > ymax)
            {
              tmp = ymax;
              ymax = ymin;
              ymin = tmp;
            }

          if (zmin > zmax)
            {
              tmp = zmax;
              zmax = zmin;
              zmin = tmp;
            }

          /* Move limits away from boundaries to avoid numerical errors in */
          /* some systems. */

          xmin = xmin + EXTRAP_L;
          xmax = xmax - EXTRAP_L;
          ymin = ymin + EXTRAP_L;
          ymax = ymax - EXTRAP_L;
          zmin = zmin + EXTRAP_L;
          zmax = zmax - EXTRAP_L;

          /* Set image size */

          xp = (long)RDB[gpl + GPL_PIX_X];
          yp = (long)RDB[gpl + GPL_PIX_Y];

          /* Get position */

          pos = RDB[gpl + GPL_POS];

          /* Boundary mode (1 = plot cell boundaries, 2 = plot material */
          /*                boundaries, 3 = plot both) */

          if ((((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0) ||
               ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)) &&
              ((long)RDB[gpl + GPL_IMP_SCALE] == 0))
            bor = 0;
          else
            bor = (long)RDB[gpl + GPL_PLOT_BOUND];

          /* Allocate memory */

          mtx = (long *)Mem(MEM_ALLOC, xp*yp, sizeof(long));

          /* Get weight window particle type */

          if ((ptr = (long)RDB[DATA_PTR_WWD0]) > VALID_PTR)
            {
              /* Check for multiple */

              if (NextItem(ptr) < VALID_PTR)
                par = (long)RDB[ptr + WWD_PARTICLE_TYPE];
              else
                par = 0;
            }
          else
            par = 0;

          /* Axis type */

          ax = (long)RDB[gpl + GPL_GEOM];
          CheckValue(FUNCTION_NAME, "ax", "", ax, 1, 3);

          /* Scale for importance plotter */

          scale = (long)RDB[gpl + GPL_IMP_SCALE];
          CheckValue(FUNCTION_NAME, "scale", "", scale, 0, 4);

          /* Get energy */

          E = RDB[gpl + GPL_IMP_E];
          CheckValue(FUNCTION_NAME, "E", "", E, -INFTY, INFTY);

          /* Get maximum and maximum importances */

          fmin = RDB[gpl + GPL_IMP_MIN];
          CheckValue(FUNCTION_NAME, "fmin", "", fmin, 0.0, 1E+200);

          fmax = RDB[gpl + GPL_IMP_MAX];
          CheckValue(FUNCTION_NAME, "fmax", "", fmax, fmin, 1E+200);

          /* Get quick plotter option */

          qp = (long)RDB[DATA_QUICK_PLOT_MODE];

          /* Plot image */

          n = PlotImage(gpl, mtx, xp, yp, xmin, xmax, ymin, ymax, zmin,
                        zmax, pos, bor, par, ax, scale, E, fmin, fmax, qp);

          /* Check for geometry errors */

          if ((n < 0) && (ini == YES))
            {
              if (n == GEOM_ERROR_NO_CELL)
                fprintf(outp, "Geometry errors in plot %s (no cell).\n",
                        fname);
              else if (n == GEOM_ERROR_MULTIPLE_CELLS)
                fprintf(outp, "Geometry errors in plot %s (overlap).\n",
                        fname);
              else
                fprintf(outp, "Geometry errors in plot %s (unknown).\n",
                        fname);
            }

          /*******************************************************************/

          /***** Draw image **************************************************/

          /* Check matrix */

          for (n = 0; n < xp; n++)
            for (m = 0; m < yp; m++)
              if ((mtx[n*yp + m] < 0) || (mtx[n*yp + m] > MAX_PLOT_COLORS - 1))
                Die(FUNCTION_NAME, "Invalid colour index %ld (%ld %ld)",
                    mtx[n*yp + m], n, m);

          /* Create image */

          im = gdImageCreate(xp, yp);

          /* Generate palette */

          for(n = 0; n < MAX_PLOT_COLORS; n++)
            palette[n] = gdImageColorAllocate(im, R[n], G[n], B[n]);

          /* Start parallel timer */

          StartTimer(TIMER_OMP_PARA);

          /* Draw image */

#ifdef OPEN_MP
#pragma omp parallel private(n, m)
#endif
          {

#ifdef OPEN_MP
#pragma omp for
#endif
            for (n = 0; n < xp; n++)
              for (m = 0; m < yp; m++)
                {
                  /* Flip y-axis in xy-plot */

                  if ((long)RDB[gpl +  GPL_GEOM] == PLOT_MODE_XY)
                    gdImageSetPixel(im, n, m, palette[mtx[n*yp + yp - m - 1]]);
                  else
                    gdImageSetPixel(im, n, m, palette[mtx[n*yp + m]]);
                }
          }

          /* Stop parallel timer */

          StopTimer(TIMER_OMP_PARA);

          /* Plot tracks */

          if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
            {
              /* Plot tracks (return value indicates if fission has */
              /* occurred during frame) */

              n = PlotTracks(gpl, im, nf, xmin, xmax, ymin, ymax,
                                zmin, zmax, xp, yp);

              /* Check flag and flash palette */

              if (n == YES)
                FlashPalette(im, R, G, B, 2);
            }

          /*******************************************************************/

          /***** Write in file and cleanup ***********************************/

          /* Open file for writing */

          if ((fp = fopen(fname, "w")) == NULL)
            Die(FUNCTION_NAME, "Unable to open file for writing");

          /* Write image (png format) */

          gdImagePng(im, fp);

          /* Free image */

          gdImageDestroy(im);

          /* Free matrix */

          Mem(MEM_FREE, mtx);

          /* Close file */

          fclose(fp);

          /*******************************************************************/

          /* Next plot */

          gpl = NextItem(gpl);
        }
    }

  if (nf > 1)
    fprintf(outp, " %3.0f%% complete\n\nOK.\n\n", 100.0);

  /***************************************************************************/

  /* Reset plotter mode */

  WDB[DATA_PLOTTER_MODE] = (double)NO;

  /* Done */

  if ((long)RDB[DATA_PART_PTR_SOURCE] < VALID_PTR)
    fprintf(outp, " 100%% complete\n\n");

  /* Check stop mode */

  if ((long)RDB[DATA_STOP_AFTER_PLOT] == STOP_AFTER_PLOT_GEOM)
    exit(-1);

  /****************************************************************************/

#endif
}

/*****************************************************************************/

/***** Flash palette to indicate fission in track plotter mode ***************/

#ifndef NO_GFX_MODE

void FlashPalette(gdImagePtr im, long *r0, long *g0, long *b0, long type)
{
  long n, r, g, b;

  /* Reset palette */

  for(n = 0; n < MAX_PLOT_COLORS; n++)
    gdImageColorDeallocate(im, n);

  /* Check type */

  if (type == 1)
    {
      /***********************************************************************/

      /***** Invert colors ***************************************************/

      for(n = 0; n < MAX_PLOT_COLORS; n++)
        {
          r = 255 - r0[n];
          g = 255 - g0[n];
          b = 255 - b0[n];

          gdImageColorAllocate(im, r, g, b);
        }

      /***********************************************************************/
    }
  else if (type == 2)
    {
      /***********************************************************************/

      /***** Brighten colors *************************************************/

      /* First 5 are unchanged */

      for(n = 0; n < 5; n++)
        gdImageColorAllocate(im, r0[n], g0[n], b0[n]);

      /* Adjust remaining */

      for(n = 5; n < MAX_PLOT_COLORS; n++)
        {
          if ((r = r0[n] + 10) > 255)
            r = 255;
          if ((g = g0[n] + 10) > 255)
            g = 255;
          if ((b = b0[n] + 10) > 255)
            b = 255;

          gdImageColorAllocate(im, r, g, b);
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid palette flash type");
}

#endif

/*****************************************************************************/
