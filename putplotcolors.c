/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : putplotcolors.c                                */
/*                                                                           */
/* Created:       2020/06/17 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Sets the RGB values for color scheme used in geometry        */
/*              plotter.                                                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PutPlotColors:"

/*****************************************************************************/

void PutPlotColors(long imp, long *R, long *G, long *B)
{
#ifndef NO_GFX_MODE

  long id, ptr, ncol, nifc, mat, n, i, m, nc0;
  long unsigned seed;
  double f;

  /* Loop over OpenMP threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Init random number sequence */

      seed = ReInitRNG(id);
      SEED[id*RNG_SZ] = seed;
    }

  /***************************************************************************/

  /***** IFC and source point animation **************************************/

  /* Reduce number of colors if IFC or source point animation */
  /* is in use */

  if ((ptr = (long)RDB[DATA_PTR_IFC0]) > VALID_PTR)
    {
      /* Use only half */

      ncol = (long)(((double)MAX_PLOT_COLORS)/2.0);

      /* Calculate number of ifc materials*/

      nifc = 0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check whether material is linked to IFC */
          /* IFC materials currently use parent material color */
          /* in order to better show the distributions */

          if ((RDB[mat + MATERIAL_USE_IFC] == (double)YES)
              && ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR))
            nifc++;

          mat = NextItem(mat);
        }

      /* Calculate number of IFC colors */

      nifc = (long)(((double)ncol - 10)/(double)nifc);
    }
  else if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
    {
      ncol = (long)(((double)MAX_PLOT_COLORS)/2.0);
      nifc = 0;
    }
  else if (((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0) &&
           ((long)RDB[DATA_TRACK_PLOT_ANIM] == YES))
    {
      ncol = MAX_PLOT_COLORS - 2*TRACK_PLOT_NCOL;
      nifc = 0;
    }
  else
    {
      ncol = MAX_PLOT_COLORS;
      nifc = 0;
    }

  /***************************************************************************/

  /***** Random colors *******************************************************/

  /* Randomize */

  for(n = 5; n < MAX_PLOT_COLORS; n++)
    {
      /* More contrast when domain decompositoin is on */

      if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
        {
          R[n] = (long)(255.0*RandF(0));
          G[n] = (long)(255.0*RandF(0));
          B[n] = (long)(255.0*RandF(0));
        }
      else
        {
          R[n] = (long)(255.0*(RandF(0) + 1.0)/2);
          G[n] = (long)(255.0*(RandF(0) + 1.0)/2);
          B[n] = (long)(255.0*(RandF(0) + 1.0)/2);
        }
    }

  /* Void colour */

  R[0] = 0;
  G[0] = 0;
  B[0] = 0;

  /* No cell error */

  R[1] = 0;
  G[1] = 255;
  B[1] = 0;

  /* Multiple cell error */

  R[2] = 255;
  G[2] = 0;
  B[2] = 0;

  /* Pointer error */

  R[3] = 255;
  G[3] = 255;
  B[3] = 0;

  /* Undefined density factor */

  R[4] = 255;
  G[4] = 0;
  B[4] = 255;

  /***************************************************************************/

  /***** Material-wise RGB values ********************************************/

  /* User-defined material colours */

  n = 5;
  i = NO;

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if material uses interface */

      if ((long)RDB[mat + MATERIAL_USE_IFC] == (double)NO)
        {
          /* Check if colour is defined */

          if ((m = (long)RDB[mat + MATERIAL_RGB]) > 0)
            {
              if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT])
                  > VALID_PTR)
                WDB[mat + MATERIAL_COLOUR_IDX] =
                  RDB[ptr + MATERIAL_COLOUR_IDX];
              else if (n < ncol)
                {
                  /* Check index */

                  if ((n < 0) || (n > MAX_PLOT_COLORS - 1))
                    Die(FUNCTION_NAME, "Index error");

                  R[n] = (long)(m/1000000.0);
                  G[n] = (long)((m - 1000000.0*R[n])/1000.0);
                  B[n] = (long)(m - 1000000.0*R[n] - 1000.0*G[n]);

                  /* Set index */

                  WDB[mat + MATERIAL_COLOUR_IDX] = (double)n;

                  /* Next colour */

                  n++;
                }
              else
                {
                  /* Use random index */

                  WDB[mat + MATERIAL_COLOUR_IDX] =
                    (double)((long)(RandF(0)*(ncol - 5))) + 5.0;

                  /* Flag for warning */

                  i = YES;
                }
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Check if all are used */

  if (i == YES)
    Note(0,
         "Number of pre-assigned material colors exceeds max (252)");

  /***************************************************************************/

  /***** Domain-specific colors in DD mode ***********************************/

  /* Check for domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    {
      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if divided */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
            WDB[mat + MATERIAL_COLOUR_IDX] =
              (double)n + RDB[mat + MATERIAL_MPI_ID];

          /* Next material */

          mat = NextItem(mat);
        }

      /* Update index */

      if (n + mpitasks < ncol)
        n = n + mpitasks;
    }

  /***************************************************************************/

  /***** Put remaining colour indexes ****************************************/

  m = n;

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if material uses interface */

      if ((long)RDB[mat + MATERIAL_USE_IFC] == (double)NO)
        {
          /* Cycle if all are used */

          if (m > ncol - 1)
            m = n;

          /* One more check */

          if (m > ncol - 1)
            m = ncol - 1;

          /* Check if index is not already given */

          if ((long)RDB[mat + MATERIAL_COLOUR_IDX] == 0)
            WDB[mat + MATERIAL_COLOUR_IDX] = (double)m++;
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Set total number of colors */

  nc0 = n;

  /***************************************************************************/

  /***** Set IFC colors ******************************************************/

  /* User-defined material colours */

  n = ncol;

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if material uses interface */

      if ((long)RDB[mat + MATERIAL_USE_IFC] == (double)YES)
        {
          /* Number of colors */

          nc0 = ncol;

          /* Check if colour is defined */

          if ((m = (long)RDB[mat + MATERIAL_RGB]) > 0)
            {
              if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) >
                  VALID_PTR)
                WDB[mat + MATERIAL_COLOUR_IDX] =
                  RDB[ptr + MATERIAL_COLOUR_IDX];
              else if (n < MAX_PLOT_COLORS)
                {
                  /* Check index */

                  if ((n < 0) || (n > MAX_PLOT_COLORS - 1))
                    Die(FUNCTION_NAME, "Index error");

                  R[n] = (long)(m/1000000.0);
                  G[n] = (long)((m - 1000000.0*R[n])/1000.0);
                  B[n] = (long)(m - 1000000.0*R[n] - 1000.0*G[n]);

                  /* Set index */

                  WDB[mat + MATERIAL_COLOUR_IDX] = (double)n;

                  /* Next colour */

                  n = n + nifc;
                }
              else
                WDB[mat + MATERIAL_COLOUR_IDX] = 0;
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Put remaining colour indexes */

  m = n;

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if material uses interface */

      if ((long)RDB[mat + MATERIAL_USE_IFC] == (double)YES)
        {
          /* Cycle if all are used */

          if (m > MAX_PLOT_COLORS - 1)
            m = n;

          /* One more check */

          if (m > MAX_PLOT_COLORS - 1)
            m = MAX_PLOT_COLORS - 1;

          /* Check if index is not already given */

          if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) >
              VALID_PTR)
            WDB[mat + MATERIAL_COLOUR_IDX] =
              RDB[ptr + MATERIAL_COLOUR_IDX];
          else if ((long)RDB[mat + MATERIAL_COLOUR_IDX] == 0)
            {
              WDB[mat + MATERIAL_COLOUR_IDX] = (double)m;
              m = m + nifc;
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Create shades */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if material uses interface */

      if ((long)RDB[mat + MATERIAL_USE_IFC] == (double)YES)
        {
          /* Get index to color */

          if ((n = (long)RDB[mat + MATERIAL_COLOUR_IDX]) < ncol)
            Die(FUNCTION_NAME, "Indexing error");

          /* Loop over shades */

          for (m = 1; m < nifc; m++)
            {
              /* Check index */

              if ((n < 0) || (n > MAX_PLOT_COLORS - 1) ||
                  (n + nifc - m < 0) || (n + nifc - m > MAX_PLOT_COLORS - 1))
                Die(FUNCTION_NAME, "Index error");

              /* Put colours */

              R[n + nifc - m] =
                R[n]*((double)(m - 1))/((double)(nifc - 1));
              G[n + nifc - m] =
                G[n]*((double)(m - 1))/((double)(nifc - 1));
              B[n + nifc - m] =
                B[n]*((double)(m - 1))/((double)(nifc - 1));
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Store values */

  WDB[DATA_PLOTTER_NIFC] = (double)nifc;
  WDB[DATA_PLOTTER_NCOL] = (double)ncol;
  WDB[DATA_PLOTTER_NC0] = (double)nc0;

  /***************************************************************************/

  /***** Adjust colors for track plotting ************************************/

  /* Check if tracks are recorded for plotting */

  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
    {
      /* Brighten up... */

      for(n = 4; n < MAX_PLOT_COLORS; n++)
        {
          /* Check index */

          if ((n < 0) || (n > MAX_PLOT_COLORS - 1))
            Die(FUNCTION_NAME, "Index error");

          if ((R[n] = (long)(R[n] + 50)) > 255)
            R[n] = 255;
          if ((G[n] = (long)(G[n] + 50)) > 255)
            G[n] = 255;
          if ((B[n] = (long)(B[n] + 50)) > 255)
            B[n] = 255;
        }

      /* Particle */

      R[0] = 0;
      G[0] = 0;
      B[0] = 0;

      /* Tail */

      if ((long)RDB[DATA_TRACK_PLOT_ANIM] == YES)
        {
          /* Photons */

          for (n = 0; n < TRACK_PLOT_NCOL; n++)
            {
              f = ((double)n/((double)TRACK_PLOT_NCOL));

              /* Green */

              R[ncol + n] = (long)(2.0*f*255.0) - 200;
              G[ncol + n] = (long)(2.0*f*255.0);
              B[ncol + n] = (long)(2.0*f*255.0) - 200;

              /* Check limits */

              if (R[ncol + n] > 255)
                R[ncol + n] = 255;
              else if (R[ncol + n] < 0)
                R[ncol + n] = 0;

              if (G[ncol + n] > 255)
                G[ncol + n] = 255;
              else if (G[ncol + n] < 0)
                G[ncol + n] = 0;

              if (B[ncol + n] > 255)
                B[ncol + n] = 255;
              else if (B[ncol + n] < 0)
                B[ncol + n] = 0;
            }

          /* Neutrons */

          for (n = 0; n < TRACK_PLOT_NCOL; n++)
            {
              f = ((double)n/((double)TRACK_PLOT_NCOL));

              /* Purple */

              R[ncol + TRACK_PLOT_NCOL + n] = (long)(2.0*f*255.0);
              G[ncol + TRACK_PLOT_NCOL + n] = (long)(2.0*f*255.0) - 200;
              B[ncol + TRACK_PLOT_NCOL + n] = (long)(2.0*f*255.0);

              /* Check limits */

              if (R[ncol + TRACK_PLOT_NCOL + n] > 255)
                R[ncol + TRACK_PLOT_NCOL + n] = 255;
              else if (R[ncol + TRACK_PLOT_NCOL + n] < 0)
                R[ncol + TRACK_PLOT_NCOL + n] = 0;

              if (G[ncol + TRACK_PLOT_NCOL + n] > 255)
                G[ncol + TRACK_PLOT_NCOL + n] = 255;
              else if (G[ncol + TRACK_PLOT_NCOL + n] < 0)
                G[ncol + TRACK_PLOT_NCOL + n] = 0;

              if (B[ncol + TRACK_PLOT_NCOL + n] > 255)
                B[ncol + TRACK_PLOT_NCOL + n] = 255;
              else if (B[ncol + TRACK_PLOT_NCOL + n] < 0)
                B[ncol + TRACK_PLOT_NCOL + n] = 0;
            }
        }
    }

  /***************************************************************************/

  /***** Set colors for importance mesh **************************************/

  /* Check mode */

  if (imp == YES)
    {
      /* Override palette */

      MakePalette(&R[nc0], &G[nc0], &B[nc0], ncol - nc0,
                  PALETTE_BLUE_RED);
    }

  /***************************************************************************/

  /***** Put colors for source point animation *******************************/

  /* Check index */

  if ((ncol < 1) || (n > MAX_PLOT_COLORS))
    Die(FUNCTION_NAME, "Index error");

  if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
    MakePalette(&R[ncol], &G[ncol], &B[ncol], ncol,
                (long)RDB[DATA_SOURCE_PT_ANIM_PALETTE]);

#endif

  /***************************************************************************/
}

/*****************************************************************************/
