/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : trackfile.c                                    */
/*                                                                           */
/* Created:       2020/01/31 (JLe)                                           */
/* Last modified: 2020/01/31 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Prints tracks in file                                        */
/*                                                                           */
/* Comments: - Toi paino pit‰‰ viel‰ katsoa ett‰ miss‰ kohtaa se             */
/*             muuttuu (pit‰‰ varmaan tehd‰ sama juttu kuin energialle).     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TrackFile:"

#define DATA_ITEMS 12

/*****************************************************************************/

void TrackFile()
{
  long id, evn, part, type, idx, sz, n;
  double x, y, z, u, v, w, E, wgt, t, r, f, dose;
  double *data[DATA_ITEMS];
  FILE *fp;

  /* T‰‰ on niist‰ Olli Hyvˆsen M&C 2017 -laskuista, ilmeisesti SCALE:ssa */
  /* MT9504 : ANSI standard (1977) gamma flux-to-dose-rate factors        */
  /* (rem/h)/(photons/cm2/s) */

  static double flux_to_dose[48][2] = {
    {1.00E-02,  0.00E+00},
    {2.00E-02,  1.88E-06},
    {3.00E-02,  7.88E-07},
    {4.50E-02,  3.95E-07},
    {6.00E-02,  2.82E-07},
    {7.00E-02,  2.60E-07},
    {7.50E-02,  2.58E-07},
    {1.00E-01,  2.69E-07},
    {1.50E-01,  3.24E-07},
    {2.00E-01,  4.36E-07},
    {2.60E-01,  5.76E-07},
    {3.00E-01,  7.07E-07},
    {4.00E-01,  8.70E-07},
    {4.50E-01,  1.03E-06},
    {5.10E-01,  1.14E-06},
    {5.12E-01,  1.20E-06},
    {6.00E-01,  1.28E-06},
    {7.00E-01,  1.44E-06},
    {8.00E-01,  1.60E-06},
    {9.00E-01,  1.76E-06},
    {1.00E+00,  1.91E-06},
    {1.20E+00,  2.11E-06},
    {1.33E+00,  2.34E-06},
    {1.44E+00,  2.49E-06},
    {1.50E+00,  2.60E-06},
    {1.57E+00,  2.68E-06},
    {1.66E+00,  2.77E-06},
    {1.80E+00,  2.91E-06},
    {2.00E+00,  3.10E-06},
    {2.15E+00,  3.29E-06},
    {2.35E+00,  3.47E-06},
    {2.50E+00,  3.65E-06},
    {2.75E+00,  3.84E-06},
    {3.00E+00,  4.08E-06},
    {3.50E+00,  4.41E-06},
    {4.00E+00,  4.83E-06},
    {4.50E+00,  5.22E-06},
    {5.00E+00,  5.43E-06},
    {5.50E+00,  6.00E-06},
    {6.00E+00,  6.37E-06},
    {6.50E+00,  6.74E-06},
    {7.00E+00,  7.11E-06},
    {7.50E+00,  7.48E-06},
    {8.00E+00,  7.84E-06},
    {1.00E+01,  8.74E-06},
    {1.20E+01,  1.02E-05},
    {1.40E+01,  1.17E-05},
    {2.00E+01,  1.43E-05}};

  /* Print */

  fprintf(outp, "Writing tracks in file...\n");

  /* Output */

  fp = fopen("tracks.dat", "w");

  /* Reset count */

  idx = 0;

  /* Loop over threads and bank (NOTE: type argument is not used in  */
  /* FromTrkBank() at the moment, so this works for photons as well) */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    while ((part = FromTrkBank(id)) > VALID_PTR)
      {
        /* Update count */

        idx++;

        /* Put particle back to bank */

        ToBank(part, id);

        /* Get size of event buffer */

        evn = (long)RDB[part + PARTICLE_PTR_EVENTS];
        CheckPointer(FUNCTION_NAME, "(evn)", DATA_ARRAY, evn);
        sz = ListSize(evn);

        if (sz == 0)
          continue;

        /* Allocate memory for data */

        for (n = 0; n < DATA_ITEMS; n++)
          data[n] = Mem(MEM_ALLOC, sz, sizeof(double));

        /* Reset energy, position and direction of motion */

        E = 0.0;
        x = 0.0;
        y = 0.0;
        z = 0.0;
        u = 0.0;
        v = 0.0;
        w = 0.0;

        /* Reset number of events */

        sz = 0;

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
                /*
                (type != TRACK_END_SURF) &&
                */
                (type != TRACK_END_WWIN) &&
                (type != TRACK_END_BC))
              {
                /* Pointer to next */

                evn = NextItem(evn);

                /* Cycle loop */

                continue;
              }

            /* Multiplications are recorded as additional starting */
            /* points. Skip to avoid zero direction vector because */
            /* the same position is recorded twice. */

            if ((type == TRACK_END_STRT) && (NextItem(evn) > VALID_PTR))
              {
                /* Next event */

                evn = NextItem(evn);

                /* Cycle loop */

                continue;
              }

            /* Calculate direction (check if last) */

            if (E > 0.0)
              {
                /* Calculate direction */

                u = x - RDB[evn + EVENT_X];
                v = y - RDB[evn + EVENT_Y];
                w = z - RDB[evn + EVENT_Z];

                /* Normalize */

                r = sqrt(u*u + v*v + w*w);
                CheckValue(FUNCTION_NAME, "r", "", r, ZERO, INFTY);

                u = u/r;
                v = v/r;
                w = w/r;
              }

            /* Get coordinates */

            x = RDB[evn + EVENT_X];
            CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);

            y = RDB[evn + EVENT_Y];
            CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);

            z = RDB[evn + EVENT_Z];
            CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);

            /* Get weight and time */

            wgt = RDB[evn + EVENT_WGT];
            CheckValue(FUNCTION_NAME, "wgt", "", wgt, -INFTY, INFTY);

            t = RDB[evn + EVENT_T];
            CheckValue(FUNCTION_NAME, "t", "", t, -INFTY, INFTY);

            /* Store data */

            data[0][sz] = t;
            data[1][sz] = x;
            data[2][sz] = y;
            data[3][sz] = z;
            data[4][sz] = u;
            data[5][sz] = v;
            data[6][sz] = w;
            data[7][sz] = E;
            data[8][sz] = wgt;

            /* Find index */

            if (E < flux_to_dose[0][0])
              dose = flux_to_dose[0][1];
            else if (E > flux_to_dose[47][0])
              dose = flux_to_dose[47][1];
            else
              {
                /* Search interval */

                for (n = 0; n < 48; n++)
                  if (flux_to_dose[n][0] > E)
                    break;

                /* Check */

                CheckValue(FUNCTION_NAME, "n", "", n, 1, 47);
                CheckValue(FUNCTION_NAME, "E", "", E, flux_to_dose[n - 1][0],
                           flux_to_dose[n][0]);

                /* Interpolate */

                f = (E - flux_to_dose[n - 1][0])/
                  (flux_to_dose[n][0] - flux_to_dose[n - 1][0]);
                CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);

                /* Put dose contribution */

                data[9][sz] = f*(flux_to_dose[n][1] - flux_to_dose[n - 1][1])
                  + flux_to_dose[n - 1][1];
              }

            /* Update energy */

            E = RDB[evn + EVENT_E];
            CheckValue(FUNCTION_NAME, "E", "", E, 0.0, INFTY);

            /* Update index */

            sz++;

            /* Next event */

            evn = NextItem(evn);
          }

        /* Print size */

        fprintf(fp, "%ld\n", sz);

        /* Print events in reverse order */

        while (sz-- > 0)
          {
            fprintf(fp, "%11.5E %12.5E %12.5E %12.5E %8.5f %8.5f %8.5f %12.5E %12.5E %12.5E\n",
                    data[0][sz], data[1][sz], data[2][sz],
                    data[3][sz], data[4][sz], data[5][sz], data[6][sz],
                    data[7][sz], data[8][sz], data[9][sz]);
          }

        /* Free allocated memory */

        for (n = 0; n < DATA_ITEMS; n++)
          Mem(MEM_FREE, data[n]);
      }

  /* Close file */

  fclose(fp);

  /* Flush bank to track plot bank */

  FlushBank();

  /* Exit OK */

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
