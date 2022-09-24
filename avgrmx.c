/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : avgrmx.c                                       */
/*                                                                           */
/* Created:       2018/11/22 (JLe)                                           */
/* Last modified: 2019/04/10 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates average importances for empty mesh cells          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AvgRMX:"

/*****************************************************************************/

void AvgRMX(long wwd)
{
  long rmx, msh0, loc0, loc1, ptr, ng, n, i, j, count;
  double sum, div, *I0, *I2, f0, f1, diff, V;
  const double *I1;

  /* Get smoothing factors */

  if ((f0 = RDB[wwd + WWD_SMOOTH_F0]) < 0.0)
    return;

  f1 = RDB[wwd + WWD_SMOOTH_F1];
  CheckValue(FUNCTION_NAME, "f1", "", f1, ZERO, INFTY);

  /* Logs */

  f0 = log10(f0);
  f1 = log10(f1);

  /* Print */

  fprintf(outp, "Applying smoothing for weight-window mesh...\n");

  /* Two loops: 1st to eliminate bad values, 2nd to smooth results */

  for (i = 0; i < 2; i++)
    {
      /* Loop over data */

      rmx = (long)RDB[DATA_PTR_RMX0];
      while (rmx > VALID_PTR)
        {
          /* Get number of energy groups */

          ng = (long)RDB[rmx + RMX_NG];
          CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

          /* Loop over iterations */

          do
            {
              /* Reset count */

              count = 0;

              /* Loop over mesh */

              msh0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              while (msh0 > VALID_PTR)
                {
                  /* Pointer to importances */

                  ptr = (long)RDB[msh0 + RMX_CELL_IMP_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  I0 = &WDB[ptr];

                  /* Use solution array for temporarily storing averages */

                  ptr = (long)RDB[msh0 + RMX_CELL_WRK_IN_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  I2 = &WDB[ptr];

                  /* Loop over energy groups */

                  for (n = 0; n < ng; n++)
                    {
                      /* Reset sum */

                      sum = 0.0;
                      div = 0.0;

                      /* Reset temporary */

                      I2[n] = 0.0;

                      /* Loop over neighbours */

                      j = 0;

                      loc0 = (long)RDB[msh0 + RMX_CELL_PTR_BOUNDS];
                      while (loc0 > VALID_PTR)
                        {
                          /* Pointer to mesh cell */

                          loc1 = (long)RDB[loc0 + RMX_CELL_BOUND_PTR_CELL];
                          CheckPointer(FUNCTION_NAME, "loc1", DATA_ARRAY, loc1);

                          /* Pointer to importances */

                          ptr = (long)RDB[loc1 + RMX_CELL_IMP_CURR];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                          I1 = &RDB[ptr];

                          /* Pointer to segment-wise solutions */

                          ptr = (long)RDB[msh0 + RMX_CELL_ADJ_SOL_IN_CURR];
                          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                          /* Volume or importance-weighting */

                          if (1 != 2)
                            V = RDB[msh0 + RMX_CELL_RVOL];
                          else
                            V = RDB[ptr + j*ng + n];

                          /* Check value */

                          CheckValue(FUNCTION_NAME, "V", "", V, 0.0, INFTY);

                          /* Update count */

                          j++;

                          /* Add to sum */

                          if (I1[n] > 0.0)
                            {
                              sum = sum + V*I1[n];
                              div = div + V;
                            }

                          /* Next neighbour */

                          loc0 = NextItem(loc0);
                        }

                      /* Check importance */

                      if ((div == 0.0) || (I0[n] == 0.0))
                        continue;
                      else
                        diff = fabs(log10(fabs(sum/div/I0[n])));

                      /* Check loop type */

                      if (i == 0)
                        {
                          /* Eliminate bad results */

                          if (diff > f0)
                            {
                              I2[n] = sum/div;
                              count++;
                            }
                        }
                      else
                        {
                          /* Smoothing */

                          if (diff > f1)
                            {
                              I2[n] = (I0[n] + sum)/(div + 1.0);
                              count++;
                            }
                        }
                    }

                  /* Next */

                  msh0 = NextItem(msh0);
                }

              /* Put smoothed data */

              msh0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              while (msh0 > VALID_PTR)
                {
                  /* Pointer to importances */

                  ptr = (long)RDB[msh0 + RMX_CELL_IMP_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  I0 = &WDB[ptr];

                  /* Use solution array for temporarily storing averages */

                  ptr = (long)RDB[msh0 + RMX_CELL_WRK_IN_CURR];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  I2 = &WDB[ptr];

                  /* Loop over energy groups */

                  for (n = 0; n < ng; n++)
                    if (I2[n] != 0.0)
                      {
                        /* Replace */

                        I0[n] = I2[n];

                        /* Reset */

                        I2[n] = 0.0;
                      }

                  /* Next */

                  msh0 = NextItem(msh0);
                }
            }
          while (count > 0);

          /* Next */

          rmx = NextItem(rmx);
        }
    }

  fprintf(outp, "OK.\n\n");

  /* Loop over data and normalize */

  rmx = (long)RDB[DATA_PTR_RMX0];
  while (rmx > VALID_PTR)
    {
      /* Re-normalize importances */

      NormalizeImp(rmx);

      /* Next */

      rmx = NextItem(rmx);
    }
}

/*****************************************************************************/
