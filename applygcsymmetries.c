/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : applygcsymmetries.c                            */
/*                                                                           */
/* Created:       2014/05/15 (JLe)                                           */
/* Last modified: 2020/07/21 (ARi)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Applies symmetry options in ADF's and pin-power form         */
/*              functions.                                                   */
/*                                                                           */
/* Comments: - Tää laskee korjaukset jotka lisätään bufferiin että saadaan   */
/*             symmetrisiä arvoja.                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ApplyGCSymmetries:"

/* Surfaces */

#define W 0
#define S 1
#define E 2
#define N 3

/* Corners */

#define NW 0
#define NE 1
#define SE 2
#define SW 3

/*****************************************************************************/

void ApplyGCSymmetries(long gcu)
{
  long ng, adf, ptr, surf, type, ns, nc, nmax, n, m, i, sym, mode;
  double *f, *f0, val, max;

  /* Check pointer to group constant universe */

  CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);

  /***************************************************************************/

  /***** Fluxes and currents in ADF's ****************************************/

  /* Reset maximum difference */

  max = 0.0;

  /* Get pointer to ADF structure and symmetry */

  if ((adf = (long)RDB[gcu + GCU_PTR_ADF]) > VALID_PTR)
    if ((sym = (long)RDB[adf + ADF_SYM]) > 0)
      {
        /* Get pointer to outer boundary surface */

        surf = (long)RDB[adf + ADF_PTR_SURF];
        CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

        /* Get surface type */

        type = (long)RDB[surf + SURFACE_TYPE];

        /* Get number of surfaces and corners */

        ns = (long)RDB[adf + ADF_NSURF];
        CheckValue(FUNCTION_NAME, "ns", "", ns, 1, 6);

        nc = (long)RDB[adf + ADF_NCORN];
        CheckValue(FUNCTION_NAME, "nc", "", nc, 0, 6);

        /* Get number of energy groups */

        ng = (long)RDB[DATA_ERG_FG_NG];
        CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

        /* Allocate memory for temporary arrays  */

        if (ns > nc)
          {
            f0 = (double *)Mem(MEM_ALLOC, ns, sizeof(double));
            f = (double *)Mem(MEM_ALLOC, ns, sizeof(double));
          }
        else
          {
            f0 = (double *)Mem(MEM_ALLOC, nc, sizeof(double));
            f = (double *)Mem(MEM_ALLOC, nc, sizeof(double));
          }

        /* Loop over parameters */

        for (i = 0; i < 15; i++)
          {
            /* Avoid compiler warning */

            ptr = -1;
            nmax = 0;
            mode = 0;

            /* Get pointer and set number of values */

            if (i == 0)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_HET_SURF_FLUX];
                nmax = ns;
                mode = 1;
              }
            else if (i == 1)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_HET_CORN_FLUX];
                nmax = nc;
                mode = 2;
              }
            else if (i == 2)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_SURF_IN_CURR];
                nmax = ns;
                mode = 1;
              }
            else if (i == 3)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_SURF_OUT_CURR];
                nmax = ns;
                mode = 1;
              }
            else if (i == 4)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_SURF_NET_CURR];
                nmax = ns;
                mode = 1;
              }
            else if (i == 5)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_MID_IN_CURR];
                nmax = ns;
                mode = 1;
              }
            else if (i == 6)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_MID_OUT_CURR];
                nmax = ns;
                mode = 1;
              }
            else if (i == 7)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_MID_NET_CURR];
                nmax = ns;
                mode = 1;
              }
            else if (i == 8)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_CORN_IN_CURR];
                nmax = nc;
                mode = 2;
              }
            else if (i == 9)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_CORN_OUT_CURR];
                nmax = nc;
                mode = 2;
              }
            else if (i == 10)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_CORN_NET_CURR];
                nmax = nc;
                mode = 2;
              }
            else if (i == 11)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_SGN_HET_SURF_FLUX];
                nmax = ns;
                mode = 1;
              }
            else if (i == 12)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_SGN_SURF_IN_CURR];
                nmax = ns;
                mode = 1;
              }
            else if (i == 13)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_SGN_SURF_OUT_CURR];
                nmax = ns;
                mode = 1;
              }
            else if (i == 14)
              {
                ptr = (long)RDB[gcu + GCU_RES_FG_DF_SGN_SURF_NET_CURR];
                nmax = ns;
                mode = 1;
              }
            else
              {
                Die(FUNCTION_NAME, "Overflow");
              }

            /* Check pointer */

            if (ptr < VALID_PTR)
              continue;

            /* Loop over energy groups */

            for (n = 0; n < ng; n++)
              {
                /* Read data */

                for (m = 0; m < nmax; m++)
                  f0[m] = BufVal(ptr, m, n);

                /* Check geometry type */

                switch (type)
                  {
                  case SURF_SQC:
                  case SURF_RECT:
                  case SURF_CUBOID:
                  case SURF_CUBE:
                    {
                      /* Infinite or truncated square prism (2D or 3D case) */
                      /* averaging only across radial directions            */

                      if (sym == 1)
                        {
                          /* Full symmetry */

                          if ((type == SURF_SQC) ||
                              (type == SURF_CUBOID) ||
                              (type == SURF_CUBE))
                            {
                              /* Average over all values if sqc */

                              val = (f0[0] + f0[1] + f0[2] + f0[3])/4.0;
                              f[0] = val - f0[0];
                              f[1] = val - f0[1];
                              f[2] = val - f0[2];
                              f[3] = val - f0[3];
                            }
                          else
                            {
                              /* Average over NS and EW and all corners */
                              /* if rect */

                              if (mode == 1)
                                {
                                  val = (f0[S] + f0[N])/2.0;
                                  f[S] = val - f0[S];
                                  f[N] = val - f0[N];
                                  val = (f0[W] + f0[E])/2.0;
                                  f[W] = val - f0[W];
                                  f[E] = val - f0[E];
                                }
                              else
                                {
                                  val = (f0[0] + f0[1] + f0[2] + f0[3])/4.0;
                                  f[0] = val - f0[0];
                                  f[1] = val - f0[1];
                                  f[2] = val - f0[2];
                                  f[3] = val - f0[3];
                                }
                            }
                        }
                      else if (sym == 2)
                        {
                          /* S1 / NS symmetry */

                          if (mode == 1)
                            {
                              val = (f0[S] + f0[N])/2.0;
                              f[W] = 0.0;
                              f[S] = val - f0[S];
                              f[E] = 0.0;
                              f[N] = val - f0[N];
                            }
                          else
                            {
                              val = (f0[NW] + f0[SW])/2.0;
                              f[NW] = val - f0[NW];
                              f[SW] = val - f0[SW];
                              val = (f0[NE] + f0[SE])/2.0;
                              f[NE] = val - f0[NE];
                              f[SE] = val - f0[SE];
                            }
                        }
                      else if (sym == 3)
                        {
                          /* C1 / SENW symmetry */

                          if (mode == 1)
                            {
                              val = (f0[W] + f0[S])/2.0;
                              f[W] = val - f0[W];
                              f[S] = val - f0[S];
                              val = (f0[E] + f0[N])/2.0;
                              f[E] = val - f0[E];
                              f[N] = val - f0[N];
                            }
                          else
                            {
                              val = (f0[SE] + f0[NW])/2.0;
                              f[SW] = 0.0;
                              f[SE] = val - f0[SE];
                              f[NW] = val - f0[NW];
                              f[NE] = 0.0;
                            }
                        }
                      else if (sym == 4)
                        {
                          /* S2 / EW symmetry */

                          if (mode == 1)
                            {
                              val = (f0[W] + f0[E])/2.0;
                              f[W] = val - f0[W];
                              f[S] = 0.0;
                              f[E] = val - f0[E];
                              f[N] = 0.0;
                            }
                          else
                            {
                              val = (f0[SW] + f0[SE])/2.0;
                              f[SW] = val - f0[SW];
                              f[SE] = val - f0[SE];
                              val = (f0[NW] + f0[NE])/2.0;
                              f[NW] = val - f0[NW];
                              f[NE] = val - f0[NE];
                            }
                        }
                      else if (sym == 5)
                        {
                          /* C2 / SWNE symmetry */

                          if (mode == 1)
                            {
                              val = (f0[W] + f0[N])/2.0;
                              f[W] = val - f0[W];
                              f[N] = val - f0[N];
                              val = (f0[E] + f0[S])/2.0;
                              f[E] = val - f0[E];
                              f[S] = val - f0[S];
                            }
                          else
                            {
                              val = (f0[SW] + f0[NE])/2.0;
                              f[SW] = val - f0[SW];
                              f[SE] = 0.0;
                              f[NW] = 0.0;
                              f[NE] = val - f0[NE];
                            }
                        }
                      else if (sym == 6)
                        {
                          /* NSEW symmetry */

                          if (mode == 1)
                            {
                              val = (f0[S] + f0[N])/2.0;
                              f[S] = val - f0[S];
                              f[N] = val - f0[N];

                              val = (f0[W] + f0[E])/2.0;
                              f[W] = val - f0[W];
                              f[E] = val - f0[E];
                            }
                          else
                            {
                              val = (f0[0] + f0[1] + f0[2] + f0[3])/4.0;
                              f[0] = val - f0[0];
                              f[1] = val - f0[1];
                              f[2] = val - f0[2];
                              f[3] = val - f0[3];
                            }
                        }
                      else if (sym != 0)
                        Die(FUNCTION_NAME, "Invalid symmetry");

                      /* Break case */

                      break;
                    }
                  case SURF_HEXYC:
                  case SURF_HEXXC:
                    {
                      /* Infinite hexagonal prism (2D case) */

                      if (sym == 1)
                        {
                          /* Full symmetry */

                          val = (f0[0] + f0[1] + f0[2] +
                                 f0[3] + f0[4] + f0[5])/6.0;
                          f[0] = val - f0[0];
                          f[1] = val - f0[1];
                          f[2] = val - f0[2];
                          f[3] = val - f0[3];
                          f[4] = val - f0[4];
                          f[5] = val - f0[5];
                        }
                      else if (sym == 2)
                        {
                          /* S1-symmetry */

                          if (mode == 1)
                            {
                              f[0] = 0.0;
                              f[3] = 0.0;
                              val = (f0[1] + f0[5])/2.0;
                              f[1] = val - f0[1];
                              f[5] = val - f0[5];
                              val = (f0[2] + f0[4])/2.0;
                              f[2] = val - f0[2];
                              f[4] = val - f0[4];
                            }
                          else
                            {
                              f[0] = 0.0;
                              val = (f0[1] + f0[2])/2.0;
                              f[1] = val - f0[1];
                              f[2] = val - f0[2];
                              f[3] = 0.0;
                              val = (f0[4] + f0[5])/2.0;
                              f[4] = val - f0[4];
                              f[5] = val - f0[5];
                            }
                        }
                      else if (sym == 3)
                        {
                          /* C1-symmetry */

                          if (mode == 1)
                            {
                              val = (f0[0] + f0[1])/2.0;
                              f[0] = val - f0[0];
                              f[1] = val - f0[1];
                              val = (f0[2] + f0[5])/2.0;
                              f[2] = val - f0[2];
                              f[5] = val - f0[5];
                              val = (f0[3] + f0[4])/2.0;
                              f[3] = val - f0[3];
                              f[4] = val - f0[4];
                            }
                          else
                            {
                              val = (f0[0] + f0[2])/2.0;
                              f[0] = val - f0[0];
                              f[2] = val - f0[2];
                              f[1] = 0.0;
                              val = (f0[3] + f0[5])/2.0;
                              f[3] = val - f0[3];
                              f[5] = val - f0[5];

                              f[4] = 0.0;
                            }
                        }
                      else if (sym == 4)
                        {
                          /* S2-symmetry */

                          if (mode == 1)
                            {
                              f[1] = 0.0;
                              f[4] = 0.0;
                              val = (f0[0] + f0[2])/2.0;
                              f[0] = val - f0[0];
                              f[2] = val - f0[2];
                              val = (f0[3] + f0[5])/2.0;
                              f[3] = val - f0[3];
                              f[5] = val - f0[5];
                            }
                          else
                            {
                              val = (f0[0] + f0[1])/2.0;
                              f[0] = val - f0[0];
                              f[1] = val - f0[1];
                              val = (f0[2] + f0[5])/2.0;
                              f[2] = val - f0[2];
                              f[5] = val - f0[5];
                              val = (f0[3] + f0[4])/2.0;
                              f[3] = val - f0[3];
                              f[4] = val - f0[4];
                            }
                        }
                      else if (sym == 5)
                        {
                          /* C2-symmetry */

                          if (mode == 1)
                            {
                              val = (f0[0] + f0[3])/2.0;
                              f[0] = val - f0[0];
                              f[3] = val - f0[3];
                              val = (f0[1] + f0[2])/2.0;
                              f[1] = val - f0[1];
                              f[2] = val - f0[2];
                              val = (f0[4] + f0[5])/2.0;
                              f[4] = val - f0[4];
                              f[5] = val - f0[5];
                            }
                          else
                            {
                              f[0] = 0.0;
                              val = (f0[1] + f0[5])/2.0;
                              f[1] = val - f0[1];
                              f[5] = val - f0[5];
                              val = (f0[2] + f0[4])/2.0;
                              f[2] = val - f0[2];
                              f[4] = val - f0[4];

                              f[3] = 0.0;
                            }
                        }
                      else if (sym == 6)
                        {
                          /* S3-symmetry */

                          if (mode == 1)
                            {
                              f[2] = 0.0;
                              f[5] = 0.0;
                              val = (f0[1] + f0[3])/2.0;
                              f[1] = val - f0[1];
                              f[3] = val - f0[3];
                              val = (f0[0] + f0[4])/2.0;
                              f[0] = val - f0[0];
                              f[4] = val - f0[4];
                            }
                          else
                            {
                              val = (f0[0] + f0[5])/2.0;
                              f[0] = val - f0[0];
                              f[5] = val - f0[5];
                              val = (f0[1] + f0[4])/2.0;
                              f[1] = val - f0[1];
                              f[4] = val - f0[4];
                              val = (f0[2] + f0[3])/2.0;
                              f[2] = val - f0[2];
                              f[3] = val - f0[3];
                            }
                        }
                      else if (sym == 7)
                        {
                          /* C3-symmetry */

                          if (mode == 1)
                            {
                              val = (f0[0] + f0[5])/2.0;
                              f[0] = val - f0[0];
                              f[5] = val - f0[5];
                              val = (f0[1] + f0[4])/2.0;
                              f[1] = val - f0[1];
                              f[4] = val - f0[4];
                              val = (f0[2] + f0[3])/2.0;
                              f[2] = val - f0[2];
                              f[3] = val - f0[3];
                            }
                          else
                            {
                              val = (f0[0] + f0[4])/2.0;
                              f[0] = val - f0[0];
                              f[4] = val - f0[4];
                              f[2] = 0.0;
                              val = (f0[1] + f0[3])/2.0;
                              f[1] = val - f0[1];
                              f[3] = val - f0[3];
                              f[5] = 0.0;
                            }
                        }
                      else if (sym != 0)
                        {
                          Die(FUNCTION_NAME, "Invalid symmetry");
                        }

                      /* Break case */

                      break;
                    }
                  }

                /* Mark scoring buffer unreduced */

                WDB[DATA_BUF_REDUCED] = (double)NO;

                /* Get maximum difference for error checking */
                /* (surface flux and currents) */

                if ((i == 0) || (i == 2) || (i == 3))
                  for (m = 0; m < nmax; m++)
                    if (fabs(f[m]/f0[m]) > max)
                      max = fabs(f[m]/f0[m]);

                /* Put data */

                for (m = 0; m < nmax; m++)
                  AddBuf(f[m], 1.0, ptr, 0, -1, m, n);

                /* Reduce buffer */

                ReduceBuffer();
              }
          }

        /* Free temporary arrays */

        Mem(MEM_FREE, f0);
        Mem(MEM_FREE, f);
      }

  /* Print warning if statistics is good enough */

  if ((long)(RDB[DATA_MICRO_CALC_BATCH_SIZE]*RDB[DATA_CYCLE_BATCH_SIZE]) >
      1000*20 && max > 0.1)
    {
      Note(0, "ADF symmetry option may be wrong");
    }

  /***************************************************************************/
}

/*****************************************************************************/
