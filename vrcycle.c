/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : vrcycle.c                                      */
/*                                                                           */
/* Created:       2018/07/03 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Runs iterations for variance reduction.                      */
/*                                                                           */
/* Comments: - Toi CON-juttu on hirveä viritelmä joka ei varmaan toimi.      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "VRCycle:"

/*****************************************************************************/

void VRCycle()
{
  long wwd, loc0, rmx, n, nbtch, pop, nmax, i, sz0, sz1, ptr;

  /* Pointer to weight window definition */

  wwd = (long)RDB[DATA_PTR_WWD0];
  CheckPointer(FUNCTION_NAME, "wwd", DATA_ARRAY, wwd);

  /* Check multiple definitions */

  if (NextItem(wwd) > VALID_PTR)
    Error(wwd, "Multipple WWD definitions in iteration");

  /* Memory count */

  WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] + (double)MemCount();

  /* Remember batch and population size */

  nbtch = (long)RDB[DATA_SRC_BATCHES];
  pop = (long)RDB[DATA_SRC_POP];

  /* Switch weight windows off if importances are defined */

  if ((long)RDB[DATA_USE_GEOM_IMP] == YES)
    WDB[DATA_USE_WEIGHT_WINDOWS] = (double)NO;

  /* Loop over iterations */

  loc0 = (long)RDB[wwd + WWD_PTR_ITER];
  while (loc0 > VALID_PTR)
    {
      /* Check iteration type */

      if ((long)RDB[loc0 + WWD_ITER_TYPE] == WWIN_ITER_USR)
        {
          /*******************************************************************/

          /***** User-defined ************************************************/

          /* Pointer to reponse matrix solution */

          rmx = (long)RDB[loc0 + WWD_ITER_PTR_RMX];
          CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

          /* Link RMX solution */

          WDB[DATA_PTR_RMX0] = (double)rmx;

          /* Put density factor */

          WDB[DATA_GLOBAL_DF] = RDB[loc0 + WWD_ITER_DF];

          /* Check if final cycle */

          if (NextItem(loc0) < VALID_PTR)
            WDB[DATA_RMTX_CALC] = (double)NO;

          /* Prepare transport cycle */

          PrepareTransportCycle();

          /* Run transport cycle */

          TransportCycle();

          /* Put type */

          WDB[wwd + WWD_TYPE] = (double)WWD_MESH_TYPE_ITER;

          /* Put scaling and normalization factors */

          WDB[wwd + WWD_MULT] = -1.0;
          WDB[wwd + WWD_POW] = -1.0;
          WDB[wwd + WWD_NORM_FACT] = 1.0;

          /* Link mesh pointers */

          WDB[wwd + WWD_PTR_MESH] = RDB[rmx + RMX_PTR_MESH];
          WDB[wwd + WWD_PTR_RMX] = (double)rmx;

          /* Write output */

          WriteWWMesh(rmx);

          /* Update iteration index */

          WDB[DATA_VR_ITER_IDX] = RDB[DATA_VR_ITER_IDX] + 1.0;

          /*******************************************************************/
        }
      else if ((long)RDB[loc0 + WWD_ITER_TYPE] == WWIN_ITER_GEO)
        {
          /*******************************************************************/

          /***** Automatically adapted (geometry based) **********************/

          /* Pointer to reponse matrix solution */

          rmx = (long)RDB[loc0 + WWD_ITER_PTR_RMX];
          CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

          /* Link RMX solution */

          WDB[DATA_PTR_RMX0] = (double)rmx;

          /* Number of mfp iterations */

          nmax = (long)RDB[loc0 + WWD_ITER_SPLIT_LOOPS];

          /* Loop over runs */

          for (n = 0; n < nmax + 1; n++)
            {
              /* Skip track sampling if not first cycle */

              if ((PrevItem(loc0) > VALID_PTR) && (n < nmax))
                continue;

              /* Print */

              if (n == 0)
                fprintf(outp, "Sampling tracks to find dense materials:\n\n");

              /* Link iteration pointer */

              WDB[rmx + RMX_PTR_ITER] = (double)loc0;

              /* Run transport cycle */

              if (n < nmax)
                ProbeVRMesh(rmx);
              else
                {
                  /* Prepare transport cycle */

                  PrepareTransportCycle();

                  /* Run simulation */

                  TransportCycle();

                  /* Apply smoothing */

                  AvgRMX(wwd);
                }

              /* Put type */

              WDB[wwd + WWD_TYPE] = (double)WWD_MESH_TYPE_ITER;

              /* Put scaling and normalization factors */

              WDB[wwd + WWD_MULT] = -1.0;
              WDB[wwd + WWD_POW] = -1.0;
              WDB[wwd + WWD_NORM_FACT] = 1.0;

              /* Link mesh pointers */

              WDB[wwd + WWD_PTR_MESH] = RDB[rmx + RMX_PTR_MESH];
              WDB[wwd + WWD_PTR_RMX] = (double)rmx;

              /* Get size before split */

              ptr = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              sz0 = ListSize(ptr);

              /* Adjust mesh (density or importance) */

              if (n < nmax)
                sz1 = AdjustRMX(rmx, 0);
              else
                sz1 = AdjustRMX(rmx, 2);

              /* Neighbours */

              while ((i = AdjustRMX(rmx, 3)) > 0)
                sz1 = sz1 + i;

              /* Memory count */

              WDB[DATA_TOT_RMX_BYTES] = RDB[DATA_TOT_RMX_BYTES]
                + (double)MemCount();

              /* Check splits and stop */

              if ((n < nmax - 1) && (sz1 == 0))
                n = nmax - 1;

              /* Print */

              if (n < nmax)
                fprintf(outp, " - %ld/%ld mesh cells split by density.\n",
                        sz1, sz0);
              else
                fprintf(outp, " - %ld/%ld mesh cells split by importance.\n",
                        sz1, sz0);

              /* Add newline */

              if (n > nmax - 2)
                fprintf(outp, "\n");
            }

          /* Write output */

          WriteWWMesh(rmx);

          /* Update iteration index */

          WDB[DATA_VR_ITER_IDX] = RDB[DATA_VR_ITER_IDX] + 1.0;

          /*******************************************************************/
        }
      else if ((long)RDB[loc0 + WWD_ITER_TYPE] == WWIN_ITER_TRK)
        {
          /*******************************************************************/

          /***** Automatically adapted (sampled tracks) **********************/

          /* Pointer to reponse matrix solution */

          rmx = (long)RDB[loc0 + WWD_ITER_PTR_RMX];
          CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

          /* Link RMX solution */

          WDB[DATA_PTR_RMX0] = (double)rmx;

          /* Number of mfp iterations */

          nmax = (long)RDB[loc0 + WWD_ITER_SPLIT_LOOPS];

          /* Loop over runs */

          for (n = 0; n < nmax + 1; n++)
            {
	      /* Put mode */

	      if (n < nmax)
		{
		  WDB[DATA_RMTX_MFP_CALC] = (double)YES;
		  WDB[DATA_SRC_BATCHES] = 1.0;
                  WDB[DATA_SRC_POP] = RDB[loc0 + WWD_ITER_SPLIT_TRACKS];
		}
	      else
		{
		  WDB[DATA_RMTX_MFP_CALC] = (double)NO;
		  WDB[DATA_SRC_BATCHES] = (double)nbtch;
                  WDB[DATA_SRC_POP] = (double)pop;
		}

              /* Prepare transport cycle or only clear statistics */

              if ((n == 0) || (n == nmax))
                PrepareTransportCycle();
              else
                ClearStat(-1);

              /* Print */

              if (n == 0)
                fprintf(outp,
                        "Running histories to find thick materials:\n\n");

              /* Link iteration pointer */

              /*
              WDB[wwd + WWD_PTR_MESH] = RDB[rmx + RMX_PTR_MESH];
	      WDB[wwd + WWD_PTR_RMX] = (double)rmx;
              */
              WDB[rmx + RMX_PTR_ITER] = (double)loc0;

              /* Run transport cycle */

              TransportCycle();

              /* Apply smoothing to RMX solution */

              if (n > nmax - 1)
                AvgRMX(wwd);

              /* Put type */

              WDB[wwd + WWD_TYPE] = (double)WWD_MESH_TYPE_ITER;

              /* Put scaling and normalization factors */

              WDB[wwd + WWD_MULT] = -1.0;
              WDB[wwd + WWD_POW] = -1.0;
              WDB[wwd + WWD_NORM_FACT] = 1.0;

              /* Link mesh pointers */

              WDB[wwd + WWD_PTR_MESH] = RDB[rmx + RMX_PTR_MESH];
              WDB[wwd + WWD_PTR_RMX] = (double)rmx;

              /* Get size before split */

              ptr = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              sz0 = ListSize(ptr);

              /* Adjust mesh (density or importance) */

              if (n < nmax)
                sz1 = AdjustRMX(rmx, 0);
              else
                sz1 = AdjustRMX(rmx, 2);

              /* Neighbours */

              while ((i = AdjustRMX(rmx, 3)) > 0)
                sz1 = sz1 + i;

              /* Memory count */

              WDB[DATA_TOT_RMX_BYTES] = RDB[DATA_TOT_RMX_BYTES]
                + (double)MemCount();

              /* Check splits and stop */

              if ((n < nmax - 1) && (sz1 == 0))
                n = nmax - 1;

              /* Print */

              if (n < nmax)
                fprintf(outp, " - %ld/%ld mesh cells split by density.\n",
                        sz1, sz0);
              else
                fprintf(outp, " - %ld/%ld mesh cells split by importance.\n",
                        sz1, sz0);

              /* Add newline */

              if (n > nmax - 2)
                fprintf(outp, "\n");
            }

          /* Write output */

          WriteWWMesh(rmx);

          /* Update iteration index */

          WDB[DATA_VR_ITER_IDX] = RDB[DATA_VR_ITER_IDX] + 1.0;

          /*******************************************************************/
        }
      else if ((long)RDB[loc0 + WWD_ITER_TYPE] == WWIN_ITER_CON)
        {
          /*******************************************************************/

          /***** Automatically adapted (contribution based) ******************/

          /* Pointer to reponse matrix solution */

          rmx = (long)RDB[loc0 + WWD_ITER_PTR_RMX];
          CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

          /* Link RMX solution */

          WDB[DATA_PTR_RMX0] = (double)rmx;

          /* Number of mfp iterations */

          nmax = 1;

          /* Loop over runs */

          for (n = 0; n < nmax + 1; n++)
            {
              /* Skip track sampling if not first cycle */

              if ((PrevItem(loc0) > VALID_PTR) && (n < nmax))
                continue;

              /* Print */

              if (n == 0)
                fprintf(outp, "Sampling tracks to find dense materials:\n\n");

              /* Link iteration pointer */

              WDB[rmx + RMX_PTR_ITER] = (double)loc0;

              /* Run transport cycle */

              if (n < nmax)
                ProbeVRMesh(rmx);
              else
                {
                  /* Prepare transport cycle */

                  PrepareTransportCycle();

                  /* Run simulation */

                  TransportCycle();

                  /* Apply smoothing */

                  AvgRMX(wwd);
                }

              /* Put type */

              WDB[wwd + WWD_TYPE] = (double)WWD_MESH_TYPE_ITER;

              /* Put scaling and normalization factors */

              WDB[wwd + WWD_MULT] = -1.0;
              WDB[wwd + WWD_POW] = -1.0;
              WDB[wwd + WWD_NORM_FACT] = 1.0;

              /* Link mesh pointers */

              WDB[wwd + WWD_PTR_MESH] = RDB[rmx + RMX_PTR_MESH];
              WDB[wwd + WWD_PTR_RMX] = (double)rmx;

              /* Get size before split */

              ptr = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              sz0 = ListSize(ptr);

              /* Adjust mesh (density or importance) */

              if (n < nmax)
                sz1 = AdjustRMX(rmx, 0);
              else
                sz1 = AdjustRMX(rmx, 4);

              /* Neighbours */

              while ((i = AdjustRMX(rmx, 3)) > 0)
                sz1 = sz1 + i;

              /* Memory count */

              WDB[DATA_TOT_RMX_BYTES] = RDB[DATA_TOT_RMX_BYTES]
                + (double)MemCount();

              /* Check splits and stop */

              if ((n < nmax - 1) && (sz1 == 0))
                n = nmax - 1;

              /* Print */

              if (n < nmax)
                fprintf(outp, " - %ld/%ld mesh cells split by density.\n",
                        sz1, sz0);
              else
                fprintf(outp, " - %ld/%ld mesh cells split by importance.\n",
                        sz1, sz0);

              /* Add newline */

              if (n > nmax - 2)
                fprintf(outp, "\n");
            }

          /* Write output */

          WriteWWMesh(rmx);

          /* Update iteration index */

          WDB[DATA_VR_ITER_IDX] = RDB[DATA_VR_ITER_IDX] + 1.0;

          /*******************************************************************/
        }
      else
        Die(FUNCTION_NAME, "Invalid iteration mode");

      /* Switch weight windows back on */

      WDB[DATA_USE_WEIGHT_WINDOWS] = (double)YES;

      /* Make sure geometry importances are off */

      WDB[DATA_USE_GEOM_IMP] = (double)NO;

      /* Plot geometry */

      GeometryPlotter(NO);

      /* Next iteration */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/
}

/*****************************************************************************/
