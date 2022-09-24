/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : coefoutput.c                                   */
/*                                                                           */
/* Created:       2014/04/15 (JLe)                                           */
/* Last modified: 2019/12/04 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Writes group constant output for coefficient calculations    */
/*                                                                           */
/* Comments: - Tähän muutettiin toi PrintCoeVals() juuri ennen 2.1.25 jakoa. */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CoefOutput:"

/*****************************************************************************/

void CoefOutput()
{
  long loc0, loc1, loc2, ptr, n, np, gcu, nv, det, i, j, nr, ZAI, MT, I, ntot;
  char tmpstr[MAX_STR], **params;
  FILE *fp;

  /* Date and time */

  time_t t = time(NULL);
  struct tm tm = *localtime(&t);

  /***************************************************************************/

  /***** Default list of output parameters ***********************************/

  char *default_params[] = {
    "INF_FLX",
    "INF_KINF",
    "INF_REP_TIME",
    "INF_PROMPT_LIFE",
    "INF_TOT",
    "INF_CAPT",
    "INF_FISS",
    "INF_NSF",
    "INF_KAPPA",
    "INF_INVV",
    "INF_NUBAR",
    "INF_ABS",
    "INF_REMXS",
    "INF_RABSXS",
    "INF_CHIT",
    "INF_CHIP",
    "INF_CHID",
    "INF_I135_YIELD",
    "INF_XE135_YIELD",
    "INF_XE135M_YIELD",
    "INF_PM149_YIELD",
    "INF_SM149_YIELD",
    "INF_I135_MICRO_ABS",
    "INF_XE135_MICRO_ABS",
    "INF_XE135M_MICRO_ABS",
    "INF_PM149_MICRO_ABS",
    "INF_SM149_MICRO_ABS",
    "INF_I135_MACRO_ABS",
    "INF_XE135_MACRO_ABS",
    "INF_XE135M_MACRO_ABS",
    "INF_PM149_MACRO_ABS",
    "INF_SM149_MACRO_ABS",
    "INF_S0",
    "INF_S1",
    "INF_S2",
    "INF_S3",
    "INF_S4",
    "INF_S5",
    "INF_S6",
    "INF_S7",
    "INF_SP0",
    "INF_SP1",
    "INF_SP2",
    "INF_SP3",
    "INF_SP4",
    "INF_SP5",
    "INF_SP6",
    "INF_SP7",
    "INF_SCATT0",
    "INF_SCATT1",
    "INF_SCATT2",
    "INF_SCATT3",
    "INF_SCATT4",
    "INF_SCATT5",
    "INF_SCATT6",
    "INF_SCATT7",
    "INF_SCATTP0",
    "INF_SCATTP1",
    "INF_SCATTP2",
    "INF_SCATTP3",
    "INF_SCATTP4",
    "INF_SCATTP5",
    "INF_SCATTP6",
    "INF_SCATTP7",
    "INF_TRANSPXS",
    "INF_DIFFCOEF",
    "B1_KINF",
    "B1_KEFF",
    "B1_REP_TIME",
    "B1_PROMPT_LIFE",
    "B1_B2",
    "B1_ERR",
    "B1_FLX",
    "B1_FISS_FLX",
    "B1_TOT",
    "B1_CAPT",
    "B1_FISS",
    "B1_NSF",
    "B1_KAPPA",
    "B1_INVV",
    "B1_NUBAR",
    "B1_ABS",
    "B1_REMXS",
    "B1_RABSXS",
    "B1_CHIT",
    "B1_CHIP",
    "B1_CHID",
    "B1_I135_YIELD",
    "B1_XE135_YIELD",
    "B1_PM149_YIELD",
    "B1_SM149_YIELD",
    "B1_I135_MICRO_ABS",
    "B1_XE135_MICRO_ABS",
    "B1_PM149_MICRO_ABS",
    "B1_SM149_MICRO_ABS",
    "B1_I135_MACRO_ABS",
    "B1_XE135_MACRO_ABS",
    "B1_PM149_MACRO_ABS",
    "B1_SM149_MACRO_ABS",
    "B1_S0",
    "B1_S1",
    "B1_S2",
    "B1_S3",
    "B1_S4",
    "B1_S5",
    "B1_S6",
    "B1_S7",
    "B1_SP0",
    "B1_SP1",
    "B1_SP2",
    "B1_SP3",
    "B1_SP4",
    "B1_SP5",
    "B1_SP6",
    "B1_SP7",
    "B1_SCATT0",
    "B1_SCATT1",
    "B1_SCATT2",
    "B1_SCATT3",
    "B1_SCATT4",
    "B1_SCATT5",
    "B1_SCATT6",
    "B1_SCATT7",
    "B1_SCATTP0",
    "B1_SCATTP1",
    "B1_SCATTP2",
    "B1_SCATTP3",
    "B1_SCATTP4",
    "B1_SCATTP5",
    "B1_SCATTP6",
    "B1_SCATTP7",
    "B1_TRANSPXS",
    "B1_DIFFCOEF",
    "CMM_DIFFCOEF",
    "CMM_TRANSPXS",
    "CMM_DIFFCOEF_X",
    "CMM_TRANSPXS_X",
    "CMM_DIFFCOEF_Y",
    "CMM_TRANSPXS_Y",
    "CMM_DIFFCOEF_Z",
    "CMM_TRANSPXS_Z",
    "TRC_DIFFCOEF",
    "TRC_TRANSPXS",
    "DF_HET_SURF_FLUX",
    "DF_HOM_SURF_FLUX",
    "DF_SURF_DF",
    "DF_HET_CORN_FLUX",
    "DF_HOM_CORN_FLUX",
    "DF_CORN_DF",
    "DF_CORN_IN_CURR",
    "DF_CORN_OUT_CURR",
    "DF_CORN_NET_CURR",
    "DF_HET_VOL_FLUX",
    "DF_HOM_VOL_FLUX",
    "DF_SURF_IN_CURR",
    "DF_SURF_OUT_CURR",
    "DF_SURF_NET_CURR",
    "DF_MID_IN_CURR",
    "DF_MID_OUT_CURR",
    "DF_MID_NET_CURR",
    "DF_SGN_SURF_IN_CURR",
    "DF_SGN_SURF_OUT_CURR",
    "DF_SGN_SURF_NET_CURR",
    "DF_SGN_HET_SURF_FLUX",
    "DF_SGN_HOM_SURF_FLUX",
    "DF_SGN_SURF_DF",
    "PPW_POW",
    "PPW_HOM_FLUX",
    "PPW_FF",
    "PPW_XYZ",
    "ALB_IN_CURR",
    "ALB_OUT_CURR",
    "ALB_TOT_ALB",
    "ALB_PART_ALB",
    "IMP_KEFF",
    "ANA_KEFF",
    "BETA_EFF",
    "LAMBDA",
    "\0"
  };

  /***************************************************************************/

  /***** Make parameter list *************************************************/

  /* Check if user-defined list is given */

  if ((ptr = (long)RDB[DATA_COEF_CALC_PTR_PARAM_LIST]) > VALID_PTR)
    {
      /* Calculate list size */

      np = 0;
      while ((long)RDB[ptr + np] > VALID_PTR)
        np++;

      /* Allocate memory */

      params = (char **)Mem(MEM_ALLOC, np + 1, sizeof(char *));

      for (n = 0; n < np + 1; n++)
        params[n] = (char *)Mem(MEM_ALLOC, MAX_STR, sizeof(char));

      /* Read values */

      for (n = 0; n < np; n++)
        strcpy(params[n], GetText(ptr++));

      /* Put terminator */

      *params[np] = '\0';
    }
  else
    {
      /* Use default list */

      params = default_params;
    }

  /***************************************************************************/

  /***** Print header data ***************************************************/

  /* Check pointer to coefficients */

  if ((loc0 = (long)RDB[DATA_PTR_COEF0]) < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check pointer to gc universes */

  if ((gcu = (long)RDB[DATA_PTR_GCU0]) < VALID_PTR)
    return;

  /* Open file for writing */

  sprintf(tmpstr, "%s.coe", GetText(DATA_PTR_INPUT_FNAME));
  if ((fp = fopen(tmpstr, "a")) == NULL)
    Error(loc0, "Unable to open file for writing");

  /* Print index, total number of runs and number of universes */

  fprintf(fp, "%ld %ld %ld %ld %ld\n", (long)RDB[DATA_COEF_CALC_RUN_IDX],
          (long)RDB[DATA_COEF_CALC_TOT_RUNS], (long)RDB[DATA_COEF_CALC_IDX],
          (long)RDB[DATA_TOT_COEF_CALC], ListSize(gcu));

  /* Pointer to branches */

  loc1 = (long)RDB[loc0 + COEF_PTR_MTX];
  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

  /* Print dimensions */

  fprintf(fp, "%ld", ListSize(loc1));

  /* Reset number of variables */

  nv = 0;

  /* Loop over matrix and print branches */

  while (loc1 > VALID_PTR)
    {
      /* Pointer to branch name */

      ptr = (long)RDB[loc1 + COEF_MTX_PTR_BRA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Print */

      fprintf(fp, " %s", GetText(ptr));

      /* Add to number of variables */

      if ((ptr = (long)RDB[loc1 + COEF_MTX_PTR_VAR]) > VALID_PTR)
        nv = nv + ListSize(ptr);

      /* Next */

      loc1 = NextItem(loc1);
    }

  /* Include history variables */

  if ((loc1 = (long)RDB[DATA_PTR_CASEMTX0]) > VALID_PTR)
    {
      /* Get pointer to branch object setting history */

      loc1 = (long)RDB[loc1 + CASEMTX_PTR_HIS];

      /* Add to number of variables */

      if ((ptr = (long)RDB[loc1 + DEP_BRA_PTR_VAR]) > VALID_PTR)
        nv = nv + ListSize(ptr);
    }

  /* Print newline */

  fprintf(fp, "\n");

  /* Print number of variables */

  fprintf(fp, "%ld ", nv + 3);

  /* Print version, date and time */

  t = time(NULL);
  tm = *localtime(&t);

  fprintf(fp, "VERSION %s DATE %d/", CODE_VERSION, tm.tm_year - 100);

  if (tm.tm_mon + 1 < 10)
    fprintf(fp, "0%d/", tm.tm_mon + 1);
  else
    fprintf(fp, "%d/", tm.tm_mon + 1);

  if (tm.tm_mday < 10)
    fprintf(fp, "0%d ", tm.tm_mday);
  else
    fprintf(fp, "%d ", tm.tm_mday);

  fprintf(fp, "TIME ");

  if (tm.tm_hour < 10)
    fprintf(fp, "0%d:", tm.tm_hour);
  else
    fprintf(fp, "%d:", tm.tm_hour);

  if (tm.tm_min < 10)
    fprintf(fp, "0%d:", tm.tm_min);
  else
    fprintf(fp, "%d:", tm.tm_min);

  if (tm.tm_sec < 10)
    fprintf(fp, "0%d", tm.tm_sec);
  else
    fprintf(fp, "%d", tm.tm_sec);

  /* Loop over matrix and print variables */

  loc1 = (long)RDB[loc0 + COEF_PTR_MTX];
  while (loc1 > VALID_PTR)
    {
      /* Loop over variables */

      ptr = (long)RDB[loc1 + COEF_MTX_PTR_VAR];
      while (ptr > VALID_PTR)
        {
          /* Print name and value */

          fprintf(fp, " %s %s", GetText(ptr + DEP_BRA_VAR_PTR_NAME),
                  GetText(ptr + DEP_BRA_VAR_PTR_VALUE));

          /* Next variable */

          ptr = NextItem(ptr);
        }

      /* Next */

      loc1 = NextItem(loc1);
    }

  /* Print history variables */

  if ((loc1 = (long)RDB[DATA_PTR_CASEMTX0]) > VALID_PTR)
    {
      /* Get pointer to branch object setting history */

      loc1 = (long)RDB[loc1 + CASEMTX_PTR_HIS];

      /* Loop over variables */

      ptr = (long)RDB[loc1 + DEP_BRA_PTR_VAR];
      while (ptr > VALID_PTR)
        {
          /* Print name and value */

          fprintf(fp, " %s %s", GetText(ptr + DEP_BRA_VAR_PTR_NAME),
                  GetText(ptr + DEP_BRA_VAR_PTR_VALUE));

          /* Next variable */

          ptr = NextItem(ptr);
        }
    }

  /* Print newline */

  fprintf(fp, "\n");

  /* Print burnup, index and total number of points */

  fprintf(fp, "%s %ld %ld\n", GetText(DATA_PTR_COEF_BU_PT),
          (long)RDB[DATA_COEF_CALC_BU_IDX], (long)RDB[DATA_TOT_COEF_BU]);

  /***************************************************************************/

  /***** Print output ********************************************************/

  /* Count number of parameters */

  np = 0;
  while (*params[np] != '\0')
    np++;

  /* Loop over gc universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Check if associated with micro-depletion data */

      if ((loc0 = (long)RDB[gcu + GCU_PTR_MDEP]) > VALID_PTR)
        nr = (long)RDB[loc0 + MDEP_N_REA] + 1;
      else
        nr = 0;

      /* Pointer to universe */

      ptr = (long)RDB[gcu + GCU_PTR_UNIV];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Print universe name and number of output parameters */

      fprintf(fp, "%s %ld\n", GetText(ptr + UNIVERSE_PTR_NAME), np + nr);

      /* Loop over parameters */

      n = 0;
      while (*params[n] != '\0')
        {
          /* Find parameter */

          loc1 = (long)RDB[gcu + GCU_PTR_FIRST_STAT];
          while ((loc1 > VALID_PTR) &&
                 (loc1 <= (long)RDB[gcu + GCU_PTR_LAST_STAT]))
            {
              /* Compare name */

              if (!strcmp(params[n], GetText(loc1 + SCORE_PTR_NAME)))
                break;

              /* Next */

              loc1 = NextItem(loc1);
            }

          /* Check pointer and reset */

          if ((loc1 < VALID_PTR) ||
              (loc1 > (long)RDB[gcu + GCU_PTR_LAST_STAT]))
            {
              /* Loop over other stats */

              loc1 = (long)RDB[DATA_PTR_SCORE0];
              while (loc1 > VALID_PTR)
                {
                  /* Compare pointer */

                  if (loc1 > (long)RDB[DATA_LAST_GLOBAL_STAT])
                    {
                      /* Reset pointer */

                      loc1 = -1;

                      /* Break loop */

                      break;
                    }

                  /* Compare name */

                  if (!strcmp(params[n], GetText(loc1 + SCORE_PTR_NAME)))
                    break;

                  /* Next */

                  loc1 = NextItem(loc1);
                }

              /* Check pointer */

              if (loc1 < VALID_PTR)
                {
                  /* Loop over detectors */

                  det = (long)RDB[DATA_PTR_DET0];
                  while (det > VALID_PTR)
                    {
                      /* Check name */

                      if (!strcmp(params[n], GetText(det + DET_PTR_NAME)))
                        {
                          /* Set pointer */

                          loc1 = (long)RDB[det + DET_PTR_STAT];

                          /* Break loop */

                          break;
                        }

                      /* Next detector */

                      det = NextItem(det);
                    }
                }

              /* Check pointer */

              if (loc1 < VALID_PTR)
                {
                  /* Print name and zero values */

                  fprintf(fp, "%s 0\n", params[n]);

                  /* Next parameter */

                  n++;

                  /* Cycle loop */

                  continue;
                }
            }

          /* Print data */

          PrintCoeVals(fp, loc1);

          /* Next parameter */

          n++;
        }

      /* Print micro-depletion data */

      if ((loc0 = (long)RDB[gcu + GCU_PTR_MDEP]) > VALID_PTR)
        {
          /* Number of reactions */

          nr = (long)RDB[loc0 + MDEP_N_REA];
          CheckValue(FUNCTION_NAME, "nr", "", nr, 1, 100000);

          /* Pointer to search keys */

          loc1 = (long)RDB[loc0 + MDEP_PTR_KEY];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Pointer to average densities */

          loc2 = RDB[gcu + GCU_RES_FG_MICRO_DEP_ADENS];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Number of energy groups */

          ntot = (long)RDB[DATA_ERG_FG_NG];
          CheckValue(FUNCTION_NAME, "ntot", "", ntot, 1, 100000);

          /* Pointer to flux */

          ptr = (long)RDB[gcu + GCU_INF_FLX];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Print */

          if ((long)RDB[DATA_COEF_CALC_INCLUDE_ERRORS] == (double)NO)
            {
              /* No errors */

              fprintf(fp, "MDEP_FLX %ld", ntot);

              for (j = 0; j < ntot; j++)
                fprintf(fp, " %12.5E", Mean(ptr, j)/RDB[loc0 + MDEP_VOLUME]);

              fprintf(fp, "\n");
            }
          else
            {
              /* Include errors */

              fprintf(fp, "MDEP_FLX %ld", 2*ntot);

              for (j = 0; j < ntot; j++)
                {
                  if (RelErr(ptr, j) < 0.00010)
                    fprintf(fp, " %11.5E %7.1E",
                            Mean(ptr, j)/RDB[loc0 + MDEP_VOLUME],
                            RelErr(ptr, j));
                  else
                    fprintf(fp, " %11.5E %7.5f",
                            Mean(ptr, j)/RDB[loc0 + MDEP_VOLUME],
                            RelErr(ptr, j));
                }

              fprintf(fp, "\n");
            }

          /* Pointer to cross sections */

          ptr = (long)RDB[gcu + GCU_RES_FG_MICRO_DEP_XS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over reactions */

          for (i = 0; i < nr; i++)
            {
              /* Get ZAI, MT and I */

              ZAI = (long)(RDB[loc1 + i]/10000.0);
              MT = (long)((RDB[loc1 + i] - (double)ZAI*10000.0)/10.0);
              I = (long)RDB[loc1 + i] - 10000*ZAI - 10*MT;

              /* Adjust I if fission */

              if (((MT == 18) || (MT == 19)) && (I > 0))
                I = I - 2;

              /* Check if errors are included and print results */

              if ((long)RDB[DATA_COEF_CALC_INCLUDE_ERRORS] == (double)NO)
                {
                  /* No errors */

                  fprintf(fp, "MDEP %ld %6ld %3ld %ld %11.5E",
                          ntot + 4, ZAI, MT, I, Mean(loc2, i));

                  for (j = 0; j < ntot; j++)
                    fprintf(fp, " %12.5E", Mean(ptr, i, j));
                }
              else
                {
                  /* Include errors */

                  if (RelErr(loc2, i) < 0.00010)
                    fprintf(fp, "MDEP %ld %6ld %3ld %ld %11.5E %7.1E",
                            2*ntot + 4, ZAI, MT, I, Mean(loc2, i),
                            RelErr(loc2, i));
                  else
                    fprintf(fp, "MDEP %ld %6ld %3ld %ld %11.5E %7.5f",
                            2*ntot + 4, ZAI, MT, I, Mean(loc2, i),
                            RelErr(loc2, i));

                  for (j = 0; j < ntot; j++)
                    {
                      if (RelErr(ptr, i, j) < 0.00010)
                        fprintf(fp, " %12.5E %7.1E", Mean(ptr, i, j),
                                RelErr(ptr, i, j));
                      else
                        fprintf(fp, " %12.5E %7.5f", Mean(ptr, i, j),
                                RelErr(ptr, i, j));
                    }
                }

              fprintf(fp, "\n");
            }
        }

      /* Next universe */

      gcu = NextItem(gcu);
    }

  /* Free allocated memory */

  if ((long)RDB[DATA_COEF_CALC_PTR_PARAM_LIST] > VALID_PTR)
    {
      for (n = 0; n < np + 1; n++)
        Mem(MEM_FREE, params[n]);

      Mem(MEM_FREE, params);
    }

  /***************************************************************************/

  /* Close file */

  fclose(fp);
}

/*****************************************************************************/
