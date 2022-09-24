/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printhistoryoutput.c                           */
/*                                                                           */
/* Created:       2011/05/18 (JLe)                                           */
/* Last modified: 2020/06/25 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Prints history of statistical variables                      */
/*                                                                           */
/* Comments: - Tää toimii nyt vaan yksiarvoisilla muuttujilla                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintHistoryOutput:"

/*****************************************************************************/

void PrintHistoryOutput()
{
  long loc0, ptr, i, n0, n1, n2, n3, n4, dim, bins[5];
  long skip, curb, cycle;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Test print option */

  if ((long)RDB[DATA_OPTI_PRINT_HIS] == NO)
    return;

  /* Check corrector step */

  if (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
      (RDB[DATA_BURN_SIE] == (double)NO))
    return;

  /* Avoid compiler warning */

  for (i = 0; i < 5; i++)
    bins[i] = 0;

  /* Set file name */

  sprintf(tmpstr, "%s_his%ld.m", GetText(DATA_PTR_INPUT_FNAME),
          (long)RDB[DATA_BURN_STEP]);

  /* File name if running SIE burnup calculation */

  if (RDB[DATA_BURN_SIE] == (double)YES)
    {
      if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
        sprintf(tmpstr, "%s_his%ldi%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                (long)RDB[DATA_BURN_STEP], (long)RDB[DATA_BURN_CI_I]);
      else
        sprintf(tmpstr, "%s_his%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                (long)RDB[DATA_BURN_STEP]);
    }

  /* Open file */

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

 /* Fix number of active and inactive cycles if using fission source passing */

  if((RDB[DATA_USE_FSP] == (double)NO) ||
     ((RDB[DATA_BURN_STEP] + RDB[DATA_SOL_REL_ITER] == 0.0)
      && (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)))
    {

      /* Set number of inactive batches */

      skip = (long)RDB[DATA_CRIT_SKIP];

      /* Set current batch number */

      curb = (long)RDB[DATA_CYCLE_IDX] + 1;
    }
  else
    {
      /* Set number of inactive batches                          */

      skip = (long)RDB[DATA_FSP_CRIT_SKIP];

      /* Fix current batch number */

      curb = (long)RDB[DATA_CYCLE_IDX] + skip - (long)RDB[DATA_CRIT_SKIP] + 1;
    }

  /* Print running time */

  ptr = (long)RDB[DATA_BTCH_TIME];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  fprintf(fp, "HIS_TIME = [\n");

  /* Loop over cycles */

  for (cycle = 0; cycle < curb; cycle++)
    {
      /* Check skip cycles */

      if (cycle == skip)
        fprintf(fp, "\n%% ----- Begin active cycles -----\n\n");

      if (cycle < skip)
        fprintf(fp, "%4ld  ", cycle + 1);
      else
        fprintf(fp, "%4ld  ", cycle - skip + 1);

      /* Fix cycle to read from in case of fission source passing */

      i = cycle + (long)RDB[DATA_CRIT_SKIP] - skip;

      fprintf(fp, "%1.5E\n", RDB[ptr + i]);
    }

  fprintf(fp, "];\n\n");

  /* Loop over statistics */

  loc0 = (long)RDB[DATA_PTR_SCORE0];
  while(loc0 > VALID_PTR)
    {
      /* Check history pointer */

      if ((long)RDB[loc0 + SCORE_PTR_HIS] > 0)
        {
          /* Get dimension */

          if ((dim = (long)RDB[loc0 + SCORE_DIM]) > 5)
            Die(FUNCTION_NAME, "Need more dimensions");

          /* Pointer to list of maximum values */

          ptr = (long)RDB[loc0 + SCORE_PTR_NMAX];

          /* Read bin sizes */

          for (i = 0; i < dim; i++)
            {
              /* Get number of bins */

              bins[i] = (long)RDB[ptr++];
            }

          /* Print variable name */

          sprintf(tmpstr, "HIS_%s", GetText(loc0 + SCORE_PTR_NAME));
          fprintf(fp, "%s = [\n", tmpstr);

          /* Loop over cycles */

          for (cycle = 0; cycle < curb; cycle++)
            {
              /* Check skip cycles */

              if (cycle == skip)
                fprintf(fp, "\n%% ----- Begin active cycles -----\n\n");

              if (cycle < skip)
                fprintf(fp, "%4ld  ", cycle + 1);
              else
                fprintf(fp, "%4ld  ", cycle - skip + 1);

              /* Fix cycle to read from in case of fission source passing */

              i = cycle + (long)RDB[DATA_CRIT_SKIP] - skip;

              /* Check number of dimensions and print data */

              if (dim == 1)
                {
                  for (n0 = 0; n0 < bins[0]; n0++)
                    fprintf(fp, "%1.5E %1.5E %1.5f  ",
                            HisVal(loc0, i, n0),
                            HisMean(loc0, i, n0),
                            HisRelErr(loc0, i, n0));
                }
              else if (dim == 2)
                {
                  for (n0 = 0; n0 < bins[0]; n0++)
                    for (n1 = 0; n1 < bins[1]; n1++)
                      fprintf(fp, "%1.5E %1.5E %1.5f  ",
                              HisVal(loc0, i, n0, n1),
                              HisMean(loc0, i, n0, n1),
                              HisRelErr(loc0, i, n0, n1));
                }
              else if (dim == 3)
                {
                  for (n0 = 0; n0 < bins[0]; n0++)
                    for (n1 = 0; n1 < bins[1]; n1++)
                      for (n2 = 0; n2 < bins[2]; n2++)
                        fprintf(fp, "%1.5E %1.5E %1.5f  ",
                                HisVal(loc0, i, n0, n1, n2),
                                HisMean(loc0, i, n0, n1, n2),
                                HisRelErr(loc0, i, n0, n1, n2));
                }
              else if (dim == 4)
                {
                  for (n0 = 0; n0 < bins[0]; n0++)
                    for (n1 = 0; n1 < bins[1]; n1++)
                      for (n2 = 0; n2 < bins[2]; n2++)
                        for (n3 = 0; n3 < bins[3]; n3++)
                          fprintf(fp, "%1.5E %1.5E %1.5f  ",
                                  HisVal(loc0, i, n0, n1, n2, n3),
                                  HisMean(loc0, i, n0, n1, n2, n3),
                                  HisRelErr(loc0, i, n0, n1, n2, n3));
                }
              else if (dim == 5)
                {
                  for (n0 = 0; n0 < bins[0]; n0++)
                    for (n1 = 0; n1 < bins[1]; n1++)
                      for (n2 = 0; n2 < bins[2]; n2++)
                        for (n3 = 0; n3 < bins[3]; n3++)
                          for (n4 = 0; n4 < bins[4]; n4++)
                            fprintf(fp, "%1.5E %1.5E %1.5f  ",
                                    HisVal(loc0, i, n0, n1, n2, n3, n4),
                                    HisMean(loc0, i, n0, n1, n2, n3, n4),
                                    HisRelErr(loc0, i, n0, n1, n2, n3, n4));
                }
              else
                Die(FUNCTION_NAME, "wtf?");

              fprintf(fp, "\n");
            }

          fprintf(fp, "];\n\n");
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Close file */

  fclose(fp);

  /***************************************************************************/
}

/*****************************************************************************/
