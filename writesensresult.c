/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writesensresult.c                              */
/*                                                                           */
/* Created:       2017/05/03 (VVa)                                           */
/* Last modified: 2018/12/05 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Writes sensitivity calculation results from one response to  */
/*              file.                                                        */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteSensResult:"

/*****************************************************************************/

void WriteSensResult(FILE *fp, char *name, long res, long nmat, long nzai,
                     long nrea, long nene, long nmu, long gen)
{
  long imat, izai, irea, iene, imu, idx, sens, skipizai, skipimat, ptr;
  char tmpstr[MAX_STR], tmpstr2[MAX_STR*2];

  /* Get pointer to the sensitivity block or return 0.0 */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  if (name == NULL)
    sprintf(tmpstr, "%s", GetText(res + SCORE_PTR_NAME));
  else
    sprintf(tmpstr, "%s", name);

  fprintf(fp, "%s = [\n", tmpstr);

  /* Open array */

  imu = 0;

  /* Get index to the "lost" zai that should be skipped */

  ptr = (long)RDB[sens + SENS_PTR_ZAI_INDICES];
  skipizai = (long)RDB[ptr + SENS_NON_ZAI_IDX];

  /* Get index to the "lost" mat that should be skipped */

  ptr = (long)RDB[sens + SENS_PTR_MAT_INDICES];
  skipimat = (long)RDB[ptr + SENS_NON_MAT_IDX];

  /* First loop is over materials (skip "not found") */

  for (imat = 0; imat < nmat; imat++)
    {
      if (imat == skipimat)
        continue;

      /* Second loop over ZAIs (skip "not found") */

      for (izai = 0; izai < nzai; izai++)
        {
          if (izai == skipizai)
            continue;

          /* Third loop over reactions */

          for (irea = 0; irea < nrea; irea++)
            {
              /* Final loop over energy bins (skip first, which is total) */

              for (iene = 1; iene < nene; iene++)
                {
                  /* Calculate bin idx */

                  idx = 1 + imat*nzai*nrea*nene*nmu + izai*nrea*nene*nmu + irea*nene*nmu + iene*nmu + imu;

                  /* Print value and estimate for relative error */

                  fprintf(fp, " %12.5E %7.1E", Mean(res, gen, idx), RelErr(res, gen, idx));

                }
            }
        }
    }

  fprintf(fp, "\n];\n\n");

  fprintf(fp, "%s = reshape(%s, [2, SENS_N_ENE, SENS_N_PERT, SENS_N_ZAI, SENS_N_MAT]);\n", tmpstr, tmpstr);
  fprintf(fp, "%s = permute(%s, [5, 4, 3, 2, 1]);\n\n", tmpstr, tmpstr);

  /* Not resolved in energy (total) */

  sprintf(tmpstr2, "%s_E_INT", tmpstr);
  fprintf(fp, "%s = [\n", tmpstr2);

  /* First loop is over materials (skip first, which is "not found") */

  for (imat = 1; imat < nmat; imat++)
    {
      /* Second loop over ZAIs (skip first, which is "not found") */

      for (izai = 1; izai < nzai; izai++)
        {
          /* Third loop over reactions */

          for (irea = 0; irea < nrea; irea++)
            {
              /* First energy bin includes all energies */

              iene = 0;

              /* Calculate bin idx */

              idx = 1 + imat*nzai*nrea*nene*nmu + izai*nrea*nene*nmu + irea*nene*nmu + iene*nmu + imu;

              /* Print value and estimate for relative error */

              fprintf(fp, " %12.5E %7.1E", Mean(res, gen, idx), RelErr(res, gen, idx));
            }
        }
    }

  fprintf(fp, "\n];\n\n");

  fprintf(fp, "%s = reshape(%s, [2, SENS_N_PERT, SENS_N_ZAI, SENS_N_MAT]);\n", tmpstr2, tmpstr2);
  fprintf(fp, "%s = permute(%s, [4, 3, 2, 1]);\n\n", tmpstr2, tmpstr2);

}
