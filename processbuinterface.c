/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processbuinterface.c                           */
/*                                                                           */
/* Created:       2020/04/20 (JLe)                                           */
/* Last modified: 2020/04/20 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reads materials from restart files for burnup interface      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessBUInterface:"

/*****************************************************************************/

void ProcessBUInterface()
{
  long sz, n, nnuc, zai, idx;
  double bu, days, days0, adens, mdens, mbu;
  char name[MAX_STR], fname[MAX_STR];
  FILE *fp;
  return;
  /* Check domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    Die(FUNCTION_NAME, "Domain decomposition in use");

  /* File name */

  if ((long)RDB[DATA_RESTART_READ_PTR_FNAME] > VALID_PTR)
    sprintf(fname, "%s", "tst_1.wrk");
  else
    sprintf(fname, "%s.wrk", "tst_1");

  /* Open file for reading */

  if ((fp = fopen(fname, "r")) == NULL)
    Error(0, "Restart file \"%s\" does not exist", fname);

  /* Print */

  fprintf(outp,
          "Reading material compositions from restart file \"%s\":\n\n",
          fname);

  /* Read loop */

  while ((sz = fread(&n, sizeof(long), 1, fp)) > 0)
    {
      /* Read name */

      if ((sz = fread(name, sizeof(char), n, fp)) == 0)
        Error(0, "Error in restart file");

      /* Put EOF */

      name[n] = '\0';

      printf("%s\n", name);

      /* Read nominal burnup and time */

      if ((sz = fread(&bu, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&days, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      /* Update index if new day */

      if (days != days0)
        idx++;

      /* Read number of nuclides, atomic and mass density and burnup */

      if ((sz = fread(&nnuc, sizeof(long), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&mdens, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      if ((sz = fread(&mbu, sizeof(double), 1, fp)) == 0)
        Error(0, "Error in restart file");

      /* Check values (mdens is not used) */

      CheckValue(FUNCTION_NAME, "nnuc", "", nnuc, 1, 10000);
      CheckValue(FUNCTION_NAME, "adens", "", adens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mdens", "", mdens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mbu", "", mbu, 0.0, 1000.0);

      /* Loop over composition */

      for (n = 0; n < nnuc; n++)
        {
          /* Read ZAI and atomic density */

          if ((sz = fread(&zai, sizeof(long), 1, fp)) == 0)
            Error(0, "Error in restart file");

          if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
            Error(0, "Error in restart file");

          /* Check if Xenon is set to zero */

          if ((long)RDB[DATA_RESTART_READ_ZERO_XE] == YES)
            if ((zai == 541350) || (zai == 541351))
              adens = 0.0;

          /* Check if Samarium is set to zero */

          if ((long)RDB[DATA_RESTART_READ_ZERO_SM] == YES)
            if (zai == 621490)
              adens = 0.0;
        }
    }

  fprintf(outp, "OK.\n\n");

  exit(0);
}

/*****************************************************************************/
