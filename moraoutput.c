/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : moraoutput.c                                   */
/*                                                                           */
/* Created:       2013/03/19 (JLe)                                           */
/* Last modified: 2013/03/19 (JLe)                                           */
/* Version:       2.1.13                                                     */
/*                                                                           */
/* Description: Writes homogenized multi-group cross sections for MORA       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MORAOutput:"

/*****************************************************************************/

void MORAOutput()
{
  long loc0, ptr, ng, nc, n, m, i;
  FILE *fp;
  char tmpstr[MAX_STR];

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_MORA0];
  while (loc0 > VALID_PTR)
    {
      /* Open file for writing */
  
      sprintf(tmpstr, "%s", GetText(loc0 + MORA_PTR_FNAME));
      
      if ((fp = fopen(tmpstr, "w")) == NULL)
        Warn(FUNCTION_NAME, "Unable to open file for writing");

      /* Get number of energy groups and cosine bins */

      ng = (long)RDB[loc0 + MORA_N_EG];
      nc = (long)RDB[loc0 + MORA_N_COS];

      /* Print group structure */

      ptr = (long)RDB[loc0 + MORA_PTR_EG];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\nENE %ld\n", ng);

      for (n = ng; n > -1; n--)
        fprintf(fp, "%12.5E\n", RDB[ptr + n]);

      /* Print cosine bins */

      fprintf(fp, "\nCOS %ld\n", nc);

      for (n = nc; n > -1; n--)
        fprintf(fp, "%12.5E\n", 2.0*((double)n)/((double)nc) - 1.0);

      /* Print flux */

      ptr = (long)RDB[loc0 + MORA_PTR_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, n), RelErr(ptr, n));

      /* Print total cross section */

      ptr = (long)RDB[loc0 + MORA_PTR_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, n), RelErr(ptr, n));
      
      /* Print capture cross section */

      ptr = (long)RDB[loc0 + MORA_PTR_CAPT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, n), RelErr(ptr, n));

      /* Print fission cross section */

      ptr = (long)RDB[loc0 + MORA_PTR_FISS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, n), RelErr(ptr, n));

      /* Print fission kappa */

      ptr = (long)RDB[loc0 + MORA_PTR_KAPPA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, n)/MEV, RelErr(ptr, n));

      /* Print prompt fission nubar */

      ptr = (long)RDB[loc0 + MORA_PTR_PNU];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, n), RelErr(ptr, n));

      /* Print delayed fission nubar */

      ptr = (long)RDB[loc0 + MORA_PTR_DNU];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, n), RelErr(ptr, n));

      /* Print prompt fission spectrum */

      ptr = (long)RDB[loc0 + MORA_PTR_CHIP];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, n), RelErr(ptr, n));

      /* Print delayed fission spectrum */

      ptr = (long)RDB[loc0 + MORA_PTR_CHID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, n), RelErr(ptr, n));

      /* Print scattering probability matrix */

      ptr = (long)RDB[loc0 + MORA_PTR_SCATTP];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        for (m = ng - 1; m > -1; m--)
          for (i = nc - 1; i > -1; i--)
            fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, m, n, i), 
                    RelErr(ptr, n, m, i));
      
      /* Print scattering weight matrix */

      ptr = (long)RDB[loc0 + MORA_PTR_SCATTW];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      fprintf(fp, "\n%s\n", GetText(ptr + SCORE_PTR_NAME));

      for (n = ng - 1; n > -1; n--)
        for (m = ng - 1; m > -1; m--)
          for (i = nc - 1; i > -1; i--)
            fprintf(fp, "%12.5E %1.5f\n", Mean(ptr, m, n, i), 
                    RelErr(ptr, n, m, i));
      
      /* Close file */
      
      fclose(fp);
      
      /* Next */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
