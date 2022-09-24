/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writeicmdata.c                                 */
/*                                                                           */
/* Created:       2013/09/15 (JLe)                                           */
/* Last modified: 2013/11/02 (JLe)                                           */
/* Version:       2.1.16                                                     */
/*                                                                           */
/* Description: Writes data for ICM calculation                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteICMData:"

/*****************************************************************************/

void WriteICMData()
{
  long idx, nseg, nsub, nmu0, nmu1, nmu2, nmua, nmus, ng0, ng1, np, icm, ptr;
  long n0, ma0, ms0, g0, n1, ma1, ms1, g1, sz, i;
  double val;
  FILE *fp;
  char tmpstr[MAX_STR];

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check if ICM data is calculated */

  if ((long)RDB[DATA_ICM_CALC] == NO)
    return;

  /* Check if in active cycles */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Open file for writing */

  if ((long)RDB[DATA_ICM_PTR_OUTFILE] > VALID_PTR)
    sprintf(tmpstr, "%s", GetText(DATA_ICM_PTR_OUTFILE));
  else
    sprintf(tmpstr, "%s.icm", GetText(DATA_PTR_INPUT_FNAME));

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Warn(FUNCTION_NAME, "Unable to open file for writing");

  /* Get number of energy groups, segments and angular bins */

  ng0 = (long)RDB[DATA_ICM_NG0];
  ng1 = (long)RDB[DATA_ICM_NG1];
  nseg = (long)RDB[DATA_ICM_NSEG];
  nsub = (long)RDB[DATA_ICM_NSUB];
  nmu0 = (long)RDB[DATA_ICM_NMU0];
  nmu1 = (long)RDB[DATA_ICM_NMU1];
  nmu2 = (long)RDB[DATA_ICM_NMU2];

  /* Calculate number of asymmetric and symmetric bins */

  nmua = nmu1;
  nmus = nmu0*nmu2;

  /* Reset index */

  idx = 1;

  /* Loop over data */

  icm = (long)RDB[DATA_PTR_ICM0];
  while (icm > VALID_PTR)
    {
      /* Get number of pins */

      np = (long)RDB[icm + ICM_NP];

      /**********************************************************************/

      /***** Header data ****************************************************/

      /* Write keyword */

      sprintf(tmpstr, "HDR");
      sz = strlen(tmpstr);
      fwrite(&sz, sizeof(long), 1, fp);
      fwrite(tmpstr, sizeof(char), sz, fp);

      /* Write id */

      sprintf(tmpstr, "%s", GetText(icm + ICM_PTR_ID));
      sz = strlen(tmpstr);
      fwrite(&sz, sizeof(long), 1, fp);
      fwrite(tmpstr, sizeof(char), sz, fp);

      /* Write input file name */

      sprintf(tmpstr, "%s", GetText(DATA_PTR_INPUT_FNAME));
      sz = strlen(tmpstr);
      fwrite(&sz, sizeof(long), 1, fp);
      fwrite(tmpstr, sizeof(char), sz, fp);

      /* write number of groups, segments, angular bins and pins */

      fwrite(&ng0, sizeof(long), 1, fp);
      fwrite(&ng1, sizeof(long), 1, fp);
      fwrite(&nseg, sizeof(long), 1, fp);
      fwrite(&nsub, sizeof(long), 1, fp);
      fwrite(&nmu0, sizeof(long), 1, fp);
      fwrite(&nmu1, sizeof(long), 1, fp);
      fwrite(&nmu2, sizeof(long), 1, fp);
      fwrite(&nmua, sizeof(long), 1, fp);
      fwrite(&nmus, sizeof(long), 1, fp);
      fwrite(&np, sizeof(long), 1, fp);

      /* Write subdivision */

      ptr = (long)RDB[DATA_ICM_PTR_SUB];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fwrite(&RDB[ptr], sizeof(double), nsub + 1, fp);

      /* Write angular bins */
      
      ptr = (long)RDB[DATA_ICM_PTR_MU0];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fwrite(&RDB[ptr], sizeof(double), nmu0 + 1, fp);

      ptr = (long)RDB[DATA_ICM_PTR_MU1];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fwrite(&RDB[ptr], sizeof(double), nmu1 + 1, fp);

      ptr = (long)RDB[DATA_ICM_PTR_MU2];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fwrite(&RDB[ptr], sizeof(double), nmu2 + 1, fp);

      /* Write groups */
      
      ptr = (long)RDB[DATA_ICM_PTR_ENE0];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fwrite(&RDB[ptr], sizeof(double), ng0 + 1, fp);

      if ((ptr = (long)RDB[DATA_ICM_PTR_ENE1]) > VALID_PTR)
        {
          ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          fwrite(&RDB[ptr], sizeof(double), ng1 + 1, fp);
        }

      /***********************************************************************/

      /***** Write coupling coefficients ("type 1") **************************/

      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_CC1]) > VALID_PTR)
        {
          /* Write keyword */

          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus*ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  for (n1 = 0; n1 < nseg; n1++)
                    for (ma1 = 0; ma1 < nmua; ma1++)
                      for (ms1 = 0; ms1 < nmus; ms1++)
                        for (g1 = 0; g1 < ng0; g1++)
                          {
                            val = Mean(ptr, n0, ma0, ms0, g0, n1, 
                                       ma1, ms1, g1);
                            fwrite(&val, sizeof(double), 1, fp);
                          }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  for (n1 = 0; n1 < nseg; n1++)
                    for (ma1 = 0; ma1 < nmua; ma1++)
                      for (ms1 = 0; ms1 < nmus; ms1++)
                        for (g1 = 0; g1 < ng0; g1++)
                          {
                            val = RelErr(ptr, n0, ma0, ms0, g0, n1, 
                                         ma1, ms0, g1);
                            fwrite(&val, sizeof(double), 1, fp);
                          }
        }
      
      /***********************************************************************/

      /***** Write coupling coefficients ("type 2") **************************/

      /* Get pointer */
      
      if ((ptr = (long)RDB[icm + ICM_RES_CC2]) > VALID_PTR)
        {
          /* Write keyword */

          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus*ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  for (n1 = 0; n1 < nseg; n1++)
                    for (ma1 = 0; ma1 < nmua; ma1++)
                      for (ms1 = 0; ms1 < nmus; ms1++)
                        for (g1 = 0; g1 < ng0; g1++)
                          {
                            val = Mean(ptr, n0, ma0, ms0, g0, n1, 
                                       ma1, ms1, g1);
                            fwrite(&val, sizeof(double), 1, fp);
                          }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  for (n1 = 0; n1 < nseg; n1++)
                    for (ma1 = 0; ma1 < nmua; ma1++)
                      for (ms1 = 0; ms1 < nmus; ms1++)
                        for (g1 = 0; g1 < ng0; g1++)
                          {
                            val = RelErr(ptr, n0, ma0, ms0, g0, n1, 
                                         ma1, ms0, g1);
                            fwrite(&val, sizeof(double), 1, fp);
                          }
        }
      
      /***********************************************************************/
      
      /***** Write assembly flux reconstruction factors ("type 1") ***********/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_AFLX1]) > VALID_PTR)
        {
          /* Write keyword */

          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus*ng1;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  for (g1 = 0; g1 < ng1; g1++)
                    {
                      val = Mean(ptr, n0, ma0, ms0, g0, g1);
                      fwrite(&val, sizeof(double), 1, fp);
                    }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  for (g1 = 0; g1 < ng1; g1++)
                    {
                      val = RelErr(ptr, n0, ma0, ms0, g0, g1);
                      fwrite(&val, sizeof(double), 1, fp);
                    }
        }

      /***********************************************************************/

      /***** Write assembly flux reconstruction factors ("type 2") ***********/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_AFLX2]) > VALID_PTR)
        {      
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus*ng1;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  for (g1 = 0; g1 < ng1; g1++)
                    {
                      val = Mean(ptr, n0, ma0, ms0, g0, g1);
                      fwrite(&val, sizeof(double), 1, fp);
                    }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  for (g1 = 0; g1 < ng1; g1++)
                    {
                      val = RelErr(ptr, n0, ma0, ms0, g0, g1);
                      fwrite(&val, sizeof(double), 1, fp);
                    }
        }

      /***********************************************************************/

      /***** Write source rate reconstruction factors ("type 1") *************/
      
      /* Get pointer */
      
      if ((ptr = (long)RDB[icm + ICM_RES_ASRC1]) > VALID_PTR)
        {      
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
      
      /***********************************************************************/
      
      /***** Write source rate reconstruction factors ("type 2") *************/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_ASRC2]) > VALID_PTR)
        {
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
  
      /***********************************************************************/
      
      /***** Write fission rate reconstruction factors ("type 1") ************/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_AFISS1]) > VALID_PTR)
        {      
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
      
      /***********************************************************************/

      /***** Write fission rate reconstruction factors ("type 2") ************/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_AFISS2]) > VALID_PTR)
        {
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
  
      /***********************************************************************/

      /***** Write absorption rate reconstruction factors ("type 1") *********/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_AABS1]) > VALID_PTR)
        {      
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
      
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
  
      /***********************************************************************/
      
      /***** Write absorption rate reconstruction factors ("type 2") *********/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_AABS2]) > VALID_PTR)
        {
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
      
      /***********************************************************************/

      /***** Write assembly power reconstruction factors ("type 1") **********/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_APOW1]) > VALID_PTR)
        {      
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
  
      /***********************************************************************/
      
      /***** Write assembly power reconstruction factors ("type 2") **********/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_APOW2]) > VALID_PTR)
        {
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
  
      /***********************************************************************/

      /***** Write leak rate reconstruction factors (type 1) *****************/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_LEAK1]) > VALID_PTR)
        {      
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
  
      /***********************************************************************/

      /***** Write leak rate reconstruction factors (type 2) *****************/
      
      /* Get pointer */

      if ((ptr = (long)RDB[icm + ICM_RES_LEAK2]) > VALID_PTR)
        {
          /* Write keyword */
          
          sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
          sz = strlen(tmpstr);
          fwrite(&sz, sizeof(long), 1, fp);
          fwrite(tmpstr, sizeof(char), sz, fp);
          
          /* Write number of values */
          
          sz = ng0*nseg*nmua*nmus;
          fwrite(&sz, sizeof(long), 1, fp);
          
          /* Write values and errors */
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = Mean(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
          
          for (n0 = 0; n0 < nseg; n0++)
            for (ma0 = 0; ma0 < nmua; ma0++)
              for (ms0 = 0; ms0 < nmus; ms0++)
                for (g0 = 0; g0 < ng0; g0++)
                  {
                    val = RelErr(ptr, n0, ma0, ms0, g0);
                    fwrite(&val, sizeof(double), 1, fp);
                  }
        }
  
      /***********************************************************************/

      /***** Write pin power reconstruction factors ("type 1") ***************/

      /* Check number of pins */

      if (np > 0)
        {
          /* Get pointer */

          if ((ptr = (long)RDB[icm + ICM_RES_PPOW1]) > VALID_PTR)
            {          
              /* Write keyword */
              
              sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
              sz = strlen(tmpstr);
              fwrite(&sz, sizeof(long), 1, fp);
              fwrite(tmpstr, sizeof(char), sz, fp);
              
              /* Write number of values */
              
              sz = ng0*nseg*nmua*nmus*np;
              fwrite(&sz, sizeof(long), 1, fp);
              
              /* Write values and errors */
              
              for (n0 = 0; n0 < nseg; n0++)
                for (ma0 = 0; ma0 < nmua; ma0++)
                  for (ms0 = 0; ms0 < nmus; ms0++)
                    for (g0 = 0; g0 < ng0; g0++)
                      for (i = 0; i < np; i++)
                        {
                          val = Mean(ptr, n0, ma0, ms0, g0, i);
                          fwrite(&val, sizeof(double), 1, fp);
                        }
              
              for (n0 = 0; n0 < nseg; n0++)
                for (ma0 = 0; ma0 < nmua; ma0++)
                  for (ms0 = 0; ms0 < nmus; ms0++)
                    for (g0 = 0; g0 < ng0; g0++)
                      for (i = 0; i < np; i++)
                        {
                          val = RelErr(ptr, n0, ma0, ms0, g0, i);
                          fwrite(&val, sizeof(double), 1, fp);
                        }
            }
        }
          
      /***********************************************************************/

      /***** Write pin power reconstruction factors ("type 2") ***************/

      /* Check number of pins */

      if (np > 0)
        {      
          /* Get pointer */
          
          if ((ptr = (long)RDB[icm + ICM_RES_PPOW2]) > VALID_PTR)
            {
              /* Write keyword */
          
              sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
              sz = strlen(tmpstr);
              fwrite(&sz, sizeof(long), 1, fp);
              fwrite(tmpstr, sizeof(char), sz, fp);
              
              /* Write number of values */
              
              sz = ng0*nseg*nmua*nmus*np;
              fwrite(&sz, sizeof(long), 1, fp);
              
              /* Write values and errors */
              
              for (n0 = 0; n0 < nseg; n0++)
                for (ma0 = 0; ma0 < nmua; ma0++)
                  for (ms0 = 0; ms0 < nmus; ms0++)
                    for (g0 = 0; g0 < ng0; g0++)
                      for (i = 0; i < np; i++)
                        {
                          val = Mean(ptr, n0, ma0, ms0, g0, i);
                          fwrite(&val, sizeof(double), 1, fp);
                        }
              
              for (n0 = 0; n0 < nseg; n0++)
                for (ma0 = 0; ma0 < nmua; ma0++)
                  for (ms0 = 0; ms0 < nmus; ms0++)
                    for (g0 = 0; g0 < ng0; g0++)
                      for (i = 0; i < np; i++)
                        {
                          val = RelErr(ptr, n0, ma0, ms0, g0, i);
                          fwrite(&val, sizeof(double), 1, fp);
                        }
            }
        }
    
      /***********************************************************************/

      /***** Write pin flux reconstruction factors ("type 1") ****************/
      
      /* Check number of pins */

      if (np > 0)
        {
          /* Get pointer */
          
          if ((ptr = (long)RDB[icm + ICM_RES_PFLX1]) > VALID_PTR)
            {          
              /* Write keyword */
              
              sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
              sz = strlen(tmpstr);
              fwrite(&sz, sizeof(long), 1, fp);
              fwrite(tmpstr, sizeof(char), sz, fp);
              
              /* Write number of values */
              
              sz = ng0*nseg*nmua*nmus*ng1*np;
              fwrite(&sz, sizeof(long), 1, fp);
              
              /* Write values and errors */
              
              for (n0 = 0; n0 < nseg; n0++)
                for (ma0 = 0; ma0 < nmua; ma0++)
                  for (ms0 = 0; ms0 < nmus; ms0++)
                    for (g0 = 0; g0 < ng0; g0++)
                      for (g1 = 0; g1 < ng1; g1++)
                        for (i = 0; i < np; i++)
                          {
                            val = Mean(ptr, n0, ma0, ms0, g0, g1, i);
                            fwrite(&val, sizeof(double), 1, fp);
                          }
              
              for (n0 = 0; n0 < nseg; n0++)
                for (ma0 = 0; ma0 < nmua; ma0++)
                  for (ms0 = 0; ms0 < nmus; ms0++)
                    for (g0 = 0; g0 < ng0; g0++)
                      for (g1 = 0; g1 < ng1; g1++)
                        for (i = 0; i < np; i++)
                          {
                            val = RelErr(ptr, n0, ma0, ms0, g0, g1, i);
                            fwrite(&val, sizeof(double), 1, fp);
                          }
            }
        }
      
      /***********************************************************************/
      
      /***** Write pin flux reconstruction factors ("type 2") ****************/

      /* Check number of pins */

      if (np > 0)
        {      
          /* Get pointer */
          
          if ((ptr = (long)RDB[icm + ICM_RES_PFLX2]) > VALID_PTR)
            {          
              /* Write keyword */
              
              sprintf(tmpstr, "%s", GetText(ptr + SCORE_PTR_NAME));
              sz = strlen(tmpstr);
              fwrite(&sz, sizeof(long), 1, fp);
              fwrite(tmpstr, sizeof(char), sz, fp);
              
              /* Write number of values */
              
              sz = ng0*nseg*nmua*nmus*ng1*np;
              fwrite(&sz, sizeof(long), 1, fp);
              
              /* Write values and errors */
              
              for (n0 = 0; n0 < nseg; n0++)
                for (ma0 = 0; ma0 < nmua; ma0++)
                  for (ms0 = 0; ms0 < nmus; ms0++)
                    for (g0 = 0; g0 < ng0; g0++)
                      for (g1 = 0; g1 < ng1; g1++)
                        for (i = 0; i < np; i++)
                          {
                            val = Mean(ptr, n0, ma0, ms0, g0, g1, i);
                            fwrite(&val, sizeof(double), 1, fp);
                          }
              
              for (n0 = 0; n0 < nseg; n0++)
                for (ma0 = 0; ma0 < nmua; ma0++)
                  for (ms0 = 0; ms0 < nmus; ms0++)
                    for (g0 = 0; g0 < ng0; g0++)
                      for (g1 = 0; g1 < ng1; g1++)
                        for (i = 0; i < np; i++)
                          {
                            val = RelErr(ptr, n0, ma0, ms0, g0, g1, i);
                            fwrite(&val, sizeof(double), 1, fp);
                          }
            }
        }
        
      /***********************************************************************/

      /* Next */

      icm = NextItem(icm);
    }

  /* Close file and exit */

  fclose(fp);
}

/*****************************************************************************/
