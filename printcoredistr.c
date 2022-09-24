/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printcoredistr.c                               */
/*                                                                           */
/* Created:       2011/05/15 (JLe)                                           */
/* Last modified: 2012/09/07 (JLe)                                           */
/* Version:       2.1.8                                                      */
/*                                                                           */
/* Description: Prints core power distribution                               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintCoreDistr:"

/*****************************************************************************/

void PrintCoreDistr()
{
  long ptr, n0, n1, n2, m0, m1, m2, m5, m10, i;
  double val, err, prof[101], cum[101];
  char outfile[MAX_STR];
  FILE *fp;

  /* Check that distribution is called */
  
  if (RDB[DATA_CORE_PDE_DEPTH] < 1.0)  
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) 
    return;
  
  /* Set file name */

  sprintf(outfile, "%s_core%ld.m", GetText(DATA_PTR_INPUT_FNAME),
          (long)RDB[DATA_BURN_STEP]);
          
  /* Open file */
  
  if ((fp = fopen(outfile, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /***************************************************************************/

  /***** Core level **********************************************************/

  if ((ptr = (long)RDB[DATA_CORE_PDE_PTR_RES0]) > VALID_PTR)
    {
      /* Reset counters */

      m0 = 0;
      m1 = 0;
      m2 = 0;
      m5 = 0;
      m10 = 0;
      
      for (i = 0; i < 100; i++)
        {
          prof[i] = 0.0;
          cum[i] = 0.0;
        }

      /* Print */

      fprintf(fp, "%% ----- Distribution at level 0 (first lattice):\n\n");
      
      fprintf(fp, "lvl0 = [\n");

      for (n0 = 0; n0 < (long)RDB[DATA_CORE_PDE_N0]; n0++)
        {
          val = Mean(ptr, n0);
          err = RelErr(ptr, n0);

          if ((val > 0.0) && (err > 0.0) && (err < 1.0))
            {
              fprintf(fp, " %3ld %1.5E %8.5f %% a\n", n0 + 1, val, err);

              m0++;

              /* Compare limits */

              if (err < 0.01)
                m1++;
              if (err < 0.02)
                m2++;
              if (err < 0.05)
                m5++;
              if (err < 0.10)
                m10++;
              
              /* Add to profile */
              
              if (err < 0.2)
                prof[(long)(err*500.0)]++;
            }
        }

      fprintf(fp, "];\n");
      
      fprintf(fp, "\n%% Number and fraction of bins with statistical errors\n");
      fprintf(fp, "%% below 0.01, 0.02, 0.05 and 0.10:\n\n");
      
      fprintf(fp, "n0 = [ %ld %ld %ld %ld ];\n", m1, m2, m5, m10);
      fprintf(fp, "f0 = [ %1.5f %1.5f %1.5f %1.5f ];\n", 
              ((double)m1)/((double)m0), ((double)m2)/((double)m0),
              ((double)m5)/((double)m0), ((double)m10)/((double)m0));

      fprintf(fp, "\n%% Distribution of statistical errors (1%% bins):\n\n");
      
      if (m0 > 0)
        {
          fprintf(fp, "prof0 = [ ");
          
          for (i = 0; i < 100; i++)
            {
              fprintf(fp, "%1.5E ", prof[i]/m0);
              
              if (i == 0)
                cum[i] = prof[i]/m0;
              else
                cum[i] = cum[i - 1] + prof[i]/m0;
            }
          
          fprintf(fp, "];\n");

          fprintf(fp, "\n%% Cumulative distribution of statistical errors ");
          fprintf(fp, "(1%% bins):\n\n");
          
          fprintf(fp, "cum0 = [ ");
          
          for (i = 0; i < 100; i++)
            fprintf(fp, "%1.5E ", cum[i]);
          
          fprintf(fp, "];\n");
        }
    }

  /***************************************************************************/

  /***** Pin level ***********************************************************/

  if ((ptr = (long)RDB[DATA_CORE_PDE_PTR_RES1]) > VALID_PTR)
    {
      /* Reset counters */

      m0 = 0;
      m1 = 0;
      m2 = 0;
      m5 = 0;
      m10 = 0;

      for (i = 0; i < 100; i++)
        {
          prof[i] = 0.0;
          cum[i] = 0.0;
        }
      
      /* Print */

      fprintf(fp, "\n%% ----- Distribution at level 1 (second lattice):\n\n");
      
      fprintf(fp, "lvl1 = [\n");
      
      for (n0 = 0; n0 < (long)RDB[DATA_CORE_PDE_N0]; n0++)
        for (n1 = 0; n1 < (long)RDB[DATA_CORE_PDE_N1]; n1++)
          {
            val = Mean(ptr, n0, n1);
            err = RelErr(ptr, n0, n1);

            if ((val > 0.0) && (err > 0.0) && (err < 1.0))
              {
                fprintf(fp, " %3ld %3ld %1.5E %8.5f %% p\n", n0 + 1, 
                        n1 + 1, val, err);

                m0++;
                
                /* Compare limits */
                
                if (err < 0.01)
                  m1++;
                if (err < 0.02)
                  m2++;
                if (err < 0.05)
                  m5++;
                if (err < 0.10)
                  m10++;

                /* Add to profile */
              
                if (err < 0.2)
                  prof[(long)(err*500.0)]++;
              }
          }

      fprintf(fp, "];\n");

      fprintf(fp, "\n%% Number and fraction of bins with statistical errors\n");
      fprintf(fp, "%% below 0.01, 0.02, 0.05 and 0.10:\n\n");

      fprintf(fp, "n1 = [ %ld %ld %ld %ld ];\n", m1, m2, m5, m10);
      fprintf(fp, "f1 = [ %1.5f %1.5f %1.5f %1.5f ];\n", 
              ((double)m1)/((double)m0), ((double)m2)/((double)m0),
              ((double)m5)/((double)m0), ((double)m10)/((double)m0));

      fprintf(fp, "\n%% Distribution of statistical errors (1%% bins):\n\n");
     
      if (m0 > 0)
        {
          fprintf(fp, "prof1 = [ ");
          
          for (i = 0; i < 100; i++)
            {
              fprintf(fp, "%1.5E ", prof[i]/m0);
              
              if (i == 0)
                cum[i] = prof[i]/m0;
              else
                cum[i] = cum[i - 1] + prof[i]/m0;
            }
          
          fprintf(fp, "];\n");

          fprintf(fp, "\n%% Cumulative distribution of statistical errors ");
          fprintf(fp, "(1%% bins):\n\n");
          
          fprintf(fp, "cum1 = [ ");
          
          for (i = 0; i < 100; i++)
            fprintf(fp, "%1.5E ", cum[i]);
          
          fprintf(fp, "];\n");
        }
    }

  /***************************************************************************/

  /***** Region level ********************************************************/
    
  if ((ptr = (long)RDB[DATA_CORE_PDE_PTR_RES2]) > VALID_PTR)
    {
      /* Reset counters */

      m0 = 0;
      m1 = 0;
      m2 = 0;
      m5 = 0;
      m10 = 0;

      for (i = 0; i < 100; i++)
        {
          prof[i] = 0.0;
          cum[i] = 0.0;
        }
      
      /* Print */

      fprintf(fp, "\n%% ----- Distribution at level 2 (axial):\n\n");

      fprintf(fp, "lvl2 = [\n");

      if ((long)RDB[DATA_CORE_PDE_N1] > 0)
        {
          for (n0 = 0; n0 < (long)RDB[DATA_CORE_PDE_N0]; n0++)
            for (n1 = 0; n1 < (long)RDB[DATA_CORE_PDE_N1]; n1++)
              for (n2 = 0; n2 < (long)RDB[DATA_CORE_PDE_N2]; n2++)
                {
                  val = Mean(ptr, n0, n1, n2);
                  err = RelErr(ptr, n0, n1, n2);
                  
                  if ((val > 0.0) && (err > 0.0) && (err < 1.0))
                    {
                      fprintf(fp, " %3ld %3ld %3ld %1.5E %8.5f %% r\n", n0 + 1, 
                              n1 + 1, n2 + 1, val, err);
                      
                      m0++;
                      
                      /* Compare limits */
                      
                      if (err < 0.01)
                        m1++;
                      if (err < 0.02)
                        m2++;
                      if (err < 0.05)
                        m5++;
                      if (err < 0.10)
                        m10++;
                      
                      /* Add to profile */
                      
                      if (err < 0.2)
                        prof[(long)(err*500.0)]++;
                    }
                }
        }
      else
        {
          for (n0 = 0; n0 < (long)RDB[DATA_CORE_PDE_N0]; n0++)
            for (n2 = 0; n2 < (long)RDB[DATA_CORE_PDE_N2]; n2++)
              {
                val = Mean(ptr, n0, 0, n2);
                err = RelErr(ptr, n0, 0, n2);
                
                if ((val > 0.0) && (err > 0.0) && (err < 1.0))
                  {
                    fprintf(fp, " %3ld   1 %3ld %1.5E %8.5f %% r\n", n0 + 1, 
                            n2 + 1, val, err);
                    
                    m0++;
                    
                    /* Compare limits */
                    
                    if (err < 0.01)
                      m1++;
                    if (err < 0.02)
                      m2++;
                    if (err < 0.05)
                      m5++;
                    if (err < 0.10)
                      m10++;
                    
                    /* Add to profile */
                    
                    if (err < 0.2)
                      prof[(long)(err*500.0)]++;
                  }
              }
        }

      fprintf(fp, "];\n");

      fprintf(fp, "\n%% Number and fraction of bins with statistical errors\n");
      fprintf(fp, "%% below 0.01, 0.02, 0.05 and 0.10:\n\n");

      fprintf(fp, "n2 = [ %ld %ld %ld %ld ];\n", m1, m2, m5, m10);
      fprintf(fp, "f2 = [ %1.5f %1.5f %1.5f %1.5f ];\n", 
              ((double)m1)/((double)m0), ((double)m2)/((double)m0),
              ((double)m5)/((double)m0), ((double)m10)/((double)m0));

      fprintf(fp, "\n%% Distribution of statistical errors (1%% bins):\n\n");

      if (m0 > 0)
        {
          fprintf(fp, "prof2 = [ ");
          
          for (i = 0; i < 100; i++)
            {
              fprintf(fp, "%1.5E ", prof[i]/m0);
              
              if (i == 0)
                cum[i] = prof[i]/m0;
              else
                cum[i] = cum[i - 1] + prof[i]/m0;
            }
          
          fprintf(fp, "];\n");

          fprintf(fp, "\n%% Cumulative distribution of statistical errors ");
          fprintf(fp, "(1%% bins):\n\n");
          
          fprintf(fp, "cum2 = [ ");
          
          for (i = 0; i < 100; i++)
            fprintf(fp, "%1.5E ", cum[i]);
          
          fprintf(fp, "];\n");
        }
    }
  
  /***************************************************************************/

  /* Close file */

  fclose(fp);

  /***************************************************************************/
}

/*****************************************************************************/
