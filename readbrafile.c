/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readbrafile.c                                  */
/*                                                                           */
/* Created:       2011/01/28 (JLe)                                           */
/* Last modified: 2016/11/17 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads energy-dependent isomeric branching ratios             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadBRAFile:"

/*****************************************************************************/

void ReadBRAFile()
{
  long pta, MAT, ZA, ZAI, LIS, NS, NWD, NXC, LFS[20], NR[20], NP[20], INTT[20];
  long mt, n, m, i, loc0, ptr, MT[20], NP3[20], NP10[20], np;
  double *E[20], *f[20], *E3[20], *xs3[20], *E10[20], *xs10[20], *Eg, *xs;
  char line[MAX_STR], *eof;
  FILE *fp;
  fpos_t pos;

  /* Check pointer */
  
  if ((pta = (long)RDB[DATA_PTR_BRADATA_FNAME_LIST]) < 1)
    return;

  /* Loop over files */

  while ((long)RDB[pta] > 0)
    {      
      /* Test format */

      WDB[DATA_DUMMY] = RDB[pta];
      TestDOSFile(GetText(DATA_DUMMY));

      /***********************************************************************/

      /***** Read multiplicities from file 9 *********************************/

      /* Open file for reading */
  
      fp = OpenDataFile(pta, "isomeric branching ratio data");
        
      /* Read data */

      while (1 != 2)
        {
          /*******************************************************************/

          /***** Read data from file *****************************************/

          /* Loop until line 1 of file 9 */
      
          do 
            eof = fgets(line, 82, fp);
          while (((line[71] != '9') || (line[78] != ' ') || (line[79] != '1')) 
                 && (eof != NULL));

          /* Check EOF */
          
          if (eof == NULL)
            break;

          /* Read ZA, LIS and NS */
          
          ZA = (long)rint(ENDFColF(1, line));
          LIS = ENDFColI(3, line);
          NS = ENDFColI(5, line);
          
          /* Check if array sizes are sufficient */

          if (NS > 20)
            Die(FUNCTION_NAME, "Increase array sizes");

          /* Get mt */
          
          mt = atoi(&line[72]);

          /* Convert all isomeric states to stae 1 */

          if (LIS > 0)
            LIS = 1;
          
          /* Calculate ZAI */

          ZAI = 10*ZA + LIS;

          /* Loop over states */
          
          for (n = 0; n < NS; n++)
            {
              /* Next line */
              
              ENDFNewLine(FUNCTION_NAME, line, fp);
              
              /* Read LFS, NR and NP */
              
              LFS[n] = ENDFColI(4, line);
              NR[n] = ENDFColI(5, line);
              NP[n] = ENDFColI(6, line);

              /* Allocate memory for data */
              
              E[n] = Mem(MEM_ALLOC, NP[n], sizeof(double));
              f[n] = Mem(MEM_ALLOC, NP[n], sizeof(double));
        
              /* Next line */
              
              ENDFNewLine(FUNCTION_NAME, line, fp);
              
              /* Check number of points */

              if (NP[n] != ENDFColI(1,line))
                Die(FUNCTION_NAME, "Mismatch in number of points");

              /* Read interpolation scheme */

              INTT[n] = ENDFColI(2, line);

              /* Next line */

              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Loop over data */

              i = 1;
              for (m = 0; m < NP[n]; m++)
                {
                  /* Read energy and fraction */

                  E[n][m] = ENDFColF(i++, line)/1E+6;
                  f[n][m] = ENDFColF(i++, line);
                  
                  /* Check values (datassa on virheitä, toi f leikataan */
                  /* ykköseksi tuolla alempana) */

                  CheckValue(FUNCTION_NAME, "E", "", E[n][m], 1E-12, 1000.0);
                  CheckValue(FUNCTION_NAME, "f", "", f[n][m], 0.0, 1.5);

                  if ((i > 6) && (m < NP[n] - 1))
                    {
                      ENDFNewLine(FUNCTION_NAME, line, fp);
                      i = 1;
                    }
                }
            }
          
          /*******************************************************************/

          /***** Store and process *******************************************/

          /* Find ground state ratio */

          for (n = 0; n < NS; n++)
            if (LFS[n] == 0)
              break;

          /* Check */

          if (n == NS)
            {
              /* No ground state found, check if single mode */

              if (NS != 1)
                Die(FUNCTION_NAME, 
                    "Multiple excited states without ground state");

              /* Calculate ground state factors */

              for (n = 0; n < NP[0]; n++)
                f[0][n] = 1.0 - f[0][n];

              /* Set index */
              
              n = 0;
            }

          /* Count points that differ from unity */

          i = 0;
          for (m = 0; m < NP[n]; m++)
            if (f[n][m] < 1.0)
              i++;

          /* Do energy cut-off */

          if (E[n][0] > RDB[DATA_NEUTRON_EMAX])
            i = 0;

          /* Check */

          if ((i > 0) && (NP[n] > 1))
            {
              /* Allocate memory for structure */
              
              loc0 = NewItem(DATA_PTR_BRA_LIST, BRA_LIST_BLOCK_SIZE);

              /* Put data */
              
              WDB[loc0 + BRA_LIST_ZAI] = (double)ZAI;
              WDB[loc0 + BRA_LIST_MT] = (double)mt;
              WDB[loc0 + BRA_LIST_FILE] = 9.0;
              WDB[loc0 + BRA_LIST_NP] = (double)NP[n];
              WDB[loc0 + BRA_LIST_INTT] = (double)INTT[n];
              WDB[loc0 + BRA_LIST_FIX_FRAC] = -1.0;
              WDB[loc0 + BRA_LIST_PTR_FNAME] = RDB[pta];

              /* Make energy grid */

              ptr = MakeEnergyGrid(NP[n], 0, 0, -1, E[n], EG_INTERP_MODE_LIN);
              WDB[loc0 + BRA_LIST_PTR_E] = (double)ptr;

              /* Allocate memory for ground state fractions */

              ptr = ReallocMem(DATA_ARRAY, NP[n]);
              WDB[loc0 + BRA_LIST_PTR_F] = (double)ptr;

              /* Copy data */

              for (m = 0; m < NP[n]; m++)
                {
                  /* Truncate to 1.0 */

                  if (f[n][m] > 1.0)
                    WDB[ptr++] = 1.0;
                  else
                    WDB[ptr++] = f[n][m];
                }
            }

          /* Free memory */
          
          for (n = 0; n < NS; n++)
            {
              Mem(MEM_FREE, E[n]);
              Mem(MEM_FREE, f[n]);
            }

          /*******************************************************************/
        }

      /* Close file */
      
      fclose(fp);
    
      /***********************************************************************/

      /***** Read cross sections from file 10 ********************************/

      /* Open file for reading */
  
      fp = OpenDataFile(pta, "isomeric branching ratio data");
        
      /* Read data */

      while (1 != 2)
        {
          /*******************************************************************/

          /***** Read data from file *****************************************/

          /* Loop until line 1 of file 1 */
      
          do 
            eof = fgets(line, 82, fp);
          while (((line[71] != '1') || (line[78] != ' ') || (line[79] != '1')) 
                 && (eof != NULL));
        
          /* Check EOF */
          
          if (eof == NULL)
            break;

          /* Get material number */

          MAT = ENDFMat(line);

          /* Read ZA */
          
          ZA = (long)rint(ENDFColF(1, line));

          /* Next line */

          ENDFNewLine(FUNCTION_NAME, line, fp);

          /* Read LIS */

          LIS = ENDFColI(4, line);
          
          /* Convert all isomeric states to stae 1 */

          if (LIS > 0)
            LIS = 1;

          /* Calculate ZAI */
          
          ZAI = 10*ZA + LIS;

          /* Skip few lines */

          ENDFNewLine(FUNCTION_NAME, line, fp);
          ENDFNewLine(FUNCTION_NAME, line, fp);

          /* Get number of comment lines and records */ 

          NWD = ENDFColI(5, line);
          NXC = ENDFColI(6, line);

          /* Skip comment lines */

          for (n = 0; n < NWD; n++)
            ENDFNewLine(FUNCTION_NAME, line, fp);

          /* Reset count */

          NS = 0;

          /* Read reaction data */

          for (n = 0; n < NXC; n++)
            {
              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Check file type and store MT (cannot handle mt 5) */

              if ((ENDFColI(3, line) == 10) && (ENDFColI(4, line) != 5))
                  MT[NS++] = ENDFColI(4, line);

              /* Check if array sizes are sufficient */

              if (NS > 20)
                Die(FUNCTION_NAME, "Increase array sizes");
            }

          /* Read data from file 3 */

          for (n = 0; n < NS; n++)
            {
              /* Remember position */

              fgetpos(fp, &pos);

              /* Loop until correct mt in file */
      
              do 
                eof = fgets(line, 82, fp);
              while (((ENDFMF(line) != 3) || (ENDFMT(line) != MT[n])) && 
                     (eof != NULL) && (MAT == ENDFMat(line)));
          
              /* Check EOF */
          
              if ((eof == NULL) || (MAT != ENDFMat(line)))
                {
                  /* File may not exist if data is not processed with reconr */

                  Note(0, "Nuclide %ld in \"%s\" has no MF 3 for MT %ld", 
                        ZAI, GetText(pta), MT[n]);

                  /* Reset number of points and pointers */

                  NP3[n] = -1;
                  E3[n] = NULL;
                  xs3[n] = NULL;

                  /* Rewind */

                  fsetpos(fp, &pos);
                  
                  /* Cycle loop */

                  continue;
                }

              /* Next line */

              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Check NR and read NP */
              
              if (ENDFColI(5, line) != 1)
                Die(FUNCTION_NAME, "NR = %ld", ENDFColI(5, line));

              NP3[n] = ENDFColI(6, line);

              /* Allocate memory for data */
              
              E3[n] = Mem(MEM_ALLOC, NP3[n], sizeof(double));
              xs3[n] = Mem(MEM_ALLOC, NP3[n], sizeof(double));

              /* Next line */
              
              ENDFNewLine(FUNCTION_NAME, line, fp);
              
              /* Check number of points */

              if (NP3[n] != ENDFColI(1,line))
                Die(FUNCTION_NAME, "Mismatch in number of points");

              /* Check interpolation scheme */

              if (ENDFColI(2, line) != 2)
                Die(FUNCTION_NAME, "INTT = %ld\n", ENDFColI(2,line));

              /* Next line */

              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Loop over data */
              
              i = 1;
              for (m = 0; m < NP3[n]; m++)
                {
                  /* Read energy and fraction */

                  E3[n][m] = ENDFColF(i++, line)/1E+6;
                  xs3[n][m] = ENDFColF(i++, line);
                  
                  /* Check values */

                  CheckValue(FUNCTION_NAME, "E", "", E3[n][m], 1E-12, 1000.0);
                  CheckValue(FUNCTION_NAME, "xs", "", xs3[n][m], 0.0, MAX_XS);

                  if ((i > 6) && (m < NP3[n] - 1))
                    {
                      ENDFNewLine(FUNCTION_NAME, line, fp);
                      i = 1;
                    }
                }
            }
        
          /* Read data from file 10 */

          for (n = 0; n < NS; n++)
            {
              /* Check if file 3 was available */

              if (NP3[n] < 0)
                {
                  /* Reset number of points and pointers */

                  NP10[n] = -1;
                  E10[n] = NULL;
                  xs10[n] = NULL;

                  /* Cycle loop */

                  continue;
                }

              /* Loop until correct mt in file */
      
              do 
                eof = fgets(line, 82, fp);
              while (((ENDFMF(line) != 10) || (ENDFMT(line) != MT[n])) && 
                     (eof != NULL));

              /* Check EOF */
          
              if (eof == NULL)
                Die(FUNCTION_NAME, 
                    "Nuclide %ld in \"%s\" has no MF 10 for MT %ld", 
                    ZAI, GetText(pta), MT[n]);

              /* Get ZA and LIS */

              ZA = (long)rint(ENDFColF(1, line));
              LIS = ENDFColI(3, line);

              /* Compare ZAI */

              if (ZAI != 10*ZA + LIS)
                Die(FUNCTION_NAME, "Mismatch in ZAI");

              /* Check NS */

              if ((m = ENDFColI(5, line)) > 4)
                Die(FUNCTION_NAME, "NS = %ld", ENDFColI(5, line));

              /* Next line */

              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Check LFS and read LFS and NP */

              if (ENDFColI(5, line) != 1)
                Die(FUNCTION_NAME, "NR = %ld", ENDFColI(5, line));
              
              LFS[n] = ENDFColI(4, line);
              NP10[n] = ENDFColI(6, line);

              /* Check that data is for ground state if multiple states */
              /* are given */

              if ((m > 1) && (LFS[n] != 0))
                Die(FUNCTION_NAME, "Multiple states without ground state");

              /* Allocate memory for data */
              
              E10[n] = Mem(MEM_ALLOC, NP10[n], sizeof(double));
              xs10[n] = Mem(MEM_ALLOC, NP10[n], sizeof(double));

              /* Next line */
              
              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Check number of points */

              if (NP10[n] != ENDFColI(1,line))
                Die(FUNCTION_NAME, "Mismatch in number of points");

              /* Check interpolation scheme */

              if (ENDFColI(2, line) != 2)
                Die(FUNCTION_NAME, "INTT = %ld\n", ENDFColI(2,line));

              /* Next line */

              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Loop over data */
              
              i = 1;
              for (m = 0; m < NP10[n]; m++)
                {
                  /* Read energy and fraction */

                  E10[n][m] = ENDFColF(i++, line)/1E+6;
                  xs10[n][m] = ENDFColF(i++, line);
                  
                  /* Check values */

                  CheckValue(FUNCTION_NAME, "E", "", E10[n][m], 1E-12, 1000.0);
                  CheckValue(FUNCTION_NAME, "xs", "", xs10[n][m], 0.0, MAX_XS);

                  if ((i > 6) && (m < NP10[n] - 1))
                    {
                      ENDFNewLine(FUNCTION_NAME, line, fp);
                      i = 1;
                    }
                }
            }

          /*******************************************************************/

          /***** Store and process *******************************************/

          /* Loop over reactions */

          for (n = 0; n < NS; n++)
            {
              /* Check if file 3 was available */

              if (NP3[n] < 0)
                continue;
                
              /* Do energy cut-off */

              if ((E3[n][0] > RDB[DATA_NEUTRON_EMAX]) || 
                  (E10[n][0] > RDB[DATA_NEUTRON_EMAX]))
                continue;
                
              /* Reset union grid pointer and number of points */

              Eg = NULL;
              np = 0;

              /* Unionize grids */

              Eg = AddPts(Eg, &np, E3[n], NP3[n]);
              Eg = AddPts(Eg, &np, E10[n], NP10[n]);

              /* Allocate memory for temporary cross section array */

              xs = Mem(MEM_ALLOC, np, sizeof(double));

              /* Allocate memory for structure */
              
              loc0 = NewItem(DATA_PTR_BRA_LIST, BRA_LIST_BLOCK_SIZE);

              /* Put data */
              
              WDB[loc0 + BRA_LIST_ZAI] = (double)ZAI;
              WDB[loc0 + BRA_LIST_MT] = (double)MT[n];
              WDB[loc0 + BRA_LIST_FILE] = 10.0;
              WDB[loc0 + BRA_LIST_NP] = (double)np;
              WDB[loc0 + BRA_LIST_INTT] = 2.0;
              WDB[loc0 + BRA_LIST_FIX_FRAC] = -1.0;
              WDB[loc0 + BRA_LIST_PTR_FNAME] = RDB[pta];

              /* Make energy grid */
              
              ptr = MakeEnergyGrid(np, 0, 0, -1, Eg, EG_INTERP_MODE_LIN);
              WDB[loc0 + BRA_LIST_PTR_E] = (double)ptr;

              ptr = ReallocMem(DATA_ARRAY, np);
              WDB[loc0 + BRA_LIST_PTR_F] = (double)ptr;

              /* Reconstruct production cross section (assume lin-lin) */

              InterpolateData(Eg, xs, np, E10[n], xs10[n], NP10[n], 0, 
                              NULL, NULL, NO);

              /* Store */

              for (m = 0; m < np; m++)
                WDB[ptr + m] = xs[m];

              /* Reconstruct total reaction cross section (assume lin-lin) */

              InterpolateData(Eg, xs, np, E3[n], xs3[n], NP3[n], 0, 
                              NULL, NULL, NO);
           
              /* Divide */

              for (m = 0; m < np; m++)
                {
                  if (xs[m] > 0.0)
                    WDB[ptr + m] = RDB[ptr + m]/xs[m];
                  else if ((LFS[n] > 0) && (RDB[ptr + m] != 0.0))
                    Die(FUNCTION_NAME, "Error in cross section");
                }

              /* Calculate ground state fraction */

              for (m = 0; m < np; m++)
                {
                  /* Truncate */
                  
                  if (RDB[ptr + m] > 1.0)
                    WDB[ptr + m] = 1.0;

                  /* Check which value was read */

                  if (LFS[n] == 0)
                    WDB[ptr + m] = RDB[ptr + m];
                  else
                    WDB[ptr + m] = 1.0 - RDB[ptr + m];
                }

              /* Free memory */

              Mem(MEM_FREE, Eg);
              Mem(MEM_FREE, xs);              
            }

          /* Free memory */
          
          for (n = 0; n < NS; n++)
            if (NP3[n] > 0)
              {
                Mem(MEM_FREE, E3[n]);
                Mem(MEM_FREE, xs3[n]);
                Mem(MEM_FREE, E10[n]);
                Mem(MEM_FREE, xs10[n]);
              }
        
          /*******************************************************************/
        }

      /* Close file */
      
      fclose(fp);

      /***********************************************************************/

      /* Next file */

      pta++;
    }

  /* Check */

  if ((loc0 = (long)RDB[DATA_PTR_BRA_LIST]) < VALID_PTR)
    return;
    
  /* Sort list */

  SortList(loc0, BRA_LIST_MT, SORT_MODE_ASCEND);
  SortList(loc0, BRA_LIST_ZAI, SORT_MODE_ASCEND);

  /* Remove duplicates */

  do
    {
      /* Reset count */

      i = 0;

      /* Loop over list */
      
      loc0 = (long)RDB[DATA_PTR_BRA_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Pointer to next */
          
          if ((ptr = NextItem(loc0)) < VALID_PTR)
            break;
          
          /* Compare ZAI and MT */
          
          if (((long)RDB[loc0 + BRA_LIST_ZAI] == 
               (long)RDB[ptr + BRA_LIST_ZAI]) &&
              ((long)RDB[loc0 + BRA_LIST_MT] == 
               (long)RDB[ptr + BRA_LIST_MT]))
            {
              /* Update counter */

              i++;

              /* Pointer to next */
              
              loc0 = NextItem(ptr);
              
              /* Remove duplicate */

              RemoveItem(ptr);
            }
          else
            loc0 = ptr;
        }
    }
  while (i > 0);

  /* Print */

  fprintf(outp, "Energy-dependent isomeric branching ratios:\n\n");

  loc0 = (long)RDB[DATA_PTR_BRA_LIST];
  while (loc0 > VALID_PTR)
    {
      fprintf(outp, " - Nuclide %s mt %ld : %s\n", 
              ZAItoIso((long)RDB[loc0 + BRA_LIST_ZAI], 1), 
              (long)RDB[loc0 + BRA_LIST_MT],
              ReactionMT((long)RDB[loc0 + BRA_LIST_MT], YES));              
      /* Next */

      loc0 = NextItem(loc0);
    }

  fprintf(outp, "\n");
}
/*****************************************************************************/
