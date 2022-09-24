/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readfissionyields.c                            */
/*                                                                           */
/* Created:       2010/09/12 (JLe)                                           */
/* Last modified: 2020/06/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reads fission yield data from ENDF format files              */
/*                                                                           */
/* Comments: - Converted from readfissyields.c in Serpent 1.1.12             */
/*           - Lisää xenon ja samarium yieldit nuklideihin                   */
/*           - Tää ei varmaan toimi jos on useampia file nameja              */
/*           - Randomization affects only independent yields                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadFissionYields:"

/*****************************************************************************/

void ReadFissionYields()
{
  double E, yield, sum, std, f;
  long n, m, i, loc0, loc1, ZA, NWD, I, LE, NFP, ptr, ZAI, Z, A, yld, pta;
  char line[256], *eof;
  FILE *fp;

  /***************************************************************************/

  /***** Neutron-induced fission  ********************************************/

  if ((pta = (long)RDB[DATA_PTR_NFYDATA_FNAME_LIST]) > 0)
    {
      /* Print */

      fprintf(outp, "Reading neutron-induced fission yields...\n");

      /* Loop over files */

      while ((long)RDB[pta] > 0)
        {
          /*******************************************************************/

          /***** Independent yields ******************************************/

          /* Test format */

          WDB[DATA_DUMMY] = RDB[pta];
          TestDOSFile(GetText(DATA_DUMMY));

          /* Open file for reading */

          fp = OpenDataFile(pta, "NFY data");

          /* Loop over endf file */

          do
            {
              /* Loop until in comment section */

              do
                eof = fgets(line, 82, fp);
              while ((atoi(&line[71]) != 1451) && (eof != NULL));

              /* Check EOF */

              if (eof == NULL)
                break;

              /* Read ZA */

              ZA = (long)rint(ENDFColF(1, line));

              /* Check that value is reasonable, if not, this is not */
              /* a fission yield file. */

              if (ZA < 89000)
                Error(0, "No fission yield data available in file \"%s\"",
                      GetText(pta));

              /* Check value */

              CheckValue(FUNCTION_NAME, "ZA", "", ZA, 89000, 120000);

              /* Read isomeric state number */

              ENDFNewLine(FUNCTION_NAME, line, fp);
              I = ENDFColI(4, line);

              /* Check value */

              CheckValue(FUNCTION_NAME, "I", "", I, 0, 1);

              /* Convert to ZAI (this is the parent, ZA and I are re-used */
              /* in the next loop). */

              ZAI = 10*ZA + I;

              /* Skip line */

              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Read header size */

              ENDFNewLine(FUNCTION_NAME, line, fp);
              NWD = ENDFColI(5, line);

              /* Skip header */

              fseek(fp, 81*(NWD + 1), SEEK_CUR);

              /* Skip remaining comment block */

              do
                ENDFNewLine(FUNCTION_NAME, line, fp);
              while (atoi(&line[71]) == 1451);

              /* Loop until fission product yield data. */

              do
                eof = fgets(line, 82, fp);
              while ((atoi(&line[71]) != 8454) && (eof != NULL));

              /* Check EOF */

              if (eof == NULL)
                Die(FUNCTION_NAME, "Error reading fission yield file \"%s\"",
                    GetText(pta));

              /* Get number of incident energies */

              LE = ENDFColI(3, line);

              /* Check value */

              CheckValue(FUNCTION_NAME, "LE", "", LE, 0, 4);

              /* Loop over energies */

              for (n = 0; n < LE; n++)
                {
                  /***********************************************************/

                  /***** Store data ******************************************/

                  /* Allocate memory for data */

                  yld = ReallocMem(ACE_ARRAY, FISSION_YIELD_BLOCK_SIZE);

                  /* Set null pointer to next */

                  ACE[yld + FISSION_YIELD_PTR_NEXT] = NULLPTR;

                  /* Check if previous exists (ei voi käyttää VALID_PTR) */

                  if ((ptr = (long)RDB[DATA_PTR_ACE_NFY_DATA]) < 0)
                    {
                      /* First definition, set pointer */

                      WDB[DATA_PTR_ACE_NFY_DATA] = (double)yld;
                    }
                  else
                    {
                      /* Find last block  (tohon ei VALID_PTR) */

                      while ((long)ACE[ptr + FISSION_YIELD_PTR_NEXT] > 0)
                        ptr = (long)ACE[ptr + FISSION_YIELD_PTR_NEXT];

                      /* Set pointer to new */

                      ACE[ptr + FISSION_YIELD_PTR_NEXT] = (double)yld;
                    }

                  /* Put parent ZAI */

                  ACE[yld + FISSION_YIELD_PARENT_ZAI] = (double)ZAI;

                  /* Put index */

                  ACE[yld + FISSION_YIELD_IDX] = (double)(n + 1);

                  /* Get energy */

                  ENDFNewLine(FUNCTION_NAME, line, fp);
                  E = 1E-6*ENDFColF(1, line);

                  /* Check value */

                  CheckValue(FUNCTION_NAME, "E", "(independent)", E, 0.0, 20.0);

                  /* Put value */

                  ACE[yld + FISSION_YIELD_E] = E;

                  /* Get number of fission products */

                  NFP = ENDFColI(6, line);

                  /* Check value */

                  CheckValue(FUNCTION_NAME, "NFP", "", NFP, 700,
                             MAX_FP_NUCLIDES);

                  /* Additional check for upper limit */

                  if (NFP > MAX_FP_NUCLIDES)
                    Die(FUNCTION_NAME,
                        "Number of FP nuclides exceeds maximum");

                  /* Put value */

                  ACE[yld + FISSION_YIELD_NFP] = (double)NFP;

                  /* Allocate memory for data */

                  loc0 = ReallocMem(ACE_ARRAY, NFP + 1);

                  /* Set pointer */

                  ACE[yld + FISSION_YIELD_PTR_DISTR] = (double)loc0;

                  /* Read first line */

                  ENDFNewLine(FUNCTION_NAME, line, fp);

                  /* Loop over data */

                  i = 1;
                  sum = 0.0;

                  for (m = 0; m < NFP; m++)
                    {
                      /* Read isotope ZA */

                      ZA = (long)ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "ZA", "", ZA, 1001, 80000);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Separate Z and A (only needed for bug check) */

                      Z = (long)((double)ZA/1000.0);
                      A = ZA - 1000*Z;

                      /* Check values */

                      CheckValue(FUNCTION_NAME, "Z", "", Z, 1, 80);
                      CheckValue(FUNCTION_NAME, "A", "", A, 1, 210);

                      /* Read isomeric state */

                      I = (long)ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "I", "", I, 0, 2);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Read yield */

                      yield = ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "yield", "",yield, 0.0, 1.0);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Check if value is perturbed */

                      if ((long)RDB[DATA_BURN_RANDOMIZE_FY] == YES)
                        {
                          /* Read standard deviation */

                          std = ENDFColF(i++, line);

                          /* Randomize */

                          if ((f = RandNorm(yield, std, 0)) > 0.0)
                            yield = f;
                        }
                      else
                        {
                          /* Skip uncertainty */

                          i++;
                        }

                      if ((i > 6) && (m < NFP - 1))
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Allocate memory */

                      loc1 = ReallocMem(ACE_ARRAY, FY_BLOCK_SIZE);

                      /* Put pointer */

                      ACE[loc0++] = (double)loc1;

                      /* Put values */

                      ACE[loc1 + FY_TGT_ZAI] = (double)(10*ZA + I);
                      ACE[loc1 + FY_INDEPENDENT_FRAC] = yield;

                      /* Add to sum */

                      sum = sum + yield;
                    }

                  /* Put null pointer */

                  ACE[loc0] = NULLPTR;

                  /* Check sum (may differ from 2 due to cut-off and ternary */
                  /* fission) */

                  CheckValue(FUNCTION_NAME, "sum", "", sum, 1.5, 2.5);

                  /* Set warning flag */

                  if ((sum < 1.9999) || (sum > 2.01))
                    WDB[DATA_WARN_NFY_SUM] = RDB[DATA_WARN_NFY_SUM] + 1.0;
                }
            }
          while (1 != 2);

          /* Close file */

          fclose(fp);

          /*******************************************************************/

          /***** Cumulative yields *******************************************/

          fp = OpenDataFile(pta, "NFY data");

          /* Pointer to data (tässä oletetaan että noi independent yieldit on */
          /* luettu ja cumulative yieldit on samassa järjestyksessä). */

          yld = (long)RDB[DATA_PTR_ACE_NFY_DATA];

          /* Loop over endf file */

          do
            {
              /* Loop until in comment section */

              do
                eof = fgets(line, 82, fp);
              while ((atoi(&line[71]) != 1451) && (eof != NULL));

              /* Check EOF */

              if (eof == NULL)
                break;

              /* Read ZA */

              ZA = (long)rint(ENDFColF(1, line));

              /* Check value */

              CheckValue(FUNCTION_NAME, "ZA", "", ZA, 89000, 120000);

              /* Read isomeric state number */

              ENDFNewLine(FUNCTION_NAME, line, fp);
              I = ENDFColI(4, line);

              /* Check value */

              CheckValue(FUNCTION_NAME, "I", "", I, 0, 1);

              /* Convert to ZAI (this is the parent, ZA and I are re-used */
              /* in the next loop). */

              ZAI = 10*ZA + I;

              /* Skip line */

              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Read header size */

              ENDFNewLine(FUNCTION_NAME, line, fp);
              NWD = ENDFColI(5, line);

              /* Skip header */

              fseek(fp, 81*(NWD + 1), SEEK_CUR);

              /* Skip remaining comment block */

              do
                ENDFNewLine(FUNCTION_NAME, line, fp);
              while (atoi(&line[71]) == 1451);

              /* Loop until fission product yield data. */

              do
                eof = fgets(line, 82, fp);
              while ((atoi(&line[71]) != 8459) && (eof != NULL));

              /* Check EOF */

              if (eof == NULL)
                Die(FUNCTION_NAME, "Error reading fission yield file \"%s\"",
                    GetText(pta));

              /* Get number of incident energies */

              LE = ENDFColI(3, line);

              /* Check value */

              CheckValue(FUNCTION_NAME, "LE", "", LE, 0, 4);

              /* Loop over energies */

              for (n = 0; n < LE; n++)
                {
                  /***********************************************************/

                  /***** Store data ******************************************/

                  /* Compare ZAI */

                  if ((long)ACE[yld + FISSION_YIELD_PARENT_ZAI] != ZAI)
                    Die(FUNCTION_NAME, "Mismatch in parent ZAI");

                  /* Get energy */

                  ENDFNewLine(FUNCTION_NAME, line, fp);
                  E = 1E-6*ENDFColF(1, line);

                  /* Compare energy */

                  if (ACE[yld + FISSION_YIELD_E] != E)
                    Die(FUNCTION_NAME, "Mismatch in energy");

                  /* Get number of fission products */

                  NFP = ENDFColI(6, line);

                  /* compare value */

                  if ((long)ACE[yld + FISSION_YIELD_NFP] != NFP)
                    Die(FUNCTION_NAME, "Mismatch in NFP");

                  /* Set pointer to data */

                  loc0 = (long)ACE[yld + FISSION_YIELD_PTR_DISTR];

                  /* Read first line */

                  ENDFNewLine(FUNCTION_NAME, line, fp);

                  /* Loop over data */

                  i = 1;
                  sum = 0.0;

                  for (m = 0; m < NFP; m++)
                    {
                      /* Read isotope ZA */

                      ZA = (long)ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "ZA", "", ZA, 1001, 80000);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Read isomeric state */

                      I = (long)ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "I", "", I, 0, 2);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Read yield */

                      yield = ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "yield", "",yield, 0.0, 1.0);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Skip uncertainty */

                      i++;

                      if ((i > 6) && (m < NFP - 1))
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Get pointer */

                      loc1 = (long)ACE[loc0++];

                      /* Compare ZAI */

                      if ((long)ACE[loc1 + FY_TGT_ZAI] != 10*ZA + I)
                        Die(FUNCTION_NAME, "Mismatch in product ZAI");

                      /* Compare yield */

                      if ((yield > 0.0) &&
                          ((long)RDB[DATA_BURN_RANDOMIZE_FY] == NO) &&
                          (ACE[loc1 + FY_INDEPENDENT_FRAC]/yield- 1.0 > 1E-4))
                        {
                          /* Check major error */

                          if (ACE[loc1 + FY_INDEPENDENT_FRAC]/yield - 1.0 >
                              0.005)
                            Die (FUNCTION_NAME,
                                 "Independent yield > cumulative (%ld %ld %E %E %E)",
                                 ZAI, 10*ZA + I,
                                 yield, ACE[loc1 + FY_INDEPENDENT_FRAC], E);

#ifdef DEBUG
                          /* Print warning */

                          Warn(FUNCTION_NAME,
                               "Independent yield > cumulative (%ld %ld %E %E %E)",
                               ZAI, 10*ZA + I,
                               yield, ACE[loc1 + FY_INDEPENDENT_FRAC], E);
#endif
                        }

                      /* Use independent fraction if cumulative is zero */

                      if (yield == 0.0)
                        yield = ACE[loc1 + FY_INDEPENDENT_FRAC];

                      /* Put yield */

                      ACE[loc1 + FY_CUMULATIVE_FRAC] = yield;
                    }

                  /* Next yield */

                  yld = (long)ACE[yld + FISSION_YIELD_PTR_NEXT];
                }
            }
          while (1 != 2);

          /* Close file */

          fclose(fp);

          /*******************************************************************/

          /* Next file */

          pta++;
        }

      /***********************************************************************/

      fprintf(outp, "OK.\n\n");
    }

  /***************************************************************************/

  /***** Spontaneous fission  ************************************************/

  if ((pta = (long)RDB[DATA_PTR_SFYDATA_FNAME_LIST]) > 0)
    {
      /* Print */

      fprintf(outp, "Reading spontaneous fission yields...\n");

      /* Loop over files */

      while ((long)RDB[pta] > 0)
        {
          /*******************************************************************/

          /***** Independent yields ******************************************/

          /* Open file for reading */

          fp = OpenDataFile(pta, "NFY data");

          /* Loop over endf file */

          do
            {
              /* Loop until in comment section */

              do
                eof = fgets(line, 82, fp);
              while ((atoi(&line[71]) != 1451) && (eof != NULL));

              /* Check EOF */

              if (eof == NULL)
                break;

              /* Read ZA */

              ZA = (long)rint(ENDFColF(1, line));

              /* Check that value is reasonable, if not, this is not */
              /* a fission yield file. */

              if (ZA < 89000)
                Error(0, "No fission yield data available in file \"%s\"",
                      GetText(pta));

              /* Check value */

              CheckValue(FUNCTION_NAME, "ZA", "", ZA, 89000, 120000);

              /* Read isomeric state number */

              ENDFNewLine(FUNCTION_NAME, line, fp);
              I = ENDFColI(4, line);

              /* Check value */

              CheckValue(FUNCTION_NAME, "I", "", I, 0, 1);

              /* Convert to ZAI (this is the parent, ZA and I are re-used */
              /* in the next loop). */

              ZAI = 10*ZA + I;

              /* Skip line */

              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Read header size */

              ENDFNewLine(FUNCTION_NAME, line, fp);
              NWD = ENDFColI(5, line);

              /* Skip header */

              fseek(fp, 81*(NWD + 1), SEEK_CUR);

              /* Skip remaining comment block */

              do
                ENDFNewLine(FUNCTION_NAME, line, fp);
              while (atoi(&line[71]) == 1451);

              /* Loop until fission product yield data. */

              do
                eof = fgets(line, 82, fp);
              while ((atoi(&line[71]) != 8454) && (eof != NULL));

              /* Check EOF */

              if (eof == NULL)
                Die(FUNCTION_NAME, "Error reading fission yield file \"%s\"",
                    GetText(pta));

              /* Get number of incident energies */

              LE = ENDFColI(3, line);

              /* Check value */

              CheckValue(FUNCTION_NAME, "LE", "", LE, 1, 1);

              /* Loop over energies */

              for (n = 0; n < LE; n++)
                {
                  /***********************************************************/

                  /***** Store data ******************************************/

                  /* Allocate memory for data */

                  yld = ReallocMem(ACE_ARRAY, FISSION_YIELD_BLOCK_SIZE);

                  /* Set null pointer to next */

                  ACE[yld + FISSION_YIELD_PTR_NEXT] = NULLPTR;

                  /* Check if previous exists (ei voi käyttää VALID_PTR) */

                  if ((ptr = (long)RDB[DATA_PTR_ACE_SFY_DATA]) < 0)
                    {
                      /* First definition, set pointer */

                      WDB[DATA_PTR_ACE_SFY_DATA] = (double)yld;
                    }
                  else
                    {
                      /* Find last block  (tohon ei VALID_PTR) */

                      while ((long)ACE[ptr + FISSION_YIELD_PTR_NEXT] > 0)
                        ptr = (long)ACE[ptr + FISSION_YIELD_PTR_NEXT];

                      /* Set pointer to new */

                      ACE[ptr + FISSION_YIELD_PTR_NEXT] = (double)yld;
                    }

                  /* Put parent ZAI */

                  ACE[yld + FISSION_YIELD_PARENT_ZAI] = (double)ZAI;

                  /* Get energy */

                  ENDFNewLine(FUNCTION_NAME, line, fp);
                  E = 1E-6*ENDFColF(1, line);

                  /* Check value */

                  CheckValue(FUNCTION_NAME, "E", "(independent)", E, 0.0, 20.0);

                  /* Put value */

                  ACE[yld + FISSION_YIELD_E] = E;

                  /* Get number of fission products */

                  NFP = ENDFColI(6, line);

                  /* Check value */

                  CheckValue(FUNCTION_NAME, "NFP", "", NFP, 700,
                             MAX_FP_NUCLIDES);

                  /* Additional check for upper limit */

                  if (NFP > MAX_FP_NUCLIDES)
                    Die(FUNCTION_NAME,
                        "Number of FP nuclides exceeds maximum");

                  /* Put value */

                  ACE[yld + FISSION_YIELD_NFP] = (double)NFP;

                  /* Allocate memory for data */

                  loc0 = ReallocMem(ACE_ARRAY, NFP + 1);

                  /* Set pointer */

                  ACE[yld + FISSION_YIELD_PTR_DISTR] = (double)loc0;

                  /* Read first line */

                  ENDFNewLine(FUNCTION_NAME, line, fp);

                  /* Loop over data */

                  i = 1;
                  sum = 0.0;

                  for (m = 0; m < NFP; m++)
                    {
                      /* Read isotope ZA */

                      ZA = (long)ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "ZA", "", ZA, 1001, 80000);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Separate Z and A (only needed for bug check) */

                      Z = (long)((double)ZA/1000.0);
                      A = ZA - 1000*Z;

                      /* Check values */

                      if ((Z < 1) || (Z > 80))
                        Die(FUNCTION_NAME, "Error in Z");

                      if ((A < 1) || (A > 210))
                        Die(FUNCTION_NAME, "Error in A");

                      /* Read isomeric state */

                      I = (long)ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "I", "", I, 0, 2);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Read yield */

                      yield = ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "yield", "",yield, 0.0, 1.0);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Check if value is perturbed */

                      if ((long)RDB[DATA_BURN_RANDOMIZE_FY] == YES)
                        {
                          /* Read standard deviation */

                          std = ENDFColF(i++, line);

                          /* Randomize */

                          if ((f = RandNorm(yield, std, 0)) > 0.0)
                            yield = f;
                        }
                      else
                        {
                          /* Skip uncertainty */

                          i++;
                        }

                      if ((i > 6) && (m < NFP - 1))
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Allocate memory */

                      loc1 = ReallocMem(ACE_ARRAY, FY_BLOCK_SIZE);

                      /* Put pointer */

                      ACE[loc0++] = (double)loc1;

                      /* Put values */

                      ACE[loc1 + FY_TGT_ZAI] = (double)(10*ZA + I);
                      ACE[loc1 + FY_INDEPENDENT_FRAC] = yield;

                      /* Add to sum */

                      sum = sum + yield;
                    }

                  /* Put null pointer */

                  ACE[loc0] = NULLPTR;

                  /* Check sum (may differ from 2 due to cut-off and ternary */
                  /* fission) */

                  CheckValue(FUNCTION_NAME, "sum", "", sum, 1.5, 2.5);

                  /* Set warning flag */

                  if ((sum < 1.9999) || (sum > 2.01))
                    WDB[DATA_WARN_SFY_SUM] = RDB[DATA_WARN_SFY_SUM] + 1.0;
                }
            }
          while (1 != 2);

          /* Close file */

          fclose(fp);

          /*******************************************************************/

          /***** Cumulative yields *******************************************/

          fp = OpenDataFile(pta, "SFY data");

          /* Pointer to data (tässä oletetaan että noi independent yieldit */
          /* on luettu ja cumulative yieldit on samassa järjestyksessä). */

          yld = (long)RDB[DATA_PTR_ACE_SFY_DATA];

          /* Loop over endf file */

          do
            {
              /* Loop until in comment section */

              do
                eof = fgets(line, 82, fp);
              while ((atoi(&line[71]) != 1451) && (eof != NULL));

              /* Check EOF */

              if (eof == NULL)
                break;

              /* Read ZA */

              ZA = (long)rint(ENDFColF(1, line));

              /* Check value */

              CheckValue(FUNCTION_NAME, "ZA", "", ZA, 89000, 120000);

              /* Read isomeric state number */

              ENDFNewLine(FUNCTION_NAME, line, fp);
              I = ENDFColI(4, line);

              /* Check value */

              CheckValue(FUNCTION_NAME, "I", "", I, 0, 1);

              /* Convert to ZAI (this is the parent, ZA and I are re-used */
              /* in the next loop). */

              ZAI = 10*ZA + I;

              /* Skip line */

              ENDFNewLine(FUNCTION_NAME, line, fp);

              /* Read header size */

              ENDFNewLine(FUNCTION_NAME, line, fp);
              NWD = ENDFColI(5, line);

              /* Skip header */

              fseek(fp, 81*(NWD + 1), SEEK_CUR);

              /* Skip remaining comment block */

              do
                ENDFNewLine(FUNCTION_NAME, line, fp);
              while (atoi(&line[71]) == 1451);

              /* Loop until fission product yield data. */

              do
                eof = fgets(line, 82, fp);
              while ((atoi(&line[71]) != 8459) && (eof != NULL));

              /* Check EOF */

              if (eof == NULL)
                Die(FUNCTION_NAME, "Error reading fission yield file \"%s\"",
                    GetText(pta));

              /* Get number of incident energies */

              LE = ENDFColI(3, line);

              /* Check value */

              CheckValue(FUNCTION_NAME, "LE", "", LE, 1, 1);

              /* Loop over energies */

              for (n = 0; n < LE; n++)
                {
                  /***********************************************************/

                  /***** Store data ******************************************/

                  /* Compare ZAI */

                  if ((long)ACE[yld + FISSION_YIELD_PARENT_ZAI] != ZAI)
                    Die(FUNCTION_NAME, "Mismatch in parent ZAI");

                  /* Get energy */

                  ENDFNewLine(FUNCTION_NAME, line, fp);
                  E = 1E-6*ENDFColF(1, line);

                  /* Compare energy */
                  /*
                    if (ACE[yld + FISSION_YIELD_E] != E)
                    Die(FUNCTION_NAME, "Mismatch in energy");
                  */
                  /* Get number of fission products */

                  NFP = ENDFColI(6, line);

                  /* compare value */

                  if ((long)ACE[yld + FISSION_YIELD_NFP] != NFP)
                    Die(FUNCTION_NAME, "Mismatch in NFP");

                  /* Set pointer to data */

                  loc0 = (long)ACE[yld + FISSION_YIELD_PTR_DISTR];

                  /* Read first line */

                  ENDFNewLine(FUNCTION_NAME, line, fp);

                  /* Loop over data */

                  i = 1;
                  sum = 0.0;

                  for (m = 0; m < NFP; m++)
                    {
                      /* Read isotope ZA */

                      ZA = (long)ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "ZA", "", ZA, 1001, 80000);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Read isomeric state */

                      I = (long)ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "I", "", I, 0, 2);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Read yield */

                      yield = ENDFColF(i++, line);

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "yield", "",yield, 0.0, 1.0);

                      if (i > 6)
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Skip uncertainty */

                      i++;

                      if ((i > 6) && (m < NFP - 1))
                        {
                          ENDFNewLine(FUNCTION_NAME, line, fp);
                          i = i - 6;
                        }

                      /* Get pointer */

                      loc1 = (long)ACE[loc0++];

                      /* Compare ZAI */

                      if ((long)ACE[loc1 + FY_TGT_ZAI] != 10*ZA + I)
                        Die(FUNCTION_NAME, "Mismatch in product ZAI");

                      /* Compare yield (ENDF/B-VI.8 and -VII datassa pieniä */
                      /* ylityksiä) */

                      if ((yield > 0.0) &&
                          ((long)RDB[DATA_BURN_RANDOMIZE_FY] == NO) &&
                          (ACE[loc1 + FY_INDEPENDENT_FRAC]/yield- 1.0 > 2E-4))
                        Warn(FUNCTION_NAME,
                             "Independent yield exceeds cumulative (%ld %E %E)",
                             ZAI, yield, ACE[loc1 + FY_INDEPENDENT_FRAC]);

                      /* Use independent fraction if cumulative is zero */

                      if (yield == 0.0)
                        yield = ACE[loc1 + FY_INDEPENDENT_FRAC];

                      /* Put yield */

                      ACE[loc1 + FY_CUMULATIVE_FRAC] = yield;
                    }

                  /* Next yield */

                  yld = (long)ACE[yld + FISSION_YIELD_PTR_NEXT];
                }
            }
          while (1 != 2);

          /* Close file */

          fclose(fp);

          /*******************************************************************/

          /* Next file */

          pta++;
        }

      fprintf(outp, "OK.\n\n");
    }

  /***************************************************************************/
}

/*****************************************************************************/
