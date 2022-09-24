/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readinput.c                                    */
/*                                                                           */
/* Created:       2010/09/21 (JLe)                                           */
/* Last modified: 2020/05/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reads input file.                                            */
/*                                                                           */
/* Comments: - If the argument type for TestParam() is PTYPE_INT, the limits */
/*             must be given in integer (or long) format.                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "surface_types.h"

#define FUNCTION_NAME "ReadInput:"

/*****************************************************************************/

void ReadInput(char *inputfile)
{
  long np, i0, i, j, k, l, n, m, loc0, loc1, loc2, loc3, ptr, type, nx, ny, nz;
  long nr, ns, ne, ni, line, r, b, g, lines, ZAI, MT, idx;
  char *input, word[MAX_STR], str[20000], **params, c;
  char fname[MAX_STR], pname[MAX_STR], tmpstr[MAX_STR];
  double val, sum, ax, ay, az, mem, xmin, xmax, ymin, ymax, zmin, zmax, lims[6];
  double mem0, mem1, r0, r1, *mlims;
  FILE *fp;

  /* Print file name */

  fprintf(outp, "Reading input file \"%s\"...\n", inputfile);

  /* Copy input file name */

  strcpy(fname, inputfile);

  /* Reset j and np to avoid compiler warnings */

  j = 0;
  np = 0;

  /* Read input file */

  input = ReadTextFile(inputfile);

  /* Calculate number of lines */

  n = 0;
  lines = 0;
  while(input[n] != '\0')
    if (input[n++] == '\n')
      lines++;

  /* Reset counters for line number calculation */

  WDB[DATA_LINE_NUM_N0] = 0.0;
  WDB[DATA_LINE_NUM_NL0] = 1.0;

  /* Avoid compiler warning */

  params = NULL;

  /***************************************************************************/

  /***** Loop over words *****************************************************/

  i0 = 0;

  while ((i = NextWord(&input[i0], word)) > 0)
    {
      /* Reset parameter name */

      *pname = '\0';

      /* update pointer */

      i0 = i0 + i;

      /* Get line number for error messages */

      line = GetLineNumber(input, i0);

      /***********************************************************************/

      /***** Read another input file *****************************************/

      if (!strcasecmp(word, "include"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 1, fname);

          /* Check that file exists */

          if ((fp = fopen(params[0], "r")) != NULL)
            fclose(fp);
          else
            {
              /* File not found */

              Error(-1, pname, fname, line,
                    "Input file \"%s\" does not exist", params[0]);
            }

          /* Remember line number */

          n = (long)RDB[DATA_LINE_NUM_N0];
          m = (long)RDB[DATA_LINE_NUM_NL0];

          /* Read input file */

          ReadInput(params[0]);

          /* Restore line number */

          WDB[DATA_LINE_NUM_N0] = (double)n;
          WDB[DATA_LINE_NUM_NL0] = (double)m;
        }

      /***********************************************************************/

      /***** Material definition *********************************************/

      else if (!strcasecmp(word, "mat"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, 4*MAX_ISOTOPES + 8,
                         fname);

          /* Get memory size */

          mem = RDB[DATA_TOTAL_BYTES];

          /* Create new item */

          loc0 = NewItem(DATA_PTR_M0, MATERIAL_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Material name */

          WDB[loc0 + MATERIAL_PTR_NAME] = (double)PutText(params[j++]);

          /* Material density */

          if (!strcmp(params[j], "sum"))
            {
              /* Set value to -infinity to calculate sum from composition */

              WDB[loc0 + MATERIAL_ADENS] = -INFTY;

              j++;
            }
          else
            {
              /* Read value */

              WDB[loc0 + MATERIAL_ADENS] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          -1000.0, 1000.0);
            }

          /* Reset temperatures and volume */

          WDB[loc0 + MATERIAL_DOPPLER_TEMP] = -1.0;
          WDB[loc0 + MATERIAL_TMS_TMIN] = INFTY;
          WDB[loc0 + MATERIAL_TMS_TMAX] = -INFTY;
          WDB[loc0 + MATERIAL_VOLUME_GIVEN] = -INFTY;

          /* Reset sum */

          sum = 0.0;

          /* Loop over parameters */

          while (j < np)
            {
              /* Check parameter */

              if (!strcmp(params[j], "tmp"))
                {
                  /***** Temperature for Doppler-breadening ******************/

                  j++;

                  /* Get temperature */

                  WDB[loc0 + MATERIAL_DOPPLER_TEMP] =
                           TestParam(pname, fname, line, params[j++],
                                     PTYPE_REAL, 0.0, 100000.0);

                  /* Set option */

                  WDB[DATA_USE_DOPPLER_PREPROCESSOR] = (double)YES;

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "tms") || !strcmp(params[j], "ettm"))
                {
                  /***** Temperature for TMS ********************************/

                  j++;

                  /* Get temperature */

                  WDB[loc0 + MATERIAL_TMS_TMIN] =
                    TestParam(pname, fname, line, params[j++],
                              PTYPE_REAL, 0.0, 100000.0);

                  /* Copy to maximum */

                  WDB[loc0 + MATERIAL_TMS_TMAX] =
                    RDB[loc0 + MATERIAL_TMS_TMIN];

                  /* Set mode */

                  WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;
                  WDB[loc0 + MATERIAL_TMS_MODE] = (double)YES;

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "rgb"))
                {
                  /***** Material colour *************************************/

                  j++;

                  /* Get r, b and g */

                  r = (long)TestParam(pname, fname, line, params[j++],
                                      PTYPE_INT, 0, 255);

                  g = (long)TestParam(pname, fname, line, params[j++],
                                      PTYPE_INT, 0, 255);

                  b = (long)TestParam(pname, fname, line, params[j++],
                                      PTYPE_INT, 0, 255);

                  /* Set color */

                  WDB[loc0 + MATERIAL_RGB] = (double)(b + 1000*g + 1000000*r);

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "vol"))
                {
                  /***** Material volume *************************************/

                  j++;

                  /* Get volume */

                  WDB[loc0 + MATERIAL_VOLUME_GIVEN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "fix"))
                {
                  /***** Default library ID and temperature*******************/

                  j++;

                  /* Get default ID and temperature */

                  WDB[loc0 + MATERIAL_DEFAULT_PTR_LIB_ID] =
                    (double)PutText(params[j++]);
                  WDB[loc0 + MATERIAL_DEFAULT_TMP] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "mass"))
                {
                  /***** Material mass ***************************************/

                  j++;

                  /* Get mass */

                  WDB[loc0 + MATERIAL_MASS_GIVEN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "burn"))
                {
                  /***** Burnable material ***********************************/

                  j++;

                  /* Set burn flag */

                  SetOption(loc0 + MATERIAL_OPTIONS, OPT_BURN_MAT);

                  /* Set burn sort flag and materials flag */

                  WDB[loc0 + MATERIAL_BURN_SORT_FLAG] = 1.0;
                  WDB[DATA_BURN_MATERIALS_FLAG] = (double)YES;

                  /* Get number of rings */

                  WDB[loc0 + MATERIAL_BURN_RINGS] =
                    TestParam(pname, fname, line, params[j++],
                              PTYPE_INT, 0, 10000000);

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "moder"))
                {
                  /***** Thermal scattering data *****************************/

                  j++;

                  /* Check number of parameters */

                  if (j > np - 3)
                    Error(loc0, "Invalid number of parameters");

                  /* Create new item (use the same structure as with */
                  /* the therm card) */

                  loc1 = NewItem(loc0 + MATERIAL_PTR_SAB, THERM_BLOCK_SIZE);

                  /* Read name */

                  WDB[loc1 + THERM_PTR_ALIAS] =
                    (double)PutText(params[j++]);

                  /* Read ZA */

                  WDB[loc1 + THERM_ZA] =
                    TestParam(pname, fname, line, params[j++],
                              PTYPE_INT, 1001, 120000);

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "tft"))
                {
                  /***** Minimum and maximum temperatures for TMS* ***********/

                  j++;

                  /* Get minimum temperature */

                  WDB[loc0 + MATERIAL_TMS_TMIN] =
                    TestParam(pname, fname, line, params[j++],
                              PTYPE_REAL, 0.0, 100000.0);

                  /* Get maximum temperature */

                  WDB[loc0 + MATERIAL_TMS_TMAX] =
                    TestParam(pname, fname, line, params[j++],
                              PTYPE_REAL, RDB[loc0 + MATERIAL_TMS_TMIN],
                              100000.0);

                  /* Set mode */

                  WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;
                  WDB[loc0 + MATERIAL_TMS_MODE] = (double)YES;

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "coeft"))
                {
                  /***** Temperature for coefficient calculation *************/

                  /* NOTE: This is for testing only */

                  j++;

                  /* Get temperature */

                  WDB[loc0 + MATERIAL_COEF_TEMP] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, 100000.0);

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "coefs"))
                {
                  /***** S(a,b)-data for coefficient calculation *************/

                  /* NOTE: This is for testing only */

                  j++;

                  /* Get temperature */

                  WDB[loc0 + MATERIAL_COEF_SAB] =
                    (double)PutText(params[j++]);

                  /***********************************************************/
                }
              else
                {
                  /***** Composition *****************************************/

                  /* Create new item */

                  loc1 = NewItem(loc0 + MATERIAL_PTR_COMP,
                                 COMPOSITION_BLOCK_SIZE);

                  /* Check number of parameters */

                  if (j > np - 2)
                    Error(loc0, "Invalid number of parameters");

                  /* Read nuclide name */

                  WDB[loc1 + COMPOSITION_PTR_NUCLIDE] =
                    (double)PutText(params[j++]);

                  /* Read fraction */

                  val = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -1E25, 1E+25);

                  WDB[loc1 + COMPOSITION_ADENS] = val;

                  /* Add to sum */

                  sum = sum + val;

                  /***********************************************************/
                }
            }

          /* Set density if sum */

          if (RDB[loc0 + MATERIAL_ADENS] == -INFTY)
            {
              /* Check sum */

              if (sum > 100.0)
                Error(loc0, "Atomic densities shoud be given in 1/barn*cm");
              else
                WDB[loc0 + MATERIAL_ADENS] = sum;
            }

          /* Put memory size */

          WDB[loc0 + MATERIAL_MEM_SIZE] = RDB[DATA_TOTAL_BYTES] - mem;
         }

      /***********************************************************************/

      /***** Depletion branches **********************************************/

      else if (!strcasecmp(word, "branch"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 1000000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_BRA0, DEP_BRA_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Branch name */

          WDB[loc0 + DEP_BRA_PTR_NAME] = (double)PutText(params[j++]);

          /* Reset normalization */

          WDB[loc0 + DEP_BRA_NORM] = -1.0;

          /* Reset zero poison flags */

          WDB[loc0 + DEP_BRA_ZERO_XE] = (double)NO;
          WDB[loc0 + DEP_BRA_ZERO_SM] = (double)NO;

          /* Avoid compiler warning */

          loc1 = -1;
          n = 0;

          /* Loop over parameters */

          while (j < np)
            {
              /* Check parameter */

              if (!strcmp(params[j], "stp"))
                {
                  /* Change in material state */

                  loc1 = NewItem(loc0 + DEP_BRA_PTR_STP,
                                 DEP_BRA_STP_BLOCK_SIZE);

                  /* Update index */

                  j++;

                  /* Read material name */

                  if (j == np)
                    Error(loc0, "Missing material name");
                  else
                    WDB[loc1 + DEP_BRA_STP_PTR_MAT] =
                      (double)PutText(params[j++]);

                  /* Read density */

                  if (j == np)
                    Error(loc0, "Missing material density");
                  else if (!strcmp(params[j], "sum"))
                    {
                      /* Set value to -infinity to calculate sum from */
                      /* composition */

                      WDB[loc1 + DEP_BRA_STP_DENSITY] = -INFTY;

                      j++;
                    }
                  else
                    WDB[loc1 + DEP_BRA_STP_DENSITY] =
                      TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                -1000.0, 1000.0);

                  /* Read temperature */

                  if (j == np)
                    Error(loc0, "Missing temperature");
                  else
                    WDB[loc1 + DEP_BRA_STP_TEMP] =
                      TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                -1.0, 100000.0);

                  /* Set flag for more parameters */

                  n = 1;
                }
              else if (!strcmp(params[j], "gcu"))
                {
                  /* Update index */

                  j++;

                  /* Change universe for group constant generation */

                  WDB[loc0 + DEP_BRA_PTR_GCU] = 0.0;
                  loc1 = NewItem(loc0 + DEP_BRA_PTR_GCU, GCU_BLOCK_SIZE);

                  /* Get universe */

                  WDB[loc1 + GCU_PTR_UNIV] = (double)PutText(params[j++]);
                }
              else if (!strcmp(params[j], "norm"))
                {
                  /* Update index */

                  j++;

                  /* Adjust normalization */

                  WDB[loc0 + DEP_BRA_NORM] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);
                }
              else if (!strcmp(params[j], "xenon"))
                {
                  /* Update index */

                  j++;

                  /* Adjust normalization */

                  WDB[loc0 + DEP_BRA_ZERO_XE] =
                    (double)(1 - TestParam(pname, fname, line,
                                           params[j++], PTYPE_LOGICAL));
                }
              else if (!strcmp(params[j], "samarium"))
                {
                  /* Update index */

                  j++;

                  /* Adjust normalization */

                  WDB[loc0 + DEP_BRA_ZERO_SM] =
                    (double)(1 - TestParam(pname, fname, line,
                                           params[j++], PTYPE_LOGICAL));
                }
              else if (!strcmp(params[j], "var"))
                {
                  /* Variable (to be passed into output) */

                  loc1 = NewItem(loc0 + DEP_BRA_PTR_VAR,
                                 DEP_BRA_VAR_BLOCK_SIZE);

                  /* Update index */

                  j++;

                  /* Read name */

                  if (j == np)
                    Error(loc0, "Missing variable name");
                  else
                    WDB[loc1 + DEP_BRA_VAR_PTR_NAME] =
                      (double)PutText(params[j++]);

                  /* Read value */

                  if (j == np)
                    Error(loc0, "Missing variable value");
                  else
                    WDB[loc1 + DEP_BRA_VAR_PTR_VALUE] =
                      (double)PutText(params[j++]);

                  /* Reset flag for more parameters */

                  n = 0;
                }
              else if (!strcmp(params[j], "repm"))
                {
                  /* Replace material */

                  loc1 = NewItem(loc0 + DEP_BRA_PTR_REPLACE_MAT,
                                 DEP_BRA_REPLACE_MAT_BLOCK_SIZE);

                  /* Update index */

                  j++;

                  /* Read first material */

                  if (j == np)
                    Error(loc0, "Missing first material");
                  else
                    WDB[loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT1] =
                      (double)PutText(params[j++]);

                  /* Read second material */

                  if (j == np)
                    Error(loc0, "Missing second material");
                  else
                    WDB[loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT2] =
                      (double)PutText(params[j++]);

                  /* Reset flag for more parameters */

                  n = 0;
                }
              else if (!strcmp(params[j], "repu"))
                {
                  /* Replace universe */

                  loc1 = NewItem(loc0 + DEP_BRA_PTR_REPLACE_UNI,
                                 DEP_BRA_REPLACE_UNI_BLOCK_SIZE);

                  /* Update index */

                  j++;

                  /* Read first universe */

                  if (j == np)
                    Error(loc0, "Missing first universe");
                  else
                    WDB[loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1] =
                      (double)PutText(params[j++]);

                  /* Read second universe */

                  if (j == np)
                    Error(loc0, "Missing second universe");
                  else
                    WDB[loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI2] =
                      (double)PutText(params[j++]);

                  /* Reset flag for more parameters */

                  n = 0;
                }
              else if (!strcmp(params[j], "tra"))
                {
                  /* Transformation */

                  loc1 = NewItem(loc0 + DEP_BRA_PTR_TRANS,
                                 DEP_BRA_TRANS_BLOCK_SIZE);

                  /* Update index */

                  j++;

                  /* Read universe */

                  if (j == np)
                    Error(loc0, "Missing universe");
                  else
                    WDB[loc1 + DEP_BRA_TRANS_PTR_UNI] =
                      (double)PutText(params[j++]);

                  /* Read second transformation */

                  if (j == np)
                    Error(loc0, "Missing transformation");
                  else
                    WDB[loc1 + DEP_BRA_TRANS_PTR_TRANS] =
                      (double)PutText(params[j++]);

                  /* Reset flag for more parameters */

                  n = 0;
                }
              else if (n == 1)
                {
                  /* Read optional parameters for stp */

                  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                  /* Allocate memory */

                  loc2 = NewItem(loc1 + DEP_BRA_STP_PTR_SAB,
                                 DEP_BRA_STP_SAB_BLOCK_SIZE);

                  /* Read therm */

                  if (j == np)
                    Error(loc0, "Missing thermal scattering data");
                  else
                    WDB[loc2 + DEP_BRA_STP_SAB_PTR_THERM] =
                      (double)PutText(params[j++]);

                  /* Read libraries */

                  if (j == np)
                    Error(loc0, "Missing first S(a,b) library");
                  else
                    WDB[loc2 + DEP_BRA_STP_SAB_PTR_LIB1] =
                      (double)PutText(params[j++]);

                  if (j == np)
                    Error(loc0, "Missing second S(a,b) library");
                  else
                    WDB[loc2 + DEP_BRA_STP_SAB_PTR_LIB2] =
                      (double)PutText(params[j++]);
                }
              else
                Error(loc0, "Input error");
            }
         }

      /***********************************************************************/

      /***** Coefficient calculation *****************************************/

      else if (!strcasecmp(word, "coef"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 1000000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_COEF0, COEF_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Number of burnup points */

          ni = (long)TestParam(pname, fname, line, params[j++],
                               PTYPE_INT, 1, 10000);
          WDB[loc0 + COEF_N_BU] = (double)ni;

          /* Allocate memory for points */

          ptr = ReallocMem(DATA_ARRAY, ni);
          WDB[loc0 + COEF_PTR_BU_PTS] = (double)ptr;

          /* Read values */

          for (n = 0; n < ni; n++)
            {
              if (j == np)
                Error(loc0, "Missing burnup points");
              else
                {
                  /* Test that point is valid */

                  TestParam(pname, fname, line, params[j], PTYPE_REAL,
                            -INFTY, INFTY);

                  /* Store as array (converted later) */

                  WDB[ptr++] = (double)PutText(params[j++]);
                }
            }

          /* Reset cumulative count */

          nr = 1;
          WDB[loc0 + COEF_N_TOT] = 1.0;

          /* Read matrix */

          while (j < np)
            {
              /* New entry */

              loc1 = NewItem(loc0 + COEF_PTR_MTX, COEF_BLOCK_SIZE);

              /* Number of branches */

              ns = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                   1, 1000);

              WDB[loc1 + COEF_MTX_N_BRA] = (double)ns;

              /* Put cumulative count */

              WDB[loc1 + COEF_MTX_N_CUMU] = (double)nr;
              nr = nr*ns;

              /* Allocate memory for branches */

              ptr = ReallocMem(DATA_ARRAY, ns);
              WDB[loc1 + COEF_MTX_PTR_BRA] = (double)ptr;

              /* Read branches */

              for (n = 0; n < ns; n++)
                {
                  if (j == np)
                    Error(loc0, "Missing branches");
                  else
                    WDB[ptr++] = (double)PutText(params[j++]);
                }

              /* Add to total */

              WDB[loc0 + COEF_N_TOT] = RDB[loc0 + COEF_N_TOT]*((double)ns);
            }

          /* Add to total */

          WDB[DATA_TOT_COEF_CALC] = RDB[DATA_TOT_COEF_CALC] + (double)nr;
          WDB[DATA_COEF_CALC_TOT_RUNS] =
            RDB[DATA_COEF_CALC_TOT_RUNS] + (double)(nr*ni);
         }

      /***********************************************************************/

      /***** History variation calculations **********************************/

      else if (!strcasecmp(word, "hisv"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 1000000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_HISV0, HISV_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;
          while (j < np)
            {
              /* New entry */

              loc1 = NewItem(loc0 + HISV_PTR_VAR, COEF_BLOCK_SIZE);

              /* Read burnup point */

              if (j == np)
                Error(loc0, "Missing burnup point");
              else
                WDB[loc1 + HISV_VAR_BU] =
                  TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                            -INFTY, 1E6);

              /* Check zero */

              if (RDB[loc1 + HISV_VAR_BU] == 0.0)
                Error(loc0, "Zero burnup not allowed as variation point");

              /* Avoid compiler warning */

              ns = -1;

              /* Number of branches */

              if (j == np)
                Error(loc0, "Missing number of variations");
              else
                ns = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                     1, 1000);

              WDB[loc1 + HISV_VAR_N_VAR] = (double)ns;

               /* Allocate memory for branches */

              ptr = ReallocMem(DATA_ARRAY, ns);
              WDB[loc1 + HISV_VAR_PTR_VAR] = (double)ptr;

              /* Read branches */

              for (n = 0; n < ns; n++)
                {
                  if (j == np)
                    Error(loc0, "Missing variations");
                  else
                    WDB[ptr++] = (double)PutText(params[j++]);
                }
            }
         }

      /***********************************************************************/

      /***** Case matrix calculation *****************************************/

      else if (!strcasecmp(word, "casematrix"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 6, 1000000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_CASEMTX0, CASEMTX_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Read name */

          WDB[loc0 + CASEMTX_PTR_NAME] = (double)PutText(params[j++]);

          /* Read number of histories */

          ni = (long)TestParam(pname, fname, line, params[j++],
                               PTYPE_INT, 1, 10000);
          WDB[loc0 + CASEMTX_N_HIS] = (double)ni;

          /* Allocate memory for array */

          ptr = ReallocMem(DATA_ARRAY, ni);
          WDB[loc0 + CASEMTX_PTR_HIS] = (double)ptr;

          /* Read values */

          for (n = 0; n < ni; n++)
            {
              if (j == np)
                Error(loc0, "Missing data");
              else
                WDB[ptr++] = (double)PutText(params[j++]);
            }

          /* Create new item */

          loc1 = NewItem(loc0 + CASEMTX_PTR_COEF, COEF_BLOCK_SIZE);

          /* Number of burnup points */

          ni = (long)TestParam(pname, fname, line, params[j++],
                               PTYPE_INT, 1, 10000);
          WDB[loc1 + COEF_N_BU] = (double)ni;

          /* Allocate memory for points */

          ptr = ReallocMem(DATA_ARRAY, ni);
          WDB[loc1 + COEF_PTR_BU_PTS] = (double)ptr;

          /* Read values */

          for (n = 0; n < ni; n++)
            {
              if (j == np)
                Error(loc0, "Missing burnup points");
              else
                {
                  /* Test that point is valid */

                  TestParam(pname, fname, line, params[j], PTYPE_REAL,
                            -INFTY, INFTY);

                  /* Store as array (converted later) */

                  WDB[ptr++] = (double)PutText(params[j++]);
                }
            }

          /* Reset cumulative count */

          nr = 1;
          WDB[loc1 + COEF_N_TOT] = 1.0;

          /* Read matrix */

          while (j < np)
            {
              /* New entry */

              loc2 = NewItem(loc1 + COEF_PTR_MTX, COEF_BLOCK_SIZE);

              /* Number of branches */

              ns = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                   1, 1000);

              WDB[loc2 + COEF_MTX_N_BRA] = (double)ns;

              /* Put cumulative count */

              WDB[loc2 + COEF_MTX_N_CUMU] = (double)nr;
              nr = nr*ns;

              /* Allocate memory for branches */

              ptr = ReallocMem(DATA_ARRAY, ns);
              WDB[loc2 + COEF_MTX_PTR_BRA] = (double)ptr;

              /* Read branches */

              for (n = 0; n < ns; n++)
                {
                  if (j == np)
                    Error(loc0, "Missing branches");
                  else
                    WDB[ptr++] = (double)PutText(params[j++]);
                }

              /* Add to total */

              WDB[loc1 + COEF_N_TOT] = RDB[loc1 + COEF_N_TOT]*((double)ns);
            }

          /* Add to total */

          WDB[loc0 + CASEMTX_TOT_COEF_CALC] =
            RDB[loc0 + CASEMTX_TOT_COEF_CALC] + (double)nr;
          WDB[loc0 + CASEMTX_COEF_CALC_TOT_RUNS] =
            RDB[loc0 + CASEMTX_COEF_CALC_TOT_RUNS] + (double)(nr*ni);
         }

      /***********************************************************************/

      /***** Mixture definition **********************************************/

      else if (!strcasecmp(word, "mix"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, 4*MAX_ISOTOPES + 8,
                         fname);

          /* Get memory size */

          mem = RDB[DATA_TOTAL_BYTES];

          /* Create new item */

          loc0 = NewItem(DATA_PTR_M0, MATERIAL_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Material name */

          WDB[loc0 + MATERIAL_PTR_NAME] = (double)PutText(params[j++]);

          /* Reset temperatures and volume */

          WDB[loc0 + MATERIAL_DOPPLER_TEMP] = -1.0;
          WDB[loc0 + MATERIAL_TMS_TMIN] = INFTY;
          WDB[loc0 + MATERIAL_TMS_TMAX] = -INFTY;
          WDB[loc0 + MATERIAL_VOLUME_GIVEN] = -INFTY;

          /* Loop over parameters */

          while (j < np)
            {
              /* Check parameter */

              if (!strcmp(params[j], "rgb"))
                {
                  /***** Mixture colour **************************************/

                  j++;

                  /* Get r, b and g */

                  r = (long)TestParam(pname, fname, line, params[j++],
                                      PTYPE_INT, 0, 255);

                  g = (long)TestParam(pname, fname, line, params[j++],
                                      PTYPE_INT, 0, 255);

                  b = (long)TestParam(pname, fname, line, params[j++],
                                      PTYPE_INT, 0, 255);

                  /* Set color */

                  WDB[loc0 + MATERIAL_RGB] = (double)(b + 1000*g + 1000000*r);

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "vol"))
                {
                  /***** Mixture volume **************************************/

                  j++;

                  /* Get volume */

                  WDB[loc0 + MATERIAL_VOLUME_GIVEN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);

                  /***********************************************************/
                }
              else if (!strcmp(params[j], "mass"))
                {
                  /***** Mixture mass ****************************************/

                  j++;

                  /* Get mass */

                  WDB[loc0 + MATERIAL_MASS_GIVEN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);

                  /***********************************************************/
                }
              else
                {
                  /***** Mixed materials *************************************/

                  /* Create new item */

                  loc1 = NewItem(loc0 + MATERIAL_PTR_MIX, MIXTURE_BLOCK_SIZE);

                  /* Check number of parameters */

                  if (j > np - 2)
                    Error(loc0, "Invalid number of parameters");

                  /* Read nuclide name */

                  WDB[loc1 + MIXTURE_PTR_MAT] =
                    (double)PutText(params[j++]);

                  /* Read fraction (read to volume fraction) */

                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -1E+25, 1E+25);

                  WDB[loc1 + MIXTURE_VFRAC] = val;

                  /***********************************************************/
                }

            }

          /* Put memory size */

          WDB[loc0 + MATERIAL_MEM_SIZE] = RDB[DATA_TOTAL_BYTES] - mem;
         }

      /***********************************************************************/

      /***** Thermal scattering data *****************************************/

      else if ((!strcasecmp(word, "therm")) ||
               (!strcasecmp(word, "thermstoch")))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          /* Maximum number of parameters = 50 (very many libraries can  */
          /* be used for otf treatment) */

          params = GetParams(word, input, &np, &i0, 2, 50, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_T0, THERM_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line ;

          /* Read data */

          j = 0;

          /* Read name */

          WDB[loc0 + THERM_PTR_ALIAS] = (double)PutText(params[j++]);

          if (np == 2)
            {
              /* No temperature data. Create one SAB nuclide */

              loc1 = NewItem(loc0 + THERM_PTR_SAB, SAB_BLOCK_SIZE);

              /* Read isotope name */

              WDB[loc1 + SAB_PTR_NAME] = (double)PutText(params[j++]);

              /* Reset temperature */

              WDB[loc0 + THERM_T] = -1.0;
              WDB[loc0 + THERM_INTERP_MODE] = THERM_INTERP_NONE;
            }

          else if (np >= 4)
            {
              /* Temperature data. Read T */

              /* T�ll� l�mp�tilalla ei ole mit��n v�li� on-the-fly */
              /* -tapauksessa, koska interface antaa l�mp�tilan */

              val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, 1E+6);
              WDB[loc0 + THERM_T] = val;

              /* Check that number of libraries makes sense */

              if (val > 0.0 && np != 4)
                Error(loc0,
                      "Invalid number of libraries for interpolation in pre-processing phase or T > 0 in on-the-fly mode");

              /* Read library names */

              for (i = 2; i < np; i++)
                {
                  /* New S(a,b) nuclide */

                  loc1 = NewItem(loc0 + THERM_PTR_SAB, SAB_BLOCK_SIZE);

                  /* Read name of nuclide  */

                  WDB[loc1 + SAB_PTR_NAME] = (double)PutText(params[j++]);
                }

              /* Set interpolation mode */

              if (RDB[loc0 + THERM_T] == 0.0)
                WDB[loc0 + THERM_INTERP_MODE] = (double)THERM_INTERP_OTF;
              else if (!strcasecmp(word, "thermstoch"))
                WDB[loc0 + THERM_INTERP_MODE] = (double)THERM_INTERP_STOCHMIX;
              else if (!strcasecmp(word, "therm"))
                WDB[loc0 + THERM_INTERP_MODE] = (double)THERM_INTERP_MAKXSF;
            }
          else
            Error(loc0, "Invalid number of parameters");

         }

      /***********************************************************************/

      /***** Surface definition **********************************************/

      else if (!strcasecmp(word, "surf"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, MAX_SURFACE_PARAMS + 2,
                             fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_S0, SURFACE_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read parameters */

          j = 0;

          /* Surface name */

          WDB[loc0 + SURFACE_PTR_NAME] = (double)PutText(params[j++]);

          /* Find type */

          for (type = 0; type < SURFACE_TYPES + 1; type++)
            if (!strcmp(params[j], surf_types[type]))
              break;

          /* Check */

          if (type == SURFACE_TYPES + 1)
            Error(loc0, "Surface type %s does not exist", params[j]);

          /* Put type */

          WDB[loc0 + SURFACE_TYPE] = (double)type + 1.0;

          /* Update counter */

          j++;

          /* Check number of parameters */

          if (np - j < surf_params[type][0])
            Error(loc0, "Not enough parameters given for surface type \"%s\"",
                  surf_types[type]);
          else if (np - j > surf_params[type][1])
            Error(loc0, "Too many parameters given for surface type \"%s\"",
                  surf_types[type]);

          /* Set number of parameters */

          WDB[loc0 + SURFACE_N_PARAMS] = (double)(np - j);

          /* Get parameters */

          if (np - j > 0)
            {
              /* Allocate memory for parameters and set pointer */

              ptr = ReallocMem(DATA_ARRAY, np - j);
              WDB[loc0 + SURFACE_PTR_PARAMS] = (double)ptr;

              /* Read parameters */

              for (n = j; n < np; n++)
                WDB[ptr++] = TestParam(pname, fname, line, params[j++],
                                        PTYPE_REAL, -INFTY, INFTY);
            }
          else
            WDB[loc0 + SURFACE_PTR_PARAMS] = NULLPTR;

          /* Calculate numeric value (NOTE: This is only to speed up */
          /* the search loops in processcells.c) */

          n = NumericStr(GetText(loc0 + SURFACE_PTR_NAME));
          WDB[loc0 + SURFACE_NAME_NUMERIC] = (double)n;
         }

      /***********************************************************************/

      /***** Cell ************************************************************/

      else if (!strcasecmp(word, "cell"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, 6 + MAX_CELL_SURFACES,
                         fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_C0, CELL_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Read name */

          if (j == np)
            Error(loc0, "Missing cell name");
          else
            WDB[loc0 + CELL_PTR_NAME] =  (double)PutText(params[j++]);

          /* Read universe */

          if (j == np)
            Error(loc0, "Missing cell universe");
          else
            WDB[loc0 + CELL_PTR_UNI] = (double)PutText(params[j++]);

          /* Material or fill */

          if (j == np)
            Error(loc0, "Missing cell material");
          else if (!strcmp(params[j], "fill"))
            {
              if (j == np - 1)
                Error(loc0, "Missing fill universe");
              else
                WDB[loc0 + CELL_PTR_FILL] = (double)PutText(params[++j]);

              WDB[loc0 + CELL_PTR_MAT] = NULLPTR;
            }
          else
            {
              WDB[loc0 + CELL_PTR_MAT] = (double)PutText(params[j]);
              WDB[loc0 + CELL_PTR_FILL] = NULLPTR;
            }

          /* Allocate memory for collision counter */

          AllocValuePair(loc0 + CELL_COL_COUNT);

          /* Update index */

          j++;

          /* Check */

          if (j == np)
            WDB[loc0 + CELL_PTR_SURF_LIST] = NULLPTR;
          else
            {
              /* Read surface list to a string */

              m = 0;
              str[0] = '\0';

              for (n = j; n < np; n++)
                {
                  /* Check string length */

                  if (strlen(str) > 19999)
                    Die(FUNCTION_NAME, "Temporary string lenght exceeded");

                  /* Check leading minus sign */

                  if (params[j][0] == '-')
                    sprintf(&str[m], " - %s ", &params[j][1]);
                  else
                    sprintf(&str[m], " %s ", params[j]);

                  /* Update counters */

                  m = strlen(str);

                  j++;
                }

              /* Put pointer */

              WDB[loc0 + CELL_PTR_SURF_LIST] = (double)PutText(str);
            }
        }

      /***********************************************************************/

      /***** Solid 3D geometry type ******************************************/

      /* NOTE: Tähän myös tuo OpenFOAM -pohjainen geometriatyyppi. */
      /* Keksitään parempi nimi myöhemmin */

      else if (!strcasecmp(word, "solid"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, MAX_INPUT_PARAMS, fname);

          /* Read data */

          j = 0;

          /* Avoid compiler warning */

          type = -1;

          /* Read type */

          if (j == np)
            Error(-1, pname, fname, line, "Missing solid geometry 3D type");
          else
            type = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                   1, 3);

          /* Read name */

          if (j == np)
            Error(-1, pname, fname, line,  "Missing universe name");
          else
            WDB[DATA_DUMMY] =  (double)PutText(params[j++]);

          /* Check type */

          if (type == SOLID_GEO_TYPE_OF)
            {
              /* Create new item */

              loc0 = NewItem(DATA_PTR_UMSH0, UMSH_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Put universe name */

              WDB[loc0 + UMSH_PTR_NAME] = RDB[DATA_DUMMY];

              /* Read background universe */

              if (j < np)
                WDB[loc0 + UMSH_PTR_BG_UNIV] = (double)PutText(params[j++]);
              else
                Error(loc0, "Missing background universe name");

              /* Create associated inteface structure */

              loc1 = NewItem(loc0 + UMSH_PTR_IFC, IFC_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc1 + PARAM_PTR_NAME] = (double)PutText(word);
              WDB[loc1 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc1 + PARAM_LINE] = (double)line;

              /* Read search mesh split limit */

              if (j < np)
                WDB[loc1 + IFC_SEARCH_MESH_ADA_SPLIT] =
                  TestParam(pname, fname, line, params[j++],
                            PTYPE_INT, 1, 10000000000);
              else
                Error(loc0, "Missing search mesh split limit");

              /* Avoid compiler warning */

              ns = -1;

              /* Read search mesh dimensions */

              if (j < np)
                ns = (long)TestParam(pname, fname, line, params[j++],
                                     PTYPE_INT, 1, 10);
              else
                Error(loc0, "Missing search mesh dimension");

              /* Allocate memory */

              ptr = ReallocMem(DATA_ARRAY, ns + 1);
              WDB[loc1 + IFC_SEARCH_MESH_ADA_PTR_SZ] = (double)ptr;

              /* Read size */

              for (n = 0; n < ns; n++)
                {
                  if (j < np)
                    WDB[ptr++] = TestParam(pname, fname, line, params[j++],
                                           PTYPE_INT, 2, 10000);
                  else
                    Error(loc0, "Missing mesh size");
                }

              /* Put terminator */

              WDB[ptr] = -1.0;

              /* Read points file */

              if (j < np)
                WDB[loc1 + IFC_PTR_OF_PFILE] =
                  (double)PutText(params[j++]);
              else
                Error(loc0, "Missing points file");

              /* Read faces file */

              if (j < np)
                WDB[loc1 + IFC_PTR_OF_FFILE] =
                  (double)PutText(params[j++]);
              else
                Error(loc0, "Missing faces file");

              /* Read owner file */

              if (j < np)
                WDB[loc1 + IFC_PTR_OF_OFILE] =
                  (double)PutText(params[j++]);
              else
                Error(loc0, "Missing owner file");

              /* Read neighbour file */

              if (j < np)
                WDB[loc1 + IFC_PTR_OF_NFILE] =
                  (double)PutText(params[j++]);
              else
                Error(loc0, "Missing neighbour file");

              /* Read materials file */

              if (j < np)
                WDB[loc1 + IFC_PTR_OF_MFILE] =
                  (double)PutText(params[j++]);
              else
                Error(loc0, "Missing materials file");
            }
          else if (type == SOLID_GEO_TYPE_STL)
            {
              /* STL type geometry */

              loc0 = NewItem(DATA_PTR_STL0, STL_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Put universe name */

              WDB[loc0 + STL_PTR_NAME] = RDB[DATA_DUMMY];

              /* Read background universe */

              if (j == np)
                Error(loc0, "Missing background universe");
              else
                WDB[loc0 + STL_PTR_BG_UNIV] =  (double)PutText(params[j++]);

              /* Read search mesh split limit */

              if (j < np)
                WDB[loc0 + STL_SEARCH_MESH_ADA_SPLIT] =
                  TestParam(pname, fname, line, params[j++],
                            PTYPE_INT, 1, 10000000000);
              else
                Error(loc0, "Missing search mesh split limit");

              /* Avoid compiler warning */

              ns = -1;

              /* Read search mesh dimensions */

              if (j < np)
                ns = (long)TestParam(pname, fname, line, params[j++],
                                     PTYPE_INT, 1, 10);
              else
                Error(loc0, "Missing search mesh dimension");

              /* Allocate memory */

              ptr = ReallocMem(DATA_ARRAY, ns + 1);
              WDB[loc0 + STL_SEARCH_MESH_PTR_SZ] = (double)ptr;

              /* Read size */

              for (n = 0; n < ns; n++)
                {
                  if (j < np)
                    WDB[ptr++] = TestParam(pname, fname, line, params[j++],
                                           PTYPE_INT, 2, 10000);
                  else
                    Error(loc0, "Missing mesh size");
                }

              /* Put terminator */

              WDB[ptr] = -1.0;

              /* Read search mode */

              if (j == np)
                Error(loc0, "Missing STL search mode");
              else
                WDB[loc0 + STL_SEARCH_MODE] =
                  TestParam(pname, fname, line, params[j++], PTYPE_INT, 1, 2);

              /* Read point merge radius */

              if (j == np)
                Error(loc0, "Missing STL point merge radius");
              else
                WDB[loc0 + STL_MERGE_RAD] =
                  TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                            -1.0, 1.0);

              /* Check index */

              if (j == np)
                Error(loc0, "Missing STL geometry files");

              /* Read data */

              while (j < np)
                {
                  /* Check word */

                  if (!strcasecmp(params[j], "file"))
                    {
                      /* Update index */

                      j++;

                      /* New file */

                      ptr = NewItem(loc0 + STL_PTR_FILES, STL_FILE_BLOCK_SIZE);

                      /* Read solid name name */

                      if (j == np)
                        Error(loc0, "Missing body name");
                      else
                        WDB[ptr + STL_FILE_PTR_SNAME] =
                          (double)PutText(params[j++]);

                      /* Read file name */

                      if (j == np)
                        Error(loc0, "Missing file name");
                      else
                        WDB[ptr + STL_FILE_PTR_FNAME] =
                          (double)PutText(params[j++]);

                      /* Read scaling factor */

                      if (j == np)
                        Error(loc0, "Missing scaling factor for file \"%s\"",
                              GetText(ptr + STL_FILE_PTR_FNAME));
                      else
                        WDB[ptr + STL_FILE_SCALING] =
                          TestParam(pname, fname, line, params[j++],
                                    PTYPE_REAL, ZERO, INFTY);

                      /* Read origin */

                      if (j > np - 3)
                        Error(loc0, "Missing origin for file \"%s\"",
                              GetText(ptr + STL_FILE_PTR_FNAME));

                      WDB[ptr + STL_FILE_X0] =
                        TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);

                      WDB[ptr + STL_FILE_Y0] =
                        TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);

                      WDB[ptr + STL_FILE_Z0] =
                        TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);
                    }
                  else if (!strcasecmp(params[j], "body"))
                    {
                      /* Update index */

                      j++;

                      /* New body */

                      ptr = NewItem(loc0 + STL_PTR_BODIES,
                                    STL_BODY_BLOCK_SIZE);

                      /* Reset material and fill pointes */

                      WDB[ptr + STL_BODY_PTR_MNAME] = NULLPTR;
                      WDB[ptr + STL_BODY_PTR_FILL] = NULLPTR;

                      /* Read solid name */

                      if (j == np)
                        Error(loc0, "Missing body name");
                      else
                        WDB[ptr + STL_BODY_PTR_BNAME] =
                          (double)PutText(params[j++]);

                      /* Read cell name */

                      if (j == np)
                        Error(loc0, "Missing cell name for body \"%s\"",
                              GetText(ptr + STL_BODY_PTR_BNAME));
                      else
                        WDB[ptr + STL_BODY_PTR_CNAME] =
                          (double)PutText(params[j++]);

                      /* Read material name */

                      if (j == np)
                        Error(loc0, "Missing material name for body \"%s\"",
                              GetText(ptr + STL_BODY_PTR_BNAME));
                      else if (!strcasecmp(params[j], "fill"))
                        {
                          /* Update index */

                          j++;

                          /* Read universe name */

                          if (j == np)
                            Error(loc0, "Missing universe for body \"%s\"",
                                  GetText(ptr + STL_BODY_PTR_BNAME));
                          else
                            WDB[ptr + STL_BODY_PTR_FILL] =
                              (double)PutText(params[j++]);
                        }
                      else
                        WDB[ptr + STL_BODY_PTR_MNAME] =
                          (double)PutText(params[j++]);
                    }
                  else
                    Error(loc0, "Entry \"%s\" should be \"file\" or \"body\"",
                          params[j]);
                }
            }
          else if (type == SOLID_GEO_TYPE_IFC)
            {

              /* Create new item */

              loc0 = NewItem(DATA_PTR_UMSH0, UMSH_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Open input file to read universe and BG universe */

              if ((fp = fopen(GetText(DATA_DUMMY), "r")) == NULL)
                Error(loc0,
                      "Multi-physics interface file \"%s\" does not exist",
                      GetText(DATA_DUMMY));

              /* Interface type */

              if (fscanf(fp, "%ld", &m) == EOF)
                Die(FUNCTION_NAME, "fscanf error");

              if(m != IFC_TYPE_OF_SOLID)
                Error(loc0,
                      "Interface type %ld not supported with UMSH solids", m);

              /* Get universe name */

              if (fscanf(fp, "%s", str) == EOF)
                Die(FUNCTION_NAME, "fscanf error");

              /* Put universe name */

              WDB[loc0 + UMSH_PTR_NAME] = (double)PutText(str);

              /* Get background universe name */

              if (fscanf(fp, "%s", str) == EOF)
                Die(FUNCTION_NAME, "fscanf error");

              /* Put BG universe name */

              WDB[loc0 + UMSH_PTR_BG_UNIV] = (double)PutText(str);

              /* Read other things in readifcofmesh.c */

              fclose(fp);

              /* Get memory size before creation of the interface */

              mem0 = RDB[DATA_ALLOC_MAIN_SIZE];

              /* Create new item */
              /* This can be updated, so it is stored to the main IFC-list */

              loc1 = NewItem(DATA_PTR_IFC0, IFC_BLOCK_SIZE);

              /* Put interface to UMSH */

              WDB[loc0 + UMSH_PTR_IFC] = (double)loc1;

              /* Put name, file name and line number for IFC */

              WDB[loc1 + PARAM_PTR_NAME] = (double)PutText(word);
              WDB[loc1 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc1 + PARAM_LINE] = (double)line;

              /* Read data */

              j = 0;

              /* Store file name */

              WDB[loc1 + IFC_PTR_INPUT_FNAME] = RDB[DATA_DUMMY];

              /* Store file name to UMSH block */

              WDB[loc0 + UMSH_PTR_FNAME] = RDB[DATA_DUMMY];

              /* Read interface */
              /* Time binning is not yet set in dynamic simulations */

              ReadInterface(loc1, (long)NO);

              /* Get memory size after creation of the interface */

              mem1 = RDB[DATA_ALLOC_MAIN_SIZE];

              /* Store ifc-memory size to interface (for MPI transfer) */

              WDB[loc1 + IFC_MEM_SIZE] = mem1 - mem0;

              /* Set flag */

              WDB[DATA_USE_DENSITY_FACTOR] = (double)YES;
            }
          else
            Die(FUNCTION_NAME, "Invalid geometry type");
        }

      /***********************************************************************/

      /***** Geometry nest ***************************************************/

      else if ((!strcasecmp(word, "nest")) || (!strcasecmp(word, "pin")) ||
               (!strcasecmp(word, "particle")))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, 10000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_NST0, NEST_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Read name */

          WDB[loc0 + NEST_PTR_NAME] = (double)PutText(params[j]);

          /* Special treatment for pin and particle */

          if (!strcmp(word, "pin"))
            type = SURF_CYL - 1;
          else if (!strcmp(word, "particle"))
            type = SURF_SPH - 1;
          else
            {
              /* Update index */

              j++;

              /* Check if surface type is given */

              for (type = 0; type < SURFACE_TYPES + 1; type++)
                if (!strcmp(params[j], surf_types[type]))
                  break;
            }

          if (type < SURFACE_TYPES)
            {
              /***************************************************************/

              /***** All regions have the same type **************************/

              /* Check */

              if (surf_params[type][0] < 1)
                Error(loc0, "Surface type %s not allowed in nests", params[j]);

              /* Put type */

              WDB[loc0 + NEST_TYPE] = (double)(type + 1);

              /* Reset pointer */

              loc1 = -1;

              /* Reset minimum radius */

              r0 = -INFTY;

              /* Loop over parameters */

              j++;

              while (j < np)
                {
                  /* Create new region */

                  loc1 = NewItem(loc0 + NEST_PTR_REGIONS, NEST_REG_BLOCK_SIZE);

                  /* Check fill entry */

                  if (!strcmp(params[j], "fill"))
                    {
                      /* Put filler universe */

                      WDB[loc1 + NEST_REG_PTR_FILL] =
                        (double)PutText(params[++j]);
                      WDB[loc1 + NEST_REG_PTR_CELL] = NULLPTR;
                    }
                  else
                    {
                      /* Put material name in cell pointer */

                      WDB[loc1 + NEST_REG_PTR_CELL] =
                        (double)PutText(params[j]);
                      WDB[loc1 + NEST_REG_PTR_FILL] = NULLPTR;
                    }

                  /* Update counter */

                  j++;

                  /* Check if last */

                  if (j < np)
                    {
                      /* Create surface */

                      loc2 = NewItem(loc1 + NEST_REG_PTR_SURF_IN,
                                     SURFACE_BLOCK_SIZE);

                      /* Put pointer */

                      WDB[loc1 + NEST_REG_PTR_SURF_IN] = (double)loc2;

                      /* Put surface type */

                      WDB[loc2 + SURFACE_TYPE] = (double)(type + 1);

                      /* Create parameter list */

                      ptr = ReallocMem(DATA_ARRAY, surf_params[type][0]);

                      /* Put pointer */

                      WDB[loc2 + SURFACE_PTR_PARAMS] = (double)ptr;

                      /* Put number of parameters */

                      WDB[loc2 + SURFACE_N_PARAMS] =
                        (double)surf_params[type][0];

                      /* Get surface parameter */

                      if ((type == SURF_PX) || (type == SURF_PY) ||
                          (type == SURF_PZ))
                        r1 = TestParam(pname, fname, line, params[j++],
                                       PTYPE_REAL, -INFTY, INFTY);
                      else
                        r1 = TestParam(pname, fname, line, params[j++],
                                       PTYPE_REAL, 0.0, INFTY);

                      /* Compare to previous */

                      if (r1 < r0)
                        Error(loc0, "regions must be nested");
                      else
                        r0 = r1;

                      /* Put surface parameter */

                      WDB[ptr + surf_params[type][0] - 1] = r1;
                    }
                }

              /* Check that last region has no surface */

              if (loc1 < 0)
                Die(FUNCTION_NAME, "Pointer error");
              else if ((long)RDB[loc1 + NEST_REG_PTR_SURF_IN] > 0)
                Error(loc0, "Last region must be unbound");

              /***************************************************************/
            }
          else
            {
              /***************************************************************/

              /***** Different type for each region **************************/

              /* Put type */

              WDB[loc0 + NEST_TYPE] = -1.0;

              /* Reset pointer */

              loc1 = -1;

              /* Loop over parameters */

              while (j < np)
                {
                  /* Create new region */

                  loc1 = NewItem(loc0 + NEST_PTR_REGIONS,
                                 NEST_REG_BLOCK_SIZE);

                  /* Check fill entry */

                  if (!strcmp(params[j], "fill"))
                    {
                      /* Put filler universe */

                      WDB[loc1 + NEST_REG_PTR_FILL] =
                        (double)PutText(params[++j]);
                      WDB[loc1 + NEST_REG_PTR_CELL] = NULLPTR;
                    }
                  else
                    {
                      /* Put material in cell pointer */

                      WDB[loc1 + NEST_REG_PTR_CELL] =
                        (double)PutText(params[j]);
                      WDB[loc1 + NEST_REG_PTR_FILL] = NULLPTR;
                    }

                  /* Update counter */

                  j++;

                  /* Check if last */

                  if (j < np)
                    {
                      /* Get surface type */

                      for (type = 0; type < SURFACE_TYPES + 1; type++)
                        if (!strcmp(params[j], surf_types[type]))
                          break;

                      /* Check */

                      if (type == SURFACE_TYPES + 1)
                        Error(loc0, "Surface type %s does not exist",
                              params[j]);
                      else if (surf_params[type][0] < 1)
                        Error(loc0, "Surface type %s not allowed in nests",
                              params[j]);

                      /* Update counter */

                      j++;

                      /* Create surface */

                      loc2 = NewItem(loc1 + NEST_REG_PTR_SURF_IN,
                                     SURFACE_BLOCK_SIZE);

                      /* Put pointer */

                      WDB[loc1 + NEST_REG_PTR_SURF_IN] = (double)loc2;

                      /* Put surface type */

                      WDB[loc2 + SURFACE_TYPE] = (double)(type + 1);

                      /* Create parameter list */

                      ptr = ReallocMem(DATA_ARRAY, surf_params[type][0]);

                      /* Put pointer */

                      WDB[loc2 + SURFACE_PTR_PARAMS] = (double)ptr;

                      /* Put number of parameters */

                      WDB[loc2 + SURFACE_N_PARAMS] =
                        (double)surf_params[type][0];

                      /* Read surface parameters */

                      for (n = 0; n < surf_params[type][0]; n++)
                        WDB[ptr++] = TestParam(pname, fname, line,
                                                params[j++], PTYPE_REAL,
                                                -INFTY, INFTY);
                    }
                }

              /* Check that last region has no surface */

              if (loc1 < 0)
                Die(FUNCTION_NAME, "Pointer error");
              else if ((long)RDB[loc1 + NEST_REG_PTR_SURF_IN] > 0)
                Error(loc0, "Last region must be unbound");

              /***************************************************************/
            }
        }

      /***********************************************************************/

      /***** Coordinate transformations **************************************/

      else if (!strcasecmp(word, "trans") || (!strcasecmp(word, "utrans")) ||
               (!strcasecmp(word, "strans")) || (!strcasecmp(word, "ftrans")) ||
               (!strcasecmp(word, "dtrans")) || (!strcasecmp(word, "ltrans")) ||
               (!strcasecmp(word, "transv")) || (!strcasecmp(word, "transa")) ||
               (!strcasecmp(word, "transb")))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, 18, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_TR0, TRANS_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Set default order */

          WDB[loc0 + TRANS_ORDER] = (double)TRANS_ORDER_TRANS_ROT;

          /* Reset burnup */

          WDB[loc0 + TRANS_BURNUP] = -INFTY;

          /* Read data */

          j = 0;

          /* Set type */

          if (!strncasecmp(pname, "trans", 5))
            {
              /* Check velocity "transv" or acceleration "transa" */

              if (!strcasecmp(word, "transv"))
                WDB[loc0 + TRANS_MOVE_TYPE] = (double)TRANS_MOVE_VEL;
              else if (!strcasecmp(word, "transa"))
                WDB[loc0 + TRANS_MOVE_TYPE] = (double)TRANS_MOVE_ACC;
              else if (!strcasecmp(word, "transb"))
                {
                  /* Read burnup */

                  WDB[loc0 + TRANS_BURNUP] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -1E9, 10000.0);
                }

              /* Check second parameter */

              if (!strcasecmp(params[j], "s"))
                {
                  WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_SURF;
                  j++;
                }
              else if (!strcasecmp(params[j], "f"))
                {
                  WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_FILL;
                  j++;
                }
              else if (!strcasecmp(params[j], "u"))
                {
                  WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_UNI;
                  j++;
                }
              else if (!strcasecmp(params[j], "d"))
                {
                  WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_DET;
                  j++;
                }
              else if (!strcasecmp(params[j], "l"))
                {
                  WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_LAT;
                  j++;
                }
              else
                WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_UNI;
            }
          else if (!strcasecmp(pname, "strans"))
            WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_SURF;
          else if (!strcasecmp(pname, "ftrans"))
            WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_FILL;
          else if (!strcasecmp(pname, "utrans"))
            WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_UNI;
          else if (!strcasecmp(pname, "dtrans"))
            WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_DET;
          else if (!strcasecmp(pname, "ltrans"))
            WDB[loc0 + TRANS_TYPE] = (double)TRANS_TYPE_LAT;
          else
            Die(FUNCTION_NAME, "Invalid transformation type");

          /* Reset level */

          WDB[loc0 + TRANS_PTR_LVL] = -1.0;

          /* Transformation name */

          WDB[loc0 + TRANS_PTR_NAME] = (double)PutText(params[j++]);

          /* Put universe, surface or detector name */

          if ((long)RDB[loc0 + TRANS_TYPE] == TRANS_TYPE_UNI)
            WDB[loc0 + TRANS_PTR_UNI] = RDB[loc0 + TRANS_PTR_NAME];
          else if ((long)RDB[loc0 + TRANS_TYPE] == TRANS_TYPE_FILL)
            WDB[loc0 + TRANS_PTR_CELL] = RDB[loc0 + TRANS_PTR_NAME];
          else if ((long)RDB[loc0 + TRANS_TYPE] == TRANS_TYPE_DET)
            WDB[loc0 + TRANS_PTR_DET] = RDB[loc0 + TRANS_PTR_NAME];
          else if ((long)RDB[loc0 + TRANS_TYPE] == TRANS_TYPE_LAT)
            WDB[loc0 + TRANS_PTR_LAT] = RDB[loc0 + TRANS_PTR_NAME];
          else
            WDB[loc0 + TRANS_PTR_SURF] = RDB[loc0 + TRANS_PTR_NAME];

          /* Read index for lattice transformation */

          if ((long)RDB[loc0 + TRANS_TYPE] == TRANS_TYPE_LAT)
            WDB[loc0 + TRANS_LAT_IDX] =
              TestParam(pname, fname, line, params[j++], PTYPE_INT,
                        1, MAX_LATTICE_ITEMS);

          /* Check if time limits of transformation are given             */
          /* NOTE: Nää rajat pätee nyt ainoastaan noihin transv ja transa */
          /*       transformaatioihin */

          if (!strncasecmp(params[j], "tlim", 4))
            {
              /* Increment j */

              j++;

              /* Check that there is at least six more parameters */

              if (np - j < 4)
                Error(loc0, "Invalid number of parameters given, was expecting "
                      "at least four parameters after \"tlim\"");

              /* Get transformation beginning time (inclusive) */

              WDB[loc0 + TRANS_T0] = TestParam(pname, fname, line, params[j++],
                                               PTYPE_REAL, 0.0, INFTY);

              /* Get transformation ending time (exclusive) */

              WDB[loc0 + TRANS_T1] = TestParam(pname, fname, line,
                                               params[j++], PTYPE_REAL,
                                               WDB[loc0 + TRANS_T0], INFTY);

              /* Get transformation time limit type */

              WDB[loc0 + TRANS_T_TYPE] = TestParam(pname, fname, line,
                                                   params[j++], PTYPE_INT,
                                                   1, 3);
            }

          /* Check if next parameter is a keyword */

          if (!strcmp(params[j], "rot"))
            {
              /* Rotation around general axis */

              /* Increment parameter */

              j++;

              /* Check remaining number of parameters */

              if (np - j != 7)
                Error(loc0, "Invalid number of parameters given");

              /* Get first point */

              WDB[loc0 + TRANS_X0] = TestParam(pname, fname, line,
                                               params[j++], PTYPE_REAL,
                                               -INFTY, INFTY);

              WDB[loc0 + TRANS_Y0] = TestParam(pname, fname, line,
                                               params[j++], PTYPE_REAL,
                                               -INFTY, INFTY);

              WDB[loc0 + TRANS_Z0] = TestParam(pname, fname, line,
                                               params[j++], PTYPE_REAL,
                                               -INFTY, INFTY);

              /* Get direction vector */

              WDB[loc0 + TRANS_ROT_AX_U] = TestParam(pname, fname, line,
                                                     params[j++], PTYPE_REAL,
                                                     -INFTY, INFTY);

              WDB[loc0 + TRANS_ROT_AX_V] = TestParam(pname, fname, line,
                                                     params[j++], PTYPE_REAL,
                                                     -INFTY, INFTY);

              WDB[loc0 + TRANS_ROT_AX_W] = TestParam(pname, fname, line,
                                                     params[j++], PTYPE_REAL,
                                                     -INFTY, INFTY);

              /* Get rotation speed */

              WDB[loc0 + TRANS_ROT_ANG] = TestParam(pname, fname, line,
                                                    params[j++], PTYPE_REAL,
                                                    -INFTY, INFTY);

              /* Set rotation flag */

              WDB[loc0 + TRANS_ROT] = (double)YES;

              /* Set explicit rotation matrix flag */

              WDB[loc0 + TRANS_EXPLI] = (double)YES;
            }
          else
            {
              /* Other transformation */

              /* Check number of remaining parameters */

              if ((np - j == 1) || (np - j == 4))
                {
                  /* Level transformation */

                  n = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                      0, MAX_GEOMETRY_LEVELS);
                  WDB[loc0 + TRANS_PTR_LVL] = (double)n;
                }

              /* Check remaining number of parameters */

              if ((np - j != 0) && (np - j != 3) && (np - j != 6) && (np - j != 12)
                  && (np - j != 7) && (np - j != 13))
                Error(loc0,"Invalid number of parameters given");

              /* Check translation */

              if (np - j > 2)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + TRANS_X0] = val;

                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + TRANS_Y0] = val;

                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + TRANS_Z0] = val;
                }

              /* Check rotation */

              if (np - j > 2)
                {
                  /* Put flag */

                  WDB[loc0 + TRANS_ROT] = (double)YES;

                  /* Check type */

                  if ((np - j == 9) || (np - j == 10))
                    {
                      /***********************************************************/

                      /***** Explicit definition *********************************/

                      /* Read data */

                      val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                      -INFTY, INFTY);
                      WDB[loc0 + TRANS_RX1] = val;

                      val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                      -INFTY, INFTY);
                      WDB[loc0 + TRANS_RX2] = val;

                      val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                      -INFTY, INFTY);
                      WDB[loc0 + TRANS_RX3] = val;

                      val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                      -INFTY, INFTY);
                      WDB[loc0 + TRANS_RX4] = val;

                      val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                      -INFTY, INFTY);
                      WDB[loc0 + TRANS_RX5] = val;

                      val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                      -INFTY, INFTY);
                      WDB[loc0 + TRANS_RX6] = val;

                      val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                      -INFTY, INFTY);
                      WDB[loc0 + TRANS_RX7] = val;

                      val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                      -INFTY, INFTY);
                      WDB[loc0 + TRANS_RX8] = val;

                      val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                      -INFTY, INFTY);
                      WDB[loc0 + TRANS_RX9] = val;

                      /* Set explicit flag */

                      WDB[loc0 + TRANS_EXPLI] = (double)YES;

                      /* Check if order is given */

                      if (j < np)
                        {
                          val = TestParam(pname, fname, line, params[j++],
                                          PTYPE_INT, 1, 2);
                          WDB[loc0 + TRANS_ORDER] = val - 1.0;
                        }

                      /***********************************************************/
                    }
                  else if ((np - j == 3) || (np - j == 4))
                    {
                      /***** Angles *********************************************/

                      /***********************************************************/

                      /* Read angles */

                      ax = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                     -360.0, 360.0);
                      ay = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                     -360.0, 360.0);
                      az = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                     -360.0, 360.0);

                      /* Convert to rad */

                      ax = ax*PI/180.0;
                      ay = ay*PI/180.0;
                      az = az*PI/180.0;

                      /* Fill matrix */

                      WDB[loc0 + TRANS_RX1] = 1.0;
                      WDB[loc0 + TRANS_RX2] = 0.0;
                      WDB[loc0 + TRANS_RX3] = 0.0;
                      WDB[loc0 + TRANS_RX4] = 0.0;
                      WDB[loc0 + TRANS_RX5] = cos(ax);
                      WDB[loc0 + TRANS_RX6] = -sin(ax);
                      WDB[loc0 + TRANS_RX7] = 0.0;
                      WDB[loc0 + TRANS_RX8] = sin(ax);
                      WDB[loc0 + TRANS_RX9] = cos(ax);

                      WDB[loc0 + TRANS_RY1] = cos(ay);
                      WDB[loc0 + TRANS_RY2] = 0.0;
                      WDB[loc0 + TRANS_RY3] = sin(ay);
                      WDB[loc0 + TRANS_RY4] = 0.0;
                      WDB[loc0 + TRANS_RY5] = 1.0;
                      WDB[loc0 + TRANS_RY6] = 0.0;
                      WDB[loc0 + TRANS_RY7] = -sin(ay);
                      WDB[loc0 + TRANS_RY8] = 0.0;
                      WDB[loc0 + TRANS_RY9] = cos(ay);

                      WDB[loc0 + TRANS_RZ1] = cos(az);
                      WDB[loc0 + TRANS_RZ2] = -sin(az);
                      WDB[loc0 + TRANS_RZ3] = 0.0;
                      WDB[loc0 + TRANS_RZ4] = sin(az);
                      WDB[loc0 + TRANS_RZ5] = cos(az);
                      WDB[loc0 + TRANS_RZ6] = 0.0;
                      WDB[loc0 + TRANS_RZ7] = 0.0;
                      WDB[loc0 + TRANS_RZ8] = 0.0;
                      WDB[loc0 + TRANS_RZ9] = 1.0;

                      /* Reset explicit flag */

                      WDB[loc0 + TRANS_EXPLI] = (double)NO;

                      /* T�� yhdist�� noi 3 erillist� matriisia yhdeksi     */
                      /* rotaatiomatriisiksi. Pit�isi teoriassa toimia,     */
                      /* mutta ei oo testattu kunnolla. Tolla if-lauseella  */
                      /* saa palautettua vanhan menetelm�n. Jos t�� toimii, */
                      /* niin ton expli-h�ss�k�n voinee poistaa sitten kun  */
                      /* on testattu (JLe / 21.2.2017 / 2.1.29) */

                      if (1 != 2)
                        {
                          WDB[loc0 + TRANS_RX1] = cos(ay)*cos(az);
                          WDB[loc0 + TRANS_RX2] = -cos(ay)*sin(az);
                          WDB[loc0 + TRANS_RX3] = sin(ay);
                          WDB[loc0 + TRANS_RX4] = cos(ax)*sin(az) +
                            cos(az)*sin(ax)*sin(ay);
                          WDB[loc0 + TRANS_RX5] = cos(ax)*cos(az) -
                            sin(ax)*sin(ay)*sin(az);
                          WDB[loc0 + TRANS_RX6] = -cos(ay)*sin(ax);
                          WDB[loc0 + TRANS_RX7] = sin(ax)*sin(az) -
                            cos(ax)*cos(az)*sin(ay);
                          WDB[loc0 + TRANS_RX8] = cos(az)*sin(ax) +
                            cos(ax)*sin(ay)*sin(az);
                          WDB[loc0 + TRANS_RX9] = cos(ax)*cos(ay);

                          WDB[loc0 + TRANS_RY1] = 0.0;
                          WDB[loc0 + TRANS_RY2] = 0.0;
                          WDB[loc0 + TRANS_RY3] = 0.0;
                          WDB[loc0 + TRANS_RY4] = 0.0;
                          WDB[loc0 + TRANS_RY5] = 0.0;
                          WDB[loc0 + TRANS_RY6] = 0.0;
                          WDB[loc0 + TRANS_RY7] = 0.0;
                          WDB[loc0 + TRANS_RY8] = 0.0;
                          WDB[loc0 + TRANS_RY9] = 0.0;

                          WDB[loc0 + TRANS_RZ1] = 0.0;
                          WDB[loc0 + TRANS_RZ2] = 0.0;
                          WDB[loc0 + TRANS_RZ3] = 0.0;
                          WDB[loc0 + TRANS_RZ4] = 0.0;
                          WDB[loc0 + TRANS_RZ5] = 0.0;
                          WDB[loc0 + TRANS_RZ6] = 0.0;
                          WDB[loc0 + TRANS_RZ7] = 0.0;
                          WDB[loc0 + TRANS_RZ8] = 0.0;
                          WDB[loc0 + TRANS_RZ9] = 0.0;

                          /* Set explicit flag */

                          WDB[loc0 + TRANS_EXPLI] = (double)YES;
                        }

                      /* Check if order is given */

                      if (j < np)
                        {
                          val = TestParam(pname, fname, line, params[j++],
                                          PTYPE_INT, 1, 2);
                          WDB[loc0 + TRANS_ORDER] = val - 1.0;
                        }
                    }
                  else
                    Error(loc0,"Number of parameters for rotation must be 3 or 9");
                }
              else
                WDB[loc0 + TRANS_ROT] = (double)NO;
            }
         }

      /***********************************************************************/

      /***** Geometry plot **************************************************/

      else if (!strcasecmp(word, "plot"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 3, 11, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_GPL0, GPL_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Init boundaries */

          WDB[loc0 + GPL_XMIN] = -INFTY;
          WDB[loc0 + GPL_XMAX] =  INFTY;
          WDB[loc0 + GPL_YMIN] = -INFTY;
          WDB[loc0 + GPL_YMAX] =  INFTY;
          WDB[loc0 + GPL_ZMIN] = -INFTY;
          WDB[loc0 + GPL_ZMAX] =  INFTY;
          WDB[loc0 + GPL_POS] =  -INFTY;

          /* Plot material boundaries by default */

          WDB[loc0 + GPL_PLOT_BOUND] = 2.0;

          /* Update counter and set index */

          WDB[DATA_N_GEOM_PLOTS] = RDB[DATA_N_GEOM_PLOTS] + 1.0;
          WDB[loc0 + GPL_IDX] = RDB[DATA_N_GEOM_PLOTS];

          /* Plot type */

          sprintf(str, "%s", params[j++]);

          /* Check if non-default type is given */

          if (strlen(str) > 1)
            {
              if (str[1] == '0')
                WDB[loc0 + GPL_PLOT_BOUND] = 0.0;
              else if (str[1] == '1')
                WDB[loc0 + GPL_PLOT_BOUND] = 1.0;
              else if (str[1] == '2')
                WDB[loc0 + GPL_PLOT_BOUND] = 2.0;
              else if (str[1] == '3')
                WDB[loc0 + GPL_PLOT_BOUND] = 3.0;
              else if ((str[1] == '4') || (str[1] == '5') || (str[1] == '6') ||
                       (str[1] == '7'))
                {
                  /* Plot material boundaries */

                  WDB[loc0 + GPL_PLOT_BOUND] = 2.0;

                  /* Linear or log scale */

                  if (str[1] == '4')
                    WDB[loc0 + GPL_IMP_SCALE] = 1.0;
                  else if (str[1] == '5')
                    WDB[loc0 + GPL_IMP_SCALE] = 2.0;
                  else if (str[1] == '6')
                    WDB[loc0 + GPL_IMP_SCALE] = 3.0;
                  else if (str[1] == '7')
                    WDB[loc0 + GPL_IMP_SCALE] = 4.0;
                  else
                    Error(loc0, "Invalid color scale %c", str[2]);

                  /* Read minimum and maximum */

                  WDB[loc0 + GPL_IMP_MIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              1E-200, 1E+200);

                  WDB[loc0 + GPL_IMP_MAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              RDB[loc0 + GPL_IMP_MIN], 1E+200);

                  /* Read energy (negative entries for photons) */

                  WDB[loc0 + GPL_IMP_E] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);
                }
            }

          /* Put size */

          if (j > np - 2)
            Error(loc0, "Missing image size");

          n = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              5, 100000);
          WDB[loc0 + GPL_PIX_X] = (double)n;

          n = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              5, 100000);
          WDB[loc0 + GPL_PIX_Y] = (double)n;

          /* Get position */

          if (j < np)
            {
              val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);
              WDB[loc0 + GPL_POS] = val;
            }

          /* Check mode */

          if ((str[0] == 'x') || (str[0] == 'X') || (str[0] == '1'))
            {
              /* YZ-plot */

              WDB[loc0 + GPL_GEOM] = (double)PLOT_MODE_YZ;

              /* Read boundaries */

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_YMIN] = val;
                }

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_YMAX] = val;
                }

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_ZMIN] = val;
                }

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_ZMAX] = val;
                }
            }

          else if ((str[0] == 'y') || (str[0] == 'Y') || (str[0] == '2'))
            {
              /* XZ-plot */

              WDB[loc0 + GPL_GEOM] = (double)PLOT_MODE_XZ;

              /* Read boundaries */

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_XMIN] = val;
                }

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_XMAX] = val;
                }

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_ZMIN] = val;
                }

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_ZMAX] = val;
                }
            }

          else if ((str[0] == 'z') || (str[0] == 'Z') || (str[0] == '3'))
            {
              /* XY-plot */

              WDB[loc0 + GPL_GEOM] = (double)PLOT_MODE_XY;

              /* Read boundaries */

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_XMIN] = val;
                }

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_XMAX] = val;
                }

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_YMIN] = val;
                }

              if (j < np)
                {
                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc0 + GPL_YMAX] = val;
                }
            }
          else
            Error(loc0, "Invalid plot type %s", str);
         }

      /***********************************************************************/

      /***** Mesh plot *******************************************************/

      else if (!strcasecmp(word, "mesh"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, 13, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_MPL0, MPL_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Init boundaries */

          WDB[loc0 + MPL_XMIN] = -INFTY;
          WDB[loc0 + MPL_XMAX] =  INFTY;
          WDB[loc0 + MPL_YMIN] = -INFTY;
          WDB[loc0 + MPL_YMAX] =  INFTY;
          WDB[loc0 + MPL_ZMIN] = -INFTY;
          WDB[loc0 + MPL_ZMAX] =  INFTY;

          /* Set default type and colormap */

          WDB[loc0 + MPL_TYPE] = MPL_TYPE_FLUXPOW;
          WDB[loc0 + MPL_COLMAP] = PALETTE_HOTCOLD;
          WDB[loc0 + MPL_COLOR_SCALE] = COLOR_SCALE_LIN;

          /* Reset number of distributions and divider flag */

          WDB[loc0 + MPL_NDIST] = 1.0;
          WDB[loc0 + MPL_DIV] = (double)NO;

          /* Reset minimum and maximum values */

          WDB[loc0 + MPL_MIN1] = -1.0;
          WDB[loc0 + MPL_MAX1] = -1.0;
          WDB[loc0 + MPL_MIN2] = -1.0;
          WDB[loc0 + MPL_MAX2] = -1.0;

          /* Read parameters */

          j = 0;

          /* Get type */

          type = (long)TestParam(pname, fname, line, params[j++],
                                 PTYPE_INT, 1, 14);

          /* Check type */

          if (type < 4)
            {
              /* Serpent 1 style flux / fission rate plot */

              WDB[loc0 + MPL_AX] = (double)type;
              WDB[loc0 + MPL_NDIST] = 2.0;
            }
          else if (type == 14)
            {
              /* Source point mesh no longer supported */

              Error(loc0, "Type 14, use detector instead");
            }
          else if (type == MPL_TYPE_FLUXTEMP)
            {
              /* Flux / temperature plot */

              WDB[loc0 + MPL_TYPE] = (double)type;
              WDB[loc0 + MPL_NDIST] = 2.0;
              WDB[loc0 + MPL_DIV] = (double)YES;

              /* Read axis */

              WDB[loc0 + MPL_AX] =
                TestParam(pname, fname, line, params[j++], PTYPE_INT, 1, 3);
            }
          else
            {
              /* Check for two-valued type */

              if (type == 4)
                WDB[loc0 + MPL_NDIST] = 2.0;

              /* Put type */

              WDB[loc0 + MPL_TYPE] = (double)type;

              /* Read color map */

              WDB[loc0 + MPL_COLMAP] =
                TestParam(pname, fname, line, params[j++], PTYPE_INT,-13, 13);

              /* Check log type */

              if ((long)RDB[loc0 + MPL_COLMAP] > 0)
                WDB[loc0 + MPL_COLOR_SCALE] = COLOR_SCALE_LIN;
              else if ((long)RDB[loc0 + MPL_COLMAP] < 0)
                {
                  WDB[loc0 + MPL_COLOR_SCALE] = COLOR_SCALE_LOG;
                  WDB[loc0 + MPL_COLMAP] = -RDB[loc0 + MPL_COLMAP];
                }
              else
                Error(loc0, "Invalid color map 0");

              /* Read detector name */

              if ((type == MPL_TYPE_DET) || (type == MPL_TYPE_DET_IMP))
                WDB[loc0 + MPL_PTR_DET] = (double)PutText(params[j++]);

              /* Put divider flag for some types */

              if ((type == MPL_TYPE_DET_IMP) || (type == MPL_TYPE_DENSITY) ||
                  (type == MPL_TYPE_DT_NEFF) || (type == MPL_TYPE_DT_GEFF))
                WDB[loc0 + MPL_DIV] = (double)YES;

              /* Read axis */

              WDB[loc0 + MPL_AX] =
                TestParam(pname, fname, line, params[j++], PTYPE_INT, 1, 4);
            }

            /* Read size */

          WDB[loc0 + MPL_NX] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT, 10, 5000);

          WDB[loc0 + MPL_NY] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT, 10, 5000);

          /* Get symmetry */

          if (j < np)
            WDB[loc0 + MPL_SYM] =
              TestParam(pname, fname, line, params[j++], PTYPE_INT, 0, 8);

          /* Check symmetry */

          if (((long)RDB[loc0 + MPL_SYM] != 0) &&
              ((long)RDB[loc0 + MPL_NX] != (long)RDB[loc0 + MPL_NY]))
            Error(loc0, "Mesh plot must be square if symmetry option is used");

          /* Get dimensions */

          if (j < np - 1)
            {
              WDB[loc0 + MPL_XMIN] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          -INFTY, INFTY);

              WDB[loc0 + MPL_XMAX] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          -INFTY, INFTY);
            }
          if (j < np - 1)
            {
              WDB[loc0 + MPL_YMIN] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          -INFTY, INFTY);

              WDB[loc0 + MPL_YMAX] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          -INFTY, INFTY);
            }
          if (j < np - 1)
            {
              WDB[loc0 + MPL_ZMIN] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          -INFTY, INFTY);

              WDB[loc0 + MPL_ZMAX] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          -INFTY, INFTY);
            }

          /* Check boundaries */

          if ((RDB[loc0 + MPL_XMAX] <= RDB[loc0 + MPL_XMIN]) ||
              (RDB[loc0 + MPL_YMAX] <= RDB[loc0 + MPL_YMIN]) ||
              (RDB[loc0 + MPL_ZMAX] <= RDB[loc0 + MPL_ZMIN]))
            Error(loc0, "Error in mesh boundaries");
         }

      /***********************************************************************/

      /***** Explicit stochastic (pebble bed) geometry ***********************/

      else if (!strcasecmp(word, "pbed"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 3, 4, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_PB0, PBED_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Read name */

          WDB[loc0 + PBED_PTR_NAME] = (double)PutText(params[j++]);

          /* Read background universe */

          WDB[loc0 + PBED_PTR_BG_UNIV] = (double)PutText(params[j++]);

          /* Read file name */

          WDB[loc0 + PBED_PTR_FNAME] = (double)PutText(params[j++]);

          /* Reset results flag */

          WDB[loc0 + PBED_CALC_RESULTS] = (double)NO;

          /* Check if results are requested */

          while (j < np)
            {
              if (!strcmp(params[j], "pow"))
                WDB[loc0 + PBED_CALC_RESULTS] = (double)YES;

              j++;
            }

          /* Add to counter */

          WDB[DATA_N_PBED] = RDB[DATA_N_PBED] + 1.0;
        }

      /***********************************************************************/

      /***** Implicit particle fuel model ************************************/

      else if (!strcasecmp(word, "disp"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 15, fname);

          Error(-1, pname, fname, line,
                "Implicit particle fuel model not supported in Serpent 2");
        }

      /***********************************************************************/

      /***** Unstructured mesh based geometry ********************************/

      else if (!strcasecmp(word, "umsh"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 15, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_UMSH0, UMSH_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Read name */

          if (j < np)
            WDB[loc0 + UMSH_PTR_NAME] = (double)PutText(params[j++]);
          else
            Error(loc0, "Missing universe name");

          /* Read background universe */

          if (j < np)
            WDB[loc0 + UMSH_PTR_BG_UNIV] = (double)PutText(params[j++]);
          else
            Error(loc0, "Missing background universe name");

          /* Create associated inteface structure */

          loc1 = NewItem(loc0 + UMSH_PTR_IFC, IFC_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc1 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc1 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc1 + PARAM_LINE] = (double)line;

          /* Read search mesh split limit */

          if (j < np)
            WDB[loc1 + IFC_SEARCH_MESH_ADA_SPLIT] =
              TestParam(pname, fname, line, params[j++],
                        PTYPE_INT, 1, 10000000000);
          else
            Error(loc0, "Missing search mesh split limit");

          /* Avoid compiler warning */

          ns = -1;

          /* Read search mesh dimensions */

          if (j < np)
            ns = (long)TestParam(pname, fname, line, params[j++],
                                 PTYPE_INT, 1, 10);
          else
            Error(loc0, "Missing search mesh dimension");

          /* Allocate memory */

          ptr = ReallocMem(DATA_ARRAY, ns + 1);
          WDB[loc1 + IFC_SEARCH_MESH_ADA_PTR_SZ] = (double)ptr;

          /* Read size */

          for (n = 0; n < ns; n++)
            {
              if (j < np)
                WDB[ptr++] = TestParam(pname, fname, line, params[j++],
                                       PTYPE_INT, 2, 10000);
              else
                Error(loc0, "Missing mesh size");
            }

          /* Put terminator */

          WDB[ptr] = -1.0;

          /* Read points file */

          if (j < np)
            WDB[loc1 + IFC_PTR_OF_PFILE] = (double)PutText(params[j++]);
          else
            Error(loc0, "Missing points file");

          /* Read faces file */

          if (j < np)
            WDB[loc1 + IFC_PTR_OF_FFILE] = (double)PutText(params[j++]);
          else
            Error(loc0, "Missing faces file");

          /* Read owner file */

          if (j < np)
            WDB[loc1 + IFC_PTR_OF_OFILE] = (double)PutText(params[j++]);
          else
            Error(loc0, "Missing owner file");

          /* Read neighbour file */

          if (j < np)
            WDB[loc1 + IFC_PTR_OF_NFILE] = (double)PutText(params[j++]);
          else
            Error(loc0, "Missing neighbour file");

          /* Read materials file */

          if (j < np)
            WDB[loc1 + IFC_PTR_OF_MFILE] = (double)PutText(params[j++]);
          else
            Error(loc0, "Missing materials file");
        }

      /***********************************************************************/

      /***** Multi-physics interface *****************************************/

      else if (!strcasecmp(word, "ifc"))
        {
          /* Update memory size */

          WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] +
            (double)MemCount();

          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 102, fname);

          /* Get memory size before creation of the interface */

          mem0 = RDB[DATA_ALLOC_MAIN_SIZE];

          /* Create new item */

          loc0 = NewItem(DATA_PTR_IFC0, IFC_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Read file name */

          WDB[loc0 + IFC_PTR_INPUT_FNAME] = (double)PutText(params[j++]);

          /* Read interface */
          /* Time binning is not yet set in dynamic simulations */

          ReadInterface(loc0, (long)NO);

          /* Read additional parameters */

          while (j < np)
            {
              if (!strcasecmp(params[j], "setmat"))
                {
                  /* Manually set materials that are linked to the interface */

                  j++;

                  /* Read number of materials */

                  if (j < np)
                    WDB[loc0 + IFC_N_MAT] = TestParam(pname, fname, line,
                                                      params[j++], PTYPE_INT, 1,
                                                      100);
                  else
                    Error(loc0, "Missing number of materials after \"setmat\"");

                  /* Get number of materials */

                  n = (long)RDB[loc0 + IFC_N_MAT];

                  /* Allocate memory for material array */

                  ptr = ReallocMem(DATA_ARRAY, n+1);

                  /* Store pointer to material array */

                  WDB[loc0 + IFC_PTR_MAT_ARR] = (double)ptr;

                  /* Read material names */

                  for (k = 0; k < n; k++)
                    {
                      if (j == np)
                        Error(loc0, "Could not read material name number %ld", k+1);

                      /* Read material name */

                      WDB[ptr++] = (double)PutText(params[j++]);
                    }

                  /* Null terminate the array */

                  WDB[ptr] = NULLPTR;
                }

              if (j < np)
                Die(FUNCTION_NAME, "Should read additional parameters\n");
            }

          /* Get memory size after creation of the interface */

          mem1 = RDB[DATA_ALLOC_MAIN_SIZE];

          /* Store ifc-memory size to interface (for MPI transfer) */

          WDB[loc0 + IFC_MEM_SIZE] = mem1 - mem0;

          /* Set flag */

          WDB[DATA_USE_DENSITY_FACTOR] = (double)YES;

          /* Update memory size */

          WDB[DATA_TOT_IFC_BYTES] = RDB[DATA_TOT_IFC_BYTES] +
            (double)MemCount();

        }

      /***********************************************************************/

      /***** General mesh based interface ************************************/

      else if (!strcasecmp(word, "dataifc"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 3, 12, fname);

          /* Get memory size before creation of the interface */

          mem0 = RDB[DATA_ALLOC_MAIN_SIZE];

          /* Create new item */

          loc0 = NewItem(DATA_PTR_DATAIFC0, DATAIFC_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Read mesh name */

          WDB[loc0 + DATAIFC_PTR_MESH] = (double)PutText(params[j++]);

          /* Read file name */

          WDB[loc0 + DATAIFC_PTR_INPUT_FNAME] = (double)PutText(params[j++]);

          /* Test file format */

          TestDOSFile(GetText(loc0 + DATAIFC_PTR_INPUT_FNAME));

          /* Open file for reading */

          if ((fp = fopen(GetText(loc0 + DATAIFC_PTR_INPUT_FNAME), "r")) == NULL)
            Error(loc0, "Data-interface file \"%s\" does not exist",
                  GetText(loc0 + DATAIFC_PTR_INPUT_FNAME));

          /* Close file after testing */

          fclose(fp);

          /* Read other data required for the type */

          if (!strcasecmp(params[j], "adens"))
            {
              j++;

              /* Store data type */

              WDB[loc0 + DATAIFC_DATA_TYPE] = DATA_TYPE_IFC_ADENS;

              /* Read and store nuclide name (ZAI) */

              WDB[loc0 + DATAIFC_PTR_NUCLIDE] = TestParam(pname, fname, line,
                                                          params[j++], PTYPE_INT, 10010, 993000);
            }
          else
            Error(loc0, "Unknown data type for data interface \"%s\"", params[j]);

          /* Material list */

          if (j < np)
            {
              /* Get number of materials */

              n = np - j;

              /* Allocate memory for material list */

              ptr = ReallocMem(DATA_ARRAY, n + 1);

              /* Store pointer to materials list */

              WDB[loc0 + DATAIFC_PTR_MAT_ARR] = (double)ptr;

              /* Read materials (linked in ProcessDataInterfaces) */

              for (; n > 0; n--)
                WDB[ptr++] = PutText(params[j++]);

              /* Put null pointer */

              WDB[ptr] = NULLPTR;
            }
        }

      /***********************************************************************/

      /***** FINIX coupling **************************************************/

      else if (!strcasecmp(word, "finrod"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 5, 7, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_FIN0, FINIX_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

#ifndef FINIX
          Error(loc0, "Source code not compiled with FINIX");
#endif
          /* Read data */

          j = 0;

          /* Read nest */

          WDB[loc0 + FINIX_PTR_UNI_NAME] = (double)PutText(params[j++]);

          /* Read rod type name */

          if (j == np)
            Error(loc0, "Missing rod type");
          else
            WDB[loc0 + FINIX_PTR_RODNAME] = (double)PutText(params[j++]);

          /* Read options name */

          if (j == np)
            Error(loc0, "Missing options");
          else
            WDB[loc0 + FINIX_PTR_OPTINAME] = (double)PutText(params[j++]);

          /* Read scenario name */

          if (j == np)
            Error(loc0, "Missing options");
          else
            WDB[loc0 + FINIX_PTR_SCENNAME] = (double)PutText(params[j++]);

          /* Read power binning name */

          if (j == np)
            Error(loc0, "Missing options");
          else
            WDB[loc0 + FINIX_PTR_POWMSH] = (double)PutText(params[j++]);

          /* Check if additional parameters are given */

          if (j < np)
            {
              if (!strcasecmp(params[j], "nrods"))
                {
                  /* Number of rods */
                  /* Power will be divided with this instead of */
                  /* nst + NEST_COUNT deduced by Serpent        */

                  j++;

                  if (j == np)
                    Error(loc0, "Missing number of rods after \"nrod\"-card");

                  /* Store number of rods */

                  WDB[loc0 + FINIX_N_RODS] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, 10000.0);

                }
              else
                {
                  /* Unknown parameter */

                  Error(loc0, "Unknown parameter \"%s\"", params[j]);
                }
            }

          /* Set ETTM mode */

          WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;

          /* Set coupled calculation flag */

          WDB[DATA_RUN_CC] = (double)YES;

        }
      else if (!strcasecmp(word, "finbin"))
        {
          /***** FINIX power binning mesh  ********************************/

          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 7, 7, fname);

          /* Allocate memory for mesh structure */

          loc0 = NewItem(DATA_PTR_FINBIN0, MESH_BLOCK_SIZE);

          /* Alloc value pair for storing index (only seek once per collision) */

          AllocValuePair(loc0 + MESH_PREV_COL_IDX);

          /* Read data */

          j = 0;

          /* Store name of power binning */

          WDB[loc0 + MESH_PTR_NAME] = (double)PutText(params[j++]);

          /* Cylindrical mesh */

          WDB[loc0 + MESH_TYPE] = (double)MESH_TYPE_CYLINDRICAL;

          /* Store flag to use lowest level coordinates in */
          /* scoring rather than universe 0 coordinates */

          WDB[loc0 + MESH_LOCAL_COORDS] = (double)YES;

          /* Check number of parameters */

          if (j > np - 6)
            Error(loc0, "Missing mesh parameters");

          /* Read parameters */

          /* Read radial binning */

          WDB[loc0 + MESH_MIN0] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      0.0, INFTY);
          WDB[loc0 + MESH_MAX0] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      RDB[loc0 + MESH_MIN0], INFTY);
          WDB[loc0 + MESH_N0] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT,
                      1, 1000000);

          /* Angular binning is fixed */

          WDB[loc0 + MESH_MIN1] = 0.0;
          WDB[loc0 + MESH_MAX1] = 2*PI;
          WDB[loc0 + MESH_N1] = 1;

          /* Read axial binning */

          WDB[loc0 + MESH_MIN2] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      -INFTY, INFTY);
          WDB[loc0 + MESH_MAX2] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      RDB[loc0 + MESH_MIN2], INFTY);
          WDB[loc0 + MESH_N2] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT,
                      1, 1000000);

          /***************************************************************/

        }

      /***** General data mesh definition  ********************************/

      else if (!strcasecmp(word, "datamesh"))
        {

          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 12, 12, fname);

          /* Allocate memory for mesh structure */

          loc0 = NewItem(DATA_PTR_MESH0, MESH_BLOCK_SIZE);

          /* Read data */

          j = 0;

          /* Store name of data mesh */

          WDB[loc0 + MESH_PTR_NAME] = (double)PutText(params[j++]);

          /* Store mesh type */

          WDB[loc0 + MESH_TYPE] = TestParam(pname, fname, line, params[j++],
                                            PTYPE_INT, MESH_TYPE_CARTESIAN,
                                            MESH_TYPE_ICYL);

          /* Store flag to use lowest level coordinates in */
          /* scoring rather than universe 0 coordinates */

          WDB[loc0 + MESH_LOCAL_COORDS] = TestParam(pname, fname, line, params[j++],
                                            PTYPE_LOGICAL);


          /* Read mesh parameters */

          /* Read binning in first direction */

          WDB[loc0 + MESH_MIN0] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      -INFTY, INFTY);
          WDB[loc0 + MESH_MAX0] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      RDB[loc0 + MESH_MIN0], INFTY);
          WDB[loc0 + MESH_N0] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT,
                      1, 1000000);

          /* Read binning in second direction */

          WDB[loc0 + MESH_MIN1] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      -INFTY, INFTY);
          WDB[loc0 + MESH_MAX1] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      RDB[loc0 + MESH_MIN1], INFTY);
          WDB[loc0 + MESH_N1] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT,
                      1, 1000000);

          /* Read binning in third direction */

          WDB[loc0 + MESH_MIN2] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      -INFTY, INFTY);
          WDB[loc0 + MESH_MAX2] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      RDB[loc0 + MESH_MIN2], INFTY);
          WDB[loc0 + MESH_N2] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT,
                      1, 1000000);

          /***************************************************************/

        }

      /***********************************************************************/

      /***** Lattice definition **********************************************/

      else if (!strcasecmp(word, "lat"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 6, MAX_LATTICE_ITEMS + 7,
                         fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_L0, LAT_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Read name */

          WDB[loc0 + LAT_PTR_NAME] = (double)PutText(params[j++]);

          /* Lattice type */

          type = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT, 1,
                                 LATTICE_TYPES);

          /* Check illegal type */

          if (type == LAT_TYPE_RND)
            Error(loc0, "Type %ld not allowed", type);

          /* Put type */

          WDB[loc0 + LAT_TYPE] = (double)type;

          /* Lattice origin (siirr� t�� noiden if-lausekkeiden sis��n?) */

          WDB[loc0 + LAT_ORIG_X0] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      -INFTY, INFTY);

          WDB[loc0 + LAT_ORIG_Y0] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      -INFTY, INFTY);

          if ((type == LAT_TYPE_CUBOID) || (type == LAT_TYPE_XPRISM) ||
              (type == LAT_TYPE_YPRISM))
            WDB[loc0 + LAT_ORIG_Z0] =
              TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                        -INFTY, INFTY);

          /* Reset pointer */

          ptr = -1;

          /* Check type */

          if ((type == LAT_TYPE_S) || (type ==LAT_TYPE_HX) ||
              (type == LAT_TYPE_HY) || (type == LAT_TYPE_XTRIAG) ||
              (type == LAT_TYPE_INFS) || (type == LAT_TYPE_INFHY) ||
              (type == LAT_TYPE_INFHX))
            {
              /***** Simple types ********************************************/

              /* Lattice size */

              if ((type == LAT_TYPE_S) || (type ==LAT_TYPE_HX) ||
                  (type == LAT_TYPE_HY) || (type == LAT_TYPE_XTRIAG))
                {
                  nx = (long)TestParam(pname, fname, line, params[j++],
                                       PTYPE_INT, 1, 200);

                  ny = (long)TestParam(pname, fname, line, params[j++],
                                       PTYPE_INT, 1, 200);
                }
              else
                {
                  nx = 1;
                  ny = 1;
                }

              /* Put values */

              WDB[loc0 + LAT_NX] = (double)nx;
              WDB[loc0 + LAT_NY] = (double)ny;
              WDB[loc0 + LAT_NTOT] = (double)(nx*ny);

              /* Lattice pitch */

              val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              ZERO, INFTY);

              WDB[loc0 + LAT_PITCH] = val;

              /* Check number of items */

              if (np - j != nx*ny)
                Error(loc0, "Invalid number of lattice elements");

              /* Allocate memory for data */

              ptr = ReallocMem(DATA_ARRAY, nx*ny + 1);

              /* Set pointer */

              WDB[loc0 + LAT_PTR_FILL] = (double)ptr;

              /* Read items */

              while (j < np)
                WDB[ptr++] = (double)PutText(params[j++]);

              /* Put null pointer */

              WDB[ptr] = NULLPTR;

              /***************************************************************/
            }
          else if (type == LAT_TYPE_ZSTACK)
            {
              /***** Vertical stack ******************************************/

              /* Number of layers */

              ns = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                   1, 1000);

              WDB[loc0 + LAT_NTOT] = (double)ns;

              /* Check number of items */

              if (np - j != 2*ns)
                Error(loc0, "Invalid number of layers");

              /* Allocate memory for data */

              loc1 = ReallocMem(DATA_ARRAY, ns + 1);
              loc2 = ReallocMem(DATA_ARRAY, ns);

              /* Set pointers */

              WDB[loc0 + LAT_PTR_FILL] = (double)loc1;
              WDB[loc0 + LAT_PTR_Z] = (double)loc2;

              /* Read data */

              for (n = 0; n < ns; n++)
                {
                  /* Z-coordinate */

                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                  WDB[loc2++] = val;

                  /* Universe */

                  WDB[loc1++] = (double)PutText(params[j++]);
                }

              /* Put null pointer */

              WDB[loc1] = NULLPTR;

              /***************************************************************/
            }
          else if (type == LAT_TYPE_CLU)
            {
              /***** Cluster type lattice ************************************/

              /* Read number of rings */

              nr = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                 1, 1000);

              WDB[loc0 + LAT_N_RINGS] = (double)nr;

              /* Loop over rings */

              for (n = 0; n < nr; n++)
                {
                  /* Create new ring */

                  loc1 = NewItem(loc0 + LAT_PTR_FILL, RING_BLOCK_SIZE);

                  /* Read number of sectors */

                  ns = (long)TestParam(pname, fname, line, params[j++],
                                       PTYPE_INT, 1, 1000);

                  WDB[loc1 + RING_N_SEC] = (double)ns;

                  /* Add to total number of elements */

                  WDB[loc0 + LAT_NTOT] = RDB[loc0 + LAT_NTOT] + (double)ns;

                  /* Read ring radius */

                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  0.0, INFTY);

                  WDB[loc1 + RING_RAD] = val;

                  /* Read tilt (convert to radians) */

                  val = TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -360.0, 360.0);

                  WDB[loc1 + RING_TILT] = PI*val/180;

                  /* Allocate memory for data */

                  ptr = ReallocMem(DATA_ARRAY, ns + 1);

                  /* Set pointer */

                  WDB[loc1 + RING_PTR_FILL] = (double)ptr;

                  /* Read items */

                  for (m = 0; m < ns; m++)
                    {
                      if (j == np)
                        Error(loc0, "Missing lattice elements");
                      else
                        WDB[ptr++] = (double)PutText(params[j++]);
                    }

                  /* Put null pointer */

                  WDB[ptr] = NULLPTR;

                  /* Compare indexes */

                  if (j > np)
                    Error(loc0, "Invalid number of lattice elements");
                }

              /* Compare indexes */

              if (j < np)
                Error(loc0, "Invalid number of lattice elements");

              /* Close ring list */

              loc1 = (long)RDB[loc0 + LAT_PTR_FILL];
              CloseList(loc1);

              /***************************************************************/
            }
          else if ((type == LAT_TYPE_CUBOID) || (type == LAT_TYPE_XPRISM) ||
                   (type == LAT_TYPE_YPRISM))
            {
              /***** 3D lattice **********************************************/

              /* Lattice size */

              nx = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                   1, 100);

              ny = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                   1, 100);

              nz = (long)TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                   1, 100);

              /* Put values */

              WDB[loc0 + LAT_NX] = (double)nx;
              WDB[loc0 + LAT_NY] = (double)ny;
              WDB[loc0 + LAT_NZ] = (double)nz;
              WDB[loc0 + LAT_NTOT] = (double)(nx*ny*nz);

              /* Lattice pitch */

              WDB[loc0 + LAT_PITCHX] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          ZERO, INFTY);

              WDB[loc0 + LAT_PITCHY] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          ZERO, INFTY);

              WDB[loc0 + LAT_PITCHZ] =
                TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                          ZERO, INFTY);

              /* Check number of items */

              if (np - j != nx*ny*nz)
                Error(loc0, "Invalid number of lattice elements");

              /* Allocate memory for data */

              ptr = ReallocMem(DATA_ARRAY, nx*ny*nz + 1);

              /* Set pointer */

              WDB[loc0 + LAT_PTR_FILL] = (double)ptr;

              /* Read items */

              while (j < np)
                WDB[ptr++] = (double)PutText(params[j++]);

              /* Put null pointer */

              WDB[ptr] = NULLPTR;

              /***************************************************************/
            }
         }

      /***********************************************************************/

      /****** Depletion history **********************************************/

      else if (!strcasecmp(word, "dep"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, 10000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_BURN_PTR_DEP, DEP_HIS_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Set burnup mode */

          WDB[DATA_BURNUP_CALCULATION_MODE] = YES;

          /* Read data */

          j = 0;

          if ((ne = np - j) < 2)
            Error(loc0, "Invalid number of parameters");

          /* Allocate memory for steps */

          loc1 = ReallocMem(DATA_ARRAY, ne);

          /* Set pointer */

          WDB[loc0 + DEP_HIS_PTR_STEPS] = (double)loc1;

          /* Reset number of steps */

          ns = 0;

          /* Loop over parameters */

          for (n = 0; n < ne; n++)
            {
              /* Check if depletion step */

              if (!strcmp(params[j], "bustep"))
                {
                  WDB[loc0 + DEP_HIS_STEP_TYPE] = DEP_STEP_BU_STEP;
                  j++;
                }
              else if (!strcmp(params[j], "daystep"))
                {
                  WDB[loc0 + DEP_HIS_STEP_TYPE] = DEP_STEP_DAY_STEP;
                  j++;
                }
              else if (!strcmp(params[j], "butot"))
                {
                  WDB[loc0 + DEP_HIS_STEP_TYPE] = DEP_STEP_BU_TOT;
                  j++;
                }
              else if (!strcmp(params[j], "daytot"))
                {
                  WDB[loc0 + DEP_HIS_STEP_TYPE] = DEP_STEP_DAY_TOT;
                  j++;
                }
              else if (!strcmp(params[j], "decstep"))
                {
                  WDB[loc0 + DEP_HIS_STEP_TYPE] = DEP_STEP_DEC_STEP;
                  j++;
                }
              else if (!strcmp(params[j], "dectot"))
                {
                  WDB[loc0 + DEP_HIS_STEP_TYPE] = DEP_STEP_DEC_TOT;
                  j++;
                }
              else if (!strcmp(params[j], "actstep"))
                {
                  WDB[loc0 + DEP_HIS_STEP_TYPE] = DEP_STEP_ACT_STEP;
                  j++;
                }
              else if (!strcmp(params[j], "acttot"))
                {
                  WDB[loc0 + DEP_HIS_STEP_TYPE] = DEP_STEP_ACT_TOT;
                  j++;
                }

              /* Check if pointer to reprocessor */

              else if (!strcmp(params[j], "pro"))
                {
                  WDB[loc0 + DEP_HIS_PTR_REPROC] =
                    (double)PutText(params[j + 1]);
                  j = j + 2;
                  n++;
                }

              /* Check if pointer to branch */

              else if (!strcmp(params[j], "bra"))
                {
                  WDB[loc0 + DEP_HIS_PTR_BRANCH] =
                    (double)PutText(params[j + 1]);
                  j = j + 2;
                  n++;
                }

              /* Should be a numberical entry */

              else
                {
                  WDB[loc1++] = TestParam(pname, fname, line, params[j++],
                                          PTYPE_REAL, ZERO, INFTY);
                  ns++;
                }
            }

          /* Set number of steps */

          if (ns == 0)
            Error(loc0, "No depletion steps provided");
          else
            WDB[loc0 + DEP_HIS_N_STEPS] = (double)ns;

          /* Add to total */

          WDB[DATA_BURN_TOT_STEPS] = RDB[DATA_BURN_TOT_STEPS] + (double)ns;
        }

      /***********************************************************************/

      /****** Depletion zone dividers ****************************************/

      else if (!strcasecmp(word, "div"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 3, 1000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_DIV0, DIV_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Init values */

          WDB[loc0 + DIV_NX] = 1.0;
          WDB[loc0 + DIV_XMIN] = -INFTY;
          WDB[loc0 + DIV_XMAX] = INFTY;
          WDB[loc0 + DIV_NY] = 1.0;
          WDB[loc0 + DIV_YMIN] = -INFTY;
          WDB[loc0 + DIV_YMAX] = INFTY;
          WDB[loc0 + DIV_NZ] = 1.0;
          WDB[loc0 + DIV_ZMIN] = -INFTY;
          WDB[loc0 + DIV_ZMAX] = INFTY;
          WDB[loc0 + DIV_NRAD] = 1.0;
          WDB[loc0 + DIV_RMIN] = 0.0;
          WDB[loc0 + DIV_RMAX] = INFTY;
          WDB[loc0 + DIV_NSEG] = 1.0;
          WDB[loc0 + DIV_SEG0] = 0.0;
          WDB[loc0 + DIV_SEP] = (double)NO;
          WDB[loc0 + DIV_LIMS_CHECK] = (double)YES;

          /* Get material pointer */

          WDB[loc0 + DIV_PTR_MAT] = (double)PutText(params[j++]);

          /* Loop over remaining parameters */

          while (j < np)
            {
              /* Check parameter type */

              strcpy(str, params[j++]);

              if (!strcmp(str, "sep"))
                {
                  /* Flag for handling each material zone separately */

                  if (j == np)
                    Error(loc0, "Missing value after \"%s\"", str);

                  /* Set level */

                  WDB[loc0 + DIV_SEP_LVL] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT, 0,
                              10000) - 1.0;

                  /* Set flag */

                  if ((long)RDB[loc0 + DIV_SEP_LVL] > -1)
                    WDB[loc0 + DIV_SEP] = (double)YES;
                  else
                    WDB[loc0 + DIV_SEP] = (double)NO;
                }
              else if (!strcmp(str, "lims"))
                {
                  /* Flag for mapping regions outside limits to  */
                  /* separate material (rather than givin error) */

                  if (j == np)
                    Error(loc0, "Missing value after \"%s\"", str);

                  /* Set level */

                  WDB[loc0 + DIV_LIMS_CHECK] =
                    TestParam(pname, fname, line, params[j++], PTYPE_LOGICAL);
                }
              else if (!strcmp(str, "subx"))
                {
                  /* x-subdivision */

                  if (j > np - 3)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Number of zones */

                  WDB[loc0 + DIV_NX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 10000);

                  /* Minimum */

                  WDB[loc0 + DIV_XMIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  /* Maximum */

                  WDB[loc0 + DIV_XMAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);
                }
              else if (!strcmp(str, "suby"))
                {
                  /* y-subdivision */

                  if (j > np - 3)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Number of yones */

                  WDB[loc0 + DIV_NY] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 10000);

                  /* Minimum */

                  WDB[loc0 + DIV_YMIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  /* Maximum */

                  WDB[loc0 + DIV_YMAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);
                }
              else if (!strcmp(str, "subz"))
                {
                  /* z-subdivision */

                  if (j > np - 3)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Number of zones */

                  WDB[loc0 + DIV_NZ] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 10000);

                  /* Minimum */

                  WDB[loc0 + DIV_ZMIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  /* Maximum */

                  WDB[loc0 + DIV_ZMAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);
                }
              else if (!strcmp(str, "subr"))
                {
                  /* Radial subdivision */

                  if (j > np - 3)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Number of zones */

                  WDB[loc0 + DIV_NRAD] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              -10000, 10000);

                  /* Check zero */

                  if ((long)RDB[loc0 + DIV_NRAD] == 0.0)
                    Error(loc0,
                          "Number of radial zones must be greater than zero");

                  /* Minimum */

                  WDB[loc0 + DIV_RMIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  /* Maximum */

                  WDB[loc0 + DIV_RMAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);
                }
              else if (!strcmp(str, "subs"))
                {
                  /* Angular segments */

                  if (j > np - 2)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Number of zones */

                  WDB[loc0 + DIV_NSEG] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 10000);

                  /* Zero angle */

                  WDB[loc0 + DIV_SEG0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, 360.0);
                }
              else if (!strcmp(str, "out"))
                {
                  /* Flag for output printing */

                  if (j == np)
                    Error(loc0, "Missing value after \"%s\"", str);

                  /* Set value */

                  WDB[loc0 + DIV_OUTPUT_FLAG] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT, 1,3);
                }
              else
                Error(loc0, "Invalid divider parameter \"%s\"", str);
            }
        }

      /***********************************************************************/

      /***** External source definition **************************************/

      else if (!strcasecmp(word, "src"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 10000000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_SRC0, SRC_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Reset weight */

          WDB[loc0 + SRC_WGT] = 1.0;

          /* Reset boundaries */

          WDB[loc0 + SRC_XMIN] = -INFTY;
          WDB[loc0 + SRC_XMAX] =  INFTY;
          WDB[loc0 + SRC_YMIN] = -INFTY;
          WDB[loc0 + SRC_YMAX] =  INFTY;
          WDB[loc0 + SRC_ZMIN] = -INFTY;
          WDB[loc0 + SRC_ZMAX] =  INFTY;
          WDB[loc0 + SRC_RMIN] =  INFTY;
          WDB[loc0 + SRC_RMAX] = -INFTY;

          /* Reset point */

          WDB[loc0 + SRC_X0] = -INFTY;
          WDB[loc0 + SRC_Y0] = -INFTY;
          WDB[loc0 + SRC_Z0] = -INFTY;

          /* Reset energy */

          WDB[loc0 + SRC_E] = -INFTY;

          /* Reset type */

          WDB[loc0 + SRC_TYPE] = -1.0;

          /* Read data */

          j = 0;

          /* Source name */

          WDB[loc0 + SRC_PTR_NAME] = (double)PutText(params[j++]);

          /* Check if type is given */

          if (j < np)
            {
              if (!strcmp(params[j], "neutron") || !strcmp(params[j], "n") ||
                  (atoi(params[j]) == PARTICLE_TYPE_NEUTRON))
                {
                  WDB[loc0 + SRC_TYPE] = (double)PARTICLE_TYPE_NEUTRON;
                  j++;
                }
              else if (!strcmp(params[j], "photon") ||
                       !strcmp(params[j], "gamma") ||
                       !strcmp(params[j], "p") ||
                       !strcmp(params[j], "g") ||
                       (atoi(params[j]) == PARTICLE_TYPE_GAMMA))
                {
                  WDB[loc0 + SRC_TYPE] = (double)PARTICLE_TYPE_GAMMA;
                  j++;
                }
            }

          /* Default type is neutron */

          if ((long)RDB[loc0 + SRC_TYPE] == -1)
            WDB[loc0 + SRC_TYPE] = (double)PARTICLE_TYPE_NEUTRON;

          /* Check particle type and set transport flag */

          if ((long)RDB[loc0 + SRC_TYPE] == PARTICLE_TYPE_NEUTRON)
            WDB[DATA_NEUTRON_TRANSPORT_MODE] = (double)YES;
          else if ((long)RDB[loc0 + SRC_TYPE] == PARTICLE_TYPE_GAMMA)
            WDB[DATA_PHOTON_TRANSPORT_MODE] = (double)YES;

          /* Loop over remaining parameters */

          while (j < np)
            {
              /* Check parameter type */

              strcpy(str, params[j++]);

              if (!strcmp(str, "se"))
                {
                  /* Energy. Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing value after \"%s\"", str);

                  /* Set value */

                  WDB[loc0 + SRC_E] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              ZERO, INFTY);
                }
              else if (!strcmp(str, "sw"))
                {
                  /* Weight. Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing value after \"%s\"", str);

                  /* Set value */

                  WDB[loc0 + SRC_WGT] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);
                }
              else if (!strcmp(str, "sm"))
                {
                  /* Material. Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing material name after \"%s\"", str);

                  /* Set value */

                  WDB[loc0 + SRC_PTR_MAT] = (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "sg"))
                {
                  /* Material for radioactive decay source */

                  if (j == np)
                    Error(loc0, "Missing material name after \"%s\"", str);
                  else
                    WDB[loc0 + SRC_PTR_RAD_SRC_MAT] =
                      (double)PutText(params[j++]);

                  /* Read mode */

                  if (j == np)
                    Error(loc0, "Missing sampling mode after \"%s\"", str);
                  else
                    WDB[loc0 + SRC_RAD_SRC_MODE] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_INT, 1, 2);

                  /* Set flag */

                  WDB[DATA_USE_DECAY_SRC] = (double)YES;
                }
              else if (!strcmp(str, "sc"))
                {
                  /* Cell. Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing cell after \"%s\"", str);

                  /* Set value */

                  WDB[loc0 + SRC_PTR_CELL] = (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "su"))
                {
                  /* Universe. Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing universe number after \"%s\"", str);

                  /* Set value */

                  WDB[loc0 + SRC_PTR_UNIV] = (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "ss"))
                {
                  /* Surface. Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing surface after \"%s\"", str);

                  /* Set value */

                  WDB[loc0 + SRC_PTR_SURF] = (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "sb"))
                {
                  /* Energy bins. Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing number of bins \"%s\"", str);

                  /* Get number of bins */

                  ne = (long)TestParam(pname, fname, line, params[j++],
                                       PTYPE_INT, 2, 1000000);

                  /* Allocate memory for structure */

                  loc1 = NewItem(loc0 + SRC_PTR_EBINS,
                                 SRC_EBIN_BLOCK_SIZE);

                  /* Put number of points */

                  WDB[loc1 + SRC_EBIN_NE] = (double)ne;

                  /* Read interpolation */

                  WDB[loc1 + SRC_EBIN_INTT] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              0, 4);

                  /* Allocate memory for data */

                  loc2 = ReallocMem(DATA_ARRAY, ne);
                  WDB[loc1 + SRC_EBIN_PTR_E] = (double)loc2;

                  loc3 = ReallocMem(DATA_ARRAY, ne);
                  WDB[loc1 + SRC_EBIN_PTR_PDF] = (double)loc3;

                  /* Loop over values */

                  for (m = 0; m < ne; m++)
                    {
                      /* Read energy */

                      if (j == np)
                        Error(loc0, "Missing energy");
                      else
                        WDB[loc2++] =
                          TestParam(pname, fname, line, params[j++],
                                    PTYPE_REAL, 0.0, 1000.0);

                      /* Read weight */

                      if (j == np)
                        Error(loc0, "Missing probability");
                      else
                        WDB[loc3++] =
                          TestParam(pname, fname, line, params[j++],
                                    PTYPE_REAL, 0.0, INFTY);
                    }
                }
              else if (!strcmp(str, "sx"))
                {
                  /* X-boundaries. Check number of parameters. */

                  if (j == np - 1)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Set values */

                  WDB[loc0 + SRC_XMIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  WDB[loc0 + SRC_XMAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  /* Swap boundaries */

                  if (RDB[loc0 + SRC_XMAX] < RDB[loc0 + SRC_XMIN])
                    {
                      val = RDB[loc0 + SRC_XMAX];
                      WDB[loc0 + SRC_XMAX] = RDB[loc0 + SRC_XMIN];
                      WDB[loc0 + SRC_XMIN] = val;
                    }
                }
              else if (!strcmp(str, "sy"))
                {
                  /* Y-boundaries. Check number of parameters. */

                  if (j == np - 1)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Set values */

                  WDB[loc0 + SRC_YMIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  WDB[loc0 + SRC_YMAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  /* Swap boundaries */

                  if (RDB[loc0 + SRC_YMAX] < RDB[loc0 + SRC_YMIN])
                    {
                      val = RDB[loc0 + SRC_YMAX];
                      WDB[loc0 + SRC_YMAX] = RDB[loc0 + SRC_YMIN];
                      WDB[loc0 + SRC_YMIN] = val;
                    }
                }
              else if (!strcmp(str, "sz"))
                {
                  /* Z-boundaries. Check number of parameters. */

                  if (j == np - 1)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Set values */

                  WDB[loc0 + SRC_ZMIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  WDB[loc0 + SRC_ZMAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  /* Swap boundaries */

                  if (RDB[loc0 + SRC_ZMAX] < RDB[loc0 + SRC_ZMIN])
                    {
                      val = RDB[loc0 + SRC_ZMAX];
                      WDB[loc0 + SRC_ZMAX] = RDB[loc0 + SRC_ZMIN];
                      WDB[loc0 + SRC_ZMIN] = val;
                    }
                }
              else if (!strcmp(str, "srad"))
                {
                  /* radial boundaries. Check number of parameters. */

                  if (j == np - 1)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Set values */

                  WDB[loc0 + SRC_RMIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);

                  WDB[loc0 + SRC_RMAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              RDB[loc0 + SRC_RMIN], INFTY);
                }
              else if (!strcmp(str, "st"))
                {
                  /* Time boundaries. Check number of parameters. */

                  if (j == np - 1)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Set values */

                  WDB[loc0 + SRC_TMIN] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);

                  WDB[loc0 + SRC_TMAX] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);

                  /* Swap boundaries */

                  if (RDB[loc0 + SRC_TMAX] < RDB[loc0 + SRC_TMIN])
                    {
                      val = RDB[loc0 + SRC_TMAX];
                      WDB[loc0 + SRC_TMAX] = RDB[loc0 + SRC_TMIN];
                      WDB[loc0 + SRC_TMIN] = val;
                    }
                }
              else if (!strcmp(str, "sp"))
                {
                  /* Point. Check number of parameters. */

                  if (j > np - 3)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Set values */

                  WDB[loc0 + SRC_X0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  WDB[loc0 + SRC_Y0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  WDB[loc0 + SRC_Z0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);
                }
              else if (!strcmp(str, "sd"))
                {
                  /* Direction. Check number of parameters. */

                  if (j > np - 3)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Set values */

                  WDB[loc0 + SRC_U0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -1.0, 1.0);

                  WDB[loc0 + SRC_V0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -1.0, 1.0);

                  WDB[loc0 + SRC_W0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -1.0, 1.0);
                }
              else if (!strcmp(str, "sr"))
                {
                  /* Reaction. Check number of parameters. */

                  if (j == np - 1)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Set values */

                  WDB[loc0 + SRC_PTR_XSDATA] = (double)PutText(params[j++]);
                  WDB[loc0 + SRC_PTR_REA] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              2, 1004);
                }
              else if (!strcmp(str, "sf"))
                {
                  /* Source file. Check number of parameters. */

                  if (j == np - 1)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Set file name */

                  WDB[loc0 + SRC_READ_PTR_FILE] = (double)PutText(params[j++]);

                  /* Get type */

                  type = (long)TestParam(pname, fname, line, params[j++],
                                         PTYPE_INT,-1, 3);

                  /* Check simple binary format */

                  if (type == 0)
                    Error(loc0, "Invalid source file type");
                  else if (type == -1)
                    {
                      WDB[loc0 + SRC_READ_BINARY] = (double)YES;
                      WDB[loc0 + SRC_READ_FILE_TYPE] = 1.0;
                    }
                  else
                    WDB[loc0 + SRC_READ_FILE_TYPE] = (double)type;
                }
              else if (!strcmp(str, "si"))
                {
                  /* User-defined routine. Check number of parameters */

                  if (j == np)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Get number of values */

                  ne = (long)TestParam(pname, fname, line, params[j++],
                                       PTYPE_INT, 0, 1000000);

                  /* Check */

                  if (ne > np - j)
                    Error(loc0, "Missing values after \"%s\"", str);

                  /* Allocate memory for structure */

                  loc1 = NewItem(loc0 + SRC_PTR_USR, SRC_USR_BLOCK_SIZE);

                  /* Set number of values */

                  WDB[loc1 + SRC_USR_NP] = (double)ne;

                  /* Allocate memory for values */

                  if (ne > 0)
                    ptr = ReallocMem(DATA_ARRAY, ne);
                  else
                    ptr = -1;

                  /* Put pointer */

                  WDB[loc1 + SRC_USR_PTR_PARAMS] = (double)ptr;

                  /* Read values */

                  for (m = 0; m < ne; m++)
                    WDB[ptr++] =
                      TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                -INFTY, INFTY);
                }
              else
                Error(loc0, "Invalid source parameter \"%s\"", str);
            }
         }

      /***********************************************************************/

      /***** Reprocessor *****************************************************/

      else if (!strcasecmp(word, "rep"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 10000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_REP0, REPROC_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* ID */

          WDB[loc0 + REPROC_PTR_NAME] = (double)PutText(params[j++]);

          /* Loop over remaining parameters */

          while (j < np)
            {
              /* Get parameter */

              strcpy(str, params[j++]);

              /* Check type */

              if (!strcmp(str, "ru"))
                {
                  /* Swap two universes, allocate memory for structure */

                  loc1 = NewItem(loc0 + REPROC_PTR_SWAP_LIST,
                                 REPROC_SWAP_BLOCK_SIZE);

                  /* Read first universe */

                  if (j == np)
                    Error(loc0, "Missing universe name after \"%s\"", str);
                  else
                    WDB[loc1 + REPROC_SWAP_PTR_UNI1] =
                      (double)PutText(params[j++]);

                  /* Read second universe */

                  if (j == np)
                    Error(loc0, "Missing universe name after \"%s\"", str);
                  else
                    WDB[loc1 + REPROC_SWAP_PTR_UNI2] =
                      (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "rm"))
                {
                  /* Remove material, allocate memory for structure */

                  loc1 = NewItem(loc0 + REPROC_PTR_REM_LIST,
                                 REPROC_REM_BLOCK_SIZE);

                  /* Read first material */

                  if (j == np)
                    Error(loc0, "Missing material name after \"%s\"", str);
                  else
                    WDB[loc1 + REPROC_REM_PTR_MAT1] =
                      (double)PutText(params[j++]);

                  /* Read second material */

                  if (j == np)
                    Error(loc0, "Missing material name after \"%s\"", str);
                  else
                    WDB[loc1 + REPROC_REM_PTR_MAT2] =
                      (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "rc"))
                {
                  /* Continuous reprocessing, allocate memory for structure */

                  loc1 = NewItem(loc0 + REPROC_PTR_CON_LIST,
                                 REPROC_CON_BLOCK_SIZE);

                  /* Read first material */

                  if (j == np)
                    Error(loc0, "Missing material name after \"%s\"", str);
                  else
                    WDB[loc1 + REPROC_CON_PTR_MAT1] =
                      (double)PutText(params[j++]);

                  /* Read second material */

                  if (j == np)
                    Error(loc0, "Missing material name after \"%s\"", str);
                  else
                    WDB[loc1 + REPROC_CON_PTR_MAT2] =
                      (double)PutText(params[j++]);

                  /* Read mass flow */

                  if (j == np)
                    Error(loc0, "Missing mass flow after \"%s\"", str);
                  else
                    WDB[loc1 + REPROC_CON_PTR_MFLOW] =
                      (double)PutText(params[j++]);

                  /* Read mode */

                  if (j == np)
                    Error(loc0, "Missing mode after \"%s\"", str);
                  else
                    WDB[loc1 + REPROC_CON_MODE] =
                      TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                0, 2);
                }
              else
                Error(loc0, "Invalid reprocessing parameter \"%s\"", str);
            }
         }

      /***********************************************************************/

      /***** Detector definition *********************************************/

      else if (!strcasecmp(word, "det"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 100000000, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_DET0, DET_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Reset volume */

          WDB[loc0 + DET_VOL] = -1.0;

          /* Reset mesh pointer */

          WDB[loc0 + DET_PTR_MESH] = NULLPTR;

          /* Read data */

          j = 0;

          /* Detector name */

          WDB[loc0 + DET_PTR_NAME] = (double)PutText(params[j++]);

          /* Check if type is given */

          if (j < np)
            {
              if (!strcmp(params[j], "neutron") || !strcmp(params[j], "n") ||
                  (atoi(params[j]) == PARTICLE_TYPE_NEUTRON))
                {
                  WDB[loc0 + DET_PARTICLE] = (double)PARTICLE_TYPE_NEUTRON;
                  j++;
                }
              else if (!strcmp(params[j], "photon") ||
                       !strcmp(params[j], "gamma") ||
                       !strcmp(params[j], "p") ||
                       !strcmp(params[j], "g") ||
                       (atoi(params[j]) == PARTICLE_TYPE_GAMMA))
                {
                  WDB[loc0 + DET_PARTICLE] = (double)PARTICLE_TYPE_GAMMA;
                  j++;
                }
            }

          /* Loop over remaining parameters */

          while (j < np)
            {
              /* Check parameter type */

              strcpy(str, params[j++]);

              if (!strcmp(str, "de"))
                {
                  /* Detector energy grid. Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing grid name after \"%s\"", str);

                  /* Set value */

                  WDB[loc0 + DET_PTR_EGRID] = (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "di"))
                {
                  /* Detector time binning Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing time bin name after \"%s\"", str);

                  /* Set value */

                  WDB[loc0 + DET_PTR_TME] = (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "dfl"))
                {
                  /* Detector flagging Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing target detector after \"%s\"", str);

                  /* Allocate memory */

                  loc1 = NewItem(loc0 + DET_PTR_FLAGGING, DET_FBIN_BLOCK_SIZE);

                  /* Get flag number */

                  if (j == np)
                    Error(loc0, "Missing flag number after \"%s\"", str);
                  else
                    WDB[loc1 + DET_FBIN_FLAG_NUMBER] =
                      TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                1, 64);

                  /* Get flag option */

                  if (j == np)
                    Error(loc0, "Missing flag option after \"%s\"", str);
                  else
                    WDB[loc1 + DET_FBIN_FLAG_OPTION] =
                      TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                -3, 3);

                  /* Check AND-logic */

                  if ((long)RDB[loc1 + DET_FBIN_FLAG_OPTION] < 0)
                    {
                      /* Set option */

                      WDB[loc1 + DET_FBIN_FLAG_AND_LOGIC] = (double)YES;

                      /* Adjust */

                      WDB[loc1 + DET_FBIN_FLAG_OPTION] =
                        -RDB[loc1 + DET_FBIN_FLAG_OPTION];
                    }
                  else
                    {
                      /* Reset option */

                      WDB[loc1 + DET_FBIN_FLAG_AND_LOGIC] = (double)NO;
                    }
                }
              else if (!strcmp(str, "dv"))
                {
                  /* Detector volume */

                  if (j == np)
                    Error(loc0, "Missing volume after \"%s\"", str);

                  WDB[loc0 + DET_VOL] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);
                }
              else if (!strcmp(str, "dfet"))
                {
                  if (j == np)
                    Error(loc0, "Missing FET type after \"%s\"", str);

                  /* FET detector! Create the FET parameters struct */

                  loc1 = NewItem(loc0 + DET_FET_PTR_PARAMS, FET_PARAMS_SIZE);

                  /* Read type */

                  WDB[loc1 + FET_PARAM_TYPE] =
                    TestParam(pname, fname, line, params[j++],
                              PTYPE_INT, FET_TYPES_START, FET_TYPES_END);

                  /* Read parameters */

                  switch (type = (long)RDB[loc1 + FET_PARAM_TYPE])
                    {
                    case FET_TYPE_CARTESIAN:
                      {
                        if (j == np)
                          Error(loc0, "Missing minimum 'x' dimension after \"%s\"", str);
                        WDB[loc1 + FET_CART_MIN_X] =
                          TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                    -INFTY, INFTY);

                        if (j == np)
                          Error(loc0, "Missing maximum 'x' dimension after \"%s\"", str);
                        WDB[loc1 + FET_CART_MAX_X] =
                          TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                    RDB[loc1 + FET_CART_MIN_X], INFTY);

                        if (j == np)
                          Error(loc0, "Missing polynomial order in dimension 'x' after \"%s\"", str);
                        WDB[loc1 + FET_CART_ORDER_X] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    -1, 100);

                        if (j == np)
                          Error(loc0, "Missing minimum 'y' dimension after \"%s\"", str);
                        WDB[loc1 + FET_CART_MIN_Y] =
                          TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                    -INFTY, INFTY);

                        if (j == np)
                          Error(loc0, "Missing maximum 'y' dimension after \"%s\"", str);
                        WDB[loc1 + FET_CART_MAX_Y] =
                          TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                    RDB[loc1 + FET_CART_MIN_Y], INFTY);

                        if (j == np)
                          Error(loc0, "Missing polynomial order in dimension 'y' after \"%s\"", str);
                        WDB[loc1 + FET_CART_ORDER_Y] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    -1, 100);

                        if (j == np)
                          Error(loc0, "Missing minimum 'z' dimension after \"%s\"", str);
                        WDB[loc1 + FET_CART_MIN_Z] =
                          TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                    -INFTY, INFTY);

                        if (j == np)
                          Error(loc0, "Missing maximum 'z' dimension after \"%s\"", str);
                        WDB[loc1 + FET_CART_MAX_Z] =
                          TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                    RDB[loc1 + FET_CART_MIN_Z], INFTY);

                        if (j == np)
                          Error(loc0, "Missing polynomial order in dimension 'z' after \"%s\"", str);
                        WDB[loc1 + FET_CART_ORDER_Z] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    -1, 100);
                        break;
                      }

                    case FET_TYPE_CYLINDRICAL:
                      {
                        if (j == np)
                          Error(loc0, "Missing radius 'r' after \"%s\"", str);
                        WDB[loc1 + FET_CYL_MAX_R] =
                          TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                    ZERO, INFTY);

                        if (j == np)
                          Error(loc0, "Missing Zernike polynomial order after \"%s\"", str);
                        WDB[loc1 + FET_CYL_ORDER_R] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    -1, 100);

                        if (j == np)
                          Error(loc0, "Missing axial start position after \"%s\"", str);
                        WDB[loc1 + FET_CYL_MIN_A] =
                          TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                    -INFTY, INFTY);

                        if (j == np)
                          Error(loc0, "Missing axial end position after \"%s\"", str);
                        WDB[loc1 + FET_CYL_MAX_A] =
                          TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                    RDB[loc1 + FET_CYL_MIN_A], INFTY);

                        if (j == np)
                          Error(loc0, "Missing axial polynomial order after \"%s\"", str);
                        WDB[loc1 + FET_CYL_ORDER_A] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    -1, 100);

                        if (j == np)
                          {
                            Warn(FUNCTION_NAME, "No orientation specified for cylindrical DET detector \"%s\""
                                 "\nAssuming axial orientation in x centered at (y, z) = (0, 0)", str);
                            WDB[loc1 + FET_CYL_ORIENTATION_A] = FET_ORIENTATION_X;
                            WDB[loc1 + FET_CYL_CENTER_A0] = 0;
                            WDB[loc1 + FET_CYL_CENTER_A1] = 0;
                          }
                        else
                          {
                            WDB[loc1 + FET_CYL_ORIENTATION_A] =
                              TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                        FET_ORIENTATIONS_START, FET_ORIENTATIONS_END);

                            switch ((long)RDB[loc1 + FET_CYL_ORIENTATION_A])
                              {
                              case FET_ORIENTATION_X:
                                {
                                  word[0] = 'y';
                                  word[1] = 'z';
                                  break;
                                }
                              case FET_ORIENTATION_Y:
                                {
                                  word[0] = 'x';
                                  word[1] = 'z';
                                  break;
                                }
                              case FET_ORIENTATION_Z:
                                {
                                  word[0] = 'x';
                                  word[1] = 'y';
                                  break;
                                }
#ifdef DEBUG /* We should never get here because the TestParam() above should have caught this */
                              default:
                                {
                                  Die(FUNCTION_NAME, "Unsupported orientation %ld", (long)RDB[loc1 + FET_CYL_ORIENTATION_A]);
                                }
#endif /* DEBUG */
                              }

                            if (j == np)
                              Error(loc0, "Missing \"%c\" coordinate of cylinder center after \"%s\"", word[0], str);
                            WDB[loc1 + FET_CYL_CENTER_A0] =
                              TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                        -INFTY, INFTY);

                            if (j == np)
                              Error(loc0, "Missing \"%c\" coordinate of cylinder center after \"%s\"", word[1], str);
                            WDB[loc1 + FET_CYL_CENTER_A1] =
                              TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                        -INFTY, INFTY);
                          }
                        break;
                      }

                    default:
                      {
                        /* Unsupported type */

                        Error(loc0, "Unsupported FET type %ld", type);
                      }
                    }
                }
              else if (!strcmp(str, "dir"))
                {
                  /* Detector direction vector */

                  /* Read vector */

                  if (j == np)
                    Error(loc0, "Missing x-component after \"%s\"", str);
                  else
                    WDB[loc0 + DET_DIRVEC_U] =
                      TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                -1.0, 1.0);

                  if (j == np)
                    Error(loc0, "Missing x-component after \"%s\"", str);
                  else
                    WDB[loc0 + DET_DIRVEC_V] =
                      TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                -1.0, 1.0);

                  if (j == np)
                    Error(loc0, "Missing x-component after \"%s\"", str);
                  else
                    WDB[loc0 + DET_DIRVEC_W] =
                      TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                -1.0, 1.0);
                }
              else if (!strcmp(str, "dx"))
                {
                  /* Cartesian mesh, x-binning */

                  if (j == np)
                    Error(loc0, "Missing minimum dimension after\"%s\"", str);

                  /* Get pointer or allocate memory for mesh structure */

                  if ((loc1 = (long)RDB[loc0 + DET_PTR_MESH]) < VALID_PTR)
                    {
                      /* Allocate memory */

                      loc1 = NewItem(loc0 + DET_PTR_MESH, MESH_BLOCK_SIZE);

                      /* Alloc memory for index (only seek once per collision) */

                      AllocValuePair(loc1 + MESH_PREV_COL_IDX);

                      /* Reset data */

                      WDB[loc1 + MESH_N0] = 1.0;
                      WDB[loc1 + MESH_MIN0] = -INFTY;
                      WDB[loc1 + MESH_MAX0] = INFTY;

                      WDB[loc1 + MESH_N1] = 1.0;
                      WDB[loc1 + MESH_MIN1] = -INFTY;
                      WDB[loc1 + MESH_MAX1] = INFTY;

                      WDB[loc1 + MESH_N2] = 1.0;
                      WDB[loc1 + MESH_MIN2] = -INFTY;
                      WDB[loc1 + MESH_MAX2] = INFTY;

                      /* Put type */

                      WDB[loc1 + MESH_TYPE] = (double)MESH_TYPE_CARTESIAN;
                    }
                  else
                    {
                      /* Check type */

                      if ((long)RDB[loc1 + MESH_TYPE] != MESH_TYPE_CARTESIAN)
                        Error(loc0, "Conflicting mesh type");
                    }

                  WDB[loc1 + MESH_MIN0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  if (j == np)
                    Error(loc0, "Missing maximum dimension after\"%s\"", str);

                  WDB[loc1 + MESH_MAX0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              RDB[loc1 + MESH_MIN0], INFTY);

                  if (j == np)
                    Error(loc0, "Missing number of bins after \"%s\"", str);

                  WDB[loc1 + MESH_N0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 1000000);
                }
              else if (!strcmp(str, "dy"))
                {
                  /* Cartesian mesh, y-binning */

                  if (j == np)
                    Error(loc0, "Missing minimum dimension after\"%s\"", str);

                  /* Get pointer or allocate memory for mesh structure */

                  if ((loc1 = (long)RDB[loc0 + DET_PTR_MESH]) < VALID_PTR)
                    {
                      /* Allocate memory */

                      loc1 = NewItem(loc0 + DET_PTR_MESH, MESH_BLOCK_SIZE);

                      /* Alloc memory for index (only seek once per collision) */

                      AllocValuePair(loc1 + MESH_PREV_COL_IDX);

                      /* Reset data */

                      WDB[loc1 + MESH_N0] = 1.0;
                      WDB[loc1 + MESH_MIN0] = -INFTY;
                      WDB[loc1 + MESH_MAX0] = INFTY;

                      WDB[loc1 + MESH_N1] = 1.0;
                      WDB[loc1 + MESH_MIN1] = -INFTY;
                      WDB[loc1 + MESH_MAX1] = INFTY;

                      WDB[loc1 + MESH_N2] = 1.0;
                      WDB[loc1 + MESH_MIN2] = -INFTY;
                      WDB[loc1 + MESH_MAX2] = INFTY;

                      /* Put type */

                      WDB[loc1 + MESH_TYPE] = (double)MESH_TYPE_CARTESIAN;
                    }
                  else
                    {
                      /* Check type */

                      if ((long)RDB[loc1 + MESH_TYPE] != MESH_TYPE_CARTESIAN)
                        Error(loc0, "Conflicting mesh type");
                    }

                  WDB[loc1 + MESH_MIN1] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  if (j == np)
                    Error(loc0, "Missing maximum dimension after\"%s\"", str);

                  WDB[loc1 + MESH_MAX1] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              RDB[loc1 + MESH_MIN1], INFTY);

                  if (j == np)
                    Error(loc0, "Missing number of bins after \"%s\"", str);

                  WDB[loc1 + MESH_N1] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 1000000);
                }
              else if (!strcmp(str, "dz"))
                {
                  /* Cartesian mesh, z-binning */

                  if (j == np)
                    Error(loc0, "Missing minimum dimension after\"%s\"", str);

                  /* Get pointer or allocate memory for mesh structure */

                  if ((loc1 = (long)RDB[loc0 + DET_PTR_MESH]) < VALID_PTR)
                    {
                      /* Allocate memory */

                      loc1 = NewItem(loc0 + DET_PTR_MESH, MESH_BLOCK_SIZE);

                      /* Alloc memory for index (only seek once per collision) */

                      AllocValuePair(loc1 + MESH_PREV_COL_IDX);

                      /* Reset data */

                      WDB[loc1 + MESH_N0] = 1.0;
                      WDB[loc1 + MESH_MIN0] = -INFTY;
                      WDB[loc1 + MESH_MAX0] = INFTY;

                      WDB[loc1 + MESH_N1] = 1.0;
                      WDB[loc1 + MESH_MIN1] = -INFTY;
                      WDB[loc1 + MESH_MAX1] = INFTY;

                      WDB[loc1 + MESH_N2] = 1.0;
                      WDB[loc1 + MESH_MIN2] = -INFTY;
                      WDB[loc1 + MESH_MAX2] = INFTY;

                      /* Put type */

                      WDB[loc1 + MESH_TYPE] = (double)MESH_TYPE_CARTESIAN;
                    }
                  else
                    {
                      /* Check type */

                      if ((long)RDB[loc1 + MESH_TYPE] != MESH_TYPE_CARTESIAN)
                        Error(loc0, "Conflicting mesh type");
                    }

                  WDB[loc1 + MESH_MIN2] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  if (j == np)
                    Error(loc0, "Missing maximum dimension after\"%s\"", str);

                  WDB[loc1 + MESH_MAX2] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              RDB[loc1 + MESH_MIN2], INFTY);

                  if (j == np)
                    Error(loc0, "Missing number of bins after \"%s\"", str);

                  WDB[loc1 + MESH_N2] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 1000000);
                }
              else if (!strcmp(str, "dn"))
                {
                  /* Curvilinear mesh */

                  if (j == np)
                    Error(loc0, "Missing type after\"%s\"", str);

                  /* Check mesh pointer */

                  if ((long)RDB[loc0 + DET_PTR_MESH] > VALID_PTR)
                    Error(loc0, "Multiple mesh types defined");

                  /* Read type */

                  n = (long)TestParam(pname, fname, line, params[j++],
                                      PTYPE_INT, -4, 4);

                  /* Check type */

                  if ((n == 1) || (n == -1))
                    {
                      /* Allocate memory for mesh structure */

                      loc1 = NewItem(loc0 + DET_PTR_MESH, MESH_BLOCK_SIZE);

                      /* Alloc memory for index (only seek once per collision) */

                      AllocValuePair(loc1 + MESH_PREV_COL_IDX);

                      /* Cylindrical mesh */

                      WDB[loc1 + MESH_TYPE] = (double)MESH_TYPE_CYLINDRICAL;

                      /* Store flag to use lowest level coordinates in */
                      /* scoring rather than universe 0 coordinates */

                      if (n == -1)
                        WDB[loc1 + MESH_LOCAL_COORDS] = (double)YES;

                      /* Check number of parameters */

                      if (j > np - 9)
                        Error(loc0, "Missing mesh parameters");

                      /* Read parameters */

                      WDB[loc1 + MESH_MIN0] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  0.0, INFTY);
                      WDB[loc1 + MESH_MAX0] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  RDB[loc1 + MESH_MIN0], INFTY);
                      WDB[loc1 + MESH_N0] =
                        TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                  1, 1000000);

                      WDB[loc1 + MESH_MIN1] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  0.0, 360.0);
                      WDB[loc1 + MESH_MAX1] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  RDB[loc1 + MESH_MIN1], 360.0);
                      WDB[loc1 + MESH_N1] =
                        TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                  1, 1000000);

                      WDB[loc1 + MESH_MIN2] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  -INFTY, INFTY);
                      WDB[loc1 + MESH_MAX2] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  RDB[loc1 + MESH_MIN2], INFTY);
                      WDB[loc1 + MESH_N2] =
                        TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                  1, 1000000);

                      /* Convert degrees to rad */

                      WDB[loc1 + MESH_MIN1] = RDB[loc1 + MESH_MIN1]*PI/180.0;
                      WDB[loc1 + MESH_MAX1] = RDB[loc1 + MESH_MAX1]*PI/180.0;
                    }
                  else if ((n == 2) || (n == -2))
                    {
                      /* Allocate memory for mesh structure */

                      loc1 = NewItem(loc0 + DET_PTR_MESH, MESH_BLOCK_SIZE);

                      /* Alloc memory for index (only seek once per collision) */

                      AllocValuePair(loc1 + MESH_PREV_COL_IDX);

                      /* Spherical mesh */

                      WDB[loc1 + MESH_TYPE] = (double)MESH_TYPE_SPHERICAL;

                      /* Store flag to use lowest level coordinates in */
                      /* scoring rather than universe 0 coordinates */

                      if (n == -2)
                        WDB[loc1 + MESH_LOCAL_COORDS] = (double)YES;

                      /* Check number of parameters */

                      if (j > np - 9)
                        Error(loc0, "Missing mesh parameters");

                      /* Read parameters */

                      WDB[loc1 + MESH_MIN0] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  0.0, INFTY);
                      WDB[loc1 + MESH_MAX0] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  RDB[loc1 + MESH_MIN0], INFTY);
                      WDB[loc1 + MESH_N0] =
                        TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                  1, 1000000);

                      WDB[loc1 + MESH_MIN1] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  0.0, 360.0);
                      WDB[loc1 + MESH_MAX1] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  RDB[loc1 + MESH_MIN1], 360.0);
                      WDB[loc1 + MESH_N1] =
                        TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                  1, 1000000);

                      WDB[loc1 + MESH_MIN2] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  0.0, 180.0);
                      WDB[loc1 + MESH_MAX2] =
                        TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                  RDB[loc1 + MESH_MIN2], 180.0);
                      WDB[loc1 + MESH_N2] =
                        TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                  1, 1000000);

                      /* Convert degrees to rad */

                      WDB[loc1 + MESH_MIN1] = RDB[loc1 + MESH_MIN1]*PI/180.0;
                      WDB[loc1 + MESH_MAX1] = RDB[loc1 + MESH_MAX1]*PI/180.0;
                      WDB[loc1 + MESH_MIN2] = RDB[loc1 + MESH_MIN2]*PI/180.0;
                      WDB[loc1 + MESH_MAX2] = RDB[loc1 + MESH_MAX2]*PI/180.0;
                    }
                  else if ((n == 3) || (n == -3))
                    {
                      /* Orthogonal mesh */

                      nx = -1;
                      ny = -1;
                      nz = -1;

                      /* Read size */

                      if (j > np - 1)
                        Error(loc0, "Missing mesh size (x)");
                      else
                        nx = (long)TestParam(pname, fname, line, params[j++],
                                             PTYPE_INT, 1, 10000);

                      if (j > np - 1)
                        Error(loc0, "Missing mesh size (y)");
                      else
                        ny = (long)TestParam(pname, fname, line, params[j++],
                                             PTYPE_INT, 1, 10000);

                      if (j > np - 1)
                        Error(loc0, "Missing mesh size (z)");
                      else
                        nz = (long)TestParam(pname, fname, line, params[j++],
                                             PTYPE_INT, 1, 10000);

                      /* Allocate memory for values */

                      mlims = (double *)Mem(MEM_ALLOC, nx + ny + nz + 3,
                                            sizeof(double));

                      /* Read limits */

                      for (m = 0; m < nx + ny + nz + 3; m++)
                        {
                          /* Read value */

                          if (j > np - 1)
                            Error(loc0, "Missing mesh boundary");
                          else
                            mlims[m] = TestParam(pname, fname, line,
                                                 params[j++], PTYPE_REAL,
                                                 -INFTY, INFTY);
                        }

                      /* Create structure */

                      loc1 = CreateMesh(MESH_TYPE_ORTHOGONAL,
                                        MESH_CONTENT_NONE, -1, nx, ny, nz,
                                        mlims, -1);

                      /* Put pointer */

                      WDB[loc0 + DET_PTR_MESH] = (double)loc1;

                      /* Store flag to use lowest level coordinates in */
                      /* scoring rather than universe 0 coordinates */

                      if (n == -3)
                        WDB[loc1 + MESH_LOCAL_COORDS] = (double)YES;

                      /* Free temporary array */

                      Mem(MEM_FREE, mlims);
                    }
                  else if ((n == 4) || (n == -4))
                    {
                      /* Unevenly-spaced cylindrical mesh */

                      nx = -1;
                      ny = -1;
                      nz = -1;

                      /* Read size */

                      if (j > np - 1)
                        Error(loc0, "Missing mesh size (r)");
                      else
                        nx = (long)TestParam(pname, fname, line, params[j++],
                                             PTYPE_INT, 1, 10000);

                      if (j > np - 1)
                        Error(loc0, "Missing mesh size (phi)");
                      else
                        ny = (long)TestParam(pname, fname, line, params[j++],
                                             PTYPE_INT, 1, 10000);

                      if (j > np - 1)
                        Error(loc0, "Missing mesh size (z)");
                      else
                        nz = (long)TestParam(pname, fname, line, params[j++],
                                             PTYPE_INT, 1, 10000);

                      /* Allocate memory for values */

                      mlims = (double *)Mem(MEM_ALLOC, nx + ny + nz + 3,
                                            sizeof(double));

                      /* Read limits */

                      for (m = 0; m < nx + ny + nz + 3; m++)
                        {
                          /* Read value */

                          if (j > np - 1)
                            Error(loc0, "Missing mesh boundary");
                          else
                            mlims[m] = TestParam(pname, fname, line,
                                                 params[j++],
                                                 PTYPE_REAL, -INFTY, INFTY);

                          /* Check */

                          if ((m < nx + 1) && (mlims[m] < 0.0))
                            Error(loc0, "Invalid radial mesh boundary %E",
                                  mlims[m]);
                          else if ((m > nx) && (m < nx + ny + 2))
                            {
                              /* Check angle */

                              if ((mlims[m] < 0.0) || (mlims[m] > 360))
                                Error(loc0,
                                      "Invalid angular mesh boundary %E",
                                      mlims[m]);
                              else
                                mlims[m] = mlims[m]*PI/180.0;
                            }
                        }

                      /* Create structure */

                      loc1 = CreateMesh(MESH_TYPE_ICYL, MESH_CONTENT_NONE,
                                        -1, nx, ny, nz, mlims, -1);

                      /* Put pointer */

                      WDB[loc0 + DET_PTR_MESH] = (double)loc1;

                      /* Store flag to use lowest level coordinates in */
                      /* scoring rather than universe 0 coordinates */

                      if (n == -4)
                        WDB[loc1 + MESH_LOCAL_COORDS] = (double)YES;

                      /* Free temporary array */

                      Mem(MEM_FREE, mlims);
                    }
                  else
                    Die(FUNCTION_NAME, "Invalid mesh type");
                }
              else if (!strcmp(str, "dh"))
                {
                  /* Hexagonal mesh mesh */

                  if (j == np)
                    Error(loc0, "Missing type after\"%s\"", str);

                  /* Check mesh pointer */

                  if ((long)RDB[loc0 + DET_PTR_MESH] > VALID_PTR)
                    Error(loc0, "Multiple mesh types defined");

                  /* Allocate memory for mesh structure */

                  loc1 = NewItem(loc0 + DET_PTR_MESH, MESH_BLOCK_SIZE);

                  /* Alloc memory for index (only seek once per collision) */

                  AllocValuePair(loc1 + MESH_PREV_COL_IDX);

                  /* Read type */

                  n = (long)TestParam(pname, fname, line, params[j++],
                                      PTYPE_INT, 2, 3);

                  /* Check type */

                  if (n == 2)
                    WDB[loc1 + MESH_TYPE] = (double)MESH_TYPE_HEXX;
                  else
                    WDB[loc1 + MESH_TYPE] = (double)MESH_TYPE_HEXY;

                  /* Check number of parameters */

                  if (j > np - 8)
                    Error(loc0, "Missing mesh parameters");

                  /* Read central coordinates */

                  WDB[loc1 + MESH_MIN0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);
                  WDB[loc1 + MESH_MIN1] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  /* Read pitch */

                  WDB[loc1 + MESH_MAX0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);

                  /* Read size */

                  WDB[loc1 + MESH_N0] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 1000000);
                  WDB[loc1 + MESH_N1] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 1000000);

                  /* Read axial binning */

                  WDB[loc1 + MESH_MIN2] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              -INFTY, INFTY);
                  WDB[loc1 + MESH_MAX2] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              RDB[loc1 + MESH_MIN2], INFTY);
                  WDB[loc1 + MESH_N2] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              1, 1000000);
                }
              else if (!strcmp(str, "dt"))
                {
                  /* Type. Check number of parameters. */

                  if (j == np)
                    Error(loc0, "Missing detector type after \"%s\"", str);

                  /* Set type */

                  WDB[loc0 + DET_TYPE] =
                    TestParam(pname, fname, line, params[j++], PTYPE_INT,
                              -6, 4);

                  /* Check if divider or multiplier */

                  if (((long)RDB[loc0 + DET_TYPE] == DETECTOR_TYPE_DIVI) ||
                      ((long)RDB[loc0 + DET_TYPE] == DETECTOR_TYPE_MULTI))
                    {
                      /* Check number of parameters. */

                      if (j == np)
                        Error(loc0,
                              "Missing divider or multiplier after \"%s\"",
                              str);

                      /* Set name */

                      WDB[loc0 + DET_PTR_MUL] = (double)PutText(params[j++]);
                    }
                }
              else if (!strcmp(str, "dr"))
                {
                  /* Detector reaction, allocate memory for structure */

                  loc1 = NewItem(loc0 + DET_PTR_RBINS, DET_RBIN_BLOCK_SIZE);

                  /* Avoid compiler warning */

                  ne = -1;

                  /* Read MT */

                  if (j == np)
                    Error(loc0, "Missing reaction mt after \"%s\"", str);
                  else
                    {
                      /* Check last character and put isomer flag */

                      c = params[j][strlen(params[j]) - 1];

                      if ((c == 'g') || (c == 'G'))
                        {
                          WDB[loc1 + DET_RBIN_RFS] = 0.0;
                          params[j][strlen(params[j]) - 1] = '\0';
                        }
                      else if ((c == 'i') || (c == 'I') ||
                               (c == 'm') || (c == 'M'))
                        {
                          WDB[loc1 + DET_RBIN_RFS] = 1.0;
                          params[j][strlen(params[j]) - 1] = '\0';
                        }
                      else
                        WDB[loc1 + DET_RBIN_RFS] = -1.0;

                      /* Put value */

                      WDB[loc1 + DET_RBIN_MT] =
                        TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                  -248, 120000);
                    }

                  /* Check for user-defined mode */

                  if ((long)RDB[loc1 + DET_RBIN_MT] == MT_USER_DEFINED)
                    {
                      /* Read response function */

                      WDB[loc1 + DET_RBIN_PTR_FUN] =
                        (double)PutText(params[j++]);
                    }
                  else
                    {
                      /* Read material name */

                      if (j == np)
                        Error(loc0, "Missing material name after \"%s\"", str);
                      else
                        WDB[loc1 + DET_RBIN_PTR_MAT] =
                          (double)PutText(params[j++]);
                    }
                }
              else if (!strcmp(str, "du"))
                {
                  /* Detector universe, allocate memory for structure */

                  loc1 = NewItem(loc0 + DET_PTR_UBINS, DET_UBIN_BLOCK_SIZE);

                  /* Read universe */

                  if (j == np)
                    Error(loc0, "Missing universe after \"%s\"", str);
                  else
                    WDB[loc1 + DET_UBIN_PTR_UNI] =
                      (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "dl"))
                {
                  /* Check previous definitions */

                  if ((long)RDB[loc0 + DET_PTR_LBINS] > VALID_PTR)
                    Error(loc0, "Multiple detector lattices not allowed");

                  /* Detector lattice, allocate memory for structure */

                  loc1 = NewItem(loc0 + DET_PTR_LBINS, DET_LBIN_BLOCK_SIZE);

                  /* Read universe */

                  if (j == np)
                    Error(loc0, "Missing lattice after \"%s\"", str);
                  else
                    WDB[loc1 + DET_LBIN_PTR_LAT] =
                      (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "dc"))
                {
                  /* Detector cell, allocate memory for structure */

                  loc1 = NewItem(loc0 + DET_PTR_CBINS, DET_CBIN_BLOCK_SIZE);

                  /* Read cell */

                  if (j == np)
                    Error(loc0, "Missing cell after \"%s\"", str);
                  else
                    WDB[loc1 + DET_CBIN_PTR_CELL] =
                      (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "dm"))
                {
                  /* Detector material, allocate memory for structure */

                  loc1 = NewItem(loc0 + DET_PTR_MBINS, DET_MBIN_BLOCK_SIZE);

                  /* Read material */

                  if (j == np)
                    Error(loc0, "Missing material after \"%s\"", str);
                  else
                    WDB[loc1 + DET_MBIN_PTR_MAT] =
                      (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "da"))
                {
                  /* Activation detector, allocate memory for structure */

                  loc1 = NewItem(loc0 + DET_PTR_ABINS, DET_ABIN_BLOCK_SIZE);

                  /* Read material */

                  if (j == np)
                    Error(loc0, "Missing material after \"%s\"", str);
                  else
                    WDB[loc1 + DET_ABIN_PTR_MAT] =
                      (double)PutText(params[j++]);

                  /* Read flux */

                  if (j == np)
                    Error(loc0, "Missing flux after \"%s\"", str);
                  else
                    WDB[loc1 + DET_ABIN_FLUX] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_REAL, -1.0, INFTY);
                }
              else if (!strcmp(str, "df"))
                {
                  /* Detector source file, read file name */

                  if (j == np)
                    Error(loc0, "Missing file name after \"%s\"", str);
                  else
                    WDB[loc0 + DET_WRITE_PTR_FILE] =
                      (double)PutText(params[j++]);

                  /* Avoid compiler warning */

                  val = 0.0;

                  /* Read fraction */

                  if (j == np)
                    Error(loc0, "Missing fraction after \"%s\"", str);
                  else
                    val = TestParam(pname, fname, line, params[j++],
                                    PTYPE_REAL, -1.0, 1.0);

                  /* Check binary flag */

                  if (val > 0.0)
                    {
                      WDB[loc0 + DET_WRITE_PROB] = val;
                      WDB[loc0 + DET_WRITE_BINARY] = (double)NO;
                    }
                  else
                    {
                      WDB[loc0 + DET_WRITE_PROB] = -val;
                      WDB[loc0 + DET_WRITE_BINARY] = (double)YES;
                    }
                }
              else if (!strcmp(str, "ds"))
                {
                  /* Surface current detector, allocate memory for structure */

                  loc1 = NewItem(loc0 + DET_PTR_SBINS, DET_SBIN_BLOCK_SIZE);

                  /* Put type */

                  WDB[loc1 + DET_SBIN_TYPE] = (double)SUPERDET_TYPE_CURRENT;

                  /* Read surface name */

                  if (j == np)
                    Error(loc0, "Missing surface after \"%s\"", str);
                  else
                    WDB[loc1 + DET_SBIN_PTR_SURF] =
                      (double)PutText(params[j++]);

                  /* Read normal direction */

                  if (j == np)
                    Error(loc0, "Missing normal direction after \"%s\"", str);
                  else
                    WDB[loc1 + DET_SBIN_SURF_NORM] =
                      TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                -2, 1);
                }
              else if (!strcmp(str, "dtl"))
                {
                  /* Super-imposet track-length volume flux detector,  */
                  /* allocate memory for structure (use same structure */
                  /* as surface current detector) */

                  loc1 = NewItem(loc0 + DET_PTR_SBINS, DET_SBIN_BLOCK_SIZE);

                  /* Put type */

                  WDB[loc1 + DET_SBIN_TYPE] = (double)SUPERDET_TYPE_TLEFLUX;

                  /* Read surface name */

                  if (j == np)
                    Error(loc0, "Missing surface after \"%s\"", str);
                  else
                    WDB[loc1 + DET_SBIN_PTR_SURF] =
                      (double)PutText(params[j++]);
                }
              else if (!strcmp(str, "dhis"))
                {
                  /* History collection for detector */

                  if (j == np)
                    Error(loc0, "Missing history indicator after \"%s\"",str);
                  else
                    WDB[loc0 + DET_WRITE_HIS] =
                      TestParam(pname, fname, line, params[j++], PTYPE_LOGICAL);

                  /* Set global option */

                  if ((long)RDB[loc0 + DET_WRITE_HIS] == YES)
                    WDB[DATA_RUN_STAT_TESTS] = (double)YES;
                }
              else if (!strcmp(str, "dg"))
                {
                  /* Delayed group flagging for detector */

                  if (j == np)
                    Error(loc0, "Missing delayed group number after \"%s\"",
                          str);
                  else
                    WDB[loc0 + DET_DELAYED_GROUP] =
                      TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                0, 8);
                }
              else if (!strcmp(str, "dumsh"))
                {
                  /* Mapping for UMSH cells */

                  loc1 = NewItem(loc0 + DET_PTR_CBINS, DET_CBIN_BLOCK_SIZE);

                  /* Read umsh universe */

                  if (j == np)
                    Error(loc0, "Missing universe name \"%s\"",str);
                  else
                    WDB[loc1 + DET_CBIN_UMSH_PTR_UMSH] =
                      (double)PutText(params[j++]);

                  /* Avoid compiler warning */

                  ni = -1;

                  /* Read number of value pairs */

                  if (j == np)
                    Error(loc0,
                          "Missing number of cell/bin pairs after \"%s\"",str);
                  else
                    ni = (long)TestParam(pname, fname, line, params[j++],
                                         PTYPE_INT, 1, 10000000000);

                  /* Allocate memory for data */

                  loc2 = ReallocMem(DATA_ARRAY, ni + 1);
                  WDB[loc1 + DET_CBIN_UMSH_PTR_CELLS] = (double)loc2;

                  loc3 = ReallocMem(DATA_ARRAY, ni + 1);
                  WDB[loc1 + DET_CBIN_UMSH_PTR_BINS] = (double)loc3;

                  /* Read data */

                  for (n = 0; n < ni; n++)
                    {
                      /* Check */

                      if (j > np - 2)
                        Error(loc0, "Incorrect number of values after \"%s\"",
                              str);

                      /* Read values */

                      WDB[loc2++] = TestParam(pname, fname, line,
                                              params[j++], PTYPE_INT,
                                              1, 10000000000);

                      WDB[loc3++] = TestParam(pname, fname, line,
                                              params[j++], PTYPE_INT,
                                              1, 10000000000);
                    }

                  /* Put terminators */

                  WDB[loc2] = -1;
                  WDB[loc3] = -1;
                }
              else
                Error(loc0, "Invalid detector parameter \"%s\"", str);
            }
         }

      /***********************************************************************/

      /***** User-defined energy grid ****************************************/

      else if (!strcasecmp(word, "ene"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 3, 1E+6, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_ENE0, ENE_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Grid name */

          WDB[loc0 + ENE_PTR_NAME] = (double)PutText(params[j++]);

          /* Get type */

          type = (long)TestParam(pname, fname, line, params[j++],
                                 PTYPE_INT, 1, 4);

          /* Put type */

          WDB[loc0 + ENE_TYPE] = (double)type;

          /* Read data */

          switch(type)
            {
            case EG_TYPE_ARB:
              {
                /* Arbitrary grid. Get number of values */

                ne = np - j;

                /* Set number of bins */

                if (ne < 2)
                  Error(loc0, "Not enough boundaries to define a grid");
                else
                  WDB[loc0 + ENE_NB] = (double)(ne - 1);

                /* Allocate memory for data */

                loc1 = ReallocMem(DATA_ARRAY, ne);

                /* Set pointer */

                WDB[loc0 + ENE_PTR_GRID] = (double)loc1;

                /* Read values */

                while(j < np)
                  WDB[loc1++] = TestParam(pname, fname, line, params[j++],
                                           PTYPE_REAL, 0.0, 1000.0);

                /* Break case */

                break;
              }
            case EG_TYPE_UNI_E:
            case EG_TYPE_UNI_L:
              {
                /* Uniform grid. Check number of parameters */

                if (np != 5)
                  Error(loc0, "Uniform grid %s must be defined by 5 values",
                        GetText(loc0 + ENE_PTR_NAME));

                /* Read number of energy bins */

                WDB[loc0 + ENE_NB] =
                  TestParam(pname, fname, line, params[j++], PTYPE_INT,
                            1, 1000000);

                /* Read minimum and maximum energy */

                WDB[loc0 + ENE_EMIN] =
                  TestParam(pname, fname, line, params[j++],
                            PTYPE_REAL, 0.0, 1000.0);

                WDB[loc0 + ENE_EMAX] =
                  TestParam(pname, fname, line, params[j++],
                            PTYPE_REAL, RDB[loc0 + ENE_EMIN], 1000.0);

                /* Check zero in logarithmic grid */

                if ((type == EG_TYPE_UNI_L) && (RDB[loc0 + ENE_EMIN] == 0.0))
                  Error(loc0, "Lower limit must be above zero");

                /* Break case */

                break;
              }
            case EG_TYPE_PREDEF:
              {
                if (j == np)
                  Error(loc0, "Missing grid name after \"%s\"", str);

                /* Predefined grid, read name */

                WDB[loc0 + ENE_PTR_PREDEF] = (double)PutText(params[j++]);

                /* Break case */

                break;
              }
            default:
              Error(loc0, "Invalid energy grid type %d", type);
            }
         }

      /***********************************************************************/

      /***** Weight window structure *****************************************/

      else if (!strcasecmp(word, "wwin"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, 1E+6, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_WWD0, WWD_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Avoid compiler warning */

          lims[0] = 0.0;
          lims[1] = 0.0;
          lims[2] = 0.0;
          lims[3] = 0.0;
          lims[4] = 0.0;
          lims[5] = 0.0;

          nx = -1;
          ny = -1;
          nz = -1;

          /* Reset parameters */

          WDB[loc0 + WWD_NORM_FACT] = -1.0;
          WDB[loc0 + WWD_POW] = -INFTY;
          WDB[loc0 + WWD_MULT] = -INFTY;
          WDB[loc0 + WWD_SRC_BIAS] = (double)NO;
          WDB[loc0 + WWD_BOUNDS_TYPE] = (double)WW_MESH_BOUNDS_AVG;
          WDB[loc0 + WWD_LIM_MIN] = -1.0;
          WDB[loc0 + WWD_LIM_MAX] = -1.0;
          WDB[loc0 + WWD_SMOOTH_F0] = -1.0;
          WDB[loc0 + WWD_SMOOTH_F1] = INFTY;

          /* Read data */

          j = 0;

          /* Name */

          WDB[loc0 + WWD_PTR_NAME] = (double)PutText(params[j++]);

          /* Loop over remaining parameters */

          while (j < np)
            {
              /* Check parameter type */

              strcpy(str, params[j++]);

              if (!strcmp(str, "wn"))
                {
                  /* Read normalization */

                  if (j > np - 1)
                    Error(loc0, "Missing renormalization factor");
                  else
                    WDB[loc0 + WWD_NORM_FACT] =
                      TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                -INFTY, INFTY);

                  /* Check if point is given */

                  if (RDB[loc0 + WWD_NORM_FACT] > 0.0)
                    {
                      if (j > np - 1)
                        Error(loc0, "Missing coordinates");
                      else
                        WDB[loc0 + WWD_NORM_X] =
                          TestParam(pname, fname, line, params[j++],
                                    PTYPE_REAL, -INFTY, INFTY);

                      if (j > np - 1)
                        Error(loc0, "Missing coordinates");
                      else
                        WDB[loc0 + WWD_NORM_Y] =
                          TestParam(pname, fname, line, params[j++],
                                    PTYPE_REAL, -INFTY, INFTY);

                      if (j > np - 1)
                        Error(loc0, "Missing coordinates");
                      else
                        WDB[loc0 + WWD_NORM_Z] =
                          TestParam(pname, fname, line, params[j++],
                                    PTYPE_REAL, -INFTY, INFTY);

                      if (j > np - 1)
                        Error(loc0, "Missing energy");
                      else
                        WDB[loc0 + WWD_NORM_E] =
                          TestParam(pname, fname, line, params[j++],
                                    PTYPE_REAL, ZERO, INFTY);
                    }
                  else if (RDB[loc0 + WWD_NORM_FACT] == 0.0)
                    Error(loc0, "Renormalization factor must be nonzero");
                }
              else if (!strcmp(str, "wt"))
                {
                  /* Read source biasing flag */

                  if (j < np)
                    WDB[loc0 + WWD_SRC_BIAS] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_LOGICAL);

                  /* Read bounds type */

                  if (j < np)
                    WDB[loc0 + WWD_BOUNDS_TYPE] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_INT, 1, 2);

                  /* Read limits */

                  if ((j < np) && (params[j][0] != 'w'))
                    WDB[loc0 + WWD_LIM_MIN] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_REAL, -1.0, INFTY);

                  if ((j < np) && (params[j][0] != 'w'))
                    WDB[loc0 + WWD_LIM_MAX] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_REAL, -1.0, INFTY);

                  /* Read smoothing factors */
                  /*
                  if ((j < np) && (params[j][0] != 'w'))
                    WDB[loc0 + WWD_SMOOTH_F0] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_REAL, 1.0, INFTY);

                  if ((j < np) && (params[j][0] != 'w'))
                    WDB[loc0 + WWD_SMOOTH_F1] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_REAL, 1.0, INFTY);
                  */
                  /* Check */

                  if ((long)RDB[loc0 + WWD_SRC_BIAS] == YES)
                    Error(loc0, "Source biasing not working in this version");
                }
              else if (!strcmp(str, "wf"))
                {
                  /* Read file name */

                  WDB[loc0 + WWD_PTR_FNAME] = (double)PutText(params[j++]);

                  /* Read type */

                  if (j > np - 1)
                    Error(loc0, "Missing file type");
                  else
                    WDB[loc0 + WWD_TYPE] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_INT, 1, 2);
                }
              else if (!strcmp(str, "wx"))
                {
                  /* Read multiplier */

                  if (j > np - 1)
                    Error(loc0, "Missing multiplier");
                  else
                    WDB[loc0 + WWD_MULT] =
                      TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                ZERO, INFTY);

                  /* Read exponent */

                  if (j > np - 1)
                    Error(loc0, "Missing exponential");
                  else
                    WDB[loc0 + WWD_POW] =
                      TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                                -INFTY, INFTY);
                }
              else if (!strcmp(str, "wi"))
                {
                  /* Avoid compiler warning */

                  type = -1;

                  /* Read iteration type */

                  if (j > np - 1)
                    Error(loc0, "Missing iteration type");
                  else
                    type = TestParam(pname, fname, line, params[j++],
                                     PTYPE_INT, 1, 4);

                  /* Read number of iterations */

                  if (j > np - 1)
                    Error(loc0, "Missing number of iterations");
                  else
                    WDB[loc0 + WWD_N_ITER] =
                      TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                1, 1000000);

                  /* Check type */

                  if ((type == WWIN_ITER_GEO) || (type == WWIN_ITER_TRK) ||
                      (type == WWIN_ITER_CON))
                    {
                      /* Create structure */

                      loc1 = NewItem(loc0 + WWD_PTR_ITER, WWD_ITER_BLOCK_SIZE);

                      /* Put type */

                      WDB[loc1 + WWD_ITER_TYPE] = (double)type;

                      /* Put density factor */

                      WDB[loc1 + WWD_ITER_DF] = 1.0;

                      /* Put split sizes (use octree mesh) */

                      WDB[loc1 + WWD_ITER_SPLIT_NX] = 2.0;
                      WDB[loc1 + WWD_ITER_SPLIT_NY] = 2.0;
                      WDB[loc1 + WWD_ITER_SPLIT_NZ] = 2.0;

                      /* Read pointer to response matrix solver */

                      if (j > np - 1)
                        Error(loc0, "Missing RMX name");
                      else
                        WDB[loc1 + WWD_ITER_PTR_RMX] =
                          PutText(params[j++]);

                      /* Split size */

                      if (j > np - 1)
                        Error(loc0, "Missing split size");
                      else
                        WDB[loc1 + WWD_ITER_SPLIT_NX] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    1, 100);

                      if (j > np - 1)
                        Error(loc0, "Missing split size");
                      else
                        WDB[loc1 + WWD_ITER_SPLIT_NY] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    1, 100);

                      if (j > np - 1)
                        Error(loc0, "Missing split size");
                      else
                        WDB[loc1 + WWD_ITER_SPLIT_NZ] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    1, 100);

                      if (j > np - 1)
                        Error(loc0, "Missing number of loops");
                      else
                        WDB[loc1 + WWD_ITER_SPLIT_LOOPS] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    0, 1000);

                      if (j > np - 1)
                        Error(loc0, "Missing number of tracks");
                      else
                        WDB[loc1 + WWD_ITER_SPLIT_TRACKS] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    1, 1000000000);

                      if (j > np - 1)
                        Error(loc0, "Missing importance split criteria");
                      else
                        WDB[loc1 + WWD_ITER_SPLIT_IMP_CRIT] =
                          TestParam(pname, fname, line, params[j++],
                                    PTYPE_REAL, 1.0, 1E18);

                      if (j > np - 1)
                        Error(loc0, "Missing maximum number of neighbours");
                      else
                        WDB[loc1 + WWD_ITER_SPLIT_MAX_NEIGH] =
                          TestParam(pname, fname, line, params[j++], PTYPE_INT,
                                    1, 1000000);

                      /* Loop over density criteria */

                      while (j < np)
                        {
                          /* Create structure */

                          loc2 = NewItem(loc1 + WWD_ITER_SPLIT_PTR_DENS,
                                         WWD_ITER_SPLIT_DENS_BLOCK_SIZE);

                          /* Read density */

                          if (j > np - 1)
                            Error(loc0, "Missing density");
                          else
                            WDB[loc2 + WWD_ITER_SPLIT_DENS_RHO] =
                              TestParam(pname, fname, line, params[j++],
                                        PTYPE_REAL, -INFTY, INFTY);

                          /* Read monimum cell size */

                          if (j > np - 1)
                            Error(loc0, "Missing minimum cell size");
                          else
                            WDB[loc2 + WWD_ITER_SPLIT_DENS_XMIN] =
                              TestParam(pname, fname, line, params[j++],
                                        PTYPE_REAL, ZERO, INFTY);

                          if (j > np - 1)
                            Error(loc0, "Missing minimum cell size");
                          else
                            WDB[loc2 + WWD_ITER_SPLIT_DENS_YMIN] =
                              TestParam(pname, fname, line, params[j++],
                                        PTYPE_REAL, ZERO, INFTY);

                          if (j > np - 1)
                            Error(loc0, "Missing minimum cell size");
                          else
                            WDB[loc2 + WWD_ITER_SPLIT_DENS_ZMIN] =
                              TestParam(pname, fname, line, params[j++],
                                        PTYPE_REAL, ZERO, INFTY);
                        }

                      /* Check pointer */

                      if ((loc2 = (long)RDB[loc1 + WWD_ITER_SPLIT_PTR_DENS])
                          < VALID_PTR)
                        Error(loc0, "No density split criteria");
                      else
                        {
                          /* Put size */

                          WDB[DATA_MAX_RMX_SPLIT_FLAGS] =
                            (double)ListSize(loc2);

                          /* Close list */

                          CloseList(loc2);

                          /* Check densities */

                          while (loc2 > VALID_PTR)
                            {
                              /* Compare densities */

                              if ((loc3 = NextItem(loc2)) > VALID_PTR)
                                if (RDB[loc2 + WWD_ITER_SPLIT_DENS_RHO]*
                                    RDB[loc3 + WWD_ITER_SPLIT_DENS_RHO] <= 0.0)
                                  Error(loc0, "Inconsistent density units");

                              /* Next */

                              loc2 = NextItem(loc2);
                            }

                          /* Sort by density */

                          loc2 = (long)RDB[loc1 + WWD_ITER_SPLIT_PTR_DENS];
                          CheckPointer(FUNCTION_NAME, "(loc2)",
                                       DATA_ARRAY, loc2);

                          if (RDB[loc2 + WWD_ITER_SPLIT_DENS_RHO] < 0.0)
                            SortList(loc2, WWD_ITER_SPLIT_DENS_RHO,
                                     SORT_MODE_ASCEND);
                          else
                            SortList(loc2, WWD_ITER_SPLIT_DENS_RHO,
                                     SORT_MODE_DESCEND);
                        }

                      /* Copy to remaining */

                      for (n = 0; n < (long)RDB[loc0 + WWD_N_ITER] - 1; n++)
                        DuplicateItem(loc1);
                    }
                  else if (type == WWIN_ITER_USR)
                    {
                      /* Loop over iterations */

                      for (n = 0; n < (long)RDB[loc0 + WWD_N_ITER]; n++)
                        {
                          /* Create structure */

                          loc1 = NewItem(loc0 + WWD_PTR_ITER,
                                         WWD_ITER_BLOCK_SIZE);

                          /* Put type */

                          WDB[loc1 + WWD_ITER_TYPE] = (double)type;

                          /* Read pointer to response matrix solver */

                          if (j > np - 1)
                            Error(loc0, "Missing RMX name");
                          else
                            WDB[loc1 + WWD_ITER_PTR_RMX] =
                              PutText(params[j++]);

                          /* Read density factor */

                          if (j > np - 1)
                            Error(loc0, "Missing density factor");
                          else
                            WDB[loc1 + WWD_ITER_DF] =
                              TestParam(pname, fname, line, params[j++],
                                        PTYPE_REAL, 0.0, 1.0);
                        }
                    }
                  else
                    Die(FUNCTION_NAME, "Invalid iteration type");

                  /* Set VR iterations on */

                  WDB[DATA_RUN_VR_ITER] = (double)YES;
                }
              else
                Error(loc0, "Invalid weight window parameter \"%s\"", str);
            }

          /* Set flag */

          WDB[DATA_USE_WEIGHT_WINDOWS] = (double)YES;
        }

      /***********************************************************************/

      /***** Response matrix solver ******************************************/

      else if (!strcasecmp(word, "wwgen"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 6, 1E+6, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_RMX0, RMX_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Avoid compiler warning */

          lims[0] = 0.0;
          lims[1] = 0.0;
          lims[2] = 0.0;
          lims[3] = 0.0;
          lims[4] = 0.0;
          lims[5] = 0.0;

          mlims = NULL;

          nx = -1;
          ny = -1;
          nz = -1;

          loc1 = -1;

          /* Read data */

          j = 0;

          /* Name */

          WDB[loc0 + RMX_PTR_NAME] = (double)PutText(params[j++]);

          /* Read convergence criterion */

          WDB[loc0 + RMX_CONV] = TestParam(pname, fname, line, params[j++],
                                           PTYPE_REAL, 0.0, 1.0);

          /* Read number of iterations */

          WDB[loc0 + RMX_N_ITER] =
            TestParam(pname, fname, line, params[j++],
                      PTYPE_INT, -1, 1000000);

          /* Check that convergence limit is provided */

          if (((long)RDB[loc0 + RMX_N_ITER] < 1) &&
              (RDB[loc0 + RMX_CONV] == 0.0))
            Error(loc0, "Convergence criterion must be provided");

          /* Read mode */

          WDB[loc0 + RMX_MODE] = TestParam(pname, fname, line, params[j++],
                                           PTYPE_INT, 1, 4);

          /* Read energy grid */

          if (!strcmp(params[j], "-1"))
            {
              WDB[loc0 + RMX_PTR_EGRID] = NULLPTR;
              j++;
            }
          else
            WDB[loc0 + RMX_PTR_EGRID] = (double)PutText(params[j++]);

          /* Read mesh type */

          type = (long)TestParam(pname, fname, line, params[j++],
                                         PTYPE_INT, -1, 8);

          /* Check */

          if (type == MESH_TYPE_CARTESIAN)
            {
              /* Check number of parameters */

              if (j > np - 9)
                Error(loc0, "Missing mesh parameters");

              /* Read limits and sizes */

              lims[0] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);
              lims[1] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, lims[0], INFTY);
              nx = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              lims[2] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);
              lims[3] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, lims[2], INFTY);
              ny = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              lims[4] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);
              lims[5] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, lims[4], INFTY);
              nz = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);
            }
          else if ((type == MESH_TYPE_HEXX) ||
                   (type == MESH_TYPE_HEXY))
            {
              /* Check number of parameters */

              if (j > np - 8)
                Error(loc0, "Missing mesh parameters");

              /* Read central coordinates */

              lims[0] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);

              lims[2] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);

              /* Test if centred */

              if ((lims[0] != 0.0) || (lims[2] != 0.0))
                Error(loc0, "Hex mesh must be centered at the origin");

              /* Read pitch */

              lims[1] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, 0.0, INFTY);

              /* Read size */

              nx = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              ny = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              /* Read axial binning */


              lims[4] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);

              lims[5] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);

              nz = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);
            }
          else if (type == MESH_TYPE_SPHERICAL)
            {
              /* Check number of parameters */

              if (j > np - 9)
                Error(loc0, "Missing mesh parameters");

              /* Read limits and sizes */

              lims[0] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, 0.0, INFTY);
              lims[1] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, lims[0], INFTY);
              nx = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              lims[2] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, 0.0, 360.0);
              lims[3] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, lims[2], 360.0);
              ny = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              lims[4] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, 0.0, 180.0);
              lims[5] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, lims[4], 180.0);
              nz = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              /* Check unsupported types */

              Error(loc0, "Spherical mesh type not supported yet");

              /* Convert degrees to rad */

              lims[2] = lims[2]*PI/180.0;
              lims[3] = lims[3]*PI/180.0;
              lims[4] = lims[4]*PI/180.0;
              lims[5] = lims[5]*PI/180.0;
            }
          else if (type == MESH_TYPE_CYLINDRICAL)
            {
              /* Check number of parameters */

              if (j > np - 9)
                Error(loc0, "Missing mesh parameters");

              /* Read limits and sizes */

              lims[0] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, 0.0, INFTY);
              lims[1] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, lims[0], INFTY);
              nx = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              lims[2] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, 0.0, 360.0);
              lims[3] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, lims[2], 360.0);
              ny = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              lims[4] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, -INFTY, INFTY);
              lims[5] = TestParam(pname, fname, line, params[j++],
                                  PTYPE_REAL, lims[4], INFTY);
              nz = (long)TestParam(pname, fname, line, params[j++],
                                   PTYPE_INT, 1, 10000);

              if (fabs(lims[2] - lims[3]) != 360.0)
                Error(loc0, "Angular binning must cover full sector");

              /* Convert degrees to rad */

              lims[2] = lims[2]*PI/180.0;
              lims[3] = lims[3]*PI/180.0;
            }
          else if (type == MESH_TYPE_ORTHOGONAL)
            {
              /* Read size */

              if (j > np - 1)
                Error(loc0, "Missing mesh size (x)");
              else
                nx = (long)TestParam(pname, fname, line, params[j++],
                                     PTYPE_INT, 1, 10000);

              if (j > np - 1)
                Error(loc0, "Missing mesh size (y)");
              else
                ny = (long)TestParam(pname, fname, line, params[j++],
                                     PTYPE_INT, 1, 10000);

              if (j > np - 1)
                Error(loc0, "Missing mesh size (z)");
              else
                nz = (long)TestParam(pname, fname, line, params[j++],
                                     PTYPE_INT, 1, 10000);

              /* Allocate memory for values */

              mlims = (double *)Mem(MEM_ALLOC, nx + ny + nz + 3,
                                    sizeof(double));

              /* Read limits */

              for (n = 0; n < nx + ny + nz + 3; n++)
                {
                  /* Read value */

                  if (j > np - 1)
                    Error(loc0, "Missing mesh boundary");
                  else
                    mlims[n] = TestParam(pname, fname, line, params[j++],
                                        PTYPE_REAL, -INFTY, INFTY);
                }
            }
          else if (type == MESH_TYPE_ICYL)
            {
              /* Read size */

              if (j > np - 1)
                Error(loc0, "Missing mesh size (r)");
              else
                nx = (long)TestParam(pname, fname, line, params[j++],
                                     PTYPE_INT, 1, 10000);

              if (j > np - 1)
                Error(loc0, "Missing mesh size (phi)");
              else
                ny = (long)TestParam(pname, fname, line, params[j++],
                                     PTYPE_INT, 1, 10000);

              if (j > np - 1)
                Error(loc0, "Missing mesh size (z)");
              else
                nz = (long)TestParam(pname, fname, line, params[j++],
                                     PTYPE_INT, 1, 10000);

              /* Allocate memory for values */

              mlims = (double *)Mem(MEM_ALLOC, nx + ny + nz + 3,
                                    sizeof(double));

              /* Read limits */

              for (n = 0; n < nx + ny + nz + 3; n++)
                {
                  /* Read value */

                  if (j > np - 1)
                    Error(loc0, "Missing mesh boundary");
                  else
                    mlims[n] = TestParam(pname, fname, line, params[j++],
                                        PTYPE_REAL, -INFTY, INFTY);

                  /* Check */

                  if ((n < nx + 1) && (mlims[n] < 0.0))
                    Error(loc0, "Invalid radial mesh boundary %E",
                          mlims[n]);
                  else if ((n > nx) && (n < nx + ny + 2))
                    {
                      /* Check angle */

                      if ((mlims[n] < 0.0) || (mlims[n] > 360))
                        Error(loc0, "Invalid angular mesh boundary %E",
                              mlims[n]);
                      else
                        mlims[n] = mlims[n]*PI/180.0;
                    }
                }

              /* Check that radial mesh starts from zero */

              if (mlims[0] != 0.0)
                Error(loc0, "First radial mesh boundary must be zero");

              /* Check angular boundaries */

              if (fabs(mlims[nx + 1] - mlims[nx + ny + 1]) != 2.0*PI)
                Error(loc0, "Angular binning must cover full sector");
            }
          else if (type == -1)
            {
              /* Check energy grid */

              if ((long)RDB[loc0 + RMX_PTR_EGRID] > VALID_PTR)
                Error(loc0, "Energy grid definition not allowed");
            }
          else
            Error(loc0, "Mesh type %ld not supported", type);

          /* Check type */

          if ((type == MESH_TYPE_ORTHOGONAL) || (type == MESH_TYPE_ICYL))
            {
              /* Create mesh */

              loc1 = CreateMesh(type, MESH_CONTENT_PTR, -1, nx, ny, nz,
                                mlims, -1);

              /* Free temporary arrays */

              Mem(MEM_FREE, mlims);
            }
          else if (type > 0)
            {
              /* Create mesh */

              loc1 = CreateMesh(type, MESH_CONTENT_PTR, -1, nx, ny, nz,
                                lims, -1);
            }

          /* Put pointer */

          WDB[loc0 + RMX_PTR_MESH] = (double)loc1;

          /* Create structures */

          idx = 0;

          for (k = 0; k < nz; k++)
            for (m = 0; m < ny; m++)
              for (n = 0; n < nx; n++)
                {
                  /* Create structure */

                  loc2 = NewItem(loc0 + RMX_PTR_MESH_DATA, RMX_CELL_BLOCK_SIZE);

                  /* Put pointer */

                  ptr = ReadMeshPtr(loc1, n, m, k);
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  WDB[ptr] = (double)loc2;

                  /* Put index and pointer */

                  WDB[loc2 + RMX_CELL_IDX] = (double)(idx);
                  WDB[loc2 + RMX_CELL_MESH_IDX] = (double)(idx++);
                  WDB[loc2 + RMX_CELL_PTR_MESH] = (double)loc1;
                }

          /* Read detectors */

          if (((ne = np - j) < 2) &&
              ((long)RDB[loc0 + RMX_MODE] != RMX_MODE_GVR))
            Error(loc0, "Missing responses");

          /* Loop over responses */

          while(j < np)
            {
              /* Allocate memory for data */

              ptr = NewItem(loc0 + RMX_PTR_DET, RMX_DET_BLOCK_SIZE);

              /* Check */

              if (np - j < 2)
                Error(loc0, "Missing response weight factor");
              else
                {
                  /* Read detector name */

                  WDB[ptr + RMX_DET_PTR_DET] = (double)PutText(params[j++]);

                  /* Read weight */

                  WDB[ptr + RMX_DET_WGT] =
                    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                              0.0, INFTY);
                }
            }

          /* Set flag */

          WDB[DATA_RMTX_CALC] = (double)YES;
        }

      /***********************************************************************/

      /***** Mass flow structure *********************************************/

      else if (!strcasecmp(word, "mflow"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 3, 1E+6, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_MFLOW0, MFLOW_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Put id */

          WDB[loc0 + MFLOW_PTR_NAME] = (double)PutText(params[j++]);

          /* Loop over values */

          while (j < np)
            {
              /* Create item */

              loc1 = NewItem(loc0 + MFLOW_PTR_DATA, MFLOW_LIST_BLOCK_SIZE);

              /* Put material id */

              WDB[loc1 + MFLOW_LIST_ZAI] = (double)PutText(params[j++]);

              /* Put rate */

              if (j < np)
                WDB[loc1 + MFLOW_LIST_RATE] =
                  TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                            -INFTY, INFTY);
              else
                Error(loc0, "Missing coefficient");
            }
        }

      /***********************************************************************/

      /***** Time bin structure **********************************************/

      else if (!strcasecmp(word, "tme"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 3, 1E+6, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_TME0, TME_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Grid name */

          WDB[loc0 + TME_PTR_NAME] = (double)PutText(params[j++]);

          /* Get type */

          type = (long)TestParam(pname, fname, line, params[j++],
                                 PTYPE_INT, 1, 3);

          /* Put type */

          WDB[loc0 + TME_TYPE] = (double)type;

          /* Read data */

          switch(type)
            {
            case TB_TYPE_ARB:
              {
                /* Arbitrary binning. Get number of values */

                ne = np - j;

                /* Set number of bins */

                WDB[loc0 + TME_NB] = (double)(ne - 1);

                /* Allocate memory for data */

                loc1 = ReallocMem(DATA_ARRAY, ne);

                /* Set pointer */

                WDB[loc0 + TME_PTR_BINS] = (double)loc1;

                /* Read values */

                while(j < np)
                  WDB[loc1++] = TestParam(pname, fname, line, params[j++],
                                           PTYPE_REAL, 0.0, 1E+9);

                /* Break case */

                break;
              }
            case TB_TYPE_UNI_T:
            case TB_TYPE_UNI_LOGT:
              {
                /* Uniform grid. Check number of parameters */

                if (np != 5)
                  Error(loc0, "Uniform binning %s must be defined by 5 values",
                        GetText(loc0 + TME_PTR_NAME));

                /* Read number of time bins */

                WDB[loc0 + TME_NB] =
                  TestParam(pname, fname, line, params[j++], PTYPE_INT,
                            1, 1000000);

                /* Read minimum and maximum time */

                WDB[loc0 + TME_TMIN] =
                  TestParam(pname, fname, line, params[j++],
                            PTYPE_REAL, 0.0, 1E+9);

                WDB[loc0 + TME_TMAX] =
                  TestParam(pname, fname, line, params[j++],
                            PTYPE_REAL, RDB[loc0 + TME_TMIN], 1E+9);

                /* Break case */

                break;
              }
            default:
              Error(loc0, "Invalid time binning type %d", type);
            }
         }

      /***********************************************************************/

      /***** User-defined response function **********************************/

      else if (!strcasecmp(word, "fun"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 3, 1E+6, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_FUN0, FUN_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Function name */

          WDB[loc0 + FUN_PTR_NAME] = (double)PutText(params[j++]);

          /* Get type */

          type = (long)TestParam(pname, fname, line, params[j++],
                                 PTYPE_INT, 1, 1);

          /* Put type */

          WDB[loc0 + FUN_TYPE] = (double)type;

          /* Read data */

          switch(type)
            {
            case 1:
              {
                /* Tabular interpolated data, get interpolation */

                WDB[loc0 + FUN_INT] =
                  TestParam(pname, fname, line, params[j++], PTYPE_INT, 1, 5);

                /* Check number of points */

                if ((np - j) % 2)
                  Error(loc0, "Odd number of values entered");

                /* Put number of points */

                ne = (np - j)/2;
                WDB[loc0 + FUN_NE] = (double)ne;

                /* Allocate memory for data */

                loc1 = ReallocMem(DATA_ARRAY, ne);
                WDB[loc0 + FUN_PTR_E] = (double)loc1;

                loc2 = ReallocMem(DATA_ARRAY, ne);
                WDB[loc0 + FUN_PTR_F] = (double)loc2;

                /* Read values */

                for (n = 0; n < ne; n++)
                  {
                    WDB[loc1++] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_REAL, 0.0, 1E+9);

                    WDB[loc2++] =
                      TestParam(pname, fname, line, params[j++],
                                PTYPE_REAL, -INFTY, INFTY);
                  }

                /* Break case */

                break;
              }
            default:
              Error(loc0, "Invalid response function type %d", type);
            }
         }

      /***********************************************************************/

      /***** RIA simulation **************************************************/

      else if (!strcasecmp(word, "ria"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 2, 2, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_RIA0, TME_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Name */

          WDB[loc0 + RIA_PTR_NAME] = (double)PutText(params[j++]);

          /* Get time binning */

          WDB[loc0 + RIA_PTR_TME] = (double)PutText(params[j++]);
         }

      /***********************************************************************/

      /***** Sensitivity calculations ****************************************/

      else if (!strcasecmp(word, "sens"))
        {
          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 100000000000, fname);

          /* Reset counter */

          j = 0;

          /* Reset number of parameters */

          k = -1;

          /* Create the sensitivity block if it doesn't already exist */

          if ((loc0 = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
            {
              loc0 = NewItem(DATA_PTR_SENS0, SENS_BLOCK_SIZE);

              /* Set default options */

              WDB[loc0 + SENS_N_MU] = 1.0;
              WDB[loc0 + SENS_N_SPAT] = 1.0;
              SetOption(loc0 + SENS_SCORE_FLAGS, SENS_SCORE_FLAG_ENE_TOT);
            }

          /* Compare name */

          if (!strcasecmp(params[j], "resp"))
            {
              /***************************************************************/

              /****** Responses for sensitivity calculation ******************/

              /* Set mode */

              WDB[DATA_SENS_MODE] = (double)SENS_MODE_BASIC;

              /* Copy parameter name */

              strcpy(pname, params[j]);

              k = j + 1;

              /* Check parameter type */

              strcpy(str, params[k++]);

              if (!strcasecmp(str, "keff"))
                {

                  /* Sensitivity of multiplication factor */

                  SetOption(loc0 + SENS_RESP_FLAGS,
                            SENS_RESP_FLAG_KEFF);

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++],
                                       PTYPE_LOGICAL);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_RESP_FLAGS,
                                  SENS_RESP_FLAG_KEFF);
                      else
                        ResetOption(loc0 + SENS_RESP_FLAGS,
                                    SENS_RESP_FLAG_KEFF);
                    }

                  /* Test if flag was set */

                  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_KEFF)
                    {
                      /* Create new repsonse block */

                      loc0 = NewItem(loc0 + SENS_PTR_RESP0, SENS_RESP_BLOCK_SIZE);

                      /* Put name, file name and line number */

                      WDB[loc0 + PARAM_PTR_NAME] =
                        (double)PutText("sens resp keff");
                      WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
                      WDB[loc0 + PARAM_LINE] = (double)line;

                      /* Put response type */

                      WDB[loc0 + SENS_RESP_TYPE] = (double)SENS_RESP_TYPE_KEFF;

                      /* Put response name */

                      WDB[loc0 + SENS_RESP_PTR_NAME] = PutText("KEFF");
                    }
                }
              else if (!strcasecmp(str, "leff"))
                {
                  /* Sensitivity of prompt generation time */

                  SetOption(loc0 + SENS_RESP_FLAGS,
                            SENS_RESP_FLAG_LEFF);

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++],
                                       PTYPE_LOGICAL);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_RESP_FLAGS,
                                  SENS_RESP_FLAG_LEFF);
                      else
                        ResetOption(loc0 + SENS_RESP_FLAGS,
                                    SENS_RESP_FLAG_LEFF);
                    }

                  /* Test if flag was set */

                  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_LEFF)
                    {
                      /* Create new repsonse block */

                      loc0 = NewItem(loc0 + SENS_PTR_RESP0, SENS_RESP_BLOCK_SIZE);

                      /* Put name, file name and line number */

                      WDB[loc0 + PARAM_PTR_NAME] =
                        (double)PutText("sens resp leff");
                      WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
                      WDB[loc0 + PARAM_LINE] = (double)line;

                      /* Put response type */

                      WDB[loc0 + SENS_RESP_TYPE] = (double)SENS_RESP_TYPE_LEFF;

                      /* Put response name */

                      WDB[loc0 + SENS_RESP_PTR_NAME] = PutText("LEFF");
                    }
                }
              else if (!strcasecmp(str, "beff"))
                {
                  /* Sensitivity of delayed neutron fraction */

                  SetOption(loc0 + SENS_RESP_FLAGS,
                            SENS_RESP_FLAG_BEFF);

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++],
                                       PTYPE_LOGICAL);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_RESP_FLAGS,
                                  SENS_RESP_FLAG_BEFF);
                      else
                        ResetOption(loc0 + SENS_RESP_FLAGS,
                                    SENS_RESP_FLAG_BEFF);
                    }

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++],
                                       PTYPE_LOGICAL);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_RESP_FLAGS,
                                  SENS_RESP_FLAG_BEFF_GROUPS);
                      else
                        ResetOption(loc0 + SENS_RESP_FLAGS,
                                    SENS_RESP_FLAG_BEFF_GROUPS);
                    }

                  /* Test if total flag was set */

                  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_BEFF)
                    {
                      /* Create new repsonse block */

                      loc1 = NewItem(loc0 + SENS_PTR_RESP0, SENS_RESP_BLOCK_SIZE);

                      /* Put name, file name and line number */

                      WDB[loc1 + PARAM_PTR_NAME] =
                        (double)PutText("sens resp beff");
                      WDB[loc1 + PARAM_PTR_FNAME] = (double)PutText(fname);
                      WDB[loc1 + PARAM_LINE] = (double)line;

                      /* Put response type */

                      WDB[loc1 + SENS_RESP_TYPE] = (double)SENS_RESP_TYPE_BEFF;

                      /* Put response name */

                      WDB[loc1 + SENS_RESP_PTR_NAME] = PutText("BEFF");
                    }

                  /* Test if group-wise flag was set */

                  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_BEFF_GROUPS)
                    {
                      for (i = 0; i < 8; i++)
                        {
                          /* Create new repsonse block */

                          loc1 = NewItem(loc0 + SENS_PTR_RESP0, SENS_RESP_BLOCK_SIZE);

                          /* Put name, file name and line number */

                          WDB[loc1 + PARAM_PTR_NAME] =
                            (double)PutText("sens resp beff");
                          WDB[loc1 + PARAM_PTR_FNAME] = (double)PutText(fname);
                          WDB[loc1 + PARAM_LINE] = (double)line;

                          /* Put response type */

                          WDB[loc1 + SENS_RESP_TYPE] = (double)SENS_RESP_TYPE_BEFF;

                          /* Put response group index */

                          WDB[loc1 + SENS_RESP_DET_BIN_IDX] = (double)(i+1);

                          /* Print response name */

                          sprintf(tmpstr, "BEFF_GROUP_%ld", i);

                          /* Put response name */

                          WDB[loc1 + SENS_RESP_PTR_NAME] = PutText(tmpstr);
                        }
                    }
                }
              else if (!strcasecmp(str, "lambda"))
                {
                  /* Sensitivity of delayed neutron fraction */

                  SetOption(loc0 + SENS_RESP_FLAGS,
                            SENS_RESP_FLAG_LAMBDA);

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++],
                                       PTYPE_LOGICAL);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_RESP_FLAGS,
                                  SENS_RESP_FLAG_LAMBDA);
                      else
                        ResetOption(loc0 + SENS_RESP_FLAGS,
                                    SENS_RESP_FLAG_LAMBDA);
                    }

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++],
                                       PTYPE_LOGICAL);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_RESP_FLAGS,
                                  SENS_RESP_FLAG_LAMBDA_GROUPS);
                      else
                        ResetOption(loc0 + SENS_RESP_FLAGS,
                                    SENS_RESP_FLAG_LAMBDA_GROUPS);
                    }

                  /* Test if total flag was set */

                  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_LAMBDA)
                    {
                      /* Create new repsonse block */

                      loc1 = NewItem(loc0 + SENS_PTR_RESP0, SENS_RESP_BLOCK_SIZE);

                      /* Put name, file name and line number */

                      WDB[loc1 + PARAM_PTR_NAME] =
                        (double)PutText("sens resp lambda");
                      WDB[loc1 + PARAM_PTR_FNAME] = (double)PutText(fname);
                      WDB[loc1 + PARAM_LINE] = (double)line;

                      /* Put response type */

                      WDB[loc1 + SENS_RESP_TYPE] = (double)SENS_RESP_TYPE_LAMBDA;

                      /* Put response name */

                      WDB[loc1 + SENS_RESP_PTR_NAME] = PutText("LAMBDA");
                    }

                  /* Test if group-wise flag was set */

                  if ((long)RDB[loc0 + SENS_RESP_FLAGS] & SENS_RESP_FLAG_LAMBDA_GROUPS)
                    {
                      for (i = 0; i < 8; i++)
                        {
                          /* Create new repsonse block */

                          loc1 = NewItem(loc0 + SENS_PTR_RESP0, SENS_RESP_BLOCK_SIZE);

                          /* Put name, file name and line number */

                          WDB[loc1 + PARAM_PTR_NAME] =
                            (double)PutText("sens resp lambda");
                          WDB[loc1 + PARAM_PTR_FNAME] = (double)PutText(fname);
                          WDB[loc1 + PARAM_LINE] = (double)line;

                          /* Put response type */

                          WDB[loc1 + SENS_RESP_TYPE] = (double)SENS_RESP_TYPE_LAMBDA;

                          /* Put response group index */

                          WDB[loc1 + SENS_RESP_DET_BIN_IDX] = (double)(i+1);

                          /* Print response name */

                          sprintf(tmpstr, "LAMBDA_GROUP_%ld", i);

                          /* Put response name */

                          WDB[loc1 + SENS_RESP_PTR_NAME] = PutText(tmpstr);
                        }
                    }
                }
              else if (!strcasecmp(str, "detratio"))
                {
                  /* Sensitivity of detector ratio */

                  /* Check that exactly three arguments are given */

                  if (np - k < 3)
                    Error(-1, pname, fname, line,
                          "Need at least 3 arguments for \"sens resp detratio\" "
                          "got %ld", np - k);

                  /* Create new detector ratio block */

                  loc0 = NewItem(loc0 + SENS_PTR_RESP0, SENS_RESP_BLOCK_SIZE);

                  /* Put name, file name and line number */

                  WDB[loc0 + PARAM_PTR_NAME] =
                    (double)PutText("sens resp detratio");
                  WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
                  WDB[loc0 + PARAM_LINE] = (double)line;

                  /* Put response type */

                  WDB[loc0 + SENS_RESP_TYPE] = (double)SENS_RESP_TYPE_RATIO;

                  /* Reset detector bin index */

                  WDB[loc0 + SENS_RESP_DET_BIN_IDX] = (double)-1;

                  /* Put response name */

                  WDB[loc0 + SENS_RESP_PTR_NAME] = PutText(params[k++]);

                  /* Put name of first detector */

                  WDB[loc0 + SENS_RESP_PTR_DET1] = PutText(params[k++]);

                  /* Put name of second detector */

                  WDB[loc0 + SENS_RESP_PTR_DET2] = PutText(params[k++]);

                  /* Put frac to score for first detector */

                  if (k < np)
                    WDB[loc0 + SENS_RESP_DET1_P] = TestParam(pname, fname, line, params[k++],
                                                             PTYPE_REAL, 0.0, 1.0);
                  else
                    WDB[loc0 + SENS_RESP_DET1_P] = 1.0;

                  /* Put frac to score for second detector */

                  if (k < np)
                    WDB[loc0 + SENS_RESP_DET2_P] = TestParam(pname, fname, line, params[k++],
                                                             PTYPE_REAL, 0.0, 1.0);
                  else
                    WDB[loc0 + SENS_RESP_DET2_P] = 1.0;

                }
              else if (!strcasecmp(str, "void"))
                {
                  /* Sensitivity of void fraction */

                  /* Check that exactly one argument is given */

                  if (np - k != 1)
                    Error(-1, pname, fname, line,
                          "Need exactly 1 arguments for \"sens resp void\" "
                          "got %ld", np - k);

                  /* Put void fraction material name */

                  WDB[loc0 + SENS_PTR_VOID_MAT] = PutText(params[k++]);

                  /* Set flag for calculating the void coef sensitivity */

                  SetOption(loc0 + SENS_RESP_FLAGS, SENS_RESP_FLAG_VOID);

                  /* Create new repsonse block */

                  loc0 = NewItem(loc0 + SENS_PTR_RESP0, SENS_RESP_BLOCK_SIZE);

                  /* Put name, file name and line number */

                  WDB[loc0 + PARAM_PTR_NAME] =
                    (double)PutText("sens resp void");
                  WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
                  WDB[loc0 + PARAM_LINE] = (double)line;

                  /* Put response type */

                  WDB[loc0 + SENS_RESP_TYPE] = (double)SENS_RESP_TYPE_VOID;

                  /* Put response name */

                  WDB[loc0 + SENS_RESP_PTR_NAME] = PutText("VOID");
                }
              else
                Error(-1, params[j], fname, line,
                      "Unknown parameter \"%s\" for \"sens resp\"", str);

              /* Set flags for recording events */

              SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_IFP);
              SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_SENS);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "pert"))
            {
              /****** Non-default sources for sensitivity calculation ********/

              /* Set mode */

              WDB[DATA_SENS_MODE] = (double)SENS_MODE_BASIC;

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check parameter type */

              strcpy(str, params[k++]);

              if (!strcasecmp(str, "chi"))
                {
                  /* Fission spectrum*/

                  SetOption(loc0 + SENS_PERT_FLAGS,
                            SENS_PERT_FLAG_CHI);

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 1);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_PERT_FLAGS,
                                  SENS_PERT_FLAG_CHI);
                      else
                        ResetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_CHI);
                    }
                }
              else if (!strcasecmp(str, "temperature"))
                {
                  /* Fission spectrum*/

                  SetOption(loc0 + SENS_PERT_FLAGS,
                            SENS_PERT_FLAG_TEMPERATURE);

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 1);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_PERT_FLAGS,
                                  SENS_PERT_FLAG_TEMPERATURE);
                      else
                        ResetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_TEMPERATURE);
                    }
                }
              else if (!strcasecmp(str, "eleg"))
                {
                  /* Legendre moments of elastic scattering secondary distribution */

                  SetOption(loc0 + SENS_PERT_FLAGS,
                            SENS_PERT_FLAG_SCATT_MOM);

                  /* Set maximum legendre moment to calculate */

                  WDB[loc0 + SENS_MAX_SCATT_MOM] = (double)3;

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 7);

                      /* Set mode */

                      if (type > 0)
                        {
                          SetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_SCATT_MOM);

                          /* Set maximum legendre moment to calculate */

                          WDB[loc0 + SENS_MAX_SCATT_MOM] = (double)type;
                        }
                      else
                        ResetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_SCATT_MOM);
                    }
                }
              else if (!strcasecmp(str, "nubar"))
                {
                  /* Total, prompt and delayed nubar */

                  SetOption(loc0 + SENS_PERT_FLAGS,
                            SENS_PERT_FLAG_NUBAR);

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 1);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_PERT_FLAGS,
                                  SENS_PERT_FLAG_NUBAR);
                      else
                        ResetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_NUBAR);
                    }
                }
              else if (!strcasecmp(str, "elamu"))
                {
                  /* Total, prompt and delayed nubar */

                  SetOption(loc0 + SENS_PERT_FLAGS,
                            SENS_PERT_FLAG_ELA_MU);

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 1);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_PERT_FLAGS,
                                  SENS_PERT_FLAG_ELA_MU);
                      else
                        ResetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_ELA_MU);
                    }

                }
              else if (!strcasecmp(str, "inlmu"))
                {
                  /* Inelastic scattering cosine */

                  SetOption(loc0 + SENS_PERT_FLAGS,
                            SENS_PERT_FLAG_INL_MU);

                  if (k < np)
                    {
                      /* Get mode */

                      type = TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 1);

                      /* Set mode */

                      if (type == 1)
                        SetOption(loc0 + SENS_PERT_FLAGS,
                                  SENS_PERT_FLAG_INL_MU);
                      else
                        ResetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_INL_MU);
                    }
                }
              else if (!strcasecmp(str, "xs"))
                {
                  /* Default nuclear cross section perturbations */

                  SetOption(loc0 + SENS_PERT_FLAGS,
                            SENS_PERT_FLAG_XS);

                  if (k < np)
                    {

                      /* Check parameter name */

                      strcpy(str, params[k++]);

                      if ((!strcasecmp(str, "all")) ||
                          (!strcasecmp(str, "allrea")))
                        {
                          /* All sum reaction modes */

                          SetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_XS);
                          ResetOption(loc0 + SENS_PERT_FLAGS,
                                      SENS_PERT_FLAG_XSMT);
                        }
                      else if (!strcasecmp(str, "allmt"))
                        {
                          /* All mt numbers */

                          SetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_XSMT);
                          ResetOption(loc0 + SENS_PERT_FLAGS,
                                      SENS_PERT_FLAG_XS);
                        }
                      else if (!strcasecmp(str, "realist"))
                        {
                          /* Set flag to identify reactions based on sum reaction mode */

                          SetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_XS);
                          ResetOption(loc0 + SENS_PERT_FLAGS,
                                      SENS_PERT_FLAG_XSMT);

                          /* Count parameters before control word */

                          n = CountSensParams(params + k, np - k);

                          /* Allocate memory for list */

                          ptr = ReallocMem(DATA_ARRAY, n + 1);

                          /* Store pointer to list */

                          WDB[loc0 + SENS_PTR_REA_LIST] = (double)ptr;

                          /* Read values */

                          for (; n > 0; n--)
                            WDB[ptr++] = PutText(params[k++]);

                          /* Put null termintor */

                          WDB[ptr] = NULLPTR;

                        }
                      else if (!strcasecmp(str, "mtlist"))
                        {
                          /* Set flag to identify reactions based on mt numbers */

                          SetOption(loc0 + SENS_PERT_FLAGS,
                                    SENS_PERT_FLAG_XSMT);
                          ResetOption(loc0 + SENS_PERT_FLAGS,
                                      SENS_PERT_FLAG_XS);

                          /* List of MT numbers */

                          /* Count parameters before control word */

                          n = CountSensParams(params + k, np - k);

                          /* Allocate memory for list */

                          ptr = ReallocMem(DATA_ARRAY, n + 1);

                          /* Store pointer to list */

                          WDB[loc0 + SENS_PTR_MT_LIST] = (double)ptr;

                          /* Read values */

                          for (; n > 0; n--)
                            WDB[ptr++] = TestParam(pname, fname, line, params[k++],
                                                   PTYPE_INT, 1, 1004);

                          /* Put null terminator */

                          WDB[ptr] = NULLPTR;

                        }
                      else if (!strcasecmp(str, "0"))
                        {
                          /* Unset the option */

                          ResetOption(loc0 + SENS_PERT_FLAGS,
                                      SENS_PERT_FLAG_XS);
                          ResetOption(loc0 + SENS_PERT_FLAGS,
                                      SENS_PERT_FLAG_XSMT);
                        }
                      else
                        {
                          /* XS based on MT number */

                          Error(-1, pname, fname, line,
                                "Unknown option \"%s\" for "
                                "\"sens pert xs\"", str);
                        }
                    }
                }
              else if (!strcasecmp(str, "zailist"))
                {
                  /* ZAIs to include in the perturbation */

                  /* Check that ZAI list has not been created yet */

                  if ((long)RDB[loc0 + SENS_PTR_ZAI_ARR] != 0)
                    Error(-1, pname, fname, line,
                          "ZAI array for sensitivity calculations has already been created "
                          "using \"sens pert zailist\".");

                  if (k < np)
                    {

                      /* Check parameter name */

                      strcpy(str, params[k]);

                      if (!strcasecmp(str, "all"))
                        {

                          /* Include all ZAIs */

                          WDB[loc0 + SENS_PTR_ZAI_ARR] = (double)-1;

                          SetOption(loc0 + SENS_SCORE_FLAGS, SENS_SCORE_FLAG_ZAI_TOT);

                        }
                      else
                        {
                          /* List of ZAIs */

                          /* Get number of ZAIs */

                          n = np - k;

                          /* Store number of ZAIs */

                          WDB[loc0 + SENS_N_ZAI] = (double)n;

                          /* Allocate memory for ZAI list */

                          ptr = ReallocMem(DATA_ARRAY, n + 1);

                          /* Store pointer to ZAI list */

                          WDB[loc0 + SENS_PTR_ZAI_ARR] = (double)ptr;

                          /* Read ZAIs */

                          for (; n > 0; n--)
                            {
                              if (!strcasecmp(params[k],"total"))
                                {
                                  /* Check for catch-all ZAI (total) */

                                  SetOption(loc0 + SENS_SCORE_FLAGS, SENS_SCORE_FLAG_ZAI_TOT);
                                  k++;

                                  /* Reduce number of total ZAI */

                                  WDB[loc0 + SENS_N_ZAI] -= 1.0;
                                }
                              else if (!strcasecmp(params[k],"sum"))
                                {
                                  /* Check for sum ZAI */

                                  SetOption(loc0 + SENS_SCORE_FLAGS, SENS_SCORE_FLAG_ZAI_SUM);
                                  k++;

                                  /* Reduce number of total ZAI */

                                  WDB[loc0 + SENS_N_ZAI] -= 1.0;
                                }
                              else
                                WDB[ptr++] = TestParam(pname, fname, line, params[k++],
                                                       PTYPE_INT, 10010, 999991);
                            }

                          /* Put null pointer */

                          WDB[ptr] = NULLPTR;
                        }
                    }
                }
              else if (!strcasecmp(str, "matlist"))
                {
                  /* Materials to include in the perturbation */

                  /* Check that material list has not been created yet */

                  if ((long)RDB[loc0 + SENS_PTR_MAT_ARR] != 0)
                    Error(-1, pname, fname, line,
                          "Material array for sensitivity calculations has already been created "
                          "using \"sens pert matlist\".");

                  if (k < np)
                    {
                      /* Check parameter name */

                      strcpy(str, params[k]);

                      if (!strcasecmp(str, "all"))
                        {

                          /* Include all materials */

                          WDB[loc0 + SENS_PTR_MAT_ARR] = (double)-1;

                          SetOption(loc0 + SENS_SCORE_FLAGS, SENS_SCORE_FLAG_MAT_TOT);

                        }
                      else
                        {
                          /* List of materials */

                          /* Get number of materials */

                          n = np - k;

                          /* Store number of materials */

                          WDB[loc0 + SENS_N_MAT] = (double)n;

                          /* Allocate memory for material list */

                          ptr = ReallocMem(DATA_ARRAY, n + 1);

                          /* Store pointer to materials list */

                          WDB[loc0 + SENS_PTR_MAT_ARR] = (double)ptr;

                          /* Read materials (linked in ProcessSensZAIs) */

                          for (; n > 0; n--)
                            {
                              if (!strcasecmp(params[k], "total"))
                                {
                                  /* Check for catch-all material (total) */

                                  SetOption(loc0 + SENS_SCORE_FLAGS, SENS_SCORE_FLAG_MAT_TOT);
                                  k++;

                                  /* Decrease number of materials included */

                                  WDB[loc0 + SENS_N_MAT] -= 1.0;
                                }
                              else if (!strcasecmp(params[k], "sum"))
                                {
                                  /* Check for sum material */

                                  SetOption(loc0 + SENS_SCORE_FLAGS, SENS_SCORE_FLAG_MAT_SUM);
                                  k++;

                                  /* Decrease number of materials included */

                                  WDB[loc0 + SENS_N_MAT] -= 1.0;
                                }
                              else

                                WDB[ptr++] = PutText(params[k++]);
                            }

                          /* Put null pointer */

                          WDB[ptr] = NULLPTR;
                        }
                    }
                }
              else if (!strcasecmp(str, "custom"))
                {
                  /* Custom perturbation (can be xGPT) */

                  if (np - k < 1)
                    Error(-1, pname, fname, line,
                          "Need at least name for \"sens resp custom\" ");

                  /* Create new custom perturbation block */

                  loc1 = NewItem(loc0 + SENS_PTR_PERT0, SENS_PERT_BLOCK_SIZE);

                  /* Put name, file name and line number */

                  WDB[loc1 + PARAM_PTR_NAME] = (double)PutText("sens pert custom");
                  WDB[loc1 + PARAM_PTR_FNAME] = (double)PutText(fname);
                  WDB[loc1 + PARAM_LINE] = (double)line;

                  /* Put perturbation type */
                  /*
                  WDB[loc1 + SENS_PERT_TYPE] = 0.0;
                  */
                  /* Put response name */

                  WDB[loc1 + SENS_PERT_PTR_NAME] = PutText(params[k++]);

                  while (k < np)
                    {
                      /* Check parameter name */

                      strcpy(str, params[k++]);

                      if (!strcasecmp(str, "efunc"))
                        {
                          /* Energy dependent function */

                          /* Get name of file */

                          WDB[loc1 + SENS_PERT_PTR_EGRID] = PutText(params[k++]);
                        }
                      else if (!strcasecmp(str, "zailist"))
                        {
                          /* List of ZAI numbers */

                          /* Count parameters before control word */

                          n = CountSensParams(params + k, np - k);

                          /* Allocate memory for ZAI list */

                          ptr = ReallocMem(DATA_ARRAY, n + 1);

                          /* Store pointer to ZAI list */

                          WDB[loc1 + SENS_PERT_PTR_ZAI_LIST] = (double)ptr;

                          /* Read ZAIs */

                          for (; n > 0; n--)
                            WDB[ptr++] = TestParam(pname, fname, line, params[k++],
                                                   PTYPE_INT, 10010, 999991);

                          /* Put null pointer */

                          WDB[ptr] = NULLPTR;
                        }
                      else if (!strcasecmp(str, "matlist"))
                        {
                          /* List of material names */

                          /* Count parameters before control word */

                          n = CountSensParams(params + k, np - k);

                          /* Allocate memory for list */

                          ptr = ReallocMem(DATA_ARRAY, n + 1);

                          /* Store pointer to list */

                          WDB[loc1 + SENS_PERT_PTR_MAT_LIST] = (double)ptr;

                          /* Read values */

                          for (; n > 0; n--)
                            WDB[ptr++] = PutText(params[k++]);

                          /* Put null pointer */

                          WDB[ptr] = NULLPTR;
                        }
                      else if (!strcasecmp(str, "realist"))
                        {
                          /* List of reaction names */

                          /* Count parameters before control word */

                          n = CountSensParams(params + k, np - k);

                          /* Allocate memory for list */

                          ptr = ReallocMem(DATA_ARRAY, n + 1);

                          /* Store pointer to list */

                          WDB[loc1 + SENS_PERT_PTR_REA_LIST] = (double)ptr;

                          /* Read values */

                          for (; n > 0; n--)
                            WDB[ptr++] = PutText(params[k++]);

                          /* Put null termintor */

                          WDB[ptr] = NULLPTR;
                        }
                      else if (!strcasecmp(str, "mtlist"))
                        {
                          /* List of MT numbers */

                          /* Count parameters before control word */

                          n = CountSensParams(params + k, np - k);

                          /* Allocate memory for list */

                          ptr = ReallocMem(DATA_ARRAY, n + 1);

                          /* Store pointer to list */

                          WDB[loc1 + SENS_PERT_PTR_MT_LIST] = (double)ptr;

                          /* Read values */

                          for (; n > 0; n--)
                            WDB[ptr++] = TestParam(pname, fname, line, params[k++],
                                                   PTYPE_INT, 1, 1004);

                          /* Put null terminator */

                          WDB[ptr] = NULLPTR;
                        }
                      else
                        {
                          Error(-1, params[j], fname, line,
                                "Unknown parameter \"%s\" for \"sens pert custom\"", str);
                        }
                    }
                }
              else
                Error(-1, params[j], fname, line,
                      "Unknown parameter \"%s\" for \"sens pert\"", str);

              /* Set flags for recording events */

              SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_IFP);
              SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_SENS);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "opt"))
            {
              /****** Options for sensitivity calculation ********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Set mode */

              WDB[DATA_SENS_MODE] = (double)SENS_MODE_BASIC;

              /* Loop over remaining parameters */

              /* Check parameter type */

              strcpy(str, params[k++]);

              if (!strcasecmp(str, "egrid"))
                {
                  /* Sens energy grid. Check number of parameters. */

                  if (k == np)
                    Error(-1, params[j], fname, line,
                          "Missing grid name after \"sens opt egrid\"");

                  /* Read energy grid name */

                  WDB[loc0 + SENS_PTR_EGRID] = (double)PutText(params[k++]);

                }
              else if (!strcasecmp(str, "dmesh"))
                {
                  /* Sens spatial mesh. Check number of parameters. */

                  if (k == np)
                    Error(-1, params[j], fname, line,
                          "Missing mesh name after \"sens opt datamesh\"");

                  /* Read energy grid name */

                  WDB[loc0 + SENS_PTR_MESH] = (double)PutText(params[k++]);

                }
              else if (!strcasecmp(str, "history"))
                {
                  /* Set option on by default is flag is not given */

                  SetOption(loc0 + SENS_RESP_FLAGS, SENS_SCORE_FLAG_HIS);

                  /* Get flag */

                  if (k < np)
                    {
                      n = TestParam(pname, fname, line, params[k++],
                                    PTYPE_LOGICAL);

                        if (n)
                          SetOption(loc0 + SENS_RESP_FLAGS, SENS_SCORE_FLAG_HIS);
                        else
                          ResetOption(loc0 + SENS_RESP_FLAGS, SENS_SCORE_FLAG_HIS);
                    }

                }
              else if (!strcasecmp(str, "latgen"))
                {
                  /* number of latent generations */

                  if (k == np)
                    Error(-1, params[j], fname, line,
                          "Missing number of latent generations after \"sens opt latgen\"");

                  /* Get number of generations */

                  n = TestParam(pname, fname, line, params[k++],
                                PTYPE_INT, 0, 100);

                  /* Put number of generations */

                  WDB[DATA_SENS_LAST_GEN] = (double)(n + 1);

                }
              else if (!strcasecmp(str, "direct"))
                {
                  /* Store events instead of direct score */

                  if (k == np)
                    Error(-1, params[j], fname, line,
                          "Missing value for \"sens opt direct\"");

                  /* Store type */

                  WDB[DATA_SENS_SCORE_TYPE] = SENS_SCORE_TYPE_DIRECT;

                  /* Get event block bank size */

                  WDB[DATA_EBLOCK_BANK_SZ] = TestParam(pname, fname, line, params[k++],
                                 PTYPE_REAL, 1e-6, INFTY);

                }
              else
                Error(-1, params[j], fname, line,
                      "Unknown parameter \"%s\" for \"sens opt\"", str);

              /* Set flags for recording events */

              SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_IFP);
              SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_SENS);
            }
          else
            {
              /***** Parameter is not defined ********************************/

              Error(-1, params[j], fname, line, "Unidentified option");

              /* NOTE: Jos toi errori poistetaan niin kaikkia parametreja ei */
              /*       lueta --> aseta k = np jotta ajo ei kaadu alempana    */
              /*       olevaan testaukseen.  */

              k = np;

              /***************************************************************/
            }
        }

      /***** Data sample  ****************************************************/

      /***** Currently samples temperature and density *****/

      else if (!strcasecmp(word, "sample"))
        {
          /* Copy parameter name */

          strcpy (pname, word);

          /* Read parameters (mesh information) */

          params = GetParams(word, input, &np, &i0, 9, 9, fname);

          /* Create new item */

          loc0 = NewItem(DATA_PTR_SAMPLE0, SAMPLE_BLOCK_SIZE);

          /* Put name, file name and line number */

          WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
          WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
          WDB[loc0 + PARAM_LINE] = (double)line;

          /* Read data */

          j = 0;

          /* Read x-binning */

          WDB[loc0 + SAMPLE_MESH_N0] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT,
                      1, 1000000);
          WDB[loc0 + SAMPLE_MESH_MIN0] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      -INFTY, INFTY);
          WDB[loc0 + SAMPLE_MESH_MAX0] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      RDB[loc0 + SAMPLE_MESH_MIN0], INFTY);

          /* Read y-binning */

          WDB[loc0 + SAMPLE_MESH_N1] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT,
                      1, 1000000);
          WDB[loc0 + SAMPLE_MESH_MIN1] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      -INFTY, INFTY);
          WDB[loc0 + SAMPLE_MESH_MAX1] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      RDB[loc0 + SAMPLE_MESH_MIN1], INFTY);

          /* Read z-binning */

          WDB[loc0 + SAMPLE_MESH_N2] =
            TestParam(pname, fname, line, params[j++], PTYPE_INT,
                      1, 1000000);
          WDB[loc0 + SAMPLE_MESH_MIN2] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      -INFTY, INFTY);
          WDB[loc0 + SAMPLE_MESH_MAX2] =
            TestParam(pname, fname, line, params[j++], PTYPE_REAL,
                      RDB[loc0 + SAMPLE_MESH_MIN2], INFTY);
         }

      /***********************************************************************/

      /***** Misc. parameters ************************************************/

      else if (!strcasecmp(word, "set"))
        {
          /* Read parameters */

          params = GetParams(word, input, &np, &i0, 1, 100000000000, fname);

          /* Reset counter */

          j = 0;

          /* Reset number of parameters */

          k = -1;

          /* Compare name */

          if (!strcasecmp(params[j], "dt"))
            {
              /****** Delta-tracking parameters ******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get probability threshold for neutron transport */

              WDB[DATA_DT_NTHRESH] =
                1.0 - TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                0.0, 1.0);

              /* Get probability threshold for photon transport */

              if (k < np)
                WDB[DATA_DT_PTHRESH] =
                  1.0 - TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                  0.0, 1.0);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "spd"))
            {
              /****** Override particle speeds *******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Neutron speed */

              if (k < np)
                WDB[DATA_NEUTRON_SPD] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, INFTY);
              else
                Error(-1, pname, fname, line, "Missing neutron speed");

              /* Photon speed */

              if (k < np)
                WDB[DATA_PHOTON_SPD] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, INFTY);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "bala"))
            {
              /****** OMP load balancing options *****************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Probability to store particles in common que */

              if (k < np)
                WDB[DATA_COMMON_QUE_LIM] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            -1, 1000000000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "stlfile"))
            {
              /****** Option to write STL search mesh in file ****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read option */

              if (k < np)
                WDB[DATA_STL_MESH_FILES] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "dd"))
            {
              /****** Domain decomposition mode ******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get mode */

              if (k < np)
                WDB[DATA_DD_MODE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 3);

              /* Check mode */

              if (((long)RDB[DATA_DD_MODE] == DD_MODE_AUTO) ||
                  ((long)RDB[DATA_DD_MODE] == DD_MODE_SECTOR))
                {
                  /* Get origin */

                  if (k < np)
                    WDB[DATA_DD_ORIG_X0] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                -INFTY, INFTY);

                  if (k < np)
                    WDB[DATA_DD_ORIG_Y0] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                -INFTY, INFTY);

                  if (k < np)
                    WDB[DATA_DD_ORIG_Z0] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                -INFTY, INFTY);

                  if (k < np)
                    WDB[DATA_DD_SECT0] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                0.0, 360.0);
                }
              else if ((long)RDB[DATA_DD_MODE] == DD_MODE_LAT)
                Error(-1, pname, fname, line, "Mode not supported yet");

              /* Switch DD on */

              if ((long)RDB[DATA_DD_MODE] != DD_MODE_NONE)
                WDB[DATA_DD_DECOMPOSE] = (double)YES;

              /* Check that source code is compiled in MPI */

#ifndef MPI

              Error(-1, pname, fname, line,
                    "Domain decomposition requires MPI mode");

#endif

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "wwb"))
            {
              /****** Weight window bounds ***********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Lower bound */

              if (k < np)
                WDB[DATA_WWD_LOWER_BOUND] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 1.0);
              else
                Error(-1, pname, fname, line, "Missing lower bound");

              /* Upper bound */

              if (k < np)
                WDB[DATA_WWD_UPPER_BOUND] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            1.0, INFTY);
              else
                Error(-1, pname, fname, line, "Missing upper bound");

              /* Roulette survival factor */

              if (k < np)
                WDB[DATA_WWD_SURVIVAL_F] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            1.0, INFTY);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "keff"))
            {
              /****** Additional scaling of fission nubar ********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read k-eff (used as scaling factor for nubar) */

              if (k < np)
                WDB[DATA_EXT_SRC_NUBAR_F] =
                  1.0/TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                1E-3, 1E3);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "maxsplit"))
            {
              /****** Limints for splitting and roulette *********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read limits */

              if (k < np)
                WDB[DATA_WWD_MAX_SPLIT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            2, 1000000000);

              if (k < np)
                WDB[DATA_WWD_MIN_ROULETTE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            ZERO, 1.0);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "rmxtest"))
            {
              /****** Test modes for RMX solver ******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read limit */

              if (k < np)
                WDB[DATA_RMX_TEST_MODE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 3);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "rmx"))
            {
              /****** Convergence criteria for RMX solver ********************/

              /* NOTE: Ei paras nimi parametrille */

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read values */

              if (k < np)
                WDB[DATA_RMX_MULT_CONVG_LIM] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 1.0);

              if (k < np)
                WDB[DATA_RMX_MULT_CONVG_PASS] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 1.0);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "cslist"))
            {
              /****** Parameters for cell search lists ***********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read limit for creating cell list */

              if (k < np)
                WDB[DATA_MAX_CELL_SEARCH_LIST] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            0, 1000000000);

              /* Read number of points for pre-sorting */

              if (k < np)
                WDB[DATA_PRESORT_NP] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            0, 1000000000);

              /* Read number of batches for pre-sorting */

              if (k < np)
                WDB[DATA_PRESORT_NB] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            0, 1000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "dfsol"))
            {
              /****** Diffusion solver options *******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get solver mode */

              if (k < np)
                WDB[DATA_HOMOFLUX_SOLVER] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 2);

              /* Diffusion coefficient used in the calculation */

              if (k < np)
                WDB[DATA_HOMOFLUX_DIFFCOEF] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 2);

              /* Get number of points for trapezoidal integration */

              if (k < np)
                WDB[DATA_ADF_TRAPZ_PT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            10, 1000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "dhfit"))
            {
              /****** Exponential function fit for decay heat ****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of functions */

              if (k < np)
                WDB[DATA_EXPO_DEC_FIT_NF] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                            1000);
              else
                Error(-1, pname, fname, line, "Missing number of terms");

              /* Get number of time points*/

              if (k < np)
                WDB[DATA_EXPO_DEC_FIT_NT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            2, 1000000);

              /* Get time limits */

              if (k < np)
                WDB[DATA_EXPO_DEC_FIT_TMIN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, INFTY);

              if (k < np)
                WDB[DATA_EXPO_DEC_FIT_TMAX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            RDB[DATA_EXPO_DEC_FIT_TMIN], INFTY);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "roulette"))
            {
              /***** Parameters for russian roulette *************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get weight limit */

              WDB[DATA_OPT_ROULETTE_W0] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          0.0, 1.0);

              /* Get probability threshold */

              if (k < np)
                WDB[DATA_OPT_ROULETTE_P0] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 1.0);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "gct"))
            {
              /***** Perform statistical tests on group constants ************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Option */

              if (k < np)
                WDB[DATA_GC_STAT_TESTS] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Set global */

              if ((long)RDB[DATA_GC_STAT_TESTS] == YES)
                WDB[DATA_RUN_STAT_TESTS] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "rnddec"))
            {
              /***** Randomize decay data ************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Option for decay constants */

              if (k < np)
                WDB[DATA_BURN_RANDOMIZE_DEC] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Option for fission yields */

              if (k < np)
                WDB[DATA_BURN_RANDOMIZE_FY] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "confi"))
            {
              /***** Set confidentiality label *******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Option */

              if (k < np)
                WDB[DATA_CONFIDENTIAL] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "tpa"))
            {
              /***** Track plot animation ************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Minimum and maximum time */

              if (k < np)
                WDB[DATA_TRACK_PLOT_TMIN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, INFTY);

              if (k < np)
                WDB[DATA_TRACK_PLOT_TMAX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            RDB[DATA_TRACK_PLOT_TMIN], INFTY);

              /* Plotted history lenght (in cm) */

              if (k < np)
                WDB[DATA_TRACK_PLOT_HIS_LENGTH] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, INFTY);

              /* Number of frames */

              if (k < np)
                WDB[DATA_TRACK_PLOT_FRAMES] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                            100000);

              /* Event buffer */

              if (k < np)
                WDB[DATA_EVENT_BANK_SZ] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                            1000000000);

              /* Set mode */

              WDB[DATA_TRACK_PLOT_ANIM] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "spa"))
            {
              /***** Source point animation **********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read palette */

              if (k < np)
                WDB[DATA_SOURCE_PT_ANIM_PALETTE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 13);
              else
                Error(-1, pname, fname, line, "Missing palette");

              /* Read fraction for colour scheme */

              if (k < np)
                WDB[DATA_SOURCE_PT_ANIM_F] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 1.0);

              /* Set mode */

              WDB[DATA_SOURCE_PT_ANIM] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "rfw"))
            {
              /***** Write binary restart file *******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Option */

              if (k < np)
                WDB[DATA_WRITE_RESTART_FILE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);
              else
                Error(-1, pname, fname, line, "Missing option");

              /* Get file name */

              if (k < np)
                WDB[DATA_RESTART_WRITE_PTR_FNAME] =
                  (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "rfr"))
            {
              /***** Read binary restart file ********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check if step is given */

              if (!strcasecmp(params[k], "idx"))
                {
                  /* Update index */

                  k++;

                  /* Read step index */

                  if (k < np)
                    WDB[DATA_RESTART_READ_IDX] = 1.0 +
                      TestParam(pname, fname, line, params[k++], PTYPE_INT,
                                0, 10000);
                  else
                    Error(-1, pname, fname, line, "Missing step index");
                }
              else if (!strcasecmp(params[k], "continue"))
                {
                  /* Update index */

                  k++;

                  /* Continue from last */

                  WDB[DATA_RESTART_READ_CONTINUE] = (double)YES;
                }
              else if (!strcasecmp(params[k], "bu"))
                {
                  /* Update index */

                  k++;

                  /* Read step index */

                  if (k < np)
                    WDB[DATA_RESTART_READ_POINT] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                ZERO, 10000);
                  else
                    Error(-1, pname, fname, line, "Missing burnup");
                }
              else if (!strcasecmp(params[k], "days"))
                {
                  /* Update index */

                  k++;

                  /* Read step index */

                  if (k < np)
                    WDB[DATA_RESTART_READ_POINT] =
                      -TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                ZERO, 1E+12);
                  else
                    Error(-1, pname, fname, line, "Missing time");
                }
              else
                {
                  /* Burnup or days */

                  if (k < np)
                    WDB[DATA_RESTART_READ_POINT] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                -1E+12, 1000.0);
                  else
                    Error(-1, pname, fname, line, "Missing burnup or time");
                }

              /* Get file name */

              if (k < np)
                WDB[DATA_RESTART_READ_PTR_FNAME] =
                  (double)PutText(params[k++]);
              else if ((long)RDB[DATA_RESTART_READ_CONTINUE] == NO)
                Error(-1, pname, fname, line, "Missing file name");

              /* Set option */

              WDB[DATA_READ_RESTART_FILE] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "delnu"))
            {
              /***** Sample delayed neutrons *********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_USE_DELNU] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "fissye"))
            {
              /***** Include energy-dependent fission yields *****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_FISSY_ENE_DEP] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "sca"))
            {
              /***************************************************************/

              /***** Source convergence acceleration *************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Number of outer iterations */

              if (k < np)
                WDB[DATA_RMX_CONVG_OUT_N] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 1000);
              else
                Error(-1, pname, fname, line,
                      "Missing number of outer iterations");

              /* Number of source batches */

              if (k < np)
                WDB[DATA_RMX_CONVG_SB] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 1000);
              else
                Error(-1, pname, fname, line,
                      "Missing number of source batches");

              /* Check existing */

              if ((long)RDB[DATA_PTR_RMX0] > VALID_PTR)
                Error(-1, pname, fname, line,
                      "Weight-window generation not allowed");

              /* Create new response-matrix item */

              loc0 = NewItem(DATA_PTR_RMX0, RMX_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(pname);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Avoid compiler warning */

              lims[0] = 0.0;
              lims[1] = 0.0;
              lims[2] = 0.0;
              lims[3] = 0.0;
              lims[4] = 0.0;
              lims[5] = 0.0;

              nx = -1;
              ny = -1;
              nz = -1;

              type = -1;

              /* Set default tolerances */

              WDB[DATA_RMX_CONVG_OUT_TOL] = 1E-2;
              WDB[loc0 + RMX_CONV] = 1E-9;
              WDB[loc0 + RMX_N_ITER] = 10000.0;

              /* Read mesh type */

              if (k < np)
                type = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 7);
              else
                Error(loc0, "Missing mesh type");

              /* Check */

              if (type == MESH_TYPE_CARTESIAN)
                {
                  /* Check number of parameters */

                  if (k > np - 9)
                    Error(loc0, "Missing mesh parameters");

                  /* Read limits and sizes */

                  lims[0] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, -INFTY, INFTY);
                  lims[1] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, lims[0], INFTY);
                  nx = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  lims[2] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, -INFTY, INFTY);
                  lims[3] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, lims[2], INFTY);
                  ny = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  lims[4] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, -INFTY, INFTY);
                  lims[5] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, lims[4], INFTY);
                  nz = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);
                }
              else if ((type == MESH_TYPE_HEXX) ||
                       (type == MESH_TYPE_HEXY))
                {
                  /* Check number of parameters */

                  if (j > np - 9)
                    Error(loc0, "Missing mesh parameters");

                  /* Read central coordinates */

                  lims[0] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, -INFTY, INFTY);

                  lims[2] = TestParam(pname, fname, line, params[k++],
                                  PTYPE_REAL, -INFTY, INFTY);

                  /* Test if centred */

                  if ((lims[0] != 0.0) || (lims[2] != 0.0))
                    Error(loc0, "Hex mesh must be centered at the origin");

                  /* Read pitch */

                  lims[1] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, 0.0, INFTY);

                  /* Read size */

                  nx = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  ny = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  /* Read axial binning */


                  lims[4] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, -INFTY, INFTY);

                  lims[5] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, -INFTY, INFTY);

                  nz = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);
                }
              else if (type == MESH_TYPE_SPHERICAL)
                {
                  /* Check number of parameters */

                  if (k > np - 9)
                    Error(loc0, "Missing mesh parameters");

                  /* Read limits and sizes */

                  lims[0] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, 0.0, INFTY);
                  lims[1] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, lims[0], INFTY);
                  nx = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  lims[2] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, 0.0, 360.0);
                  lims[3] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, lims[2], 360.0);
                  ny = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  lims[4] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, 0.0, 180.0);
                  lims[5] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, lims[4], 180.0);
                  nz = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  /* Check unsupported types */

                  Error(loc0, "Spherical mesh type not supported yet");

                  /* Convert degrees to rad */

                  lims[2] = lims[2]*PI/180.0;
                  lims[3] = lims[3]*PI/180.0;
                  lims[4] = lims[4]*PI/180.0;
                  lims[5] = lims[5]*PI/180.0;
                }
              else if (type == MESH_TYPE_CYLINDRICAL)
                {
                  /* Check number of parameters */

                  if (k > np - 9)
                    Error(loc0, "Missing mesh parameters");

                  /* Read limits and sizes */

                  lims[0] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, 0.0, INFTY);
                  lims[1] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, lims[0], INFTY);
                  nx = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  lims[2] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, 0.0, 360.0);
                  lims[3] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, lims[2], 360.0);
                  ny = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  lims[4] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, -INFTY, INFTY);
                  lims[5] = TestParam(pname, fname, line, params[k++],
                                      PTYPE_REAL, lims[4], INFTY);
                  nz = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  if (fabs(lims[2] - lims[3]) != 360.0)
                    Error(loc0, "Angular binning must cover full sector");

                  /* Convert degrees to rad */

                  lims[2] = lims[2]*PI/180.0;
                  lims[3] = lims[3]*PI/180.0;
                }
              else
                Error(loc0, "Mesh type %ld not supported", type);

              /* Create mesh */

              loc1 = CreateMesh(type, MESH_CONTENT_PTR, -1, nx, ny, nz,
                                lims, -1);
              WDB[loc0 + RMX_PTR_MESH] = (double)loc1;

              /* Create structures */

              idx = 0;

              for (l = 0; l < nz; l++)
                for (m = 0; m < ny; m++)
                  for (n = 0; n < nx; n++)
                    {
                      /* Create structure */

                      loc2 = NewItem(loc0 + RMX_PTR_MESH_DATA,
                                     RMX_CELL_BLOCK_SIZE);

                      /* Put pointer */

                      ptr = ReadMeshPtr(loc1, n, m, l);
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      WDB[ptr] = (double)loc2;

                      /* Put index and mesh pointer */

                      WDB[loc2 + RMX_CELL_IDX] = (double)(idx);
                      WDB[loc2 + RMX_CELL_MESH_IDX] = (double)(idx++);
                      WDB[loc2 + RMX_CELL_PTR_MESH] = (double)loc1;
                    }

              /* Read submesh size */

              if (k < np)
                WDB[loc0 + RMX_CONVG_SUBMESH_NX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 100);
              else
                WDB[loc0 + RMX_CONVG_SUBMESH_NX] = 1.0;

              if (k < np)
                WDB[loc0 + RMX_CONVG_SUBMESH_NY] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 100);
              else
                WDB[loc0 + RMX_CONVG_SUBMESH_NY] = 1.0;

              if (k < np)
                WDB[loc0 + RMX_CONVG_SUBMESH_NZ] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 100);
              else
                WDB[loc0 + RMX_CONVG_SUBMESH_NZ] = 1.0;

              /* Check hexagonal */

              if ((type == MESH_TYPE_HEXX) ||
                  (type == MESH_TYPE_HEXY))
                {
                  /* Check that size matches */

                  if (RDB[loc0 + RMX_CONVG_SUBMESH_NX] !=
                      RDB[loc0 + RMX_CONVG_SUBMESH_NY])
                    Error(loc0, "Mismatch in hexagonal submesh size");

                  /* Read pitch */

                  if (k < np)
                    WDB[loc0 + RMX_CONVG_SUBMESH_P] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                ZERO, INFTY);
                  else
                    WDB[loc0 + RMX_CONVG_SUBMESH_P] = lims[1];
                }

              /* Read convergence criterion for response-matrix calculation */

              if (k < np)
                WDB[loc0 + RMX_CONV] = TestParam(pname, fname, line,
                                                 params[k++],
                                                 PTYPE_REAL, 0.0, 1E-3);

              /* Read convergence criterion for outer iterations */

              if (k < np)
                WDB[DATA_RMX_CONVG_OUT_TOL] = TestParam(pname, fname, line,
                                                        params[k++],
                                                        PTYPE_REAL, 0.0, 1E-2);

              /* Switch modes on */

              WDB[DATA_RMX_CONVG_ACC] = (double)YES;
              WDB[DATA_RMTX_CALC] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "wie"))
            {
              /***** Wieland method ******************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Avoid compiler warning */

              val = -1;

              /* Read parameter */

              if (k < np)
                val = TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                -0.999999, INFTY);
              else
                Error(-1, pname, fname, line, "Missing value");

              /* Check */

              if (val < 0.0)
                {
                  /* Put probability */

                  WDB[DATA_WIELANDT_P] = -val;

                  /* Set mode */

                  WDB[DATA_WIELANDT_MODE] = WIELANDT_MODE_FIX_P;
                }
              else
                {
                  /* Check value */

                  if ((val < 0.5) || (val > 3.0))
                    Error(-1, pname, fname, line,
                          "Weilandt k-eff should be between 0.5 and 3");

                  /* Put keff */

                  WDB[DATA_WIELANDT_KEFF] = val;

                  /* Set mode */

                  WDB[DATA_WIELANDT_MODE] = WIELANDT_MODE_FIX_K;
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "voidc"))
            {
              /***** Remove void cells ***************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_REMOVE_VOID_CELLS] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "lost"))
            {
              /***** Maximum collisions in undefined positions ***************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_MAX_UNDEF_POS] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, -1,
                            1000000000);

              /* Set mode */

              WDB[DATA_IGNORE_UNDEFINED_CELLS] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "mbtch"))
            {
              /***** Batch size for MPI transfer *****************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_OPTI_MPI_BATCH_SIZE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1000,
                            10000000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "poi"))
            {
              /***** Fission product poison xs *******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_OPTI_POISON_CALC] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Read volume */

              if (k < np)
                WDB[DATA_POISON_XS_VOL] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, INFTY);

              /* Separate treatment for Xe-135m */

              if (k < np)
                WDB[DATA_OPTI_POISON_CALC_XE135M] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "eddi"))
            {
              /***** Eddington factors ***************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_OPTI_EDDINGTON_CALC] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ba"))
            {
              /****** Nuclides treated separately as burnable absorbers ******/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check if material list is given */

              if (k == np)
                Error(-1, pname, fname, line, "Missing values");

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, np - k + 1);

              /* Set pointer */

              WDB[DATA_PTR_BA_CHAIN] = (double)loc0;

              /* Read ZAI's */

              while (k < np)
                WDB[loc0++] = TestParam(pname, fname, line, params[k++],
                                        PTYPE_INT, 10010, 840000);

              /* Put null terminator */

              WDB[loc0] = -1;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "xenon"))
            {
              /****** Xenon equilibrium calculation **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              if (k < np)
                WDB[DATA_XENON_EQUILIBRIUM_MODE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Check if material list is given */

              if (k < np)
                {
                  /* Allocate memory */

                  loc0 = ReallocMem(DATA_ARRAY, np - k + 1);

                  /* Set pointer */

                  WDB[DATA_PTR_XENON_MAT_LIST] = (double)loc0;

                  /* Read materials */

                  while (k < np)
                    WDB[loc0++] = (double)PutText(params[k++]);

                  /* Put null pointer */

                  WDB[loc0] = NULLPTR;
                }
              else if ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] == NO)
                WDB[DATA_XENON_EQUILIBRIUM_MODE] = -1.0;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "samarium"))
            {
              /****** Samarium equilibrium calculation ***********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              if (k < np)
                WDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Check if material list is given */

              if (k < np)
                {
                  /* Allocate memory */

                  loc0 = ReallocMem(DATA_ARRAY, np - k + 1);

                  /* Set pointer */

                  WDB[DATA_PTR_SAMARIUM_MAT_LIST] = (double)loc0;

                  /* Read materials */

                  while (k < np)
                    WDB[loc0++] = (double)PutText(params[k++]);

                  /* Put null pointer */

                  WDB[loc0] = NULLPTR;
                }
              else if ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] == NO)
                WDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] = -1.0;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "nphys"))
            {
              /***** Neutron reaction sampling *******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Fission */

              if (k < np)
                WDB[DATA_NPHYS_SAMPLE_FISS] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Capture */

              if (k < np)
                WDB[DATA_NPHYS_SAMPLE_CAPT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Scattering */

              if (k < np)
                WDB[DATA_NPHYS_SAMPLE_SCATT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ngamma"))
            {
              /***** Photon production from neutron reactions ****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Photon production */

              if (k < np)
                WDB[DATA_PHOTON_PRODUCTION] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 2);

              /* Mimimum weight */

              if (k < np)
                WDB[DATA_PHOTON_IMPL_WMIN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, INFTY);

              /* Maximum number of emitted photons */

              if (k < np)
                WDB[DATA_PHOTON_IMPL_NMAX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 1000000);

              /* Check and include data */

              if ((long)RDB[DATA_PHOTON_PRODUCTION] != NO)
                WDB[DATA_INCLUDE_PHOT_PROD_XS] = (double)YES;

              /* Set both transport modes on */

              WDB[DATA_NEUTRON_TRANSPORT_MODE] = (double)YES;
              WDB[DATA_PHOTON_TRANSPORT_MODE] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "depout"))
            {
              /***** Divided materials in burnup output **********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_BURN_MAT_OUTPUT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 3);

              /* Print interval */

              if (k < np)
                WDB[DATA_BURN_PRINT_STEP] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 2);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "mcvol"))
            {
              /****** Monte Carlo volume calculation *************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of points if not given as command-line parameter */

              WDB[DATA_VOLUME_MC_NMAX] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          -1E12, 1E+12);

              /* Set calculation on */

              WDB[DATA_VOLUME_MC] = YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "root"))
            {
              /****** Set root universe **************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get universe name */

              WDB[DATA_PTR_ROOT_UNIVERSE] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "csw"))
            {
              /****** Write criticality source in file ***********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get file name */

              WDB[DATA_PTR_CRIT_SRC_DET] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "gsw"))
            {
              /****** Write (n,gamma) source in file *************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get file name */

              WDB[DATA_PTR_NGAMMA_SRC_DET] = (double)PutText(params[k++]);

              /* Check option */

              if (k < np)
                WDB[DATA_NGAMMA_SRC_SIM] =
                  TestParam(pname, fname, line, params[k++],
                            PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "tle"))
            {
              /****** Track-length estimate of neutron flux ******************/

              Die(FUNCTION_NAME, "Is this doing anything?");

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Delta-tracking switch */

              if (k < np)
                {
                  n = (long)TestParam(pname, fname, line, params[k++],
                                      PTYPE_LOGICAL);

                  if (n == YES)
                    WDB[DATA_OPT_USE_DT] = (double)NO;
                  else
                    WDB[DATA_OPT_USE_DT] = (double)YES;
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "seed"))
            {
              /****** Random number seed *************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Seed */

              if (k < np)
                {
                  /* Get value */

                  parent_seed = (unsigned long)atol(params[k++]);

                  /* Check */

                  if (parent_seed < 1)
                    Error(-1, pname, fname, line, "RNG seed must be positive");

                  /* Write seed file */

                  sprintf(str, "%s.seed", GetText(DATA_PTR_INPUT_FNAME));

                  if ((fp = fopen(str, "w")) == NULL)
                    Error(-1, pname, fname, line,
                          "Unable to open seed file for writing");

                  fprintf(fp, "%lu\n", parent_seed);

                  fclose(fp);
                }

              /* Read starting batch */

              if (k < np)
                WDB[DATA_NHIST_BATCH0] =
                  (double)TestParam(pname, fname, line, params[k++], PTYPE_INT,
                                    -100000000000, 100000000000);

              /* Set replay option */

              WDB[DATA_OPTI_REPLAY] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "title"))
            {
              /****** Problem title ******************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Title */

              if (k < np)
                WDB[DATA_PTR_TITLE] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "acelib"))
            {
              /***** ACE directory file **************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Allocate memory */

              ptr = ReallocMem(DATA_ARRAY, np - k + 1);

              /* Put pointer */

              WDB[DATA_PTR_ACEDATA_FNAME_LIST] = (double)ptr;

              /* Loop over values */

              while (k < np)
                WDB[ptr++] = (double)PutText(params[k++]);

              /* Put null pointer */

              WDB[ptr] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "coverxlib"))
            {
              /***** COVERX data file ****************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Allocate memory */

              ptr = ReallocMem(DATA_ARRAY, np - k + 1);

              /* Put pointer */

              WDB[DATA_PTR_COVERXDATA_FNAME_LIST] = (double)ptr;

              /* Loop over values */

              while (k < np)
                WDB[ptr++] = (double)PutText(params[k++]);

              /* Put null pointer */

              WDB[ptr] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "declib"))
            {
              /***** Decay data file *****************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Allocate memory */

              ptr = ReallocMem(DATA_ARRAY, np - k + 1);

              /* Put pointer */

              WDB[DATA_PTR_DECDATA_FNAME_LIST] = (double)ptr;

              /* Loop over values */

              while (k < np)
                WDB[ptr++] = (double)PutText(params[k++]);

              /* Put null pointer */

              WDB[ptr] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "nfylib"))
            {
              /***** Neutron-induced fission yields **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Allocate memory */

              ptr = ReallocMem(DATA_ARRAY, np - k + 1);

              /* Put pointer */

              WDB[DATA_PTR_NFYDATA_FNAME_LIST] = (double)ptr;

              /* Loop over values */

              while (k < np)
                WDB[ptr++] = (double)PutText(params[k++]);

              /* Put null pointer */

              WDB[ptr] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "sfylib"))
            {
              /***** Spontaneous fission yields ******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Allocate memory */

              ptr = ReallocMem(DATA_ARRAY, np - k + 1);

              /* Put pointer */

              WDB[DATA_PTR_SFYDATA_FNAME_LIST] = (double)ptr;

              /* Loop over values */

              while (k < np)
                WDB[ptr++] = (double)PutText(params[k++]);

              /* Put null pointer */

              WDB[ptr] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "bralib"))
            {
              /***** Isomeric branching ratios *******************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Allocate memory */

              ptr = ReallocMem(DATA_ARRAY, np - k + 1);

              /* Put pointer */

              WDB[DATA_PTR_BRADATA_FNAME_LIST] = (double)ptr;

              /* Loop over values */

              while (k < np)
                WDB[ptr++] = (double)PutText(params[k++]);

              /* Put null pointer */

              WDB[ptr] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "egrid"))
            {
              /****** Energy grid parameters *********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Reconstruction tolerance */

              if (k < np)
                WDB[DATA_ERG_TOL] = TestParam(pname, fname, line, params[k++],
                                               PTYPE_REAL, 0.0, 0.1);

              /* Minimum energy */

              if (k < np)
                WDB[DATA_NEUTRON_EMIN] = TestParam(pname, fname, line,
                                                    params[k++], PTYPE_REAL,
                                                    0.0, 1000.0);

              /* Maximum energy */

              if (k < np)
                WDB[DATA_NEUTRON_EMAX] = TestParam(pname, fname, line,
                                                    params[k++], PTYPE_REAL,
                                                    RDB[DATA_NEUTRON_EMIN],
                                                    1000.0);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "opti"))
            {
              /****** Optimization mode **************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_OPTI_MODE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 4);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "arr"))
            {
              /****** Analog reaction rate calculation ***********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode for neutrons */

              if (k < np)
                WDB[DATA_ANA_RR_NCALC] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 2);

              /* Mode for photons */

              if (k < np)
                WDB[DATA_ANA_RR_PCALC] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 1);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "repro"))
            {
              /****** Reproducibility mode ***********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get mode */

              n = (long)TestParam(pname, fname, line, params[k++],
                                  PTYPE_INT, 0, 2);

              /* Check */

              if (n == 0)
                {
                  WDB[DATA_OPTI_OMP_REPRODUCIBILITY] = (double)NO;
                  WDB[DATA_OPTI_MPI_REPRODUCIBILITY] = (double)NO;
                }
              else if (n == 1)
                {
                  WDB[DATA_OPTI_OMP_REPRODUCIBILITY] = (double)YES;
                  WDB[DATA_OPTI_MPI_REPRODUCIBILITY] = (double)NO;
                }
              else if (n == 2)
                {
                  WDB[DATA_OPTI_OMP_REPRODUCIBILITY] = (double)YES;
                  WDB[DATA_OPTI_MPI_REPRODUCIBILITY] = (double)YES;
                }
              else
                Error(-1, pname, fname, line, "Invalid reproducibility mode");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "U235H"))
            {
              /****** U-235 fission heating value ****************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get value */

              WDB[DATA_NORM_U235_FISSE] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL, 0.0,
                          INFTY);

              /* Convert units */

              WDB[DATA_NORM_U235_FISSE] = RDB[DATA_NORM_U235_FISSE]*MEV;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ures"))
            {
              /****** Unresolved resonance probability tables ****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              if (k < np)
                WDB[DATA_USE_URES] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Reset pointer */

              loc0 = -1;

              /* Loop over parameters */

              while (k < np)
                {
                  /* Check if parameter is mode (not processed in Serpent 2) */

                  if ((strlen(params[k]) == 1) && ((*params[k] == '1') ||
                                                   (*params[k] == '2') ||
                                                   (*params[k] == '3')))

                    k++;

                  /* Check if parameter is nuclide */

                  else if (params[k][strlen(params[k]) - 1] == 'c')
                    {
                      /* Check pointer */

                      if (loc0 < 0)
                        {
                          /* Allocate memory */

                          loc0 = ReallocMem(DATA_ARRAY, np + 1);

                          /* Set pointer */

                          WDB[DATA_URES_PTR_USE_LIST] = (double)loc0;
                        }

                      /* Read nuclide */

                      WDB[loc0++] = (double)PutText(params[k++]);
                    }

                  /* Check if parameter is dilution cut-off */

                  else if ((atof(params[k]) >= 0.0) && (atof(params[k]) <= 1.0))
                    WDB[DATA_URES_DILU_CUT] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                0.0, 1.0);

                  /* Error */

                  else
                    Error(-1, pname, fname, line,
                          "Given value \"%s\" is not mode, cut-off or nuclide name", params[k]);
                }

              /* Set null pointer */

              if (loc0 > 0)
                WDB[loc0] = NULLPTR;

              /* Skip remaining parameters */

              k = np;

              /* TODO: lisää moodit, cut-offit ja muut */

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "fissh"))
            {
              /****** Nuclide-wise fission energies **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of input parameters */

              ni = np - 1;

              /* Check count */

              if (ni % 2)
                Error(-1, pname, fname, line, "Odd number of values entered");

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni + 1);

              /* Set pointer */

              WDB[DATA_NORM_PTR_FISSH] = (double)loc0;

              /* Read data */

              for (n = 0; n < ni/2; n++)
                {
                  /* Read ZAI */

                  WDB[loc0++] = TestParam(pname, fname, line, params[k++],
                                          PTYPE_INT, 90000, 1100000);

                  /* Read energy in MeV */

                  WDB[loc0++] = TestParam(pname, fname, line, params[k++],
                                          PTYPE_REAL, ZERO, INFTY);
                }

              /* Put null terminator */

              WDB[loc0] = -1.0;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "transmurea"))
            {
              /****** Transmutation reactions included in burnup mtx *********/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of input parameters */

              ni = np - 1;

              /* Check count */

              if (ni % 2)
                Error(-1, pname, fname, line, "Odd number of values entered");

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni + 1);

              /* Set pointer */

              WDB[DATA_DEP_TRANSMU_REA_LIST] = (double)loc0;

              /* Read data */

              for (n = 0; n < ni/2; n++)
                {
                  /* Read nuclide */

                  WDB[loc0++] = PutText(params[k++]);

                  /* Read mt */

                  WDB[loc0++] = PutText(params[k++]);
                }

              /* Put null terminator */

              WDB[loc0] = -1.0;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "isobra"))
            {
              /****** Fixed branching ratios to isomeric states **************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of input parameters */

              ni = np - 1;

              /* Check count */

              if (ni % 3)
                Error(-1, pname, fname, line,
                      "Incorrect number of values entered");

              /* Read data */

              for (n = 0; n < ni/3; n++)
                {
                  /* Create structure */

                  loc0 = NewItem(DATA_PTR_FIX_BRA0, FIX_BRA_BLOCK_SIZE);

                  /* Read ZAI */

                  WDB[loc0 + FIX_BRA_ZAI] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_INT, 10010, 1000000);

                  /* Read MT */

                  WDB[loc0 + FIX_BRA_MT] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_INT, 4, 199);

                  /* Read fraction */

                  WDB[loc0 + FIX_BRA_FRAC] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_REAL, 0.0, 1.0);
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "imp"))
            {
              /****** Importances ********************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Avoid compiler warning */

              type = -1;

              /* Get type */

              if (params[k][0] == 'e')
                {
                  /* Update index */

                  k++;

                  /* Energy importances, get power */

                  if (k < np)
                    WDB[DATA_ENE_IMP_POW] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                -INFTY, INFTY);

                  /* Set flag */

                  WDB[DATA_ENE_IMP_CALC] = (double)YES;
                }
              else
                {
                  /* Check type */

                  if (params[k][0] == 'c')
                    type = 1;
                  else if (params[k][0] == 'm')
                    type = 2;
                  else if (params[k][0] == 'u')
                    type = 3;
                  else
                    Error(-1, pname, fname, line, "Invalid type %s", params[k]);

                  k++;

                  /* Get number of input parameters */

                  ni = np - k;

                  /* Check count */

                  if (ni % 2)
                    Error(-1, pname, fname, line,
                          "Incorrect number of values entered");

                  /* Read data */

                  for (n = 0; n < ni/2; n++)
                    {
                      /* Create structure */

                      loc0 = NewItem(DATA_PTR_IMP0, IMP_BLOCK_SIZE);

                      /* Put name, file name and line number */

                      WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(pname);
                      WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
                      WDB[loc0 + PARAM_LINE] = (double)line;

                      /* Read cell */

                      if (type == 1)
                        WDB[loc0 + IMP_PTR_CELL] = PutText(params[k++]);
                      else if (type == 2)
                        WDB[loc0 + IMP_PTR_MAT] = PutText(params[k++]);
                      else if (type == 3)
                        WDB[loc0 + IMP_PTR_UNI] = PutText(params[k++]);
                      else
                        Error(-1, pname, fname, line, "Invalid type");

                      /* Read importance */

                      WDB[loc0 + IMP_I] =
                        TestParam(pname, fname, line, params[k++],
                                  PTYPE_REAL, 0.0, INFTY);
                    }

                  /* Set importance flag */

                  WDB[DATA_USE_GEOM_IMP] = (double)YES;
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "decomp"))
            {
              /****** Elemental decomposition ********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get option */

              if (k < np)
                WDB[DATA_ELEM_DECOMP] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Get number of input parameters */

              if ((ni = np - k) > 0)
                {
                  /* Allocate memory */

                  ptr = ReallocMem(DATA_ARRAY, ni + 1);
                  WDB[DATA_ELEM_DECOMP_PTR_LIST] = (double)ptr;

                  /* Read values */

                  while (k < np)
                    WDB[ptr++] = PutText(params[k++]);

                  /* Put null terminator */

                  WDB[ptr] = NULLPTR;
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "xscalc"))
            {
              /***** Burnup xs mode ******************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              if (k < np)
                {
                  n = (long)TestParam(pname, fname, line, params[k++],
                                      PTYPE_INT, 1, 3);
                  if (n == 1)
                    WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)NO;
                  else if (n == 2)
                    WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)YES;
                  else
                    Error(-1, pname, fname, line,
                          "Serpent 2 doesn't support xscalc = 3");
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "iter"))
            {
              /***************************************************************/

              /***** K-eff iteration *****************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check mode */

              if (k > np - 1)
                Error(-1, pname, fname, line, "Missing iteration mode");
              else if (!strcmp(params[k], "usr"))
                {
                  /* User-defined iteration (useriter.c subroutine) */

                  k++;

                  /* Set mode */

                  WDB[DATA_ITER_MODE] = (double)ITER_MODE_USER;

                  /* Read number of cycles */

                  WDB[DATA_ITER_NCYC] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT,
                              0, 10000);

                  /* Read k-eff */

                  WDB[DATA_ITER_KEFF] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              0.5, 2.0);

                  /* Number of paramers */

                  if ((ns = np - k) > 0)
                    {
                      /* Put number of input parametres */

                      WDB[DATA_ITER_USR_N_PARAM] = (double)ns;

                      /* Allocate memory for pointer array */

                      ptr = ReallocMem(DATA_ARRAY, ns);
                      WDB[DATA_ITER_USR_PTR_PARAM] = (double)ptr;

                      /* Loop over parameters */

                      for (j = 0; j < ns; j++)
                        WDB[ptr++] = PutText(params[k++]);
                    }
                }
              else if (!strcmp(params[k], "alb"))
                {
                  /* Albedo iteration */

                  k++;

                  /* Set mode */

                  WDB[DATA_ITER_MODE] = (double)ITER_MODE_ALBEDO;

                  /* Stop tracks at outer boundary */

                  WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

                  /* Read number of cycles */

                  if (k < np)
                    WDB[DATA_ITER_NCYC] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT,
                                1, 10000);

                  /* Read k-eff */

                  if (k < np)
                    WDB[DATA_ITER_KEFF] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                0.5, 2.0);

                  /* Read x-, y- and z-factors */

                  if (k < np)
                    WDB[DATA_ITER_ALB_F1] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                0.0, 1.0);

                  if (k < np)
                    WDB[DATA_ITER_ALB_F2] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                0.0, 1.0);

                  if (k < np)
                    WDB[DATA_ITER_ALB_F3] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                0.0, 1.0);
                }
              else if (!strcmp(params[k], "nuc"))
                {
                  /* Nuclide concentration iteration */

                  k++;

                  /* Set mode */

                  WDB[DATA_ITER_MODE] = (double)ITER_MODE_NUCLIDE;

                  /* Check number of parameters */

                  if (np - k < 4)
                    Error(-1, pname, fname, line, "Too few arguments for "
                          "\"iter nuc\", need at least 4, got %ld", np - k);

                  /* Read number of cycles */

                  WDB[DATA_ITER_NCYC] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT,
                              1, 10000);

                  /* Read k-eff */

                  WDB[DATA_ITER_KEFF] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              0.5, 2.0);

                  /* Create new item (dummy) */

                  loc0 = NewItem(DATA_PTR_ITER_NUC0, ITER_NUCLIDE_BLOCK_SIZE);

                  /* Put name, file name and line number */

                  WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(pname);
                  WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
                  WDB[loc0 + PARAM_LINE] = (double)line;

                  /* Get number of ZAIs */

                  n = TestParam(pname, fname, line, params[k++],
                                PTYPE_INT, 1, 100);

                  /* Allocate memory for ZAI list */

                  ptr = ReallocMem(DATA_ARRAY, n + 1);

                  /* Store ZAI list temporarily to first block */

                  WDB[loc0 + ITER_NUCLIDE_ZAI] = ptr;

                  /* Check that the number of remaining parameters is at least n */

                  if (np - k < n)
                    Error(-1, pname, fname, line,
                          "Should read %ld ZAIs for \"set nuc\", but only %ld "
                          "parameters remain.", n, np - k);

                  /* Read ZAIs */

                  for (; n > 0; n--)
                    WDB[ptr++] = TestParam(pname, fname, line, params[k++],
                                           PTYPE_INT, 10010, 999991);

                  /* Put null pointer */

                  WDB[ptr] = NULLPTR;

                  /* Check if material list is given */

                  if (k < np)
                    {
                      /* Get number of materials */

                      n = TestParam(pname, fname, line, params[k++],
                                    PTYPE_INT, 1, 100);

                      /* Allocate memory */

                      ptr = ReallocMem(DATA_ARRAY, n + 1);

                      /* Set pointer */

                      WDB[loc0 + ITER_NUCLIDE_PTR_MATERIAL_LIST] = (double)ptr;

                      /* Check that the number of remaining parameters is at least n */

                      if (np - k < n)
                        Error(-1, pname, fname, line,
                              "Should read %ld materials for \"set nuc\", but "
                              "only %ld parameters remain.", n, np - k);

                      /* Read materials */

                      for (; n > 0; n--)
                        WDB[ptr++] = (double)PutText(params[k++]);

                      /* Put null pointer */

                      WDB[ptr] = NULLPTR;
                    }

                  /* Put initial iteration value */

                  WDB[DATA_ITER_VAL] = 1.0;

                }
              else
                Error(-1, pname, fname, line, "Invalid iteration mode \"%s\"",
                      params[k]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "impl"))
            {
              /****** Implicit treatment of reaction modes *******************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Implicit capture */

              if (k < np)
                WDB[DATA_OPT_IMPL_CAPT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Implicit (n,xn) */

              if (k < np)
                WDB[DATA_OPT_IMPL_NXN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Implicit fission */

              if (k < np)
                {
                  /* Get nubar */

                  WDB[DATA_OPT_IMPL_FISS_NUBAR] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              0.0, 100.0);

                  /* Set flag */

                  if ((long)RDB[DATA_OPT_IMPL_FISS_NUBAR] > 0.0)
                    WDB[DATA_OPT_IMPL_FISS] = (double)YES;
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "dbrc"))
            {
              /****** Doppler-broadening rejection correction ****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Minimum energy */

              WDB[DATA_DBRC_EMIN] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          1E-8, 1.0);
              WDB[DATA_DBRC_EMAX] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          RDB[DATA_DBRC_EMIN], 1.0);

              /* Get number of isotopes */

              ni = np - k;

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni + 1);

              /* Set pointer */

              WDB[DATA_PTR_DBRC] = (double)loc0;

              /* Read data */

              for (n = 0; n < ni; n++)
                WDB[loc0++] = (double)PutText(params[k++]);

              /* Set null */

              WDB[loc0] = NULLPTR;

              /* Set flag */

              WDB[DATA_USE_DBRC] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "qparam_dbrc"))
            {
              /***** Q-parameter for DBRC ************************************/

              /* Qparameter (confidence level) for temperature majorants */
              /* to be used with DBRC*/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read parameter */

              if (k < np)
                {
#ifndef TRADMAJO
                  /* In case of "revisited" majorant, values around */
                  /* 1E-5 - 1E-4 are reasonal. */

                  WDB[DATA_QPARAM_DBRC] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_REAL, 1E-10, 1.0);
#else
                  /* In case of "traditional majorant", values around */
                  /* 2.0-4.0 are reasonal (the meaning of QPARAM is   */
                  /* different) */

                  WDB[DATA_QPARAM_DBRC] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_REAL, 1.0, 10.0);
#endif
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "qparam_tms"))
            {
              /***** Q-parameter for TMS *************************************/

              /* Qparameter (confidence level) for temperature majorants */
              /* to be used with TMS*/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read parameter */

              if (k < np)
                {
#ifndef TRADMAJO
                  /* In case of "revisited" majorant, values around */
                  /* 1E-5 - 1E-4 are reasonal. */

                  WDB[DATA_QPARAM_TMS] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_REAL, 1E-10, 1.0);

#else
                  /* In case of "traditional majorant", values around */
                  /* 2.5 - 4.0 are reasonal (the meaning of QPARAM is */
                  /* different) */

                  WDB[DATA_QPARAM_TMS] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_REAL, 1.0, 10.0);
#endif
                }
              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "bc"))
            {
              /****** Boundary conditions ************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check for single or multiple boundary conditions */

              if ((np == 2) || (np == 3))
                {
                  /* Condition type */

                  if (k < np)
                    {
                      if (!strcasecmp(params[k], "black"))
                        WDB[DATA_GEOM_BC0] = (double)BC_BLACK;
                      else if (!strcasecmp(params[k], "reflective"))
                        WDB[DATA_GEOM_BC0] = (double)BC_REFLECTIVE;
                      else if (!strcasecmp(params[k], "periodic"))
                        WDB[DATA_GEOM_BC0] = (double)BC_PERIODIC;
                      else
                        WDB[DATA_GEOM_BC0] =
                          TestParam(pname, fname, line, params[k], PTYPE_INT,
                                    1, 3);
                      k++;
                    }

                  /* Albedo */

                  if (k < np)
                    WDB[DATA_GEOM_ALBEDO1] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                0.0, 100.0);

                  /* Set partial conditions */

                  WDB[DATA_GEOM_BC1] = RDB[DATA_GEOM_BC0];
                  WDB[DATA_GEOM_BC2] = RDB[DATA_GEOM_BC0];
                  WDB[DATA_GEOM_BC3] = RDB[DATA_GEOM_BC0];
                }
              else if ((np == 4) || (np == 5))
                {
                  /* Condition in x-direction */

                  if (k < np)
                    {
                      if (!strcasecmp(params[k], "black"))
                        WDB[DATA_GEOM_BC1] = (double)BC_BLACK;
                      else if (!strcasecmp(params[k], "reflective"))
                        WDB[DATA_GEOM_BC1] = (double)BC_REFLECTIVE;
                      else if (!strcasecmp(params[k], "periodic"))
                        WDB[DATA_GEOM_BC1] = (double)BC_PERIODIC;
                      else
                        WDB[DATA_GEOM_BC1] =
                          TestParam(pname, fname, line, params[k], PTYPE_INT,
                                    1, 3);

                      k++;
                    }
                  else
                    Error(-1, pname, fname, line,
                          "Missing boundary condition in x-direction");

                  /* Condition in y-direction */

                  if (k < np)
                    {
                      if (!strcasecmp(params[k], "black"))
                        WDB[DATA_GEOM_BC2] = (double)BC_BLACK;
                      else if (!strcasecmp(params[k], "reflective"))
                        WDB[DATA_GEOM_BC2] = (double)BC_REFLECTIVE;
                      else if (!strcasecmp(params[k], "periodic"))
                        WDB[DATA_GEOM_BC2] = (double)BC_PERIODIC;
                      else
                        WDB[DATA_GEOM_BC2] =
                          TestParam(pname, fname, line, params[k], PTYPE_INT,
                                    1, 3);

                      k++;
                    }
                  else
                    Error(-1, pname, fname, line,
                          "Missing boundary condition in y-direction");

                  /* Condition in z-direction */

                  if (k < np)
                    {
                      if (!strcasecmp(params[k], "black"))
                        WDB[DATA_GEOM_BC3] = (double)BC_BLACK;
                      else if (!strcasecmp(params[k], "reflective"))
                        WDB[DATA_GEOM_BC3] = (double)BC_REFLECTIVE;
                      else if (!strcasecmp(params[k], "periodic"))
                        WDB[DATA_GEOM_BC3] = (double)BC_PERIODIC;
                      else
                        WDB[DATA_GEOM_BC3] =
                          TestParam(pname, fname, line, params[k], PTYPE_INT,
                                    1, 3);

                      k++;
                    }
                  else
                    Error(-1, pname, fname, line,
                          "Missing boundary condition in z-direction");

                  /* Albedo */

                  if (k < np)
                    WDB[DATA_GEOM_ALBEDO1] =
                      TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                0.0, 1.0);

                  /* Reset common condition */

                  if (((long)RDB[DATA_GEOM_BC1] != BC_BLACK) ||
                      ((long)RDB[DATA_GEOM_BC2] != BC_BLACK) ||
                      ((long)RDB[DATA_GEOM_BC3] != BC_BLACK))
                    WDB[DATA_GEOM_BC0] = -1.0;
                }

              /* Set other two albedos */

              WDB[DATA_GEOM_ALBEDO2] = RDB[DATA_GEOM_ALBEDO1];
              WDB[DATA_GEOM_ALBEDO3] = RDB[DATA_GEOM_ALBEDO1];

              /***************************************************************/
            }
          else if ((!strcasecmp(params[j], "dff")) ||
                   (!strcasecmp(params[j], "dfr")))
            Die(FUNCTION_NAME, "Poistit naa: dff, dfr");
          else if (!strcasecmp(params[j], "adf"))
            {
              /***** Discontinuity factors ***********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Create new item */

              loc0 = NewItem(DATA_PTR_ADF0, ADF_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(pname);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Get universe */

              if (k < np)
                WDB[loc0 + ADF_PTR_GCU] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing universe name");

              /* Get surface */

              if (k < np)
                WDB[loc0 + ADF_PTR_SURF] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing surface name");

              /* Reset enforce option */

              WDB[loc0 + ADF_ENFORCE_REF] = (double)NO;

              /* Read symmetry */

              if (k < np)
                {
                  /* Check type */

                  if ((!strcasecmp(params[k], "no")) ||
                      (!strcasecmp(params[k], "none")))
                    {
                      n = 0;
                      k++;
                    }
                  else if (!strcasecmp(params[k], "full"))
                    {
                      n = 1;
                      k++;
                    }
                  else if ((!strcasecmp(params[k], "s1")) ||
                      (!strcasecmp(params[k], "ns")) ||
                      (!strcasecmp(params[k], "sn")))
                    {
                      n = 2;
                      k++;
                    }
                  else if ((!strcasecmp(params[k], "c1")) ||
                      (!strcasecmp(params[k], "nesw")) ||
                      (!strcasecmp(params[k], "swne")))
                    {
                      n = 3;
                      k++;
                    }
                  else if ((!strcasecmp(params[k], "s2")) ||
                           (!strcasecmp(params[k], "ew")) ||
                           (!strcasecmp(params[k], "we")))
                    {
                      n = 4;
                      k++;
                    }
                  else if ((!strcasecmp(params[k], "c2")) ||
                           (!strcasecmp(params[k], "nwse")) ||
                           (!strcasecmp(params[k], "senw")))
                    {
                      n = 5;
                      k++;
                    }
                  else if ((!strcasecmp(params[k], "s3")) ||
                           (!strcasecmp(params[k], "nsew")) ||
                           (!strcasecmp(params[k], "ewns")))
                    {
                      n = 6;
                      k++;
                    }
                  else if (!strcasecmp(params[k], "c3"))
                    {
                      n = 7;
                      k++;
                    }
                  else
                    n = (long)TestParam(pname, fname, line, params[k++],
                                        PTYPE_INT, 0, 7);

                  /* Put symmetry */

                  WDB[loc0 + ADF_SYM] = (double)n;
                }
              else
                Error(-1, pname, fname, line, "Missing symmetry");

              /* Read enforce flag */

              if (k < np)
                WDB[loc0 + ADF_ENFORCE_REF] =
                  TestParam(pname, fname, line, params[k++],
                            PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "alb"))
            {
              /***** Albedos (group constant) ********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Create new item */

              loc0 = NewItem(DATA_PTR_ALB0, ALB_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(pname);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Get universe */

              if (k < np)
                WDB[loc0 + ALB_PTR_GCU] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing universe name");

              /* Get surface */

              if (k < np)
                WDB[loc0 + ALB_PTR_SURF] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing surface name");

              /* Get direction */

              if (k < np)
                {
                  /* Get parameter */

                  n = (long)TestParam(pname, fname, line, params[k++],
                                      PTYPE_INT, -1, 1);

                  /* Check */

                  if (n == -1)
                    {
                      /* Inward direction */

                      WDB[loc0 + ALB_DIR] = 0;
                    }
                  else if (n == 0)
                    {
                      /* Not allowed */

                      Error(-1, pname, fname, line,
                            "Direction must be -1 (= inward) or 1 (= outward)");
                    }
                  else
                    {
                      /* Outward direction */

                      WDB[loc0 + ALB_DIR] = 1;
                    }
                }
              else
                Error(-1, pname, fname, line,
                      "Missing direction (-1 = inward, 1 = outward)");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "mdep"))
            {
              /****** Cross sections for micro-depletion *********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Create new item */

              loc0 = NewItem(DATA_PTR_MDEP0, MDEP_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(pname);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Get universe */

              if (k < np)
                WDB[loc0 + MDEP_PTR_GCU] = (double)PutText(params[k++]);
              else
                Error(loc0, "Missing universe name");

              /* Avoid compiler warning */

              ni = -1;

              /* Read volume */

              if (k < np)
                WDB[loc0 + MDEP_VOLUME] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, INFTY);
              else
                Error(loc0, "Missing number of materials");


              /* Read number of materials */

              if (k < np)
                ni = TestParam(pname, fname, line, params[k++], PTYPE_INT,
                               0, 100000);
              else
                Error(loc0, "Missing number of materials");

              /* Check */

              if (ni == 0)
                WDB[loc0 + MDEP_PTR_MAT] = NULLPTR;
              else if (ni > np - k)
                Error(loc0, "Missing material entries");
              else
                {
                  /* Allocate memory for data */

                  ptr = ReallocMem(DATA_ARRAY, ni + 1);
                  WDB[loc0 + MDEP_PTR_MAT] = (double)ptr;

                  /* Read data */

                  for (n = 0; n < ni; n++)
                    WDB[ptr++] = (double)PutText(params[k++]);

                  /* Put null terminator */

                  WDB[ptr] = -1.0;
                }

              /* Get number of remaining entries */

              ni = np - k;

              /* Check */

              if (ni == 0)
                Error(loc0, "Missing reactions entries");
              if (ni % 2)
                Error(loc0, "Invalid number of reaction entries");

              /* Allocate memory for data */

              ptr = ReallocMem(DATA_ARRAY, ni/2 + 1);
              WDB[loc0 + MDEP_PTR_KEY] = (double)ptr;

              /* Put number of values */

              WDB[loc0 + MDEP_N_REA] = (double)(ni/2);

              /* Read data */

              for (n = 0; n < ni; n = n + 2)
                {
                  /* Get ZAI */

                  ZAI = TestParam(pname, fname, line, params[k++],
                                  PTYPE_INT, 10010, 992550);

                  /* Check last character of MT */

                  c = params[k][strlen(params[k]) - 1];

                  /* Put flag */

                  if ((c == 'g') || (c == 'G'))
                    {
                      m = 1;
                      params[k][strlen(params[k]) - 1] = '\0';
                    }
                  else if ((c == 'i') || (c == 'I') ||
                           (c == 'm') || (c == 'M'))
                    {
                      m = 2;
                      params[k][strlen(params[k]) - 1] = '\0';
                    }
                  else
                    m = 0;

                  /* Get mt */

                  MT = TestParam(pname, fname, line, params[k++],
                                 PTYPE_INT, 2, 999);

                  /* Check fission modes and convert */

                  if ((MT > 180) && (MT < 190))
                    {
                      m = MT - 180 + 2;
                      MT = 18;
                    }
                  else if ((MT > 190) && (MT < 200))
                    {
                      m = MT - 190 + 2;
                      MT = 19;
                    }

                  /* Store search key */

                  WDB[ptr++] = (double)(10000*ZAI + 10*MT + m);
                }

              /* Put null terminator */

              WDB[ptr] = -1.0;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "icm"))
            {
              /***** Data for interface current method ***********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if (k > np - 1)
                Error(-1, pname, fname, line, "Missing parameters");

              while (k < np)
                {
                  /* Check */

                  if (!strcasecmp(params[k], "erg"))
                    {
                      /* Read ICM energy grid */

                      if (++k < np)
                        WDB[DATA_ICM_PTR_ENE0] = (double)PutText(params[k++]);
                      else
                        Error(-1, pname, fname, line,
                              "Missing ICM energy grid");

                      /* Read reconstruction energy grid */

                      if (k < np)
                        WDB[DATA_ICM_PTR_ENE1] = (double)PutText(params[k++]);
                      else
                        Error(-1, pname, fname, line,
                              "Missing reconstruction energy grid");
                    }
                  else if (!strcasecmp(params[k], "sub"))
                    {
                      /* Read segmentation */

                      if (++k < np)
                        WDB[DATA_ICM_NSUB] =
                          TestParam(pname, fname, line, params[k++], PTYPE_INT,
                                    1, 50);
                      else
                        Error(-1, pname, fname, line, "Missing segmentation");

                      /* Allocate memory for data */

                      ns = (long)RDB[DATA_ICM_NSUB];
                      ptr = ReallocMem(DATA_ARRAY, ns + 1);
                      WDB[DATA_ICM_PTR_SUB] = (double)ptr;

                      /* First bin */

                      WDB[ptr++] = 0.0;

                      /* Read values */

                      for (n = 0; n < ns - 1; n++)
                        WDB[ptr++] = TestParam(pname, fname, line, params[k++],
                                               PTYPE_REAL, 0.0, 1.0);

                      /* Last bin */

                      WDB[ptr++] = 1.0;
                    }
                  else if (!strcasecmp(params[k], "mu0"))
                    {
                      /* Read perpendicular bins */

                      if (++k < np)
                        WDB[DATA_ICM_NMU0] =
                          TestParam(pname, fname, line, params[k++], PTYPE_INT,
                                    1, 50);
                      else
                        Error(-1, pname, fname, line, "Missing mu bins");

                      /* Allocate memory for data */

                      ns = (long)RDB[DATA_ICM_NMU0];
                      ptr = ReallocMem(DATA_ARRAY, ns + 1);
                      WDB[DATA_ICM_PTR_MU0] = (double)ptr;

                      /* First bin */

                      WDB[ptr++] = 0.0;

                      /* Read values */

                      for (n = 0; n < ns - 1; n++)
                        WDB[ptr++] = TestParam(pname, fname, line, params[k++],
                                               PTYPE_REAL, 0.0, 1.0);

                      /* Last bin */

                      WDB[ptr++] = 1.0;
                    }
                  else if (!strcasecmp(params[k], "mu1"))
                    {
                      /* Read parallel bins */

                      if (++k < np)
                        WDB[DATA_ICM_NMU1] =
                          TestParam(pname, fname, line, params[k++], PTYPE_INT,
                                    1, 50);
                      else
                        Error(-1, pname, fname, line, "Missing mu bins");

                      /* Allocate memory for data */

                      ns = (long)RDB[DATA_ICM_NMU1];
                      ptr = ReallocMem(DATA_ARRAY, ns + 1);
                      WDB[DATA_ICM_PTR_MU1] = (double)ptr;

                      /* First bin */

                      WDB[ptr++] = -1.0;

                      /* Read values */

                      for (n = 0; n < ns - 1; n++)
                        WDB[ptr++] = TestParam(pname, fname, line, params[k++],
                                               PTYPE_REAL, -1.0, 1.0);

                      /* Last bin */

                      WDB[ptr++] = 1.0;
                    }
                  else if (!strcasecmp(params[k], "mu2"))
                    {
                      /* Read axial bins */

                      if (++k < np)
                        WDB[DATA_ICM_NMU2] =
                          TestParam(pname, fname, line, params[k++], PTYPE_INT,
                                    1, 50);
                      else
                        Error(-1, pname, fname, line, "Missing mu bins");

                      /* Allocate memory for data */

                      ns = (long)RDB[DATA_ICM_NMU2];
                      ptr = ReallocMem(DATA_ARRAY, ns + 1);
                      WDB[DATA_ICM_PTR_MU2] = (double)ptr;

                      /* First bin */

                      WDB[ptr++] = -1.0;

                      /* Read values */

                      for (n = 0; n < ns - 1; n++)
                        WDB[ptr++] = TestParam(pname, fname, line, params[k++],
                                               PTYPE_REAL, -1.0, 1.0);

                      /* Last bin */

                      WDB[ptr++] = 1.0;
                    }
                  else
                    {
                      /* Create new item */

                      loc0 = NewItem(DATA_PTR_ICM0, ICM_BLOCK_SIZE);

                      /* Put name, file name and line number */

                      WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(pname);
                      WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
                      WDB[loc0 + PARAM_LINE] = (double)line;

                      /* Read ID */

                      WDB[loc0 + ICM_PTR_ID] = (double)PutText(params[k++]);

                      /* Read surface */

                      if (k < np)
                        WDB[loc0 + ICM_PTR_SURF] =
                          (double)PutText(params[k++]);
                      else
                        Error(-1, pname, fname, line, "Missing surface name");

                      /* Read lattice */

                      if (k < np)
                        WDB[loc0 + ICM_PTR_LAT] = (double)PutText(params[k++]);

                      /* Switch calculation on */

                      WDB[DATA_ICM_CALC] = (double)YES;
                    }
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ppw"))
            {
              /***** Pin power distributions *********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Create new item */

              loc0 = NewItem(DATA_PTR_PPW0, PPW_BLOCK_SIZE);

              /* Get universe */

              if (k < np)
                WDB[loc0 + PPW_PTR_GCU] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing universe name");

              /* Get lattice name */

              if (k < np)
                WDB[loc0 + PPW_PTR_LAT] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing lattice name");

              /* TODO: Lisää symmetria */

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "trc"))
            {
              /***** Material-wise transport correction **********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Create new item */

              loc0 = NewItem(DATA_PTR_TRC0, TRANSP_CORR_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(pname);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Get material name */

              if (k < np)
                WDB[loc0 + TRANSP_CORR_PTR_MAT] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing material name");

              /* Get file name */

              if (k < np)
                WDB[loc0 + TRANSP_CORR_FNAME] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing file name");

              /* Get minimum energy */

              if (k < np)
                WDB[loc0 + TRANSP_CORR_EMIN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, INFTY);
              else
                Error(-1, pname, fname, line, "Missing minimum energy");

              /* Read nuclides */

              if ((ni = np - k) > 0)
                {
                  /* Allocate memory for ZAI */

                  ptr = ReallocMem(DATA_ARRAY, ni + 1);
                  WDB[loc0 + TRANSP_CORR_PTR_ISO] = (double)ptr;

                  /* Read data */

                  for (n = 0; n < ni; n++)
                    WDB[ptr++] = (double)TestParam(pname, fname, line,
                                                   params[k++], PTYPE_INT,
                                                   10010, 992550);
                  /* Put null terminator */

                  WDB[ptr] = -1.0;
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "csm"))
            {
              /***** Super-imposed cell search mesh **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Create new item */

              loc0 = NewItem(DATA_PTR_CSM0, CELL_MESH_BLOCK_SIZE);

              /* Get universe */

              if (k < np)
                WDB[loc0 + CELL_MESH_PTR_UNI] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing universe name");

              /* Get mesh type */

              if (k < np)
                WDB[loc0 + CELL_MESH_TYPE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 5);
              else
                Error(-1, pname, fname, line, "Missing mesh type");

              /* Check number of remaining parameters */

              if (np - k != 9)
                Error(-1, pname, fname, line,
                      "Invalid number of parameters given");

              /* Read min, max and size (times three) */

              WDB[loc0 + CELL_MESH_MIN0] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          -INFTY, INFTY);

              WDB[loc0 + CELL_MESH_MAX0] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          -INFTY, INFTY);

              WDB[loc0 + CELL_MESH_N0] =
                TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 1000);

              WDB[loc0 + CELL_MESH_MIN1] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          -INFTY, INFTY);

              WDB[loc0 + CELL_MESH_MAX1] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          -INFTY, INFTY);

              WDB[loc0 + CELL_MESH_N1] =
                TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 1000);

              WDB[loc0 + CELL_MESH_MIN2] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          -INFTY, INFTY);

              WDB[loc0 + CELL_MESH_MAX2] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                          -INFTY, INFTY);

              WDB[loc0 + CELL_MESH_N2] =
                TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 1000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "mora"))
            {
              /***** Cross sections for MORA calculation *********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Create new item */

              loc0 = NewItem(DATA_PTR_MORA0, MORA_BLOCK_SIZE);

              /* Get file name */

              if (k < np)
                WDB[loc0 + MORA_PTR_FNAME] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing file name");

              /* Get universe */

              if (k < np)
                WDB[loc0 + MORA_PTR_UNIV] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing universe name");

              /* Get group structure */

              if (k < np)
                WDB[loc0 + MORA_PTR_EG] = (double)PutText(params[k++]);
              else
                Error(-1, pname, fname, line, "Missing group structure");

              /* Get number of cosine bins */

              if (k < np)
                WDB[loc0 + MORA_N_COS] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 128);
              else
                Error(-1, pname, fname, line, "Missing number of cosine bins");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "usym"))
            {
              /***** Universe symmetry ***************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Create new item */

              loc0 = NewItem(DATA_PTR_SYM0, SYMMETRY_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(params[j]);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Check number of parameters */

              if ((np == 8) || (np == 9))
                {
                  /* Serpent 2 style, read Universe name */

                  WDB[loc0 + SYMMETRY_PTR_UNI] = (double)PutText(params[k++]);

                  /* Read axis */

                  WDB[loc0 + SYMMETRY_AXIS] = (double)PutText(params[k++]);

                  /* Read mode */

                  WDB[loc0 + SYMMETRY_BC] = (double)PutText(params[k++]);

                  /* Read symmetry origin */

                  WDB[loc0 + SYMMETRY_X0] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_REAL, -INFTY, INFTY);

                  WDB[loc0 + SYMMETRY_Y0] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_REAL, -INFTY, INFTY);

                  /* Read angles */

                  WDB[loc0 + SYMMETRY_THETA0] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL, 0.0,
                              360.0);

                  WDB[loc0 + SYMMETRY_ROT] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL, 0.0,
                              180.0);

                  /* Read coordinate transformation option */

                  if (k < np)
                    WDB[loc0 + SYMMETRY_COORD_TRANS] =
                      TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);
                }
              else if ((np == 5) || (np == 3))
                {
                  /* Serpent 1 style, read Universe name */

                  WDB[loc0 + SYMMETRY_PTR_UNI] =(double)PutText(params[k++]);

                  /* Read symmetry (allow quqdrant only for now) */

                  WDB[loc0 + SYMMETRY_SYM] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_INT, 4, 4);

                  /* Read symmetry origin */

                  if (k < np)
                    WDB[loc0 + SYMMETRY_X0] =
                      TestParam(pname, fname, line, params[k++],
                                PTYPE_REAL, -INFTY, INFTY);

                  if (k < np)
                    WDB[loc0 + SYMMETRY_Y0] =
                      TestParam(pname, fname, line, params[k++],
                                PTYPE_REAL, -INFTY, INFTY);
                }
              else
                Error(-1, pname, fname, line,
                      "Invalid number of input values");

              /* Check (this is due to bug in code) */

              if ((RDB[loc0 + SYMMETRY_X0] != 0.0) || (RDB[loc0 + SYMMETRY_Y0]))
                Error(-1, pname, fname, line,
                      "Symmetries don't work unless centeret at (0,0)");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "fmtx"))
            {
              /***** Fission matrix ******************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Create new item */

              loc0 = ReallocMem(DATA_ARRAY, FMTX_BLOCK_SIZE);
              WDB[DATA_PTR_FMTX] = (double)loc0;

              /* Reset pointers and level */

              WDB[loc0 + FMTX_PTR_MAT] = NULLPTR;
              WDB[loc0 + FMTX_PTR_UNI] = NULLPTR;
              WDB[loc0 + FMTX_PTR_MESH] = NULLPTR;
              WDB[loc0 + FMTX_LVL] = -1;

              /* Read type */

              WDB[DATA_FMTX_TYPE] =
                TestParam(pname, fname, line, params[k++],
                          PTYPE_INT, 1, 4);

              /* Check type */

              if ((long)RDB[DATA_FMTX_TYPE] == FISSION_MATRIX_TYPE_MAT)
                {
                  /* Read materials */

                  if (k < np)
                    {
                      /* Allocate memory */

                      loc1 = ReallocMem(DATA_ARRAY, np - k + 1);
                      WDB[loc0 + FMTX_PTR_MAT] = (double)loc1;

                      /* Read materials */

                      while(k < np)
                        WDB[loc1++] = (double)PutText(params[k++]);

                      /* Put null terminator */

                      WDB[loc1] = NULLPTR;
                    }
                }
              else if ((long)RDB[DATA_FMTX_TYPE] == FISSION_MATRIX_TYPE_UNI)
                {
                  /* Read universes */

                  if (k < np)
                    {
                      /* Allocate memory */

                      loc1 = ReallocMem(DATA_ARRAY, np - k + 1);
                      WDB[loc0 + FMTX_PTR_UNI] = (double)loc1;

                      /* Read materials */

                      while(k < np)
                        WDB[loc1++] = (double)PutText(params[k++]);

                      /* Put null terminator */

                      WDB[loc1] = NULLPTR;
                    }
                }
              else if ((long)RDB[DATA_FMTX_TYPE] == FISSION_MATRIX_TYPE_LVL)
                {
                  /* Read level */

                  WDB[loc0 + FMTX_LVL] =
                    TestParam(pname, fname, line, params[k++],
                              PTYPE_INT, 0, MAX_GEOMETRY_LEVELS);
                }
              else if ((long)RDB[DATA_FMTX_TYPE] == FISSION_MATRIX_TYPE_XYZ)
                {
                  /* Check parameters */

                  if (np - k != 9)
                    Error(-1, pname, fname, line,
                          "Invalid number of parameters");

                  /* Read values */

                  xmin = TestParam(pname, fname, line, params[k++],
                                   PTYPE_REAL, -INFTY, INFTY);
                  xmax = TestParam(pname, fname, line, params[k++],
                                   PTYPE_REAL, xmin, INFTY);
                  nx = (long)TestParam(pname, fname, line, params[k++],
                                      PTYPE_INT, 1, 1000000);

                  ymin = TestParam(pname, fname, line, params[k++],
                                   PTYPE_REAL, -INFTY, INFTY);
                  ymax = TestParam(pname, fname, line, params[k++],
                                   PTYPE_REAL, ymin, INFTY);
                  ny = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 1000000);

                  zmin = TestParam(pname, fname, line, params[k++],
                                   PTYPE_REAL, -INFTY, INFTY);
                  zmax = TestParam(pname, fname, line, params[k++],
                                   PTYPE_REAL, zmin, INFTY);
                  nz = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 1000000);

                  /* Put mesh variables */

                  lims[0] = xmin;
                  lims[1] = xmax;
                  lims[2] = ymin;
                  lims[3] = ymax;
                  lims[4] = zmin;
                  lims[5] = zmax;

                  /* Create mesh */

                  ptr = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_RES, -1,
                                   nx, ny, nz, lims, -1);

                  /* Put pointer */

                  WDB[loc0 + FMTX_PTR_MESH] = (double)ptr;
                }
              else
                Error(-1, pname, fname, line, "Invalid fission matrix type");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "xsplot"))
            {
              /****** Cross section data plot ********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Resolution */

              if (k < np)
                WDB[DATA_XSPLOT_NE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 10,
                            10000);

              /* boundaries */

              if (k < np)
                WDB[DATA_XSPLOT_EMIN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, INFTY);

              if (k < np)
                WDB[DATA_XSPLOT_EMAX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            RDB[DATA_XSPLOT_EMIN], INFTY);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "cfe") ||
                   !strcasecmp(params[j], "minxs"))
            {
              /****** Minimum macroscopic xs for CFE ************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mean distance for neutrons */

              if (k < np)
                WDB[DATA_CFE_N_MIN_L] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            -1.0, INFTY);

              /* Mean time for neutrons */

              if (k < np)
                WDB[DATA_CFE_N_MIN_T] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            -1.0, INFTY);

              /* Mean distance for photons */

              if (k < np)
                WDB[DATA_CFE_G_MIN_L] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            -1.0, INFTY);

              /* Mean time for photons */

              if (k < np)
                WDB[DATA_CFE_G_MIN_T] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            -1.0, INFTY);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "inftrk"))
            {
              /****** Handling of infinite tracking loops *******************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Maximum number of loops for neutrons */

              if (k < np)
                WDB[DATA_NEUTRON_MAX_TRACK_LOOP] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 10000000);

              /* Error flag for neutrons */

              if (k < np)
                WDB[DATA_NEUTRON_MAX_TRACK_LOOP_ERR] =
                  TestParam(pname, fname, line, params[k++],
                            PTYPE_LOGICAL);

              /* Maximum number of loops for photons */

              if (k < np)
                WDB[DATA_PHOTON_MAX_TRACK_LOOP] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 10000000);

              /* Error flag for photons */

              if (k < np)
                WDB[DATA_PHOTON_MAX_TRACK_LOOP_ERR] =
                  TestParam(pname, fname, line, params[k++],
                            PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "nps"))
            {
              /****** Run parameters (source mode) **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Number of neutron histories */

              if (k < np)
                {
                  /* Use atol, because TestParam cannot handle long ints */

                  n = atol(params[k++]);

                  /* Check value */

                  if (n < 1)
                    Error(-1, pname, fname, line,
                          "Minimum number of source particles is 1");
                  else
                    WDB[DATA_SRC_POP] = (double)n;
                }

              /* Number of batches */

              if (k < np)
                WDB[DATA_SRC_BATCHES] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 100000000000);
              else
                {
                  /* Set default to 200 */

                  WDB[DATA_SRC_BATCHES] = 200;
                }

              /* Divide number of histories by batches */

              WDB[DATA_SRC_POP] =
                (double)((long)(RDB[DATA_SRC_POP]/RDB[DATA_SRC_BATCHES]));

              /* Check that batch size is reasonable */

              if ((long)RDB[DATA_SRC_POP] < 1)
                Error(-1, pname, fname, line,
                      "Insufficient source size (integer value expected)");

              /* Read pointer to time binning */

              if (k < np)
                WDB[DATA_DYN_PTR_TIME_BINS] = (double)PutText(params[k++]);

              /* Set source mode */

              WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_SRC;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "pop"))
            {
              /****** Run parameters (criticality mode) *********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Neutron population */

              if (k < np)
                WDB[DATA_CRIT_POP] =
                  TestParam(pname, fname, line,  params[k++], PTYPE_INT,
                            10, 100000000000);

              /* Number of active cycles */

              if (k < np)
                WDB[DATA_CRIT_CYCLES] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 100000000000);

              /* Number of inactive cycles */

              if (k < np)
                WDB[DATA_CRIT_SKIP] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1,  100000000000);

              /* Initial guess for k-eff */

              if (k < np)
                WDB[DATA_CYCLE_KEFF] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.5, 3.0);

              /* Batching interval */

              if (k < np)
                {
                  /* Get number of batches */

                  n = (long)TestParam(pname, fname, line, params[k++],
                                      PTYPE_INT, 1,
                                      (long)RDB[DATA_CRIT_CYCLES]);

                  /* Check */

                  if ((long)RDB[DATA_CRIT_CYCLES] % n)
                    Error(-1, pname, fname, line,
                          "Number of cycles must be divisible with the number of batches");

                  /* Put value */

                  WDB[DATA_BATCH_INTERVAL] = (double)n;
                }

              /* Number of parallel eigenvalue calculations */

              if (k < np)
                WDB[DATA_N_POP_EIG] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 10000);

              /* Set criticality mode */

              WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_CRIT;

              /* Set neutron transport flag */

              WDB[DATA_NEUTRON_TRANSPORT_MODE] = (double)YES;

             /* Store original population size for corrector iteration */

              WDB[DATA_BURN_CI_ORIG_NBATCH] = RDB[DATA_CRIT_POP];
              WDB[DATA_BURN_CI_ORIG_CYCLES] = RDB[DATA_CRIT_CYCLES];
              WDB[DATA_BURN_CI_ORIG_SKIP]   = RDB[DATA_CRIT_SKIP];

              /* Check number of eigenvalue simulations */

              if (((long)RDB[DATA_CRIT_POP]) % ((long)RDB[DATA_N_POP_EIG]))
                Error(-1, pname, fname, line,
                      "Population size must be divisible with number of eigenvalue simulations");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "gpop"))
            {
              /****** Growing population *************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get maximum number of histories */

              if (k < np)
                {
                  /* Use atol, because TestParam cannot handle long ints */

                  n = atol(params[k++]);

                  /* Check value */

                  if (n < 1)
                    Error(-1, pname, fname, line,
                          "Minimum number of source particles is 1");
                  else
                    WDB[DATA_GROW_POP_MAX_HIS] = (double)n;
                }
              else
                Error(-1, pname, fname, line,
                      "Missing maximum number of histories");

              /* Get maximum population */

              if (k < np)
                WDB[DATA_GROW_POP_MAX_POP] =
                  (double)TestParam(pname, fname, line, params[k++], PTYPE_INT,
                                    1, 1000000000);

              /* Set parameter c */

              if (k < np)
                WDB[DATA_GROW_POP_C] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, 1e12);

              /* Set gpop mode on */

              WDB[DATA_GROW_POP_SIM] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "nfg"))
            {
              /****** Few-group parameters ***********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Serpent 1 style default 1-, 2- or 4-group structure */

              if ((np == 2) && (!strcasecmp(params[k], "1") ||
                                !strcasecmp(params[k], "2") ||
                                !strcasecmp(params[k], "4")))
                {
                  /* Reset pointer to predefined */

                  WDB[DATA_ERG_FG_PTR_PREDEF] = NULLPTR;

                  /* Get number of groups */

                  ne = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 1, 10000);

                  /* Set value */

                  WDB[DATA_ERG_FG_NG] = (double)ne;

                  /* Allocate memory */

                  loc0 = ReallocMem(DATA_ARRAY, ne + 1);

                  /* Set pointer */

                  WDB[DATA_ERG_FG_PTR_GRID] = (double)loc0;

                  /* Set boundaries */

                  if (ne == 1)
                    {
                      WDB[loc0++] = 0.0;
                      WDB[loc0++] = INFTY;
                    }
                  else if (ne == 2)
                    {
                      WDB[loc0++] = 0.0;
                      WDB[loc0++] = 0.625E-6;
                      WDB[loc0++] = INFTY;
                    }
                  else if (ne == 4)
                    {
                      WDB[loc0++] = 0.0;
                      WDB[loc0++] = 0.625E-6;
                      WDB[loc0++] = 5.5E-3;
                      WDB[loc0++] = 0.821;
                      WDB[loc0++] = INFTY;
                    }
                  else
                    Die(FUNCTION_NAME, "Error in boundaries");
                }

              /* Single entry is interpreted as a pre-defined group structure */

              else if (np == 2)
                WDB[DATA_ERG_FG_PTR_PREDEF] = (double)PutText(params[k++]);

              /* Serpent 1-style structure with given boundaries */

              else if (k < np)
                {
                  /* Reset pointer to predefined */

                  WDB[DATA_ERG_FG_PTR_PREDEF] = NULLPTR;

                  /* Get number of groups */

                  ne = (long)TestParam(pname, fname, line, params[k++],
                                       PTYPE_INT, 2, 10000);

                  /* Set value */

                  WDB[DATA_ERG_FG_NG] = (double)ne;

                  /* Check if structure is given */

                  if (np > 2)
                    {
                      /* Check number of given parameters */

                      if (ne != np - 1)
                        Error(-1, pname, fname, line,
                              "Invalid number of values in few-group structure");

                      /* Allocate memory */

                      loc0 = ReallocMem(DATA_ARRAY, ne + 1);

                      /* Set pointer */

                      WDB[DATA_ERG_FG_PTR_GRID] = (double)loc0;

                      /* Set first value */

                      WDB[loc0++] = 0.0;

                      /* Read boundaries */

                      for (n = 0; n < ne - 1 ; n++)
                        WDB[loc0++] =
                          TestParam(pname, fname, line, params[k++],
                                    PTYPE_REAL, 0.0, 1000.0);

                      /* Set last value */

                      WDB[loc0++] = INFTY;
                    }
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "gcu"))
            {
              /***** Universes for group constant generation *****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check exclude mode */

              if (atoi(params[k]) == -1)
                {
                  /* Put negative pointer */

                  WDB[DATA_PTR_GCU0] = -1.0;

                  /* Skip remaining parameters */

                  k = np;
                }
              else
                {
                  /* Get number of universes */

                  ni = np - k;

                  /* Loop over universes */

                  for (n = 0; n < ni; n++)
                    {
                      /* Allocate memory */

                      loc0 = NewItem(DATA_PTR_GCU0, GCU_BLOCK_SIZE);

                      /* Get universe */

                      WDB[loc0 + GCU_PTR_UNIV] = (double)PutText(params[k++]);
                    }
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "cpd"))
            {
              /****** Core power distribution ********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Depth */

              if (k < np)
                WDB[DATA_CORE_PDE_DEPTH] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 3);

              /* Reset Z-bins */

              WDB[DATA_CORE_PDE_NZ] = 1.0;
              WDB[DATA_CORE_PDE_ZMIN] = -INFTY;
              WDB[DATA_CORE_PDE_ZMAX] = INFTY;

              /* Number of Z-bins */

              if (k < np)
                WDB[DATA_CORE_PDE_NZ] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                            1000);

              /* Minimum z */

              if (k < np)
                WDB[DATA_CORE_PDE_ZMIN] =
                  TestParam(pname, fname, line, params[k++],
                                    PTYPE_REAL, -INFTY, INFTY);

              /* Maximum z */

              if (k < np)
                WDB[DATA_CORE_PDE_ZMAX] =
                  TestParam(pname, fname, line, params[k++],
                            PTYPE_REAL, RDB[DATA_CORE_PDE_ZMIN], INFTY);

              /* Read levels */

              if (k < np)
                WDB[DATA_CORE_PDE_L1] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, -1,
                            1000);

              if (k < np)
                WDB[DATA_CORE_PDE_L2] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, -1,
                            1000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ufs"))
            {
              /****** Uniform fission source method **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Reset size and dimensions */

              WDB[DATA_UFS_NX] = 1.0;
              WDB[DATA_UFS_NY] = 1.0;
              WDB[DATA_UFS_NZ] = 1.0;

              WDB[DATA_UFS_XMIN] = -INFTY;
              WDB[DATA_UFS_XMAX] = INFTY;
              WDB[DATA_UFS_YMIN] = -INFTY;
              WDB[DATA_UFS_YMAX] = INFTY;
              WDB[DATA_UFS_ZMIN] = -INFTY;
              WDB[DATA_UFS_ZMAX] = INFTY;

              /* Get mode */

              if (k < np)
                WDB[DATA_UFS_MODE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 4);

              /* Get order */

              if (k < np)
                WDB[DATA_UFS_ORDER] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL, -5.0,
                            5.0);

              /* Check number of parameters */

              if (np == 6)
                {
                  /* Get mesh size */

                  WDB[DATA_UFS_NX] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                                1000);
                  WDB[DATA_UFS_NY] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                              1000);
                  WDB[DATA_UFS_NZ] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                              1000);
                }
              else if (np == 7)
                {
                  /* Get lattice name */

                  WDB[DATA_UFS_PTR_LAT] = (double)PutText(params[k++]);

                  /* Get axial size and dimensions */

                  WDB[DATA_UFS_NZ] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                              1000);
                  WDB[DATA_UFS_ZMIN] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              -INFTY, INFTY);
                  WDB[DATA_UFS_ZMAX] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              -INFTY, INFTY);
                }
              else if (np == 12)
                {
                  /* Get mesh size and dimensions */

                  WDB[DATA_UFS_NX] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                              1000);
                  WDB[DATA_UFS_XMIN] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              -INFTY, INFTY);
                  WDB[DATA_UFS_XMAX] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              -INFTY, INFTY);

                  WDB[DATA_UFS_NY] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                              1000);
                  WDB[DATA_UFS_YMIN] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              -INFTY, INFTY);
                  WDB[DATA_UFS_YMAX] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              -INFTY, INFTY);

                  WDB[DATA_UFS_NZ] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                              1000);
                  WDB[DATA_UFS_ZMIN] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              -INFTY, INFTY);
                  WDB[DATA_UFS_ZMAX] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              -INFTY, INFTY);
                }
              else if ((long)RDB[DATA_UFS_MODE] != UFS_MODE_RMX)
                Error(-1, params[j], fname, line,
                      "Invalid number of parameters given");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "entr"))
            {
              /****** Fission source entropy *********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Number of x bins */

              if (k < np)
                WDB[DATA_ENTROPY_NX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 1000);

              /* Number of y bins */

              if (k < np)
                WDB[DATA_ENTROPY_NY] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 1000);

              /* Number of z bins */

              if (k < np)
                WDB[DATA_ENTROPY_NZ] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 1000);

              /* Limits in x-direction */

              if (k < np)
                WDB[DATA_ENTROPY_XMIN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            -INFTY, INFTY);
              if (k < np)
                WDB[DATA_ENTROPY_XMAX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            RDB[DATA_ENTROPY_XMIN], INFTY);

              /* Limits in y-direction */

              if (k < np)
                WDB[DATA_ENTROPY_YMIN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            -INFTY, INFTY);
              if (k < np)
                WDB[DATA_ENTROPY_YMAX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            RDB[DATA_ENTROPY_YMIN], INFTY);

              /* Limits in z-direction */

              if (k < np)
                WDB[DATA_ENTROPY_ZMIN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            -INFTY, INFTY);
              if (k < np)
                WDB[DATA_ENTROPY_ZMAX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            RDB[DATA_ENTROPY_ZMIN], INFTY);

              /***************************************************************/
            }
          else if (
                   !strcasecmp(params[j], "power") ||
                   !strcasecmp(params[j], "powdens") ||
                   !strcasecmp(params[j], "genrate") ||
                   !strcasecmp(params[j], "fissrate") ||
                   !strcasecmp(params[j], "absrate") ||
                   !strcasecmp(params[j], "lossrate") ||
                   !strcasecmp(params[j], "flux") ||
                   !strcasecmp(params[j], "sfrate") ||
                   !strcasecmp(params[j], "acti") ||
                   !strcasecmp(params[j], "srcrate"))
            {
              /****** Normalization ******************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j;

              /* Check number of parameters */

              if (np < 2)
                Error(-1, params[j], fname, line,
                      "Not enough parameters for normalization");

              /* Create new block */

              loc0 = NewItem(DATA_PTR_NORM, NORM_BLOCK_SIZE);

              /* Put name, file name and line number */

              WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(params[j]);
              WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
              WDB[loc0 + PARAM_LINE] = (double)line;

              /* Reset coefficients */

              WDB[loc0 + NORM_POWER] = -1.0;
              WDB[loc0 + NORM_POWDENS] = -1.0;
              WDB[loc0 + NORM_GENRATE] = -1.0;
              WDB[loc0 + NORM_FISSRATE] = -1.0;
              WDB[loc0 + NORM_ABSRATE] = -1.0;
              WDB[loc0 + NORM_LOSSRATE] = -1.0;
              WDB[loc0 + NORM_FLUX] = -1.0;
              WDB[loc0 + NORM_SRCRATE] = -1.0;
              WDB[loc0 + NORM_ACTI_ZAI] = -1.0;
              WDB[loc0 + NORM_SFRATE] = -1.0;

              /* Read value */

              val = TestParam(pname, fname, line, params[k + 1], PTYPE_REAL,
                              -1.0, INFTY);

              /* Check type and store value */

              if (!strcasecmp(params[k], "power"))
                WDB[loc0 + NORM_POWER] = val;
              else if (!strcasecmp(params[k], "powdens"))
                WDB[loc0 + NORM_POWDENS] = val;
              else if (!strcasecmp(params[k], "flux"))
                WDB[loc0 + NORM_FLUX] = val;
              else if (!strcasecmp(params[k], "genrate"))
                WDB[loc0 + NORM_GENRATE] = val;
              else if (!strcasecmp(params[k], "fissrate"))
                WDB[loc0 + NORM_FISSRATE] = val;
              else if (!strcasecmp(params[k], "absrate"))
                WDB[loc0 + NORM_ABSRATE] = val;
              else if (!strcasecmp(params[k], "lossrate"))
                WDB[loc0 + NORM_LOSSRATE] = val;
              else if (!strcasecmp(params[k], "srcrate"))
                WDB[loc0 + NORM_SRCRATE] = val;
              else if (!strcasecmp(params[k], "sfrate"))
                WDB[loc0 + NORM_SFRATE] = val;
              else if (!strcasecmp(params[k], "acti"))
                {
                  /* Check ZAI */

                  if ((val != -1) && ((val < 10010) || (val > 1120000)))
                    Error(loc0, "Invalid ZAI %s", params[k]);
                  else
                    WDB[loc0 + NORM_ACTI_ZAI] = val;

                  /* Get activity */

                  WDB[loc0 + NORM_ACTI_A] =
                    TestParam(pname, fname, line, params[k + 2], PTYPE_REAL,
                              ZERO, INFTY);

                  /* Increase count */

                  k++;

                  /* Initialize source rate */

                  WDB[loc0 + NORM_SRCRATE] = 1.0;
                }
              else
                Die(FUNCTION_NAME, "Error!");

              k = k + 2;

              /* Read material name */

              if (k < np)
                WDB[loc0 + NORM_PTR_MAT] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "pcc"))
            {
              /****** Predictor corrector method *****************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if (np < 2)
                Error(-1, params[j], fname, line,
                      "Predictor-corrector mode not given");

              /* Check mode */

              if (!strcasecmp(params[k], "ce") ||
                  !strcasecmp(params[k], "0"))
                {
                  /* CE (old "Euler's method") */

                  WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_CONSTANT;
                  WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_NONE;
                  WDB[DATA_BURN_PRED_NSS] = 1.0;
                  WDB[DATA_BURN_CORR_NSS] = -1.0;

                  /* Update index */

                  k++;

                  /* Check if substeps are given */

                  if (k < np)
                    Error(-1, pname, fname, line,
                          "Sub-steps not allowed in this mode");
                }

              else if (!strcasecmp(params[k], "celi") ||
                       !strcasecmp(params[k], "1"))
                {
                  /* CE/LI (old "predictor-corrector method") */

                  WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_CONSTANT;
                  WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_LINEAR;
                  WDB[DATA_BURN_PRED_NSS] = 1.0;
                  WDB[DATA_BURN_CORR_NSS] = 1.0;

                  /* Update index */

                  k++;

                  /* number of corrector substeps */

                  if (k < np)
                    WDB[DATA_BURN_CORR_NSS] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                                1000);

                  /* Check number of entries */

                  if (k < np)
                    Error(-1, pname, fname, line,
                          "Second sub-step entry not allowed");
                }

              else if (!strcasecmp(params[k], "le") ||
                       !strcasecmp(params[k], "2"))
                {
                  /* LE */

                  WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_LINEAR;
                  WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_NONE;
                  WDB[DATA_BURN_PRED_NSS] = 1.0;
                  WDB[DATA_BURN_CORR_NSS] = -1.0;

                  /* Update index */

                  k++;

                  /* number of predictor substeps */

                  if (k < np)
                    WDB[DATA_BURN_PRED_NSS] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                                1000);

                  /* Check number of entries */

                  if (k < np)
                    Error(-1, pname, fname, line,
                          "Second sub-step entry not allowed");
                }

              else if (!strcasecmp(params[k], "leli") ||
                       !strcasecmp(params[k], "3"))
                {
                  /* LE/LI */

                  WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_LINEAR;
                  WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_LINEAR;
                  WDB[DATA_BURN_PRED_NSS] = 1.0;
                  WDB[DATA_BURN_CORR_NSS] = 1.0;

                  /* Update index */

                  k++;

                  /* number of predictor substeps */

                  if (k < np)
                    WDB[DATA_BURN_PRED_NSS] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                                1000);

                  /* number of corrector substeps */

                  if (k < np)
                    WDB[DATA_BURN_CORR_NSS] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                                1000);
                }

              else if (!strcasecmp(params[k], "leqi") ||
                       !strcasecmp(params[k], "4"))
                {
                  /* LE/QI */

                  WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_LINEAR;
                  WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_QUADRATIC;
                  WDB[DATA_BURN_PRED_NSS] = 1.0;
                  WDB[DATA_BURN_CORR_NSS] = 1.0;

                  /* Update index */

                  k++;

                  /* number of predictor substeps */

                  if (k < np)
                    WDB[DATA_BURN_PRED_NSS] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                                1000);

                  /* number of corrector substeps */

                  if (k < np)
                    WDB[DATA_BURN_CORR_NSS] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                                1000);
                }
              else if (!strcasecmp(params[k], "cece") ||
                       !strcasecmp(params[k], "6"))
                {
                  /* CE/CE (corrector ce is constant backwards extrapolation)*/

                  WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_CONSTANT;
                  WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_CONSTANT;
                  WDB[DATA_BURN_PRED_NSS] = 1.0;
                  WDB[DATA_BURN_CORR_NSS] = 1.0;

                  /* Update index */

                  k++;

                }
              else if (!strcasecmp(params[k], "5"))
                {
                  /* Custom mode, update index */

                  k++;

                  /* Read predictor interpolation */

                  val = TestParam(pname, fname, line, params[k++], PTYPE_INT,
                                  0, 1);

                  if (val == 0)
                    WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_CONSTANT;
                  else
                    WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_LINEAR;

                  /* interpolation order on corrector */

                  if (k < np)
                    {
                      val = TestParam(pname, fname, line, params[k++],
                                      PTYPE_INT, 1, 2);

                      if(val == 1)
                        WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_LINEAR;
                      else
                        WDB[DATA_BURN_CORR_TYPE] =
                          (double)CORR_TYPE_QUADRATIC;
                    }

                  /* number of predictor substeps */

                  if (k < np)
                    WDB[DATA_BURN_PRED_NSS] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                                1000);

                  /* Number of corrector substeps */

                  if (k < np)
                    WDB[DATA_BURN_CORR_NSS] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                                1000);
                }
              else
                Error(-1, pname, fname, line,
                      "Invalid predictor-corrector mode %s", params[k]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ci"))
            {
              /****** Corrector iteration (for dufek/hybrid methods) *********/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if (np != 3)
                Error(-1, params[j], fname, line,
                      "Corrector iteration parameters missing");

              /* Check mode */

              if (!strcasecmp(params[k], "inner") ||
                  !strcasecmp(params[k], "0"))
                {
                  /* Inner iteration mode */

                  WDB[DATA_BURN_CI_TYPE] = (double)CI_TYPE_INNER;

                  /* Update index */

                  k++;

                }
              else if (!strcasecmp(params[k], "outer") ||
                       !strcasecmp(params[k], "1"))
                {
                  /* Outer iteration mode */

                  WDB[DATA_BURN_CI_TYPE] = (double)CI_TYPE_OUTER;

                  /* Update index */

                  k++;

                }
              else
                Error(-1, params[j], fname, line,
                      "Unknown corrector iteration mode %s", params[k]);

              /* Maximum number of iterations (ATM this is also the min.) */

              WDB[DATA_BURN_CI_MAXI] =
                TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                          1000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "sie"))
            {
              /****** Stochastic Implicit Euler ******/

             /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if (np != 2)
                Error(-1, params[j], fname, line,
                      "Stochastic Implicit Euler parameters missing");

              /* Set mode */

              WDB[DATA_BURN_SIE] = (double)YES;

              /* Inner corrector iteration mode */

              WDB[DATA_BURN_CI_TYPE] = (double)CI_TYPE_INNER;

              /* Maximum number of iterations (ATM this is also the min.) */

              WDB[DATA_BURN_CI_MAXI] =
                TestParam(pname, fname, line, params[k++], PTYPE_REAL, -1000.0,
                          1000.0);

              /* Constant extrapolation on predictor and corrector results */
              /* in Stochastic Implicit Euler */

              WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_CONSTANT;
              WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_CONSTANT;
              WDB[DATA_BURN_PRED_NSS] = 1.0;
              WDB[DATA_BURN_CORR_NSS] = 1.0;

              /* Check if tolerance given */

              if(RDB[DATA_BURN_CI_MAXI] < 0)
                {

                  WDB[DATA_BURN_CI_TOLER] = -RDB[DATA_BURN_CI_MAXI];

                  WDB[DATA_BURN_CI_MAXI] = 500;

                }
              else
                {
                  /* Number of iterations given */

                  WDB[DATA_BURN_CI_TOLER] = ZERO;

                }
            }
          else if (!strcasecmp(params[j], "cpop"))
            {
              /****** alternate population size in corrector iter. *********/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if ((np < 4) || (np > 5))
                Error(-1, params[j], fname, line,
                      "Corrector alt. pop. must have 3 or 4 parameters");

              /* Neutron population */

              WDB[DATA_BURN_CI_NBATCH] =
                TestParam(pname, fname, line,  params[k++], PTYPE_INT,
                          10, 100000000000);

              /* Number of active cycles */

              WDB[DATA_BURN_CI_CYCLES] =
                TestParam(pname, fname, line, params[k++], PTYPE_INT,
                          1, 100000000000);

              /* Number of inactive cycles */

              WDB[DATA_BURN_CI_SKIP] =
                TestParam(pname, fname, line, params[k++], PTYPE_INT,
                          0,  100000000000);

              if (k < np)
                {
                  /* Number of inactive cycles on further inner iterations */

                  WDB[DATA_BURN_CI_SKIP_2] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT,
                              0,  100000000000);

                }
              else
                {
                  /* Number of inactive cycles on further inner iterations */
                  /* same as for first inner iteration */

                  WDB[DATA_BURN_CI_SKIP_2] = RDB[DATA_BURN_CI_SKIP];

                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "inventory"))
            {
              /****** Output for burnup calculation **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of isotopes */

              ni = np - k;

              /* Read data */

              for (n = 0; n < ni; n++)
                {
                  /* Check top option */

                  if (!strcasecmp(params[k], "top"))
                    {
                      /* Update counter */

                      k++;

                      /* Get number of nuclides */

                      ns = (long)TestParam(pname, fname, line, params[k++],
                                           PTYPE_INT, 1, 100000);

                      /* Check type */

                      if (!strcasecmp(params[k], "mass"))
                        WDB[DATA_BURN_INV_TOP_MASS] = (double)ns;
                      else if (!strcasecmp(params[k], "activity"))
                        WDB[DATA_BURN_INV_TOP_ACTIVITY] = (double)ns;
                      else if (!strcasecmp(params[k], "sf"))
                        WDB[DATA_BURN_INV_TOP_SF] = (double)ns;
                      else if (!strcasecmp(params[k], "gsrc"))
                        WDB[DATA_BURN_INV_TOP_GSRC] = (double)ns;
                      else if ((!strcasecmp(params[k], "decheat")) ||
                               (!strcasecmp(params[k], "dh")))
                        WDB[DATA_BURN_INV_TOP_DECAY_HEAT] = (double)ns;
                      else if (!strcasecmp(params[k], "ingtox"))
                        WDB[DATA_BURN_INV_TOP_ING_TOX] = (double)ns;
                      else if (!strcasecmp(params[k], "inhtox"))
                        WDB[DATA_BURN_INV_TOP_INH_TOX] = (double)ns;
                      else if (!strcasecmp(params[k], "all"))
                        {
                          WDB[DATA_BURN_INV_TOP_MASS] = (double)ns;
                          WDB[DATA_BURN_INV_TOP_ACTIVITY] = (double)ns;
                          WDB[DATA_BURN_INV_TOP_SF] = (double)ns;
                          WDB[DATA_BURN_INV_TOP_GSRC] = (double)ns;
                          WDB[DATA_BURN_INV_TOP_DECAY_HEAT] = (double)ns;
                          WDB[DATA_BURN_INV_TOP_ING_TOX] = (double)ns;
                          WDB[DATA_BURN_INV_TOP_INH_TOX] = (double)ns;
                        }
                      else
                        Error(-1, params[j], fname, line,
                            "Unknown quantity %s for top inventory",
                            params[k]);

                      /* Update counters */

                      n = n + 2;
                      k++;
                    }
                  else
                    {
                      /* Create new block */

                      loc0 = NewItem(DATA_BURN_PTR_INVENTORY,
                                     INVENTORY_BLOCK_SIZE);
                      WDB[loc0 + INVENTORY_PTR_NAME]
                        = (double)PutText(params[k]);
                      WDB[loc0 + INVENTORY_PTR_ENTRY]
                        = (double)PutText(params[k++]);
                    }
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "blockdt"))
            {
              /****** List of material where DT is never used ****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of parameters */

              ni = np - k;

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni + 1);
              WDB[DATA_DT_PTR_BLOCK_LIST] = (double)loc0;

              /* Read data */

              for (n = 0; n < ni; n++)
                WDB[loc0++] = (double)PutText(params[k++]);

              /* Put null terminator */

              WDB[loc0] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "otfburn"))
            {
              /****** List of nuclides included in OTF-burnup ****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of parameters */

              ni = np - k;

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni + 1);
              WDB[DATA_OTF_BURN_NUC_LIST] = (double)loc0;

              /* Read data */

              for (n = 0; n < ni; n++)
                WDB[loc0++] = (double)PutText(params[k++]);

              /* Put null terminator */

              WDB[loc0] = NULLPTR;

              /* Put mode on */

              WDB[DATA_OTF_BURN_MODE] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "coefpara"))
            {
              /****** List of output parameters in coefficient calculations **/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get statistical error flag */

              n = (long)TestParam(pname, fname, line, params[k++],
                                  PTYPE_INT, 0, 1);

              /* Check */

              if (n == 0)
                WDB[DATA_COEF_CALC_INCLUDE_ERRORS] = (double)NO;
              else
                WDB[DATA_COEF_CALC_INCLUDE_ERRORS] = (double)YES;

              /* Get number of parameters */

              ni = np - k;

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni + 1);
              WDB[DATA_COEF_CALC_PTR_PARAM_LIST] = (double)loc0;

              /* Read data */

              for (n = 0; n < ni; n++)
                WDB[loc0++] = (double)PutText(params[k++]);

              /* Put null terminator */

              WDB[loc0] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "forcedt"))
            {
              /****** List of material where DT is always used ***************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of parameters */

              ni = np - k;

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni);
              WDB[DATA_DT_PTR_FORCE_LIST] = (double)loc0;

              /* Read data */

              for (n = 0; n < ni; n++)
                WDB[loc0++] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "runtme"))
            {
              /****** Maximum transport cycle running time *******************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get value */

              WDB[DATA_MAX_TRANSPORT_RUNTIME] =
                TestParam(pname, fname, line, params[k++],
                          PTYPE_REAL, -1.0, INFTY);

              /* Convert to seconds */

              WDB[DATA_MAX_TRANSPORT_RUNTIME] =
                60.0*RDB[DATA_MAX_TRANSPORT_RUNTIME];

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "mvol"))
            {
              /****** Material volumes ***************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of entries */

              ni = np - k;

              /* Check */

              if (ni%3)
                Error(-1, params[j], fname, line,
                      "Invalid number of parameters");
              else
                ni = (long)((double)ni/3.0);

              /* Loop over entries */

              for (n = 0; n < ni; n++)
                {
                  /* Create new block */

                  loc0 = NewItem(DATA_PTR_MVOL0, MVOL_BLOCK_SIZE);

                  /* Put data */

                  WDB[loc0 + MVOL_PTR_MAT] = (double)PutText(params[k++]);
                  WDB[loc0 + MVOL_REG_IDX] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT,
                              0, 100000000000);
                  WDB[loc0 + MVOL_VOL] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              0.0, INFTY);
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "his"))
            {
              /***** Output file control *************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Print history file */

              if (k < np)
                WDB[DATA_OPTI_PRINT_HIS] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Switch fission source entropy calculation on */

              WDB[DATA_OPTI_ENTROPY_CALC] = (double)YES;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "tcut"))
            {
              /***** Time cut-off ********************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get limit */

              if (k < np)
                WDB[DATA_TIME_CUT_TMAX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, INFTY);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ecut"))
            {
              /***** Energy cut-off ******************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get limit for neutrons */

              if (k < np)
                {
                  val = TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                                  -1.0, INFTY);

                  if (val >= 0.0)
                    WDB[DATA_NEUTRON_ECUT] = val;
                }

              /* Get limit for photons */

              if (k < np)
                WDB[DATA_PHOTON_ECUT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, RDB[DATA_PHOTON_EMAX]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ecutmat"))
            {
              /****** Photon energy cut-offs in materials ********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of entries */

              ni = np - k;

              /* Check */

              if ((ni == 0) || (ni%2))
                Error(-1, pname, fname, line,
                      "Invalid number of parameters");
              else
                ni = (long)((double)ni/2.0);

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni+1);
              loc1 = ReallocMem(DATA_ARRAY, ni+1);

              /* Set material and energy pointers */

              WDB[DATA_PHOTON_ECUTMAT_MAT] = (double)loc0;
              WDB[DATA_PHOTON_ECUTMAT_E] = (double)loc1;

              /* Loop over entries */

              for (n = 0; n < ni; n++)
                {
                  /* Put data */
                  WDB[loc0++] = (double)PutText(params[k++]);
                  WDB[loc1++] = TestParam(pname, fname, line, params[k++],
                     PTYPE_REAL, RDB[DATA_PHOTON_EMIN], RDB[DATA_PHOTON_EMAX]);
                }

              /* Put null terminator */

              WDB[loc0] = NULLPTR;
              WDB[loc1] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ecutdens"))
            {
              /****** Photon energy cut-offs for mass densities **************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of entries */

              ni = np - k;

              /* Check */

              if ((ni == 0) || (ni%2))
                Error(-1, pname, fname, line,
                      "Invalid number of parameters");
              else
                ni = (long)((double)ni/2.0);

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni+1);
              loc1 = ReallocMem(DATA_ARRAY, ni+1);

              /* Set density and energy pointers */

              WDB[DATA_PHOTON_ECUTDENS_DENS] = (double)loc0;
              WDB[DATA_PHOTON_ECUTDENS_E] = (double)loc1;

              /* Loop over entries */

              for (n = 0; n < ni; n++)
                {
                  /* Put data */
                  WDB[loc0++] = TestParam(pname, fname, line, params[k++],
                     PTYPE_REAL, ZERO, INFTY);
                  WDB[loc1++] = TestParam(pname, fname, line, params[k++],
                     PTYPE_REAL, RDB[DATA_PHOTON_EMIN], RDB[DATA_PHOTON_EMAX]);

                  /* Check density */
                  if ((n > 0) && (WDB[loc0-1] <= WDB[loc0-2]))
                    Error(-1, pname, fname, line, "Densities are not in "
                                                  "ascending order");

                  /* Check energy */
                  if ((n > 0) && (WDB[loc1-1] <= WDB[loc1-2]))
                    Error(-1, pname, fname, line, "Cut-off energies are not in "
                                                  "ascending order");
                }

              /* Put null terminator */

              WDB[loc0] = NULLPTR;
              WDB[loc1] = NULLPTR;

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "mfpcut"))
            {
              /***** Photon mean free path cut-off ***************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get limit */

              if (k < np)
                WDB[DATA_PHOTON_MFPCUT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            ZERO, INFTY);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "mfpcutmat"))
            {
              /****** Photon mean free path cut-offs in materials ************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of entries */

              ni = np - k;

              /* Check */

              if ((ni == 0) || (ni%2))
                Error(-1, pname, fname, line,
                      "Invalid number of parameters");
              else
                ni = (long)((double)ni/2.0);

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni+1);
              loc1 = ReallocMem(DATA_ARRAY, ni+1);

              /* Set material and mean free path pointers */

              WDB[DATA_PHOTON_MFPCUTMAT_MAT] = (double)loc0;
              WDB[DATA_PHOTON_MFPCUTMAT_MFP] = (double)loc1;

              /* Loop over entries */

              for (n = 0; n < ni; n++)
                {
                  /* Put data */
                  WDB[loc0++] = (double)PutText(params[k++]);
                  WDB[loc1++] = TestParam(pname, fname, line, params[k++],
                     PTYPE_REAL, ZERO, INFTY);
                }

              /* Put null terminator */

              WDB[loc0] = NULLPTR;
              WDB[loc1] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "mfpcutdens"))
            {
              /****** Photon mean free path cut-offs for mass densities ******/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of entries */

              ni = np - k;

              /* Check */

              if ((ni == 0) || (ni%2))
                Error(-1, pname, fname, line,
                      "Invalid number of parameters");
              else
                ni = (long)((double)ni/2.0);

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni+1);
              loc1 = ReallocMem(DATA_ARRAY, ni+1);

              /* Set density and mean free path pointers */

              WDB[DATA_PHOTON_MFPCUTDENS_DENS] = (double)loc0;
              WDB[DATA_PHOTON_MFPCUTDENS_MFP] = (double)loc1;

              /* Loop over entries */

              for (n = 0; n < ni; n++)
                {
                  /* Put data */
                  WDB[loc0++] = TestParam(pname, fname, line, params[k++],
                     PTYPE_REAL, ZERO, INFTY);
                  WDB[loc1++] = TestParam(pname, fname, line, params[k++],
                     PTYPE_REAL, ZERO, INFTY);

                  /* Check density */
                  if ((n > 0) && (WDB[loc0-1] <= WDB[loc0-2]))
                    Error(-1, pname, fname, line, "Densities are not in "
                                                  "ascending order");

                  /* Check energy */
                  if ((n > 0) && (WDB[loc1-1] <= WDB[loc1-2]))
                    Error(-1, pname, fname, line, "Mean free path cut-offs "
                                                 "are not in ascending order");
                }

              /* Put null terminator */

              WDB[loc0] = NULLPTR;
              WDB[loc1] = NULLPTR;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "gcut"))
            {
              /***** Generation cut-off **************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get limit */

              if (k < np)
                {
                  WDB[DATA_GEN_CUT] =
                    TestParam(pname, fname, line, params[k++], PTYPE_INT,
                              1, MAX_GENERATIONS);

                  WDB[DATA_MAX_PROMPT_CHAIN_LENGTH] = RDB[DATA_GEN_CUT];
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "fpcut"))
            {
              /***** Fission product yield cut-off ***************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get limit */

              if (k < np)
                WDB[DATA_DEP_FP_YIELD_CUTOFF] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,  0.0,
                            2.0);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "printm"))
            {
              /***** Print compositions flag *********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              if (k < np)
                WDB[DATA_BURN_PRINT_COMP] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Limit */

              if (k < np)
                WDB[DATA_BURN_PRINT_COMP_LIM] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL, 0.0,
                            1.0);

              /* Print a note */

              fprintf(outp, "\nNOTE: The \"set printm\" feature is outdated ");
              fprintf(outp, "and not recommended to be\n      used anymore. ");
              fprintf(outp, "For restart calculations use the binary restart");
              fprintf(outp, "\n      files (\"set rfw\") and (\"set rfr\"):");
              fprintf(outp, "\n\n      http://serpent.vtt.fi/mediawiki/index.php/Input_syntax_manual#set_rfw\n\n");
              fprintf(outp, "      Isotopic data is best obtained");
              fprintf(outp, " from the _dep.m output.\n\n");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "depmtx"))
            {
              /***** Print burnup matrix flag ********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              if (k < np)
                WDB[DATA_BURN_PRINT_DEPMTX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "outp"))
            {
              /***** Output print interval ***********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              if (k < np)
                WDB[DATA_PRINT_INTERVAL] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1,
                            10000000);
              else
                Error(-1, pname, fname, line, "Missing print interval");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "micro"))
            {
              /****** Micro-group structure **********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Energy grid name */

              if (k < np)
                WDB[DATA_MICRO_PTR_EGRID] = (double)PutText(params[k++]);

              /* Get micro-group calculation option */

              if (k < np)
                if (TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 2)
                    == 2)
                  WDB[DATA_MICRO_CALC_BATCH_SIZE] = -1.0;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "gspec"))
            {
              /****** Energy group structure for gamma spectra ***************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Energy grid name */

              if (k < np)
                WDB[DATA_GSPEC_PTR_EGRID] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "cmm"))
            {
              /****** Switch CMM calculation on or off ***********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get option */

              if (k < np)
                WDB[DATA_CMM_CALC] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "dix"))
            {
              /****** Switch double-indexing on or off ***********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get option */

              if (k < np)
                WDB[DATA_OPTI_DIX] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "fum"))
            {
              /****** Critical spectrum calculation **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Set calculation on */

              WDB[DATA_B1_CALC] = (double)YES;

              /* Energy grid name */

              if (k < np)
                WDB[DATA_MICRO_PTR_EGRID] = (double)PutText(params[k++]);

              /* Get micro-group calculation option */

              if (k < np)
                if (TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 2)
                    == 2)
                  WDB[DATA_MICRO_CALC_BATCH_SIZE] = -1.0;

              /* Check remaining */

              if (k < np)
                {
                  /* Check mode */

                  if ((params[k][0] == 'o') || (params[k][0] == 'O') ||
                      (params[k][0] == '0'))
                    {
                      /* old method */

                      WDB[DATA_B1_MODE] = (double)CRIT_SPECTRUM_OLD;
                      k++;
                    }
                  else if ((params[k][0] == 'b') || (params[k][0] == 'B') ||
                      (params[k][0] == '1'))
                    {
                      /* B1 method */

                      WDB[DATA_B1_MODE] = (double)CRIT_SPECTRUM_B1;
                      k++;
                    }
                  else if ((params[k][0] == 'p') || (params[k][0] == 'P') ||
                           (params[k][0] == '2'))
                    {
                      /* P1 method */

                      WDB[DATA_B1_MODE] = (double)CRIT_SPECTRUM_P1;
                      k++;
                    }
                  else if ((params[k][0] == 'f') || (params[k][0] == 'F') ||
                           (params[k][0] == '3'))
                    {
                      /* "Casmo" method */

                      WDB[DATA_B1_MODE] = (double)CRIT_SPECTRUM_FM;
                      k++;

                      /* Get diffusion coefficient */

                      if (k < np)
                        WDB[DATA_B1_FM_DIFF] =
                          (double)TestParam(pname, fname, line, params[k++],
                                            PTYPE_INT, 1, 3);
                    }
                  else
                    Error(-1, params[j], fname, line, "Invalid mode: %s",
                          params[k]);
                }

              /* Get error limit */

              if (k < np)
                WDB[DATA_B1_ERR_LIMIT] =
                  TestParam(pname, fname, line, params[k++],
                            PTYPE_REAL, 0.0, 1.0);

              /* Get target keff */

              if (k < np)
                WDB[DATA_B1_KEFF_TGT] =
                  TestParam(pname, fname, line, params[k++],
                            PTYPE_REAL, ZERO, 4.0);

              /* Get maximum number of iterations */

              if (k < np)
                WDB[DATA_B1_MAX_ITER] =
                  TestParam(pname, fname, line, params[k++],
                            PTYPE_INT, 2, 10000);

              /* Get initial guess for B2 */

              if (k < np)
                WDB[DATA_B1_B2_INIT] =
                  TestParam(pname, fname, line, params[k++],
                            PTYPE_REAL, ZERO, INFTY);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "bumode"))
            {
              /***** Burnup calculation mode *********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Read mode */

              if (!strcasecmp(params[k], "tta"))
                WDB[DATA_BURN_BUMODE] = (double)BUMODE_TTA;
              else if (!strcasecmp(params[k], "cram"))
                WDB[DATA_BURN_BUMODE] = (double)BUMODE_CRAM;
              else
                {
                  /* Read mode */

                  WDB[DATA_BURN_BUMODE] =
                    TestParam(pname, fname, line, params[k], PTYPE_INT, 1, 3);

                  /* Convert mode 3 to 1 for compatibility with Serpent 1 */

                  if (RDB[DATA_BURN_BUMODE] == 3.0)
                    WDB[DATA_BURN_BUMODE] = (double)BUMODE_TTA;
                }

              /* Update index */

              k++;

              /* Check if cram order is given */

              if (k < np)
                {
                  /* Check mode */

                  if ((long)RDB[DATA_BURN_BUMODE] == BUMODE_CRAM)
                    {
                      /* Get order */

                      n = (long)TestParam(pname, fname, line, params[k++],
                                          PTYPE_INT, 2, 16);

                      /* Check */

                      if (n % 2)
                        Error(0, "Invalid CRAM order must be divisible by 2");
                      else
                        WDB[DATA_BURN_CRAM_K] = (double)n;
                    }
                }

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ttacut"))
            {
              /***** TTA cut-off *********************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get limit */

              if (k < np)
                WDB[DATA_DEP_TTA_CUTOFF] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL, 0.0,
                            1.0);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "bunorm"))
            {
              /***** Normalization in burnup calculation *********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_NORM_BURN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 1, 3);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "nbuf"))
            {
              /***** Neutron buffer factor ***********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Factor */

              if (k < np)
                WDB[DATA_PART_NBUF_FACTOR] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            1.0, 1000000000.0);

              /* Event bank */

              if (k < np)
                WDB[DATA_EVENT_BANK_SZ] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 10000000000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "gbuf"))
            {
              /***** Photon buffer factor ************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Factor */

              if (k < np)
                WDB[DATA_PART_GBUF_FACTOR] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            1.0, 1000000000.0);

              /* Event bank */

              if (k < np)
                WDB[DATA_EVENT_BANK_SZ] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            10, 10000000000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "pbuf"))
            {
              /***** Precursor buffer factor *********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_PART_PBUF_FACTOR] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            1.0, 1000000000.0);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "sfbuf"))
            {
              /***** Source file buffer **************************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_SRC_FILE_BUF_SIZE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 1000000000);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "shbuf"))
            {
              /***** Shared scoring buffer ***********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Option for BUF array */

              if (k < np)
                WDB[DATA_OPTI_SHARED_BUF] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /* Option for RES2 array */

              if (k < np)
                WDB[DATA_OPTI_SHARED_RES2] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ppid"))
            {
              /***** User given parent process id *****************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_PPID] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 0,
                            1000000);

              /* Set coupled calculation mode on */

              WDB[DATA_RUN_CC] = (double)YES;

              /* Set signaling mode to file */

              if(RDB[DATA_CC_SIG_MODE] == (double)SIG_MODE_NONE)
                WDB[DATA_CC_SIG_MODE] = (double)SIG_MODE_POSIX;
              else
                Error(-1, params[j], fname, line,
                      "Multiple signalling modes defined (set ppid and set comfile)");

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "pport"))
            {
              /***** User given parent process port ***************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_PPORT] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT, 0,
                            1000000);

              /* Set coupled calculation mode on */

              WDB[DATA_RUN_CC] = (double)YES;

              /* Set signaling mode to file */

              if(RDB[DATA_CC_SIG_MODE] == (double)SIG_MODE_NONE)
                WDB[DATA_CC_SIG_MODE] = (double)SIG_MODE_SOCKET;
              else
                Error(-1, params[j], fname, line,
                      "Multiple signalling modes defined (set ppid, set comfile, set pport)");

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "fsp"))
            {
              /***** Fission source passing between iterations / dep steps ****/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if(k < np)
                WDB[DATA_USE_FSP] = TestParam(pname, fname, line, params[k++],
                                              PTYPE_LOGICAL);

              /* Read alternate number of inactive cycles */

              if (RDB[DATA_USE_FSP] == (double)YES)
                {
                  if(k == np)
                    Error(-1, params[j], fname, line,
                          "Missing number of critical cycles with fission source passing");
                  else
                    WDB[DATA_FSP_CRIT_SKIP] =
                      TestParam(pname, fname, line, params[k++], PTYPE_INT, 0, 100000000);

                  /* Store number of inactive cycles so that it can be */
                  /* alternated with SIE burnup scheme correctors      */

                  WDB[DATA_BURN_CI_ORIG_FSP_SKIP] = RDB[DATA_FSP_CRIT_SKIP];
                }

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "ekn"))
            {
              /***** Klein-Nishina threshold energy ***************************/

              /* (scattering function used below the threshold) */

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_PHOTON_EKN] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL, -INFTY,
                            INFTY);

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "cdop"))
            {
              /***** Doppler broadening of Compton photons ********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_PHOTON_USE_DOPPLER] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "fallback_cp_2129"))
           {
             /***** Fallback option for Compton profiles *********************/

             /* Copy parameter name */

             strcpy (pname, params[j]);

             k = j + 1;

             /* Mode */

             if (k < np)
               WDB[DATA_PHOTON_CP_FALLBACK_2129] =
                 TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

             /****************************************************************/
           }
          else if (!strcasecmp(params[j], "ttb"))
            {
              /***** Thick-target bremsstralung for electrons and positrons **/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_PHOTON_USE_TTB] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "ttbpm"))
            {
              /***** Bremsstralung model for positrons ************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_PHOTON_TTBPM] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "pdatadir"))
            {
              /***** Directory for photon transport data **********************/

              k = j + 1;

              /* Check the last character */

              if (params[k][strlen(params[k]) - 1] != '/')
                strcat(params[k], "/");

              /* Check the first character and set the path accordingly */

              if (params[k][0] == '/') {
                /* Assuming absolute path */
                WDB[DATA_PHOTON_DATA_DIR] = (double)PutText(params[k]);
              }
              else {
                if ((getenv("SERPENT_DATA")) != NULL) {
                  /* Assuming relative path with respect to the environment
                   * variable */
                  sprintf(tmpstr, "%s/%s", getenv("SERPENT_DATA"), params[k]);
                  WDB[DATA_PHOTON_DATA_DIR] = (double)PutText(tmpstr);
                }
                else {
                  /* The environment variable is not set - try absolute path */
                  sprintf(tmpstr, "/%s", params[k]);
                  WDB[DATA_PHOTON_DATA_DIR] = (double)PutText(tmpstr);
                }
              }

              k++;

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "cea"))
            {
              /***** Compton electron angular distribution model **************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_PHOTON_COMP_EANG] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "elmee"))
            {
              /****** List of materials for elmee card ***********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of entries */

              ni = np - k;

              /* Check */

              if ((ni == 0) || (ni%2))
                Error(-1, pname, fname, line,
                      "Invalid number of parameters");
              else
                ni = (long)((double)ni/2.0);

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni+1);
              loc1 = ReallocMem(DATA_ARRAY, ni+1);

              WDB[DATA_ELECTRON_MEE_MAT] = (double)loc0;
              WDB[DATA_ELECTRON_MEE] = (double)loc1;

              /* Loop over entries */

              for (n = 0; n < ni; n++)
                {
                  /* Put data */
                  WDB[loc0++] = (double)PutText(params[k++]);
                  WDB[loc1++] = TestParam(pname, fname, line, params[k++],
                      PTYPE_REAL, 0.0, 1.0);  /* 1 MeV maximum */
                }

              /* Put null terminator */

              WDB[loc0] = NULLPTR;
              WDB[loc1] = NULLPTR;

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "elcond"))
            {
              /****** List of materials for elcond card **********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of entries */

              ni = np - k;

              /* Check */

              if ((ni == 0) || (ni%2))
                Error(-1, pname, fname, line,
                      "Invalid number of parameters");
              else
                ni = (long)((double)ni/2.0);

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni+1);
              loc1 = ReallocMem(DATA_ARRAY, ni+1);

              WDB[DATA_ELECTRON_COND_MAT] = (double)loc0;
              WDB[DATA_ELECTRON_COND] = (double)loc1;

              /* Loop over entries */

              for (n = 0; n < ni; n++)
                {
                  /* Put data */
                  WDB[loc0++] = (double)PutText(params[k++]);
                  WDB[loc1++] = TestParam(pname, fname, line, params[k++],
                      PTYPE_INT, 0, 2);
                }

              /* Put null terminator */

              WDB[loc0] = NULLPTR;
              WDB[loc1] = NULLPTR;

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "elgas"))
            {
              /****** List of materials for elgas card ***********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of entries */

              ni = np - k;

              /* Check the number of entries */

              if ((ni == 0) || (ni%2))
                Error(-1, pname, fname, line,
                      "Invalid number of parameters");
              else
                ni = (long)((double)ni/2.0);

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni+1);
              loc1 = ReallocMem(DATA_ARRAY, ni+1);

              WDB[DATA_ELECTRON_GAS_MAT] = (double)loc0;
              WDB[DATA_ELECTRON_GAS] = (double)loc1;

              /* Loop over entries */

              for (n = 0; n < ni; n++)
                {
                  /* Put data */
                  WDB[loc0++] = (double)PutText(params[k++]);
                  WDB[loc1++] = TestParam(pname, fname, line, params[k++],
                      PTYPE_LOGICAL);
                }

              /* Put null terminator */

              WDB[loc0] = NULLPTR;
              WDB[loc1] = NULLPTR;

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "printelsp"))
            {
              /***** Print stopping power data flag **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                WDB[DATA_ELECTRON_PRINT_SP] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "elspn"))
            {
              /***** Stopping power energy grid size *************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check that elspe has not been set */

              if (RDB[DATA_ELECTRON_SP_E] > VALID_PTR)
                Error(-1, pname, fname, line,
                      "\"%s\" can't be used if \"elspe\" has been set",
                      pname);

              /* Mode */

              if (k < np)
                WDB[DATA_ELECTRON_SP_N] = TestParam(pname, fname, line,
                                            params[k++], PTYPE_INT, 50, 5000);

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "elspe"))
            {
              /***** Stopping power energy grid ******************************/
              /* NOTE: This option is used for debugging only */

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get number of entries */

              ni = np - k;

              /* Check */

              if ((ni < 5) || (ni > 5000))
                Error(-1, pname, fname, line,
                      "Invalid number of parameters (should be between 5 "
                      "and 5000)");

              WDB[DATA_ELECTRON_SP_N] = (double)ni;

              /* Allocate memory */

              loc0 = ReallocMem(DATA_ARRAY, ni+1);

              WDB[DATA_ELECTRON_SP_E] = (double)loc0;

              /* Loop over entries */

              for (n = 0; n < ni; n++) {
                  WDB[loc0++] =  TestParam(pname, fname, line, params[k++],
                      PTYPE_REAL, RDB[DATA_PHOTON_EMIN], 1.0e4);
              }

              /* Put null terminator */

              WDB[loc0] = NULLPTR;


              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "comfile"))
            {
              /***** User given communication file names **********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get com in file */

              if (k < np)
                {
                  WDB[DATA_PTR_COM_IN] = (double)PutText(params[k++]);
                }
              else
                Error(-1, params[j], fname, line,
                      "Missing communications filename to Serpent");

              /* Get com out file */

              if (k < np)
                {
                  WDB[DATA_PTR_COM_OUT] = (double)PutText(params[k++]);
                }
              else
                Error(-1, params[j], fname, line,
                      "Missing communications filename from Serpent");


              /* Reset infile */

              if ((fp = fopen(GetText(DATA_PTR_COM_IN), "w")) == NULL)
                Error(-1, "Cannot open in communication file \"%s\" for writing",
                      GetText(DATA_PTR_COM_IN));
              fprintf(fp, "-1\n");
              fclose(fp);

              /* Reset outfile */

              if ((fp = fopen(GetText(DATA_PTR_COM_OUT), "w")) == NULL)
                Error(-1, "Cannot open out communication file \"%s\" for writing",
                      GetText(DATA_PTR_COM_OUT));
              fprintf(fp, "-1\n");
              fclose(fp);

              /* Set coupled calculation mode on */

              WDB[DATA_RUN_CC] = (double)YES;

              /* Set signaling mode to file */

              if(RDB[DATA_CC_SIG_MODE] == (double)SIG_MODE_NONE)
                WDB[DATA_CC_SIG_MODE] = (double)SIG_MODE_FILE;
              else
                Error(-1, params[j], fname, line,
                      "Multiple signalling modes defined (set ppid and set comfile)");

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "relfactor"))
            {
              /***** Solution relaxation factor for coupled calculation *******/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                {
                  WDB[DATA_SOL_REL_FACT] = TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 1.0);
                }
              else
                Error(-1, params[j], fname, line,
                      "Missing solution relaxation factor");

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "memfrac"))
            {
              /***** Solution relaxation factor for coupled calculation *******/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                {
                  WDB[DATA_CPU_MEM_FRAC] = TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 1.0);
                }
              else
                Error(-1, params[j], fname, line,
                      "Missing value after memfrac");

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "finrodfile"))
            {
              /***** FINIX rod definition file ********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Store path */

              WDB[DATA_PTR_FINROD_FNAME] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "finscenariofile"))
            {
              /***** FINIX scenario definition file **************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Store path */

              WDB[DATA_PTR_FINSCEN_FNAME] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "finoptionsfile"))
            {
              /***** FINIX options definition file ***************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Store path */

              WDB[DATA_PTR_FINOPTI_FNAME] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "fininitfile"))
            {
              /***** FINIX options definition file ***************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Store path */

              WDB[DATA_PTR_FININIT_FNAME] = (double)PutText(params[k++]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ccmaxiter"))
            {
              /***** Solution relaxation factor for coupled calculation *******/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                {
                  WDB[DATA_SOL_REL_MAX_ITER] = TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 1000000);
                }
              else
                Error(-1, params[j], fname, line,
                      "Missing value after memfrac");

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "ccmaxpop"))
            {
              /***** Solution relaxation factor for coupled calculation *******/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                {
                  WDB[DATA_SOL_REL_MAX_POP] = TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, INFTY);
                }
              else
                Error(-1, params[j], fname, line,
                      "Missing value after memfrac");

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "savesrc"))
            {

              /** File name for precursor concentrations **/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if (k == np)
                Error(-1, params[j], fname, line, "Invalid number of parameters for savesrc");

              if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
                {
                  /* Allocate memory for precursor detector block if not already created */
                  /* This might have been created with set dynsrc */

                  /* Create precursor detector */

                  loc0 = NewItem(DATA_PTR_PREC_DET, PRECDET_BLOCK_SIZE);

                  WDB[loc0 + PARAM_PTR_NAME] = (double)PutText(word);
                  WDB[loc0 + PARAM_PTR_FNAME] = (double)PutText(fname);
                  WDB[loc0 + PARAM_LINE] = (double)line;

                  /* Initialize meshing */

                  WDB[loc0 + PRECDET_N0] = 1.0;
                  WDB[loc0 + PRECDET_N1] = 1.0;
                  WDB[loc0 + PRECDET_N2] = 1.0;

                  /* Initialize fractions of points to save */

                  WDB[loc0 + PRECDET_SAVE_FRAC_LIVE] = 1.0;
                  WDB[loc0 + PRECDET_SAVE_FRAC_PREC] = 1.0;
                }

              /* Store output file name base */
              WDB[loc0 + PRECDET_PTR_OUT_FNAME] = (double)PutText(params[k++]);

              /* Fraction of live neutrons to store */
              if (k < np)
                WDB[loc0 + PRECDET_SAVE_FRAC_LIVE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 1.0);

              /* Fraction of precursors to store */
              if (k < np)
                WDB[loc0 + PRECDET_SAVE_FRAC_PREC] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 1.0);

              /* Read mesh dimensions if given (default is 1x1x1) */
              /* Mesh dimension in x-direction */
              if (k < np)
                WDB[loc0 + PRECDET_N0] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 10000);

              /* Mesh dimension in y-direction */
              if (k < np)
                WDB[loc0 + PRECDET_N1] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 10000);

              /* Mesh dimension in z-direction */
              if (k < np)
                WDB[loc0 + PRECDET_N2] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 10000);

              /****************************************************************/
            }
          else if (!strcasecmp(params[j], "dynsrc"))
            {
              /****************************************************************/

              /** File name for live+delayed neutron source **/
              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if (k == np)
                Error(-1, params[j], fname, line, "Invalid number of parameters for dynsrc");

              /* Set neutron transport mode on */

              WDB[DATA_NEUTRON_TRANSPORT_MODE] = (double)YES;

              if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
                {
                  /* Allocate memory for precursor detectors if not already created */
                  /* This might have been created with set savesrc */

                  /* Create precursor detector */

                  loc0 = NewItem(DATA_PTR_PREC_DET, PRECDET_BLOCK_SIZE);

                  /* Initialize meshing */

                  WDB[loc0 + PRECDET_N0] = 1.0;
                  WDB[loc0 + PRECDET_N1] = 1.0;
                  WDB[loc0 + PRECDET_N2] = 1.0;

                  /* Initialize fractions of points to save */

                  WDB[loc0 + PRECDET_SAVE_FRAC_LIVE] = 1.0;
                  WDB[loc0 + PRECDET_SAVE_FRAC_PREC] = 1.0;

                }

              /* Store input file name base */

              WDB[loc0 + PRECDET_PTR_IN_FNAME] = (double)PutText(params[k++]);

              /* Test that source files exist */

              sprintf(tmpstr, "%s.main", GetText(loc0 + PRECDET_PTR_IN_FNAME));

              if ((fp = fopen(tmpstr, "r")) != NULL)
                  fclose(fp);
              else
                {
                  /* File not found */

                  Error(-1, pname, fname, line,
                        "Dynamic source file \"%s\" does not exist", tmpstr);
                }

              sprintf(tmpstr, "%s.prec", GetText(loc0 + PRECDET_PTR_IN_FNAME));

              if ((fp = fopen(tmpstr, "r")) != NULL)
                  fclose(fp);
              else
                {
                  /* File not found */

                  Error(-1, pname, fname, line,
                        "Dynamic source file \"%s\" does not exist", tmpstr);
                }

              sprintf(tmpstr, "%s.live", GetText(loc0 + PRECDET_PTR_IN_FNAME));

              if ((fp = fopen(tmpstr, "r")) != NULL)
                  fclose(fp);
              else
                {
                  /* File not found */

                  Error(-1, pname, fname, line,
                        "Dynamic source file \"%s\" does not exist", tmpstr);
                }

              /* Set precursor transport mode on if requested */
              /* Default is 0 == use regular precursor mesh */

              if (k < np)
                {
                  WDB[DATA_PRECURSOR_TRANSPORT_MODE] =
                    TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

                  /* Switch to internal PREC_MODE_FLAGS */

                  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == 0.0)
                    WDB[DATA_PRECURSOR_TRANSPORT_MODE] = PREC_MODE_MESH;
                  else
                    WDB[DATA_PRECURSOR_TRANSPORT_MODE] = PREC_MODE_POINT;
                }

              /* Set precursor buffer factor */

              if ((RDB[DATA_PRECURSOR_TRANSPORT_MODE] != PREC_MODE_NONE) &&
                  (WDB[DATA_PART_PBUF_FACTOR] < 1.0))
                WDB[DATA_PART_PBUF_FACTOR] = 1.0;

              /* Set delayed neutron emission off at fission, but include */
              /* delayed nu in total nu */

              WDB[DATA_USE_DELNU] = DELNU_MODE_PREC;

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "precsrcf"))
            {
              /***** Precursor number of source points factor ****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                {
                  WDB[DATA_PREC_SRC_FACT] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              0.0, 1000000.0);
                }
              else
                Error(-1, params[j], fname, line,
                      "Missing precursor source point factor");

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "precthresh"))
            {
              /***** Threshold for saving additional precursors in transient **/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                {
                  WDB[DATA_PREC_STORE_TRESH] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              -100.0, 100.0);
                }
              else
                Error(-1, params[j], fname, line,
                      "Missing precursor storing treshold factor");

              /* Check that the threshold was nonzero */

              if (fabs(RDB[DATA_PREC_STORE_TRESH]) < ZERO)
                Error(-1, params[j], fname, line,
                      "Need nonzero precursor threshold. Was given %E.",
                      RDB[DATA_PREC_STORE_TRESH]);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "ifp"))
            {
              /****** Iterated fission probability chain length **************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get mode */

              if (k < np)
                WDB[DATA_IFP_CHAIN_LENGTH] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            1, 100);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "transtime"))
            {
              /****** Time to use for moving transformations *****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Get mode */

              if (k < np)
                WDB[DATA_TDEP_TRANS_TYPE] =
                  TestParam(pname, fname, line, params[k++], PTYPE_INT,
                            0, 1);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "transnorm"))
            {
              /***** Additional normalization factor for transients **********/
              /* Default is 1.0 */

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Mode */

              if (k < np)
                {
                  WDB[DATA_TRANS_NORM_FACT] =
                    TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                              ZERO, 1e12);
                }
              else
                Error(-1, params[j], fname, line,
                      "Missing transient normalization factor");
              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "edepmode"))
            {
              /***** Energy deposition mode **********************************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if (k == np)
                Error(-1, params[j], fname, line,
                      "Invalid number of parameters for edepmode");

              /* Store energy deposition mode */

              WDB[DATA_EDEP_MODE] = TestParam(pname, fname, line, params[k++],
                                              PTYPE_INT, 0, 3);

              /* Set flag to read kerma from ace file if needed */

              if ((long)RDB[DATA_EDEP_MODE] > EDEP_MODE_MT458)
                WDB[DATA_INCLUDE_HEAT_PROD_XS] = (double)YES;

              /* Additional options */

              if (k < np)
                WDB[DATA_EDEP_CAPT_E] =
                  TestParam(pname, fname, line, params[k++], PTYPE_REAL,
                            0.0, 50.0);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "edepdel"))
            {
              /***** Delayed components for energy deposition ****************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if (k == np)
                Error(-1, params[j], fname, line,
                      "Invalid number of parameters for edepdel");

              WDB[DATA_EDEP_DELAYED] = TestParam(pname, fname, line,
                                                 params[k++], PTYPE_LOGICAL);

              if (k < np)
                WDB[DATA_EDEP_LOCAL_EGD] =
                  TestParam(pname, fname, line, params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "edepkcorr"))
            {
              /***** k-correction for energy deposition **********************/

              /* Copy parameter name */

              strcpy (pname, params[j]);

              k = j + 1;

              /* Check number of parameters */

              if (k == np)
                Error(-1, params[j], fname, line,
                      "Invalid number of parameters for edepkcorr");

              WDB[DATA_EDEP_KEFF_CORR] = TestParam(pname, fname, line,
                                                   params[k++], PTYPE_LOGICAL);

              /***************************************************************/
            }
          else if (!strcasecmp(params[j], "abs"))
            {
              /***** Unsupported Serpent 1 parameters ************************/

              /* These may seriously affect the results */

              Error(-1, params[j], fname, line,
                    "Parameter is not supported in Serpent 2");

              /***************************************************************/
            }
          else if ((!strcasecmp(params[j], "sym")) ||
                   (!strcasecmp(params[j], "remxs")) ||
                   (!strcasecmp(params[j], "stabcut")) ||
                   (!strcasecmp(params[j], "xsfcut")) ||
                   (!strcasecmp(params[j], "outfile")))

            {
              /***** Serpent 1 parameters ************************************/

              /* Less serious, print warning */

              Note(0,
                   "Serpent 1 input option %s on line %ld of file %s not used",
                   params[j], line, fname);

              /* Skip values */

              k = np;

              /***************************************************************/
            }
          else
            {
              /***** Parameter is not defined ********************************/

              Error(-1, params[j], fname, line, "Unidentified option");

              /* NOTE: Jos toi errori poistetaan niin kaikkia parametreja ei */
              /*       lueta --> aseta k = np jotta ajo ei kaadu alempana    */
              /*       olevaan testaukseen.  */
              /*
              k = np;
              */
              /***************************************************************/
            }

          /* Check number of parameters */

          if (k != np)
            Error(-1, params[j], fname, line, "Invalid number of parameters");
        }

      /* Free parameter list (can be null if empty input file is read) */

      if (params != NULL)
        {
          if (np > 0)
            for (n = 0; n < np + 1; n++)
              Mem(MEM_FREE, params[n]);

          Mem(MEM_FREE, params);
        }
    }

  /**************************************************************************/

  /* Free memory */

  Mem(MEM_FREE, input);
}

/*****************************************************************************/
