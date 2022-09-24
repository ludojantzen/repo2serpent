/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcfe.c                                    */
/*                                                                           */
/* Created:       2018/01/26 (BWe)                                           */
/* Last modified: 2018/11/14 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads FE-dependent multi-physics interfaces                  */
/*                                                                           */
/* Comments:    Based on readifcfunc.c                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCFE:"

/* Local definitions */

long ReadParameters(long, long, long, FILE *, char *, char *);
void SwapStrings(char **, char **);

/*****************************************************************************/

void ReadIFCFE(long loc0, long update)
{
  long bufferSize, inParams, outParams, loc1, val, nmat, i, ptr;
  long feType, interfaceType, makeOutput, copyInputSpecs = NO, n, npIn;
  double *coefficients, value;
  char tmpstr[MAX_STR];
  char *const fileName = GetText(loc0 + IFC_PTR_INPUT_FNAME);
  FILE *fp;

  if (update == YES)
    {
      /* Get the existing pointers */

      inParams = (long)RDB[loc0 + IFC_FET_INPUT_PARAMS_PTR];
      outParams = (long)RDB[loc0 + IFC_FET_OUTPUT_PARAMS_PTR];
    }
  else
    {
      /* Allocate memory */

      inParams = ReallocMem(DATA_ARRAY, FET_PARAMS_SIZE);
      outParams = ReallocMem(DATA_ARRAY, FET_PARAMS_SIZE);

      /* Store pointers */

      WDB[loc0 + IFC_FET_INPUT_PARAMS_PTR] = (double)inParams;
      WDB[loc0 + IFC_FET_OUTPUT_PARAMS_PTR] = (double)outParams;
    }

#ifndef DEBUG
  /* The filename and interface type are implicitly available, otherwise how  */
  /* would we get here from `ReadInterface()`?                                */

  fp = fopen(fileName, "r");
  val = fscanf(fp, "%ld", &interfaceType);
#else

  /* Avoid compiler warning */

  val = 0;

  /* Check anyways if we are in debug mode */

  if ((fp = fopen(fileName, "r")) == NULL)
    Error(loc0, "Multi-physics interface file \"%s\" does not exist", fileName);

  /* Read interface type */

  if (fscanf(fp, "%ld", &interfaceType) == EOF)
    Error(loc0, "No interface type provided in \"%s\"", fileName);
#endif /* DEBUG */

  /* Read FE type */

  if (fscanf(fp, "%ld", &feType) == EOF)
    Error(loc0, "Missing FET type after interface type \"%ld\" in \"%s\"",
          interfaceType, fileName);

  /* Using only one material */

  nmat = 1;

  if (!update)
    {
      /* Allocate memory for material array */

      ptr = ReallocMem(DATA_ARRAY, nmat + 1);

      /* Store pointer to material array */

      WDB[loc0 + IFC_PTR_MAT_ARR] = (double)ptr;
    }
  else
    {
      /* Get material array pointer from memory */

      ptr = (long)RDB[loc0 + IFC_PTR_MAT_ARR];
    }

  /* Read material names */

  for (i = 0; i < nmat; i++)
    {

      /* Read interface material */

      if (fscanf(fp, "%s", tmpstr) == EOF)
        Error(loc0, "Missing interface material after FE type \"%ld\" in \"%s\"",
              feType, fileName);

      /* Store material name */

      if (update == NO)
        WDB[ptr + i] = (double)PutText(tmpstr);
    }

  /* Null terminate material name array */

  if (update == NO)
    WDB[ptr + nmat] = (double)NULLPTR;

  /* Read output flag */

  if (fscanf(fp, "%ld", &makeOutput) == EOF)
    Error(loc0, "Missing output flag after interface material \"%s\" in \"%s\"",
          tmpstr, fileName);

  WDB[loc0 + IFC_CALC_OUTPUT] = (double)makeOutput;

  /* Read output file name and axial binning for output */

  if (makeOutput == YES)
    {
      /* Read name of the outermost material in the pin being scored            */
      /* This follows a thought made by VVa in the Serpent forum:               */
      /* http://ttuki.vtt.fi/serpent/viewtopic.php?f=24&t=2086&sid=2cc96988a4e8ed2c5940a1c6ea2de47c#p5347 */

      /* Read material name */

      if (fscanf(fp, "%s", tmpstr) == EOF)
        Die(FUNCTION_NAME, "Could not read material name from %s",
            GetText(loc0 + IFC_PTR_INPUT_FNAME));

      /* Store material name */

      if (update == NO)
        WDB[loc0 + IFC_PTR_OUTPUT_MAT] = (double)PutText(tmpstr);

      /* Read output file name */

      if (fscanf(fp, "%s", tmpstr) == EOF)
        Error(loc0,
              "Missing output file name after FET power scoring material \"%s\" in \"%s\"",
              tmpstr, fileName);

      WDB[loc0 + IFC_PTR_OUTPUT_FNAME] = (double)PutText(tmpstr);

      /* Loop over previous and check for duplicate file names */

      loc1 = PrevItem(loc0);
      while (loc1 > VALID_PTR)
        {
          /* Compare file names */

          if ((long)RDB[loc1 + IFC_PTR_OUTPUT_FNAME] > VALID_PTR)
            if (CompareStr(loc0 + IFC_PTR_OUTPUT_FNAME,
                           loc1 + IFC_PTR_OUTPUT_FNAME))
              Error(loc0,
                    "Output file name \"%s\" already used by \"%s\"",
                    GetText(loc0 + IFC_PTR_OUTPUT_FNAME),
                    GetText(loc1 + IFC_PTR_INPUT_FNAME));

          /* Pointer to previous */

          loc1 = PrevItem(loc1);
        }

      /* Should the input FE specifications be copied for the output specs? */

      if (fscanf(fp, "%ld", &copyInputSpecs) == EOF)
        Error(loc0,
              "Missing input-options-copy flag after output file name \"%s\" in \"%s\"",
              tmpstr, fileName);

      if (copyInputSpecs == NO)
        ReadParameters(loc0, feType, outParams, fp, fileName, "IFC_output_FET_parameters");
    }


  /*******************************************************************/

  /***** Functional expansion coefficients input *********************/

  /* Get parameters */

  npIn = ReadParameters(loc0, feType, inParams, fp, fileName, "IFC_input_FET_parameters");

  /* Get the input coefficients array */

  if ((update == NO) ||
      ((update == YES) && (npIn > (long)RDB[loc0 + IFC_FET_INPUT_COEFS_SIZE])))
    {
#ifdef DEBUG
      if ((update == YES) && ((long)RDB[loc0 + IFC_FET_INPUT_COEFS_SIZE] > 0))
        Warn(FUNCTION_NAME, "Resizing the FET input coefficients array; the previous array is being discarded but no deallocated");
#endif /* DEBUG */
      bufferSize = ceil(npIn * 1.5);
      loc1 = ReallocMem(DATA_ARRAY, bufferSize * sizeof (double));
      WDB[loc0 + IFC_FET_INPUT_COEFS_SIZE] = (double)bufferSize;
      WDB[loc0 + IFC_FET_PTR_INPUT_COEFS_ARRAY] = (double)loc1;
    }
  else
    {
      loc1 = (long)RDB[loc0 + IFC_FET_PTR_INPUT_COEFS_ARRAY];
      memset(&WDB[loc1],
             0.0,
             (long)RDB[loc0 + IFC_FET_INPUT_COEFS_SIZE] * sizeof (double));
    }

  /* Get the array of coefficients directly */

  coefficients = &WDB[loc1];

  /* Copy everything that needs to be copied if the output parameters are     */
  /* specified to match the input parameters.                                 */

  if (copyInputSpecs == YES)
    memcpy(&WDB[outParams + FET_PARAMS_COPY_START],
           &RDB[inParams + FET_PARAMS_COPY_START],
           (FET_PARAMS_COPY_END - FET_PARAMS_COPY_START) * sizeof (double));

  /* Ensure that the cache arrays are sized properly for FET calculations. */

  AllocFETCache(inParams, FET_USE_COEF, -1, update);
  if (makeOutput == YES)
    AllocFETCache(outParams, FET_GENERATE_COEF, -1, update);

  for (n = 0; n < npIn; ++n)
    {
      /* Get value */

      if (fscanf(fp, "%lf", &value) == EOF)
        Error(loc0,
              "No enough FET coefficients provided. Only %ld of %ld were found.",
              n, npIn);

      /* Put the value */

      coefficients[n] = value;
    }

  /* Reset maximum density to avoid error in processinterface.c */
  /* (value is not used for anything) */

  WDB[loc0 + IFC_MAX_DENSITY] = 0.0;

  /*******************************************************************/

  /* Close file */

  fclose(fp);
}

/*****************************************************************************/

long ReadParameters(long loc0, long feType, long ptr, FILE *fp,
                    char *fileName, char *parameterName)
{
  long np;
  char stringsArray[2][MAX_STR], c0, c1;
  char *pStr, *vStr;
  double n;

  pStr = stringsArray[0];
  vStr = stringsArray[1];

  /* Avoid compiler warnings */

  np = 0;
  c0 = 'x';
  c1 = 'y';

  WDB[ptr + FET_PARAM_TYPE] = feType;

  switch (feType)
    {
    case FET_TYPE_CARTESIAN:
      {
        /* Read the orders and dimensions */

        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing minimum 'x' dimension in \"%s\"", fileName);
        WDB[ptr + FET_CART_MIN_X] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, -INFTY, INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing maximum 'x' dimension after \"%s\" in \"%s\"", pStr, fileName);
        WDB[ptr + FET_CART_MAX_X] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, RDB[ptr + FET_CART_MIN_X], INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing polynomial order in dimension 'x' after \"%s\" in \"%s\"", pStr, fileName);
        n = TestParam(parameterName, fileName, -1, vStr, PTYPE_INT, -1, 100);
        WDB[ptr + FET_CART_ORDER_X] = n < 0 ? FET_DEFAULT_ORDER : n;
        WDB[ptr + FET_CART_NCOEF_X] = RDB[ptr + FET_CART_ORDER_X] + 1;

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing minimum 'y' dimension after \"%s\" in \"%s\"", pStr, fileName);
        WDB[ptr + FET_CART_MIN_Y] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, -INFTY, INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing maximum 'y' dimension after \"%s\" in \"%s\"", pStr, fileName);
        WDB[ptr + FET_CART_MAX_Y] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, RDB[ptr + FET_CART_MIN_Y], INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing polynomial order in dimension 'y' after \"%s\" in \"%s\"", pStr, fileName);
        n = TestParam(parameterName, fileName, -1, vStr, PTYPE_INT, -1, 100);
        WDB[ptr + FET_CART_ORDER_Y] = n < 0 ? FET_DEFAULT_ORDER : n;
        WDB[ptr + FET_CART_NCOEF_Y] = RDB[ptr + FET_CART_ORDER_Y] + 1;

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing minimum 'z' dimension after \"%s\" in \"%s\"", pStr, fileName);
        WDB[ptr + FET_CART_MIN_Z] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, -INFTY, INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing maximum 'z' dimension after \"%s\" in \"%s\"", pStr, fileName);
        WDB[ptr + FET_CART_MAX_Z] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, RDB[ptr + FET_CART_MIN_Z], INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing polynomial order in dimension 'z' after \"%s\" in \"%s\"", pStr, fileName);
        n = TestParam(parameterName, fileName, -1, vStr, PTYPE_INT, -1, 100);
        WDB[ptr + FET_CART_ORDER_Z] = n < 0 ? FET_DEFAULT_ORDER : n;
        WDB[ptr + FET_CART_NCOEF_Z] = RDB[ptr + FET_CART_ORDER_Z] + 1;

        /* Calculate the number of coefficients */

        np = RDB[ptr + FET_CART_NCOEF_X]
          * RDB[ptr + FET_CART_NCOEF_Y]
          * RDB[ptr + FET_CART_NCOEF_Z];

        break;
      }

    case FET_TYPE_CYLINDRICAL:
      {
        /* Read the orders and dimensions */

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing radius 'r' in \"%s\"", fileName);
        WDB[ptr + FET_CYL_MAX_R] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, -INFTY, INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing Zernike polynomial order after \"%s\" in \"%s\"", pStr, fileName);
        n = TestParam(parameterName, fileName, -1, vStr, PTYPE_INT, -1, 100);
        np = (long)(n < 0 ? FET_DEFAULT_ORDER : n);
        WDB[ptr + FET_CYL_ORDER_R] = np;
        WDB[ptr + FET_CYL_NCOEF_R] = (double)((np + 1) * (np + 2) / 2);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing axial start position after \"%s\" in \"%s\"", pStr, fileName);
        WDB[ptr + FET_CYL_MIN_A] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, -INFTY, INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing axial end position after \"%s\" in \"%s\"", pStr, fileName);
        WDB[ptr + FET_CYL_MAX_A] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, RDB[ptr + FET_CYL_MIN_A], INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing axial polynomial order after \"%s\" in \"%s\"", pStr, fileName);
        n = TestParam(parameterName, fileName, -1, vStr, PTYPE_INT, -1, 100);
        WDB[ptr + FET_CYL_ORDER_A] = n < 0 ? FET_DEFAULT_ORDER : n;
        WDB[ptr + FET_CYL_NCOEF_A] = RDB[ptr + FET_CYL_ORDER_A] + 1;

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing axial orientation after \"%s\" in \"%s\"", pStr, fileName);
        WDB[ptr + FET_CYL_ORIENTATION_A] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_INT, FET_ORIENTATIONS_START, FET_ORIENTATIONS_END);

        switch ((long)RDB[ptr + FET_CYL_ORIENTATION_A])
          {
          case FET_ORIENTATION_X:
            {
              c0 = 'y';
              c1 = 'z';

              break;
            }

          case FET_ORIENTATION_Y:
            {
              c0 = 'x';
              c1 = 'z';

              break;
            }

          case FET_ORIENTATION_Z:
            {
              c0 = 'x';
              c1 = 'y';

              break;
            }

#ifdef DEBUG /* We should never get here because the TestParam() above should have caught this */
          default:
            {
              Die(FUNCTION_NAME, "Unsupported orientation \"%ld\"  in \"%s\"",
                  (long)RDB[ptr + FET_CYL_ORIENTATION_A], fileName);
            }
#endif /* DEBUG */
          }

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing \"%c\" coordinate of cylinder center after \"%s\" in \"%s\"", c0, pStr, fileName);
        WDB[ptr + FET_CYL_CENTER_A0] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, -INFTY, INFTY);

        SwapStrings(&pStr, &vStr);
        if (fscanf(fp, "%s", vStr) == EOF)
          Error(loc0, "Missing \"%c\" coordinate of cylinder center after \"%s\" in \"%s\"", c1, pStr, fileName);
        WDB[ptr + FET_CYL_CENTER_A1] =
          TestParam(parameterName, fileName, -1, vStr,
                    PTYPE_REAL, -INFTY, INFTY);

        /* Calculate the number of coefficients */

        np = RDB[ptr + FET_CYL_NCOEF_R] * RDB[ptr + FET_CYL_NCOEF_A];

        break;
      }

    default:
      {
        Die(FUNCTION_NAME, "Unsupported FE type \"%ld\" in \"%s\"", feType, fileName);
      }
    }

  return np;
}

void SwapStrings(char **string1, char **string2)
{
  char *temp = *string1;
  *string1 = *string2;
  *string2 = temp;
}
