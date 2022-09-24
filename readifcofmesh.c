/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcofmesh.c                                */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2018/11/08 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads OpenFOAM multi-physics interfaces                      */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*             - Tiheysjakauman yksiköt luetaan nyt dimensions -vektorista   */
/*               (2.12.2014 / 2.1.23 / JLe)                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCOFMesh:"

/*****************************************************************************/

void ReadIFCOFMesh(long ifc, long update)
{
  long loc1, ptr, type, nc, n, nmat;
  long i, j;
  long dflist, tmplist;
  double d0, T0;
  char tmpstr[MAX_STR], pfile[MAX_STR], ffile[MAX_STR], ofile[MAX_STR];
  char nfile[MAX_STR], rfile[MAX_STR], tfile[MAX_STR], mapfile[MAX_STR];
  char matfile[MAX_STR];
  FILE *fp;

  /* Open file for reading */

  if ((fp = fopen(GetText(ifc + IFC_PTR_INPUT_FNAME), "r")) == NULL)
    Error(ifc, "Multi-physics interface file \"%s\" does not exist",
          GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Read interface type */

  if (fscanf(fp, "%ld", &type) == EOF)
    Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

  /* Read material name or uni & BG uni */

  if(type != IFC_TYPE_OF_SOLID)
    {
      if (fscanf(fp, "%s", tmpstr) == EOF)
        Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

      /* Allocate memory for the material array and store material name */

      if (update == NO)
        {
          /* IFC files always contain only one material name */

          nmat = 1;

          /* Store number of materials */

          WDB[ifc + IFC_N_MAT] = (double)nmat;

          /* Allocate memory for material array */

          ptr = ReallocMem(DATA_ARRAY, nmat + 1);

          /* Store pointer to material array */

          WDB[ifc + IFC_PTR_MAT_ARR] = (double)ptr;

          /* Store material name */

          WDB[ptr++] = (double)PutText(tmpstr);

          /* Null terminate the array */

          WDB[ptr] = NULLPTR;
        }
    }
  else
    {
      /* Read universe name (stored in readinput) */

      if (fscanf(fp, "%s", tmpstr) == EOF)
        Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

      /* Read BG universe name (stored in readinput) */

      if (fscanf(fp, "%s", tmpstr) == EOF)
        Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

    }

  /* Read output flag */

  if (fscanf(fp, "%ld", &n) == EOF)
    Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

  /* Check output flag value */

  CheckValue(FUNCTION_NAME, "n", "Output flag", n, -IFC_PRINT_ERROR_BOTH, IFC_PRINT_ERROR_BOTH);

  /* Read output file name for output */

  if (n != NO)
    {
      /* Set flag */

      WDB[ifc + IFC_CALC_OUTPUT] = (double)YES;

      /* Check negative flag (for scoring energy deposited in all materials )*/

      if (n < 0)
        {
          /* Set flag */

          WDB[ifc + IFC_SCORE_ALL_MATERIALS] = (double)YES;

          /* Remove minus sign */

          n = -1;
        }

      /* Set error printing flag */

      if (n != YES)
        WDB[ifc + IFC_PRINT_ERROR] = (double)n;

      /* Read only file name */

      if (fscanf(fp, "%s", tmpstr) == EOF)
        Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

      WDB[ifc + IFC_PTR_OUTPUT_FNAME] = (double)PutText(tmpstr);

      /* Loop over previous and check duplicate file names */

      loc1 = PrevItem(ifc);
      while (loc1 > VALID_PTR)
        {
          /* Compare file names */

          if ((long)RDB[loc1 + IFC_PTR_OUTPUT_FNAME] > VALID_PTR)
            if (CompareStr(ifc + IFC_PTR_OUTPUT_FNAME,
                           loc1 + IFC_PTR_OUTPUT_FNAME))
              Error(ifc,
                    "Duplicate output file name with distribution \"%s\"",
                    GetText(loc1 + IFC_PTR_INPUT_FNAME));

          /* Pointer to previous */

          loc1 = PrevItem(loc1);
        }
    }
  else
    {
      /* Reset flag */

      WDB[ifc + IFC_CALC_OUTPUT] = (double)NO;
    }

  /*******************************************************************/

  /***** OpenFoam mesh ***********************************************/

  /* Read nominal temperature and density */

  if (fscanf(fp, "%lf %lf", &d0, &T0) == EOF)
    Error(ifc, "Missing nominal temperature or density");

  /* Store nominal density and temperature */

  WDB[ifc + IFC_NOMINAL_DENSITY] = d0;
  WDB[ifc + IFC_NOMINAL_TEMPERATURE] = T0;

  /* Read search mesh size */

  if (fscanf(fp, "%ld", &i) == EOF)
    Die(FUNCTION_NAME, "fscanf error (8)");
  else if (i < 2)
    Error(ifc, "Error in search mesh split criterion: Value given was %ld, must be >= 2", i);
  else
    {
      if(update == NO)
        WDB[ifc + IFC_SEARCH_MESH_ADA_SPLIT] = (double)i;
    }

  if (fscanf(fp, "%ld", &i) == EOF)
    Die(FUNCTION_NAME, "fscanf error (9)");
  else if ((i < 1) || (i > 10))
    Error(ifc, "Error in search mesh depth");

  /* Allocate memory for search mesh size */

  if (update == NO)
    {
      ptr = ReallocMem(DATA_ARRAY, i + 1);
      WDB[ifc + IFC_SEARCH_MESH_ADA_PTR_SZ] = (double)ptr;
    }
  else
    ptr = -1;

  /* Read sizes */

  for (n = 0; n < i; n++)
    {
      if (fscanf(fp, "%ld", &j) == EOF)
        Die(FUNCTION_NAME, "fscanf error (10)");
      else if (j < 1)
        Error(ifc, "Error in search mesh size");
      else
        {
          if(update == NO)
            WDB[ptr++] = (double)j;
        }
    }

  /* Put terminator */

  if(update == NO)
    WDB[ptr] = -1.0;

  /* Read file names (points, faces, owner and neighbour) */

  if (fscanf(fp, "%s", pfile) == EOF)
    Error(ifc, "Missing path to points file");
  else
    TestDOSFile(pfile);

  WDB[ifc + IFC_PTR_OF_PFILE] = (double)PutText(pfile);

  if (fscanf(fp, "%s", ffile) == EOF)
    Error(ifc, "Missing path to faces file");
  else
    TestDOSFile(ffile);

  WDB[ifc + IFC_PTR_OF_FFILE] = (double)PutText(ffile);

  if (fscanf(fp, "%s", ofile) == EOF)
    Error(ifc, "Missing path to owner file");
  else
    TestDOSFile(ofile);

  WDB[ifc + IFC_PTR_OF_OFILE] = (double)PutText(ofile);

  if (fscanf(fp, "%s", nfile) == EOF)
    Error(ifc, "Missing path to neighbour file");
  else
    TestDOSFile(nfile);

  WDB[ifc + IFC_PTR_OF_NFILE] = (double)PutText(nfile);

  if ((type == IFC_TYPE_OF_MAT) || (type == IFC_TYPE_OF_SOLID))
    {

      /* Read material file */

      if (fscanf(fp, "%s", matfile) == EOF)
        Error(ifc, "Missing path to material file");
      else
        TestDOSFile(matfile);

      WDB[ifc + IFC_PTR_OF_MFILE] = (double)PutText(matfile);

    }

  /* Read density file */

  if (fscanf(fp, "%s", rfile) == EOF)
    Error(ifc, "Missing path to density file");
  else
    TestDOSFile(rfile);

  WDB[ifc + IFC_PTR_OF_RFILE] = (double)PutText(rfile);

  /* Read type flag */

  if (strcmp(rfile, "-1"))
    {
      if (fscanf(fp, "%ld", &n) == EOF)
        Error(ifc, "Missing density distribution type flag");
      else if (n != 1)
        Error(ifc, "Invalid density distribution type flag");

      /* Store type flag */

      WDB[ifc + IFC_DENSITY_TYPE] = n;
    }

  /* Read temperature file */

  if (fscanf(fp, "%s", tfile) == EOF)
    Error(ifc, "Missing path to temperature file");
  else
    TestDOSFile(tfile);

  WDB[ifc + IFC_PTR_OF_TFILE] = (double)PutText(tfile);

  /* Read type flag */

  if (strcmp(tfile, "-1"))
    {
      if (fscanf(fp, "%ld", &n) == EOF)
        Error(ifc, "Missing temperature distribution type flag");
      else if ((n != 1) && (n != 2))
        Error(ifc, "Invalid temperature distribution type flag");

      /* Store type flag */

      WDB[ifc + IFC_TEMPERATURE_TYPE] = n;
    }

  /* Read output mapping file */

  if ((long)RDB[ifc + IFC_CALC_OUTPUT] == YES)
    {
      if (fscanf(fp, "%s", mapfile) == EOF)
        Error(ifc, "Missing path to output map file");
      else
        TestDOSFile(mapfile);

      WDB[ifc + IFC_PTR_OF_MAPFILE] = (double)PutText(mapfile);
    }
  else
    WDB[ifc + IFC_PTR_OF_MAPFILE] = (double)PutText("-1");

  /* Close file */

  fclose(fp);

  /*******************************************************************/

  /***** Read basic mesh-data (no fields) ****************************/

  if (update == NO)
    ReadOFMesh(ifc);

  /*******************************************************************/

  /***** Read density distribution ***********************************/

  /* Check if density distribution is given */

  if (strcmp(rfile, "-1"))
    ReadOFDensities(ifc, update);
  else
    {
      /* Put nominal values */

      WDB[ifc + IFC_MAX_DENSITY] = d0;
      WDB[ifc + IFC_MIN_DENSITY] = d0;

      /* Create density factor list with nominal values */
      /* this way the density factor list always exists */

      if (update == NO)
        {
          /* Get number of cells */

          nc = (long)RDB[ifc + IFC_NC_PARENTS];

          /* Allocate memory for list */

          dflist = ReallocMem(DATA_ARRAY, nc);

          /* Store pointer to list */

          WDB[ifc + IFC_PTR_DF_LIST] = (double)dflist;

          /* Set cell densities to nominal value */

          for (n = 0; n < nc; n++)
            {
              WDB[dflist + n] = d0;
            }
        }
    }

  /*******************************************************************/

  /***** Read temperature distribution *******************************/

  /* Check if temperature distribution is given */

  if (strcmp(tfile, "-1"))
    ReadOFTemperatures(ifc, update);
  else
    {
      /* Put nominal values */

      WDB[ifc + IFC_MAX_TEMP] = T0;
      WDB[ifc + IFC_MIN_TEMP] = T0;

      /* Create temperature list with nominal values */
      /* this way the temperature list always exists */

      if (update == NO)
        {
          /* Get number of cells */

          nc = (long)RDB[ifc + IFC_NC_PARENTS];

          /* Allocate memory for list */

          tmplist = ReallocMem(DATA_ARRAY, nc);

          /* Store pointer to list */

          WDB[ifc + IFC_PTR_TMP_LIST] = (double)tmplist;

          /* Set cell temperatures */

          for (n = 0; n < nc; n++)
            {
              WDB[tmplist + n] = T0;
            }
        }
    }

  /* Switch type to a regular tet mesh */

  WDB[ifc + IFC_TYPE] = (double)IFC_TYPE_TET_MESH;

  /*******************************************************************/

  fprintf(outp, "\n");

  /*******************************************************************/

}

/*****************************************************************************/
