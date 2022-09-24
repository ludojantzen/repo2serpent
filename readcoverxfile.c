/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readcoverxfile.c                               */
/*                                                                           */
/* Created:       2010/06/06 (VVa)                                           */
/* Last modified: 2018/09/17 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads data from COVERX format multigroup covariance library  */
/*              file into COV block                                          */
/*                                                                           */
/* Comments: - Currently reads all nuclides and all reactions.               */
/*           - The COVERX format is described in:                            */
/*             B. T. Rearden and M. A. Jessee, Eds., "SCALE Code System",    */
/*             ORNL/TM-2005/39, Version 6.2.1, Oak Ridge National            */
/*             Laboratory, Oak Ridge, Tennessee (2016).                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadCOVERXFile:"

/*
#define COVERX_DEBUG_PRINT
*/
long COVERXreadBOF(FILE *, FILE *, long, char *, char *, long *, long *);
long COVERXreadBOFASCII(FILE *, FILE *, long, char *, char *, long *, long *);
long COVERXreadFC(FILE *, FILE *, long [7], long);
long COVERXreadFCASCII(FILE *, FILE *, long [7]);
long COVERXreadFD(FILE *, FILE *, long, long, char *, long);
long COVERXreadFDASCII(FILE *, FILE *, long, long, char *);
long COVERXreadEGrid(FILE *, FILE *, long, double *, long);
long COVERXreadEGridASCII(FILE *, FILE *, long, double *);
long COVERXreadMMTControl(FILE *, FILE *, long, long);
long COVERXreadMMTControlASCII(FILE *, FILE *, long);
long COVERXreadMMTData(FILE *, FILE *, long, long, long);
long COVERXreadMMTDataASCII(FILE *, FILE *, long, long);
long COVERXreadMatrixControl(FILE *, FILE *, long [5], long);
long COVERXreadMatrixControlASCII(FILE *, FILE *, long [5]);
long COVERXreadMBC(FILE *, FILE *, long, long, long *, long *, long);
long COVERXreadMBCASCII(FILE *, FILE *, long, long, long *, long *);
long COVERXreadMatrix(FILE *, FILE *, long, long *,
                      long *, double *, long);
long COVERXreadMatrixASCII(FILE *, FILE *, long, long *,
                           long *, double *);
long COVERXskipMatrix(FILE *, FILE *, long, long, long);
long COVERXskipMatrixASCII(FILE *, FILE *, long, long);
long SCALEIDtoZAI(long);

/*****************************************************************************/

void ReadCOVERXFile(long nameptr)
{
  char tmpstr[MAX_STR], tmpstr2[MAX_STR], fname[MAX_STR];
  double *covdata, *grid;
  long neutrongrid, gammagrid, loc0, ptr, nuc, included, asciiflag;
  long nholl, ngroup, nngrup, nggrup, mult, nmmp;
  long ntype, nmtrix, mat1, mat2, mt1, mt2, nblock, *bandwidth, *diagpos;
  long i, j, k, j2, k2, n1, n2, tmplong, filecontrol[7], ZAI1, ZAI2;
  long matrixcontrol[5];
  FILE *fp, *dbgp;

#ifdef COVERX_DEBUG_PRINT
  dbgp = outp;
#else
  dbgp = fopen("/dev/null", "w");
#endif

  /* Test format */
  /*
  TestDOSFile(GetText(nameptr));
  */
  /* Open file for reading */

  fp = OpenDataFile(nameptr, "COVERX data file");

  /* Print name to string */

  sprintf(fname, "%s", GetText(nameptr));

  /* Print out note to terminal */

  fprintf(outp, "Reading COVERX-format file \"%s\":\n", fname);

  /* This seems to be the correct value */

  mult = 6;

  /* Initialize some values */

  neutrongrid = 0;
  gammagrid = 0;
  asciiflag = NO;

  /******************************/
  /* Read beginning of the file */
  /******************************/

  if (COVERXreadBOF(fp, dbgp, mult, tmpstr, tmpstr2, &tmplong, &asciiflag))
    Die(FUNCTION_NAME, "Error in beginning of file of COVERX file \"%s\"", fname);

  fprintf(outp, " - BOF data: \"%s\" \"%s\" version %ld\n", tmpstr, tmpstr2, tmplong);

  /*************************/
  /* Read the file control */
  /*************************/

  if (COVERXreadFC(fp, dbgp, filecontrol, asciiflag))
    Die(FUNCTION_NAME, "Error in file control of COVERX file \"%s\"", fname);

  /* Store the seven values to well named variables */

  ngroup = filecontrol[0];
  nngrup = filecontrol[1];
  nggrup = filecontrol[2];
  ntype  = filecontrol[3];
  nmmp   = filecontrol[4];
  nmtrix = filecontrol[5];
  nholl  = filecontrol[6];

  fprintf(outp, " - %-3ld energy groups total (%-3ld neutron, %-3ld gamma)\n",
          ngroup, nngrup, nggrup);

  /* Check type */

  if ((ntype < 1) || (ntype > 3))
    Die(FUNCTION_NAME, "Invalid data type %ld", ntype);
  else if (ntype == 1)
    {
      fprintf(outp, " - Data type is covariance matrix (standard deviation)\n");
      Die(FUNCTION_NAME, "Absolute covariance matrices are not yet supported.\n");
    }
  else if (ntype == 2)
    {
      fprintf(outp, " - Data type is relative covariance matrix (relative standard deviation)\n");
    }
  else if (ntype == 3)
    {
      fprintf(outp, " - Data type is correlation matrix (standard deviation)\n");
      Die(FUNCTION_NAME, "Correlation matrices are not yet supported.\n");
    }
  /* Print out number of data */

  fprintf(outp, " - File contains data for %ld material-reaction pairs in %ld matrices\n", nmmp, nmtrix);

  /* Check whether the file contains both neutron and photon groups */
  /* (reading hasn't been implemented yet) */

  if (nngrup*nggrup != 0)
    Die(FUNCTION_NAME, "COVERX format file \"%s\" contains both neutron and "
        "gamma groups, reading of such files has not been implemented yet.", fname);

  /*****************************/
  /* Read the file description */
  /*****************************/

  if (COVERXreadFD(fp, dbgp, mult, nholl, tmpstr, asciiflag))
    Die(FUNCTION_NAME, "Error in the file description of COVERX file \"%s\"", fname);

  fprintf(outp, " - File description: %s\n", tmpstr);

  /*************************************/
  /* Read the neutron group boundaries */
  /*************************************/

  if (nngrup > 0)
    {
      /* Allocate memory for the grid */

      grid = Mem(MEM_ALLOC, nngrup+1, sizeof(double));

      /* Read the energy grid */

      if (COVERXreadEGrid(fp, dbgp, nngrup, grid, asciiflag))
        Die(FUNCTION_NAME, "Error in the neutron group boundaries of COVERX file \"%s\"", fname);

      /* Print out information about the group structure */

      fprintf(outp, " - Neutron group structure between %.5E and %.5E MeV\n",
              grid[0], grid[nngrup]);

      /* Create the energy grid for the neutron groups */
      /* Keep track of pointer for now */

      neutrongrid = MakeEnergyGrid(nngrup+1, 0,0,-1, grid, EG_INTERP_MODE_LIN);

      /* Free the temporary grid allocated by COVERXreadEGrid */

      Mem(MEM_FREE, grid);

    }

  /***********************************/
  /* Read the gamma group boundaries */
  /***********************************/

  if (nggrup > 0)
    {
      /* Allocate memory for the grid */

      grid = Mem(MEM_ALLOC, nggrup+1, sizeof(double));

      /* Read the energy grid */

      if (COVERXreadEGrid(fp, dbgp, nggrup, grid, asciiflag))
        Die(FUNCTION_NAME, "Error in the gamma group boundaries of COVERX file \"%s\"", fname);

      /* Print out information about the group structure */

      fprintf(outp, " - Neutron group structure between %.5E and %.5E MeV\n",
              grid[0], grid[nggrup]);

      /* Create the energy grid for the neutron groups */
      /* Keep track of pointer for now */

      gammagrid = MakeEnergyGrid(nggrup+1, 0,0,-1, grid, EG_INTERP_MODE_LIN);

      /* Free the temporary grid allocated by COVERXreadEGrid */

      Mem(MEM_FREE, grid);
    }

  /**************************************/
  /* Read the material-reaction control */
  /* currently just skips the data      */
  /**************************************/

  if (COVERXreadMMTControl(fp, dbgp, nmmp, asciiflag))
    Die(FUNCTION_NAME, "Error in the material-reaction control of COVERX file \"%s\"", fname);

  /********************************************************************/
  /* Loop over all material-reaction pairs to read XS and error files */
  /********************************************************************/

  if (COVERXreadMMTData(fp, dbgp, nmmp, ngroup, asciiflag))
    Die(FUNCTION_NAME, "Error in the material-reaction XS and error data of COVERX file \"%s\"", fname);

  /*********************************************/
  /* Loop over all matrices to read everything */
  /*********************************************/

  fprintf(dbgp, "\n *** Reading all matrices:\n");

  /* Reset number of included data */

  included = 0;

  /* Allocate memory for working arrays */

  bandwidth = Mem(MEM_ALLOC, ngroup, sizeof(long));
  diagpos = Mem(MEM_ALLOC, ngroup, sizeof(long));

  for (i = 0; i < nmtrix; i++)
    {
      /***************************/
      /* Read the matrix control */
      /***************************/

      if (COVERXreadMatrixControl(fp, dbgp, matrixcontrol, asciiflag))
        Die(FUNCTION_NAME, "Error in matrix control for matrix %ld/%ld "
            "of COVERX file \"%s\"", i+1, nmtrix, fname);


      /* Store the seven values to well named variables */

      mat1   = matrixcontrol[0];
      mt1    = matrixcontrol[1];
      mat2   = matrixcontrol[2];
      mt2    = matrixcontrol[3];
      nblock = matrixcontrol[4];

      /* Check the number of blocks, I'm not sure how the data will be */
      /* organized in the case of multiple blocks */

      if (nblock != 1)
        Die(FUNCTION_NAME, "Number of blocks in the matrix is not 1 (it is %ld)"
                ", reading of such matrices has not yet been implemented.");

      /***************************************************/
      /* Figure out whether we want to store this matrix */
      /***************************************************/

      /* Store it only if both mat-mt pairs are included in the simulation */

      n1 = 0;
      n2 = 0;

      /* Convert SCALE IDs to ZAI */

      ZAI1 = SCALEIDtoZAI(mat1);
      ZAI2 = SCALEIDtoZAI(mat2);

      /* Skip nuclides that couldn't be translated to ZAI*/

      if (ZAI1*ZAI2 == 0)
        Die(FUNCTION_NAME, "This should not happen");

      /* Others should be simply ZA with I=0*/
      /* Find first material */

      nuc = (long)RDB[DATA_PTR_NUC0];

      while (nuc > VALID_PTR)
        {
          /* Check if this nuclide corresponds to mat1 */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == ZAI1)
            n1 = 1;

          /* Check if this nuclide corresponds to mat2 */

          if ((long)RDB[nuc + NUCLIDE_ZAI] == ZAI2)
            n2 = 1;

          /* Check if already found both */

          if (n1*n2 == 1)
            break;

          /* Compare next nuclide */

          nuc = NextItem(nuc);
        }

      /* Skip nuclides that are not part of the calculation */

      if (n1*n2 == 0)
        {
          /*
          fprintf(outp, "Skipping (%ld/mt %ld, %ld/mt %ld) (not present)\n", ZAI1, mt1,
              ZAI2, mt2);
          */
          if (COVERXskipMatrix(fp, dbgp, ngroup, nblock, asciiflag))
            Die(FUNCTION_NAME, "Error in matrix data for matrix %ld/%ld "
                "of COVERX file \"%s\"", i+1, nmtrix, fname);

          continue;
        }

      /* Increment the number of included material-reaction pair pairs */

      included++;

      /* Create the matrix */

      loc0 = NewItem(DATA_PTR_COVMTX0, COVMTX_BLOCK_SIZE);

      /* Store information about the nuclides and MTs */

      WDB[loc0 + COVMTX_ZAI1] = (double)ZAI1;
      WDB[loc0 + COVMTX_ZAI2] = (double)ZAI2;
      WDB[loc0 + COVMTX_MT1] = (double)mt1;
      WDB[loc0 + COVMTX_MT2] = (double)mt2;

      /* Store number of energy groups */

      WDB[loc0 + COVMTX_NG] = ngroup;

      /**************************************/
      /* Read block control for this matrix */
      /**************************************/

      if (COVERXreadMBC(fp, dbgp, ngroup, nblock, bandwidth, diagpos, asciiflag))
        Die(FUNCTION_NAME, "Error in block control for matrix %ld/%ld "
            "of COVERX file \"%s\"", i+1, nmtrix, fname);


      fprintf(dbgp, "\nReading block control for matrix:\n");

      /************************************/
      /* Read matrix data for this matrix */
      /************************************/

      /* This needs to be zeroed out after each matrix */
      /* since not all elements will be set */

      covdata = Mem(MEM_ALLOC, ngroup*ngroup, sizeof(double));

      /* Read data to the temporary matrix */

      if (COVERXreadMatrix(fp, dbgp, ngroup, bandwidth, diagpos, covdata, asciiflag))
        Die(FUNCTION_NAME, "Error in matrix data for matrix %ld/%ld "
            "of COVERX file \"%s\"", i+1, nmtrix, fname);

      /* Allocate memory for the data in the covariance matrix */

      ptr = ReallocMem(DATA_ARRAY, ngroup*ngroup);

      /* Store pointer */

      WDB[loc0 + COVMTX_PTR_DATA] = ptr;

      /* Store pointer to energy grid */

      WDB[loc0 + COVMTX_PTR_EGRID] = neutrongrid;

      /* Store data in the reverse group order    */
      /* stored data goes from low energy to high */

      for (j = 0; j < ngroup; j++)
        {
          for (k = 0; k < ngroup; k++)
            {
              j2 = ngroup-1-j;
              k2 = ngroup-1-k;
              WDB[ptr + j*ngroup + k] = covdata[j2*ngroup + k2];
            }
        }

      /* Print out matrix shape for debugging */

      /*
        for (j = 0; j < ngroup; j++)
          {
            for (k = 0; k < ngroup; k++)
              {
                j2 = ngroup-1-j;
                k2 = ngroup-1-k;
                if (covdata[j2*ngroup + k2] == 0)
                  fprintf(outp, " ");
                else
                  fprintf(outp, "*");
              }
            fprintf(outp, "\n");
          }
      */
      Mem(MEM_FREE, covdata);
    }

  /* Free the working arrays */

  Mem(MEM_FREE, bandwidth);
  Mem(MEM_FREE, diagpos);

  fprintf(outp, " - Included covariance data for %ld (nuclide-mt, nuclide-mt) pairs\n",
          included);

  fprintf(dbgp, "\nAll read\n");
  fprintf(outp, "\n");

  /* Check that we cannot read anything else from the file */

  if (asciiflag == YES)
    {
      tmpstr[0] = '\0';

      if (fscanf(fp, "%s", tmpstr) != 0)
        {
          if (tmpstr[0] != '\0')
            Warn(FUNCTION_NAME, "Additional data at the end of a COVERX-format file or "
                 "something wrong with the reading routine.");
        }

      if (feof(fp) == 0)
        Warn(FUNCTION_NAME, "Additional data at the end of a COVERX-format file or "
             "something wrong with the reading routine.");
    }
  else if (fread(tmpstr, sizeof(char), 1, fp) != 0)
    Warn(FUNCTION_NAME, "Additional data at the end of a COVERX-format file or "
         "something wrong with the reading routine.");

  fclose(fp);

#ifdef COVERX_DEBUG_PRINT

#else
  fclose(dbgp);
#endif

  return;
}

/*****************************************************************************/

long COVERXreadBOF(FILE *fp, FILE *dbgp, long mult, char *hname, char *huse, long *ivers, long *asciiflag)
{
  int dum;
  long expect;
  char tmpstr[3];

  /* Expected size of block */

  expect = sizeof(char)*6*3+sizeof(int);

  /* Check ASCII */

  if (fread(&tmpstr, sizeof(char), 3, fp) != 3)
    return 1;

  /* Go back to beginning */

  rewind(fp);

  if ((tmpstr[0] == ' ')  && (tmpstr[1] == '0')  && (tmpstr[2] == 'v'))
    {
      /* File in ASCII format */

      *asciiflag = YES;

      return COVERXreadBOFASCII(fp, dbgp, mult, hname, huse, ivers, asciiflag);
    }
  else
    {
      /* File not in ASCII format */

      *asciiflag = NO;
    }

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadBOF", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Read name (coverx) */

  if (fread(hname, sizeof(char), mult, fp) != mult)
    return 1;

  /* Null terminate */
  hname[6] = '\0';

  /* Read first part of user id */

  if (fread(huse, sizeof(char), mult, fp) != mult)
    return 1;

  /* Read second part of user id */

  if (fread(huse+mult, sizeof(char), mult, fp) != mult)
    return 1;

  /* Null terminate */
  huse[12] = '\0';

  /* Read file version number */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  *ivers = (long)dum;

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadBOF", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  return 0;
}

long COVERXreadBOFASCII(FILE *fp, FILE *dbgp, long mult, char *hname, char *huse, long *ivers, long *asciiflag)
{
  int dum;
  char tmpstr1[MAX_STR], tmpstr2[MAX_STR];

  /* Read " 0v coverx " */

  if (fread(hname, sizeof(char), 11, fp) != 11)
    return 1;

  /* Terminate string */

  hname[11] = '\0';

  /* Read use (the length of this may depend on library?) */

  if (fscanf(fp, "%s %s", tmpstr1, tmpstr2) != 2)
    return 1;

  /* Catenate the use */

  sprintf(huse, "%s %s", tmpstr1, tmpstr2);

  /* Read file version number */

  if (fscanf(fp, "%d", &dum) != 1)
    fprintf(outp, "Third %s\n", huse);

  *ivers = (long)dum;

  return 0;
}

/*****************************************************************************/

long COVERXreadFC(FILE *fp, FILE *dbgp, long filecontrol[7], long asciiflag)
{
  int dum;
  long expect, i;

  if (asciiflag == YES)
    return COVERXreadFCASCII(fp, dbgp, filecontrol);

  /* Expected size of block */

  expect = 7*sizeof(int);

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadFC", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Read data to array */

  for (i = 0; i < 7; i++)
    {
      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;

      /* Cast to long and store */

      filecontrol[i] = (long)dum;
    }

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadFC", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  return 0;
}

/*****************************************************************************/

long COVERXreadFCASCII(FILE *fp, FILE *dbgp, long filecontrol[7])
{
  int dum;
  long i;
  char tmpstr[MAX_STR], expect[3] = "1d\0";

  /* Read block number */

  if (fscanf(fp, "%s", tmpstr) != 1)
    return 1;

  /* Check block number */

  if (strcasecmp(tmpstr, expect) != 0)
    Die("COVERXreadFC", "Unexpected upcoming block number %s, expected %s.",
        tmpstr, expect);

  /* Read data to array */

  for (i = 0; i < 7; i++)
    {
      if (fscanf(fp, "%d", &dum) != 1)
        return 1;

      /* Cast to long and store */

      filecontrol[i] = (long)dum;
    }

  return 0;
}

/*****************************************************************************/

long COVERXreadFD(FILE *fp, FILE *dbgp, long mult, long nholl, char *description, long asciiflag)
{
  int dum;
  long expect, i;

  if (asciiflag == YES)
    return COVERXreadFDASCII(fp, dbgp, mult, nholl, description);

  /* If nholl is zero, no description block is included */

  if (nholl == 0)
    {
      description[0] = '\0';
      return 0;
    }

  /* Expected size of block */

  expect = nholl*mult*sizeof(char);

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadFD", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Read data to description */

  for (i = 0; i < nholl; i++)
    {
      if (fread(description + i*mult, sizeof(char), mult, fp) != mult)
        return 1;

    }

  /* Null terminate the string */

  description[i*mult] = '\0';

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadFD", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  return 0;
}

/*****************************************************************************/

long COVERXreadFDASCII(FILE *fp, FILE *dbgp, long mult, long nholl, char *description)
{
  char tmpstr[MAX_STR], expect[3] = "2d\0";

  /* Read block number */

  if (fscanf(fp, "%s", tmpstr) != 1)
    return 1;

  /* Check block number */

  if (strcasecmp(tmpstr, expect) != 0)
    Die("COVERXreadFC", "Unexpected upcoming block number %s, expected %s.",
        tmpstr, expect);

  /* If nholl is zero, no description block is included */

  if (nholl == 0)
    {
      description[0] = '\0';
      return 0;
    }

  /* Read data to description */

  if (fread(description, sizeof(char), nholl*mult, fp) != nholl*mult)
    return 1;

  /* Null terminate the string */

  description[nholl*mult] = '\0';

  return 0;
}

/*****************************************************************************/

long COVERXreadEGrid(FILE *fp, FILE *dbgp, long ngroup, double *grid, long asciiflag)
{
  int dum;
  double *tmpgrid;
  float tmpfloat;
  long expect, i;

  if (asciiflag == YES)
    return COVERXreadEGridASCII(fp, dbgp, ngroup, grid);

  /* Allocate memory for a temporary grid used for reversal*/

  tmpgrid = Mem(MEM_ALLOC, ngroup+1, sizeof(double));

  /* Expected size of block */

  expect = (ngroup+1)*sizeof(float);

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadEGrid", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Read data to temporary grid */

  for (i = 0; i < ngroup+1; i++)
    {
      /* Read block data */

      if (fread(&tmpfloat, sizeof(float), 1, fp) != 1)
        return 1;

      /* Store to temporary grid */

      tmpgrid[i] = (double)tmpfloat;
    }

  /* Read size of block (end) */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadEGrid", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  /* Reverse the energy grid and convert from eV to MeV*/

  for (i = 0; i < ngroup+1; i++)
    grid[i] = tmpgrid[ngroup-i]*1e-6;

  /* Free the temporary grid */

  Mem(MEM_FREE, tmpgrid);

  return 0;
}

/*****************************************************************************/

long COVERXreadEGridASCII(FILE *fp, FILE *dbgp, long ngroup, double *grid)
{
  double *tmpgrid;
  float tmpfloat;
  long i;
  char tmpstr[MAX_STR], expect[3] = "3d\0";

  /* Read block number */

  if (fscanf(fp, "%s", tmpstr) != 1)
    return 1;

  /* Check block number */

  if (strcasecmp(tmpstr, expect) != 0)
    Die("COVERXreadFC", "Unexpected upcoming block number %s, expected %s.",
        tmpstr, expect);

  /* Allocate memory for a temporary grid used for reversal*/

  tmpgrid = Mem(MEM_ALLOC, ngroup+1, sizeof(double));

  /* Read data to temporary grid */

  for (i = 0; i < ngroup+1; i++)
    {
      /* Read block data */

      if (fscanf(fp, "%f", &tmpfloat) != 1)
        return 1;

      /* Store to temporary grid */

      tmpgrid[i] = (double)tmpfloat;
    }

  /* Reverse the energy grid and convert from eV to MeV*/

  for (i = 0; i < ngroup+1; i++)
    grid[i] = tmpgrid[ngroup-i]*1e-6;

  /* Free the temporary grid */

  Mem(MEM_FREE, tmpgrid);

  return 0;
}

/*****************************************************************************/

long COVERXreadMMTControl(FILE *fp, FILE *dbgp, long nmmp, long asciiflag)
{
  int dum;
  long expect, i;

  if (asciiflag == YES)
    return COVERXreadMMTControlASCII(fp, dbgp, nmmp);

  /* Expected size of block */

  expect = 3*nmmp*sizeof(int);

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadMMTControl", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Read and skip data (MATID, MTID, MTWGT (weighting of condensation)) */

  for (i = 0; i < 3*nmmp; i++)
    if (fread(&dum, sizeof(int), 1, fp) != 1)
      return 1;

  /* Read size of block (end) */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadMMTControl", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  return 0;
}

/*****************************************************************************/

long COVERXreadMMTControlASCII(FILE *fp, FILE *dbgp, long nmmp)
{
  long i, j;
  char tmpstr[MAX_STR], expect[3] = "5d\0";
  char word[7], tmpchar;

  /* Read block number */

  if (fscanf(fp, "%s", tmpstr) != 1)
    return 1;

  /* Check block number */

  if (strcasecmp(tmpstr, expect) != 0)
    Die("COVERXreadMMTControlASCII", "Unexpected upcoming block number %s, expected %s.",
        tmpstr, expect);

  /* Skip one space */

  if (fread(word, sizeof(char), 1, fp) != 1)
        return 1;

  /* Read and skip data (MATID, MTID, MTWGT (weighting of condensation)) */

  for (i = 0; i < 3*nmmp; i++)
    {
      for (j = 0; j < 6; j++)
        {
          if (fread(&tmpchar, sizeof(char), 1, fp) != 1)
            return 1;

          /* Skip newlines */

          if (tmpchar == '\n')
            j--;
          else
            word[j] = tmpchar;

        }

      /* Null terminate */

      word[6] = '\0';
    }

  return 0;
}

/*****************************************************************************/

long COVERXreadMMTData(FILE *fp, FILE *dbgp, long nmmp, long ngroup, long asciiflag)
{
  int dum;
  float tmpfloat;
  long expect, i, j;

  if (asciiflag == YES)
    return COVERXreadMMTDataASCII(fp, dbgp, nmmp, ngroup);

  for (i = 0; i < nmmp; i++)
    {
      /* Expected size of block */

      expect = ngroup*2*sizeof(float);

      /* Read size of block */

      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;

      /* Check block size */

      if (dum != expect)
        Die("COVERXreadMMTData", "Unexpected upcoming block size %ld, expected %ld.",
            dum, expect);

      /* Read and skip data (XSdata and error data) */

      for (j = 0; j < ngroup*2; j++)
        if (fread(&tmpfloat, sizeof(float), 1, fp) != 1)
          return 1;

      /* Read size of block (end) */

      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;

      /* Check block size */

      if (dum != expect)
        Die("COVERXreadMMTData", "Unexpected past block size %ld, expected %ld.",
            dum, expect);
    }

  return 0;
}

/*****************************************************************************/

long COVERXreadMMTDataASCII(FILE *fp, FILE *dbgp, long nmmp, long ngroup)
{
  float tmpfloat;
  long i, j;
  char tmpstr[MAX_STR], expect[3] = "6d\0";

  for (i = 0; i < nmmp; i++)
    {
      /* Read block number */

      if (fscanf(fp, "%s", tmpstr) != 1)
        return 1;

      /* Check block number */

      if (strcasecmp(tmpstr, expect) != 0)
        Die("COVERXreadMMTDataASCII", "Unexpected upcoming block number %s, expected %s.",
            tmpstr, expect);

      /* Read and skip data (XSdata and error data) */

      for (j = 0; j < ngroup*2; j++)
        if (fscanf(fp, "%f", &tmpfloat) != 1)
          return 1;

    }

  return 0;
}

/*****************************************************************************/

long COVERXreadMatrixControl(FILE *fp, FILE *dbgp, long matrixcontrol[5], long asciiflag)
{
  int dum;
  long expect, i;

  if (asciiflag == YES)
    return COVERXreadMatrixControlASCII(fp, dbgp, matrixcontrol);

  /* Expected size of block */

  expect = 5*sizeof(int);

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadMatrixControl", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Read data to array */

  for (i = 0; i < 5; i++)
    {
      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;

      /* Cast to long and store */

      matrixcontrol[i] = (long)dum;
    }

  /* Read size of block (end) */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadMatrixControl", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  return 0;
}

/*****************************************************************************/

long COVERXreadMatrixControlASCII(FILE *fp, FILE *dbgp, long matrixcontrol[5])
{
  int dum;
  long i, j;
  char tmpchar, tmpstr[MAX_STR], word[7], expect[3] = "7d\0";

  /* Read block number */

  if (fscanf(fp, "%s", tmpstr) != 1)
    return 1;

  /* Check block number */

  if (strcasecmp(tmpstr, expect) != 0)
    Die("COVERXreadMatrixControlASCII", "Unexpected upcoming block number %s, expected %s.",
        tmpstr, expect);

  /* Skip one space */

  if (fread(word, sizeof(char), 1, fp) != 1)
        return 1;

  /* Read data to array */

  for (i = 0; i < 5; i++)
    {
      for (j = 0; j < 6; j++)
        {
          if (fread(&tmpchar, sizeof(char), 1, fp) != 1)
            return 1;

          /* Skip newlines */

          if (tmpchar == '\n')
            j--;
          else
            word[j] = tmpchar;

        }

      /* Null terminate */

      word[6] = '\0';

      /* Scan integer from word */

      if (sscanf(word, "%d", &dum) != 1)
        return 1;

      /* Cast to long and store */

      matrixcontrol[i] = (long)dum;
    }

  return 0;
}

/*****************************************************************************/

long COVERXreadMBC(FILE *fp, FILE *dbgp, long ngroup, long nblock,
                   long *bandwidth, long *diagpos, long asciiflag)
{
  int dum;
  long expect, i;

  if (asciiflag == YES)
    return COVERXreadMBCASCII(fp, dbgp, ngroup, nblock, bandwidth, diagpos);

  /* Expected size of block */

  expect = (2*ngroup+nblock)*sizeof(int);

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadMBC", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Read data to array */

  for (i = 0; i < ngroup; i++)
    {
      /* Read JBAND(i) -- bandwith for group i*/

      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;
      bandwidth[i] = (long)dum;

      /* Read IJJ(i) -- position of diagonal element for group i*/

      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;
      diagpos[i] = (long)dum;
    }

  for (i = 0; i < nblock; i++)
    {
      /* Read LGRP(i) -- number of groups in block i */

      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;

    }

  /* Read size of block (end) */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadMBC", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  return 0;
}

/*****************************************************************************/

long COVERXreadMBCASCII(FILE *fp, FILE *dbgp, long ngroup, long nblock,
                   long *bandwidth, long *diagpos)
{
  int dum;
  long i;
  char tmpstr[MAX_STR], expect[3] = "8d\0";

  /* Read block number */

  if (fscanf(fp, "%s", tmpstr) != 1)
    return 1;

  /* Check block number */

  if (strcasecmp(tmpstr, expect) != 0)
    Die("COVERXreadMBCASCII", "Unexpected upcoming block number %s, expected %s.",
        tmpstr, expect);

  /* Read data to array */

  for (i = 0; i < ngroup; i++)
    {
      /* Read JBAND(i) -- bandwith for group i*/

      if (fscanf(fp, "%d", &dum) != 1)
        return 1;

      /* Cast to long and store */

      bandwidth[i] = (long)dum;

      /* Read IJJ(i) -- position of diagonal element for group i*/

      if (fscanf(fp, "%d", &dum) != 1)
        return 1;

      /* Cast to long and store */

      diagpos[i] = (long)dum;
    }

  for (i = 0; i < nblock; i++)
    {
      /* Read LGRP(i) -- number of groups in block i */

      if (fscanf(fp, "%d", &dum) != 1)
        return 1;

    }

  return 0;
}

/*****************************************************************************/

long COVERXreadMatrix(FILE *fp, FILE *dbgp, long ngroup, long *bandwidth,
                      long *diagpos, double *covdata, long asciiflag)
{
  int dum;
  float tmpfloat;
  long expect, i, j, offset;

  if (asciiflag == YES)
    return COVERXreadMatrixASCII(fp, dbgp, ngroup, bandwidth, diagpos, covdata);

  /* Calculate the expected size of block */

  expect = 0;

  for (i = 0; i < ngroup; i++)
    expect += bandwidth[i];

  expect *= sizeof(int);

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadMatrix", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Loop over all groups to read the rows */

  for (i = 0; i < ngroup; i++)
    {
      /* Calculate offset to be used for this group */
      /* FORTRAN indexing starts from 1 */

      offset = i - (diagpos[i] - 1);

      /* Read nonzero data for this group */

      for (j = 0; j < bandwidth[i]; j++)
        {
          if(fread(&tmpfloat, sizeof(float), 1, fp) != 1)
            return 1;

          covdata[(i*ngroup) + (j + offset)] = (double)tmpfloat;
        }
    }

  /* Read size of block (end) */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXreadMatrix", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  return 0;
}

/*****************************************************************************/

long COVERXreadMatrixASCII(FILE *fp, FILE *dbgp, long ngroup, long *bandwidth,
                      long *diagpos, double *covdata)
{
  float tmpfloat;
  long i, j, offset;
  char tmpstr[MAX_STR], expect[3] = "9d\0";

  /* Read block number */

  if (fscanf(fp, "%s", tmpstr) != 1)
    return 1;

  /* Check block number */

  if (strcasecmp(tmpstr, expect) != 0)
    Die("COVERXreadMatrixASCII", "Unexpected upcoming block number %s, expected %s.",
        tmpstr, expect);

  /* Loop over all groups to read the rows */

  for (i = 0; i < ngroup; i++)
    {
      /* Calculate offset to be used for this group */
      /* FORTRAN indexing starts from 1 */

      offset = i - (diagpos[i] - 1);

      /* Read nonzero data for this group */

      for (j = 0; j < bandwidth[i]; j++)
        {
          if(fscanf(fp, "%f", &tmpfloat) != 1)
            return 1;

          covdata[(i*ngroup) + (j + offset)] = (double)tmpfloat;
        }
    }

  return 0;
}

/*****************************************************************************/

long COVERXskipMatrix(FILE *fp, FILE *dbgp, long ngroup, long nblock, long asciiflag)
{
  int dum;
  float tmpfloat;
  long expect, nonzeros, i;

  if (asciiflag == YES)
    return COVERXskipMatrixASCII(fp, dbgp, ngroup, nblock);

  /************************/
  /* Matrix block control */
  /************************/

  /* Expected size of block */

  expect = (2*ngroup+nblock)*sizeof(int);

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXskipMatrix", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Reset number of nonzero elements */

  nonzeros = 0;

  /* Read data to array */

  for (i = 0; i < ngroup; i++)
    {
      /* Read JBAND(i) -- bandwith for group i*/

      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;
      nonzeros += (long)dum;

      /* Read IJJ(i) -- position of diagonal element for group i*/

      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;
    }

  for (i = 0; i < nblock; i++)
    {
      /* Read LGRP(i) -- number of groups in block i */

      if (fread(&dum, sizeof(int), 1, fp) != 1)
        return 1;
    }

  /* Read size of block (end) */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXskipMatrix", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  /*****************************************************/
  /* Move on to reading and discarding the matrix data */
  /*****************************************************/

  /* Calculate the expected size of block */

  expect = nonzeros*sizeof(int);

  /* Read size of block */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXskipMatrix", "Unexpected upcoming block size %ld, expected %ld.",
        dum, expect);

  /* Loop over all groups to read the rows */

  for (i = 0; i < nonzeros; i++)
    if(fread(&tmpfloat, sizeof(float), 1, fp) != 1)
      return 1;

  /* Read size of block (end) */

  if (fread(&dum, sizeof(int), 1, fp) != 1)
    return 1;

  /* Check block size */

  if (dum != expect)
    Die("COVERXskipMatrix", "Unexpected past block size %ld, expected %ld.",
        dum, expect);

  return 0;
}

/*****************************************************************************/

long COVERXskipMatrixASCII(FILE *fp, FILE *dbgp, long ngroup, long nblock)
{
  int dum;
  float tmpfloat;
  long nonzeros, i;
  char tmpstr[MAX_STR], expect[3] = "9d\0";

  /************************/
  /* Matrix block control */
  /************************/

  sprintf(expect, "8d");

  /* Read block number */

  if (fscanf(fp, "%s", tmpstr) != 1)
    return 1;

  /* Check block number */

  if (strcasecmp(tmpstr, expect) != 0)
    Die("COVERXreadMatrixASCII", "Unexpected upcoming block number %s, expected %s.",
        tmpstr, expect);

  /* Reset number of nonzero elements */

  nonzeros = 0;

  /* Read data to array */

  for (i = 0; i < ngroup; i++)
    {
      /* Read JBAND(i) -- bandwith for group i*/

      if (fscanf(fp, "%d", &dum) != 1)
        return 1;
      nonzeros += (long)dum;

      /* Read IJJ(i) -- position of diagonal element for group i*/

      if (fscanf(fp, "%d", &dum) != 1)
        return 1;
    }

  for (i = 0; i < nblock; i++)
    {
      /* Read LGRP(i) -- number of groups in block i */

      if (fscanf(fp, "%d", &dum) != 1)
        return 1;
    }

  /*****************************************************/
  /* Move on to reading and discarding the matrix data */
  /*****************************************************/

  sprintf(expect, "9d");

  /* Read block number */

  if (fscanf(fp, "%s", tmpstr) != 1)
    return 1;

  /* Check block number */

  if (strcasecmp(tmpstr, expect) != 0)
    Die("COVERXreadMatrixASCII", "Unexpected upcoming block number %s, expected %s.",
        tmpstr, expect);

  /* Loop over all groups to read the rows */

  for (i = 0; i < nonzeros; i++)
    if(fscanf(fp, "%f", &tmpfloat) != 1)
      return 1;

  return 0;
}

/*****************************************************************************/

long SCALEIDtoZAI(long ID)
{
  long i, Z, A, ZA;
  long metastableIDs[] = {
    1047110,
    1095242,
    1095244,
    1048115,
    1027058,
    1099254,
    1067166,
    1061148,
    1052127,
    1052129,
    0
  };
  long boundIDs[] = {
    1013027, /* albound (?) */
    5004009, /* be-beo */
    3004009, /* bebound (?)  */
    4001002, /* d-cryo_ortho */
    5001002, /* d-cryo_para  */
    8001002, /* dfreegas     */
    1026000, /* febound (?)  */
    3006000, /* graphite     */
    6001001, /* h-benzene    */
    5006000, /* c-benzene (?), listed as h-benzene/duplicate to c */
    4001001, /* h-cryo_ortho */
    5001001, /* h-cryo_para  */
    1001001, /* h-liquid_ch4 */
    9001001, /* h-poly       */
    2001001, /* h-solid_ch4  */
    7001001, /* h-zrh2       */
    8001001, /* hfreegas     */
    5008016, /* o-beo        */
    1008016, /* o-uo2        */
    14728,   /* si-28 in SiO2 (old format?) */
    14729,   /* si-29 in SiO2 (old format?) */
    14730,   /* si-30 in SiO2 (old format?) */
    1014028, /* si-28 in SiO2 (old format?) */
    1014029, /* si-29 in SiO2 (old format?) */
    1014030, /* si-30 in SiO2 (old format?) */
    1092235, /* u-uo2        */
    1040090, /* zr-90-zr5h8  */
    1040091, /* zr-91-zr5h8  */
    1040092, /* zr-92-zr5h8  */
    1040093, /* zr-93-zr5h8  */
    1040094, /* zr-94-zr5h8  */
    1040095, /* zr-95-zr5h8  */
    1040096, /* zr-96-zr5h8  */
    0
  };

  /* Conversion based on SCALE 6.2 covariance libraries */
  /* (Tbl. 10.2.1 in SCALE 6.2.1 manual) */

  if (ID < 1000000)
    {
      /* Get Z and A */

      Z = (long)((double)ID/100.0);
      A = ID - 100*Z;

      /* Some silicon entries are bound scatterers identified by A */

      if (A < 299)
        {
          /* Non-bound non-metastable state, should be ZA directly */

          /* Additional check */

          if ((ID < 1001) || (ID > 100255))
            Die("SCALEIDtoZAI", "Cannot translate SCALE ID %ld to ZAI", ID);

          /* Non-metastable state -> I is 0 */

          return ID*10;
        }
    }

  /****************************************************/
  /* Check if the ID is in the list of metastable IDs */
  /****************************************************/

  i = 0;
  while (metastableIDs[i] != 0)
    {
      /* Check match */

      if (ID == metastableIDs[i])
        break;

      /* Compare next */

      i++;
    }

  /* Check if found */

  if (metastableIDs[i] != 0)
    {
      /* Remove millions (complex flag) from ID */

      ID = ID - 1000000*(long)((double)ID/1000000);

      /* Add I */

      return ID*10+1;
    }

  /***********************************************/
  /* Check if the ID is in the list of bound IDs */
  /***********************************************/

  i = 0;
  while (boundIDs[i] != 0)
    {
      /* Check match */

      if (ID == boundIDs[i])
        break;

      /* Compare next */

      i++;
    }

  /* Check if found */

  if (boundIDs[i] != 0)
    {
      /* Remove millions (complex flag) from ID (if present) */

      ZA = ID - 1000000*(long)((double)ID/1000000);

      /* Get Z and A */

      Z = (long)((double)ZA/100.0);
      A = ZA - 100*Z;

      /* Check for a larger than normal A (Si in SiO2) */

      if (A > 299)
        {
          /* Check that A is between 700 and 800 */

          if ((A < 700) || (A > 800))
            Die("SCALEIDtoZAI", "Could not translate ID %ld", ID);

          /* Need to remove 700 from A */

          A = A - 700;

        }

      /* Add some strange I for the bound scatterers */

      return Z*1000+A*10+9;
    }

  /* For now, let's die on unknown IDs */

  Die("SCALEIDtoZAI", "Could not translate ID %ld to ZAI", ID);

  return 0;
}
