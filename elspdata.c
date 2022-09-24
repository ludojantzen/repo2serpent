/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : elspdata.c                                     */
/*                                                                           */
/* Created:       2016/07/09 (TKa)                                           */
/* Last modified: 2017/03/08 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Functions for reading electron stoppping power data          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

/*****************************************************************************/

void ElSPDataRead(ElSPData *SPData) {
  /* Reads electron stopping power data into ElectronSPData struct.
   *
   * Some assumptions of the data:
   * - Units:
   *    * Electron kinetic energy [MeV]
   *    * Stopping powers [MeV cm^2/g]
   * - At least three columns must exist in each table (energy, col SP, rad Sp)
   * - Maximum number of columns is five
   * - The columns are (if available):
   *     energy | col. SP | rad. SP | tot. SP | density effect parameter
   * */
  long i, j, Z, Zidx, nZ, nrows, nrowstmp, ncolumnstmp, ncolumns, firstelem;
  const int nbuf = 500;  /* NOTE: same as linebuf length */
  char fname[MAX_STR], linebuf[500], *buf0, *buf1;
  static char * const FUNCTION_NAME = "ElSPDataRead:";
  double Eprev;
  double *rowdata;
  FILE *fp;

  fprintf(outp, " - reading electron stopping power data...\n");

  if (!SPData)
    Die(FUNCTION_NAME, "NULL pointer struct");

  /* Avoid compiler warning */
  nrows = ncolumns = 0;

  /* Initialize the struct to zero */
  memset(SPData, 0, sizeof(ElSPData));

  /* Set file name */
  sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
          GetText(DATA_ELECTRON_SP_FNAME));

  /* Open the file */
  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);


  /* Find the number of elements in the data */
  nZ = 0;
  while (fgets(linebuf, nbuf, fp))
    if ((sscanf(linebuf, "Element %ld\n", &Z) == 1))
      nZ++;
  rewind(fp);

  /* Put the number of elements */
  SPData->sz = nZ;

  /* Allocate memory */
  SPData->map = (long *)Mem(MEM_ALLOC, SPData->sz, sizeof(long));
  SPData->SPcole = (double **)Mem(MEM_ALLOC, SPData->sz, sizeof(double*));
  SPData->SPrade = (double **)Mem(MEM_ALLOC, SPData->sz, sizeof(double*));
  SPData->SPtote = (double **)Mem(MEM_ALLOC, SPData->sz, sizeof(double*));
  SPData->delta = (double **)Mem(MEM_ALLOC, SPData->sz, sizeof(double*));

  firstelem = 1;
  Eprev = -INFTY;
  rowdata = NULL;

  /* Loop over elements */
  for (Zidx = 0; Zidx < nZ; Zidx++) {

    /* Check the element is correct */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s for Zidx=%ld", fname, Zidx);

    if (sscanf(linebuf, "Element %ld\n", &Z) != 1)
      Die(FUNCTION_NAME, "Can't find the element Z=%ld", Z);

    if ((Z < 1) || (Z > PHOTON_ZMAX))
      Die(FUNCTION_NAME, "Z=%ld out of bounds [1, %ld]", Z, PHOTON_ZMAX);

    /* Map Z and the index */
    SPData->map[Zidx] = Z;

    /* Read the table size */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s", fname);

    if (sscanf(linebuf, "Ndata %ld %ld", &nrowstmp, &ncolumnstmp) != 2)
      Die(FUNCTION_NAME, "Can't read the table size on line %s", linebuf);


    if (firstelem) {
      /* Set sizes */
      nrows = nrowstmp;
      ncolumns = ncolumnstmp;

      /* Check the number of columns */
      if ((ncolumns < 3) || (ncolumns > 5))
        Die(FUNCTION_NAME, "Number of columns %ld out of bounds [%ld, %ld] "
                           "for Z = %ld", ncolumns, 3, 5, Z);

      SPData->nE = nrows;

      /* Allocate memory for the energy grid */
      SPData->E = (double *)Mem(MEM_ALLOC, SPData->nE, sizeof(double));

      rowdata = (double *)Mem(MEM_ALLOC, ncolumns, sizeof(double));
    }
    else {
      /* Check the read values */
      if (nrowstmp != nrows)
        Die(FUNCTION_NAME, "Number of energies is not the same for all "
                           "the elements (%ld for Z = %ld)",
            nrowstmp, Z);
      if (ncolumnstmp != ncolumns)
        Die(FUNCTION_NAME, "Number of columns is not the same for all "
                           "the elements (%ld for Z = %ld)",
            ncolumnstmp, Z);
    }


    /* Allocate memory */
    SPData->SPcole[Zidx] = (double *)Mem(MEM_ALLOC, SPData->nE, sizeof(double));
    SPData->SPrade[Zidx] = (double *)Mem(MEM_ALLOC, SPData->nE, sizeof(double));
    SPData->SPtote[Zidx] = NULL;
    SPData->delta[Zidx] = NULL;

    /* Allocate memory for additional data if needed */
    if (ncolumns >= 4)
      SPData->SPtote[Zidx] = (double *)Mem(MEM_ALLOC, SPData->nE, sizeof(double));
    if (ncolumns == 5)
      SPData->delta[Zidx] = (double *)Mem(MEM_ALLOC, SPData->nE, sizeof(double));


    /* Loop over energy (=read rows) */
    for (i = 0; i < nrows; i++) {
      if (!fgets(linebuf, nbuf, fp))
        Die(FUNCTION_NAME, "Can't read the file %s", fname);

      /* Read the row data */
      buf0 = linebuf;
      for (j = 0; j < ncolumns; j++) {
        rowdata[j] = strtod(buf0, &buf1);
        if (buf0 == buf1)
          Die(FUNCTION_NAME, "Can't read the data for Z=%ld, column %ld", Z, j);
        buf0 = buf1;
      }

      /* Put / check energy */
      if (firstelem) {
        SPData->E[i] = rowdata[0];

        /* Check that the energy array is unique and in ascending order */
        if (SPData->E[i] <= Eprev)
          Die(FUNCTION_NAME, "Energy array is not unique or sorted for Z=%ld", Z);
        Eprev = SPData->E[i];
      }
      else if (rowdata[0] != SPData->E[i]) {
        Die(FUNCTION_NAME, "Energy grid differs for Z=%ld", Z);
      }

      /* Set stopping powers */
      SPData->SPcole[Zidx][i] = rowdata[1];
      SPData->SPrade[Zidx][i] = rowdata[2];

      /* Set additional data if available */
      if (ncolumns >= 4)
        SPData->SPtote[Zidx][i] = rowdata[3];
      if (ncolumns== 5)
        SPData->delta[Zidx][i] = rowdata[4];
    }

    firstelem = 0;
  }


  /* Close the file */
  fclose(fp);

  if (!rowdata)
    Die(FUNCTION_NAME, "Data was not found or read from the file %s", fname);

  /* Check the minimum and maximum energy */
  if (SPData->E[0] > RDB[DATA_PHOTON_EMIN])
    Die(FUNCTION_NAME, "Energy %E above the minimum %E",
        SPData->E[0], RDB[DATA_PHOTON_EMIN]);

  if (SPData->E[SPData->nE-1] < RDB[DATA_PHOTON_EMAX])
    Die(FUNCTION_NAME, "Energy %E below the maximum %E",
        SPData->E[SPData->nE-1], RDB[DATA_PHOTON_EMAX]);

  /***************************************************************************/

  /* Free memory */
  Mem(MEM_FREE, rowdata);
}

/*****************************************************************************/


/*****************************************************************************/

void ElSPDataAlloc(ElSPData *SPData, long sz, long nE) {
  /* Allocates memory for the elements of a ElSPData struct
   * Returns
   * */
  long i;
  static char * const FUNCTION_NAME = "ElSPDataAlloc:";

  if (!SPData)
    Die(FUNCTION_NAME, "NULL pointer");

  /* Initialize the struct to zero */
  memset(SPData, 0, sizeof(ElSPData));

  SPData->sz = sz;
  SPData->nE = nE;
  SPData->E = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
  SPData->map = (long *)Mem(MEM_ALLOC, sz, sizeof(long));
  SPData->SPcole = (double **)Mem(MEM_ALLOC, sz, sizeof(double*));
  SPData->SPrade = (double **)Mem(MEM_ALLOC, sz, sizeof(double*));
  SPData->SPtote = (double **)Mem(MEM_ALLOC, sz, sizeof(double*));
  SPData->SPcolp = (double **)Mem(MEM_ALLOC, sz, sizeof(double*));
  SPData->SPradp = (double **)Mem(MEM_ALLOC, sz, sizeof(double*));
  SPData->SPtotp = (double **)Mem(MEM_ALLOC, sz, sizeof(double*));
  SPData->delta = (double **)Mem(MEM_ALLOC, sz, sizeof(double*));

  for (i = 0; i < sz; i++) {
    SPData->SPcole[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double*));
    SPData->SPrade[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double*));
    SPData->SPtote[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double*));
    SPData->SPcolp[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double*));
    SPData->SPradp[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double*));
    SPData->SPtotp[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double*));
    SPData->delta[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double*));
  }

}

/*****************************************************************************/


/*****************************************************************************/

long ElSPDataGetIdx(ElSPData *SPData, long value) {
  /* Finds the value from the map and returns the index.
   * Returns -1, if value not found.
   * */
  long i;
  static char * const FUNCTION_NAME = "ElSPDataGetIdx:";

  if (value == 0)
    Die(FUNCTION_NAME, "Zero value not supported");

  if (!SPData || !SPData->map) {
    Warn(FUNCTION_NAME, "NULL pointer");
    return -1;
  }

  for (i = 0; i < SPData->sz; i++)
    if (SPData->map[i] == value)
      return i;

  return -1;
}

/*****************************************************************************/


/*****************************************************************************/

void ElSPDataFree(ElSPData *SPData) {
  /* Frees memory allocated for ElectronSPData struct
   * */
  long i;
  static char * const FUNCTION_NAME = "ElSPDataFree:";

  if (!SPData) {
    Warn(FUNCTION_NAME, "NULL pointer");
    return;
  }

  if (SPData->SPcole) {
    for (i = 0; i < SPData->sz; i++)
      Mem(MEM_FREE, SPData->SPcole[i]);
    Mem(MEM_FREE, SPData->SPcole);
  }

  if (SPData->SPrade) {
    for (i = 0; i < SPData->sz; i++)
      Mem(MEM_FREE, SPData->SPrade[i]);
    Mem(MEM_FREE, SPData->SPrade);
  }

  if (SPData->SPtote) {
    for (i = 0; i < SPData->sz; i++)
      Mem(MEM_FREE, SPData->SPtote[i]);
    Mem(MEM_FREE, SPData->SPtote);
  }

  if (SPData->SPcolp) {
    for (i = 0; i < SPData->sz; i++)
      Mem(MEM_FREE, SPData->SPcolp[i]);
    Mem(MEM_FREE, SPData->SPcolp);
  }

  if (SPData->SPradp) {
    for (i = 0; i < SPData->sz; i++)
      Mem(MEM_FREE, SPData->SPradp[i]);
    Mem(MEM_FREE, SPData->SPradp);
  }

  if (SPData->SPtotp) {
    for (i = 0; i < SPData->sz; i++)
      Mem(MEM_FREE, SPData->SPtotp[i]);
    Mem(MEM_FREE, SPData->SPtotp);
  }

  if (SPData->delta) {
    for (i = 0; i < SPData->sz; i++)
      Mem(MEM_FREE, SPData->delta[i]);
    Mem(MEM_FREE, SPData->delta);
  }

  Mem(MEM_FREE, SPData->E);
  Mem(MEM_FREE, SPData->map);

}

/*****************************************************************************/
