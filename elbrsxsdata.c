/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : elbrsxsdata.c                                  */
/*                                                                           */
/* Created:       2016/07/09 (TKa)                                           */
/* Last modified: 2017/03/08 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Functions for reading scaled electron bremsstrahlung cross   */
/*              section data                                                 */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

/*****************************************************************************/

void ElBrSXSDataRead(ElBrSXSData *SXSData) {
  /* Reads scaled bremsstrahlung data.
   * */
  const int nbuf = 500;  /* NOTE: same as linebuf length */
  long i, j, Z, Zidx, Zread, nEtmp, putE;
  double Etmp;
  char fname[MAX_STR], linebuf[500];
  static char * const FUNCTION_NAME = "ElBrSXSDataRead:";
  FILE *fp;

  fprintf(outp, " - reading bremsstrahlung cross section data...\n");

  if (!SXSData)
    Die(FUNCTION_NAME, "NULL pointer struct");

  /* Initialize the struct to zero */
  memset(SXSData, 0, sizeof(ElBrSXSData));

  /* Set the number of elements */
  SXSData->sz = PHOTON_ZMAX;

  /* Set file name */
  sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
          GetText(DATA_ELECTRON_BR_FNAME));

  /* Open the file */
  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  /* Read number of kappa */
  if (!fgets(linebuf, nbuf, fp))
    Die(FUNCTION_NAME, "Can't read the file %s", fname);

  if (sscanf(linebuf, "Kappa %ld", &SXSData->nkappa) == 0)
    Die(FUNCTION_NAME, "Can't read the size of kappa");

  /* Allocate memory */
  SXSData->kappa = (double *)Mem(MEM_ALLOC, SXSData->nkappa, sizeof(double));
  SXSData->map = (long *)Mem(MEM_ALLOC, SXSData->sz, sizeof(long));

  /* Read kappa */
  for (i = 0; i < SXSData->nkappa; i++) {
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s", fname);
    if (sscanf(linebuf, "%lf", & SXSData->kappa[i]) == 0)
      Die(FUNCTION_NAME, "Can't read kappa");
  }

  SXSData->SXSe = (double ***)Mem(MEM_ALLOC, SXSData->sz, sizeof(double**));
  SXSData->SXSeT = (double ***)Mem(MEM_ALLOC, SXSData->sz, sizeof(double**));

  /* Loop over elements */
  putE = 1;

  for (Zidx = 0; Zidx < SXSData->sz; Zidx++) {

    /* Data index for the element */
    Z = Zidx + 1;

    /* Map Z and the index */
    SXSData->map[Zidx] = Z;

    /* Check the element is correct */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s for Z=%ld", fname, Z);

    if ((sscanf(linebuf, "Element %ld\n", &Zread) != 1) || (Zread != Z))
      Die(FUNCTION_NAME, "Can't find the element Z=%ld", Z);


    /* Read the data array size */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s for Z=%ld", fname, Z);

    if (sscanf(linebuf, "NTe %ld", &nEtmp) != 1)
      Die(FUNCTION_NAME, "Can't find the array size for Z=%ld", Z);


    /* Allocate memory for data arrays */
    if (putE) {
      SXSData->nE = nEtmp;
      SXSData->E = (double *)Mem(MEM_ALLOC, SXSData->nE, sizeof(double));
    }
    else if (nEtmp != SXSData->nE) {
      Die(FUNCTION_NAME, "Energy grid differs for Z=%ld", Z);
    }

    SXSData->SXSe[Zidx] = (double **)Mem(MEM_ALLOC, SXSData->nE, sizeof(double*));

    for (i = 0; i < SXSData->nE; i++) {

      /* Read the electron kinetic energy */
      if (!fgets(linebuf, nbuf, fp))
        Die(FUNCTION_NAME, "Can't read the file %s", fname);

      if (sscanf(linebuf, "Te %lf", &Etmp) != 1)
        Die(FUNCTION_NAME, "Can't find the electron energy for Z=%ld", Z);

      Etmp *= 1.0e-6;  /* eV to MeV */

      /* Put / check energy */
      if (putE)
        SXSData->E[i] = Etmp;
      else if (Etmp != SXSData->E[i])
        Die(FUNCTION_NAME, "Energy grid differs for Z=%ld", Z);


      /* Allocate memory */
      SXSData->SXSe[Zidx][i] = (double *)Mem(MEM_ALLOC, SXSData->nkappa, sizeof(double));

      for (j = 0; j < SXSData->nkappa; j++) {

        /* Scaled cross section */
        if (!fgets(linebuf, nbuf, fp))
          Die(FUNCTION_NAME, "Can't read the file %s", fname);

        if (sscanf(linebuf, "%lf", &SXSData->SXSe[Zidx][i][j]) != 1)
          Die(FUNCTION_NAME, "Can't find the scaled cross section for Z=%ld", Z);

        SXSData->SXSe[Zidx][i][j] *= 1.0e-3;   /* millibarn to barn */
      }
    }

    /* Set the transpose */
    SXSData->SXSeT[Zidx] = (double **)Mem(MEM_ALLOC, SXSData->nkappa, sizeof(double*));
    for (j = 0; j < SXSData->nkappa; j++) {
      SXSData->SXSeT[Zidx][j] = (double *)Mem(MEM_ALLOC, SXSData->nE, sizeof(double));
      for (i = 0; i < SXSData->nE; i++)
        SXSData->SXSeT[Zidx][j][i] = SXSData->SXSe[Zidx][i][j];
    }

    putE = 0;

  }

  /* Check the minimum and maximum energy */
  if (SXSData->E[0] > RDB[DATA_PHOTON_EMIN])
    Die(FUNCTION_NAME, "Minimum energy %E in the SXS data above the DATA_PHOTON_EMIN = %E",
        SXSData->E[0], RDB[DATA_PHOTON_EMIN]);

  if (SXSData->E[SXSData->nE-1] < RDB[DATA_PHOTON_EMAX])
    Die(FUNCTION_NAME, "Maximum energy %E in the SXS data below DATA_PHOTON_EMAX = %E",
        SXSData->E[SXSData->nE-1], RDB[DATA_PHOTON_EMAX]);

  /* Close the file */
  fclose(fp);
}

/*****************************************************************************/


/*****************************************************************************/

void ElBrSXSDataAlloc(ElBrSXSData *SXSData, long sz, long nE,long nkappa) {
  /* Allocates memory for the elements of ElBrSXSData struct.
   */
  long i, j;
  static char * const FUNCTION_NAME = "ElBrSXSDataAlloc:";

  if (!SXSData)
    Die(FUNCTION_NAME, "NULL pointer");

  /* Initialize the struct to zero */
  memset(SXSData, 0, sizeof(ElBrSXSData));

  SXSData->sz = sz;
  SXSData->nE = nE;
  SXSData->nkappa = nkappa;
  SXSData->E = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
  SXSData->kappa = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));
  SXSData->map = (long *)Mem(MEM_ALLOC, sz, sizeof(long));
  SXSData->SXSe = (double ***)Mem(MEM_ALLOC, sz, sizeof(double**));
  SXSData->SXSeT = (double ***)Mem(MEM_ALLOC, sz, sizeof(double**));
  SXSData->SXSp = (double ***)Mem(MEM_ALLOC, sz, sizeof(double**));
  SXSData->SXSpT = (double ***)Mem(MEM_ALLOC, sz, sizeof(double**));

  for (i = 0; i < sz; i++) {
    SXSData->SXSe[i] = (double **)Mem(MEM_ALLOC, nE, sizeof(double*));
    SXSData->SXSeT[i] = (double **)Mem(MEM_ALLOC, nkappa, sizeof(double*));
    SXSData->SXSp[i] = (double **)Mem(MEM_ALLOC, nE, sizeof(double*));
    SXSData->SXSpT[i] = (double **)Mem(MEM_ALLOC, nkappa, sizeof(double*));

    for (j = 0; j < nE; j++) {
      SXSData->SXSe[i][j] = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));
      SXSData->SXSp[i][j] = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));
    }

    for (j = 0; j < nkappa; j++) {
      SXSData->SXSeT[i][j] = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
      SXSData->SXSpT[i][j] = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
    }
  }

}

/*****************************************************************************/


/*****************************************************************************/

long ElBrSXSDataGetIdx(ElBrSXSData *SXSData, long value) {
  /* Finds the value from the map and returns the index.
   * Returns -1, if value not found.
   */
  long i;
  static char * const FUNCTION_NAME = "ElBrSXSDataGetIdx:";

  if (value == 0)
    Die(FUNCTION_NAME, "Zero value not supported");

  if (!SXSData || !SXSData->map) {
    Warn(FUNCTION_NAME, "NULL pointer");
    return -1;
  }

  for (i = 0; i < SXSData->sz; i++)
    if (SXSData->map[i] == value)
      return i;

  return -1;
}

/*****************************************************************************/


/*****************************************************************************/

void ElBrSXSDataFree(ElBrSXSData *SXSData) {
  /* Frees memory allocated for ElBrSXSData struct
   *
   * TODO: check NULLs
   */
  long i, j;
  static char * const FUNCTION_NAME = "ElBrSXSDataFree:";

  if (!SXSData) {
    Warn(FUNCTION_NAME, "NULL pointer");
    return;
  }


  if (SXSData->SXSe) {
    for (i = 0; i < SXSData->sz; i++) {
      for (j = 0; j < SXSData->nE; j++)
        Mem(MEM_FREE, SXSData->SXSe[i][j]);
      Mem(MEM_FREE, SXSData->SXSe[i]);
    }
    Mem(MEM_FREE, SXSData->SXSe);
  }

  if (SXSData->SXSeT) {
    for (i = 0; i < SXSData->sz; i++) {
      for (j = 0; j < SXSData->nkappa; j++)
        Mem(MEM_FREE, SXSData->SXSeT[i][j]);
      Mem(MEM_FREE, SXSData->SXSeT[i]);
    }
    Mem(MEM_FREE, SXSData->SXSeT);
  }

  if (SXSData->SXSp) {
    for (i = 0; i < SXSData->sz; i++) {
      for (j = 0; j < SXSData->nE; j++)
        Mem(MEM_FREE, SXSData->SXSp[i][j]);
      Mem(MEM_FREE, SXSData->SXSp[i]);
    }
    Mem(MEM_FREE, SXSData->SXSp);
  }

  if (SXSData->SXSpT) {
    for (i = 0; i < SXSData->sz; i++) {
      for (j = 0; j < SXSData->nkappa; j++)
        Mem(MEM_FREE, SXSData->SXSpT[i][j]);
      Mem(MEM_FREE, SXSData->SXSpT[i]);
    }
    Mem(MEM_FREE, SXSData->SXSpT);
  }

  Mem(MEM_FREE, SXSData->map);
  Mem(MEM_FREE, SXSData->E);
  Mem(MEM_FREE, SXSData->kappa);
}

/*****************************************************************************/
