/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : deffcordata.c                                  */
/*                                                                           */
/* Created:       2016/07/09 (TKa)                                           */
/* Last modified: 2017/02/25 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Functions for reading density effect correction              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

/*****************************************************************************/

void DeffcorDataRead(DeffcorData *deffcorData) {
  /* Reads density effect correction data into DeffcorData struct.
   * */
  const int nbuf = 5000;  /* NOTE: same as linebuf length */
  long i, Z, Zread, Zidx;
  char fname[MAX_STR], linebuf[5000];
  static char * const FUNCTION_NAME = "DeffcorDataRead:";
  DeffcorElemData *deffcorZ;
  FILE *fp;

  fprintf(outp, " - reading density effect correction data...\n");

  if (!deffcorData)
    Die(FUNCTION_NAME, "NULL pointer struct");

  /* Initialize the struct to zero */
  memset(deffcorData, 0, sizeof(DeffcorData));


  /* Set file name */
  sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
          GetText(DATA_ELECTRON_DEFFCOR_FNAME));

  /* Open the file */
  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  /* Skip the header */
  while (fgets(linebuf, nbuf, fp) && (linebuf[0] == '%'));

  /* Read the number of elements */
  if (sscanf(linebuf, "Nelement %ld\n", &deffcorData->sz) != 1)
    Die(FUNCTION_NAME, "Can't find the number of elements in file %s", fname);

  CheckValue(FUNCTION_NAME, "Nelement", "", deffcorData->sz, PHOTON_ZMAX, 120);

  
  /* Allocate memory */
  deffcorData->elemData = Mem(MEM_ALLOC, deffcorData->sz, sizeof(DeffcorElemData));
  deffcorData->map = (long *)Mem(MEM_ALLOC, deffcorData->sz, sizeof(long));

  /* Loop over elements */
  for (Zidx = 0; Zidx < deffcorData->sz; Zidx++) {

    Z = Zidx + 1;

    /* Map Z and the index */
    deffcorData->map[Zidx] = Z;

    /* Data struct for the element */
    deffcorZ = &(deffcorData->elemData[Zidx]);

    /* Set atomic number */
    deffcorZ->Z = Z;

    /* Check the element is correct */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s", fname);

    if ((sscanf(linebuf, "Element %ld\n", &Zread) != 1) || (Zread != Z))
      Die(FUNCTION_NAME, "Can't find the element Z=%ld", Z);


    /* Read the mean excitation energy I */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s", fname);

    if (sscanf(linebuf, "I %lf\n", &deffcorZ->I) != 1)
      Die(FUNCTION_NAME, "Can't find I for the element Z=%ld", Z);

    CheckValue(FUNCTION_NAME, "I", "", deffcorZ->I, 0, INFTY);

    /* ev to MeV */
    deffcorZ->I *= 1.0e-6;

    /* Read the number of shells */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s", fname);

    if (sscanf(linebuf, "NSS %ld\n", &deffcorZ->nss) != 1)
      Die(FUNCTION_NAME, "Can't find NSS for the element Z=%ld", Z);

    CheckValue(FUNCTION_NAME, "nss", "", deffcorZ->nss, 1, 30);


    /* Allocate memory */
    deffcorZ->gsconf = (long *)Mem(MEM_ALLOC, deffcorZ->nss, sizeof(long));
    deffcorZ->ebi = (double *)Mem(MEM_ALLOC, deffcorZ->nss, sizeof(double));


    /* Read the ground state configuration and binding energies */
    for (i = 0; i < deffcorZ->nss; i++) {
      if (!fgets(linebuf, nbuf, fp))
        Die(FUNCTION_NAME, "Can't read the file %s", fname);

      if (sscanf(linebuf, "%ld %lf\n", &deffcorZ->gsconf[i], &deffcorZ->ebi[i]) != 2)
        Die(FUNCTION_NAME, "Can't find shell data for the element Z=%ld", Z);

      CheckValue(FUNCTION_NAME, "gsconf", "", deffcorZ->gsconf[i], -15, 15);
      CheckValue(FUNCTION_NAME, "ebi", "", deffcorZ->ebi[i], 0, INFTY);

      deffcorZ->ebi[i] *= 1.0e-6;  /* eV to MeV */
    }
  }

  /* Close the file */
  fclose(fp);

  /***************************************************************************/

}

/*****************************************************************************/


/*****************************************************************************/

long DeffcorDataGetIdx(DeffcorData *deffcorData, long value) {
  /* Finds the value from the map and returns the index.
   * Returns -1, if value not found.
   * */
  long i;
  static char * const FUNCTION_NAME = "DeffcorDataGetIdx:";

  if (value == 0)
    Die(FUNCTION_NAME, "Zero value not supported");

  if (!deffcorData || !deffcorData->map) {
    Warn(FUNCTION_NAME, "NULL pointer");
    return -1;
  }

  for (i = 0; i < deffcorData->sz; i++)
    if (deffcorData->map[i] == value)
      return i;

  return -1;
}

/*****************************************************************************/


/*****************************************************************************/

void DeffcorDataFree(DeffcorData *deffcorData) {
  /* Frees memory allocated for DeffcorData struct.
   * */
  long i;
  static char * const FUNCTION_NAME = "DeffcorDataFree:";

  if (!deffcorData) {
    Warn(FUNCTION_NAME, "NULL pointer");
    return;
  }

  /* Loop over element-wise data */
  for (i = 0; i < deffcorData->sz; i++) {
    if (deffcorData->elemData[i].gsconf)
      Mem(MEM_FREE, deffcorData->elemData[i].gsconf);
    if (deffcorData->elemData[i].ebi)
      Mem(MEM_FREE, deffcorData->elemData[i].ebi);
  }
  Mem(MEM_FREE, deffcorData->elemData);
  Mem(MEM_FREE, deffcorData->map);
}

/*****************************************************************************/
