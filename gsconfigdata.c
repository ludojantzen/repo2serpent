/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : gsconfigdata.c                                 */
/*                                                                           */
/* Created:       2016/07/09 (TKa)                                           */
/* Last modified: 2017/04/24 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Functions for reading ground state configuration             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

/*****************************************************************************/

void GSconfigDataRead(GSconfigData *gsconfig) {
  /* Reads ground state configuration data into GSconfigData struct.
   * */
  const int nbuf = 5000;  /* NOTE: same as linebuf length */
  long i, Z, Zread, Zidx;
  char fname[MAX_STR], linebuf[5000];
  static char * const FUNCTION_NAME = "GSconfigDataRead:";
  GSconfigElemData *gsconfigZ;
  FILE *fp;

  fprintf(outp, " - reading ground state configuration data...\n");

  if (!gsconfig)
    Die(FUNCTION_NAME, "NULL pointer struct");

  /* Initialize the struct to zero */
  memset(gsconfig, 0, sizeof(GSconfigData));


  /* Set file name */
  sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
          GetText(DATA_ELECTRON_GSCONFIG_FNAME));

  /* Open the file */
  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  /* Skip the header */
  while (fgets(linebuf, nbuf, fp) && (linebuf[0] == '%'));

  /* Read the number of elements */
  if (sscanf(linebuf, "Nelement %ld\n", &gsconfig->sz) != 1)
    Die(FUNCTION_NAME, "Can't find the number of elements in file %s", fname);

  CheckValue(FUNCTION_NAME, "Nelement", "", (double)gsconfig->sz, PHOTON_ZMAX,
             120);

  
  /* Allocate memory */
  gsconfig->elemData = Mem(MEM_ALLOC, gsconfig->sz, sizeof(GSconfigElemData));
  gsconfig->map = (long *)Mem(MEM_ALLOC, gsconfig->sz, sizeof(long));

  /* Loop over elements */
  for (Zidx = 0; Zidx < gsconfig->sz; Zidx++) {

    Z = Zidx + 1;

    /* Map Z and the index */
    gsconfig->map[Zidx] = Z;

    /* Data struct for the element */
    gsconfigZ = &(gsconfig->elemData[Zidx]);

    /* Set atomic number */
    gsconfigZ->Z = Z;

    /* Check the element is correct */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s", fname);

    if ((sscanf(linebuf, "Element %ld\n", &Zread) != 1) || (Zread != Z))
      Die(FUNCTION_NAME, "Can't find the element Z=%ld", Z);


    /* Read the mean excitation energy I */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s", fname);

    if (sscanf(linebuf, "I %lf\n", &gsconfigZ->I) != 1)
      Die(FUNCTION_NAME, "Can't find I for the element Z=%ld", Z);

    CheckValue(FUNCTION_NAME, "I", "", gsconfigZ->I, 0, INFTY);

    /* ev to MeV */
    gsconfigZ->I *= 1.0e-6;

    /* Read the number of shells */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read the file %s", fname);

    if (sscanf(linebuf, "NSS %ld\n", &gsconfigZ->nss) != 1)
      Die(FUNCTION_NAME, "Can't find NSS for the element Z=%ld", Z);

    CheckValue(FUNCTION_NAME, "nss", "", (double)gsconfigZ->nss, 1, 30);


    /* Allocate memory */
    gsconfigZ->gsconf = (long *)Mem(MEM_ALLOC, gsconfigZ->nss, sizeof(long));
    gsconfigZ->ebi = (double *)Mem(MEM_ALLOC, gsconfigZ->nss, sizeof(double));


    /* Read the ground state configuration and binding energies */
    for (i = 0; i < gsconfigZ->nss; i++) {
      if (!fgets(linebuf, nbuf, fp))
        Die(FUNCTION_NAME, "Can't read the file %s", fname);

      if (sscanf(linebuf, "%ld %lf\n", &gsconfigZ->gsconf[i],
                 &gsconfigZ->ebi[i]) != 2)
        Die(FUNCTION_NAME, "Can't find shell data for the element Z=%ld", Z);

      CheckValue(FUNCTION_NAME, "gsconf", "", (double)gsconfigZ->gsconf[i],
                 -15., 15.);
      CheckValue(FUNCTION_NAME, "ebi", "", gsconfigZ->ebi[i], 0, INFTY);

      gsconfigZ->ebi[i] *= 1.0e-6;  /* eV to MeV */
    }
  }

  /* Close the file */
  fclose(fp);

  /***************************************************************************/

}

/*****************************************************************************/


/*****************************************************************************/

long GSconfigDataGetIdx(GSconfigData *gsconfig, long value) {
  /* Finds the value from the map and returns the index.
   * Returns -1, if value not found.
   * */
  long i;
  static char * const FUNCTION_NAME = "GSconfigDataGetIdx:";

  if (value == 0)
    Die(FUNCTION_NAME, "Zero value not supported");

  if (!gsconfig || !gsconfig->map) {
    Warn(FUNCTION_NAME, "NULL pointer");
    return -1;
  }

  for (i = 0; i < gsconfig->sz; i++)
    if (gsconfig->map[i] == value)
      return i;

  return -1;
}

/*****************************************************************************/


/*****************************************************************************/

void GSconfigDataFree(GSconfigData *gsconfig) {
  /* Frees memory allocated for GSconfigData struct.
   * */
  long i;
  static char * const FUNCTION_NAME = "GSconfigDataFree:";

  if (!gsconfig) {
    Warn(FUNCTION_NAME, "NULL pointer");
    return;
  }

  /* Loop over element-wise data */
  for (i = 0; i < gsconfig->sz; i++) {
    if (gsconfig->elemData[i].gsconf)
      Mem(MEM_FREE, gsconfig->elemData[i].gsconf);
    if (gsconfig->elemData[i].ebi)
      Mem(MEM_FREE, gsconfig->elemData[i].ebi);
  }
  Mem(MEM_FREE, gsconfig->elemData);
  Mem(MEM_FREE, gsconfig->map);
}

/*****************************************************************************/
