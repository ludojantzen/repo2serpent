/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processrayleigh.c                              */
/*                                                                           */
/* Created:       2014/11/21 (TKa)                                           */
/* Last modified: 2017/03/08 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Processes Rayleigh scattering data                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessRayleigh:"


/*****************************************************************************/

void ProcessRayleigh(long loc0, long nuc) {

  long i, idx, Ndata, ptr;
  const int Nbuf = 500; /* NOTE: same as linebuf length */
  double x2max, c;
  double *x2, *cohff2;
  static const double zt = 6.5053004697083e+03; /* (1e6/(h*c)*angstrom)^2 */
  const double f0 = 1.0e-6;
  char fname[MAX_STR], linebuf[500], strZ[5], tmpstr[100];
  FILE *fp;

  /* Check pointers */
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  sprintf(strZ, "%ld", (long)RDB[nuc + NUCLIDE_Z]);

  /* Set file name */
  sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
          GetText(DATA_PHOTON_COH_FNAME));

  if (!(fp = fopen(fname, "r")))
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);


  while (fgets(linebuf, Nbuf, fp)){

    if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
       (strcmp(tmpstr, strZ) == 0)) {

      if (!fgets(linebuf, Nbuf, fp))
         Die(FUNCTION_NAME, "Can't read incoherent data");

      sscanf(linebuf, "Ndata %ld\n", &Ndata);
      CheckValue(FUNCTION_NAME, "Ndata", "", (double)Ndata, 3.0, INFTY);

      /* Temporary arrays for reading data */
      x2 = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));
      cohff2 = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));

      /* Read data */
      for (i = 0; i < Ndata; i++) {
        if (fgets(linebuf, Nbuf, fp) != NULL) {
          sscanf(linebuf, "%lf %lf\n", &x2[i], &cohff2[i]);

          /* Check limits */
          CheckValue(FUNCTION_NAME, "x", "", x2[i], 0.0, INFTY);
          CheckValue(FUNCTION_NAME, "cohff", "", cohff2[i], 0.0, RDB[nuc + NUCLIDE_Z]);

          x2[i] *= x2[i];
          cohff2[i] *= cohff2[i];
        }
        else
          Die(FUNCTION_NAME, "Can't read coherent data");
      }

      /* Find the maximum squared momentum transfer corresponding to
       * DATA_PHOTON_EMAX and update the array size*/
      x2max = zt*RDB[DATA_PHOTON_EMAX]*RDB[DATA_PHOTON_EMAX];
      idx = SearchArray(x2, x2max, Ndata);
      if (idx == -1)
        Die(FUNCTION_NAME, "Maximum squared momentum transfer could not found");
      Ndata = idx + 1;

      /* Store array size */
      WDB[loc0 + PHOTON_DIST_RAYL_N] = (double)Ndata;

      /* NOTE: We set the first element of the squared momentum transfer
       * array to be non-zero, so that we can use log-log interpolation.
       * The introduced error is negligible. TODO: Or should the first
       * interval be interpolated/sampled linearly? */
      x2[0] = f0*x2[1];

      /* Squared momentum transfers */
      ptr = ReallocMem(DATA_ARRAY, Ndata);
      WDB[loc0 + PHOTON_DIST_RAYL_X2] = (double)ptr;
      for (i = 0; i < Ndata; i++)
        WDB[ptr++] = x2[i];

      /* Log squared momentum transfers */
      ptr = ReallocMem(DATA_ARRAY, Ndata);
      WDB[loc0 + PHOTON_DIST_RAYL_LX2] = (double)ptr;
      for (i = 0; i < Ndata; i++)
        WDB[ptr++] = log(x2[i]);

      /* Store squared form factors */
      ptr = ReallocMem(DATA_ARRAY, Ndata);
      WDB[loc0 + PHOTON_DIST_RAYL_FF2] = (double)ptr;
      for (i = 0; i < Ndata; i++)
        WDB[ptr++] = cohff2[i];

      /* Caculate and store the interpolation constants */
      ptr = ReallocMem(DATA_ARRAY, Ndata-1);
      WDB[loc0 + PHOTON_DIST_RAYL_C] = (double)ptr;
      for (i = 0; i < Ndata-1; i++) {
        c = log(cohff2[i+1]/cohff2[i])/log(x2[i+1]/x2[i]);
        if (c == -1)
          Die(FUNCTION_NAME, "Interpolation constant c = -1");
        WDB[ptr++] = c;
      }

      /* Integrated squared form factors */
      ptr = ReallocMem(DATA_ARRAY, Ndata);
      WDB[loc0 + PHOTON_DIST_RAYL_FF2INT] = (double)ptr;
      TrapzRealCum(x2, cohff2, &WDB[ptr], Ndata, 5);

      /* Free memory */
      Mem(MEM_FREE, x2);
      Mem(MEM_FREE, cohff2);

      /* Exit loop */
      break;

    }
  }

  /* Close file */
  fclose(fp);

}

/*****************************************************************************/
