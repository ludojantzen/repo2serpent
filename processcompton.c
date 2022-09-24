/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processcompton.c                               */
/*                                                                           */
/* Created:       2014/11/21 (TKa)                                           */
/* Last modified: 2018/02/21 (TKa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes Compton data                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

/* TODO: 
 * - Read Compton-profile data only if RDB[DATA_PHOTON_USE_DOPPLER] == 1
 * 
 * */

void ProcessComptonExtrapolation(double pz_end, double J_end, double I_end,
                                 double *a, double *b);

/*****************************************************************************/

void ProcessCompton(long loc0, long nuc) {

  static char * const FUNCTION_NAME = "ProcessCompton:";
  long i, j, Ndata, ptr0, ptr1, Nss, Npz, idx, pzminidx, Npzlin, Npzlog,
      idx_log_start;
  const long Npzd = 31;
  const int Nbuf = 500;  /* NOTE: same as linebuf length */
  const double pzminabs = 1.0/FS_CONST;
  double cpintinf, dummy, Elim, tmp, cpint100, pzlogstart;
  double cpint_overest_max, cpint_norm;
  double *pzdarr, *pzarr, *xt, *Sxt, *exta, *extb, *exti, *intmin, *pzlin,
      *pzlog;
  double **cpdmatrix, **cpmatrix, **cpintmatrix;
  char fname[MAX_STR], linebuf[500], strZ[5], tmpstr[100], *buf0, *buf1;
  FILE *fp;
  const long Npzincl = 15;
  const double pzincl[15] = {12., 17., 25., 35., 45., 50., 55., 65., 70., 
                             80., 90., 110., 120., 130., 1.0/FS_CONST};
  double pzdmod[4];
  double cpdmod[4];
  
  /* Check pointers */
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);


  sprintf(strZ, "%ld", (long)RDB[nuc + NUCLIDE_Z]);

  /***** Read incohrent scattering functions *********************************/

  xt = NULL;
  Sxt = NULL;

  /* Set file name */
  sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
          GetText(DATA_PHOTON_INCOH_FNAME));

  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  while (fgets(linebuf, Nbuf, fp)) {

    if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
       (strcmp(tmpstr, strZ) == 0)) {

      if (fgets(linebuf, Nbuf, fp)) {
        if (sscanf(linebuf, "Ndata %ld\n", &Ndata) == 0)
          Die(FUNCTION_NAME, "Can't read incoherent scattering functions: "
                             "Ndata");

        CheckValue(FUNCTION_NAME, "Ndata", "", (double)Ndata, 3.0, INFTY);

        /* Store array size */
        WDB[loc0 + PHOTON_DIST_COMP_NISF] = (double)Ndata;

        /* Allocate memory for momentum transfers */
        xt = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));

        /* Allocate memory for scattering functions */
        Sxt = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));

        /* Read and store data */
        for (i = 0; i < Ndata; i++) {
          if (fgets(linebuf, Nbuf, fp) != NULL) {
            if (sscanf(linebuf, "%lf %lf\n", &xt[i], &Sxt[i]) != 2)
              Die(FUNCTION_NAME, "Can't read incoherent scattering functions "
                                 "data: xt, Sxt");
          }
          else
            Die(FUNCTION_NAME, "Can't read incoherent data");
        }

      }
      else
        Die(FUNCTION_NAME, "Can't read incoherent data");

      /* Exit loop */
      break;
    }
  }
  fclose(fp);

  if (!xt || !Sxt)
    Die(FUNCTION_NAME, "Incoherent scattering function data not found");

  /***************************************************************************/


  /***** Process incoherent scattering functions *****************************/

  /* Momentum transfers */
  ptr0 = ReallocMem(DATA_ARRAY, Ndata);
  WDB[loc0 + PHOTON_DIST_COMP_ISFX] = (double)ptr0;

  for (i = 0; i < Ndata; i++)
    WDB[ptr0++] = xt[i];


  /* Log-momentum transfers */
  ptr0 = ReallocMem(DATA_ARRAY, Ndata);
  WDB[loc0 + PHOTON_DIST_COMP_LISFX] = (double)ptr0;

  /* Scattering functions */
  ptr1 = ReallocMem(DATA_ARRAY, Ndata);
  WDB[loc0 + PHOTON_DIST_COMP_LISF] = (double)ptr1;

  /* Check the first elements */
  if (xt[0] == 0.0 || Sxt[0] == 0.0) {
    /* Don't use log */
    WDB[ptr0++] = xt[0];
    WDB[ptr1++] = Sxt[0];
  }
  else {
    WDB[ptr0++] = log(xt[0]);
    WDB[ptr1++] = log(Sxt[0]);
  }

  /* Store the rest as log-log */
  for (i = 1; i < Ndata; i++) {
    WDB[ptr0++] = log(xt[i]);
    WDB[ptr1++] = log(Sxt[i]);
  }

  /* Free memory */
  Mem(MEM_FREE, xt);
  Mem(MEM_FREE, Sxt);

  /***************************************************************************/


  /***** Read and process Compton profile data  ******************************/
  
  /* Check the energy limit - below the limit all the binding energies 
   * are not in descending order (at least above Z=56), and the 
   * sampling method will fail. */
  Elim = 0.0002;

  if (RDB[DATA_PHOTON_EMIN] < Elim)
    Die(FUNCTION_NAME, "Compton scattering model doesn't work below %f MeV",
        Elim);

  Nss = 0;  

  /* Set file name */
  sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
          GetText(DATA_PHOTON_CP_FNAME));

  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  while (fgets(linebuf, Nbuf, fp)) {
    if ((sscanf(linebuf, "#S %s %*s\n", tmpstr) == 1)
        && (strcmp(strZ, tmpstr) == 0)) {

      /* Read and store number of subshells */
      if (fgets(linebuf, Nbuf, fp)) {
        if (strncmp(linebuf, "#N", strlen("#N")))
          Die(FUNCTION_NAME, "Can't read Compton profile data: #N missing");
        buf0 = &linebuf[strlen("#N")];
        Nss = strtol(buf0, &buf0, 10) - 2;
        WDB[loc0 + PHOTON_DIST_COMP_NSS] = (double)Nss;
      }
      else
        Die(FUNCTION_NAME, "Can't read Compton profile data");

      if (Nss == 0)
        Die(FUNCTION_NAME, "Couldn't find number of subshells. WTF?");

      /* Read and store number of electrons in subshells */
      if (fgets(linebuf, Nbuf, fp)) {
        if (strncmp(linebuf, "#UOCCUP", strlen("#UOCCUP")))
          Die(FUNCTION_NAME, "Can't read Compton profile data: #UOCCUP missing");
        buf0 = &linebuf[strlen("#UOCCUP")];

        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_ELN] = (double)ptr0;

        for (i = 0; i < Nss; i++) {
          WDB[ptr0++] = strtod(buf0, &buf1);
          if (buf0 == buf1)
            Die(FUNCTION_NAME, "Can't read the electron numbers for Z=%s",
                strZ);
          buf0 = buf1;
        }
      }
      else
        Die(FUNCTION_NAME, "Can't read Compton profile data");


      /* Read and store the binding energies */
      /* NOTE: These binding energies differ from the ENDF binding energies! */
      if (fgets(linebuf, Nbuf, fp)) {
        if (strncmp(linebuf, "#UBIND", strlen("#UBIND")))
          Die(FUNCTION_NAME, "Can't read Compton profile data: #UBIND missing");
        buf0 = &linebuf[strlen("#UBIND")];

        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_EBI] = (double)ptr0;

        for (i = 0; i < Nss; i++) {
          WDB[ptr0++] = strtod(buf0, &buf1) / 1.0e6; /* from eV to MeV*/
          if (buf0 == buf1)
            Die(FUNCTION_NAME, "Can't read the binding energies for Z=%s",
                strZ);
          buf0 = buf1;
        }
      }
      else
        Die(FUNCTION_NAME, "Can't read Compton profile data");


      /* Skip the header */
      if (fgets(linebuf, Nbuf, fp)) {
        if (strncmp(linebuf, "#L", strlen("#L")))
          Die(FUNCTION_NAME, "Can't read Compton profile data: #L missing");
      }
      else
        Die(FUNCTION_NAME, "Can't read Compton profile data");


      /* Allocate memory for cpmatrixd and pzarrd */
      cpdmatrix = (double **)Mem(MEM_ALLOC, Nss, sizeof(double*));
      for (i = 0; i < Nss; i++)
        cpdmatrix[i] = (double *)Mem(MEM_ALLOC, Npzd, sizeof(double));
      pzdarr = (double *)Mem(MEM_ALLOC, Npzd, sizeof(double));

      /* Read data to pzarr and cpmatrix */
      for (j = 0; j < Npzd; j++) {
        if (fgets(linebuf, Nbuf, fp)) {

          /* pz */
          buf0 = linebuf;
          pzdarr[j] = strtod(buf0, &buf1);
          if (buf0 == buf1)
            Die(FUNCTION_NAME, "Can't read pz for Z=%s", strZ);
          buf0 = buf1;

          /* Skip the sum of Compton profiles */
          dummy = strtod(buf0, &buf1);
          dummy = dummy;               /* Avoid compiler warning */
          if (buf0 == buf1)
            Die(FUNCTION_NAME, "Can't read the sum of Compton profiles for "
                               "Z=%s", strZ);
          buf0 = buf1;

          /* Compton profiles */
          for (i = 0; i < Nss; i++) {
            cpdmatrix[i][j] = strtod(buf0, &buf1);
            if (buf0 == buf1)
              Die(FUNCTION_NAME, "Can't read Compton profiles for Z=%s", strZ);
            buf0 = buf1;
          }
        }
        else
          Die(FUNCTION_NAME, "Can't read Compton profile data");
      }


      /************************************************************************/


      if (RDB[DATA_PHOTON_CP_FALLBACK_2129] == (double)YES) {

        /* Use Compton profile model from 2.1.29 (will be removed after 2.1.30) */

        /***** Interpolation and extrapolation of the Compton profile data *****/

        Npz = Npzd + Npzincl;

        cpmatrix = (double **)Mem(MEM_ALLOC, Nss, sizeof(double*));
        for (i = 0; i < Nss; i++)
          cpmatrix[i] = (double *)Mem(MEM_ALLOC, Npz, sizeof(double));
        pzarr = (double *)Mem(MEM_ALLOC, Npz, sizeof(double));

        /* Include interpolation and extrapolation points to pzarr  */
        for (i = 0; i < Npzd; i++)
          pzarr[i] = pzdarr[i];

        for (i = 0; i < Npzincl; i++)
          pzarr[Npzd + i] = pzincl[i];

        /* Sort the array */
        SortArray(pzarr, Npz);

        /* Check for duplicates */
        for (i = 0; i < Npz - 1; i++)
          if (pzarr[i] == pzarr[i+1])
            Die(FUNCTION_NAME, "Compton profile data pz contains duplicate "
                               "values");

        for (j = 0; j < Npz; j++) {

          if (pzarr[j] < pzdarr[Npzd-1]) {

            /* Interpolate log-lin */

            if (pzarr[j] == pzarr[Npz - 1])
              idx = Npz - 2;
            else if ((idx = SearchArray(pzdarr, pzarr[j], Npzd)) == -1)
              Die(FUNCTION_NAME, "Compton profile data pz not found: %.5E",
                  pzarr[j]);

            for (i = 0; i < Nss; i++) {
              if (pzdarr[idx] == 0.0)
                cpmatrix[i][j] = cpdmatrix[i][idx];
              else
                cpmatrix[i][j] = ENDFInterp(4, pzarr[j], pzdarr[idx],
                                            pzdarr[idx+1], cpdmatrix[i][idx],
                                            cpdmatrix[i][idx+1]);
            }
          }
          else {

            /* Extrapolate log-lin */

            for (i = 0; i < Nss; i++) {
              cpmatrix[i][j] = cpdmatrix[i][Npzd-1] *
                  pow(cpdmatrix[i][Npzd-1]/cpdmatrix[i][Npzd-2],
                  (pzarr[j] - pzdarr[Npzd-1])/(pzdarr[Npzd-1] - pzdarr[Npzd-2]));
            }
          }
        }

        /* Check that all the values are positive */
        for (j = 0; j < Npz; j++) {
          for (i = 0; i < Nss; i++) {
            if (cpmatrix[i][j] < 0.0)
              Die(FUNCTION_NAME, "Negative value in the Compton profile data:\n "
                                 "cpmatrix[%ld][%ld] = %.f", i, j, cpmatrix[i][j]);
          }
        }

        /* Find the index of pzminabs if it's included in pzincl */
        pzminidx = -1;
        for (i = 0; i < Npzincl; i++) {
          if (pzincl[i] == pzminabs) {
            for (j = 0; j < Npz; j++) {
              if (pzarr[j] == pzminabs) {
                pzminidx = j;
                break;
              }
            }
            break;
          }
        }

        if (pzminidx == -1)
          Die(FUNCTION_NAME, "pzminabs not found");

        /***********************************************************************/

        /***** Normalize and integrate Compton profiles ************************/

        cpintmatrix = (double **)Mem(MEM_ALLOC, Nss, sizeof(double*));
        for (i = 0; i < Nss; i++)
          cpintmatrix[i] = (double *)Mem(MEM_ALLOC, Npz, sizeof(double));
        exta = (double *)Mem(MEM_ALLOC, Nss, sizeof(double));

        for (i = 0; i < Nss; i++) {

          /* Numerical integration using trapezoidal rule */
          cpintmatrix[i][0] = 0.0;
          for (j = 1; j < Npz; j++)
            cpintmatrix[i][j] = cpintmatrix[i][j-1] + 0.5*(pzarr[j] - pzarr[j-1])
                                *(cpmatrix[i][j] + cpmatrix[i][j-1]);

          /* Extrapolation coefficient */
          exta[i] = log(cpdmatrix[i][Npzd-1]/cpdmatrix[i][Npzd-2])
              /(pzdarr[Npzd-1] - pzdarr[Npzd-2]);

          /* Check that exta is negative */
          if (exta[i] >= 0)
            Die(FUNCTION_NAME, "Positive extrapolation coefficient a in "
                               "exp(a*pz)");

          /* Integral from -infty to infty */
          cpintinf = 2.0*(cpintmatrix[i][Npz-1] - cpmatrix[i][Npz-1]/exta[i]);

          /* Check the integral */
          if ((cpintinf < 0.95) || (cpintinf > 1.06))
            Die(FUNCTION_NAME, "Compton profile integral out of limits %f",
                cpintinf);

          /* Normalize Compton profiles and the integrals */
          for (j = 0; j < Npz; j++) {
            cpmatrix[i][j] = cpmatrix[i][j]/cpintinf;
            cpintmatrix[i][j] = cpintmatrix[i][j]/cpintinf;
          }
        }

        /*********************************************************************/

        /***** Store data ****************************************************/

        /* Store array size */
        WDB[loc0 + PHOTON_DIST_COMP_NCP] = (double)Npz;

        /* Store pz */
        ptr0 = ReallocMem(DATA_ARRAY, Npz);
        WDB[loc0 + PHOTON_DIST_COMP_CPPZ] = (double)ptr0;
        for (j = 0; j < Npz; j++)
          WDB[ptr0++] = pzarr[j];

        /* Store the index of pzminabs */
        WDB[loc0 + PHOTON_DIST_COMP_PZMINIDX] = (double)pzminidx;

        /* Store compton profiles */
        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CP] = (double)ptr0;

        for (i = 0; i < Nss; i++) {
          ptr1 = ReallocMem(DATA_ARRAY, Npz);
          WDB[ptr0++] = (double)ptr1;
          for (j = 0; j < Npz; j++)
            WDB[ptr1++] = cpmatrix[i][j];
        }

        /* Store integrated Compton profiles (=cdf) */
        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CPINT] = (double)ptr0;

        for (i = 0; i < Nss; i++) {
          ptr1 = ReallocMem(DATA_ARRAY, Npz);
          WDB[ptr0++] = (double)ptr1;
          for (j = 0; j < Npz; j++)
            WDB[ptr1++] = cpintmatrix[i][j];
        }

        /* Interpolation coefficients */
        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CPA] = (double)ptr0;

        for (i = 0; i < Nss; i++) {
          ptr1 = ReallocMem(DATA_ARRAY, Npz-1);
          WDB[ptr0++] = (double)ptr1;
          for (j = 0; j < Npz-1; j++)
            WDB[ptr1++] = (cpmatrix[i][j+1] - cpmatrix[i][j])
                          /(pzarr[j+1] - pzarr[j]);
        }

        /* Store the extrapolation coefficients */
        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CPEXTA] = (double)ptr0;
        for (i = 0; i < Nss; i++)
          WDB[ptr0++] = exta[i];

        /* Store electron shell number cdf */
        ptr0 = ReallocMem(DATA_ARRAY, Nss + 1);
        WDB[loc0 + PHOTON_DIST_COMP_ELNCDF] = (double)ptr0;
        WDB[ptr0] = 0.0;
        for (i = 0; i < Nss; i++)
          WDB[ptr0+1+i] = RDB[ptr0+i]
              + RDB[(long)RDB[loc0 + PHOTON_DIST_COMP_ELN] + i];

        /* Check the last element */
        if ((long)RDB[nuc + NUCLIDE_Z] != (long)WDB[ptr0+Nss])
          Die(FUNCTION_NAME, "incorrect electron number cdf");

        /* Free memory */
        for (i = 0; i < Nss; i++) {
          Mem(MEM_FREE, cpdmatrix[i]);
          Mem(MEM_FREE, cpmatrix[i]);
          Mem(MEM_FREE, cpintmatrix[i]);
        }
        Mem(MEM_FREE, cpdmatrix);
        Mem(MEM_FREE, cpmatrix);
        Mem(MEM_FREE, cpintmatrix);
        Mem(MEM_FREE, pzdarr);
        Mem(MEM_FREE, pzarr);
        Mem(MEM_FREE, exta);

        /* Exit loop */
        break;

        /*********************************************************************/
      }
      else {
        
        Npz = 0;
        idx_log_start = 4;
        pzlogstart = pzdarr[idx_log_start];
        Npzlin = 20;
        Npzlog = 500;
        Npz = 0;

        /* Create pz arrays for the linear and log regions */
        pzlin = MakeArray(0.0, pzdarr[idx_log_start], Npzlin, 1);
        pzlog = MakeArray(pzdarr[idx_log_start], pzdarr[Npzd-1], Npzlog, 2);

        /* Set the exact values at the ends to avoid numerical problems */
        pzlin[Npzlin-1] = pzdarr[idx_log_start];
        pzlog[0] = pzdarr[idx_log_start];
        pzlog[Npzlog-1] = pzdarr[Npzd-1];

        /* Merge arrays */
        pzarr = AddPts(NULL, &Npz, pzdarr, Npzd);
        pzarr = AddPts(pzarr, &Npz, pzlin, Npzlin);
        pzarr = AddPts(pzarr, &Npz, pzlog, Npzlog);

        Mem(MEM_FREE, pzlin);
        Mem(MEM_FREE, pzlog);

        /***** Interpolation and extrapolation of the Compton profile data *****/

        /* Allocate memory for Compton profiles */

        cpmatrix = (double **)Mem(MEM_ALLOC, Nss, sizeof(double*));
        for (i = 0; i < Nss; i++)
          cpmatrix[i] = (double *)Mem(MEM_ALLOC, Npz, sizeof(double));


        pzdmod[0] = pzdarr[Npzd - 3];
        pzdmod[1] = pzdarr[Npzd - 2];
        pzdmod[2] = pzdarr[Npzd - 1];
        pzdmod[3] = 1.1*pzdmod[2];

        /* Loop over shell profiles */
        for (i = 0; i < Nss; i++) {

          /* Loop over pz */
          for (j = 0; j < Npz; j++) {

            if (pzarr[j] < pzlogstart) {

              /* Interpolation on lin-lin scale */

              cpmatrix[i][j] = LagrangeInterpCubic(pzdarr, cpdmatrix[i], Npzd,
                                                   pzarr[j], 2);

            }
            else if (pzarr[j] <= pzdarr[Npzd - 2]) {

              /* Interpolation on log-log scale */

              cpmatrix[i][j] = LagrangeInterpCubic(pzdarr, cpdmatrix[i], Npzd,
                                                   pzarr[j], 5);
            }
            else if (pzarr[j] > pzdarr[Npzd - 2]) {

              /* Last interval. Include an extra point to avoid oscillation. */

              cpdmod[0] = cpdmatrix[i][Npzd - 3];
              cpdmod[1] = cpdmatrix[i][Npzd - 2];
              cpdmod[2] = cpdmatrix[i][Npzd - 1];

              /* Extrapolate the Compton profile beyond the last data point */
              cpdmod[3] = ENDFInterp(5, pzdmod[3], pzdmod[1], pzdmod[2],
                  cpdmod[1], cpdmod[2]);

              cpmatrix[i][j] = LagrangeInterpCubic(pzdmod, cpdmod, 4,
                                                   pzarr[j], 5);
            }

            /* Adjust the value if adjacent Compton profiles are close to each
             * other to avoid numerical problems in the sampling method. */
            if ((j > 0) &&
                (fabs(1.0 - cpmatrix[i][j]/cpmatrix[i][j-1]) < 1e-13))
              cpmatrix[i][j] = cpmatrix[i][j-1];
          }

        }

        /* Check that all the values are positive */
        for (j = 0; j < Npz; j++) {
          for (i = 0; i < Nss; i++) {
            if (cpmatrix[i][j] < 0.0)
              Die(FUNCTION_NAME, "Negative value in the Compton profile data:\n "
                                 "cpmatrix[%ld][%ld] = %E", i, j, cpmatrix[i][j]);
          }
        }


        if (pzarr[Npz-1] >= pzminabs)
          Die(FUNCTION_NAME, "Maximum pz %f >= 1/alpha (=%f)", pzarr[Npz-1],
              pzminabs);

        /*********************************************************************/

        /***** Integrate and normalize Compton profiles **********************/

        cpintmatrix = (double **)Mem(MEM_ALLOC, Nss, sizeof(double*));
        for (i = 0; i < Nss; i++)
          cpintmatrix[i] = (double *)Mem(MEM_ALLOC, Npz, sizeof(double));
        exta = (double *)Mem(MEM_ALLOC, Nss, sizeof(double));
        extb = (double *)Mem(MEM_ALLOC, Nss, sizeof(double));
        exti = (double *)Mem(MEM_ALLOC, Nss, sizeof(double));
        intmin = (double *)Mem(MEM_ALLOC, Nss, sizeof(double));


        for (i = 0; i < Nss; i++) {

          /* Numerical integration using trapezoidal rule */
          cpintmatrix[i][0] = 0.0;
          for (j = 1; j < Npz; j++)
            cpintmatrix[i][j] = cpintmatrix[i][j-1] + 0.5*(pzarr[j] - pzarr[j-1])
                                *(cpmatrix[i][j] + cpmatrix[i][j-1]);

          /* Integral between 0 and the maximum in the data (100) */
          cpint100 = cpintmatrix[i][Npz-1];

          /* Check the integral */
          if (cpint100 > 0.505)
            Die(FUNCTION_NAME, "Compton profile integral %f out of limits",
                cpint100);


          if (cpint100 > 0.5) {
            /* The integral is overestimated, and it is set to be slightly
             * below the maximum 0.5 so that problems in sampling pz are
             * avoided. (In rare cases, pz is limited between -137 and
             * -pz_max, in which case pz can't be sampled if the integral
             * between 0 and pz_max is exactly 0.5.) Note that ?? becomes
             * very large, if cpint_overest_max is close to 0.5. */

            cpint_overest_max = 0.4999;
            cpint_norm = 2.0*cpint100*0.5/cpint_overest_max;

            for (j = 0; j < Npz; j++) {
              cpmatrix[i][j] = cpmatrix[i][j]/cpint_norm;
              cpintmatrix[i][j] = cpintmatrix[i][j]/cpint_norm;
            }
          }

          /* Solve extrapolation coefficients */
          ProcessComptonExtrapolation(pzarr[Npz-1], cpmatrix[i][Npz-1],
              cpintmatrix[i][Npz-1], &exta[i], &extb[i]);

          /* Store integration coefficient */
          exti[i] = exp(-(1.0 + extb[i]*pzarr[Npz-1])
              *(1.0 + extb[i]*pzarr[Npz-1]));

          /* Calculate integral between 0 and ~137 */
          tmp = 1.0 + extb[i]*pzminabs;
          intmin[i] = 0.5*(1.0 - exta[i]/extb[i]*exp(-tmp*tmp));
        }


        /***********************************************************************/


        /***** Store data ******************************************************/

        /* Store array size */
        WDB[loc0 + PHOTON_DIST_COMP_NCP] = (double)Npz;

        /* Store pz */
        ptr0 = ReallocMem(DATA_ARRAY, Npz);
        WDB[loc0 + PHOTON_DIST_COMP_CPPZ] = (double)ptr0;
        for (j = 0; j < Npz; j++)
          WDB[ptr0++] = pzarr[j];

        /* Store the index of pzminabs */
        WDB[loc0 + PHOTON_DIST_COMP_PZMINIDX] = -1; /* TODO: Remove in the future */

        /* Store compton profiles */
        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CP] = (double)ptr0;

        for (i = 0; i < Nss; i++) {
          ptr1 = ReallocMem(DATA_ARRAY, Npz);
          WDB[ptr0++] = (double)ptr1;
          for (j = 0; j < Npz; j++)
            WDB[ptr1++] = cpmatrix[i][j];
        }

        /* Store integrated Compton profiles (=cdf) */
        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CPINT] = (double)ptr0;

        for (i = 0; i < Nss; i++) {
          ptr1 = ReallocMem(DATA_ARRAY, Npz);
          WDB[ptr0++] = (double)ptr1;
          for (j = 0; j < Npz; j++)
            WDB[ptr1++] = cpintmatrix[i][j];
        }

        /* Interpolation coefficients */
        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CPA] = (double)ptr0;

        for (i = 0; i < Nss; i++) {
          ptr1 = ReallocMem(DATA_ARRAY, Npz-1);
          WDB[ptr0++] = (double)ptr1;
          for (j = 0; j < Npz-1; j++)
            WDB[ptr1++] = (cpmatrix[i][j+1] - cpmatrix[i][j])
                          /(pzarr[j+1] - pzarr[j]);
        }


        /* Store the extrapolation coefficients */
        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CPEXTA] = (double)ptr0;
        for (i = 0; i < Nss; i++)
          WDB[ptr0++] = exta[i];

        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CPEXTB] = (double)ptr0;
        for (i = 0; i < Nss; i++)
          WDB[ptr0++] = extb[i];

        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CPEXTI] = (double)ptr0;
        for (i = 0; i < Nss; i++)
          WDB[ptr0++] = exti[i];

        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_COMP_CPINTMIN] = (double)ptr0;
        for (i = 0; i < Nss; i++)
          WDB[ptr0++] = intmin[i];


        /* Store electron shell number cdf */
        ptr0 = ReallocMem(DATA_ARRAY, Nss + 1);
        WDB[loc0 + PHOTON_DIST_COMP_ELNCDF] = (double)ptr0;
        WDB[ptr0] = 0.0;
        for (i = 0; i < Nss; i++)
          WDB[ptr0+1+i] = RDB[ptr0+i]
              + RDB[(long)RDB[loc0 + PHOTON_DIST_COMP_ELN] + i];

        /* Check the last element */
        if ((long)RDB[nuc + NUCLIDE_Z] != (long)WDB[ptr0+Nss])
          Die(FUNCTION_NAME, "incorrect electron number cdf");

        /* Free memory */
        for (i = 0; i < Nss; i++) {
          Mem(MEM_FREE, cpdmatrix[i]);
          Mem(MEM_FREE, cpmatrix[i]);
          Mem(MEM_FREE, cpintmatrix[i]);
        }
        Mem(MEM_FREE, cpdmatrix);
        Mem(MEM_FREE, cpmatrix);
        Mem(MEM_FREE, cpintmatrix);
        Mem(MEM_FREE, pzdarr);
        Mem(MEM_FREE, pzarr);
        Mem(MEM_FREE, exta);
        Mem(MEM_FREE, extb);
        Mem(MEM_FREE, exti);
        Mem(MEM_FREE, intmin);

        /* Exit loop */
        break;

      }

    }
  }

  fclose(fp);

  /***************************************************************************/

}

/*****************************************************************************/


/*****************************************************************************/

void ProcessComptonExtrapolation(double pz_end, double J_end, double I_end,
                                 double *a, double *b) {
  /*
   * Solves coefficients a and b for the extrapolated Compton profile function
   *
   *   J = a*(1 + b*|pz|)*exp(-(1 + b*|pz|)^2).
   *
   * Coefficients are solved by using the continuity requirement at pz_end,
   * and the integral
   *
   * Arguments:
   *
   *  pz_end  :  The last value of pz array in the Compton profile data.
   *  J_end   :  Value of the Compton profile at pz_end.
   *  I_end   :  Integral of the Compton profile between 0 and pz_end.
   *
   * */
  static char * const FUNCTION_NAME = "ProcessComptonExtrapolation:";
  double D, Iextrap;

  /* Check arguments */
  if (pz_end <= 0.0)
    Die(FUNCTION_NAME, "Non-positive argument pz_end %E", pz_end);
  if (J_end < 0.0)
    Die(FUNCTION_NAME, "Negative argument J_end %E", J_end);
  if (I_end >= 0.5)
    Die(FUNCTION_NAME, "Argument I_end above 0.5");

  /* Calculate and check discriminant */
  D = 1.0 + 4.0*pz_end*J_end/(1.0 - 2.0*I_end);

  if (D < 0.0)
    Die(FUNCTION_NAME, "Negative discriminant");

  /* Calculate b and a. Positive solution for b is used so that the derivative
   * of J is negative at pz_end. */
  *b = (-1.0 + sqrt(D))/(2.0*pz_end);
  *a = *b*(1.0 - 2.0*I_end)*exp((1.0 + *b*pz_end)*(1.0 + *b*pz_end));

  /* Integral between pz_end and infinity */
  Iextrap = *a/(2.0**b)*exp(-(1.0 + *b*pz_end)*(1.0 + *b*pz_end));

  /* Check the integral */

  if (I_end + Iextrap != 0.5) {

    if (fabs(I_end + Iextrap - 0.5) > 1.e-10) {
      Die(FUNCTION_NAME, "Incorrect Compton profile integral %.15E (error 1)",
          I_end + Iextrap);
    }
    else {
      /* Adjust the integral */
      Iextrap = fabs(I_end - 0.5);

      if (I_end + Iextrap != 0.5)
        Die(FUNCTION_NAME, "Incorrect Compton profile integral %.15E "
                           "(error 2)", I_end + Iextrap);
    }
  }
}

/*****************************************************************************/
