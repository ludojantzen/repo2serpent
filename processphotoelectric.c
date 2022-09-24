/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processphotoelectric.c                         */
/*                                                                           */
/* Created:       2014/11/21 (TKa)                                           */
/* Last modified: 2018/01/30 (TKa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Processes photoelectric data                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"


/*****************************************************************************/

void ProcessPhotoelectric(long loc0, long nuc) {

  static char * const FUNCTION_NAME = "ProcessPhotoelectric:";
  long ptrEb, ptrE, ptrp, ptrp2, ptrlp, ptrlp2, i, j, k, Ndata, Nssd,
      MTnumber, Nxstotd, idxss0, Nxstot, Nss, ptridxss, idxtot0, idx_tot_max,
      NMTd, flag_total_found, idx_ss;
  const int nbuf = 500;  /* NOTE: same as linebuf length */
  long *Nxsssd;
  double Etmp, xstmp;
  double *Etotd, *xstotd, *Eread, *xsread, *psssum;
  double **Essd, **xsssd;
  long fallback_2129 = 0;
  static long fallback_warning = 1;
  char fname[MAX_STR], linebuf[500], strZ[5], tmpstr[100];
  long MTset[600] = {0};
  FILE *fp;


#ifdef PRINT_PE_DATA
  /* Print photoelectric effect data for debugging */
  FILE *fp_out;
  FILE *fp_out2;

  if (!(fp_out = fopen("pe_data.out", "a")))
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  if (!(fp_out2 = fopen("pe_psssum.out", "a")))
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);
#endif

  /* Check pointers */
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Initialize variables */
  Nssd = Nxstotd = 0;
  Nxsssd = NULL;
  Essd = xsssd = NULL;
  Etotd = xstotd = NULL;

  /* Set element string */
  sprintf(strZ, "%ld", (long)RDB[nuc + NUCLIDE_Z]);


  /***** Read photoelectric xs data ******************************************/

  /* Set file name */
  sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
          GetText(DATA_PHOTON_PEXS_FNAME));

  if (!(fp = fopen(fname, "r"))) {

    fallback_2129 = 1;
    if (fallback_warning) {
      Warn(FUNCTION_NAME, "Unable to open file %s for reading (2.1.30 "
           "photoelectric xs data). Previous data format will be used if "
           "available.\nNOTE: Support for the old format will be removed "
           "in 2.1.31.", fname);
      fallback_warning = 0;
    }
  }

  if (!fallback_2129) {

    /* Skip the header */
    while (fgets(linebuf, nbuf, fp) && (linebuf[0] == '%'));

    /* Find the element */
    while ((sscanf(linebuf, "Element %s\n", tmpstr) != 1) ||
           (strcmp(tmpstr, strZ) != 0)) {
      if (!fgets(linebuf, nbuf, fp))
        Die(FUNCTION_NAME, "Can't find Z=%s in photoelectric data file: %s\n",
            strZ, fname);
    }

    /* Read number of MTs */
    if (!fgets(linebuf, nbuf, fp))
      Die(FUNCTION_NAME, "Can't read photoelectric data, file: %s, line:\n"
                         " %s", fname, linebuf);

    if (sscanf(linebuf, "NMT %ld\n", &NMTd) == 0)
      Die(FUNCTION_NAME, "Can't read photoelectric data: NMT");

    if ((NMTd == 1) || (NMTd < 0))
      Die(FUNCTION_NAME, "Error in file %s: incorrect NMT %ld", fname, NMTd);

    Nssd = 0;
    idx_ss = 0;
    flag_total_found = 0;
    Eread = xsread = NULL;


    if (NMTd > 1) {

      /* Set number of subshells */
      Nssd = NMTd - 1;
      CheckValue(FUNCTION_NAME, "NSS", "", (double)Nssd, 1, PHOTON_NSS_MAX);

      /* Subshell energies, cross sections and array sizes */
      Essd = (double **)Mem(MEM_ALLOC, Nssd, sizeof(double*));
      xsssd = (double **)Mem(MEM_ALLOC, Nssd, sizeof(double*));
      Nxsssd = (long *)Mem(MEM_ALLOC, Nssd, sizeof(long));

      /* Loop over MTs */
      for (i = 0; i < NMTd; i++) {

        /* Read MT number */
        if (!fgets(linebuf, nbuf, fp))
          Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n"
                             " line: %s", fname, linebuf);
        if (sscanf(linebuf, "MT %ld\n", &MTnumber) == 0)
          Die(FUNCTION_NAME, "Can't read photoelectric data: MT");

        /* Read number of data */
        if (!fgets(linebuf, nbuf, fp))
          Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n"
                             " line: %s", fname, linebuf);
        if (sscanf(linebuf, "Ndata %ld\n", &Ndata) == 0)
          Die(FUNCTION_NAME, "Can't read photoelectric data: Ndata");

        CheckValue(FUNCTION_NAME, "MTnumber", "", (double)Ndata, 3.0, INFTY);


        if (MTnumber == 522) {
          /* Total cross section data */

          Etotd = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));
          xstotd = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));
          Nxstotd = Ndata;
          Eread = Etotd;
          xsread = xstotd;
          flag_total_found = 1;
        }
        else if ((MTnumber >= 534) && (MTnumber <= 572) ) {
          /* Subshell data */

          Essd[idx_ss] = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));
          xsssd[idx_ss] = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));
          Nxsssd[idx_ss] = Ndata;
          Eread = Essd[idx_ss];
          xsread = xsssd[idx_ss];
          idx_ss++;

          /* Check the number of shells */
          if (idx_ss > Nssd)
            Die(FUNCTION_NAME, "Mismatch in the number of subshells for Z=%s "
                " in data file %s (error 1)", strZ, fname);
        }
        else
          Die(FUNCTION_NAME, "Incorrect MT %ld in the data file %s", MTnumber,
              fname);

        /* Check duplicate MTs */
        if (MTset[MTnumber])
          Die(FUNCTION_NAME, "Duplicate MT %ld in the data file %s", MTnumber,
              fname);
        else
          MTset[MTnumber] = 1;

        /* Read data */
        for (j = 0; j < Ndata; j++) {

          if (!fgets(linebuf, nbuf, fp))
            Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n"
                               " line: %s", fname, linebuf);

          if (sscanf(linebuf, "%lf %lf\n", &Etmp, &xstmp) != 2)
            Die(FUNCTION_NAME, "Can't read photoelectric data: Etmp, xstmp");

          /* Check data */
          if (Etmp <= 0.0)
            Die(FUNCTION_NAME, "Non-positive energy %E in subshell xs data",
                Etmp);

          if (xstmp <= 0.0)
            Die(FUNCTION_NAME, "Non-positive xs %E in subshell xs data",
                xstmp);

          if ((j > 0) && (Etmp < Eread[j-1]))
            Die(FUNCTION_NAME, "Shell energy array not in ascending order");

          Eread[j] = Etmp*1.0e-6; /* ev to MeV */
          xsread[j] = xstmp;
        }
      }
    }

    fclose(fp);

    /* Check the number of shells */
    if (idx_ss != Nssd)
      Die(FUNCTION_NAME, "Mismatch in the number of subshells for Z=%s "
          " in data file %s (error 2)", strZ, fname);
  }

  /***************************************************************************/


  if (fallback_2129) {
    /* Read photoelectric xs data in the format used in Serpent
     * versions 2.1.24 - 2.1.29 */

    /***** Read photoelectric shell xs data **********************************/

    /* Set file name */
    sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
            GetText(DATA_PHOTON_PESS_FNAME));

    if (!(fp = fopen(fname, "r")))
      Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

    while (fgets(linebuf, nbuf, fp)){

      if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
         (strcmp(tmpstr, strZ) == 0)) {

        if (!fgets(linebuf, nbuf, fp))
          Die(FUNCTION_NAME, "Can't read photoelectric data, file: %s , line:\n"
                             " %s", fname, linebuf);

        if (sscanf(linebuf, "NSS %ld\n", &Nssd) == 0)
          Die(FUNCTION_NAME, "Can't read photoelectric data: NSS");

        CheckValue(FUNCTION_NAME, "NSS", "", (double)Nssd, 0, INFTY);

        if (Nssd > 0) {

          Essd = (double **)Mem(MEM_ALLOC, Nssd, sizeof(double*));
          xsssd = (double **)Mem(MEM_ALLOC, Nssd, sizeof(double*));
          Nxsssd = (long *)Mem(MEM_ALLOC, Nssd, sizeof(long));

          /* Loop over shells */
          for (i = 0; i < Nssd; i++) {

            /* Read MT number */
            if (!fgets(linebuf, nbuf, fp))
              Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n"
                                 " line: %s", fname, linebuf);
            if (sscanf(linebuf, "MT %ld\n", &MTnumber) == 0)
              Die(FUNCTION_NAME, "Can't read photoelectric data: MT");

            CheckValue(FUNCTION_NAME, "MTnumber", "", (double)MTnumber,
                       534.0, 599.0);

            /* Read number of data */
            if (!fgets(linebuf, nbuf, fp))
              Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n"
                                 " line: %s", fname, linebuf);
            if (sscanf(linebuf, "Ndata %ld\n", &Ndata) == 0)
              Die(FUNCTION_NAME, "Can't read photoelectric data: Ndata");

            CheckValue(FUNCTION_NAME, "MTnumber", "", (double)Ndata, 3.0, INFTY);

            /* Number of data points */
            Nxsssd[i] = Ndata;

            Essd[i] = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));
            xsssd[i] = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));

            /* Read data */
            for (j = 0; j < Ndata; j++) {

              if (!fgets(linebuf, nbuf, fp))
                Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n"
                                   " line: %s", fname, linebuf);

              if (sscanf(linebuf, "%lf %lf\n", &Etmp, &xstmp) != 2)
                Die(FUNCTION_NAME, "Can't read photoelectric data: Etmp, xstmp");

              /* Check data */
              if (Etmp <= 0.0)
                Die(FUNCTION_NAME, "Non-positive energy %E in subshell xs data",
                    Etmp);
              if (xstmp <= 0.0)
                Die(FUNCTION_NAME, "Non-positive xs %E in subshell xs data",
                    xstmp);
              if ((j > 0) && (Etmp < Essd[i][j-1]))
                Die(FUNCTION_NAME, "Shell energy array not in ascending order");

              Essd[i][j] = Etmp/1.0e6;  /* From ev to MeV */
              xsssd[i][j] = xstmp;
            }
          }
        }

        /* Exit loop */
        break;

      }
    }

    fclose(fp);

    /*************************************************************************/


    /***** Read total photoelectric xs ***************************************/

    /* Here the total xs is the sum of the interpolated xs read above */

    /* Set file name */
    sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
            GetText(DATA_PHOTON_PETOT_FNAME));

    if (!(fp = fopen(fname, "r")))
      Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

    while (fgets(linebuf, nbuf, fp)) {

      if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
         (strcmp(tmpstr, strZ) == 0)) {

        if (!fgets(linebuf, nbuf, fp))
          Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n"
                             " line: %s", fname, linebuf);

        if (sscanf(linebuf, "Ndata %ld\n", &Nxstotd) == 0)
          Die(FUNCTION_NAME, "Can't read photoelectric data: Ndata");

        CheckValue(FUNCTION_NAME, "Ndata", "", (double)Nxstotd, 0, INFTY);

        Etotd = (double *)Mem(MEM_ALLOC, Nxstotd, sizeof(double));
        xstotd = (double *)Mem(MEM_ALLOC, Nxstotd, sizeof(double));

        flag_total_found = 1;

        /* Read data */
        for (j = 0; j < Nxstotd; j++) {

          if (!fgets(linebuf, nbuf, fp))
            Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n"
                               " line: %s", fname, linebuf);

          if (sscanf(linebuf, "%lf %lf\n", &Etmp, &xstmp) != 2)
            Die(FUNCTION_NAME, "Can't read photoelectric data: Etmp, xstmp");

          /* Check data */
          if (Etmp <= 0.0)
            Die(FUNCTION_NAME, "Non-positive energy %E in total xs data", Etmp);
          if (xstmp <= 0.0)
            Die(FUNCTION_NAME, "Non-positive xs %E in total xs data", xstmp);
          if ((j > 0) && (Etmp < Etotd[j-1]))
            Die(FUNCTION_NAME, "Shell energy array not in ascending order");

          Etotd[j] = Etmp/1.0e6;    /* From ev to MeV */
          xstotd[j] = xstmp;
        }

        /* Exit loop */
        break;

      }
    }

    fclose(fp);

    /* Check that total xs data was found */
    if ((Nssd > 0) && (Nxstotd == 0))
      Die(FUNCTION_NAME, "Total xs data not found for Z=%s", strZ);

    /*************************************************************************/

  }

  /***************************************************************************/


  /***** Process data ********************************************************/

  /* Find shells with binding energies above DATA_PHOTON_EMIN */
  for (Nss = 0; Nss < Nssd; Nss++)
    if (Essd[Nss][0] < RDB[DATA_PHOTON_EMIN])
      break;


  if (Nss == 0) {
    /* No shell xs data found (Nssd==0) or DATA_PHOTON_EMIN above
     * the highest binding energy */
    WDB[loc0 + PHOTON_DIST_PE_NSS] = 0;
  }
  else {

    if (!flag_total_found)
      Die("Total cross section data not found for Z=&s in data file %s",
          strZ, fname);

    /* Check that the binding energies are in descending order
     * NOTE: This fails when RDB[DATA_PHOTON_EMIN] < 1.5e-4 MeV */
    for (i = 0; i < Nss-1; i++)
      if (Essd[i][0] < Essd[i+1][0])
        Die(FUNCTION_NAME, "Binding energies are not in descending order for "
                           "Z=%s", strZ);

    /* Find the minimum energy index of total xs data */
    idxtot0 = -1;
    for (j = 0; j < Nxstotd; j++) {
      if (Etotd[j] == Essd[Nss-1][0]) {
        idxtot0 = j;
        break;
      }
    }

    if (idxtot0 < 0)
      Die(FUNCTION_NAME, "Minimum not found idxtot0=%ld, Z=%s %E", idxtot0,
          strZ, Etotd[idxtot0]);

    /* Find the maximum energy index of total xs data */
    idx_tot_max = 0;
    for (j = idxtot0; j < Nxstotd; j++) {
      if (Etotd[j] > RDB[DATA_PHOTON_EMAX]) {
        idx_tot_max = j;
        break;
      }
    }

    if ((idx_tot_max == idxtot0) || (Etotd[idx_tot_max] < RDB[DATA_PHOTON_EMAX]))
      Die(FUNCTION_NAME, "Maximum energy not found in total xs data");

    /* Energy grid size */
    Nxstot = idx_tot_max - idxtot0 + 1;

    /* Store grid size */
    WDB[loc0 + PHOTON_DIST_PE_NE] = (double)Nxstot;

    /* Store energy array */
    ptrE = ReallocMem(DATA_ARRAY, Nxstot);
    WDB[loc0 + PHOTON_DIST_PE_E] = (double)ptrE;
    for (j = idxtot0; j <= idx_tot_max; j++)
      WDB[ptrE++] = log(Etotd[j]);

    /* Store number of shells for which the binding energy is above the
     * minimum energy */
    WDB[loc0 + PHOTON_DIST_PE_NSS] = (double)Nss;

    /* Shell probability array (TODO: voi poistaa, ei käytetä missään) */
    ptrp = ReallocMem(DATA_ARRAY, Nss);
    WDB[loc0 + PHOTON_DIST_PE_PSS] = (double)ptrp;

    /* Shell log-probability array */
    ptrlp = ReallocMem(DATA_ARRAY, Nss);
    WDB[loc0 + PHOTON_DIST_PE_LPSS] = (double)ptrlp;

    /* Location indexes of shell log-probabilities on the energy grid */
    ptridxss = ReallocMem(DATA_ARRAY, Nss);
    WDB[loc0 + PHOTON_DIST_PE_IDXSS] = (double)ptridxss;

    /* Binding energies  */
    ptrEb = ReallocMem(DATA_ARRAY, Nss);
    WDB[loc0 + PHOTON_DIST_PE_EBI] = (double)ptrEb;

    /* Total probability array (used for checking only) */
    psssum = (double *)Mem(MEM_ALLOC, Nxstotd, sizeof(double));
    memset(psssum, 0, Nxstotd*sizeof(double));

    /* Loop over shells and interpolate shell xs on total xs energy grid */
    for (i = 0; i < Nss; i++) {

      /* Find the start index of the shell */
      idxss0 = idx_tot_max;
      for (j = 0; j <= idx_tot_max; j++) {
        if (Etotd[j] == Essd[i][0]) {

          /* Start index found */

          if (j == 0) {
            /* First energy value in the grid, which is stored only once in
             * the data */
            idxss0 = j;
          }
          else {
            /* The energy should be stored twice in the data */
            if (Etotd[j] != Etotd[j+1])
              Die(FUNCTION_NAME, "Energy %E not stored twice in total xs data "
                                 "for shell=%ld, Z=%s", Etotd[j], i, strZ);
            idxss0 = j+1;
          }

          break;
        }
      }

      if (idxss0 >= idx_tot_max)
        Die(FUNCTION_NAME, "Energy %E not found in total xs data for "
                           "shell=%ld Z=%s", Essd[i][0], i, strZ);

      if(idxss0 - idxtot0 < 0)
        Die(FUNCTION_NAME, "Negative shell location index for shell=%ld Z=%s",
            Essd[i][0], i, strZ);

      /* Store shell location index */
      WDB[ptridxss++] = (double)(idxss0 - idxtot0);

      /* Store binding energy which is assumed to be the first
       * element of the data */
      WDB[ptrEb++] = Essd[i][0];

      /* Shell probability array */
      ptrp2 = ReallocMem(DATA_ARRAY, idx_tot_max - idxss0 + 1);
      WDB[ptrp++] = (double)ptrp2;

      /* Shell log-probability array */
      ptrlp2 = ReallocMem(DATA_ARRAY, idx_tot_max - idxss0 + 1);
      WDB[ptrlp++] = (double)ptrlp2;

      /* Interpolate shell cross sections on total energy gid.
       * NOTE: InterpolateData() is not used due to co-incident energy
       * points. */
      k = 0;

      for (j = idxss0; j <= idx_tot_max; j++) {

        /* Find the index from shell energy array */
        for (; k < Nxsssd[i]-1; k++)
          if (Etotd[j] < Essd[i][k+1])
            break;

        if (k == Nxsssd[i]-1) {

          /* Check the last energy and index */
          if ((Etotd[j] < RDB[DATA_PHOTON_EMAX]) || (j != idx_tot_max))
            Die(FUNCTION_NAME, "Last energy %E below maximum %E", Etotd[j],
                RDB[DATA_PHOTON_EMAX]);

          if (Etotd[j] == Essd[i][k])
            k--;
          else if (Etotd[j] > Essd[i][k])
            Die(FUNCTION_NAME, "Total xs energy %E above maximum shell xs "
                               "energy %E", Etotd[j], Essd[i][k]);
        }

        /* Interpolate shell xs */
        xstmp = ENDFInterp(5, Etotd[j], Essd[i][k], Essd[i][k+1], xsssd[i][k],
                           xsssd[i][k+1]);

        /* Store log-probability */
        WDB[ptrp2++] = xstmp/xstotd[j];
        WDB[ptrlp2++] = log(xstmp/xstotd[j]);

        /* Increase total probability */
        psssum[j] += xstmp/xstotd[j];

#ifdef PRINT_PE_DATA
        fprintf(fp_out, "%E %E\n", Etotd[j], xstmp/xstotd[j]);
#endif

      }
    }

    /* Loop over the energy grid and check total probabilities */
    for (j = 0; j <= idx_tot_max; j++) {

#ifdef PRINT_PE_DATA
      fprintf(fp_out2, "%E %E %E\n", Etotd[j], psssum[j], psssum[j] - 1.0);
#endif

      /* NOTE: This fails due to numerical accuracy when RDB[DATA_PHOTON_EMIN]
       * is small (below 1e-5 MeV or so). */
      if (psssum[j] > 1.0)
        Die(FUNCTION_NAME, "The sum of shell probabilities %.14E is above 1 "
                           "at %E MeV for Z=%s", psssum[j], Etotd[j], strZ);
    }

    Mem(MEM_FREE, psssum);
  }



  /* Free memory */
  if (Nssd > 0) {
    for (i = 0; i < Nssd; i++) {
      Mem(MEM_FREE, Essd[i]);
      Mem(MEM_FREE, xsssd[i]);
    }
    Mem(MEM_FREE, Essd);
    Mem(MEM_FREE, xsssd);
    Mem(MEM_FREE, Nxsssd);
  }
  if (Nxstotd > 0) {
    Mem(MEM_FREE, Etotd);
    Mem(MEM_FREE, xstotd);
  }


#ifdef PRINT_PE_DATA
  fclose(fp_out);
  fclose(fp_out2);
#endif

  /***************************************************************************/

}

/*****************************************************************************/


/*****************************************************************************/

void ProcessPhotoelectricFluorescenceCDF(double *E, long nE) {
  /* Calculates and stores CDFs of shell probabilities at fluorescence
   * energies to enable faster sampling in the function Photoelectric().
   * The probabilities are calculated as in Photoelectric().
   * */
  static char * const FUNCTION_NAME = "ProcessPhotoelectricIncludeAREnergies:";
  long nuc, rea, ptd, ptr, Nssd, i, j, ssmin, idx, NEd, idxmin, idxE, szE,
      ptrpeare, ptrpearcdf, ptrpearcdf2;
  double lE, if0, cdfss;
  const double *Ebd, *idxssd, *lpd, *lpss, *lEarrd;

  if (nE == 0)
    Die(FUNCTION_NAME, "nradtrZ is zero");

  /* Loop over nuclides */
  nuc = (long)RDB[DATA_PTR_NUC0];

  while (nuc > VALID_PTR) {

    /* Only photon nuclides are included */
    if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON) {
      /* Next nuclide */
      nuc = NextItem(nuc);
      continue;
    }

    rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

    while (rea > VALID_PTR) {

      if ((long)RDB[rea + REACTION_MT] != 522) {
        rea = NextItem(rea);
        continue;
      }

      /* Pointer to photon reaction data */
      ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
      CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

      /* Size of the shell probability array */
      Nssd = (long)RDB[ptd + PHOTON_DIST_PE_NSS];

      /* Skip if photoelectric data is not found */
      if (Nssd == 0) {
        rea = NextItem(rea);
        continue;
      }

      /* Transition energies */
      ptrpeare = ReallocMem(DATA_ARRAY, nE);
      WDB[ptd + PHOTON_DIST_PE_AR_E] = (double)ptrpeare;

      /* Shell probabilities */
      ptrpearcdf = ReallocMem(DATA_ARRAY, nE);
      WDB[ptd + PHOTON_DIST_PE_AR_CDF] = (double)ptrpearcdf;


      ptr = (long)RDB[ptd + PHOTON_DIST_PE_EBI];
      CheckPointer(FUNCTION_NAME, "(EBI)", DATA_ARRAY, ptr);
      Ebd = &RDB[ptr];

      ptr = (long)RDB[ptd + PHOTON_DIST_PE_IDXSS];
      CheckPointer(FUNCTION_NAME, "(IDXSS)", DATA_ARRAY, ptr);
      idxssd = &RDB[ptr];

      ptr = (long)RDB[ptd + PHOTON_DIST_PE_LPSS];
      CheckPointer(FUNCTION_NAME, "(LPSS)", DATA_ARRAY, ptr);
      lpd = &RDB[ptr];

      NEd = (long)RDB[ptd + PHOTON_DIST_PE_NE];

      ptr = (long)RDB[ptd + PHOTON_DIST_PE_E];
      CheckPointer(FUNCTION_NAME, "(E)", DATA_ARRAY, ptr);
      lEarrd = &RDB[ptr];



      /* Loop over transition energies */
      for (i = 0; i < nE; i++) {

        WDB[ptrpeare++] = E[i];

        if (E[i] < RDB[DATA_PHOTON_EMIN])
          Die(FUNCTION_NAME, "Photon energy %E below the lower limit %E",
              E[i], RDB[DATA_PHOTON_EMIN]);

        if (E[i] < Ebd[Nssd-1]) {
          WDB[ptrpearcdf++] = NULLPTR;
          continue;
        }

        /* Find the minimum subhell (K=0, L1=1, L2=2, etc) */
        for (ssmin = 0; ssmin < Nssd; ssmin++)
          if (E[i] >= Ebd[ssmin])
            break;

        /* Find energy interval from the energy grid corresponding to the
         *  minimum subshell index */
        idxmin = (long)idxssd[ssmin];

        if (ssmin > 0)
          szE = (long)idxssd[ssmin-1] - idxmin + 1;
        else
          szE = NEd - idxmin;

        lE = log(E[i]);

        /* Avoid numerical problems */
        if (fabs(lE/lEarrd[idxmin] - 1.0) <= 1.0e-14) {
          /* Transition energy is equal to the binding energy */
          idxE = idxmin;
        }
        else {
          idxE = SearchArray(&lEarrd[idxmin], lE, szE) + idxmin;
          CheckValue(FUNCTION_NAME, "lE", "", lE, lEarrd[idxE], lEarrd[idxE+1]);
        }

        if (idxE < idxmin)
          Die(FUNCTION_NAME, "Energy %E not found in the subshell xs data",
              E[i]);

        /* Interpolation factor */
        if0 = (lE - lEarrd[idxE])/(lEarrd[idxE+1] - lEarrd[idxE]);

        /* Allocate memory for shell probabilities */
        ptrpearcdf2 = ReallocMem(DATA_ARRAY, Nssd - ssmin);
        WDB[ptrpearcdf++] = (double)ptrpearcdf2;

        /* Reset variables to avoid dangling pointers */
        Ebd = &RDB[(long)RDB[ptd + PHOTON_DIST_PE_EBI]];
        idxssd = &RDB[(long)RDB[ptd + PHOTON_DIST_PE_IDXSS]];
        lpd = &RDB[(long)RDB[ptd + PHOTON_DIST_PE_LPSS]];
        lEarrd = &RDB[(long)RDB[ptd + PHOTON_DIST_PE_E]];

        cdfss = 0.0;

        /* Loop over shells */
        for (j = ssmin; j < Nssd; j++) {

          /* Log-probability array of the shell */
          lpss = &RDB[(long)lpd[j]];

          /* Probability array index corresponding to the energy index */
          idx = idxE - (long)idxssd[j];

          /* Log-log interpolation of the shell probability */
          cdfss += exp(lpss[idx] + if0*(lpss[idx+1] - lpss[idx]));
          WDB[ptrpearcdf2++] = cdfss;
        }

        if (cdfss > 1.0)
          Die(FUNCTION_NAME, "Total photoelectric probability exceeds 1");
      }

      /* Next reaction */
      rea = NextItem(rea);
    }

    /* Next nuclide */
    nuc = NextItem(nuc);
  }

}

/*****************************************************************************/
