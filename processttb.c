/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processttb.c                                   */
/*                                                                           */
/* Created:       2014/08/11 (TKa)                                           */
/* Last modified: 2017/03/08 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Processes thick-target bremsstrahlung data                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

/* Local function definitions */
static void TTBprint2Struct(FILE *, char *, char *, double **, long, long);


#ifdef TTB_LOGINTEGRALS
static double LogIntegral(const double *, const double *, long, long, long);
static void CumLogIntegral(const double *, const double *, double *, long, long, long);
static void CumLogLogIntegral(const double *x, const double *lx, const double *y, double *F, long n);
#endif
/*****************************************************************************/


/*****************************************************************************/

void ProcessTTB(ElBrSXSData *SXSData, ElSPData *SPData) {

  /* SXSData and SPData contain the material data */

  static char * const FUNCTION_NAME = "ProcessTTB:";
  long i, j, loc0, ptr, ptr1, mat, nE, nkappa, idxSP, idxSXS, nint, idxkj,
      ptrlEk;
  double kappaj, SXSe1, SXSp1, c, tmp;
  double *beta2, *integrande, *integrandp, *Yne, *Ynp, *Eeprint;
  double **pdfBrEe, **pdfBrEp, **cdfBrEe, **cdfBrEp, **coefftmp;
  const double *Ee, *lEe, *Ek, *kappa, *SPtote, *SPtotp;
  const  double **SXSeT, **SXSpT;
  char fname[MAX_STR];
  char matstr[MAX_STR];
  FILE *fp;
  long printttbdata = 0;  /* Use only for debugging */

  fprintf(outp, " - processing TTB data...\n");

  /* Open file for printing the TTB data if needed */
  if (printttbdata) {
    sprintf(fname, "%s_ttbdata.m", GetText(DATA_PTR_INPUT_FNAME));
    if ((fp = fopen(fname, "w")) == NULL)
      Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);   
  }

  /***** Set energy grids ****************************************************/

  /* Set separate pointers for electron and photon energy grids for clarifying
   * purposes. Note that the energy grids must be identical. */

  /* Grid size (used for both photons and electrons) */
  nE = (long)RDB[DATA_ELECTRON_SP_N];

  /* Electron energy grid */
  ptr = (long)RDB[DATA_ELECTRON_SP_E];
  CheckPointer(FUNCTION_NAME, "(LE)", DATA_ARRAY, ptr);
  Ee = &RDB[ptr];

  /* Electron log-energy grid */
  ptr = (long)RDB[DATA_ELECTRON_SP_LE];
  CheckPointer(FUNCTION_NAME, "(LE)", DATA_ARRAY, ptr);
  lEe = &RDB[ptr];

  /* Photon energy grid */
  ptr = (long)RDB[DATA_ELECTRON_SP_E];
  CheckPointer(FUNCTION_NAME, "(LE)", DATA_ARRAY, ptr);
  Ek = &RDB[ptr];

  /* Photon energy per electron energy used in the scaled cross section data */
  nkappa = SXSData->nkappa;
  kappa = (const double *)SXSData->kappa;

  /* Check the energy grid */
  for (i = 0; i < nE; i++) {
    if (Ee[i] != SXSData->E[i])
      Die(FUNCTION_NAME, "TTB energy grid differs from SXSData energy grid");
    if (Ee[i] != SPData->E[i])
      Die(FUNCTION_NAME, "TTB energy grid differs from SPData energy grid");
  }

  /* Find the index corresponding to the energy threshold for linear
   * interpolation of the TTB number yield (linear interpolation used
   * is used below the energy) */
  WDB[DATA_PHOTON_TTBLINEIDX] = -1;
  for (i = 0; i < nE; i++) {
    if (Ee[i] > RDB[DATA_PHOTON_TTBLINE]) {
      WDB[DATA_PHOTON_TTBLINEIDX] = (double)i+1;
      break;
    }
  }

  if (RDB[DATA_PHOTON_TTBLINEIDX] < 0)
    Die(FUNCTION_NAME, "Energy threshold ttbline = %E not found in the energy "
                       "grid", RDB[DATA_PHOTON_TTBLINE]);

  /***************************************************************************/

  /***** Allocate memory *****************************************************/

  beta2 = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
  integrande = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
  integrandp = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
  Yne = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
  Ynp = (double *)Mem(MEM_ALLOC, nE, sizeof(double));

  coefftmp = (double **)Mem(MEM_ALLOC, nE, sizeof(double*));
  pdfBrEe = (double **)Mem(MEM_ALLOC, nE, sizeof(double*));
  pdfBrEp = (double **)Mem(MEM_ALLOC, nE, sizeof(double*));
  cdfBrEe = (double **)Mem(MEM_ALLOC, nE, sizeof(double*));
  cdfBrEp = (double **)Mem(MEM_ALLOC, nE, sizeof(double*));

  for (i = 0; i < nE; i++) {
    beta2[i] = Ee[i]*(Ee[i] + 2.0*E_RESTMASS)
               /((Ee[i] + E_RESTMASS)*(Ee[i] + E_RESTMASS));
    coefftmp[i] = (double *)Mem(MEM_ALLOC, 4, sizeof(double));
    pdfBrEe[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
    pdfBrEp[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
    cdfBrEe[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
    cdfBrEp[i] = (double *)Mem(MEM_ALLOC, nE, sizeof(double));
  }

  /***************************************************************************/

  /***************************************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR) {

    /* Skip parent materials */

    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
      /* Pointer to next */
      mat = NextItem(mat);

      /* Cycle loop */
      continue;
    }

    /* Reset variables to avoid dangling pointers */
    Ee = &RDB[(long)RDB[DATA_ELECTRON_SP_E]];
    lEe = &RDB[(long)RDB[DATA_ELECTRON_SP_LE]];
    Ek = &RDB[(long)RDB[DATA_ELECTRON_SP_E]];

    /* Get data indexes */
    if ((idxSP = ElSPDataGetIdx(SPData, mat)) == -1)
      Die(FUNCTION_NAME, "Material %s not found in SPData",
          GetText(mat + MATERIAL_PTR_NAME));

    if ((idxSXS = ElBrSXSDataGetIdx(SXSData, mat)) == -1)
      Die(FUNCTION_NAME, "Material %s not found in SXSData",
          GetText(mat + MATERIAL_PTR_NAME));

    SXSeT = (const double **)SXSData->SXSeT[idxSXS];
    SXSpT = (const double **)SXSData->SXSpT[idxSXS];
    SPtote = (const double *)SPData->SPtote[idxSP];
    SPtotp = (const double *)SPData->SPtotp[idxSP];

    /* Loop over the photon energy */
    for (i = 0; i < nE-1; i++) {

      nint = 0;
      idxkj = nkappa - 1;

      /* Loop over the electron energy - constructs the integrand for the pdf */
      for (j = i; j < nE; j++) {

        kappaj = Ek[i]/Ee[j];

        /* Find the interval for kappaj */
        while ((kappaj <= kappa[idxkj]) && (--idxkj >= 0));

        /* Check the index */
        if (idxkj < 0)
          Die(FUNCTION_NAME, "idxkj below zero");
        if (idxkj > nkappa-2)
          Die(FUNCTION_NAME, "idxkj above the maximum %ld", kappa[nkappa-2]);

        /* Check kappa */
        if (kappaj < kappa[idxkj])
          Die(FUNCTION_NAME, "kappaj %E below lower limit %E", kappaj, kappa[idxkj]);
        if (kappaj > kappa[idxkj+1])
          Die(FUNCTION_NAME, "kappaj %E above upper limit %E", kappaj, kappa[idxkj+1]);

        /* Interpolate the SXS (lin-lin) */
        SXSe1 = ENDFInterp(2, kappaj, kappa[idxkj], kappa[idxkj+1],
                           SXSeT[idxkj][j], SXSeT[idxkj+1][j]);
        SXSp1 = ENDFInterp(2, kappaj, kappa[idxkj], kappa[idxkj+1],
                           SXSpT[idxkj][j], SXSpT[idxkj+1][j]);

        /* Set the integrand */
        integrande[nint] = SXSe1/(beta2[j]*SPtote[j]*Ek[i]);
        integrandp[nint] = SXSp1/(beta2[j]*SPtotp[j]*Ek[i]);

        nint++;
      }

      /* Integrate! */
      if ((nint > 2) && (nint <= nE)) {
        /* Cubic spline integration on log-energy and linear-integrand
         * scale (mode 3) (NOTE: Interpolation mode 2 (lin-lin) seems to
         * be as good as 3) */
        CSplineCumIntegral0(&Ee[i], integrande, &pdfBrEe[i][i], nint, 0., 0.,
                            3, &lEe[i], coefftmp);
        CSplineCumIntegral0(&Ee[i], integrandp, &pdfBrEp[i][i], nint, 0., 0.,
                            3, &lEe[i], coefftmp);
      }
      else if (nint == 2) {
        /* The last points are integrated using trapezoidal rule on log-log-scale */
        pdfBrEe[i][i+1] = TrapzReal(&Ee[i], integrande, 2, NULL, 5);
        pdfBrEp[i][i+1] = TrapzReal(&Ee[i], integrandp, 2, NULL, 5);
      }
      else {
        Die(FUNCTION_NAME, "Incorrect number of integration points: %ld", nint);
      }
    }

    /* Set transposes */
    for (i = 0; i < nE; i++) {
      for (j = i; j < nE; j++) {
        tmp = pdfBrEe[j][i];
        pdfBrEe[j][i] = pdfBrEe[i][j];
        pdfBrEe[i][j] = tmp;

        tmp = pdfBrEp[j][i];
        pdfBrEp[j][i] = pdfBrEp[i][j];
        pdfBrEp[i][j] = tmp;
      }
    }

    /* The last value of the PDF, which corresponds to Ek=Te, is set to a
     * small non-zero value to enable logarithmic interpolation between all
     * points. */
    for (i = 1; i < nE; i++) {
      pdfBrEe[i][i] = 1.e-6*pdfBrEe[i][i-1];
      pdfBrEp[i][i] = 1.e-6*pdfBrEp[i][i-1];
    }

    /* Integrate the bremsstrahlung CDFs from the PDFs (log-log) */
    /* NOTE: The cdf is zero when i = 0 */
    for (i = 1; i < nE; i++) {
      TrapzRealCum(Ek, pdfBrEe[i], cdfBrEe[i], i+1, 5);
      TrapzRealCum(Ek, pdfBrEp[i], cdfBrEp[i], i+1, 5);
    }

    /* Set the photon number yields */
    for (i = 1; i < nE; i++) {
      Yne[i] = cdfBrEe[i][i];
      Ynp[i] = cdfBrEp[i][i];
    }

    /* Set a small non-zero first value at the cutoff energy */
    Yne[0] = 1.e-6*Yne[1];
    Ynp[0] = 1.e-6*Ynp[1];


    /***** Set data **********************************************************/

    /* Allocate memory for TTB block */
    loc0 = NewItem(mat + MATERIAL_PTR_TTB, TTB_BLOCK_SIZE);

    /* Grid size */
    WDB[loc0 + TTB_NE] = (double)nE;

    /* Energy array pointer */
    WDB[loc0 + TTB_E] = RDB[DATA_ELECTRON_SP_E];

    /* Log energy array pointer */
    WDB[loc0 + TTB_LE] = RDB[DATA_ELECTRON_SP_LE];

    /* Bremsstrahlung pdf for electrons */
    ptr = ReallocMem(DATA_ARRAY, nE);
    WDB[loc0 + TTB_BREPDF] = (double)ptr;
    for (i = 0; i < nE; i++) {
      ptr1 = ReallocMem(DATA_ARRAY, i+1);
      WDB[ptr++] = (double)ptr1;
      for (j = 0; j < i+1; j++)
        WDB[ptr1++] = pdfBrEe[i][j];
    }

    /* Bremsstrahlung cdf for electrons */
    ptr = ReallocMem(DATA_ARRAY, nE);
    WDB[loc0 + TTB_BRECDF] = (double)ptr;
    for (i = 0; i < nE; i++) {
      ptr1 = ReallocMem(DATA_ARRAY, i+1);
      WDB[ptr++] = (double)ptr1;
      for (j = 0; j < i+1; j++)
        WDB[ptr1++] = cdfBrEe[i][j];
    }

    /* Photon number yield for electrons */
    ptr = ReallocMem(DATA_ARRAY, nE);
    WDB[loc0 + TTB_YE] = (double)ptr;
    for (i = 0; i < nE; i++)
      WDB[ptr++] = Yne[i];

    /* log photon number yield for electrons above the first energy */
    ptr = ReallocMem(DATA_ARRAY, nE);
    WDB[loc0 + TTB_LYE] = (double)ptr;
    for (i = 0; i < nE; i++)
      WDB[ptr++] = log(Yne[i]);

    /* Interpolation constant */
    ptr = ReallocMem(DATA_ARRAY, nE);
    WDB[loc0 + TTB_INTCE] = (double)ptr;
    for (i = 0; i < nE; i++) {
      ptr1 = ReallocMem(DATA_ARRAY, i+1);
      WDB[ptr++] = (double)ptr1;
      ptrlEk = (long)RDB[DATA_ELECTRON_SP_LE];
      for (j = 0; j < i; j++) {
        c = log(pdfBrEe[i][j+1]/pdfBrEe[i][j])
            / (RDB[ptrlEk + j + 1] - RDB[ptrlEk + j]);

        if (c == -1)
          Die(FUNCTION_NAME, "Interpolation constant c = -1");
        if (isnan(c))
          Die(FUNCTION_NAME, "Interpolation constant c = NaN");

        WDB[ptr1++] = c;
      }
    }

    /* Set separate bremsstrahlung data for positrons */
    if ((long)RDB[DATA_PHOTON_TTBPM] == YES) {

      /* Bremsstrahlung pdf for positrons */
      ptr = ReallocMem(DATA_ARRAY, nE);
      WDB[loc0 + TTB_BRPPDF] = (double)ptr;
      for (i = 0; i < nE; i++) {
        ptr1 = ReallocMem(DATA_ARRAY, i+1);
        WDB[ptr++] = (double)ptr1;
        for (j = 0; j < i+1; j++)
          WDB[ptr1++] = pdfBrEp[i][j];
      }

      /* Bremsstrahlung cdf for positrons */
      ptr = ReallocMem(DATA_ARRAY, nE);
      WDB[loc0 + TTB_BRPCDF] = (double)ptr;
      for (i = 0; i < nE; i++) {
        ptr1 = ReallocMem(DATA_ARRAY, i+1);
        WDB[ptr++] = (double)ptr1;
        for (j = 0; j < i+1; j++)
          WDB[ptr1++] = cdfBrEp[i][j];
      }

      /* Photon number yield for positrons */
      ptr = ReallocMem(DATA_ARRAY, nE);
      WDB[loc0 + TTB_YP] = (double)ptr;
      for (i = 0; i < nE; i++)
        WDB[ptr++] = Ynp[i];

      /* log photon number yield for positrons */
      ptr = ReallocMem(DATA_ARRAY, nE);
      WDB[loc0 + TTB_LYP] = (double)ptr;
      for (i = 0; i < nE; i++)
        WDB[ptr++] = log(Ynp[i]);


      /* Interpolation constant */
      ptr = ReallocMem(DATA_ARRAY, nE);
      WDB[loc0 + TTB_INTCP] = (double)ptr;
      for (i = 0; i < nE; i++) {
        ptr1 = ReallocMem(DATA_ARRAY, i+1);
        WDB[ptr++] = (double)ptr1;
        ptrlEk = (long)RDB[DATA_ELECTRON_SP_LE];
        for (j = 0; j < i; j++) {
          c = log(pdfBrEp[i][j+1]/pdfBrEp[i][j])
              / (RDB[ptrlEk + j + 1] - RDB[ptrlEk + j]);

          if (c == -1)
            Die(FUNCTION_NAME, "Interpolation constant c = -1");
          if (isnan(c))
            Die(FUNCTION_NAME, "Interpolation constant c = NaN");

          WDB[ptr1++] = c;
        }
      }
    }

    /*************************************************************************/

    /* Print TTB data */
    if (printttbdata) {
      sprintf(matstr, "%s", GetText(mat + MATERIAL_PTR_NAME));
      Eeprint = (double *)Ee;
      TTBprint2Struct(fp, matstr, "E", &Eeprint, 1, nE);
      TTBprint2Struct(fp, matstr, "Yne", &Yne, 1, nE);
      TTBprint2Struct(fp, matstr, "Ynp", &Ynp, 1, nE);
      TTBprint2Struct(fp, matstr, "pdfBrEe", pdfBrEe, nE, nE);
      TTBprint2Struct(fp, matstr, "cdfBrEe", cdfBrEe, nE, nE);
      TTBprint2Struct(fp, matstr, "pdfBrEp", pdfBrEp, nE, nE);
      TTBprint2Struct(fp, matstr, "cdfBrEp", cdfBrEp, nE, nE);
    }


    /* Next material */
    mat = NextItem(mat);
  }

  /***************************************************************************/

  /* Close ttb data file */
  if (printttbdata)
    fclose(fp);


  /* Free memory */
  for (i = 0; i < nE; i++) {
    Mem(MEM_FREE, pdfBrEe[i]);
    Mem(MEM_FREE, cdfBrEe[i]);
    Mem(MEM_FREE, pdfBrEp[i]);
    Mem(MEM_FREE, cdfBrEp[i]);
    Mem(MEM_FREE, coefftmp[i]);
  }

  Mem(MEM_FREE, pdfBrEe);
  Mem(MEM_FREE, cdfBrEe);
  Mem(MEM_FREE, pdfBrEp);
  Mem(MEM_FREE, cdfBrEp);
  Mem(MEM_FREE, coefftmp);
  Mem(MEM_FREE, integrande);
  Mem(MEM_FREE, integrandp);
  Mem(MEM_FREE, Yne);
  Mem(MEM_FREE, Ynp);
  Mem(MEM_FREE, beta2);

}
/*****************************************************************************/


/*****************************************************************************/

static void TTBprint2Struct(FILE *fp, char *strMat, char *strVar,
                            double **A, long n, long m) {
  /* Prints data into Matlab struct.
   * */
  static char * const FUNCTION_NAME = "TTBprint2Struct:";
  long i, j;
  char *strStruct = "ttb_data";

  if (n < 1 || m < 1)
    Die(FUNCTION_NAME, "n < 1 or m < 1");
  if (m == 1)
    Die(FUNCTION_NAME, "n = 1 not supported");

  if (n == 1) {
    fprintf(fp, "%s.%s.%s = [ ", strStruct, strMat, strVar);
    for (j = 0; j < m; j++)
      fprintf(fp, "%.6E ", (*A)[j]);
    fprintf(fp, "];\n");
  }
  else {
    fprintf(fp, "%s.%s.%s = [...\n", strStruct, strMat, strVar);
    for (i = 0; i < n; i++) {
      fprintf(fp, " [ ");
      for (j = 0; j < m; j++)
        fprintf(fp, "%.6E ", A[i][j]);
      fprintf(fp, "]; ...\n");
    }
    fprintf(fp, " ];\n");
  }

}

/*****************************************************************************/


#ifdef TTB_LOGINTEGRALS
/*****************************************************************************/

static double LogIntegral(const double *x, const double *y, long N, long mode,
                          long islogx) {
  /* Helper function */
  static char * const FUNCTION_NAME = "LogIntegral:";
  long ninteg, i;
  double x1, x2, y1, y2, lxr, c, F;

  if (N < 2)
    Die(FUNCTION_NAME, "LogIntegral: Number of elements less than two");

  ninteg = N - 1;
  F = 0.;

  if (mode == 0) {
    /* Integral of log-log interpolation on linear scale */
    for (i = 0; i < ninteg; i++) {
      x1 = x[i];
      x2 = x[i+1];
      y1 = y[i];
      y2 = y[i+1];


//      if (x1 == x2)
//        Die(FUNCTION_NAME, "x1 == x2");
//      if (x1 == 0.)
//        Die(FUNCTION_NAME, "x1 == 0");
//      if (x2 == 0.)
//        Die(FUNCTION_NAME, "x2 == 0");
//      if (y1 == 0.)
//        Die(FUNCTION_NAME, "y1 == 0");
//      if (y2 == 0.)
//        Die(FUNCTION_NAME, "y2 == 0");

      if (islogx)
        lxr = x2 - x1;
      else
        lxr = log(x2/x1);

      c = log(y2/y1)/lxr;

      if (c == -1) {
        F += y1*x1*lxr;
      }
      else {
//        F = F + y1/(a + 1.0)*(pow(x2/x1, a)*x2 - x1);
//        double tmp1 = y1/(a + 1.0)*(pow(x2/x1, a)*x2 - x1);
//        double tmp2 = y1/(a + 1.0)*(exp(a*log(x2/x1))*x2 - x1);
//        double tmp3 = y1/(a + 1.0)*(y2/y1*x2 - x1);
//        printf(" %.12E %.12E %.12E\n", tmp1, tmp2, tmp3);
//        F += y1/(c + 1.0)*(exp(c*log(x2/x1))*x2 - x1); /* TODO: Also to ProcessRayleigh()? */
        F += (y2*x2 - y1*x1)/(c + 1.0); /* TODO: Also to ProcessRayleigh()? */
      }
    }
  }
  else if (mode == 1) {
    /* Integral of log-lin interpolation on linear scale */
    for (i = 0; i < ninteg; i++) {
      x1 = x[i];
      x2 = x[i+1];
      y1 = y[i];
      y2 = y[i+1];

      if (x1 == x2)
        Die(FUNCTION_NAME, "x1 == x2");
      if (x1 == 0)
        Die(FUNCTION_NAME, "x1 == 0");
      if (x2 == 0)
        Die(FUNCTION_NAME, "x2 == 0");

      F += y1*(x2 - x1) + (y2 - y1)*(x2 + (x1 - x2)/(log(x2/x1)));
    }
  }
  else
    Die(FUNCTION_NAME, "unknown rule");

  return F;

}

/*****************************************************************************/


/*****************************************************************************/

static void CumLogLogIntegral(const double *x, const double *lx, const double *y,
                            double *F, long n) {
  /* TODO: description */

  static char * const FUNCTION_NAME = "CumLogLogIntegral:";
  long i;
  double c, lxr;

  F[0] = 0;

  for (i = 1; i < n; i++) {

    /* Check arrays */
    if (x[i-1] == x[i])
      Die(FUNCTION_NAME, "x1 == x2");
    if (x[i-1] == 0.)
      Die(FUNCTION_NAME, "x1 == 0");
    if (x[i] == 0.)
      Die(FUNCTION_NAME, "x2 == 0");
    if (y[i-1] == 0.)
      Die(FUNCTION_NAME, "y1 == 0");
    if (y[i] == 0.)
      Die(FUNCTION_NAME, "y2 == 0");

    lxr = lx[i] - lx[i-1];
    c = log(y[i]/y[i-1])/lxr;

    if (c == -1)
      F[i]= F[i-1] + y[i-1]*x[i-1]*lxr;
    else
      F[i]= F[i-1] + (y[i]*x[i] - y[i-1]*x[i-1])/(c + 1.0); /* TODO: Also to ProcessRayleigh()? */
  }
}

/*****************************************************************************/


/*****************************************************************************/

//static void CumLogIntegral(const double *x, const double *y, double *F,
//                             long N, long mode, long islogx) {
//  /* Helper function */

//  long i;
////  double xtmp[2];
////  double ytmp[2];

////  if (N < 2)
////    Die(FUNCTION_NAME, "CumLogIntegral: Number of elements less than two");

//  F[0] = 0;

//  for (i = 1; i < N; i++) {
////    xtmp[0] = x[i-1];
////    xtmp[1] = x[i];
////    ytmp[0] = y[i-1];
////    ytmp[1] = y[i];
////    F[i] = F[i-1] + LogIntegral(xtmp, ytmp, 2, mode, islogx);
//    F[i] = F[i-1] + LogIntegral(&x[i-1], &y[i-1], 2, mode, islogx);
//  }

//}

/*****************************************************************************/
#endif
