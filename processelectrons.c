/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processelectrons.c                             */
/*                                                                           */
/* Created:       2016/06/07 (TKa)                                           */
/* Last modified: 2017/03/10 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Processes electron and positron interaction data             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"


/* Local function definitions */
static void ProcessElectronsInput();
static void IntegrateSPrad(ElBrSXSData *, long, double *, double, long);
static double ElIcompESTAR(const long *, const double *, const double *, long,
                           long);
static void ElDeffcorSternheimer(double, double, const double *, const double *,
                                 long, const double *, long, double *);
static double PositronRp(double, double);
static void ProcessSPcolBethe(const double *, long, double, double, double *,
                              double *, double *);
static void RadiativeYield(const double *, const double *,
                           const double *, long,  double *);
static void PrintSPDataMat(FILE *, ElSPData *);


#ifdef TEST_EL_FUTURE
static void BrMFP(ElBrSXSData *SXSData, long idx, double Ecut, double multiplier,
                  double *mfp);
static void RangeCSDA(double *E, double *Stot, double *R, long n);
static double LogIntegral(const double *x, const double *y, long N, long mode);
#endif

/*****************************************************************************/

void ProcessElectrons() {

  static char * const FUNCTION_NAME = "ProcessElectrons:";
  long mat, iso, nuc, loc0, elptr, Zidx, Z, i, j, idxSXS, idxSP, idxDC, iZ,
       nmat, idxZ, ssidx, nZ, nsstot, idxMat, negocc, nE, dieflag, ptr;
  double eldens, adensZ2, Z2eff, Emin, Emax, Rp;
  long *Zmat, *Zset;
  double *E, *IZ, *ebi, *fel, *adens, *mfracmdens;
  ElBrSXSData SXSData;
  ElSPData SPData;
  GSconfigData gsconfigData;
  ElBrSXSData SXSDataZ;
  ElSPData SPDataZ;
  ElBrSXSData SXSDataMat;
  ElSPData SPDataMat;
  char fname[MAX_STR];
  FILE *fp;


  fprintf(outp, "Processing electron data:\n\n");


  /***** Allocate memory for electron block for all the materials ************/
  /* TODO: Tähän joku if(elmat)? */

  /* Loop over materials */
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR) {

    /* Skip parent materials */
    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
      mat = NextItem(mat);
      continue;
    }

    NewItem(mat + MATERIAL_PTR_EL, EL_BLOCK_SIZE);

    /* Pointer to electron data */
    elptr = (long)RDB[mat + MATERIAL_PTR_EL];
    CheckPointer(FUNCTION_NAME, "(elptr)", DATA_ARRAY, elptr);

    /* Set default values */

    /* Mean excitation energy: negative value -> calculated here */
    WDB[elptr + EL_MEE] = -1.0;

    /* Material conductivity: 2 -> specified here */
    WDB[elptr + EL_COND] = 2;

    /* Material phase : condensed (affects only mixtures) */
    WDB[elptr + EL_GAS] = NO;

    mat = NextItem(mat);
  }

  /***************************************************************************/


  /***** Read electron data **************************************************/

  /* Scaled bremsstrahlung cross section data */
  ElBrSXSDataRead(&SXSData);

  /* Stopping power data */
  if (!(long)RDB[DATA_ELECTRON_SP_RAD_CALC])
    ElSPDataRead(&SPData);

  /* Ground state configuration data */
  GSconfigDataRead(&gsconfigData);

  /***************************************************************************/


  /* Process input options */
  ProcessElectronsInput();


  /***** Create electron energy grid *****************************************/

  dieflag = 0;

  if (RDB[DATA_ELECTRON_SP_E] > VALID_PTR) {

    /* User-defined energy grid used for debugging (program terminates after
     * printing stopping powers) */
    dieflag = 1;

    nE = (long)RDB[DATA_ELECTRON_SP_N];
    E = (double *)Mem(MEM_ALLOC, nE, sizeof(double));

    /* Check the energy grid is in ascending orderd */
    E[0] = RDB[(long)RDB[DATA_ELECTRON_SP_E]];
    for (i = 1; i < nE; i++) {
      E[i] = RDB[(long)RDB[DATA_ELECTRON_SP_E] + i];
      if (E[i] <= E[i-1])
        Die(FUNCTION_NAME, "Stopping power energy grid not in ascending order");
    }
  }
  else {
    /* Set logarithmic energy grid.
     * NOTE:
     * - TTB model assumes that the energy array is uniformly distributed
     *   in log. Electron/positron energy search must be changed if
     *   uniform log is not used.
     * - The global minimum RDB[DATA_PHOTON_EMIN] is used instead of a
     *   material-dependent value which would be too complicated.
     */
    nE = (long)RDB[DATA_ELECTRON_SP_N];
    Emin = RDB[DATA_PHOTON_EMIN];
    Emax = RDB[DATA_PHOTON_EMAX];
    E = MakeArray(Emin, Emax, nE, 2);

    /* Store the energy grid */
    loc0 = ReallocMem(DATA_ARRAY, nE);
    WDB[DATA_ELECTRON_SP_E] = (double)loc0;
    for (i = 0; i < nE; i++)
      WDB[loc0++] = E[i];
  }

  /* Store log energy grid */
  loc0 = ReallocMem(DATA_ARRAY, nE);
  WDB[DATA_ELECTRON_SP_LE] = (double)loc0;
  for (i = 0; i < nE; i++)
    WDB[loc0++] = log(E[i]);

  /***************************************************************************/


  /***** Process element data ************************************************/

  fprintf(outp, " - processing element data...\n");

  /* Set structs for element-wise data */

  /* Scaled cross sections */
  ElBrSXSDataAlloc(&SXSDataZ, SXSData.sz, nE, SXSData.nkappa);
  memcpy(SXSDataZ.kappa, SXSData.kappa, SXSData.nkappa*sizeof(double));
  memcpy(SXSDataZ.E, E, nE*sizeof(double));

  /* Stopping powers */
  ElSPDataAlloc(&SPDataZ, PHOTON_ZMAX, nE);
  memcpy(SPDataZ.E, E, nE*sizeof(double));


  /* Loop over nuclides */
  idxZ = 0;
  nuc = (long)RDB[DATA_PTR_NUC0];

  while (nuc > VALID_PTR) {

    /* Only photon nuclides are included */
    if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON) {
      /* Next nuclide */
      nuc = NextItem(nuc);
      continue;
    }

    Z = (long)RDB[nuc + NUCLIDE_Z];

    /* Get SXS data index */
    if ((idxSXS = ElBrSXSDataGetIdx(&SXSData, Z)) == -1)
      Die(FUNCTION_NAME, "SXSData: Z=%ld not found", Z);

    /* Set mapping for the data */
    SXSDataZ.map[idxZ] = Z;
    SPDataZ.map[idxZ] = Z;

    /* Interpolate scaled cross section */
    for (i = 0; i < SXSData.nkappa; i++) {

      /* Interpolate the transpose (log(E),linear SXS) */
      CSplineInterpolate0(SXSData.E, SXSData.SXSeT[idxSXS][i], SXSData.nE,
                          0.0, 0.0, SXSDataZ.E, SXSDataZ.SXSeT[idxZ][i],
                          SXSDataZ.nE, 3);

      for (j = 0; j < SXSDataZ.nE; j++)
        SXSDataZ.SXSe[idxZ][j][i] = SXSDataZ.SXSeT[idxZ][i][j];
    }

    /* Compute radiative stopping powers */
    if ((long)RDB[DATA_ELECTRON_SP_RAD_CALC]) {
      /* Integrate from the SXS data */
      IntegrateSPrad(&SXSDataZ, idxZ, SPDataZ.SPrade[idxZ],
                     (double)(Z*Z)*N_AVOGADRO/RDB[nuc + NUCLIDE_AW], 2);
    }
    else {
      /* Get SP data index */
      if ((idxSP = ElSPDataGetIdx(&SPData, Z)) == -1)
        Die(FUNCTION_NAME, "SPData: Z=%ld not found", Z);

      /* Cubic spline interpolation (log-log) of the SP data */
      CSplineInterpolate0(SPData.E, SPData.SPrade[idxSP], SPData.nE, 0., 0.,
                          SPDataZ.E, SPDataZ.SPrade[idxZ], SPDataZ.nE, 5);
    }

    /* Next nuclide */
    nuc = NextItem(nuc);
    idxZ++;
  }

  /***************************************************************************/


  /***** Process material data ***********************************************/

  fprintf(outp, " - processing material data...\n");

  /* Calculate the total number of materials */
  nmat = 0;
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR) {
    /* Skip parent materials */
    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
      mat = NextItem(mat);
      continue;
    }
    nmat++;
    mat = NextItem(mat);
  }

  if (nmat == 0)
    Die(FUNCTION_NAME, "No materials found");


  /* Structs for material data */
  ElBrSXSDataAlloc(&SXSDataMat, nmat, SXSDataZ.nE, SXSDataZ.nkappa);
  memcpy(SXSDataMat.kappa, SXSDataZ.kappa, SXSDataZ.nkappa*sizeof(double));
  memcpy(SXSDataMat.E, SXSDataZ.E, SXSDataZ.nE*sizeof(double));

  ElSPDataAlloc(&SPDataMat, nmat, SPDataZ.nE);
  memcpy(SPDataMat.E, SPDataZ.E, SPDataZ.nE*sizeof(double));

  /* Flag array for elements */
  Zset = (long *)Mem(MEM_ALLOC, PHOTON_ZMAX, sizeof(long));

  /* Atomic numbers found in material */
  Zmat = (long *)Mem(MEM_ALLOC, PHOTON_ZMAX, sizeof(long));

  /* Atomic densities of elements found in material */
  adens = (double *)Mem(MEM_ALLOC, PHOTON_ZMAX, sizeof(double));

  /* Mass fraction multiplied by mass density for elements found in material */
  mfracmdens = (double *)Mem(MEM_ALLOC, PHOTON_ZMAX, sizeof(double));


  /* Loop over materials */

  idxMat = 0;
  mat = (long)RDB[DATA_PTR_M0];

  while (mat > VALID_PTR) {

    /* Skip parent materials */
    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
      mat = NextItem(mat);
      continue;
    }

    /* Pointer to electron data */
    elptr = (long)RDB[mat + MATERIAL_PTR_EL];
    CheckPointer(FUNCTION_NAME, "(elptr)", DATA_ARRAY, elptr);
    
    /* Initialize arrays */
    for (i = 0; i < PHOTON_ZMAX; i++)
      Zset[i] = -1;

    memset(Zmat, 0, PHOTON_ZMAX*sizeof(long));
    memset(adens, 0, PHOTON_ZMAX*sizeof(double));
    memset(mfracmdens, 0, PHOTON_ZMAX*sizeof(double));

    /* Set mapping for the material data */
    SXSDataMat.map[idxMat] = mat;
    SPDataMat.map[idxMat] = mat;

    nZ = 0;         /* number of elements in the material  */
    nsstot = 0;     /* total number of shells in the material */
    eldens = 0.0;   /* electron density in the material */
    negocc = 0;     /* flag for negative occupation numbers */
    Z2eff = 0.0;    /* "effective" Z squared */


    /***** Loop over material composition and set elemental data *************/

    iso = (long)RDB[mat + MATERIAL_PTR_COMP];
    CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

    while (iso > VALID_PTR) {

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* 4.3.2017 / 2.1.29 / JLe: palamalaskun materiaaleissa listan 1.   */
      /* nuklidia käytetään tallentamaan informaatiota hävinneestä atomi- */
      /* tiheydestä. Se pitää hypätä yli */

      if ((long)RDB[nuc + NUCLIDE_ZAI] == -1) {
        /* Skip nuclide */
        iso = NextItem(iso);
        continue;
      }

      /* JLe (23.2.2017 / 2.2.29): tää on lisätty kytkettyä moodia varten. */
      /* materiaalikoostumukset on annettu neutronivaikutusaloina, ja      */
      /* fotonidataan pitää mennä erillisten pointterien kautta.           */
      /* HUOM: jos materiaalissa on saman alkuaineen kahta isotooppia niin */
      /* sitä vastaava fotonidata tulee käsiteltyä kahteen kertaan. */

      if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_DATA]) > VALID_PTR)
        nuc = ptr;

      /* Only photon nuclides are included */
      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON) {
        /* Next */
        iso = NextItem(iso);
        continue;
      }

      Z = (long)RDB[nuc + NUCLIDE_Z];

      /* Check if the element is new */
      if (Zset[Z-1] == -1) {
        Zmat[nZ] = Z;
        Zset[Z-1] = nZ;
        nZ++;

        /* Check the number of elements */
        if (nZ > PHOTON_ZMAX)
          Die(FUNCTION_NAME,  "Number of elements in material %s is above the "
              "maximum %ld",  GetText(mat + MATERIAL_PTR_NAME), PHOTON_ZMAX);
      }
      
      /* Index for Z */
      Zidx = Zset[Z-1];
      
      /* Atomic density */ 
      adens[Zidx] += RDB[iso + COMPOSITION_ADENS];
        
      /* Mass fraction multiplied by mass density */
      mfracmdens[Zidx] += RDB[nuc + NUCLIDE_AW]*RDB[iso + COMPOSITION_ADENS]
                          /N_AVOGADRO;

      /* Next */
      iso = NextItem(iso);
    }
    
    /* Check the number of elements */
    if (nZ == 0)
      Die(FUNCTION_NAME, "No elements found in material %s",
          GetText(mat + MATERIAL_PTR_NAME));

    if (nZ > PHOTON_ZMAX)
      Die(FUNCTION_NAME, "Number of elements %ld above the maximum %ld in "
          "material %s", nZ, PHOTON_ZMAX, GetText(mat + MATERIAL_PTR_NAME));

    /*************************************************************************/
    
    
    /***** Loop over the elements found in the material **********************/

    for (iZ = 0; iZ < nZ; iZ++) {

      Z = Zmat[iZ];

      /* Get data indexes */
      if ((idxSXS = ElBrSXSDataGetIdx(&SXSDataZ, Z)) == -1)
        Die(FUNCTION_NAME, "SXSDataZ: Z=%ld not found", Z);
      if ((idxSP = ElSPDataGetIdx(&SPDataZ, Z)) == -1)
        Die(FUNCTION_NAME, "SPDataZ: Z=%ld not found", Z);
      if ((idxDC = GSconfigDataGetIdx(&gsconfigData, Z)) == -1)
        Die(FUNCTION_NAME, "gsconfigData: Z=%ld not found", Z);

      /* Set variables */
      eldens += (double)Z*adens[iZ];
      nsstot += gsconfigData.elemData[idxDC].nss;
      adensZ2 = adens[iZ]*(double)(Z*Z);
      Z2eff += adensZ2;

      /* Set electron data */
      for (i = 0; i < SPDataZ.nE; i++) {

        /* NOTE: Assuming equal energy grid in SPDataZ and SXSDataZ */

        /* Radiative stopping power (additivity approximation) */
        SPDataMat.SPrade[idxMat][i] += mfracmdens[iZ]*SPDataZ.SPrade[idxSP][i];

        /* Scaled bremsstrahlung data (additivity approximation) */
        for (j = 0; j < SXSDataZ.nkappa; j++)
          SXSDataMat.SXSe[idxMat][i][j] += adensZ2*SXSDataZ.SXSe[idxSXS][i][j];
      }

      /* Check if negative occupation number exist */
      if (negocc == 0) {
        for (j = 0; j < gsconfigData.elemData[idxDC].nss; j++) {
          if (gsconfigData.elemData[idxDC].gsconf[j] < 0) {
            negocc = 1;
            break;
          }
        }
      }
      
    }

    /* Check the number of shells */
    if (nsstot == 0)
      Die(FUNCTION_NAME, "No shells found in material %s",
          GetText(mat + MATERIAL_PTR_NAME));

    if (nsstot > PHOTON_NSS_MAX*PHOTON_ZMAX)
      Die(FUNCTION_NAME, "Number of shells %ld above the maximum %ld in "
          "material %s", nsstot, PHOTON_NSS_MAX*PHOTON_ZMAX,
          GetText(mat + MATERIAL_PTR_NAME));

    /*************************************************************************/


    /* Set positron data */
    Z2eff /= RDB[mat + MATERIAL_ADENS];  /* effective Z*Z */

    for (i = 0; i < SPDataZ.nE; i++) {
      Rp = PositronRp(SPDataZ.E[i], Z2eff);
      SPDataMat.SPradp[idxMat][i] = Rp*SPDataMat.SPrade[idxMat][i];

      for (j = 0; j < SXSDataZ.nkappa; j++) {
        SXSDataMat.SXSeT[idxMat][j][i] = SXSDataMat.SXSe[idxMat][i][j];
        SXSDataMat.SXSp[idxMat][i][j] = Rp*SXSDataMat.SXSe[idxMat][i][j];
        SXSDataMat.SXSpT[idxMat][j][i] = SXSDataMat.SXSp[idxMat][i][j];
      }
    }


    /* Correct the material phase if needed */
    /* TODO: Warn if both "elmee" and "elgas" card are used for this
     * material, because "elgas" won't have any effect. */
    if (((long)RDB[elptr + EL_GAS] == 1) && (nZ == 1)) {
      Warn(FUNCTION_NAME,
           "Overriding the gas phase of the material %s set by the user:\n"
           "In Serpent, the gas phase doesn't affect the mean excitation energy of a\n"
           "single element material. The material is set to be a \"non-gas\".",
           GetText(mat + MATERIAL_PTR_NAME));
      WDB[elptr + EL_GAS] = 0;
    }


    /***** Set the conductive state of the material **************************/

    switch ((long)RDB[elptr + EL_COND]) {
      case 0: {
        /* Material is a non-conductor - no actions needed */
        break;
      }
      case 1: {
        /* Check if the conductive state if applicable */
        if (negocc == 0) {
          Warn(FUNCTION_NAME,
               "Overriding the conductive state of the material %s set by the user:\n"
               "No conduction electrons found, the material is set to be a non-conductor.",
               GetText(mat + MATERIAL_PTR_NAME));
          WDB[elptr + EL_COND] = 0;
        }
        break;
      }
      case 2: {
        /* Default behaviour, if the user has not set the conductive state: */
        /*  - a single element material is a conductor if conduction electrons
         *    were found, otherwise the material is a non-conductor
         *  - a multielement material is a non-conductor */
        if (nZ == 1) {
          if (negocc == 1)
            WDB[elptr + EL_COND] = 1;
          else
            WDB[elptr + EL_COND] = 0;
        }
        else {
          WDB[elptr + EL_COND] = 0;
        }

        break;
      }
      default: {
        Die(FUNCTION_NAME, "Unsupported value for conductivity.");
        break;
      }
    }

    /*************************************************************************/


    /***** Loop over the elements and set material data **********************/

    /* Allocate memory */
    IZ = (double *)Mem(MEM_ALLOC, nZ, sizeof(double));
    ebi = (double *)Mem(MEM_ALLOC, nsstot, sizeof(double));
    fel = (double *)Mem(MEM_ALLOC, nsstot, sizeof(double));

    ssidx = 0;
    for (i = 0; i < nZ; i++) {

      Zidx = Zmat[i] - 1;

      /* Check that Z is correct in the gsconfigData */
      if (Zmat[i] != gsconfigData.elemData[Zidx].Z)
        Die(FUNCTION_NAME, "Incorrect atomic number Z=%ld in data struct",
            gsconfigData.elemData[Zidx].Z);

      /* Mean excitation energies */
      IZ[i] = gsconfigData.elemData[Zidx].I;

      /* Loop over shells */
      for (j = 0; j < gsconfigData.elemData[Zidx].nss; j++) {

        /* Binding energy */
        ebi[ssidx] = gsconfigData.elemData[Zidx].ebi[j];

        /* Electron fraction:
         *  - if the material is a conductor, a negative value indicates a
         *    conduction electron*/
        if ((long)RDB[elptr + EL_COND] == YES)
          fel[ssidx] = (double)gsconfigData.elemData[Zidx].gsconf[j]
                        *adens[i]/eldens;
        else
          fel[ssidx] = fabs((double)gsconfigData.elemData[Zidx].gsconf[j])
                        *adens[i]/eldens;

        ssidx++;
      }
    }

    /*************************************************************************/


    /***** Set the mean excitation energy ***********************************/

    if (RDB[elptr + EL_MEE] < 0.0) {
      if (nZ > 1) {
        /* Calculate the mean excitation energy for a compound */
        WDB[elptr + EL_MEE] = ElIcompESTAR(Zmat, adens, IZ, nZ,
                                           (long)RDB[elptr + EL_GAS]);
      }
      else {
        /* Single-element material */
        WDB[elptr + EL_MEE] = IZ[0];
      }
    }
    /*************************************************************************/


    /* Calculate the density effect correction using the Sternheimer's
     * method */
    ElDeffcorSternheimer(eldens, RDB[elptr + EL_MEE], fel, ebi, nsstot,
                         SPDataMat.E, SPDataMat.nE, SPDataMat.delta[idxMat]);

    /* Calculate collision stopping powers using Bethe's formula */
    ProcessSPcolBethe(SPDataMat.E, SPDataMat.nE, WDB[elptr + EL_MEE], eldens,
                      SPDataMat.delta[idxMat], SPDataMat.SPcole[idxMat],
                      SPDataMat.SPcolp[idxMat]);


    /* Set the total stopping power */
    for (i = 0; i < SPDataMat.nE; i++) {
      SPDataMat.SPtote[idxMat][i] = SPDataMat.SPrade[idxMat][i]
                                    + SPDataMat.SPcole[idxMat][i];
      SPDataMat.SPtotp[idxMat][i] = SPDataMat.SPradp[idxMat][i]
                                    + SPDataMat.SPcolp[idxMat][i];
    }

    /* Free memory */
    Mem(MEM_FREE, IZ);
    Mem(MEM_FREE, ebi);
    Mem(MEM_FREE, fel);

    /* Next */
    mat = NextItem(mat);
    idxMat++;
  }

  /***************************************************************************/


  /* Process TTB data */
  if ((long)RDB[DATA_PHOTON_USE_TTB] == YES)
    ProcessTTB(&SXSDataMat, &SPDataMat);


  /* Print stopping power data */
  if ((long)RDB[DATA_ELECTRON_PRINT_SP] == YES) {

    /* Set file name */
    sprintf(fname, "%s_elsp.m", GetText(DATA_PTR_INPUT_FNAME));

    /* Open the file */
    if ((fp = fopen(fname, "w")) == NULL)
      Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

    /* Print */
    PrintSPDataMat(fp, &SPDataMat);

    /* Close the file */
    fclose(fp);
  }


  /* Free memory */
  Mem(MEM_FREE, E);
  Mem(MEM_FREE, Zset);
  Mem(MEM_FREE, Zmat);
  Mem(MEM_FREE, adens);
  Mem(MEM_FREE, mfracmdens);

  if (!(long)RDB[DATA_ELECTRON_SP_RAD_CALC])
    ElSPDataFree(&SPData);

  ElBrSXSDataFree(&SXSData);
  GSconfigDataFree(&gsconfigData);
  ElSPDataFree(&SPDataZ);
  ElBrSXSDataFree(&SXSDataZ);
  ElSPDataFree(&SPDataMat);
  ElBrSXSDataFree(&SXSDataMat);

  /* Terminate the program if the energy grid was set by the user */
  if (dieflag)
    Die(FUNCTION_NAME, "User-defined energy grid (set elspe) is used only for "
                       "debugging.");

  fprintf(outp, "OK.\n\n");


#ifdef TTB_UNIT_TEST
  /* TODO: Omaan tiedostoon, ja TTBenergy fiksusti mukaan */

  /* Set random number seed */

  time_t tseed;
  time(&tseed);
  parent_seed = (unsigned long)tseed;
  long id = 0;

  /* Initialize RNG */
  unsigned long seed = ReInitRNG(0);
  SEED[id*RNG_SZ] = seed;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];

  while (mat > VALID_PTR) {

    /* Skip parent materials */
    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
      mat = NextItem(mat);
      continue;
    }

    long idxout;
    double lEe;
    long nbr;
    long ptd, ptr, nEe, idx, part;
    double Yn, Edep;
    const double *Eed, *lEed;

    /* Pointer to TTB data */
    ptd = (long)RDB[mat + MATERIAL_PTR_TTB];
    CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

    /* Energy grid size */
    nEe = (long)RDB[ptd + TTB_NE];

    /* Electron energy grid */
    ptr = (long)RDB[ptd + TTB_E];
    CheckPointer(FUNCTION_NAME, "(E)", DATA_ARRAY, ptr);
    Eed = &RDB[ptr];

    /* Logarithmic energy grid */
    ptr = (long)RDB[ptd + TTB_LE];
    CheckPointer(FUNCTION_NAME, "(LE)", DATA_ARRAY, ptr);
    lEed = &RDB[ptr];

    /* Loop over  */

    for (i = 0; i < nEe; i++) {

      /* Without log-energy (log energy calculated in TTByield()) */

      /* Test electron number yield */
      TTByield(mat, Eed[i], 0, &lEe, MT_ELECTRON_PP_EL, &nbr, &idxout, id);

      /* Force TTB photon production */
      if (i > 0)
        nbr += 10;

//      part = FromStack(PARTICLE_TYPE_GAMMA, id);
part = 0;

      TTBenergy(mat, part, Eed[i], 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                1.0, 0, MT_ELECTRON_PP_EL,
                idxout, nbr, lEe, &Edep, id);

//      long partTTB = FromQue(id);
//      ToStack(part, id);

      /* Test positron number yield */
      TTByield(mat, Eed[i], 0, &lEe, MT_ELECTRON_PP_POS, &nbr, &idxout, id);

      /* Force TTB photon production */
      if (i > 0)
        nbr += 10;

      TTBenergy(mat, part, Eed[i], 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                1.0, 0, MT_ELECTRON_PP_POS,
                idxout, nbr, lEe, &Edep, id);

      /* With log-energy */
      lEe = lEed[i];

      /* Test electron number yield */
      TTByield(mat, Eed[i], 1, &lEe, MT_ELECTRON_PP_EL, &nbr, &idxout, id);

      /* Test positron number yield */
      TTByield(mat, Eed[i], 1, &lEe, MT_ELECTRON_PP_POS, &nbr, &idxout, id);


      /* Increase the energy slightly */
      if (i < nEe - 1) {
        double Ee = Eed[i]*(1.0 + 1.0e-15);
        TTByield(mat, Ee, 0, &lEe, MT_ELECTRON_PP_EL, &nbr, &idxout, id);
        TTBenergy(mat, part, Ee, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                  1.0, 0, MT_ELECTRON_PP_EL,
                  idxout, nbr, lEe, &Edep, id);
      }
      /* Decrease the energy slightly */
      if (i > 0) {
        double Ee = Eed[i]*(1.0 - 1.0e-15);
        TTByield(mat, Ee, 0, &lEe, MT_ELECTRON_PP_EL, &nbr, &idxout, id);
        TTBenergy(mat, part, Ee, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                  1.0, 0, MT_ELECTRON_PP_EL,
                  idxout, nbr, lEe, &Edep, id);
      }

//printf("unit test %ld\n", nbr);
    }

    /* Next material */
    mat = NextItem(mat);
  }
#endif

}

/*****************************************************************************/


/*****************************************************************************/

void ProcessElectronsInput() {
  /* Processes user input.
   * */
  static char * const FUNCTION_NAME = "ProcessElectronsInput:";
  long ptr0, ptr1, mat;

  fprintf(outp, " - processing input options...\n");

  /***** Materials given in elmee card ***************************************/

  /* Loop over materials given in the card */
  ptr0 = (long)RDB[DATA_ELECTRON_MEE_MAT];
  ptr1 = (long)RDB[DATA_ELECTRON_MEE];

  if ((ptr0 > VALID_PTR) && (ptr1 > VALID_PTR)) {

    while ((long)RDB[ptr0] > VALID_PTR) {

      /* Check pointer */
      if ((long)RDB[ptr1] < 0)
        Die(FUNCTION_NAME, "Pointer error: elmee for material %s",
            GetText(ptr0));

      /* Loop over materials */
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR) {

        /* Skip parent materials */
        if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
          mat = NextItem(mat);
          continue;
        }

        /* Compare material name */
        if (CompareStr(ptr0, mat + MATERIAL_PTR_NAME)) {
          /* Set the value */
          WDB[(long)RDB[mat + MATERIAL_PTR_EL] + EL_MEE] = RDB[ptr1];
          break;
        }

        mat = NextItem(mat);
      }

      /* Check pointer */
      if (mat < VALID_PTR)
        Note(0, "Material %s given in elmee card not found", GetText(ptr0));

      ptr0++;
      ptr1++;
    }
  }

  /***************************************************************************/

  /***** Materials given in elcond card **************************************/

  /* Loop over materials given in the card */
  ptr0 = (long)RDB[DATA_ELECTRON_COND_MAT];
  ptr1 = (long)RDB[DATA_ELECTRON_COND];

  if ((ptr0 > VALID_PTR) && (ptr1 > VALID_PTR)) {

    while ((long)RDB[ptr0] > VALID_PTR) {

      /* Check pointer */
      if ((long)RDB[ptr1] < 0)
        Die(FUNCTION_NAME, "Pointer error: elcond for material %s",
            GetText(ptr0));

      /* Loop over materials */
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR) {

        /* Skip parent materials */
        if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
          mat = NextItem(mat);
          continue;
        }

        /* Compare material name */
        if (CompareStr(ptr0, mat + MATERIAL_PTR_NAME)) {
          /* Set the value */
          WDB[(long)RDB[mat + MATERIAL_PTR_EL] + EL_COND] = RDB[ptr1];
          break;
        }

        mat = NextItem(mat);
      }

      /* Check pointer */
      if (mat < VALID_PTR)
        Note(0, "Material %s given in elcond card not found", GetText(ptr0));

      ptr0++;
      ptr1++;
    }
  }

  /***************************************************************************/

  /***** Materials given in elgas card ***************************************/

  /* Loop over materials given in the card */
  ptr0 = (long)RDB[DATA_ELECTRON_GAS_MAT];
  ptr1 = (long)RDB[DATA_ELECTRON_GAS];

  if ((ptr0 > VALID_PTR) && (ptr1 > VALID_PTR)) {

    while ((long)RDB[ptr0] > VALID_PTR) {

      /* Check pointer */
      if ((long)RDB[ptr1] < 0)
        Die(FUNCTION_NAME, "Pointer error: elgas for material %s",
            GetText(ptr0));

      /* Loop over materials */
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR) {

        /* Skip parent materials */
        if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
          mat = NextItem(mat);
          continue;
        }

        /* Compare material name */
        if (CompareStr(ptr0, mat + MATERIAL_PTR_NAME)) {
          /* Set the value */
          WDB[(long)RDB[mat + MATERIAL_PTR_EL] + EL_GAS] = RDB[ptr1];
          break;
        }

        mat = NextItem(mat);
      }

      /* Check pointer */
      if (mat < VALID_PTR)
        Note(0, "Material %s given in elgas card not found", GetText(ptr0));

      ptr0++;
      ptr1++;
    }
  }

  /***************************************************************************/

}

/*****************************************************************************/


/*****************************************************************************/

static void IntegrateSPrad(ElBrSXSData *SXSData, long idx, double *SPrad,
                           double multiplier, long mode) {
  /* Calculates radiative stopping powers from the bremsstrahlung data.
   * */
  static char * const FUNCTION_NAME = "IntegrateSPrad:";
  long i, j;
  double integral;
  double beta2, kappa0, SXS0;
  long idxk0;


  if (mode == 0) {
    /* 0 -> Ecut */

    idxk0 = SXSData->nkappa - 1;

    for (i = 0; i < SXSData->nE; i++) {
      integral = 0.0;
      kappa0 = RDB[DATA_PHOTON_EMIN]/SXSData->E[i];

      /* Find the interval for kappa0 */
      while ((kappa0 <= SXSData->kappa[idxk0]) && (--idxk0 >= 0));

      if (kappa0 < SXSData->kappa[idxk0])
        Die(FUNCTION_NAME, "kappa0 %f below lower limit", kappa0);
      if (kappa0 > SXSData->kappa[idxk0+1])
        Die(FUNCTION_NAME, "kappa0 %f above upper limit", kappa0);

      /* Interpolate SXS (lin-lin) */
      SXS0 = ENDFInterp(2, kappa0, SXSData->kappa[idxk0],
                        SXSData->kappa[idxk0+1], SXSData->SXSeT[idx][idxk0][i],
                        SXSData->SXSeT[idx][idxk0+1][i]);

      for (j = 1; j <= idxk0; j++) {
        integral += 0.5*(SXSData->SXSe[idx][i][j] + SXSData->SXSe[idx][i][j-1])
                    *(SXSData->kappa[j] - SXSData->kappa[j-1]);
      }

      /* Last interval */
      integral += 0.5*(SXS0 + SXSData->SXSe[idx][i][idxk0])
                  *(kappa0 - SXSData->kappa[idxk0]);

      beta2 = SXSData->E[i]*(SXSData->E[i] + 2.0*E_RESTMASS)
              /((SXSData->E[i] + E_RESTMASS)*(SXSData->E[i] + E_RESTMASS));
      SPrad[i] = multiplier/beta2*SXSData->E[i]*integral;
    }

  }
  else if (mode == 1) {
    /* Ecut -> E */

    idxk0 = SXSData->nkappa - 1;

    for (i = 0; i < SXSData->nE; i++) {
      integral = 0.0;
      kappa0 = RDB[DATA_PHOTON_EMIN]/SXSData->E[i];

      /* Find the interval for kappa0 */
      while ((kappa0 <= SXSData->kappa[idxk0]) && (--idxk0 >= 0));

      if (kappa0 < SXSData->kappa[idxk0])
        Die(FUNCTION_NAME, "kappa0 %f below lower limit", kappa0);
      if (kappa0 > SXSData->kappa[idxk0+1])
        Die(FUNCTION_NAME, "kappa0 %f above upper limit", kappa0);

      /* Interpolate SXS (lin-lin) */
      SXS0 = ENDFInterp(2, kappa0, SXSData->kappa[idxk0],
                        SXSData->kappa[idxk0+1], SXSData->SXSeT[idx][idxk0][i],
                        SXSData->SXSeT[idx][idxk0+1][i]);

      /* Trapezoid rule */
      integral += 0.5*(SXSData->SXSe[idx][i][idxk0+1] + SXS0)
                  *(SXSData->kappa[idxk0+1] - kappa0);

      /* The rest */
      for (j = idxk0 + 2; j < SXSData->nkappa; j++) {
        integral += 0.5*(SXSData->SXSe[idx][i][j] + SXSData->SXSe[idx][i][j-1])
                    *(SXSData->kappa[j] - SXSData->kappa[j-1]);
      }

      beta2 = SXSData->E[i]*(SXSData->E[i] + 2.0*E_RESTMASS)
              /((SXSData->E[i] + E_RESTMASS)*(SXSData->E[i] + E_RESTMASS));
      SPrad[i] = multiplier/beta2*SXSData->E[i]*integral;
    }

    /* Set small non-zero first value to allow log-log interpolation */
    SPrad[0] = 1.0e-8*SPrad[1];

  }
  else if (mode == 2) {
    /* 0 -> E */

    for (i = 0; i < SXSData->nE; i++) {
      integral = 0.0;

      /* Trapezoid rule */
      for (j = 1; j < SXSData->nkappa; j++) {
        integral += 0.5*(SXSData->SXSe[idx][i][j] + SXSData->SXSe[idx][i][j-1])
                    *(SXSData->kappa[j] - SXSData->kappa[j-1]);
      }

      beta2 = SXSData->E[i]*(SXSData->E[i] + 2.0*E_RESTMASS)
              /((SXSData->E[i] + E_RESTMASS)*(SXSData->E[i] + E_RESTMASS));
      SPrad[i] = multiplier/beta2*SXSData->E[i]*integral;
    }
  }

}

/*****************************************************************************/


/*****************************************************************************/

static double ElIcompESTAR(const long *Z, const double *adens,
                           const double *Ielements, long nZ, long isgas) {
  /* Calculates and returns the mean excitation energy I of a compound or
   * mixture using the Bragg's additivity rule. The I-values of the atomic
   * constituents given here are from [1], and most of them can also be found
   * in [2].
   *
   * [1] M. J. Berger, J. Coursey, M. Zucker, and J. Chang, "Stopping-power
   *     and range tables for electrons, protons, and helium ions," NIST,
   *     http://www.nist.gov/pml/data/star/index.cfm
   * [2] M. J. Berger, "Electron Stopping Powers for Transport Calculations,"
   *     in Monte Carlo Transport of Electrons and Photons, volume 38 of Ettore
   *     Majorana International Science Series, p. 57–80, Springer US, 1988.
   * */
  long i;
  double Ii, Imat, eldens;

  Imat = eldens = 0.0;

  for (i = 0; i < nZ; i++) {

    /* Set the I-values (in MeV) of the atomic constituents. */

    if (isgas) {
        /* Gas phase */
        switch (Z[i]) {
          case 1:  /* Hydrogen */
              Ii = 19.2e-6;
              break;
          case 3:  /* WTF-value from ESTAR: Lithium */
              Ii = 34.e-6;
              break;
          case 4:  /* WTF-value from ESTAR: Beryllium */
              Ii = 38.6e-6;
              break;
          case 5:  /* WTF-value from ESTAR: Boron */
              Ii = 49.e-6;
              break;
          case 6:  /* Carbon */
              Ii = 70.e-6;
              break;
          case 7:  /* Nitrogen */
              Ii = 82.e-6;
              break;
          case 8:  /* Oxygen */
              Ii = 97.e-6;
              break;
          default:
            Ii = Ielements[i];
        }
    }
    else {
        /* Condensed material (liquid or solid) */
        switch (Z[i]) {
          case 1:  /* Hydrogen */
              Ii = 19.2e-6;
              break;
          case 6:  /* Carbon */
              Ii = 81.e-6;
              break;
          case 7:  /* Nitrogen */
              Ii = 82.e-6;
              break;
          case 8:  /* Oxygen */
              Ii = 106.e-6;
              break;
          case 9:  /* Fluorine */
              Ii = 112.e-6;
              break;
          case 17:  /* Chlorine */
              Ii = 180.e-6;
              break;
          case 35:  /* Bromine (condensed value) */
              Ii = 357.e-6;
              break;
          default:
              Ii = 1.13*Ielements[i];
        }
    }

    /* Use the Bragg's additivity rule */
    Imat += (double)Z[i]*adens[i]*log(Ii);
    eldens += (double)Z[i]*adens[i];
  }

  /*Imat = exp(Imat/ZperAmat);*/
  Imat = exp(Imat/eldens);

  return Imat;
}

/*****************************************************************************/


/*****************************************************************************/

static void ElDeffcorSternheimer(double eldens, double Imat, const double *fel,
                                 const double *ebi, long nss, const double *E,
                                 long nE, double *delta) {
  /* Solves the density effect correction using the Sternheimer's method.
   *
   * See:
   * R. M. Sternheimer, S. M. Seltzer, and M. J. Berger,
   * "Density effect for the ionization loss of charged particles in various
   * substances," Phys. Rev. B, 26, 6067-6076 (1982)
   * */
  static char * const FUNCTION_NAME = "ElDeffcorSternheimer:";
  static const double plc = 28.8159382677053e-6; /* Plasma energy constant */
  static const double tol = 1.0e-7;   /* Tolerance for the Newton's method */
  static const long itermax = 1e5;
  long i, j, iter;
  double lImat, Epl2, mu, muold, mu2, gmu, dgmu, l, l2, lold, beta2,
      betath2, gl, dgl, fiEff, EiL, gmuEff, l2f, gli, beta2i;
  double fmeb[PHOTON_NSS_MAX*PHOTON_ZMAX];
  double li2[PHOTON_NSS_MAX*PHOTON_ZMAX];

  /***** Set variables *******************************************************/

  /* Plasma energy squared */
  Epl2 = plc*plc*eldens/N_AVOGADRO;

  /* Effective number of conduction electrons */
  fiEff = 0.0;
  for (i = 0; i < nss; i++) {
    if (fel[i] < 0.0) /* Conduction electron */
      fiEff += fabs(fel[i]);
  }

  lImat = log(Imat);

  if (fiEff == 0.0)
    gmuEff = 0.0;
  else
    gmuEff = fiEff*log(fiEff*Epl2);

  /***************************************************************************/

  /***** Solve the parameter mu using the Newton's method ********************/

  mu = 0.5; /* NOTE: The initial guess must be positive */
  iter = 0;

  do {
    /* Function and the derivative */
    gmu = 0.0;
    dgmu = 0.0;

    /* Loop over shells */
    for (i = 0; i < nss; i++) {

      /* Skip conduction electrons */
      if (fel[i] < 0.0)
        continue;

      EiL = ebi[i]*ebi[i]*mu*mu + 2.0/3.0*fel[i]*Epl2;

      /* Function at mu */
      gmu += fel[i]*log(EiL);

      /* Derivative at mu */
      dgmu += ebi[i]*ebi[i]*fel[i]*mu/EiL;
    }

    /* Include conduction electrons */
    if (fiEff != 0.0)
      gmu += gmuEff;

    gmu = gmu/2.0 - lImat;

    /* Update the solution */
    muold = mu;
    mu = muold - gmu/dgmu;

    /* Too large initial guess can lead into negative mu, decrease the value */
    if (mu < 0.0)
      mu = muold/2.0;

    /* Limit iterations */
    if (++iter > itermax)
      Die(FUNCTION_NAME, "SolveDelta mu: maximum number of iterations exceeded");

  } while (fabs(mu/muold - 1.0) > tol);

  mu2 = mu*mu;

  /***************************************************************************/

  /* Set parameters */
  for (i = 0; i < nss; i++) {
    fmeb[i] = mu2*ebi[i]*ebi[i]/Epl2;
    li2[i] = fmeb[i] + 2.0/3.0*fel[i];
  }

  /* Calculate threshold beta squared */
  betath2 = 0.0;
  if (fiEff == 0.0) {
    for (i = 0; i < nss; i++)
      betath2 += fel[i]/fmeb[i];
    betath2 = 1.0/(1.0 + betath2);
  }

  /***** Loop over electron energies *****************************************/

  /* Initial guess used for the first energy (must be positive!) */
  l = 0.01;

  for (j = 0; j < nE; j++) {

    /* v^2/c^2 */
    beta2 = E[j]*(E[j] + 2.0*E_RESTMASS)/((E[j] + E_RESTMASS)
                                          *(E[j] + E_RESTMASS));
    beta2i = 1.0/beta2;
    iter = 0;

    /* l and delta are zero when fiEff is zero and beta < betath
     * (see PhysRevB.26.6067) */
    if (beta2 < betath2) {
      delta[j] = 0.0;
      continue;
    }
    else {
      /* 1.07*l seems to be a good initial guess when logarithmic
       * energy grid is used */
        l = 1.07*l;
    }

    /***** Solve parameter l using the Newton's method with mu solved above **/

    do {
      gl = 0.0;
      dgl = 0.0;
      l2 = l*l;

      /* Loop over shells */
      for (i = 0; i < nss; i++) {

        /* Skip conduction electrons */
        if (fel[i] < 0.0)
          continue;

        l2f = 1.0/(fmeb[i] + l2);
        gli = fel[i]*l2f;

        /* Function at l */
        gl += gli;

        /* Derivative at l */
        dgl += gli*l2f;
      }

      /* Include conduction electrons */
      if (fiEff != 0.0) {
        gli = fiEff/l2;
        gl += gli;
        dgl += gli/l2;
      }

      gl += 1.0 - beta2i;
      dgl *= -2.0*l;

      /* Update the solution */
      lold = l;
      l = lold - gl/dgl;

      /* Too large initial guess can lead into negative l, decrease the value */
      if (l < 0.0)
        l = lold/2.0;

      /* Limit iterations */
      if (++iter > itermax)
        Die(FUNCTION_NAME, "SolveDelta l: maximum number of iterations exceeded");

    } while (fabs(l - lold) > tol*lold);

    l2 = l*l;

    /*************************************************************************/

    /***** Calculate the density effect correction ***************************/

    delta[j] = 0.0;

    /* Loop over shells */
    for (i = 0; i < nss; i++) {

      /* Skip conduction electrons */
      if (fel[i] < 0.0)
        continue;

      delta[j] += fel[i]*log((li2[i] + l2)/li2[i]);
    }

    /* Include conduction electrons */
    if (fiEff != 0.0)
      delta[j] += fiEff*log((fiEff + l2)/fiEff);

    delta[j] -= l2*(1.0 - beta2);

    /*************************************************************************/  
  }

  /***************************************************************************/

}

/*****************************************************************************/


/*****************************************************************************/

static double PositronRp(double E, double Z2eff) {
  /* Calculates the scaling factor Rp for positron bremsstrahlung DCS.
   * Rp is defined as the approximative ratio of the radiative stopping
   * powers for positrons and electrons.
   *   Z2eff : "effective" squared atomic number
   *
   * Source: F. Salvat et al., "PENELOPE-2011: A Code System for Monte Carlo
   * Simulation of Electron and Photon Transport" (2011)
   * */
  double tp, Rp;
  tp = log(1.0 + 1.0e6/Z2eff*E/E_RESTMASS);
  Rp = 1.0 - exp(tp*(-1.2359e-1 + tp*(6.1274e-2 + tp*(-3.1516e-2 + tp*(7.7446e-3
              + tp*(-1.0595e-3 + tp*(7.0568e-5 - 1.8080e-6*tp)))))));
  return Rp;
}

/*****************************************************************************/


/*****************************************************************************/

void ProcessSPcolBethe(const double *E, long nE, double Imat, double eldens,
                       double *delta, double *SPcole, double *SPcolp) {
  /* Calculates collision stopping powers (in MeV/cm) for electrons and
   * positrons with Bethe's formula.
   * */
  static char * const FUNCTION_NAME = "ProcessSPcolBethe:";
  static const double re = 2.8179403267e-13;  /* classical electron radius [cm] */
  long i;
  double tau, beta2, Fe, Fp, SPconst, SPtmp;

  if (!delta)
    Die(FUNCTION_NAME, "Density effect correction array not found");


  SPconst = 2.0*PI*re*re*E_RESTMASS*eldens/BARN;

  /* Loop over electron energy */
  for (i = 0; i < nE; i++) {

    tau = E[i]/E_RESTMASS;
    beta2 = E[i]*(E[i] + 2.0*E_RESTMASS)/((E[i] + E_RESTMASS)
                                          *(E[i] + E_RESTMASS));

    /* Variable Fe for electrons and Fp for positrons */
    Fe = (1.0 - beta2)*(1.0 + tau*tau/8.0 - (2.0*tau + 1.0)*LOG2);
    Fp = 2.0*LOG2 - (beta2/12.0)*
         (23.0 + (14.0 + (10.0 + 4.0/(tau + 2.0))/(tau + 2.0))/(tau + 2.0));

    /* Electron and positron collision stopping power at the particle energy */
    SPtmp = SPconst/beta2*(log(E[i]*E[i]/(Imat*Imat)) + log(1.0 + tau/2.0)
                           - delta[i]);
    SPcole[i] = SPtmp + SPconst/beta2*Fe;
    SPcolp[i] = SPtmp + SPconst/beta2*Fp;
  }

}

/*****************************************************************************/


/*****************************************************************************/

static void RadiativeYield(const double *E, const double *SPrad,
                           const double *SPtot, long nE,  double *Yrad) {
  /* Integrates radiative yield.
   * */
  long i;
  double *radPerTot;

  /* Allocate memory */
  radPerTot = (double *)Mem(MEM_ALLOC, nE, sizeof(double));

  /* Loop over energies */
  for (i = 0; i < nE; i++)
    radPerTot[i] = SPrad[i]/SPtot[i];

  /* Calculate yield (log-log) */
  TrapzRealCum(E, radPerTot, Yrad, nE, 5);

  /* Loop over energies */
  for (i = 0; i < nE; i++)
    Yrad[i] /= E[i];

  Mem(MEM_FREE, radPerTot);
}

/*****************************************************************************/


/*****************************************************************************/

static void PrintSPDataMat(FILE *fp, ElSPData *SPDataMat) {
  /* Prints stopping power data in Matlab format file.
   * */
  static char * const FUNCTION_NAME = "PrintSPDataMat";
  long i, j, mat, mat0, matfound, elptr;
  double rho;
  double *Yrade, *Yradp;

  fprintf(outp, " - printing electron stopping power data...\n");

  Yrade = (double *)Mem(MEM_ALLOC, SPDataMat->nE, sizeof(double));
  Yradp = (double *)Mem(MEM_ALLOC, SPDataMat->nE, sizeof(double));

  /* Print header */
  fprintf(fp, "%% Stopping powers and radiative yields for electrons and positrons\n");
  fprintf(fp, "%%  - Serpent version        : %s\n", CODE_VERSION);
  fprintf(fp, "%%\n");
  fprintf(fp, "%%  Units:\n");
  fprintf(fp, "%%  - mass density           : g/cm^3\n");
  fprintf(fp, "%%  - mean excitation energy : MeV\n");
  fprintf(fp, "%%  - energy                 : MeV\n");
  fprintf(fp, "%%  - stopping power (SP)    : MeV cm^2/g\n");
  fprintf(fp, "%%\n");
  fprintf(fp, "%%  Notes:\n");
  fprintf(fp, "%%  - \"delta\" is the density effect correction (dimensionless)");
  fprintf(fp, "\n\n");


  for (i = 0; i < SPDataMat->sz; i++) {

    /* Material pointer */
    mat = SPDataMat->map[i];

    /* Check that the material exists */
    matfound = 0;
    mat0 = (long)RDB[DATA_PTR_M0];
    while (mat0 > VALID_PTR) {
      if (mat0 == mat) {
        matfound = 1;
        break;
      }
      mat0 = NextItem(mat0);
    }
    if (!matfound)
      Die(FUNCTION_NAME, "Material pointer %ld not found", mat);

    /* Calculate radiative yields */
    RadiativeYield(SPDataMat->E, SPDataMat->SPrade[i], SPDataMat->SPtote[i],
                   SPDataMat->nE, Yrade);
    RadiativeYield(SPDataMat->E, SPDataMat->SPradp[i], SPDataMat->SPtotp[i],
                   SPDataMat->nE, Yradp);


    /* Pointer to electron data */
    elptr = (long)RDB[mat + MATERIAL_PTR_EL];

    /* Mass density */
    rho = (double)RDB[mat + MATERIAL_MDENS];

    /* Print labels */
    fprintf(fp, "%% Material : %s\n", GetText(mat + MATERIAL_PTR_NAME));
    fprintf(fp, "%% - mass density            : %E\n", rho);
    fprintf(fp, "%% - mean excitation energy  : %E\n", RDB[elptr + EL_MEE]);
    fprintf(fp, "%% - is conductor            : %ld\n", (long)RDB[elptr + EL_COND]);
    fprintf(fp, "%% - is gas                  : %ld\n", (long)RDB[elptr + EL_GAS]);
    fprintf(fp, "%%                ______________________Electrons_______________________"
                "   ______________________Positrons_______________________\n");
    fprintf(fp, "%%       energy         SP col        SP rad        SP tot     rad yield"
                "         SP col        SP rad       SP tot      rad yield          delta\n");

    /* Print material variable name */
    fprintf(fp, "ELSP_%s = [\n", GetText(mat + MATERIAL_PTR_NAME));

    /* Loop over energy and print data */
    for (j = 0; j < SPDataMat->nE; j++) {

      fprintf(fp, "  %8E   %8E  %8E  %8E  %8E   %8E  %8E  %8E  %8E   %8E\n",
              SPDataMat->E[j],
              SPDataMat->SPcole[i][j]/rho, SPDataMat->SPrade[i][j]/rho,
              SPDataMat->SPtote[i][j]/rho, Yrade[j],
              SPDataMat->SPcolp[i][j]/rho, SPDataMat->SPradp[i][j]/rho,
              SPDataMat->SPtotp[i][j]/rho, Yradp[j],
              SPDataMat->delta[i][j]);
    }

    fprintf(fp, "];\n\n");
  }


  Mem(MEM_FREE, Yrade);
  Mem(MEM_FREE, Yradp);
}

/*****************************************************************************/


#ifdef TEST_EL_FUTURE
/*****************************************************************************/

static void BrMFP(ElBrSXSData *SXSData, long idx, double Ecut, double multiplier,
                  double *mfp) {

  static char * const FUNCTION_NAME = "BrMFP";
  long i, j, idxkc, mode, i0;
  double kappac, Ei, SXSc, integral, beta2;
  double *kappai;
  double *integrand;

  mode = 1;   /* log-log */

  integrand = (double *)Mem(MEM_ALLOC, SXSData->nE, sizeof(double));
  kappai = (double *)Mem(MEM_ALLOC, SXSData->nkappa, sizeof(double));

  /* if Ecut == SXSData->E[0] ? */
  mfp[0] = 0;
  i0 = 1;

  idxkc = SXSData->nkappa - 1;

  for (i = i0; i < SXSData->nE; i++) {
    Ei = SXSData->E[i];

    kappac = Ecut/Ei;


    /* Find the interval for kappac */
    while ((kappac <= SXSData->kappa[idxkc]) && (--idxkc >= 0));

    if (kappac < SXSData->kappa[idxkc])
      Die(FUNCTION_NAME, "kappa0 %f below lower limit", kappac);
    if (kappac > SXSData->kappa[idxkc+1])
      Die(FUNCTION_NAME, "kappa0 %f above upper limit", kappac);


    /* Interpolate SXS (lin-lin) */
    SXSc = ENDFInterp(2, kappac, SXSData->kappa[idxkc], SXSData->kappa[idxkc+1],
                       SXSData->SXSeT[idx][idxkc][i], SXSData->SXSeT[idx][idxkc+1][i]);

    kappai[idxkc] = kappac;
    integrand[idxkc] = SXSc/kappac;

    for (j = idxkc+1; j < SXSData->nkappa; j++) {
      kappai[j] = SXSData->kappa[j];
      integrand[j] = SXSData->SXSe[idx][i][j]/SXSData->kappa[j];
    }

    integral = LogIntegral(&kappai[idxkc], &integrand[idxkc],
                           SXSData->nkappa - idxkc, mode);

    beta2 = Ei*(Ei + 2.0*E_RESTMASS)/((Ei + E_RESTMASS)*(Ei + E_RESTMASS));
    mfp[i] = multiplier/beta2*integral;

  }

  Mem(MEM_FREE, integrand);
  Mem(MEM_FREE, kappai);


}

/*****************************************************************************/


/*****************************************************************************/

static void RangeCSDA(double *E, double *Stot, double *R, long n) {
  long i;
  double *integrand;

  integrand = (double *)Mem(MEM_ALLOC, n, sizeof(double));

  for (i = 0; i < n; i++)
    integrand[i] = 1.0/Stot[i];

  CumLogLogIntegral(E, NULL, integrand, R, n);

  for (i = 0; i < n; i++)
    printf("%f  %f\n", E[i], R[i]);

  Mem(MEM_FREE, integrand);
}

/*****************************************************************************/

#endif
