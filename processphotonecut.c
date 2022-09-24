/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processphotonecut.c                            */
/*                                                                           */
/* Created:       2016/12/21 (TKa)                                           */
/* Last modified: 2017/04/10 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Processes material-wise photon energy cutoffs                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessPhotonEcut:"

/*****************************************************************************/

void ProcessPhotonEcut() {

  long mat, ptr0, ptr1, nuc, iso, rea, reatot, ptd, ptr, nss, nssmat, bsflag,
      i, maxwidth, lenstr, warnflag;
  double E1, E2, E3, mfp1, mfp2, mfp3, mfpmat, dEb3, Ecutadj;
  double const *Ebd;
  double Ebmat[PHOTON_ZMAX*PHOTON_NSS_MAX];
  const double Emaxmfp = 1.0;  /* Maximum cutoff energy when mfp cutoff is used */
  const double Eminmfp = 1.0000001*RDB[DATA_PHOTON_EMIN];
  const double bstol = 1.0E-5; /* Tolerance used for solving the cutoff energy */
  const double dEb = 1.0E-4;   /* Adjustment factor for binding energies */

  fprintf(outp, "Setting photon energy cutoffs...\n");


  /***** Set energy and mean free path cutoffs given in the input cards *****/

  /* The cutoffs are set in the following order:
   * - ecut       : energy cutoff applied for all materials
   * - mfpcut     : mfp cutoff applied for all materials
   * - ecutdens   : list of mass density and energy cutoff pairs, each cutoff
   *                is applied for materials with densities above the given
   *                density
   * - mfpcutdens : list of mass density and mfp cutoff pairs, each cutoff is
   *                applied for materials with densities above the given
   *                density
   * - ecutmat    : list of material name and energy cutoff pairs, each cutoff
   *                is applied for the given material
   * - mfpcutmat  : list of material name and mfp cutoff pairs, each cutoff
   *                is applied for the given material
   *
   * Note that the mean free paths that are first set here are used for
   * solving the corresponding approximate energy cutoffs. The actual mean
   * free paths corresponding to the cutoff energies (given or solved) are
   * set at the end. */

  /***** Initialize material energy cutoffs **********************************/

  /* Loop over materials */
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR) {

    /* Skip parent materials */
    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
      mat = NextItem(mat);
      continue;
    }

    WDB[mat + MATERIAL_PHOTON_ECUT] = RDB[DATA_PHOTON_ECUT];

    mat = NextItem(mat);
  }

  /***************************************************************************/


  /***** Initialize material mfp cutoffs *************************************/

  mat = (long)RDB[DATA_PTR_M0];

  while (mat > VALID_PTR) {

    /* Skip parent materials */
    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
      mat = NextItem(mat);
      continue;
    }

    WDB[mat + MATERIAL_PHOTON_MFPCUT] = RDB[DATA_PHOTON_MFPCUT];

    mat = NextItem(mat);
  }

  /***************************************************************************/


  /***** Process materials given in ecutdens card ****************************/

  /* Loop over densities given in the card */
  ptr0 = (long)RDB[DATA_PHOTON_ECUTDENS_DENS];
  ptr1 = (long)RDB[DATA_PHOTON_ECUTDENS_E];

  warnflag = 1;

  if ((ptr0 > VALID_PTR) && (ptr1 > VALID_PTR)) {

    while (RDB[ptr0] > 0.0) {

      /* Check cutoff */
      if ((long)RDB[ptr1] < 0)
        Die(FUNCTION_NAME, "Energy cutoff error in ecutdens for material %s",
            GetText(ptr0));

      /* Loop over materials */
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR) {

        /* Skip parent materials */
        if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
          mat = NextItem(mat);
          continue;
        }

        /* Compare material density to the given value */
        if (RDB[mat + MATERIAL_MDENS] >= RDB[ptr0]) {

          /* Check the given energy */
          if (RDB[ptr1] < RDB[DATA_PHOTON_EMIN]) {

            /* This should happen only when cross section minimum
             *  is above the grid minimum, see ProcessNuclides()) */

            if (warnflag)
              Warn(FUNCTION_NAME, "Energy %E MeV given in the ecutdens card is\n"
                 "below the minimum %E MeV.\n"
                 "The energy is set to %E.",
                 RDB[ptr1], RDB[DATA_PHOTON_EMIN], RDB[DATA_PHOTON_EMIN]);

            /* Set the value */
            WDB[mat + MATERIAL_PHOTON_ECUT] = RDB[DATA_PHOTON_EMIN];

            warnflag = 0;
          }
          else {
            /* Set the value */
            WDB[mat + MATERIAL_PHOTON_ECUT] = RDB[ptr1];
          }

          /* Set a flag for not using mfp to solve energy */
          WDB[mat + MATERIAL_PHOTON_MFPCUT] = -INFTY;
        }

        mat = NextItem(mat);
      }

      ptr0++;
      ptr1++;
    }
  }

  /***************************************************************************/


  /***** Set mfp cutoffs given in mfpcutdens card ****************************/

  ptr0 = (long)RDB[DATA_PHOTON_MFPCUTDENS_DENS];
  ptr1 = (long)RDB[DATA_PHOTON_MFPCUTDENS_MFP];

  if ((ptr0 > VALID_PTR) && (ptr1 > VALID_PTR)) {

    while (RDB[ptr0] > 0.0) {

      /* Check cutoff */
      if ((long)RDB[ptr1] < 0)
        Die(FUNCTION_NAME, "Error in mfp given in the mfpcutdens card "
                           "for material %s", GetText(ptr0));

      /* Loop over materials */
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR) {

        /* Skip parent materials */
        if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
          mat = NextItem(mat);
          continue;
        }

        /* Compare material density */
        if (RDB[mat + MATERIAL_MDENS] >= RDB[ptr0]) {
          WDB[mat + MATERIAL_PHOTON_MFPCUT] = RDB[ptr1];
        }

        mat = NextItem(mat);
      }

      ptr0++;
      ptr1++;
    }
  }

  /***************************************************************************/


  /***** Process materials given in ecutmat card *****************************/

  /* Loop over materials given in the card */
  ptr0 = (long)RDB[DATA_PHOTON_ECUTMAT_MAT];
  ptr1 = (long)RDB[DATA_PHOTON_ECUTMAT_E];

  if ((ptr0 > VALID_PTR) && (ptr1 > VALID_PTR)) {

    while ((long)RDB[ptr0] > VALID_PTR) {

      /* Check pointer */
      if ((long)RDB[ptr1] < 0)
        Die(FUNCTION_NAME, "Pointer error: ecut for material %s",
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

          /* Check the given energy (needed only when cross section minimum
           * is above the grid minimum, see ProcessNuclides()) */
          if (RDB[ptr1] < RDB[DATA_PHOTON_EMIN]) {

            Warn(FUNCTION_NAME,
                 "Energy %E MeV given in the ecutmat card for\n"
                 "material %s is below the minimum %E MeV.\n"
                 "The material cutoff energy is set to %E.", RDB[ptr1],
                 GetText(mat + MATERIAL_PTR_NAME), RDB[DATA_PHOTON_EMIN],
                 RDB[DATA_PHOTON_EMIN]);

            /* Set the value */
            WDB[mat + MATERIAL_PHOTON_ECUT] = RDB[DATA_PHOTON_EMIN];
          }
          else {
            /* Set the value */
            WDB[mat + MATERIAL_PHOTON_ECUT] = RDB[ptr1];
          }

          /* Set a flag for not using mfp to solve energy */
          WDB[mat + MATERIAL_PHOTON_MFPCUT] = -INFTY;

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


  /***** Set mfp cutoffs given in mfpcutmat card *****************************/

  ptr0 = (long)RDB[DATA_PHOTON_MFPCUTMAT_MAT];
  ptr1 = (long)RDB[DATA_PHOTON_MFPCUTMAT_MFP];

  if ((ptr0 > VALID_PTR) && (ptr1 > VALID_PTR)) {

    while ((long)RDB[ptr0] > VALID_PTR) {

      /* Check pointer */
      if ((long)RDB[ptr1] < 0)
        Die(FUNCTION_NAME, "Pointer error: mfpcutmat for material %s",
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
          WDB[mat + MATERIAL_PHOTON_MFPCUT] = RDB[ptr1];
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

  /***************************************************************************/


  /***** Convert mean free path cutoffs to energy cutoffs ********************/

  /* Loop over materials */
  mat = (long)RDB[DATA_PTR_M0];

  while (mat > VALID_PTR) {

    /* Skip parent materials */
    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
      mat = NextItem(mat);
      continue;
    }

    mfpmat = RDB[mat + MATERIAL_PHOTON_MFPCUT];

    /* Check if the mfp cutoff was not set */
    if (mfpmat < ZERO) {
      mat = NextItem(mat);
      continue;
    }

    /***** Set binding energy array for the material *************************/

    /* Binding energies given in the photoelectric data are used. */

    nssmat = 0;
    memset(Ebmat, 0, PHOTON_ZMAX*PHOTON_NSS_MAX*sizeof(double));

    /* Loop over composition */

    iso = (long)RDB[mat + MATERIAL_PTR_COMP];

    while (iso > VALID_PTR) {

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];

      /* Only photon nuclides are included */
      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON) {
        /* Next */
        iso = NextItem(iso);
        continue;
      }

      /* Find photoelectric data */
      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while ((rea > VALID_PTR) && ((long)RDB[rea + REACTION_MT] != 522))
        rea = NextItem(rea);

      if ((long)RDB[rea + REACTION_MT] != 522)
        Die(FUNCTION_NAME, "Photoelectric data (mt 522) not found for "
                           "material %s", GetText(mat + MATERIAL_PTR_NAME));

      ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
      CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

      /* Get number of shells */
      nss = (long)RDB[ptd + PHOTON_DIST_PE_NSS];

      /* No binding energies above the minimum cutoff */
      if (nss == 0) {
        iso = NextItem(iso);
        continue;
      }

      if (nss > PHOTON_NSS_MAX)
        Die(FUNCTION_NAME, "Number of shells %ld above the maximum %d", nss,
            PHOTON_NSS_MAX);

      /* Binding energies */
      ptr = (long)RDB[ptd + PHOTON_DIST_PE_EBI];
      CheckPointer(FUNCTION_NAME, "(EBI)", DATA_ARRAY, ptr);
      Ebd = &RDB[ptr];

      /* Store binding energies */
      for (i = 0; i < nss; i++)
        Ebmat[nssmat++] = Ebd[i];

      /* Next */
      iso = NextItem(iso);
    }

    /* Sort binding energies */
    if (nssmat > 0)
      SortArray(Ebmat, nssmat);

    /*************************************************************************/


    /***** Find the energy interval corresponding to the mfp *****************/

    /* The mfp is assumed to be monotonously increasing between the binding
     * energies and between the highest binding energy and Emax.
     * Note the usage of the adjustment factor dEb. */

    reatot = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
    bsflag = 1;   /* Flag for binary search */

    /* Mfp at the lower and upper limit */
    mfp1 = 1.0/PhotonMacroXS(reatot, Eminmfp, 0);
    mfp3 = 1.0/PhotonMacroXS(reatot, Emaxmfp, 0);


    if (mfpmat < mfp1) {

      /* mfp below the lower limit */

      Warn(FUNCTION_NAME, "Mean free path %E cm is below the material "
                          "minimum \n%E cm at %E MeV in material %s. "
                          "\nThe material cutoff energy is set to %E MeV.",
           mfpmat, mfp1, RDB[DATA_PHOTON_EMIN], RDB[DATA_PHOTON_EMIN],
           GetText(mat + MATERIAL_PTR_NAME));
      bsflag = 0;
      E2 = RDB[DATA_PHOTON_EMIN];
    }
    else if (mfpmat > mfp3) {

      /* mfp above the upper limit */

      Warn(FUNCTION_NAME, "Mean free path %E cm is above the "
                          "maximum \n%E cm at %E MeV in material %s. "
                          "\nThe material cutoff energy is set to %E MeV.",
           mfpmat, mfp3, Emaxmfp, Emaxmfp,
           GetText(mat + MATERIAL_PTR_NAME));
      bsflag = 0;
      E2 = Emaxmfp;
    }
    else if (nssmat == 0) {

      /* No binding energies found - use Emin and Emax*/

      E1 = Eminmfp;
      E3 = Emaxmfp;
    }
    else if ((mfp1 = 1.0/PhotonMacroXS(reatot, (1.0 + dEb)*Ebmat[nssmat-1], 0))
             < mfpmat) {

      /* The cutoff energy is between the highest binding energy and Emax */

      E1 = (1.0 + dEb)*Ebmat[nssmat-1];
      E3 = Emaxmfp;
    }
    else {

      /* Find the energy interval between the binding energies. */

      i = 0;
      E1 = Eminmfp;

      while (1) {

        /* Set the upper energy limit */
        dEb3 = dEb;
        do {
          E3 = (1.0 - dEb3)*Ebmat[i];
          dEb3 *= 0.9;
        } while (E3 < E1);

        /* Calculate mfp at the lower and upper limit */
        mfp1 = 1.0/PhotonMacroXS(reatot, E1, 0);
        mfp3 = 1.0/PhotonMacroXS(reatot, E3, 0);

        /* Test the material mfp */
        if (mfpmat < mfp1) {
          Warn(FUNCTION_NAME, "Warning 1:\nEnergy interval not found for mean"
               " free path %f cm in material %s. \nThe material cutoff energy"
               " is set to %E MeV.", mfpmat, GetText(mat + MATERIAL_PTR_NAME),
               RDB[DATA_PHOTON_EMIN]);
          bsflag = 0;
          E2 = RDB[DATA_PHOTON_EMIN];
          break;
        }
        else if (mfpmat > mfp3) {
          /* Next interval */
          E1 = (1.0 + dEb)*Ebmat[i];
          i++;

          if (i > nssmat) {
            Warn(FUNCTION_NAME, "Warning 2:\nEnergy interval not found for "
                 "mean free path %f cm in material %s. \nThe material cutoff "
                 "energy is set to %E MeV.", mfpmat,
                 GetText(mat + MATERIAL_PTR_NAME), RDB[DATA_PHOTON_EMIN]);
            bsflag = 0;
            E2 = RDB[DATA_PHOTON_EMIN];
            break;
          }
        }
        else {
          /* Interval found */
          break;
        }
      }

    }

    /*************************************************************************/


    /***** Perform binary search *********************************************/

    if (bsflag) {
      do {
        E2 = 0.5*(E3 + E1);
        mfp2 = 1.0/PhotonMacroXS(reatot, E2, 0);

        if (mfp2 < mfpmat) {
          E1 = E2;
          mfp1 = mfp2;
        }
        else {
          E3 = E2;
          mfp3 = mfp2;
        }

      } while(fabs(mfp1/mfp3 - 1.0) > bstol);

      E2 = 0.5*(E3 + E1);

      /* Check energy */

      if ((E2 < E1) || (E2 > E3)) {
        Warn(FUNCTION_NAME, "Binary search failed to find the mean "
             "free path %f cm in material %s. \nThe material cutoff energy "
             "is set to %E MeV.", mfpmat, GetText(mat + MATERIAL_PTR_NAME),
             RDB[DATA_PHOTON_EMIN]);
        E2 = RDB[DATA_PHOTON_EMIN];
      }
    }

    /*************************************************************************/

    /* Set the cutoff energy */
    WDB[mat + MATERIAL_PHOTON_ECUT] = E2;

    /* Next material */
    mat = NextItem(mat);
  }

  /***************************************************************************/


  /***** Set mean free path cutoffs corresponding to the energy cutoffs ******/

  /* Find the longest material name for output printing */
  maxwidth = 8;
  mat = (long)RDB[DATA_PTR_M0];

  while (mat > VALID_PTR) {
    lenstr = strlen(GetText(mat + MATERIAL_PTR_NAME));
    if (lenstr > maxwidth)
      maxwidth = lenstr;
    mat = NextItem(mat);
  }
  maxwidth++;

  /* Print header */
  fprintf(outp, "%*s  %19s  %19s\n", (int)maxwidth, "Material",
                "Cutoff energy (MeV)", "Mean free path (cm)");

  /* Loop over materials */
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR) {

    /* Skip parent materials */
    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT) {
      mat = NextItem(mat);
      continue;
    }

    Ecutadj = RDB[mat + MATERIAL_PHOTON_ECUT];

    /* Adjust the cutoff if needed so that PhotonMacroXS() won't fail */
    if (Ecutadj == RDB[DATA_PHOTON_EMIN])
      Ecutadj *= 1.00000001;

    /* Set the material mfp cutoff */
    reatot = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
    WDB[mat + MATERIAL_PHOTON_MFPCUT] = 1.0/PhotonMacroXS(reatot, Ecutadj, 0);

    /* Print the material, energy and the mfp */
    fprintf(outp, "%*s  %19E  %19E\n", (int)maxwidth,
            GetText(mat + MATERIAL_PTR_NAME), RDB[mat + MATERIAL_PHOTON_ECUT],
            RDB[mat + MATERIAL_PHOTON_MFPCUT]);

    /* Next material */
    mat = NextItem(mat);
  }

  /***************************************************************************/

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
