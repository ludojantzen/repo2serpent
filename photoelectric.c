/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : photoelectric.c                                */
/*                                                                           */
/* Created:       2011/04/15 (JLe)                                           */
/* Last modified: 2018/01/30 (TKa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Handles photoelectric effect for photons                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Photoelectric:"

/*****************************************************************************/

void Photoelectric(long mat, long rea, long part, double E, double x, double y,
                   double z, double u, double v, double w, double  wgt,
                   double  t, long id) {

  long ptd, ptr, Nssd, i, ss, ssmin, idx, NlEd, idxmin, idxE, szE, nbr,
      lflag, imax, idx_ar;
  double lE, Ee, lEe, beta, r1, ue, ve, we, mue, if0, cdfss, EdepAR,
      EdepTTB, Edep, a, k, gamma, fr_max, mue_max;
  const double *Ebd, *idxssd, *lpd, *lpss, *lEarrd, *pd, *cdfss_ar;

  /* Initialize values and avoid compiler warnings */
  NlEd = lflag = 0;
  Ee = EdepAR = EdepTTB = 0.0;
  Ebd = idxssd = lpd = lEarrd = pd = NULL;

  /* Pointer to photon reaction data */
  ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Size of the shell probability array */
  Nssd = (long)RDB[ptd + PHOTON_DIST_PE_NSS];

  if (Nssd > 0) {
    /* Binding energies */
    ptr = (long)RDB[ptd + PHOTON_DIST_PE_EBI];
    CheckPointer(FUNCTION_NAME, "(EBI)", DATA_ARRAY, ptr);
    Ebd = &RDB[ptr];
  }


  if (Nssd == 0 || E <= Ebd[Nssd-1]) {
    /* Photon energy is smaller than any ionization energy above
     * DATA_PHOTON_EMIN.
     * Photoelectron is assumed to obtain all the photon energy. */
    Ee = E;
  }
  else {

    /* Find the minimum subhell (K=0, L1=1, L2=2, etc) */
    for (ssmin = 0; ssmin < Nssd; ssmin++)
      if (E >= Ebd[ssmin])
        break;

    /* Check photon type and fluorescence index flag */

    if (((long)RDB[part + PARTICLE_PHOTON_TYPE] == PHOTON_TYPE_FLUOR)
            && ((idx_ar = (long)RDB[part + PARTICLE_PE_AR_IDX]) >= 0)) {
      /* Separate treatment for fluorescence photons:
       * Shell is sampled using the precomputed CDF. */

      /* Check pointers */
      CheckPointer(FUNCTION_NAME, "(AR_E)", DATA_ARRAY,
                   (long)RDB[ptd + PHOTON_DIST_PE_AR_E]);
      CheckPointer(FUNCTION_NAME, "(AR_E)", DATA_ARRAY,
                   (long)RDB[ptd + PHOTON_DIST_PE_AR_CDF]);
      CheckPointer(FUNCTION_NAME, "(AR_E)", DATA_ARRAY,
                   (long)RDB[(long)RDB[ptd + PHOTON_DIST_PE_AR_CDF] + idx_ar]);

      cdfss_ar = &RDB[(long)RDB[(long)RDB[ptd + PHOTON_DIST_PE_AR_CDF] + idx_ar]];

      /* Check the photon energy is equal to the table energy */
      CheckValue(FUNCTION_NAME, "Fluorescence energy", "", E,
                 RDB[(long)RDB[ptd + PHOTON_DIST_PE_AR_E] + idx_ar],
                 RDB[(long)RDB[ptd + PHOTON_DIST_PE_AR_E] + idx_ar]);

      /* Set log-energy */
      lE = RDB[part + PARTICLE_LE];

#ifdef DEBUG
      /* Check the log-energy */
      if (fabs(lE/log(E) - 1) > 1.e-11)
        Die(FUNCTION_NAME, "Stored AR log-energy %.15E differs from the true value %.15E", lE, log(E));
#endif

      /* Sample the shell */
      r1 = RandF(id);
      ss = -1;
      imax = Nssd - ssmin;

      for (i = 0; i < imax; i++) {
        if (r1 <= cdfss_ar[i]) {
          ss = i + ssmin;
          break;
        }

#ifdef DEBUG
        if ((cdfss_ar[i] < 0.0) || (cdfss_ar[i] > 1.0))
          Die(FUNCTION_NAME, "Shell CDF for fluorescence energies out of "
                             "bounds (value = %E)", cdfss_ar[i]);
        if ((i > 0) && (cdfss_ar[i] < cdfss_ar[i-1]))
          Die(FUNCTION_NAME, "Shell CDF for fluorescence energies is "
                             "not in ascending order.");
#endif
      }
    }
    else {
      /* Normal treatment */

      ptr = (long)RDB[ptd + PHOTON_DIST_PE_IDXSS];
      CheckPointer(FUNCTION_NAME, "(IDXSS)", DATA_ARRAY, ptr);
      idxssd = &RDB[ptr];

      ptr = (long)RDB[ptd + PHOTON_DIST_PE_LPSS];
      CheckPointer(FUNCTION_NAME, "(LPSS)", DATA_ARRAY, ptr);
      lpd = &RDB[ptr];

      NlEd = (long)RDB[ptd + PHOTON_DIST_PE_NE];

      ptr = (long)RDB[ptd + PHOTON_DIST_PE_E];
      CheckPointer(FUNCTION_NAME, "(E)", DATA_ARRAY, ptr);
      lEarrd = &RDB[ptr];

      lE = log(E);

      /* Find energy interval from the energy grid corresponding to the minimum
       * subshell index */
      idxmin = (long)idxssd[ssmin];

      if (ssmin > 0)
        szE = (long)idxssd[ssmin-1] - idxmin + 1;
      else
        szE = NlEd - idxmin;

      idxE = SearchArray(&lEarrd[idxmin], lE, szE) + idxmin;

      /* NOTE / TODO: This might fail if the energy is equal to the binding
       * energy (test this) */
      if (idxE < idxmin)
        Die(FUNCTION_NAME, "Energy not found in the subshell xs data");

      /* Interpolation factor */
      if0 = (lE - lEarrd[idxE])/(lEarrd[idxE+1] - lEarrd[idxE]);

      CheckValue(FUNCTION_NAME, "lE", "", lE, lEarrd[idxE], lEarrd[idxE+1]);

      /* Sample the shell number */

      cdfss = 0.0;
      ss = -1;
      r1 = RandF(id);

      for (i = ssmin; i < Nssd; i++) {

        /* Log-probability array of the shell */
        lpss = &RDB[(long)lpd[i]];

        /* Probability array index corresponding to the energy index */
        idx = idxE - (long)idxssd[i];

        /* Log-log interpolation of the shell probability */
        cdfss += exp(lpss[idx] + if0*(lpss[idx+1] - lpss[idx]));

        if (r1 <= cdfss) {
          ss = i;
          break;
        }
      }

      /* NOTE: This check might fail due to inaccuracies in the shell
       * cross sections. */
      CheckValue(FUNCTION_NAME, "cdfss", "\nThis check might fail due to "
                 "inaccuracies in the shell cross sections.", cdfss, 0.0, 1.0);
    }


    if (ss == -1) {
      /* Photoelectric effect occurred with an electron for which
       * ionization energy is smaller than DATA_PHOTON_EMIN.
       * Photoelectron is assumed to obtain all the photon energy. */
      Ee = E;
      lEe = lE;
      lflag = 1;
    }
    else {

      /* Set the electron energy */
      Ee = E - Ebd[ss];

      /* Atomic relaxation, store locally deposited energy */
      AtomicRelaxation(mat, rea, part, ss, x, y, z, wgt, t, &EdepAR, id);

      CheckValue(FUNCTION_NAME, "EdAR", "", EdepAR, -INFTY, Ebd[ss]);
    }
  }

  CheckValue(FUNCTION_NAME, "Ee", "", Ee, 0.0, E);


  if ((Ee < RDB[mat + MATERIAL_PHOTON_ECUT]) ||
      ((long)RDB[DATA_PHOTON_USE_TTB] == NO)) {
    /* Electron energy is deposited locally */
    Edep = Ee + EdepAR;
  }
  else {

    /* Sample the number of TTB photons */
    TTByield(mat, Ee, lflag, &lEe, MT_ELECTRON_PE, &nbr, &idx, id);


    if (nbr > 0) {

      /***** Sample the direction of the photoelectron ***********************/

      beta = sqrt(Ee*(Ee + 2.0*E_RESTMASS))/(Ee + E_RESTMASS);
      gamma = Ee/E_RESTMASS + 1.0;
      k = 0.5*gamma*(gamma - 1.0)*(gamma - 2.0);

      if (Ee < 0.5) {
        /* Use the low-energy sampling method */

        if (k < 0)
          fr_max = 1.0 + k*(1.0 - beta);
        else
          fr_max = 1.0 + k*(1.0 + beta);

        do {
          /* Rejection loop */
          do {
            r1 = RandF(id);
          } while (4.0*(1.0 - r1)*r1 < RandF(id));

          /* Cosine of the direction angle */
          mue = (2.0*r1 + beta - 1.0)/(2.0*beta*r1 - beta + 1.0);

        } while (RandF(id)*fr_max > 1.0 + k*(1.0 - beta*mue));

      }
      else {
        /* Use the high-energy sampling method */

        a = 1.0 + k - beta*k;
        mue_max = (2.0*beta - a)/(beta*(2.0 - a));

        if ((mue_max > 1.0) || (mue_max < -1.0))
            fr_max = 2.0*(1.0 + k*(1.0 + beta))/((1.0 + beta)*(1.0 + beta));
        else
            fr_max = a*a/(4.0*beta*(1.0 - beta));

        do {

          do {
            /* Rejection sampling loop */

            r1 = RandF(id);
            mue = (2.0*r1 + beta - 1.0)/(2.0*beta*r1 - beta + 1.0);

          } while (RandF(id)*2.0 > 1.0 + mue);

          a = 1.0 - beta*mue;

        } while (RandF(id)*fr_max*a*a > (1.0 - mue)*(1.0 + k*a));

      }

      CheckValue(FUNCTION_NAME, "mue", "", mue, -1.0, 1.0);

      /* Electron direction cosines */
      ue = u;
      ve = v;
      we = w;

      /* Sanity check for mu and direction vectors (for NAN's etc.) */
      CheckValue(FUNCTION_NAME, "mue", "", mue, -1.01, 1.01);
      CheckValue(FUNCTION_NAME, "ue", "", ue, -1.01, 1.01);
      CheckValue(FUNCTION_NAME, "ve", "", ve, -1.01, 1.01);
      CheckValue(FUNCTION_NAME, "we", "", we, -1.01, 1.01);

      /* Rotate direction cosines */
      AziRot(mue, &ue, &ve, &we, id);

      /***********************************************************************/

      /* Use the TTB-approximation for the electron */
      TTBenergy(mat, part, Ee, x, y, z, ue, ve, we, wgt, t, MT_ELECTRON_PE,
           idx, nbr, lEe, &EdepTTB, id);

      CheckValue(FUNCTION_NAME, "EdTTB", "", EdepTTB, -INFTY, Ee);

      /* Total deposited energy */
      Edep = EdepTTB + EdepAR;

    }
    else {
      /* No bremsstrahlung photons were created, the whole electron energy is
       * deposited locally */
      Edep = Ee + EdepAR;
    }

  }


  /* Put particle back in stack */

  ToStack(part, id);

  /* Score pulse-height detector */

  PulseDet(part, mat, Edep, x, y, z, wgt, id);

  /* Score photon heat detector */

  ScorePhotonHeat(part, mat, Edep, x, y, z, E, t, wgt, id);

}

/*****************************************************************************/
