/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ttb.c                                          */
/*                                                                           */
/* Created:       2014/06/15 (TKa)                                           */
/* Last modified: 2017/04/24 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Thick-target bremsstrahlung approximation for electrons and  */
/*              positrons                                                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

/*****************************************************************************/

void TTByield(long mat, double Ee, long lflag, double *lEe, long mt, long *nbr,
              long *idxout, long id) {
  /* Samples the number of TTB photons emitted by an electron or positron.
   * */
  static char * const FUNCTION_NAME = "TTByield:";
  long ptd, ptr, nEe, idx;
  double Yn;
  const double *Eed, *lEed, *Ynd, *lYnd;
  static const double errtol = 1.e-13;   /* Tolerance for index error */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Check material cutoff energy TODO: -> CheckValue() */
  if (Ee < RDB[mat + MATERIAL_PHOTON_ECUT])
    Die(FUNCTION_NAME, "Electron energy %E below the cutoff %E (mt = %ld).",
       Ee, RDB[mat + MATERIAL_PHOTON_ECUT], mt);

  /* Check upper limit TODO: -> CheckValue() */
  if (Ee > RDB[DATA_PHOTON_EMAX])
    Die(FUNCTION_NAME, "Electron energy %E above the photon maximum "
        "%E (mt = %ld).", Ee, RDB[DATA_PHOTON_EMAX], mt);

  /* Check TTB mode */
  if ((long)RDB[DATA_PHOTON_USE_TTB] == NO)
    Die(FUNCTION_NAME, "The function was called although TTB is switched "
                        "off (mt = %ld).", mt);

  ptd = (long)RDB[mat + MATERIAL_PTR_TTB];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Energy grid size */
  nEe = (long)RDB[ptd + TTB_NE];

  /* Electron (and bremsstrahlung photon) energy grid */
  ptr = (long)RDB[ptd + TTB_E];
  CheckPointer(FUNCTION_NAME, "(E)", DATA_ARRAY, ptr);
  Eed = &RDB[ptr];

  /* Logarithmic energy grid */
  ptr = (long)RDB[ptd + TTB_LE];
  CheckPointer(FUNCTION_NAME, "(LE)", DATA_ARRAY, ptr);
  lEed = &RDB[ptr];

  /* Select positron or electron data */

  if (((long)RDB[DATA_PHOTON_TTBPM] == YES) &&
      (mt == MT_ELECTRON_PP_POS)) {

    /* Positron */

    ptr = (long)RDB[ptd + TTB_YP];
    CheckPointer(FUNCTION_NAME, "(YP)", DATA_ARRAY, ptr);
    Ynd = &RDB[ptr];

    ptr = (long)RDB[ptd + TTB_LYP];
    CheckPointer(FUNCTION_NAME, "(LYP)", DATA_ARRAY, ptr);
    lYnd = &RDB[ptr];
  }
  else {

    /* Electron */

    ptr = (long)RDB[ptd + TTB_YE];
    CheckPointer(FUNCTION_NAME, "(YE)", DATA_ARRAY, ptr);
    Ynd = &RDB[ptr];

    ptr = (long)RDB[ptd + TTB_LYE];
    CheckPointer(FUNCTION_NAME, "(LYE)", DATA_ARRAY, ptr);
    lYnd = &RDB[ptr];
  }

#ifdef DEBUG
  /* Check the log-energy */
  if (lflag && (fabs((*lEe - log(Ee))/log(Ee)) > errtol))
    Die(FUNCTION_NAME, "Log-energy argument %.15E differs from the true value"
                       " %.15E (mt = %ld).", *lEe, log(Ee), mt);
#endif

  /* Find energy interval */

  if ((Ee > RDB[DATA_PHOTON_TTBLINE]) || (lflag)) {

    if (!lflag)
      *lEe = log(Ee);

    /* Find log energy (NOTE: assuming log-interpolated energy array) */
    idx = (long)((*lEe - lEed[0])/(lEed[1] - lEed[0]));

    /* Check the index */
    if ((idx < 0) || (idx >= nEe))
      Die(FUNCTION_NAME, "Energy %E not found, index = %ld (mt = %ld)", Ee,
          idx, mt);

    if (idx == nEe - 1)
      idx--;

    /* Check energy limits, the index can be wrong due to floating point accuracy */
    if (Ee < Eed[idx]) {
      if ((Eed[idx] - Ee)/Eed[idx] < errtol)
        idx--;
      else
        Die(FUNCTION_NAME, "Energy %E not found in the interval [%E, %E], "
                           "index = %ld (mt = %ld)", Ee, Eed[idx], Eed[idx+1],
                           idx, mt);
    }
    else if (Ee > Eed[idx+1]) {
      if ((Ee - Eed[idx+1])/Eed[idx+1] < errtol)
        idx++;
      else
        Die(FUNCTION_NAME, "Energy %E not found in the interval [%E, %E], "
                           "index = %ld (mt = %ld)", Ee, Eed[idx], Eed[idx+1],
                           idx, mt);
    }
  }
  else {

    idx = SearchArray(Eed, Ee, (long)RDB[DATA_PHOTON_TTBLINEIDX]);

    if (idx == -1)
      Die(FUNCTION_NAME, "Energy %E not found, index = %ld (mt = %ld)", Ee,
          idx, mt);
  }

  /* Set the data index */
  *idxout = idx;

  /* Interpolate photon number yield */
  if (Ee <= RDB[DATA_PHOTON_TTBLINE]) {
    /* Linear interpolation */
    Yn = Ynd[idx] + (Ynd[idx+1] - Ynd[idx])*(Ee - Eed[idx])
        /(Eed[idx+1] - Eed[idx]);
  }
  else {
    /* Linear interpolation on log-log scale */
    Yn = exp(lYnd[idx] + (lYnd[idx+1] - lYnd[idx])*(*lEe - lEed[idx])
         /(lEed[idx+1] - lEed[idx]));
  }

  /* Sample number of photons.
   * NOTE: The number of photons must not be scored here because the material
   * cutoff energy is applied in TTBenergy(). */
  *nbr = (long)(Yn + RandF(id));

  /* Set zero number of photons at the lowest energy */
  if ((idx == 0) && (Ee == Eed[0]))
    *nbr = 0;

  /* Set log energy if photons are to be created */
  if ((*nbr > 0) && (Ee <= RDB[DATA_PHOTON_TTBLINE]) && !lflag)
    *lEe = log(Ee);

  CheckValue(FUNCTION_NAME, "nbr", "", (double)(*nbr), 0.0, INFTY);
}

/*****************************************************************************/


/*****************************************************************************/

void TTBenergy(long mat, long part, double Ee, double x, double y, double z,
               double u, double v, double w, double wgt, double t, long mt,
               long idx, long nbr, double lEe, double *Edep, long id) {
  /* Samples the TTB photon energies.
   * */
  static char * const FUNCTION_NAME = "TTBenergy:";
  long ptd, ptr, newp, brcdfptr, brpdfptr, intcptr, ncdf, i, nbrs;
  double  sumEk, rcdf, Ek, cdfmax;
  const double *Eed, *lEed, *brcdf, *brpdf, *intc;

  CheckValue(FUNCTION_NAME, "nbr", "", (double)nbr, 0, INFTY);
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

#ifdef DEBUG
  /* Check energy and TTB mode */
  if (Ee < RDB[mat + MATERIAL_PHOTON_ECUT])
    Die(FUNCTION_NAME, "Electron energy %E below the cutoff %E (mt = %ld)",
       Ee, RDB[mat + MATERIAL_PHOTON_ECUT], mt);

  if ((long)RDB[DATA_PHOTON_USE_TTB] == NO)
    Die(FUNCTION_NAME, "The function was called although TTB is switched "
            "off (mt = %ld)", mt);
#endif

  if (nbr == 0) {
    /* Electron/positron energy is deposited locally */
    *Edep = Ee;
    return;
  }

  ptd = (long)RDB[mat + MATERIAL_PTR_TTB];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Electron (and bremsstrahlung photon) energy grid */
  ptr = (long)RDB[ptd + TTB_E];
  CheckPointer(FUNCTION_NAME, "(E)", DATA_ARRAY, ptr);
  Eed = &RDB[ptr];

  /* Logarithmic energy grid */
  ptr = (long)RDB[ptd + TTB_LE];
  CheckPointer(FUNCTION_NAME, "(LE)", DATA_ARRAY, ptr);
  lEed = &RDB[ptr];

  /* Select positron or electron data */
  if (((long)RDB[DATA_PHOTON_TTBPM] == YES) &&
      (mt == MT_ELECTRON_PP_POS)) {

    /* Positron */

    brcdfptr = (long)RDB[ptd + TTB_BRPCDF];
    CheckPointer(FUNCTION_NAME, "(BRPCDF)", DATA_ARRAY, brcdfptr);

    brpdfptr = (long)RDB[ptd + TTB_BRPPDF];
    CheckPointer(FUNCTION_NAME, "(BRPPDF)", DATA_ARRAY, brpdfptr);

    intcptr = (long)RDB[ptd + TTB_INTCP];
    CheckPointer(FUNCTION_NAME, "(INTCP)", DATA_ARRAY, intcptr);
  }
  else {

    /* Electron */

    brcdfptr = (long)RDB[ptd + TTB_BRECDF];
    CheckPointer(FUNCTION_NAME, "(BRECDF)", DATA_ARRAY, brcdfptr);

    brpdfptr = (long)RDB[ptd + TTB_BREPDF];
    CheckPointer(FUNCTION_NAME, "(BREPDF)", DATA_ARRAY, brpdfptr);

    intcptr = (long)RDB[ptd + TTB_INTCE];
    CheckPointer(FUNCTION_NAME, "(INTCE)", DATA_ARRAY, intcptr);
  }

  /* Set bremsstrahlung energy cdf and pdf, the grid is selected using
   * interpolation weights */
  if ((lEed[idx] + RandF(id)*(lEed[idx+1] - lEed[idx]) <= lEe) || (idx == 0)) {
    /* Upper energy */
    brcdf = &RDB[(long)RDB[brcdfptr + idx + 1]];
    brpdf = &RDB[(long)RDB[brpdfptr + idx + 1]];
    intc = &RDB[(long)RDB[intcptr + idx + 1]];
    ncdf = idx + 2; /* +2 due to the index change */

    /* Interpolate maximum cdf (ProcessTTB checks that intc[idx] != -1) */
    cdfmax = brcdf[idx] + Eed[idx]*brpdf[idx]/(1.0 + intc[idx])*
             (exp((1.0 + intc[idx])*(lEe - lEed[idx])) - 1.0);
  }
  else {
    /* Lower energy (NOTE: the maximum photon energy will be below the
     * electron energy) */
    brcdf = &RDB[(long)RDB[brcdfptr + idx]];
    brpdf = &RDB[(long)RDB[brpdfptr + idx]];
    intc = &RDB[(long)RDB[intcptr + idx]];
    ncdf = idx + 1;
    cdfmax = brcdf[idx];
  }

  /***** Sample photon energies **********************************************/

  nbrs = 0;
  sumEk = 0.0;

  for (i = 0; i < nbr; i++) {

    rcdf = RandF(id)*cdfmax;
    idx = SearchArray(brcdf, rcdf, ncdf);

    /* Check index */
    if (idx == -1) {
      /* Check if the value is equal to the maximum */
      if (rcdf == brcdf[ncdf-1])
        idx = ncdf - 2;
      else
        Die(FUNCTION_NAME, "rcdf not found (Ee: %E cdfmax: %E rcdf: %E "
                       "brcdf[ncdf-1]: %E, mt = %ld)", Ee, cdfmax, rcdf,
                       brcdf[ncdf-1], mt);
    }

    /* Check limits to avoid numerical problems */
    if (rcdf == brcdf[idx]) {
      Ek = Eed[idx];
    }
    else if (rcdf == brcdf[idx+1]) {
      Ek = Eed[idx+1];
    }
    else {
      /* Calculate photon energy */
      Ek = Eed[idx]*pow((1.0 + intc[idx])*(rcdf - brcdf[idx])
                      /(brpdf[idx]*Eed[idx]) + 1.0, 1.0/(1.0 + intc[idx]));
    }

    /* Check photon energy */
    CheckValue(FUNCTION_NAME, "photon energy", "", Ek, Eed[idx], Eed[idx+1]);
    CheckValue(FUNCTION_NAME, "photon energy", "", Ek, RDB[DATA_PHOTON_EMIN], Ee);

    /* Apply material cutoff energy */
    if (Ek > RDB[mat + MATERIAL_PHOTON_ECUT]) {

      sumEk += Ek;
      nbrs++;

      /* Create new photon */
      newp = DuplicateParticle(part, id);

      WDB[newp + PARTICLE_X] = x;
      WDB[newp + PARTICLE_Y] = y;
      WDB[newp + PARTICLE_Z] = z;
      WDB[newp + PARTICLE_WGT] = wgt;
      WDB[newp + PARTICLE_T] = t;
      WDB[newp + PARTICLE_E] = Ek;
      WDB[newp + PARTICLE_PTR_MAT] = (double)mat;

      /* Approximation: photon direction is equal to the electron direction */
      WDB[newp + PARTICLE_U] = u;
      WDB[newp + PARTICLE_V] = v;
      WDB[newp + PARTICLE_W] = w;

      /* Set photon emission type */

      WDB[newp + PARTICLE_PHOTON_TYPE] = PHOTON_TYPE_TTB;

      /* Score particle balance (JLe 27.1.2017 / 2.1.28) */

      ptr = (long)RDB[RES_G_BALA_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_SRC_TTB, 0);
      AddBuf(wgt, 1.0, ptr, id, -1, BALA_G_SRC_TTB, 1);
      AddBuf(wgt*Ek, 1.0, ptr, id, -1, BALA_G_SRC_TTB, 2);

      ToQue(newp, id);
    }
  }

  /***************************************************************************/

  /* Set locally deposited energy */
  /* NOTE: The deposited energy can be negative at high electron energies.
   * However, the average deposited energy with large number of histories is
   * positive. */
  *Edep = Ee - sumEk;

}

/*****************************************************************************/


/*****************************************************************************/

void TTB(long mat, long part, double Ee, double x, double y, double z,
         double u, double v, double w, double wgt, double t, long mt,
         double *Edep, long id) {
  /* Wrapper function for TTB photon production combining TTByield() and
   * TTBenergy().
   * */
  double lEe;
  long nbr, idx;
  lEe = 0.0;
  nbr = idx = 0;

  /* Sample number of photons */
  TTByield(mat, Ee, 0, &lEe, mt, &nbr, &idx, id);

  /* Sample photon energies */
  TTBenergy(mat, part, Ee, x, y, z, u, v, w, wgt, t, mt,
       idx, nbr, lEe, Edep, id);
}

/*****************************************************************************/

