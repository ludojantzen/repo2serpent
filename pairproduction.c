/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : pairproduction.c                               */
/*                                                                           */
/* Created:       2011/04/15 (JLe)                                           */
/* Last modified: 2018/01/30 (TKa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Handles pair production of photons                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

/* Local function definitions */
static double PairProductionSGHmu(double, long);


/*****************************************************************************/

void PairProduction(long mat, long rea, long part, double E, double x,
                    double y, double z, double u0, double v0, double w0,
                    double wgt, double t, long id) {
  static char * const FUNCTION_NAME = "PairProduction:";
  long ptd, ptr;
  double E2, lE, lE2, eps, eps2, epsmin, epsmax, tepsmin, stepsmin, xseps,
      xsepsmax, F0, G0, G, G2, F0fc, phi1, phi2, epsE, Ee, Ep, mue, mup, ue,
      ve, we, up, vp, wp, a, b, Edepe, Edepp, Ed;
  const double *F0c, *mdxsfc;

  /***************************************************************************/

  /* Check reaction pointer */
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to photon reaction data */
  ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Coefficient array F0 */
  ptr = (long)RDB[ptd + PHOTON_DIST_PP_F0];
  CheckPointer(FUNCTION_NAME, "(F0)", DATA_ARRAY, ptr);
  F0c = &RDB[ptr];

  /* Fit coefficients for the maximum differential xs */
  if (E <= 200.0) {
    ptr = (long)RDB[ptd + PHOTON_DIST_PP_MDXSFC1];
    CheckPointer(FUNCTION_NAME, "(MDXSFC1)", DATA_ARRAY, ptr);
    mdxsfc = &RDB[ptr];
  }
  else {
    ptr = (long)RDB[ptd + PHOTON_DIST_PP_MDXSFC2];
    CheckPointer(FUNCTION_NAME, "(MDXSFC2)", DATA_ARRAY, ptr);
    mdxsfc = &RDB[ptr];
  }

  /* Check and correct photon energy if needed (E can be below 2*E_RESTMASS
   * if different electron rest mass is used in xs data)
   * TODO: Edep-korjaus */
  if (E < 2.0*E_RESTMASS) {
    Warn(FUNCTION_NAME, "Photon energy %.10E is below pair production threshold "
                        "energy\n%.10E. Photon energy is set equal to the threshold.",
         E, 2.0*E_RESTMASS);
    E = 2.0*E_RESTMASS;
  }

  /***************************************************************************/


  if ((long)RDB[DATA_PHOTON_USE_TTB] == YES) {

    /***** Sample the electron and positron energy ***************************/

    eps = 0.0;
    epsmin = E_RESTMASS/E;
    epsmax = 1.0 - epsmin;

    /* Sample the electron reduced energy eps */

    if (E < 2.0) {
      /* Sample eps from uniform distribution */
      eps = epsmin + RandF(id)*(epsmax - epsmin);
    }
    else if (E <= 1.0e6) {
      /* Sample eps from the differential xs */

      tepsmin = 2.0*epsmin;
      stepsmin = sqrt(tepsmin);

      F0 = F0c[0]*stepsmin + tepsmin*(F0c[1] + F0c[2]*stepsmin
           + F0c[3]*tepsmin);
      F0fc = F0 + RDB[ptd + PHOTON_DIST_PP_FC];
      G0 = RDB[ptd + PHOTON_DIST_PP_G0]*epsmin;

      /* Rational function fit for the maximum of the differential xs */
      if (E <= 200.0) {
        E2 = E*E;
        xsepsmax = (mdxsfc[0]*E2 + mdxsfc[1]*E + mdxsfc[2])
                   /(E2 + mdxsfc[3]*E + mdxsfc[4]);
      }
      else {
        lE = log(E);
        lE2 = lE*lE;
        xsepsmax = (mdxsfc[0]*lE2 + mdxsfc[1]*lE + mdxsfc[2])
                   /(lE2 + mdxsfc[3]*lE + mdxsfc[4]);
      }

      /* Rejection sampling loop */
      do {

        /* Sample from a uniform distribution */
        eps = epsmin + RandF(id)*(epsmax - epsmin);
        eps2 = eps*eps;
        G = G0/(eps - eps2);

        if (G < 1.0) {
          G2 = G*G;
          phi1 = 20.867 - 3.242*G + 0.625*G2;
          phi2 = 20.209 - 1.930*G - 0.086*G2;
          xseps = (2.0*(eps2 - eps) + 1.0)*(phi1 + F0fc) +
              2.0/3.0*(eps - eps2)*(phi2 + F0fc);
        }
        else {
          phi1 = 21.12 - 4.184*log(G + 0.952);
          xseps = (4.0/3.0*(eps2 - eps) + 1.0)*(phi1 + F0fc);
        }

      } while (xsepsmax*RandF(id) > xseps);

    }
    else {
      Die(FUNCTION_NAME, "Photon energy exceeds the pair production limit");
    }

    CheckValue(FUNCTION_NAME, "eps", "", eps, 0, 1.0);

    /* Calculate the electron and positron kinetic energies */
    epsE = eps*E;
    Ee = epsE - E_RESTMASS;
    Ep = E - epsE - E_RESTMASS;

    /* Check energies */
    if (Ee < 0.0)
      Die(FUNCTION_NAME, "Negative electron energy %E (photon energy : %.15E)",
          Ee, E);
    if (Ep < 0.0)
      Die(FUNCTION_NAME, "Negative positron energy %E (photon energy : %.15E)",
          Ep, E);

    /*************************************************************************/

    /* Sample electron direction */
    mue = PairProductionSGHmu(Ee, id);

    /* Initialize electron direction cosines */
    ue = u0;
    ve = v0;
    we = w0;

    /* Sanity check for mu and direction vectors (for NAN's etc.) */
    CheckValue(FUNCTION_NAME, "mue", "", mue, -1.0, 1.0);
    CheckValue(FUNCTION_NAME, "ue", "", ue, -1.01, 1.01);
    CheckValue(FUNCTION_NAME, "ve", "", ve, -1.01, 1.01);
    CheckValue(FUNCTION_NAME, "we", "", we, -1.01, 1.01);

    /* Rotate */
    AziRot(mue, &ue, &ve, &we, id);

    /* Sample positron direction */
    mup = PairProductionSGHmu(Ep, id);

    CheckValue(FUNCTION_NAME, "mup", "", mup, -1.0, 1.0);

    /* Calculate positron direction cosines */
    a = sqrt((1.0 - mup*mup)/(1.0 - mue*mue));
    b = mup + a*mue;
    up = b*u0 - a*ue;
    vp = b*v0 - a*ve;
    wp = b*w0 - a*we;

    CheckValue(FUNCTION_NAME, "positron direction cosines", "",
               up*up + vp*vp + wp*wp - 1.0, -1E-4, 1E-4);

    /* TTB approximation for the electron */
    if (Ee > RDB[mat + MATERIAL_PHOTON_ECUT])
      TTB(mat, part, Ee, x, y, z, ue, ve, we, wgt, t, MT_ELECTRON_PP_EL,
          &Edepe, id);
    else
      Edepe = Ee;

    /* TTB approximation for the positron */
    if (Ep > RDB[mat + MATERIAL_PHOTON_ECUT])
      TTB(mat, part, Ep, x, y, z, up, vp, wp, wgt, t, MT_ELECTRON_PP_POS,
          &Edepp, id);
    else
      Edepp = Ep;

    CheckValue(FUNCTION_NAME, "Edepe", "", Edepe, -INFTY, Ee);
    CheckValue(FUNCTION_NAME, "Edepp", "", Edepp, -INFTY, Ep);

  }
  else {
    /* Electron and positron energy are deposited locally, set Ede and Edp
     * for convenience */
    Edepe = E - 2.0*E_RESTMASS;
    Edepp = 0.0;

    /* Set direction cosines to zero (they are not used) */
    vp = up = wp = 0.0;

    /* Set zero positron energy (not used for anything) */
    Ep = 0.0;
  }

  /* Positron annihilation (TODO: katkaisu PosAnnihilation():iin) */
  if (RDB[mat + MATERIAL_PHOTON_ECUT] <= E_RESTMASS) {
    PosAnnihilation(mat, part, Ep, x, y, z, up, vp, wp, wgt, t, id);
  }
  else {
    /* Add annihilation photons to deposited energy */
    Edepp += 2.0*E_RESTMASS;
  }

  /* Put incident photon back in stack */
  ToStack(part, id);

  /* Total deposited energy */
  Ed = Edepe + Edepp;

  /* Score pulse-height detector */

  PulseDet(part, mat, Ed, x, y, z, wgt, id);

  /* Score photon heat detector */

  ScorePhotonHeat(part, mat, Ed, x, y, z, E, t, wgt, id);

}

/*****************************************************************************/

static double PairProductionSGHmu(double E, long id) {
  /* Samples the cosine of the polar emission angle of an electron (or
   * positron) from the leading order term of the Sauter-Gluckstern-Hull
   * distribution.
   * */
  double r, beta, mu;

  beta = sqrt(E*(E + 2.0*E_RESTMASS))/(E + E_RESTMASS);
  do {
    r = RandF(id);
    mu = (2.0*r + beta - 1.0)/(2.0*beta*r - beta + 1.0);
  } while (1.0 - mu*mu == 0.0);

  return mu;
}

/*****************************************************************************/

