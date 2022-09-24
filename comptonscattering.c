/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : comptonscattering.c                            */
/*                                                                           */
/* Created:       2011/04/15 (JLe)                                           */
/* Last modified: 2018/01/30 (TKa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Handles incoherent Compton scattering of photons             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ComptonScattering:"

/* TODO:
 * - High-energy Compton? Is Klein-Nishina (with or without scattering
 *   function) adequate then?
 * */

/* Local function definitions */
void PartialPivotGauss(double[3][3], double[3], double[3]);


/*****************************************************************************/

void ComptonScattering(long mat, long rea, long part, double *E, double x,
                       double y, double z, double *u, double *v, double *w,
                       double wgt, double t, long id) {

  long ptd, ptr, Nxtd, Nssd, lo, lomax, Npzd, i, ss, eli, Nssar,
      ptrelncdfar, pznegative, pzminidx, ssmin, nbr, idx;
  double mu, Ek, Ek0, Emax, Ee0, Ee, lSxt, xtmax, lSxtmax, xt,
      pzmaxd, pzmax, tmp, cpintmax, cpint, pzmaxabs, Ebss, reln, a1, a2,
      pz, pz2, rcpint, Ediscr, mue, mue0, ue0, ve0, we0, ue, ve, we, uq, vq,
      wq, u0, v0, w0, Edep, EdepTTB, EdepAR, q, pe0c, pec, cbeta, cdelta,
      calpha, nu, pz2nu2, lEe, cpi_sgn;
  const double *xtd, *lxtd, *lSxtd, *elnd, *Ebd, *cpd, *cpintd, *cpad, *pzd,
      *cpssd, *cpintssd, *cpassd, *extad, *elncdf, *ebidar, *extbd, *extid,
      *cpintmind;
  static const double zeta = 57.0320106155672; /* 1e6/(sqrt(2)*h*c)*angstrom */
  static const double fscerm = FS_CONST*E_RESTMASS;
  static const double warntol = 1.e-10;      /* Warning tolerance for Doppler broadened energy */
  double cdopss[PHOTON_NSS_MAX][4];
  double Auvw[3][3];
  double buvw[3], uvw[3];

  /***************************************************************************/

  /* Check reaction pointer */
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to photon reaction data */
  ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Size of the incoherent scattering function data */
  Nxtd = (long)RDB[ptd + PHOTON_DIST_COMP_NISF];

  /* Momentum transfers (x = sin(theta/2)/lambda) */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_ISFX];
  CheckPointer(FUNCTION_NAME, "(ISFX)", DATA_ARRAY, ptr);
  xtd = &RDB[ptr];

  /* Log-momentum transfers */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_LISFX];
  CheckPointer(FUNCTION_NAME, "(LISFX)", DATA_ARRAY, ptr);
  lxtd = &RDB[ptr];

  /* Incoherent scattering functions */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_LISF];
  CheckPointer(FUNCTION_NAME, "(LISF)", DATA_ARRAY, ptr);
  lSxtd = &RDB[ptr];

  /* Size of Compton profile array */
  Npzd = (long)RDB[ptd + PHOTON_DIST_COMP_NCP];

  /* Index of the minimum pz in the data */
  pzminidx = (long)RDB[ptd + PHOTON_DIST_COMP_PZMINIDX];

  /* pz (projection of the intial electron momentum on momentum transfer) */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_CPPZ];
  CheckPointer(FUNCTION_NAME, "(CPPZ)", DATA_ARRAY, ptr);
  pzd = &RDB[ptr];

  /* Compton profiles */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_CP];
  CheckPointer(FUNCTION_NAME, "(CP)", DATA_ARRAY, ptr);
  cpd = &RDB[ptr];

  /* Integrated Compton profiles */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_CPINT];
  CheckPointer(FUNCTION_NAME, "(CPINT)", DATA_ARRAY, ptr);
  cpintd = &RDB[ptr];

  /* Interpolation coefficients */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_CPA];
  CheckPointer(FUNCTION_NAME, "(CPA)", DATA_ARRAY, ptr);
  cpad = &RDB[ptr];

  /* Extrapolation coefficients */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_CPEXTA];
  CheckPointer(FUNCTION_NAME, "(CPEXTA)", DATA_ARRAY, ptr);
  extad = &RDB[ptr];

  if (RDB[DATA_PHOTON_CP_FALLBACK_2129] == (double)NO) {

    /* The second extrapolation coefficients */
    ptr = (long)RDB[ptd + PHOTON_DIST_COMP_CPEXTB];
    CheckPointer(FUNCTION_NAME, "(CPEXTB)", DATA_ARRAY, ptr);
    extbd = &RDB[ptr];

    /* Integration coefficients */
    ptr = (long)RDB[ptd + PHOTON_DIST_COMP_CPEXTI];
    CheckPointer(FUNCTION_NAME, "(CPEXTI)", DATA_ARRAY, ptr);
    extid = &RDB[ptr];

    /* TODO: tän vois ehkä nimetä paremmin? */
    /* Compton profile integrals between 0 and ~137 */
    ptr = (long)RDB[ptd + PHOTON_DIST_COMP_CPINTMIN];
    CheckPointer(FUNCTION_NAME, "(CPINTMIN)", DATA_ARRAY, ptr);
    cpintmind = &RDB[ptr];
  }
  else {
    extbd = NULL;
    extid = NULL;
    cpintmind = NULL;
  }

  /* Number of subshells */
  Nssd = (long)RDB[ptd + PHOTON_DIST_COMP_NSS];
  CheckValue(FUNCTION_NAME, "NSS", "", (double)Nssd, 1, PHOTON_NSS_MAX);

  /* Electron binding energies */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_EBI];
  CheckPointer(FUNCTION_NAME, "(EBI)", DATA_ARRAY, ptr);
  Ebd = &RDB[ptr];

  /* Number of electrons per subshell */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_ELN];
  CheckPointer(FUNCTION_NAME, "(ELN)", DATA_ARRAY, ptr);
  elnd = &RDB[ptr];

  /* Electron shell number cdf */
  ptr = (long)RDB[ptd + PHOTON_DIST_COMP_ELNCDF];
  CheckPointer(FUNCTION_NAME, "(ELNCDF)", DATA_ARRAY, ptr);
  elncdf = &RDB[ptr];

  /***************************************************************************/

  /***************************************************************************/

  /* Store initial photon energy */
  Ek0 = *E;

  /* Total energy of the electron before the scattering (approximately) */
  Ee0 = E_RESTMASS;

  /* The binding energy of the electron (zero when Doppler broadening is not
   * used) */
  Ebss = 0.0;

  if (RDB[DATA_PHOTON_CP_FALLBACK_2129] == (double)YES) {
    /* Check the minimum pz in the data */
    CheckValue(FUNCTION_NAME, "minimum pz", "", pzd[pzminidx], Ee0/fscerm,
               Ee0/fscerm);
  }

  /* Avoid compiler warnings */
  ss = 0;
  pz = cpint = 0.0;

  if (Ek0 > RDB[DATA_PHOTON_EKN]) {
    /* The photon energy and scattering cosine are sampled from Klein-Nishina
     * formula. */
    KleinNishina(Ek0, &Ek, &mu, id);
  }

  else {

    /* Use incoherent scattering function with Klein-Nishina formula */

    /*************************************************************************/
    lomax = -10;
    xtmax = zeta*Ek0*SQRT2;

    /* Calculate the maximum incoherent scattering function */

    if ((xtmax < xtd[1]) && ((lxtd[0] == 0.0) || (lSxtd[0] == 0.0))) {
      /* Lin-lin linear interpolation, when the first element is zero */
      /* NOTE: The first elements are stored normally, the second ones as
       * log-log */
      lSxtmax = log(lSxtd[0] + (exp(lSxtd[1]) - lSxtd[0])*(xtmax - lxtd[0])
               /(exp(lxtd[1]) - lxtd[0]));
    }
    else if (xtmax >= xtd[Nxtd-1]) {
      lSxtmax = lxtd[Nxtd-1];
      lomax = Nxtd-2;
    }
    else {
      /* Find the lower boundary, the first non-log element is excluded */
      lomax = SearchArray(&xtd[1], xtmax, Nxtd - 1) + 1;

      if (lomax == 0)
        Die(FUNCTION_NAME, "xtmax = %.5E not found", xtmax);

      /* Log-log linear interpolation (log maximum) */
      lSxtmax = lSxtd[lomax] + (lSxtd[lomax+1] - lSxtd[lomax])
                *(log(xtmax) - lxtd[lomax])/(lxtd[lomax+1] - lxtd[lomax]);
    }

    /* Rejection loop */

    do {
      /* Sample scattering cosine from Klein-Nishina formula */
      KleinNishina(Ek0, &Ek, &mu, id);

      /* Check scattering cosine */
      CheckValue(FUNCTION_NAME, "mu", "", mu, -1.0, 1.0);

      /***** Use the incoherent form factor as a rejection function **/
      xt = zeta*Ek0*sqrt(1.0 - mu);

      if ((xt < xtd[1]) && ((lxtd[0] == 0.0) || (lSxtd[0] == 0.0))) {
        lSxt = log(lSxtd[0] + (exp(lSxtd[1]) - lSxtd[0])*(xt - lxtd[0])
                 /(exp(lxtd[1]) - lxtd[0]));
      }
      else {
        /* Find the lower boundary, exclude the first non-log element */
        lo = SearchArray(&xtd[1], xt, lomax + 1) + 1;

        if (lo == 0)
          Die(FUNCTION_NAME, "xt = %E not found", xt);

        CheckValue(FUNCTION_NAME, "xt", "", xt, xtd[lo], xtd[lo+1]);

        /* Interpolate the scattering function (log) */
        lSxt = lSxtd[lo] + (lSxtd[lo+1] - lSxtd[lo])
               *(log(xt) - lxtd[lo])/(lxtd[lo+1] - lxtd[lo]);
      }

    } while (RandF(id) > exp(lSxt - lSxtmax));

    nu = Ek0/Ek;

    /*************************************************************************/


    /*************************************************************************/

    if ((long)RDB[DATA_PHOTON_USE_DOPPLER]) {

      /* Doppler broadening of the photon energy, i.e., the initial momentum
       * and the binding energy of the electron are taken into account */

      /* Set variables */
      a1 = 1.0 - mu*mu;
      a2 = 1.0 + nu*(nu - 2.0*mu);
      memset(cdopss, 0, PHOTON_NSS_MAX*4*sizeof(double));

      /* Maximum pz in the data */
      pzmaxd = pzd[Npzd - 1];

      /* Find the minimum shell number */
      for (ssmin = 0; ssmin < Nssd; ssmin++)
        if (Ebd[ssmin] < Ek0)
          break;

      if (RDB[DATA_PHOTON_CP_FALLBACK_2129] == (double)YES) {

        /* Use Compton profile model from 2.1.29 (will be removed after 2.1.30) */

      /***********************************************************************/

      while (1) {

        /***** Sample electron shell *****************************************/

        while (1) {

          /* Sample the shell according to the number of electrons per shell */
          reln = RandF(id)*(elncdf[Nssd] - elncdf[ssmin]) + elncdf[ssmin];
          ss = SearchArray(&elncdf[ssmin], reln, Nssd + 1 - ssmin) + ssmin;

          CheckValue(FUNCTION_NAME, "ss", "", (double)ss, (double)ssmin,
                     (double)Npzd - 1.0);
          CheckValue(FUNCTION_NAME, "ss", "PHOTON_NSS_MAX", (double)ss,
                     0, PHOTON_NSS_MAX);

          /* Set subshell data */
          cpssd = &RDB[(long)cpd[ss]];
          cpintssd = &RDB[(long)cpintd[ss]];
          cpassd = &RDB[(long)cpad[ss]];

          if (!(long)cdopss[ss][0]) {

            /* Calculate the integral */

            /* Shell data flag */
            cdopss[ss][0] = 1.0;

            /* Calculate the maximum pz of the subshell */
            tmp = Ek0*(Ek0 - Ebd[ss])*(1.0 - mu);
            pzmax = (tmp - Ebd[ss]*Ee0)/(sqrt(2.0*tmp + Ebd[ss]*Ebd[ss])
                                         *fscerm);
            cdopss[ss][1] = pzmax;

            /* Calculate the integral between [0, pzmax] */
            if (pzmax >= pzmaxd) {
              /* Integral of the exponential extrapolation */
              cpintmax = cpintssd[Npzd-1] + cpssd[Npzd-1]/extad[ss]
                         *(exp(extad[ss]*(pzmax - pzmaxd)) - 1.0);
            }
            else {
              /* Linear interpolation */
              pzmaxabs = fabs(pzmax);

              /* Get the lover boundary of pzmax */
              if (pzmaxabs == pzmaxd)
                lo = Npzd - 2;
              else
                lo = SearchArray(pzd, pzmaxabs, Npzd);

              CheckValue(FUNCTION_NAME, "lo", "", (double)lo, 0,
                         (double)Npzd - 2.0);

              /* Calculate the integral */
              cpintmax = cpintssd[lo] + (0.5*cpassd[lo]*(pzmaxabs - pzd[lo])
                         + cpssd[lo])*(pzmaxabs - pzd[lo]);
            }

            /* Calculate the integral between pzmin and pzmax */
            if (pzmax < 0.0)
              cpint = cpintssd[pzminidx] - cpintmax;
            else
              cpint = cpintssd[pzminidx] + cpintmax;

            cdopss[ss][2] = cpintmax;
            cdopss[ss][3] = cpint;

            CheckValue(FUNCTION_NAME, "cpint", "", cpint, 0.0, 1.0);

            /* Sampling not needed for single shell atoms */
            if (Nssd == 1)
              break;

          }
          else {
            /* The integral was already calculated */

            pzmax = cdopss[ss][1];
            cpintmax = cdopss[ss][2];
            cpint = cdopss[ss][3];
          }

          /* Accept or reject the shell */
          if (RandF(id) <= cpint)
            break;
        }

        /* Binding energy of the sampled shell */
        Ebss = Ebd[ss];

        /* Maximum scattered photon energy */
        Emax = Ek0 - Ebss;

        /*********************************************************************/


        /***** Sample pz and solve the corresponding photon energy ***********/

        while (1) {

          /***** Sample pz from the Compton profile of the sampled subshell **/

          pznegative = 0;

          if (pzmax > 0.0) {
            rcpint = RandF(id)*cpint;
            if (rcpint < cpintssd[pzminidx]) {
              /* pz between [pzmin, 0] */
              pznegative = 1;
            }
            else {
              /* pz between [0, pzmax] */
              rcpint = rcpint - cpintssd[pzminidx];
              pznegative = 0;
            }
          }
          else {
            /* pz between [pzmin, pzmax], both are negative */
            rcpint = cpintmax + RandF(id)*cpint;
            pznegative = 1;
          }

          if (rcpint > cpintssd[Npzd-1]) {
            /* Tail */
            lo = Npzd - 1;

            if (pznegative)
              Die(FUNCTION_NAME, "negative pz");

            pz = pzmaxd + log(extad[ss]*(rcpint - cpintssd[lo])
                              /cpssd[lo] + 1.0)/extad[ss];

            CheckValue(FUNCTION_NAME, "pz", "", pz, pzmaxd, INFTY);
          }
          else {
            /* Get the lower boundary */
            lo = SearchArray(cpintssd, rcpint, Npzd);

            CheckValue(FUNCTION_NAME, "lo", ": rcpint not found", (double)lo, 0,
                       (double)Npzd - 2.0);

            /* Calculate pz */
            if (cpassd[lo] == 0.0) {
              /* The slope is zero (note: cpss[lo] is always greater than zero)*/
              pz = pzd[lo] + (rcpint - cpintssd[lo]) / cpssd[lo];
            }
            else {
              /* Use the integral of the linearly interpolated pdf */
              pz = pzd[lo] + (sqrt(cpssd[lo]*cpssd[lo] + 2.0*cpassd[lo]
                   *(rcpint - cpintssd[lo])) - cpssd[lo])/cpassd[lo];
            }

            CheckValue(FUNCTION_NAME, "pz", "", pz, pzd[lo], pzd[lo+1]);
          }

          /* Set the sign */
          if (pznegative)
            pz = -pz;

          /*******************************************************************/

          /***** Solve the energy ********************************************/

          pz *= FS_CONST;
          pz2 = pz*pz;
          Ediscr = a2 - pz2*a1;
          pz2nu2 = pz2 - nu*nu;

          if ((Ediscr < 0.0) || (pz2nu2 == 0.0)) {
            /* Sample new pz to avoid negative discriminant or division by zero.
             * Ediscr can be negative due to floating point accuracy. These are
             * rare events. */
            Warn(FUNCTION_NAME, "Energy sampling failed (Ediscr: %E, pz2nu2: "
                                "%E)\n Sampling new energy...", Ediscr, pz2nu2);
            continue;
          }

          /* Photon energy */
          Ek = Ek0*(pz2*mu - nu - pz*sqrt(Ediscr))/pz2nu2;

          /* Check energy */
          if (Ek > Emax) {
            if ((Ek - Emax)/Emax < warntol) {
              /* Photon energy can be above the maximum due to floating point
               * accuracy. */
              continue;
            }
            else {
              Warn(FUNCTION_NAME, "Sampled photon energy %.15E above the "
                                  "maximum %.15E.\n Sampling new energy... "
                                  "\n %f %ld",
                                   Ek, Emax, pzmax, pznegative);
              continue;
            }
          }
          else if (Ek < 0.0) {
            Warn(FUNCTION_NAME, "Sampled photon energy %.15E below zero.\n"
                                "Sampling new energy... ", Ek);
            continue;
          }

#ifdef DEBUG
          if (isnan(Ek))
            Die(FUNCTION_NAME, "NaN photon energy in Doppler broadening.");
#endif

          /* Exit the pz loop */
          break;

          /*******************************************************************/
        }

        /*********************************************************************/

        /* Rejection sampling (approximation to Ribberfors DDCS) */
        if (RandF(id)*Ek0 <= Ek)
          break; 
      }

      /***********************************************************************/

      }
      else  {

      /***********************************************************************/


      while (1) {

        /***** Sample electron shell *****************************************/

        while (1) {

          /* Sample the shell according to the number of electrons per shell */
          reln = RandF(id)*(elncdf[Nssd] - elncdf[ssmin]) + elncdf[ssmin];
          ss = SearchArray(&elncdf[ssmin], reln, Nssd + 1 - ssmin) + ssmin;

          CheckValue(FUNCTION_NAME, "ss", "", (double)ss, (double)ssmin,
                     (double)Npzd - 1.0);
          CheckValue(FUNCTION_NAME, "ss", "PHOTON_NSS_MAX", (double)ss,
                     0, PHOTON_NSS_MAX);

          /* Set subshell data */
          cpssd = &RDB[(long)cpd[ss]];
          cpintssd = &RDB[(long)cpintd[ss]];
          cpassd = &RDB[(long)cpad[ss]];

          if (!(long)cdopss[ss][0]) {

            /* Calculate the integral */

            /* Shell data flag */
            cdopss[ss][0] = 1.0;

            /* Calculate the maximum pz of the subshell */
            tmp = Ek0*(Ek0 - Ebd[ss])*(1.0 - mu);
            pzmax = (tmp - Ebd[ss]*Ee0)/(sqrt(2.0*tmp + Ebd[ss]*Ebd[ss])
                                         *fscerm);
            cdopss[ss][1] = pzmax;
            pzmaxabs = fabs(pzmax);

            /* Calculate the integral between [0, pzmax] */
            if (pzmaxabs >= pzmaxd) {

              /* Integral of the exponential extrapolation */

              tmp = 1.0 + extbd[ss]*pzmaxabs;
              cpintmax = 0.5*(1.0 - extad[ss]/extbd[ss]*exp(-tmp*tmp));
            }
            else {
              /* Linear interpolation */

              /* Get the lover boundary of pzmax */
              if (pzmaxabs == pzmaxd)
                lo = Npzd - 2;
              else
                lo = SearchArray(pzd, pzmaxabs, Npzd);

              CheckValue(FUNCTION_NAME, "lo", "", (double)lo, 0,
                         (double)Npzd - 2.0);

              /* Calculate the integral */
              cpintmax = cpintssd[lo] + (0.5*cpassd[lo]*(pzmaxabs - pzd[lo])
                         + cpssd[lo])*(pzmaxabs - pzd[lo]);
            }

            /* Calculate the integral between pzmin and pzmax */
            if (pzmax < 0.0)
              cpint = cpintmind[ss] - cpintmax;
            else
              cpint = cpintmind[ss] + cpintmax;

            cdopss[ss][2] = cpintmax;
            cdopss[ss][3] = cpint;

            CheckValue(FUNCTION_NAME, "cpint", "", cpint, 0.0, 1.0);

            /* Sampling not needed for single shell atoms */
            if (Nssd == 1)
              break;

          }
          else {
            /* The integral was already calculated */
            pzmax = cdopss[ss][1];
            cpintmax = cdopss[ss][2];
            cpint = cdopss[ss][3];
            pzmaxabs = fabs(pzmax);
          }

          /* Accept or reject the shell */
          if (RandF(id) <= cpint)
            break;
        }

        /* Binding energy of the sampled shell */
        Ebss = Ebd[ss];

        /* Maximum scattered photon energy */
        Emax = Ek0 - Ebss;

        /*********************************************************************/


        /***** Sample pz and solve the corresponding photon energy ***********/

        while (1) {

          /***** Sample pz from the Compton profile of the sampled subshell **/

          pznegative = 0;

          if (pzmax > 0.0) {
            rcpint = RandF(id)*cpint;

            if (rcpint < cpintmind[ss]) {
              /* pz between [pzmin, 0] */

              rcpint = cpintmind[ss] - rcpint;
              pznegative = 1;
            }
            else {
              /* pz between [0, pzmax] */
              rcpint = rcpint - cpintmind[ss];
              pznegative = 0;
            }
          }
          else {
            /* pz between [pzmin, pzmax], both are negative */
            rcpint = cpintmax + RandF(id)*cpint;
            pznegative = 1;
          }


          if (rcpint > cpintssd[Npzd-1]) {

            /* Extrapolation of the Compton profile tail */

            lo = Npzd - 1;

            pz = (sqrt(-log(extbd[ss]/extad[ss]*(1.0 - 2.0*rcpint))) - 1.0)
                /extbd[ss];

            if (pznegative)
              CheckValue(FUNCTION_NAME, "pz", "(negative extrapolation)",
                         pz, pzmaxd, 1.0/FS_CONST);
            else
              CheckValue(FUNCTION_NAME, "pz", "(positive extrapolation)",
                         pz, pzmaxd, pzmaxabs);
          }
          else {

            /* Get the lower boundary */
            lo = SearchArray(cpintssd, rcpint, Npzd);

            CheckValue(FUNCTION_NAME, "lo", "(rcpint not found)", (double)lo, 0,
                       (double)Npzd - 2.0);

            /* Calculate pz */
            if (cpassd[lo] == 0.0) {
              /* The slope is zero (note: cpss[lo] is always greater than zero)*/
              pz = pzd[lo] + (rcpint - cpintssd[lo]) / cpssd[lo];
            }
            else {
              /* Use the integral of the linearly interpolated pdf */

              /* Set the sign of the interpolation coefficient */
              if (cpassd[lo] < 0.0)
                cpi_sgn = -1;
              else
                cpi_sgn = 1;

              tmp = cpssd[lo]/cpassd[lo];
              pz = pzd[lo] + (cpi_sgn*sqrt((cpssd[lo]*tmp + 2.0*(rcpint
                              - cpintssd[lo]))/cpassd[lo]) - tmp);
            }

            CheckValue(FUNCTION_NAME, "pz", "", pz, pzd[lo], pzd[lo+1]);
          }

          /* Set the sign of pz */
          if (pznegative)
            pz = -pz;

          /*******************************************************************/

          /***** Solve the energy ********************************************/

          pz *= FS_CONST;
          pz2 = pz*pz;
          Ediscr = a2 - pz2*a1;
          pz2nu2 = pz2 - nu*nu;

          if ((Ediscr < 0.0) || (pz2nu2 == 0.0)) {
            /* Sample new pz to avoid negative discriminant or division by zero.
             * Ediscr can be negative due to floating point accuracy. These are
             * rare events. */
            Warn(FUNCTION_NAME, "Energy sampling failed (Ediscr: %E, pz2nu2: "
                                "%E)\n Sampling new energy...", Ediscr, pz2nu2);
            continue;
          }

          /* Photon energy */
          Ek = Ek0*(pz2*mu - nu - pz*sqrt(Ediscr))/pz2nu2;

          /* Check energy */
          if (Ek > Emax) {
            if ((Ek - Emax)/Emax < warntol) {
              /* Photon energy can be above the maximum due to floating point
               * accuracy. */
              continue;
            }
            else {
              Warn(FUNCTION_NAME, "Sampled photon energy %.15E above the "
                                  "maximum %.15E (pzmax:%E pz:%E negative:%ld)."
                                  "\nSampling new energy...\n",
                                   Ek, Emax, pzmax, pz/FS_CONST, pznegative);
              continue;
            }
          }
          else if (Ek < 0.0) {
            Warn(FUNCTION_NAME, "Sampled photon energy %.15E below zero.\n"
                                "Sampling new energy... ", Ek);
            continue;
          }

#ifdef DEBUG
          if (isnan(Ek))
            Die(FUNCTION_NAME, "NaN photon energy in Doppler broadening.");
#endif

          /* Exit the pz loop */
          break;

          /*******************************************************************/
        }

        /*********************************************************************/

        /* Rejection sampling (approximation to the Ribberfors DDCS) */
        if (RandF(id)*Ek0 <= Ek)
          break;
      }

      /***********************************************************************/

      }

    }

    /*************************************************************************/
  }

  /***************************************************************************/

  /***** Set photon energy and direction *************************************/

  /* Check energy */
  CheckValue(FUNCTION_NAME, "E", "", Ek, 0.0, Ek0 - Ebss);

  /* Set the photon energy */
  *E = Ek;

  /* Store direction cosines of incident photon */
  u0 = *u;
  v0 = *v;
  w0 = *w;

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u", "", u0, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v", "", v0, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w", "", w0, -1.01, 1.01);

  /* Rotate direction cosines of the photon */
  AziRot(mu, u, v, w, id);

  /***************************************************************************/


  /***** Calculate electron energy and direction *****************************/

  /* Set the electron kinetic energy */
  Ee = Ek0 + Ee0 - Ek - Ebss - E_RESTMASS;

  CheckValue(FUNCTION_NAME, "Ee", "", Ee, 0.0, Ek0 - Ebss);

  if ((Ee < RDB[mat + MATERIAL_PHOTON_ECUT]) ||
      ((long)RDB[DATA_PHOTON_USE_TTB] == NO)) {
    /* Electron energy is deposited locally */
    EdepTTB = Ee;
  }
  else {

    /* Sample the number of TTB photons */
    TTByield(mat, Ee, 0, &lEe, MT_ELECTRON_COMPTON, &nbr, &idx, id);

    if (nbr > 0) {

      /* Magnitude and direction cosines of momentum transfer vector */
      q = sqrt(Ek0*Ek0 + Ek*Ek - 2.0*Ek0*Ek*mu);

      if (q < ZERO)
        Die(FUNCTION_NAME, "Zero momentum transfer vector magnitude");

      uq = (Ek0*u0 - Ek*(*u))/q;
      vq = (Ek0*v0 - Ek*(*v))/q;
      wq = (Ek0*w0 - Ek*(*w))/q;


      /* Electron direction */

      if ((long)RDB[DATA_PHOTON_COMP_EANG] &&
          (long)RDB[DATA_PHOTON_USE_DOPPLER]) {

        /* Use the "pe sampling method" described in Toni's thesis */

        /* TODO:
         * - faster sampling of the pre-collision electron direction
         * - rotation matrix would probably be faster than PartialPivotGauss */

        CheckValue(FUNCTION_NAME, "momentum transfer vector direction cosines",
                   "", uq*uq + vq*vq + wq*wq - 1.0, -1E-4, 1E-4);

        pz *= E_RESTMASS;

        /* Sample the diection of the pre-collision electron */
        do {
          IsotropicDirection(&ue0, &ve0, &we0, id);
          cdelta = ue0*uq + ve0*vq + we0*wq;
          pe0c = -pz/cdelta;
        } while (pe0c < 0);


        mue0 = ue0*u0 + ve0*v0 + we0*w0;
        pec = sqrt(q*q + pe0c*pe0c - 2*pz*q);
        mue = (Ek0 - Ek*mu + pe0c*mue0)/pec;
        calpha = (q*cdelta + pe0c)/pec;
        cbeta = (Ek0*mue + pe0c*calpha - pec)/Ek;

        if (mue < -1 || mue > 1)
          Die(FUNCTION_NAME, "mue out of bounds");

        if (calpha < -1 || calpha > 1)
          Die(FUNCTION_NAME, "calpha out of bounds");

        if (cbeta < -1 || cbeta > 1)
          Die(FUNCTION_NAME, "cbeta out of bounds");

        /* Set the matrix */
        Auvw[0][0] = ue0;
        Auvw[0][1] = ve0;
        Auvw[0][2] = we0;
        Auvw[1][0] = u0;
        Auvw[1][1] = v0;
        Auvw[1][2] = w0;
        Auvw[2][0] = *u;
        Auvw[2][1] = *v;
        Auvw[2][2] = *w;

        buvw[0] = calpha;
        buvw[1] = mue;
        buvw[2] = cbeta;

        /* Solve electron direction cosines */
        PartialPivotGauss(Auvw, buvw, uvw);
        ue = uvw[0];
        ve = uvw[1];
        we = uvw[2];

        /* TODO: tarkista kulmat ym */
      }
      else {
        /* Approximation: electron travels in the direction of the momentum
         * transfer vector. This is equal to the free-electron scattering
         * angle when Doppler broadening is not used. */
        ue = uq;
        ve = vq;
        we = wq;
      }

      CheckValue(FUNCTION_NAME, "electron direction cosines", "",
                 ue*ue + ve*ve + we*we - 1.0, -1E-4, 1E-4);

      /* TTB-approximation for the Compton electron */
      TTBenergy(mat, part, Ee, x, y, z, ue, ve, we, wgt, t,
                MT_ELECTRON_COMPTON, idx, nbr, lEe, &EdepTTB, id);
    }
    else {
      /* No bremsstrahlung photons were created, the whole electron energy is
       * deposited locally */
      EdepTTB = Ee;
    }

  }

  CheckValue(FUNCTION_NAME, "EdTTB", "", EdepTTB, -INFTY, Ee);

  /***************************************************************************/


  /***** Atomic relaxation ***************************************************/

  if (Ebss < RDB[mat + MATERIAL_PHOTON_ECUT]) {
    /* Binding energy below cutoff, energy is deposited locally */
    EdepAR = Ebss;
  }
  else {

    /* Pointer to atomic relaxation data */
    ptr = (long)RDB[(long)RDB[rea + REACTION_PTR_NUCLIDE] + NUCLIDE_PTR_RELAX];
    CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, ptr);

    /* Number of subshells in the atomic relaxation data */
    Nssar = (long)RDB[ptr + RELAX_NSS];
    CheckValue(FUNCTION_NAME, "Nssd", "", (double)Nssar, 1, PHOTON_NSS_MAX);

    /* Pointer to shell eletron number CDF of the relaxation model */
    ptrelncdfar = (long)RDB[ptr + RELAX_ELNCDF];
    CheckPointer(FUNCTION_NAME, "(ELNCDF)", DATA_ARRAY, ptrelncdfar);

    /* Check that the relaxation and Compton electron numbers agree */
    if (Nssar <= ss) {
      ss = -1;
    }
    else if ((long)elncdf[ss + 1] != (long)RDB[ptrelncdfar + ss + 1]) {

      /* Sample the electron number starting from 1 (note: elncdf[0] = 0) */
      eli = (long)(elncdf[ss] + ceil(RandF(id)*elnd[ss]));

      /* Find the subshell in the relaxation data which corresponds to the
       * electron number */
      ss = -1;
      for (i = 0; i < Nssar; i++) {
        if ((long)RDB[ptrelncdfar + i + 1] >= eli) {
          ss = i;
          break;
        }
      }
    }

    if (ss == -1) {
      /* Relaxtion data not available for the subshell, energy is deposited
       * locally */
      /* NOTE: This can also happen when different binding energies are used
       * here and in AtomicRelaxation(). */
      EdepAR = Ebss;
    }
    else {
      /* Atomic relaxation */
      AtomicRelaxation(mat, rea, part, ss, x, y, z, wgt, t, &EdepAR, id);

      /* Due to the different binding energies used by AtomicRelaxation(), the
       * deposited energy is corrected so that the total energy is conserved.
       * EdepAR can be negative, but it is a rare event. */
      ebidar = &RDB[(long)RDB[ptr + RELAX_EBI]];
      CheckPointer(FUNCTION_NAME, "(RELAX_EBI)", DATA_ARRAY,
                   (long)RDB[ptr + RELAX_EBI]);
      EdepAR = EdepAR - ebidar[ss] + Ebss;
    }
  }

  /***************************************************************************/

  /* Total deposited energy */
  Edep = EdepTTB + EdepAR;

  /* Score pulse-height detector */

  PulseDet(part, mat, Edep, x, y, z, wgt, id);

  /* Score photon heat detector */

  ScorePhotonHeat(part, mat, Edep, x, y, z, Ek0, t, wgt, id);

}

/*****************************************************************************/


/*****************************************************************************/
void PartialPivotGauss(double A[3][3], double b[3], double x[3]) {
  /* Solves linear equation Ax=b with Gaussian elimination using scaled
   * partial pivoting.
   * NOTE: The elements of A and b are changed! */
  long i, j, k, tmp;
  long idx[3];
  double smax, aij, r, rmax, sum;
  double s[3];
  long n = 3;

  /* Calculate a scale factor array */
  for (i = 0; i < n; i++) {
    smax = 0.0;
    for (j = 0; j < n; j++) {
      aij = fabs(A[i][j]);
      if (aij > smax)
        smax = aij;
    }
    s[i] = smax;

    /* Initialize the index array */
    idx[i] = i;
  }

  /* Forward elimination */
  for (i = 0; i < n - 1; i++) {

    rmax = 0;
    k = i;

    /* Find pivot equation */
    for (j = i; j < n; j++) {
      r = fabs(A[idx[j]][i])/s[idx[j]];
      if (r > rmax) {
        rmax = r;
        k = j;
      }
    }

    /* Update the index array */
    tmp = idx[k];
    idx[k] = idx[i];
    idx[i] = tmp;

    /* Elimination */
    for (j = i+1; j < n; j++) {
      r = A[idx[j]][i]/A[idx[i]][i];
      A[idx[j]][i] = r;
      for (k = i+1; k < n; k++)
        A[idx[j]][k] = A[idx[j]][k] - r*A[idx[i]][k];

      /* Calculate b */
      b[idx[j]] = b[idx[j]] - A[idx[j]][i]*b[idx[i]];
    }
  }

  /* Back substitution */
  x[n-1] = b[idx[n-1]]/A[idx[n-1]][n-1];
  for (i = n-2; i >= 0; i--) {
    sum = b[idx[i]];
    for (j = i + 1; j < n; j++)
      sum = sum - A[idx[i]][j]*x[j];
    x[i] = sum/A[idx[i]][i];
  }

}
/*****************************************************************************/
