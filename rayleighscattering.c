/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : rayleighscattering.c                           */
/*                                                                           */
/* Created:       2014/06/18 (TKa)                                           */
/* Last modified: 2017/02/09 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Handles coherent Rayleigh scattering for photons             */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RayleighScattering:"

/*****************************************************************************/

void RayleighScattering(long rea, double E, double *u, double *v, double *w,
                        long id)
{
  long ptr, ptd, Nd, idx0, idx1;
  double x2, x2max, lx2max, ff2intmax, mu, rff2intd;
  static const double zt = 6.5053004697083e+03; /* (1e6/(h*c)*angstrom)^2 */
  const double *x2d, *lx2d, *ff2d, *ff2intd, *cd;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to photon reaction data */
  ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Size of the data arrays */
  Nd = (long)RDB[ptd + PHOTON_DIST_RAYL_N];

  /* Squared momentum transfer (= x^2) */
  ptr = (long)RDB[ptd + PHOTON_DIST_RAYL_X2];
  CheckPointer(FUNCTION_NAME, "(X2)", DATA_ARRAY, ptr);
  x2d = &RDB[ptr];

  /* Log squared momentum transfer */
  ptr = (long)RDB[ptd + PHOTON_DIST_RAYL_LX2];
  CheckPointer(FUNCTION_NAME, "(LX2)", DATA_ARRAY, ptr);
  lx2d = &RDB[ptr];

  /* Pointer to squared form factors */
  ptr = (long)RDB[ptd + PHOTON_DIST_RAYL_FF2];
  CheckPointer(FUNCTION_NAME, "(FF2D)", DATA_ARRAY, ptr);
  ff2d = &RDB[ptr];

  /* Pointer to integrated squared form factor (i.e. cumulative integral) */
  ptr = (long)RDB[ptd + PHOTON_DIST_RAYL_FF2INT];
  CheckPointer(FUNCTION_NAME, "(FF2INT)", DATA_ARRAY, ptr);
  ff2intd = &RDB[ptr];

  /* Interpolation constants */
  ptr = (long)RDB[ptd + PHOTON_DIST_RAYL_C];
  CheckPointer(FUNCTION_NAME, "(C)", DATA_ARRAY, ptr);
  cd = &RDB[ptr];
  

  /* Calculate the maximum x^2 */
  x2max = zt*E*E;
  lx2max = log(x2max);

  /***** Calculate the maximum integrated squared form factor ****************/

  if (x2max < x2d[Nd-1]) {

    /* Find the maximum from the log momentum transfer array */
    idx0 = SearchArray(lx2d, lx2max, Nd);

    if (idx0 == -1)
      Die(FUNCTION_NAME, "Maximum squared form factor smaller than the minimum "
                         "in the data");

    /* NOTE: ProcessRayleigh() checks that cd[idx0] != -1 */
    ff2intmax = ff2intd[idx0] + ff2d[idx0]*x2d[idx0]/(1.0 + cd[idx0])
                *(exp((1.0 + cd[idx0])*(lx2max - lx2d[idx0])) - 1.0);

    CheckValue(FUNCTION_NAME, "ff2intmax", "", ff2intmax, ff2intd[idx0],
               ff2intd[idx0+1]);
  }
  else {
    ff2intmax = ff2intd[Nd-1];
    idx0 = Nd - 2;
  }

  /***************************************************************************/


  /***** Sample mu using rejection sampling **********************************/

  do {
    /* Generate a random number between 0 and ff2intmax */
    rff2intd = RandF(id)*ff2intmax;

    /* Search interval, note the upper limit idx0 + 2 */
    idx1 = SearchArray(ff2intd, rff2intd, idx0 + 2);

    if (idx1 == -1)
      Die(FUNCTION_NAME, "rff2intd not found");

    CheckValue(FUNCTION_NAME, "idx1", "", (double)idx1, 0.0, (double)idx0);

    x2 = x2d[idx1]*pow((1.0 + cd[idx1])*(rff2intd - ff2intd[idx1])
                       /(ff2d[idx1]*x2d[idx1]) + 1.0, 1.0/(1.0 + cd[idx1]));

    CheckValue(FUNCTION_NAME, "x2", "", x2, x2d[idx1], x2d[idx1+1]);

    mu = 1.0 - 2.0*x2/x2max;

  } while (1.0 + mu*mu < 2.0*RandF(id));

  /***************************************************************************/

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u", "", *u, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v", "", *v, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w", "", *w, -1.01, 1.01);

  /* Rotate direction cosines */
  AziRot(mu, u, v, w, id);

}

/*****************************************************************************/
