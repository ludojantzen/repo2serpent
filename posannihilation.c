/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : posannihilation.c                              */
/*                                                                           */
/* Created:       2016/09/26 (TKa)                                           */
/* Last modified: 2017/03/10 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Positron annihilation in flight and at rest                  */
/*                                                                           */
/* Comments: In flight annihilation cannot be used with the current TTB      */
/*           method.                                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"


/* Local function definitions */
void PosAnnihilationInflight(double Ep, double up, double vp, double wp,
                             double *Eg1, double uvw1[3], double *Eg2,
                             double uvw2[3], long id);
void PosAnnihilationAtrest(double *Eg1, double uvw1[3], double *Eg2,
                           double uvw2[3], long id);

/*****************************************************************************/

void PosAnnihilation(long mat, long part, double Ep, double x, double y,
                     double z, double up, double vp, double wp, double wgt,
                     double t, long id) {
  /* NOTE: part is a photon pointer!  */
  static char * const FUNCTION_NAME = "PosAnnihilation:";
  long new1, new2, ptr;
  double Eg1, Eg2;
  double uvw1[3];
  double uvw2[3];
  double EcutAn = INFTY;

  /* Choose annihilation model */

  if (Ep > EcutAn) {
    /* In flight */
    Die(FUNCTION_NAME, "Positron in-flight annihilation not supported");
    /* TODO for electron transport: Energy deposition */
    PosAnnihilationInflight(Ep, up, vp, wp, &Eg1, uvw1, &Eg2, uvw2, id);

    /* TODO / NOTE: Annihilaatio-fotoneille pitää tehdä energiakatkaisu tässä
     * funktiossa, ja se pitää ottaa pois PairProduction():sta. */
  }
  else {
    /* At rest */
    /* TODO for electron transport: Deposit the residual
     * energy (Ed=Ep-Eg1-Eg2)? */
    PosAnnihilationAtrest(&Eg1, uvw1, &Eg2, uvw2, id);

    if (Eg1 != Eg2)
      Die(FUNCTION_NAME, "Different annihilation photon energies");
  }

  /* Sanity checks */
  CheckValue(FUNCTION_NAME, "u1", "", uvw1[0], -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v1", "", uvw1[1], -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w1", "", uvw1[2], -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u2", "", uvw2[0], -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v2", "", uvw2[1], -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w2", "", uvw2[2], -1.01, 1.01);


  /* Make two new photons */
  new1 = DuplicateParticle(part, id);
  new2 = DuplicateParticle(part, id);

  /* Put variables */
  WDB[new1 + PARTICLE_X] = x;
  WDB[new1 + PARTICLE_Y] = y;
  WDB[new1 + PARTICLE_Z] = z;
  WDB[new1 + PARTICLE_WGT] = wgt;
  WDB[new1 + PARTICLE_T] = t;
  WDB[new1 + PARTICLE_E] = Eg1;
  WDB[new1 + PARTICLE_U] = uvw1[0];
  WDB[new1 + PARTICLE_V] = uvw1[1];
  WDB[new1 + PARTICLE_W] = uvw1[2];
  WDB[new1 + PARTICLE_PTR_MAT] = (double)mat;

  /* Set photon emission type */

  WDB[new1 + PARTICLE_PHOTON_TYPE] = PHOTON_TYPE_ANNIH;

  /* Put variables */
  WDB[new2 + PARTICLE_X] = x;
  WDB[new2 + PARTICLE_Y] = y;
  WDB[new2 + PARTICLE_Z] = z;
  WDB[new2 + PARTICLE_WGT] = wgt;
  WDB[new2 + PARTICLE_T] = t;
  WDB[new2 + PARTICLE_E] = Eg2;
  WDB[new2 + PARTICLE_U] = uvw2[0];
  WDB[new2 + PARTICLE_V] = uvw2[1];
  WDB[new2 + PARTICLE_W] = uvw2[2];
  WDB[new2 + PARTICLE_PTR_MAT] = (double)mat;

  /* Set photon emission type */

  WDB[new2 + PARTICLE_PHOTON_TYPE] = PHOTON_TYPE_ANNIH;

  /* Score particle balance */

  ptr = (long)RDB[RES_G_BALA_SRC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf(2.0, 1.0, ptr, id, -1, BALA_G_SRC_ANNIH, 0);
  AddBuf(2.0*wgt, 1.0, ptr, id, -1, BALA_G_SRC_ANNIH, 1);
  AddBuf(wgt*(Eg1+Eg2), 1.0, ptr, id, -1, BALA_G_SRC_ANNIH, 2);

  /* Put photons in queue */
  ToQue(new1, id);
  ToQue(new2, id);
}

/*****************************************************************************/

void PosAnnihilationInflight(double Ep, double up, double vp, double wp,
                             double *Eg1, double uvw1[3], double *Eg2,
                             double uvw2[3], long id) {
  /* Positron annihilation in flight.
   *
   * Samples the two-photon annihilation differential cross section given
   * by Heitler. Two rejection sampling methods are used, method 1 when the
   * positron energy is below Elim and method 2 above it. It was tested
   * (Serpent 2.1.26 with the default compiler options) that the method 2
   * becomes faster at about 8 MeV. At this energy, the efficiency of the
   * method 1 is only about 0.23, whereas the method 2 has the efficiency of
   * 0.97. However, the method 2 uses log() and exp() functions, which
   * explains the speed difference.
   * */
  static char * const FUNCTION_NAME = "PosAnnihilationInflight:";
  static const double Elim = 8.0;
  double gamma, gamma2, zeta, zeta_min, zeta_max, kappa, sg21, g13, Cninv,
         tmp, mu1, mu2, a, b;

  /* Set variables */
  gamma = 1.0 + Ep/E_RESTMASS;
  gamma2 = gamma*gamma;
  sg21 = sqrt(gamma2 - 1.0);
  zeta_min = 1.0/(gamma + 1.0 + sg21);
  zeta_max = 0.5;
  kappa = (gamma2 + 4.0*gamma + 1.0);

  /* Sample zeta = Eg1/(Ep+2*E_RESTMASS) */

  if (Ep < Elim) {
    /* Sampling method 1: */
    /* Rejection sampling using a uniform rejection function */

    g13 = (gamma + 1.0);
    g13 *= 2.0*g13*g13;

    do {
      zeta = zeta_min + (zeta_max - zeta_min)*RandF(id);
      tmp = 1.0/(zeta*(1.0 - zeta));
    } while (RandF(id)*g13 > tmp*(kappa + 2.0 - tmp));
  }
  else {
    /* Sampling method 2: */
    /* Rejection sampling using a rejection function
     *  f(zeta) = Cn/(zeta*(1-zeta)) */

    /* Inverse of the normalization constant Cn */
    Cninv = log(1.0/zeta_min - 1.0);

    do {
      tmp = exp(RandF(id)*Cninv);
      zeta = 1.0/(1.0 + tmp);
    } while (RandF(id)*(kappa - 2.0) > kappa + 2.0 - (1.0 + tmp)/(1.0 - zeta));
  }

  /* Check zeta */
  if (zeta < zeta_min)
    Die(FUNCTION_NAME, "Variable zeta = %f below %f", zeta_min);
  if (zeta > zeta_max)
    Die(FUNCTION_NAME, "Variable zeta = %f above %f", zeta_max);

  /* Set photon energies (Eg1 < Eg2) */
  *Eg1 = zeta*(Ep + 2.0*E_RESTMASS);
  *Eg2 = (1.0 - zeta)*(Ep + 2.0*E_RESTMASS);

  /* Set the polar angle cosines */
  mu1 = (gamma + 1.0 - 1.0/zeta)/sg21;
  mu2 = (gamma + 1.0 - 1.0/(1.0 - zeta))/sg21;

  /* Initialize the direction cosines of the photon 1 */
  uvw1[0] = up;
  uvw1[1] = vp;
  uvw1[2] = wp;

  /* Sanity checks */
  CheckValue(FUNCTION_NAME, "mu1", "", mu1, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "mu2", "", mu2, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "up", "", up, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "vp", "", vp, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "wp", "", wp, -1.01, 1.01);

  /* Rotate */
  AziRot(mu1, uvw1, uvw1+1, uvw1+2, id);

  /* Calculate the direction cosines of the photon 2 */
  /* NOTE: The variable 'a' can also be written as
   *  a = sqrt((1.0 - mu2*mu2)/(1.0 - mu1*mu1))
   * but it reduces neatly to zeta/(1-zeta). */
  a = zeta/(1.0 - zeta);
  b = mu2 + a*mu1;
  uvw2[0] = b*up - a*uvw1[0];
  uvw2[1] = b*vp - a*uvw1[1];
  uvw2[2] = b*wp - a*uvw1[2];

  CheckValue(FUNCTION_NAME, "photon number 2 direction cosines", "",
             uvw2[0]*uvw2[0] + uvw2[1]*uvw2[1] + uvw2[2]*uvw2[2] - 1.0,
             -1E-4, 1E-4);

}
/*****************************************************************************/


/*****************************************************************************/

void PosAnnihilationAtrest(double *Eg1, double uvw1[3], double *Eg2,
                           double uvw2[3], long id) {
  /* Positron annihilation at rest.
   * */

  /* Put energies */
  *Eg1 = E_RESTMASS;
  *Eg2 = E_RESTMASS;

  /* Sample the direction of the first photon isotropically */
  IsotropicDirection(uvw1, uvw1+1, uvw1+2, id);

  /* The second photon goes to the opposite direction */
  uvw2[0] = -uvw1[0];
  uvw2[1] = -uvw1[1];
  uvw2[2] = -uvw1[2];
}

/*****************************************************************************/
