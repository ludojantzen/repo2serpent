/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : inelasticscattering.c                          */
/*                                                                           */
/* Created:       2011/02/28 (JLe)                                           */
/* Last modified: 2017/06/06 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Handles inelastic scattering reactions for neutrons          */
/*                                                                           */
/* Comments: - From Serpent 1.1.14 (2.9.2010)                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InelasticScattering:"

/*****************************************************************************/

void InelasticScattering(long rea, double *E, double *u, double *v, double *w,
                         long id)
{
  double mu, E0, ac, Ecm;
  long erg, law;

  /***************************************************************************/

  /***** Get initial values and check ****************************************/

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Check initial energy */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);

  /* Get pointer to energy distribution */

  erg = (long)RDB[rea + REACTION_PTR_ERG];

  /* NOTE: ProcessEDistributions():issa on alussa mt-tarkistus, joka   */
  /* poistaa myöhemmmin lisättyjä reaktioita (7.2.2016 / 2.1.25 / JLe) */

  if (erg < VALID_PTR)
    Die(FUNCTION_NAME, "null pointer mt = %ld\n",
        (long)RDB[rea + REACTION_MT]);

  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Get energy distribution type */

  law = (long)RDB[erg + ERG_LAW];

  /***************************************************************************/

  /***** Remember some values before the collision ***************************/

  /* Initial energy */

  E0 = *E;

  /***************************************************************************/

  /***** Sample energy and scattering angle **********************************/

  if ((law == 44) || (law == 61) || (law == 66))
    {
      /* Combined energy-angle -distributions */

      SampleENDFLaw(rea, -1, E0, E, &mu, id);
    }
  else
    {
      /* Sample energy and angle from separate distributions */

      SampleENDFLaw(rea, -1, E0, E, NULL, id);
      mu = SampleMu(rea, -1, E0, NULL, NULL, id);
    }

  /* Check that scattering has occured */

  if (*E == E0)
    return;

  /* Check values (Sampled energy can be higher for some scattering laws) */

  CheckValue(FUNCTION_NAME, "mu", "(2)", mu, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "E", "(2)", *E, 0.0, INFTY);

  /****************************************************************************/

  /***** Make final adjustments ***********************************************/

  /* Convert values to L-frame if sampled in C-frame. */

  if ((long)RDB[rea + REACTION_TY] < 0)
    {
      /* Get target awr constant */

      ac = RDB[rea + REACTION_AWR] + 1.0;

      /* Final energy in C-frame */

      Ecm = *E;

      /* Get target awr constant */

      ac = RDB[rea + REACTION_AWR] + 1.0;

      /* Energy from C- to L-frame */

      *E = Ecm + (E0 + 2.0*mu*ac*sqrt(E0*Ecm))/(ac*ac);

      /* Cosine from C- to L-frame */

      mu = mu*sqrt(Ecm/(*E)) + sqrt(E0/(*E))/ac;

      /* Adjust values if > +- 1.0 */

      if (mu > 1.0)
        mu = 1.0;
      else if (mu < -1.0)
        mu = -1.0;
    }

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u", "", *u, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v", "", *v, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w", "", *w, -1.01, 1.01);

  /* Rotate direction cosines around a random azimuthal angle */

  AziRot(mu, u, v, w, id);
}

/*****************************************************************************/
