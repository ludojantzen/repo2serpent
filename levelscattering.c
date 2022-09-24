/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : levelscattering.c                              */
/*                                                                           */
/* Created:       2011/02/28 (JLe)                                           */
/* Last modified: 2020/05/14 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Handles inelastic level scattering for neutrons              */
/*                                                                           */
/* Comments: - From Serpent 1.1.13 (10.8.2010)                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "LevelScattering:"

/*****************************************************************************/

void LevelScattering(long rea, double *E, double *u, double *v, double *w,
                     long id)
{
  double awr, muc, E0, Q, C;
  double V, Vx, Vy, Vz, Vt, Vtx, Vty, Vtz, Px, Py, Pz;
  double Vcx, Vcy, Vcz, u0, v0, w0;
  long nuc;

  /* Check that scattering is in C-frame */

  if ((long)RDB[rea + REACTION_TY] != -1)
    {
      /* L-frame scattering, use general inelastic routine */

      InelasticScattering(rea, E, u, v, w, id);

      /* Exit subroutine */

      return;
    }

  /***************************************************************************/

  /***** Get initial values and check ****************************************/

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to nuclide */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Check initial energy */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);

  /* Get atomic weight ratio */

  awr = RDB[nuc + NUCLIDE_AWR];
  CheckValue(FUNCTION_NAME, "awr", "", awr, 0.9, 300.0);

  /* Get reaction Q-value */

  Q = RDB[rea + REACTION_Q];
  CheckValue(FUNCTION_NAME, "Q", "", Q, -20.0, 10.0);

  /* Initial neutron velocity */

  V = sqrt(*E);

  Vx = *u*V;
  Vy = *v*V;
  Vz = *w*V;

  /***************************************************************************/

  /***** Remember some values before the collision ***************************/

  /* Total energy */

  E0 = *E;

  /* Direction cosines */

  u0 = *u;
  v0 = *v;
  w0 = *w;

  /* Total momentum */

  Px = Vx;
  Py = Vy;
  Pz = Vz;

  /* Check */

  if (Px*Px + Py*Py + Pz*Pz < ZERO)
    Die(FUNCTION_NAME, "Error in momentum");

  /***************************************************************************/

  /***** Transformation to C-frame *******************************************/

  /* Calculate velocity of centre-of-mass */

  Vcx = Vx/(awr + 1.0);
  Vcy = Vy/(awr + 1.0);
  Vcz = Vz/(awr + 1.0);

  /* Neutron velocities in C-frame */

  Vx = Vx - Vcx;
  Vy = Vy - Vcy;
  Vz = Vz - Vcz;

  V = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

  /* Target velocities in C-frame */

  Vtx = -Vcx;
  Vty = -Vcy;
  Vtz = -Vcz;

  Vt = sqrt(Vtx*Vtx + Vty*Vty + Vtz*Vtz);

  /***************************************************************************/

  /***** Scattering in C-frame ***********************************************/

  /* Sample scattering cosine in C-frame. NOTE: Varmista että toi energia */
  /* tosiaan on L-framessa. Muutos C-frameen kertoimella awr/(1.0 + awr). */

  muc = SampleMu(rea, -1, E0, NULL, NULL, id);

  /* Calculate direction cosines */

  *u = Vx/V;
  *v = Vy/V;
  *w = Vz/V;

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "muc", "", muc, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u", "", *u, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v", "", *v, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w", "", *w, -1.01, 1.01);

  /* Rotate */

  AziRot(muc, u, v, w, id);

  /* Calculate change in neutron speed. */

  if ((C = V*V + Q*awr/(awr + 1.0)) > 0.0)
    V = sqrt(C);
  else
    {
      /* Set original direction cosines (energy has not yet been changed) */

      *u = u0;
      *v = v0;
      *w = w0;

      /* Print warning */

#ifdef DEBUG

      Warn(FUNCTION_NAME, "Energy below threshold (%s, E = %E, Q = %E)",
           GetText(nuc + NUCLIDE_PTR_NAME), *E, Q);

#endif

      /* Ignore reaction */

      return;
    }

  /* Calculate change in target speed */

  Vt = sqrt(Vt*Vt + Q/(awr*(awr + 1.0)));

  /* Velocities after collision. */

  Vx = *u*V;
  Vy = *v*V;
  Vz = *w*V;

  Vtx = -*u*Vt;
  Vty = -*v*Vt;
  Vtz = -*w*Vt;

  /***************************************************************************/

  /***** Transformation back to L-frame **************************************/

  /* Neutron velocities in L-frame */

  Vx = Vx + Vcx;
  Vy = Vy + Vcy;
  Vz = Vz + Vcz;

  V = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

  /* Target velocities in L-frame */

  Vtx = Vtx + Vcx;
  Vty = Vty + Vcy;
  Vtz = Vtz + Vcz;

  Vt = sqrt(Vtx*Vtx + Vty*Vty + Vtz*Vtz);

  /* Set neutron energy */

  *E = (Vx*Vx + Vy*Vy + Vz*Vz);

  /* Set direction cosines */

  *u = Vx/V;
  *v = Vy/V;
  *w = Vz/V;

  /***************************************************************************/

  /***** Check conservation of energy and momentum ***************************/

#ifdef DEBUG

  /* Check relative change in total energy */

  E0 = 1.0 - (-Q + *E + awr*Vt*Vt)/E0;
  CheckValue(FUNCTION_NAME, "E", "", E0, -1E-2, 1E-2);

  /* Check relative change in momentum components */

  if (Px != 0.0)
    Px = 1.0 - (Vx + awr*Vtx)/Px;

  if (Py != 0.0)
    Py = 1.0 - (Vy + awr*Vty)/Py;

  if (Pz != 0.0)
    Pz = 1.0 - (Vz + awr*Vtz)/Pz;

  CheckValue(FUNCTION_NAME, "Px", "", Px, -1E-4, 1E-4);
  CheckValue(FUNCTION_NAME, "Py", "", Py, -1E-4, 1E-4);
  CheckValue(FUNCTION_NAME, "Pz", "", Pz, -1E-4, 1E-4);

#endif

  /***************************************************************************/
}

/*****************************************************************************/
