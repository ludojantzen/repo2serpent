/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoresurf.c                                    */
/*                                                                           */
/* Created:       2014/08/14 (JLe)                                           */
/* Last modified: 2014/08/15 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Handles all surface based tallying.                          */
/*                                                                           */
/* Comments: - The track used to check surface crossings refers to a track   */
/*             between two points where the direction of the neutron may     */
/*             have changed, i.e. physical collisions and outer boundary     */
/*             crossings.                                                    */
/*                                                                           */
/*           - Toi s:n arvo voi mennä numeriikan takia nollaksi jos paikka-  */
/*             koordinaatit on tosi isoja ja radan pituus pieni --> FAIL.    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreSurf:"

/*****************************************************************************/

void ScoreSurf(long part, double *x0, double *y0, double *z0, double x, 
               double y, double z, double u, double v, double w, double E, 
               double wgt, double t, long id)
{
  double s;

  /* Check if tracks are stopped at outer boundary */
                  
  if ((long)RDB[DATA_STOP_AT_BOUNDARY] == NO)
    return;

  /* Calculate length of track */

  s = sqrt((x - *x0)*(x - *x0) + (y - *y0)*(y - *y0) + (z - *z0)*(z - *z0));
  CheckValue(FUNCTION_NAME, "s", "", s, 0.0, INFTY);

  /* Extrapolate (tää tarvitaan geometrian kanssa yhtenevien */
  /* pintojen takia). */

  s = s + EXTRAP_L;

  /* Score super-imposed detectors */

  SuperDet(part, *x0, *y0, *z0, u, v, w, s, E, t, wgt, id);

  /* For EDo (18.1.2014) */
  
  DiffCoefED(2, *x0, *y0, *z0, u, v, w, s, E, wgt, id);
  
  /* Score discontinuity factors */
  
  ScoreDF(*x0, *y0, *z0, u, v, w, s, E, wgt, id);

  /* Score albedo surface currents */

  ScoreAlb(part, *x0, *y0, *z0, u, v, w, s, E, wgt, id);

  /* Score ICM surface currents */

  ScoreICMTrk(part, *x0, *y0, *z0, u, v, w, s, E, wgt, id);

  /* Set new initial coordinates for track */

  *x0 = x;
  *y0 = y;
  *z0 = z;
}

/*****************************************************************************/
