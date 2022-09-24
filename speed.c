/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : speed.c                                        */
/*                                                                           */
/* Created:       2012/10/11 (JLe)                                           */
/* Last modified: 2019/03/29 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates speed from energy                                 */
/*                                                                           */
/* Comments: - Speed can be set by user e.g. for track plot animations.      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Speed:"

/*****************************************************************************/

double Speed(long type, double E)
{
  double spd;
  
  /* Avoid compiler warning */

  spd = -1.0;

  /* Check type */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      /* Check override */

      if ((spd = RDB[DATA_NEUTRON_SPD]) < 0.0)
        {
          /* Calculate relativistic speed */
          
          spd = E/NEUTRON_E0 + 1.0;
          spd = SPD_C*sqrt(1.0 - 1.0/(spd*spd));
        }
    }
  else if (type == PARTICLE_TYPE_GAMMA)
    {
      /* Check override and set photon speed */

      if ((spd = RDB[DATA_PHOTON_SPD]) < 0.0)
        spd = SPD_C;
    }
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Return speed */

  return spd;
}

/*****************************************************************************/
