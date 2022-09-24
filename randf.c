/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : randf.c                                        */
/*                                                                           */
/* Created:       2010/11/22 (JLe)                                           */
/* Last modified: 2015/04/04 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Samples a uniformly distribute random number on the unit     */
/*              interval.                                                    */
/*                                                                           */
/* Comments: - Calls Rand64(), which is based on an algorithm that produces  */
/*             the same random number sequence for a particle history in     */
/*             both serial and parallel modes.                               */
/*                                                                           */
/*           - Should be used only during the transport cycle, for other     */
/*             purposes use C-function drand48().                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RandF:"

/*****************************************************************************/

double RandF(long id)
{
  unsigned long seed;
  double f;

  /* Get seed */

  seed = SEED[id*RNG_SZ];
  CheckValue(FUNCTION_NAME, "seed", "", seed, 1, INFTY);

  /* Sample rng */
  
  seed *= 2862933555777941757;
  seed += 12345;

  /* Conversion to floating point number in interval [0,1) */
  
  f = (double)(seed >> 12);
  f /= 0x0010000000000000;
  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);

  /* Store seed */

  SEED[id*RNG_SZ] = seed;

  /* Return value */

  return f;
}

/*****************************************************************************/
