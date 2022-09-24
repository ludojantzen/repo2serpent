/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resettimer.c                                   */
/*                                                                           */
/* Created:       2010/11/14 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Resets a wall-clock timer                                    */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ResetTimer:"

/*****************************************************************************/

void ResetTimer(long i)
{
  /* Check index */

  if ((i < 0) || (i > TOT_TIMERS))
    Die(FUNCTION_NAME, "Invalid timer index %d", i);

  /* Check Open MP thread number */

  if (OMP_THREAD_NUM != 0)
    return;

  /* Reset values */

  timer[i].t0 = 0.0;
  timer[i].t = 0.0;
  timer[i].cpu_t0 = 0.0;
  timer[i].cpu_t = 0.0;
}

/*****************************************************************************/
