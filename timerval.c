/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : timerval.c                                     */
/*                                                                           */
/* Created:       2010/11/14 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns wall-clock time from timer                           */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "TimerVal:"

/*****************************************************************************/

double TimerVal(long i)
{
  struct timeb tb;
  double t;

  /* Check index */

  if ((i < 0) || (i > TOT_TIMERS))
    Die(FUNCTION_NAME, "Invalid timer index %d", i);

  /* Return total time if timer is stopped */

  if (timer[i].on == NO)
    return timer[i].t;

  /* Get current wall-clock time */
  
  ftime(&tb);
  t = (double)tb.time + (double)tb.millitm/1000.0;
  
  /* Return time */

  return timer[i].t + t - timer[i].t0;
}

/*****************************************************************************/
