/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stoptimer.c                                    */
/*                                                                           */
/* Created:       2010/11/14 (JLe)                                           */
/* Last modified: 2016/10/04 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Stops a wall-clock timer                                     */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "StopTimer:"

/*****************************************************************************/

void StopTimer(long i)
{
  struct timeb tb;
  double t;

  /* Check index */

  if ((i < 0) || (i > TOT_TIMERS))
    Die(FUNCTION_NAME, "Invalid timer index %d", i);

  /* Check Open MP thread number */

  if (OMP_THREAD_NUM != 0)
    return;

  /* Set timer of */

  if (timer[i].on == YES)
    timer[i].on = NO;
  else
    Die(FUNCTION_NAME, "timer %d is already off", i);

  /* Get current wall-clock time */
  
  ftime(&tb);
  t = (double)tb.time + (double)tb.millitm/1000.0;
  
  /* Add time to counter */

  timer[i].t = timer[i].t + t - timer[i].t0;

  /* Get current CPU time */
  
  t = (double)clock();
  
  /* Add time to counter */

  timer[i].cpu_t = timer[i].cpu_t + t - timer[i].cpu_t0;
}

/*****************************************************************************/
