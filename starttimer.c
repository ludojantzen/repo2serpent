/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : starttimer.c                                   */
/*                                                                           */
/* Created:       2010/11/14 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Starts a wall-clock timer                                    */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "StartTimer:"

/*****************************************************************************/

void StartTimer(long i)
{
  struct timeb tb;
  double t;

  /* Check index */

  if ((i < 0) || (i > TOT_TIMERS))
    Die(FUNCTION_NAME, "Invalid timer index %d", i);

  /* Check Open MP thread number */

  if (OMP_THREAD_NUM != 0)
    return;

  /* Set timer on */

  if (timer[i].on == NO)
    timer[i].on = YES;
  else
    Die(FUNCTION_NAME, "timer %d is already on", i);

  /* Get current wall-clock time */
  
  ftime(&tb);
  t = (double)tb.time + (double)tb.millitm/1000.0;
  
  /* Set value */

  timer[i].t0 = t;

  /* Get current CPU time */

  t = (double)clock();
  
  /* Set value */

  timer[i].cpu_t0 = t;
}

/*****************************************************************************/
