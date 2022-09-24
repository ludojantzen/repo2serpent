/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : clearbuf.c                                     */
/*                                                                           */
/* Created:       2010/11/10 (JLe)                                           */
/* Last modified: 2018/09/13 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Clears scoring buffer(s) in all bins                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ClearBuf:"

/*****************************************************************************/

void ClearBuf()
{
  long sz, max, i, nseg;

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "BUF array not ready for access");

  /* Get buffer segment and data size */

  sz = (long)RDB[DATA_REAL_BUF_SIZE];
  max = (long)RDB[DATA_ALLOC_BUF_SIZE];

  /* Number of segments */

  if ((long)RDB[DATA_OPTI_SHARED_BUF] == YES)
    nseg = 1;
  else
    nseg = (long)RDB[DATA_OMP_MAX_THREADS];

  /* Start parallel timer */

  StartTimer(TIMER_OMP_PARA);

  /* Loop over segments and reset data */

#ifdef OPEN_MP
#pragma omp parallel private(i)
#endif
  {

#ifdef OPEN_MP
#pragma omp for
#endif

    for (i = 0; i < nseg; i++)
      memset(&BUF[i*sz], 0.0, max*sizeof(double));
  }

  /* Stop parallel timer */

  StopTimer(TIMER_OMP_PARA);

  /* Reset reduced flag */

  WDB[DATA_BUF_REDUCED] = (double)NO;

  /* Reset reduced flag for RES2 array (miksi tää tehdään?) */

  WDB[DATA_RES2_REDUCED] = (double)NO;
}

/*****************************************************************************/
