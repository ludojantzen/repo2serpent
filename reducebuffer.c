/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reducebuffer.c                                 */
/*                                                                           */
/* Created:       2010/11/12 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Reduces data from OpenMP distributed buffer to thread 0      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReduceBuffer:"

/*****************************************************************************/

void ReduceBuffer()
{
  long sz, i, n, max, nseg;
  double val;

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "BUF array not ready for access");

  /* Check and set reduced flag */
  
  if ((long)RDB[DATA_BUF_REDUCED] == YES)
    return;
  else
    WDB[DATA_BUF_REDUCED] = (double)YES;

  /* Return if shared buffer */

  if ((long)RDB[DATA_OPTI_SHARED_BUF] == YES)
    return;

  /* Number of segments */
  
  if ((nseg = (long)RDB[DATA_OMP_MAX_THREADS]) < 2)
    return;
 
  /* Get buffer segment and data size */

  sz = (long)RDB[DATA_REAL_BUF_SIZE];
  max = (long)RDB[DATA_ALLOC_BUF_SIZE];
 
  /* Copy data */

#ifdef OPEN_MP
#pragma omp parallel private(val, n, i)
#endif
  {
    
#ifdef OPEN_MP
#pragma omp for
#endif          
    
    for (n = 0; n < max; n++)
      {
        val = 0.0;
        
        for (i = 1; i < nseg; i++)
          BUF[n] += BUF[i*sz + n];
        

      }
  }
  
  /* Clear other buffers */

#ifdef OPEN_MP
#pragma omp parallel private(i)
#endif
  {
    
#ifdef OPEN_MP
#pragma omp for
#endif          
    
    for (i = 1; i < nseg; i++)
      memset(&BUF[i*sz], 0.0, max*sizeof(double));
  }
  
  /****************************************************************************/
}

/*****************************************************************************/
