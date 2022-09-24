/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : omptestcomp.c                                  */
/*                                                                           */
/* Created:       2019/02/26 (JLe)                                           */
/* Last modified: 2019/02/26 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reset thread-wise completed flags                            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "OMPTestComp:"

/*****************************************************************************/

long OMPTestComp()
{
  long ptr, id, ok;

  /* Get pointer */

  ptr = (long)RDB[DATA_PTR_OMP_COMPLETED];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

#ifdef OPEN_MP
#pragma omp critical (comp)
#endif
  {
    /* Reset return value */

    ok = YES;

  /* Loop over threads */
    
    for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
      if (RDB[ptr + id] > 0.0)
        {
          /* Set return value */

          ok = NO;

          /* Break loop */

          break;
        }
  }

  /* Return value */

  return ok;
}

/*****************************************************************************/
