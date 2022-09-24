/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ompresetcomp.c                                 */
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

#define FUNCTION_NAME "OMPResetComp:"

/*****************************************************************************/

void OMPResetComp(long id)
{
  long ptr;

  /* Get pointer */

  ptr = (long)RDB[DATA_PTR_OMP_COMPLETED];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

#ifdef OPEN_MP
#pragma omp critical (comp)
#endif
  {
    /* Check if thread has been given */

    if (id > -1)
      WDB[ptr + id] = 0.0;
    else
      {
        /* Loop over threads */
        
        for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
          WDB[ptr + id] = 0.0;
      }
  }
}

/*****************************************************************************/
