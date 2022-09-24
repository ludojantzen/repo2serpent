/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ompsetcomp.c                                   */
/*                                                                           */
/* Created:       2019/02/26 (JLe)                                           */
/* Last modified: 2019/02/26 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Set thread-wise completed flag                               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "OMPSetComp:"

/*****************************************************************************/

void OMPSetComp(long id)
{
  long ptr;

  /* Get pointer */

  ptr = (long)RDB[DATA_PTR_OMP_COMPLETED];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

#ifdef OPEN_MP
#pragma omp critical (comp)
#endif
  {
    /* Put flag */
    
    WDB[ptr + id] = 1.0;
  }
}

/*****************************************************************************/
