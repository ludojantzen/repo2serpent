/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addprivateres.c                                */
/*                                                                           */
/* Created:       2011/11/11 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Adds value in RES2 data block                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddPrivateRes:"

/*****************************************************************************/

void AddPrivateRes(long ptr, double val, long id)
{
  long sz;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "RES2 array not ready for access");
  
  /* Check if shared or private */

  if ((long)RDB[DATA_OPTI_SHARED_RES2] == YES)
    {
      /* Shared array, put data */

#ifdef OPEN_MP
#pragma omp atomic
#endif
      RES2[ptr] += val;
    }
  else
    {
      /* Get size of data block */
  
      sz = (long)RDB[DATA_REAL_RES2_SIZE];
  
      /* Put value */
  
      RES2[ptr + id*sz] = RES2[ptr + id*sz] + val;
    }
}

/*****************************************************************************/
