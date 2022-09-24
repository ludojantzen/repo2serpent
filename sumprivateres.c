/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sumprivateres.c                                */
/*                                                                           */
/* Created:       2011/11/15 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Calculates sum of values from RES2 block                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SumPrivateRes:"

/*****************************************************************************/

double SumPrivateRes(long ptr)
{
  long i, sz;
  double sum;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "RES2 array not ready for access");
  
  /* Check if shared or private */

  if ((long)RDB[DATA_OPTI_SHARED_RES2] == YES)
    {
      /* Shared array, get data */
      
      sum = RES2[ptr];
    }
  else
    {
      /* Get size of data block */
  
      sz = (long)RDB[DATA_REAL_RES2_SIZE];

      /* Reset sum */

      sum = 0.0;
  
      /* Loop over OpenMP threads and add to sum */
  
      for (i = 0; i < (long)RDB[DATA_OMP_MAX_THREADS]; i++)
        sum = sum + RES2[ptr + i*sz];
    }

  /* Return sum */

  return sum;
}

/*****************************************************************************/
