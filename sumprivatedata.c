/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sumprivatedata.c                               */
/*                                                                           */
/* Created:       2011/11/10 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Calculates sum of values from PRIVA block                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SumPrivateData:"

/*****************************************************************************/

double SumPrivateData(long ptr)
{
  long i, sz;
  double sum;

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "PRIVA array not ready for access");

  /* Get size of data block */

  sz = (long)RDB[DATA_REAL_PRIVA_SIZE];

  /* Check pointer */

  if ((ptr < 0) || (ptr > sz - 1))
    Die(FUNCTION_NAME, "Pointer error");

  /* Reset sum */

  sum = 0.0;  

  /* Loop over OpenMP threads and add to sum */

  for (i = 0; i < (long)RDB[DATA_OMP_MAX_THREADS]; i++)
    sum = sum + PRIVA[ptr + i*sz];

  /* Return sum */

  return sum;
}

/*****************************************************************************/
