/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : clearprivatedata.c                             */
/*                                                                           */
/* Created:       2011/11/10 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Sets all values to zero in PRIVA block                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ClearPrivateData:"

/*****************************************************************************/

void ClearPrivateData(long ptr)
{
  long i, sz;

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "PRIVA array not ready for access");

  /* Get size of data block */

  sz = (long)RDB[DATA_REAL_PRIVA_SIZE];

  /* Check pointer */

  if ((ptr < 0) || (ptr > sz - 1))
    Die(FUNCTION_NAME, "Pointer error");

  /* Loop over OpenMP threads and set values to zero */

  for (i = 0; i < (long)RDB[DATA_OMP_MAX_THREADS]; i++)
    PRIVA[ptr + i*sz] = 0.0;
}

/*****************************************************************************/
