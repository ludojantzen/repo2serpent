/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cleaerprivateres.c                             */
/*                                                                           */
/* Created:       2011/11/11 (JLe)                                           */
/* Last modified: 2018/03/14 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Clears values in RES2 and RES3 data blocks                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ClearPrivateRes:"

/*****************************************************************************/

void ClearPrivateRes()
{
  long sz, max, i, nseg;

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "RES2 array not ready for access");

  /* Get buffer segment and data size */

  sz = (long)RDB[DATA_REAL_RES2_SIZE];
  max = (long)RDB[DATA_ALLOC_RES2_SIZE];

  /* Number of segments */
  
  if ((long)RDB[DATA_OPTI_SHARED_RES2] == YES)
    nseg = 1;
  else
    nseg = (long)RDB[DATA_OMP_MAX_THREADS];

  /* Loop over segments */

  for (i = 0; i < nseg; i++)
    memset(&RES2[i*sz], 0.0, max*sizeof(double));

  /* Check domain decomposition */
  
  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    {
      /* Get buffer size */

      max = (long)RDB[DATA_ALLOC_RES3_SIZE];

      /* Reset data */

      memset(RES3, 0.0, max*sizeof(double));
    }
}

/*****************************************************************************/
