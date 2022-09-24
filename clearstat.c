/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : clearstat.c                                    */
/*                                                                           */
/* Created:       2010/11/10 (JLe)                                           */
/* Last modified: 2013/01/07 (JLe)                                           */
/* Version:       2.1.12                                                     */
/*                                                                           */
/* Description: Clears statistics in all bins                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ClearStat:"

/*****************************************************************************/

void ClearStat(long ptr)
{
  long nmax, n, stp;

  /* Check mode */

  if (ptr < 0)
    {
      /***********************************************************************/

      /***** Clear all statistics ********************************************/

      /* Loop over scores */

      ptr = (long)RDB[DATA_PTR_SCORE0];
      while(ptr > VALID_PTR)
        {
          /* Clear buffer */

          ClearStat(ptr);

          /* Next */
          
          ptr = NextItem(ptr);
        }

      /* Reduce RES2 array */

      ReducePrivateRes();

      /* Collect data for variance reduction routines */

      CollectVRMeshData();

      /* Clear private results */

      ClearPrivateRes();

      /* Clear buffer - tää pitää tehdä että batching */
      /* menee oikein 7.1.2013 / 2.1.12 (JLe) */ 

      ClearBuf();

      /* Reset batching counter */

      WDB[DATA_BATCH_COUNT] = 0.0;

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Clear single bin ************************************************/

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "", DATA_ARRAY, ptr);
      
      /* Get size */

      nmax = (long)RDB[ptr + SCORE_STAT_SIZE];

      if ((long)RDB[ptr + SCORE_DIM] < 1)
        Die(FUNCTION_NAME, "Tälle pitää tehdä jotain");

      /* Loop over values */
      
      for (n = 0; n < nmax; n++)
        {
          /* Get pointer to statistics */
          
          stp = (long)RDB[ptr + SCORE_PTR_DATA] + n*STAT_BLOCK_SIZE;      

          /* Clear values */

          memset(&RES1[stp], 0.0, STAT_BLOCK_SIZE*sizeof(double));          
        }

      /* TODO: Clear history */
      
      /***********************************************************************/
    }
}

/*****************************************************************************/
