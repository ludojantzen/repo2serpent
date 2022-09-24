/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorepb.c                                      */
/*                                                                           */
/* Created:       2011/11/20 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Scores power distributions for explicit geometry             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScorePB:"

/*****************************************************************************/

void ScorePB(double fissE, double wgt, long id)
{
  long pbd, pbl, ptr, ncol, idx;

  /* Check pointer */

  if ((pbd = (long)RDB[DATA_PTR_PB0]) < 0)
    return;

  /* Check fission energy */

  if (fissE < ZERO)
    return;

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Loop over geometries */

  while (pbd > 0)
    {
      /* Get pointer to pebble where collision occured */

      if ((ptr = (long)RDB[pbd + PBED_PTR_POW]) > VALID_PTR)
        if ((pbl = (long)TestValuePair(pbd + PBED_PTR_COL_PEBBLE, 
                                       (double)ncol, id)) > VALID_PTR)
          {
            /* Get pebble index */
            
            idx = (long)RDB[pbl + PEBBLE_IDX];
            
            /* Score */
            
            AddBuf1D(fissE, wgt, ptr, id, idx);
          }
      
      /* Next */
      
      pbd = NextItem(pbd);
    }
}

/*****************************************************************************/
