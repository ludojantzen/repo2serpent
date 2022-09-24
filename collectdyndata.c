/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectdyndata.c                               */
/*                                                                           */
/* Created:       2012/11/11 (JLe)                                           */
/* Last modified: 2019/03/16 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Collects interval data from dynamic simulation               */

/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectDynData:"

/*****************************************************************************/

void  CollectDynData()
{
  long ptr, tb, id;
  double wgt0, wgt1, dt;

  /* Check time dependence */

  if ((RDB[DATA_TIME_CUT_TMAX] == INFTY) || (RDB[DATA_DYN_DT] > ZERO))
    return;

  /* Get time bin index and interval length */
  
  tb = (long)RDB[DATA_DYN_TB];
  dt = RDB[DATA_TIME_CUT_TMAX] - RDB[DATA_TIME_CUT_TMIN];
    
  /* Get previous weight */

  wgt0 = RDB[DATA_DYN_WGT0];

  /* Reset new weight */

  wgt1 = 0.0;

  /* Loop over threads */
  
  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get pointer to bank */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to first item after dummy */

      ptr = NextItem(ptr);

      /* Loop */
      
      while (ptr > VALID_PTR)
        {
          /* Add to weight */

          wgt1 = wgt1 + RDB[ptr + PARTICLE_WGT];
          
          /* Next */

          ptr = NextItem(ptr);
        }
    }

  /* Score population */

  ptr = (long)RDB[RES_DYN_POP];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddStat(wgt1/RDB[DATA_SRC_POP], ptr, tb);

  /* Score period */

  ptr = (long)RDB[RES_DYN_PERIOD];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  if (wgt0 != wgt1)
    AddStat(dt/log(wgt1/wgt0), ptr, tb);
  else
    AddStat(0.9999999*INFTY, ptr, tb);

  /* Put new weight */

  WDB[DATA_DYN_WGT0] = wgt1;
}

/*****************************************************************************/
