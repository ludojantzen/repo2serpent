/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorepinpower.c                                */
/*                                                                           */
/* Created:       2013/08/05 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Scores power distributions for pin-power reconstruction      */
/*                                                                           */
/* Comments: - NOTE: noissa micro-group -muuttujissa on macro-group jako     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScorePinPower:"

/*****************************************************************************/

void ScorePinPower(long mat, long part, double E, double fissE, double wgt, 
                   long id)
{
  long ppw, ptr, lat, gcu, uni, ncol, n, ntot, np, ng;
  double x, y, z;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Check if active cycle */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Check if calculation is required */
  
  if ((ppw = (long)RDB[DATA_PTR_PPW0]) < VALID_PTR)
    return;

  /* Check fission energy */

  if (fissE == 0.0)
    return;

  /* Number of energy groups */
  
  ntot = (long)RDB[DATA_ERG_FG_NG];

  /* Get pointer to few-group structure */
  
  ptr = (long)RDB[DATA_ERG_FG_PTR_GRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Get few-group index */
              
  if ((ng = GridSearch(ptr, E)) > -1)
    {
      /* Convert index */
      
      ng = ntot - ng - 1;
      CheckValue(FUNCTION_NAME, "ng", "", ng, 0, ntot - 1);
    }
  else
    return;

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Loop over pin-power distributions */

  while (ppw > VALID_PTR)
    {
      /* Pointer to gcu universe */
      
      gcu = (long)RDB[ppw + PPW_PTR_GCU];
      CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);
      
      /* Pointer to universe */
      
      uni = (long)RDB[gcu + GCU_PTR_UNIV];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Number of pins */
      
      np = (long)RDB[ppw + PPW_NP];

      /* Compare collision number */
      
      if (TestValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, id) > -1)
        {
          /* Pointer to lattice */

          lat = (long)RDB[ppw + PPW_PTR_LAT];
          CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);
  
          /* Get index */
      
          if ((n = (long)TestValuePair(lat + LAT_COL_CELL_IDX, 
                                       (double)ncol, id)) > -1)
            {
              /* Score power (group-wise and total) */

              ptr = (long)RDB[gcu + GCU_MICRO_PPW_POW];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + n*(ntot + 1) + ng, fissE*wgt, id);
              AddPrivateRes(ptr + n*(ntot + 1) + ntot, fissE*wgt, id);

              /* Get coordinates */
              
              ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_X];
              CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
              x = GetPrivateData(ptr, id);
              
              ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Y];
              CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
              y = GetPrivateData(ptr, id);
              
              ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Z];
              CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
              z = GetPrivateData(ptr, id);
              
              /* Score coordinates (group-wise and total) */

              ptr = (long)RDB[gcu + GCU_MICRO_PPW_XYZ];
              CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + 0*np*(ntot + 1) + n*(ntot + 1) + ng, 
                            fissE*wgt*x, id);
              AddPrivateRes(ptr + 1*np*(ntot + 1) + n*(ntot + 1) + ng, 
                            fissE*wgt*y, id);
              AddPrivateRes(ptr + 2*np*(ntot + 1) + n*(ntot + 1) + ng, 
                            fissE*wgt*z, id);
              AddPrivateRes(ptr + 0*np*(ntot + 1) + n*(ntot + 1) + ntot, 
                            fissE*wgt*x, id);
              AddPrivateRes(ptr + 1*np*(ntot + 1) + n*(ntot + 1) + ntot, 
                            fissE*wgt*y, id);
              AddPrivateRes(ptr + 2*np*(ntot + 1) + n*(ntot + 1) + ntot, 
                            fissE*wgt*z, id);
            }
        }
      
      /* Next */

      ppw = NextItem(ppw);
    }
}

/*****************************************************************************/
