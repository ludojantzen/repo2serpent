/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorecmm.c                                     */
/*                                                                           */
/* Created:       2014/06/25 (JLe)                                           */
/* Last modified: 2019/02/15 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Updates coordinates and records scores for cumulative        */
/*              migration method diffusion coefficients                      */
/*                                                                           */
/* Comments: - Vanha implementaatio poistettiin 15.2.2018 / 2.1.31 (JLe)     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreCMM:"

/*****************************************************************************/

/* Return YES if coordinates of last entering into group need to be updated,
   otherwise return NO. On input E1 will have the resulting neutron energy
   if scattering has occured. If E1 < 0.0, either capture, fission or energy
   or weight cut has occured. If E1 < -2.0 then energy or weight cut has
   occured. */

long ScoreCMM(double dx0, double dy0, double dz0, double dx, double dy,
              double dz, double E0, double E1, double wgt0, double wgt1, 
              long id)
{
  long ptr, ncol, gcu, g0, g1, uni, ng, upd;
  double dwgt;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return (long)NO;

  /* Check if active cycle */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return (long)NO;

  /* Check that CMM calculation is on */

  if ((long)RDB[DATA_CMM_CALC] == NO)
    return (long)NO;

  /* Default to not to update coordinates */

  upd = (long)NO;

  /* Get pointer to microgroup energy grid */
      
  ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);        
  
  /* Number of groups */
      
  ng = (long)RDB[ptr + ENERGY_GRID_NE] - 1;
 
  /* Get group index */
      
  if ((g0 = GridSearch(ptr, E0)) > -1)
    {
      /* Convert index */
      
      g0 = ng - g0 - 1;
      CheckValue(FUNCTION_NAME, "ng2", "", g0, 0, ng - 1);
    }
  else
    {
      /* Energy is out of bounds */

      return (long)NO;
    }

  /* Get out-going few-group index */

  if (E1 < ZERO)
    {
      /* Capture, fission or energy or weight cut (not scattering) */

      g1 = -2;

      if (E1 < -2.0)
        {
          /* Energy or weight cut */

          g1 = -3;
        }
    }
  else if ((g1 = GridSearch(ptr, E1)) > -1)
    {
      /* Convert index */

      g1 = ng - g1 - 1;
      CheckValue(FUNCTION_NAME, "g1", "", g1, 0, ng - 1);
    }
  else
    {
      /* Energy is out of bounds */

      return (long)NO;
    }

  /* Check if group-to-group scattering */

  if (g0 == g1)
    {
      /* Check if weight does not change */

      if (fabs(wgt0 - wgt1) < ZERO)
        return (long)NO;
    }

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Check for multiple levels */

  if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
    {
      /* Single level, get pointer */

      if ((gcu = (long)TestValuePair(DATA_GCU_PTR_UNI, (double)ncol, id))
          < VALID_PTR)
        return (long)NO;
    }
  else
    {
      /* Multiple levels, get pointer to list */

      gcu = (long)RDB[DATA_PTR_GCU0];
      CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);
    }

  /* Loop over universes */

  while (gcu > VALID_PTR)
    {
      /* Check multi-level mode */

      if ((long)RDB[DATA_MULTI_LEVEL_GCU] == YES)
        {
          /* Pointer to universe */

          uni = (long)RDB[gcu + GCU_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Check collision */

          if (TestValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, id) < 0.0)
            {
              /* Next universe */

              gcu = NextItem(gcu);

              /* Cycle loop */

              continue;
            }
        }

      /* Get stat pointers */

      ptr = (long)RDB[gcu + GCU_BUF_CMM_CUMU_R2];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

      /* If weight changes, score with the weight change if neutron stays
         in the group, otherwise score with the original weight. In case
         of energy or weight cut, score with the original weight */
      
      if (fabs(wgt0 - wgt1) > ZERO && ((g1 == -2) || (g0 == g1)))
        dwgt = wgt0 - wgt1;
      else
        {
          dwgt = wgt0;
          
          /* Update coordinates */
          
          upd = (long)YES;
        }
      
      /* Store change in position */
      
      AddPrivateRes(ptr + ng + g0, dwgt*(dx*dx + dy*dy + dz*dz -
                                         dx0*dx0 - dy0*dy0 - dz0*dz0), id);
      AddPrivateRes(ptr + 2*ng + g0, dwgt*(dx*dx - dx0*dx0), id);
      AddPrivateRes(ptr + 3*ng + g0, dwgt*(dy*dy - dy0*dy0), id);
      AddPrivateRes(ptr + 4*ng + g0, dwgt*(dz*dz - dz0*dz0), id);

      /* Next universe */

      gcu = NextItem(gcu);
    }

  /* Exit OK */

  return upd;
}

/******************************************************************************/
