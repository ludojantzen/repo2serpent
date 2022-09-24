/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : leak.c                                         */
/*                                                                           */
/* Created:       2011/03/12 (JLe)                                           */
/* Last modified: 2017/08/15 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Handles leakage                                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Leak:"

/*****************************************************************************/

void Leak(long part, double x, double y, double z, double u, double v, 
          double w, double E, double wgt, long id)
{
  long ptr, type, n0, ma0, ms0, ng0, icm;
  double wgt0;

  /* Check pointer */
      
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  
  /* Put particle back in stack */

  ToStack(part, id);

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Score total leak rate */

  if (type == PARTICLE_TYPE_NEUTRON)
    ptr = (long)RDB[RES_TOT_NEUTRON_LEAKRATE];
  else
    ptr = (long)RDB[RES_TOT_PHOTON_LEAKRATE];

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(1.0, wgt, ptr, id, 0);

  /* Score particle balance */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      ptr = (long)RDB[RES_N_BALA_LOSS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf(1.0, 1.0, ptr, id, -1, BALA_N_LOSS_LEAK, 0);
      AddBuf(wgt, 1.0, ptr, id, -1, BALA_N_LOSS_LEAK, 1);
    }
  else if (type == PARTICLE_TYPE_GAMMA)
    {
      ptr = (long)RDB[RES_G_BALA_LOSS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_LOSS_LEAK, 0);
      AddBuf(wgt, 1.0, ptr, id, -1, BALA_G_LOSS_LEAK, 1);
    }

  /* Check active cycle and corrector step */

  if ((RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]) ||
      (type == PARTICLE_TYPE_GAMMA) ||
      ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
    return;

  /* Leak rates for six-factor formula */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      if (E < 0.625E-6)
        ptr = (long)RDB[RES_SIX_FF_LT];
      else
        ptr = (long)RDB[RES_SIX_FF_LF];
      
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, wgt, ptr, id, 0);
    }

  /***************************************************************************/
  
  /***** ICM leak rate *******************************************************/

  /* Get universe pointer */
  
  if ((long)RDB[DATA_ICM_CALC] == YES)
    if ((icm = (long)RDB[part + PARTICLE_ICM_PTR_ICM]) > VALID_PTR)
      {
        /* Get data (toi nollallinen muuttuja on jo käytössä) */
        
        n0 = (long)RDB[part + PARTICLE_ICM_IDX];
        ma0 = (long)RDB[part + PARTICLE_ICM_MUA];
        ms0 = (long)RDB[part + PARTICLE_ICM_MUS];
        ng0 = (long)RDB[part + PARTICLE_ICM_G];
        wgt0 = RDB[part + PARTICLE_ICM_WGT];
              
        /* Check */
        
        if (n0 > -1)
          {
            ptr = (long)RDB[icm + ICM_RES_LEAK1];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            AddBuf(1.0, wgt, ptr, id, -1, n0, ma0, ms0, ng0);
            
            ptr = (long)RDB[icm + ICM_RES_LEAK2];
            CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            AddBuf(wgt0, wgt, ptr, id, -1, n0, ma0, ms0, ng0);
          }
      }

  /***************************************************************************/
}

/*****************************************************************************/
