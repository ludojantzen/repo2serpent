/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : geoimportance.c                                */
/*                                                                           */
/* Created:       2018/11/20 (JLe)                                           */
/* Last modified: 2019/02/20 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Handles splitting by geometry importance                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GeoImportance:"

/*****************************************************************************/

long GeoImportance(long trk, long part, long cell, double x, double y, 
                   double z, double u, double v, double w, double E, 
                   double *wgt, double t, long id)
{
  long np, new, ptr, type;
  double I0, I1, r;

  /* Check if importances are used */

  if ((long)RDB[DATA_USE_GEOM_IMP] == NO)
    return trk;

  /* Weight windows override importances */

  if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == YES)
    return trk;

  /* Check cell pointer */

  if (cell < VALID_PTR)
    return trk;

  /* Get previous importance */

  I0 = RDB[part + PARTICLE_PREV_I];

  /* Get new importance */

  I1 = RDB[cell + CELL_IMP];
  CheckValue(FUNCTION_NAME, "I1", "", I1, 0.0, INFTY);

  /* Compare */

  if (I0 == I1)
    return trk;

  /* Store importance */

  WDB[part + PARTICLE_PREV_I] = I1;

  /* Check */

  if (I0 <= 0.0)
    {
      /* Exit */

      return trk;
    }

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Calculate ratio */

  r = I1/I0;
  CheckValue(FUNCTION_NAME, "r", "", r, ZERO, INFTY);

  /* Check */

  if (r < 1.0)
    {
      /***********************************************************************/
      
      /***** Play russian roulette *******************************************/

      /* Russian roulette */
      
      if (RandF(id) < r)
        {
          /* Score particle balance */

          if (type == PARTICLE_TYPE_NEUTRON)
            {
              ptr = (long)RDB[RES_N_BALA_SRC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(*wgt*(1.0/r - 1.0), 1.0, ptr, id, -1, BALA_N_SRC_VR, 1);
            }
          else if (type == PARTICLE_TYPE_GAMMA)
            {
              ptr = (long)RDB[RES_G_BALA_SRC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(*wgt*(1.0/r - 1.0), 1.0, ptr, id, -1, BALA_G_SRC_VR, 1);
            }
          
          /* Increase weight */

          *wgt = *wgt/r;
        }
      else
        {
          /* Put particle back in stack */
          
          if (part > VALID_PTR)
            ToStack(part, id);

          /* Score cut-off */

          if (type == PARTICLE_TYPE_NEUTRON)
            ptr = (long)RDB[RES_TOT_NEUTRON_CUTRATE];
          else
            ptr = (long)RDB[RES_TOT_PHOTON_CUTRATE];

          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, *wgt, ptr, id, 0);              

          /* Score particle balance */

          if (type == PARTICLE_TYPE_NEUTRON)
            {
              ptr = (long)RDB[RES_N_BALA_LOSS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(1.0, 1.0, ptr, id, -1, BALA_N_LOSS_CUT, 0);
              AddBuf(*wgt, 1.0, ptr, id, -1, BALA_N_LOSS_CUT, 1);
            }
          else if (type == PARTICLE_TYPE_GAMMA)
            {
              ptr = (long)RDB[RES_G_BALA_LOSS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(1.0, 1.0, ptr, id, -1, BALA_G_LOSS_CUT, 0);
              AddBuf(*wgt, 1.0, ptr, id, -1, BALA_G_LOSS_CUT, 1);
            }
          
          /* Return weight cut-off */
          
          return TRACK_END_WCUT;
        }
         
      /***********************************************************************/
    }
  else 
    {
      /***********************************************************************/
      
      /***** Split history ***************************************************/

      /* Truncate to integer value */
      
      np = (long)r;
      
      /* Sample additional */
      
      if (RandF(id) < r - (double)np)
        np++;
      
      /* Set new weight */
      
      *wgt = *wgt/((double)np);
      
      /* Remove original */
      
      np--;
      
      /* Check number of splits */
      
      if (np > 0)
        {
          /* Check excessive splitting */
          
          if (np > (long)RDB[DATA_WWD_MAX_SPLIT])
            {
              /* Check particle type */
              
              if (type == PARTICLE_TYPE_NEUTRON)
                Note(0, "Excessive neutron splitting at importance boundary");
              else
                Note(0, "Excessive photon splitting at importance boundary");
              
              /* Truncate */
              
              np = (long)RDB[DATA_WWD_MAX_SPLIT];              
            }
          
          /* Duplicate incident neutron */
      
          new = DuplicateParticle(part, id);
      
          /* Put state variables */
          
          WDB[new + PARTICLE_X] = x;
          WDB[new + PARTICLE_Y] = y;
          WDB[new + PARTICLE_Z] = z;
          
          WDB[new + PARTICLE_U] = u;
          WDB[new + PARTICLE_V] = v;
          WDB[new + PARTICLE_W] = w;
          
          WDB[new + PARTICLE_E] = E;
          WDB[new + PARTICLE_WGT] = *wgt;
          WDB[new + PARTICLE_T] = t;
          
          /* Set multiplicity */
          
          WDB[new + PARTICLE_MULTIPLICITY] = (double)(np - 1);
          
          /* Put particle in que */
          
          ToQue(new, id);
          
          /* Score particle balance */
          
          if (type == PARTICLE_TYPE_NEUTRON)
            {
              ptr = (long)RDB[RES_N_BALA_SRC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf((double)np, 1.0, ptr, id, -1, BALA_N_SRC_VR, 0);
            }
          else if (type == PARTICLE_TYPE_GAMMA)
            {
              ptr = (long)RDB[RES_G_BALA_SRC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf((double)np, 1.0, ptr, id, -1, BALA_G_SRC_VR, 0);
            }
        }
      
      /***********************************************************************/
    }
  
  /* Exit subroutine */

  return trk;
}

/****************************************************************************/
