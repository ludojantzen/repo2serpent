/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : weightwindow.c                                 */
/*                                                                           */
/* Created:       2012/04/13 (JLe)                                           */
/* Last modified: 2019/09/27 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Handles particle hit in weight window                        */
/*                                                                           */
/* Comments: - HUOM: painokatkaistuja historioita ei lasketa mukaan          */
/*             analogisiin detektoreihin sen vuoksi että osa hiukkasista     */
/*             saa ruletissa suuremman painon, ja keskimääräinen paino       */
/*             säilyy. Simulaation kannalta painokatkaisu toimii siis eri    */
/*             tavalla kuin esim. energiakatkaisu.                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WeightWindow:"

/*****************************************************************************/

long WeightWindow(long trk, long part, long type, double x, double y, double z,
                  double u, double v, double w, double E, double *wgt,
                  double t, long itype, long id)
{
  long np, new, ptr, cell;
  double min, max, p, a;

  /* Check if weight windows are used */

  if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO)
    return trk;

  /* Get importance */

  if ((p = WWImportance(type, x, y, z, u, v, w, E, itype)) == 0.0)
    {
      /* Check if geometry importances are used */

      if ((long)RDB[DATA_USE_GEOM_IMP] == NO)
        return trk;

      /* Check cell pointer */

      if ((cell = WhereAmI(x, y, z, u, v, w, id)) < VALID_PTR)
        return trk;

      /* Get cell importance */

      if ((p = RDB[cell + CELL_IMP]) <= 0.0)
        return trk;

      /* Cell importances and weight windows should not be used */
      /* simultaneously */

      Die(FUNCTION_NAME, "Should not be here");
    }

  /* Boundaries */

  min = RDB[DATA_WWD_LOWER_BOUND]/p;
  max = RDB[DATA_WWD_UPPER_BOUND]/p;

  /* Compare weight to boundaries */

  if (*wgt < min)
    {
      /***********************************************************************/

      /***** Play russian roulette *******************************************/

      /* Calculate survival probability */

      p = *wgt/min;

      /* Multiply by user-given factor */

      p = p*RDB[DATA_WWD_SURVIVAL_F];

      /* Compare to minimum */

      if (p < RDB[DATA_WWD_MIN_ROULETTE])
        p = RDB[DATA_WWD_MIN_ROULETTE];

      /* Russian roulette */

      if (RandF(id) < p)
        {
          /* Score particle balance */

          if (type == PARTICLE_TYPE_NEUTRON)
            {
              ptr = (long)RDB[RES_N_BALA_SRC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(*wgt*(1.0/p - 1.0), 1.0, ptr, id, -1, BALA_N_SRC_VR, 1);
            }
          else if (type == PARTICLE_TYPE_GAMMA)
            {
              ptr = (long)RDB[RES_G_BALA_SRC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(*wgt*(1.0/p - 1.0), 1.0, ptr, id, -1, BALA_G_SRC_VR, 1);
            }

          /* Score weight balance */

          ptr = (long)RDB[RES_WW_BALA_ROULETTE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(*wgt*(1.0/p - 1.0), 1.0, ptr, id, 2 - type);

          /* Increase weight */

          *wgt = *wgt/p;
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

          /* Score weight balance */

          ptr = (long)RDB[RES_WW_BALA_ROULETTE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(-(*wgt), 1.0, ptr, id, 2 - type);

          /* Return weight cut-off */

          return TRACK_END_WCUT;
        }

      /***********************************************************************/
    }
  else if (*wgt > max)
    {
      /***********************************************************************/

      /***** Split history ***************************************************/

      /* Two methods */

      if (1 == 2)
        {
          /* Calculate number of splits */

          np = (long)(*wgt/max);

          /* Divide weight */

          *wgt = *wgt/((double)(np + 1));
        }
      else
        {
          /* Calculate expected number of particles */

          a = *wgt/max;

          /* Truncate to integer value */

          np = (long)a;

          /* Sample additional */

          if (RandF(id) < a - (double)np)
            np++;

          /* Check excessive splitting */

          if (np > (long)RDB[DATA_WWD_MAX_SPLIT])
            {
              /* Check particle type */
              /*
              if (type == PARTICLE_TYPE_NEUTRON)
                Note(0,
                     "Excessive neutron splitting at weight window boundary");
              else
                Note(0,
                     "Excessive photon splitting at weight window boundary");
              */
              /* Truncate */

              np = (long)RDB[DATA_WWD_MAX_SPLIT];
            }

          /* Score weight balance */

          if (type == PARTICLE_TYPE_NEUTRON)
            {
              ptr = (long)RDB[RES_N_BALA_SRC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(np*max - *wgt, 1.0, ptr, id, -1, BALA_N_SRC_VR, 1);
            }
          else if (type == PARTICLE_TYPE_GAMMA)
            {
              ptr = (long)RDB[RES_G_BALA_SRC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(np*max - *wgt, 1.0, ptr, id, -1, BALA_G_SRC_VR, 1);
            }

          /* Score weight balance */

          ptr = (long)RDB[RES_WW_BALA_SPLIT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(np*max - *wgt, 1.0, ptr, id, 2 - type);

          /* Set new weight */

          *wgt = max;

          /* Remove original */

          np--;
        }

      /* Check number of splits */

      if (np > 0)
        {
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

  /* Check */

  if (RDB[DATA_WWD_SURVIVAL_F] == 1.0)
    if ((*wgt < 0.99999*min) || (*wgt > 1.00001*max))
      Die(FUNCTION_NAME, "Weight is out of bounds : %E %E %E", min,
          *wgt, max);

  /* Exit subroutine */

  return trk;
}

/****************************************************************************/
