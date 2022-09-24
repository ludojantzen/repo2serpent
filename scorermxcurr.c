/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorermxcurr.c                                 */
/*                                                                           */
/* Created:       2016/04/23 (JLe)                                           */
/* Last modified: 2020/03/04 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Scores interface currents for response matrix solver         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreRMXCurr:"

/*****************************************************************************/

void ScoreRMXCurr(long part, double x, double y, double z, double E,
                  double wgt, long id)
{
  long rmx, msh, new, old, nmax, i, j, n, m, ptr, ng, idx, ptr0, ptr1, mode;

  /* Check if response matrix calculation is on */

  if ((long)RDB[DATA_RMTX_CALC] == NO)
    return;

  /* Check active cycle and corrector step */

  if ((RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]) ||
      ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
    return;

  /* Pointer to weight window structure */

  rmx = (long)RDB[DATA_PTR_RMX0];
  CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

  /* Check particle type */

  if ((long)RDB[part + PARTICLE_TYPE] != (long)RDB[rmx + RMX_PARTICLE_TYPE])
    return;

  /* Pointer to mesh */

  if ((msh = (long)RDB[rmx + RMX_PTR_MESH]) < VALID_PTR)
    Die(FUNCTION_NAME, "No mesh");

  /* Get mode (ei määritetty convergence acceleration -moodissa) */

  mode = (long)RDB[rmx + RMX_MODE];
  CheckValue(FUNCTION_NAME, "mode", "", mode, 0, 4);

  /* Find mesh cell */

  if ((new = MeshPtr(msh, x, y, z)) < VALID_PTR)
    Error(rmx, "Mesh must cover the entire geometry");

  /* Pointer to structure */

  new = (long)RDB[new];
  CheckPointer(FUNCTION_NAME, "(new)", DATA_ARRAY, new);

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Get old mesh cell pointer */

  old = (long)RDB[part + PARTICLE_ICM_PTR_ICM];
  CheckPointer(FUNCTION_NAME, "(old)", DATA_ARRAY, old);

  /* Get maximum number of neighbours */

  nmax = (long)RDB[old + RMX_CELL_MTX_SIZE];
  CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /* Get pointer to few-group structure */

  if ((ptr = (long)RDB[rmx + RMX_PTR_EGRID]) > VALID_PTR)
    {
      /* Get group index */

      if ((m = GridSearch(ptr, E)) < 0)
        Error(rmx, "Particle energy %1.2E beyond energy boundaries", E);
    }
  else
    m = 0;

  /* Get entry group */

  n = (long)RDB[part + PARTICLE_ICM_G];
  CheckValue(FUNCTION_NAME, "n", "", n, 0, ng - 1);

  /* Check value */

  CheckValue(FUNCTION_NAME, "m", "", m, 0, ng - 1);

  /* Get entry index */

  i = (long)RDB[part + PARTICLE_ICM_IDX];
  CheckValue(FUNCTION_NAME, "i", "", i, -2, nmax);

  /* Find neighbour cell for old: i -> j */

  ptr = (long)RDB[old + RMX_CELL_PTR_BOUNDS];
  for (j = 0; j < nmax; j++)
    {
      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

      /* Compare */

      if ((long)RDB[ptr + RMX_CELL_BOUND_PTR_CELL] == new)
        break;

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Check mode */

  if (mode != RMX_MODE_WWG)
    {
      /* Check if found */

      if (j < nmax)
        {
          /* Check if particle came from source or another boundary */

          if (i == -1)
            {
              /* Score source current vector */

              ptr = (long)RDB[old + RMX_CELL_MC_CURR_SRCC_OUT];
              CheckPointer(FUNCTION_NAME, "(ptr2)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + j*ng*ng + n*ng + m, wgt, id);
            }
          else if (i > -1)
            {
              /* Score passing currents */

              ptr = (long)RDB[old + RMX_CELL_MC_CURR_THROUGH];
              CheckPointer(FUNCTION_NAME, "(ptr3)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + i*nmax*ng*ng + n*nmax*ng + j*ng + m, wgt, id);
            }

          /* Score net outward current (-2 meinaa failurea) */

          if (i > -2)
            {
              ptr = (long)RDB[old + RMX_CELL_MC_CURR_NET_OUT];
              CheckPointer(FUNCTION_NAME, "(ptr5)", RES2_ARRAY, ptr);
              AddPrivateRes(ptr + j*ng + m, wgt, id);
            }
        }
      else if (new != old)
        {
          /* Score failure */

          ptr = (long)RDB[RES_RMX_CURR_SEARCH_FAIL];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, 1.0, ptr, id, 0);
        }

      /* Score counter */

      ptr = (long)RDB[RES_RMX_CURR_SEARCH_FAIL];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, 1.0, ptr, id, 1);
    }

  /* Check if boundary was not crossed */

  if (new == old)
    {
      /* This may happen if weight windows are used with different mesh    */
      /* structure (wwdis.c returns the shortest distance to both meshes). */
      /* Ingore for now, but better option would be to set up a flag that  */
      /* indicates the type of mesh boundary where the track was stopped.  */

      if (((long)RDB[DATA_RMX_CONVG_ACC] == YES) ||
          ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == YES))
        return;

      /* Check that point is not on boundary (NOTE: tää testi ei tuu */
      /* toimimaan adaptiivisen meshin kanssa) */

      if (MeshIndex(msh, x, y, z, 1E-6) == -2)
        Error(rmx, "Source point is on or too close to mesh boundary");
      else
        Die(FUNCTION_NAME, "Boundary not crossed %E %E %E : %f %f", x, y, z,
            sqrt(x*x + y*y), PolarAngle(x, y)*180/PI);
    }

  /* Get maximum number of neighbours */

  nmax = (long)RDB[new + RMX_CELL_MTX_SIZE];
  CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

  /* Find neighbour cell for new: j -> i */

  ptr = (long)RDB[new + RMX_CELL_PTR_BOUNDS];
  for (i = 0; i < nmax; i++)
    {
      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(ptr4)", DATA_ARRAY, ptr);

      /* Compare */

      if ((long)RDB[ptr + RMX_CELL_BOUND_PTR_CELL] == old)
        break;

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Store pointers and variables */

  WDB[part + PARTICLE_ICM_PTR_ICM] = (double)new;

  /* Tätä ei käytetä RMX-laskennassa? */
  /*
  WDB[part + PARTICLE_ICM_WGT] = wgt;
  */

  WDB[part + PARTICLE_ICM_G] = (double)m;

  if (i < nmax)
    {
      /* Score net inward current (HUOM! toi on naapuriin menevä virta, */
      /* sen takia ryhmä m eikä n). */

      if ((long)RDB[rmx + RMX_MODE] != RMX_MODE_WWG)
        {
          /* RMX modes */

          ptr = (long)RDB[new + RMX_CELL_MC_CURR_NET_IN];
          CheckPointer(FUNCTION_NAME, "(ptr5)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + i*ng + m, wgt, id);
        }
      else
        {
          /* WWG modes */

          idx = (long)RDB[part + PARTICLE_ICM_EVENT_IDX];
          CheckValue(FUNCTION_NAME, "idx", "", idx, 1, MAX_RMX_BUFF);

          /* Score currents */

          ptr0 = (long)RDB[new + RMX_CELL_MC_WWG_CURR0];
          CheckPointer(FUNCTION_NAME, "(ptr0)", RES2_ARRAY, ptr0);

          ptr1 = (long)RDB[new + RMX_CELL_MC_WWG_CURR1];
          CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr1);

          idx = StoreRMXEvent(rmx, ptr0 + i*ng + m, ptr1 + i*ng + m, wgt,
                              idx, id);

          /* Put event index */

          WDB[part + PARTICLE_ICM_EVENT_IDX] = (double)idx;
        }

      /* Put new index */

      WDB[part + PARTICLE_ICM_IDX] = (double)i;
    }
  else
    {
      /* Put -2 to indicate index fail */

      WDB[part + PARTICLE_ICM_IDX] = -2.0;
    }
}

/*****************************************************************************/
