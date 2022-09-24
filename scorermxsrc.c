/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorermxsrc.c                                  */
/*                                                                           */
/* Created:       2016/04/23 (JLe)                                           */
/* Last modified: 2020/03/04 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Scores source term for response matrix solver                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreRMXSrc:"

/*****************************************************************************/

void ScoreRMXSrc(long part, double x, double y, double z, double E, double wgt,
                 long id)
{
  long rmx, msh, loc0, ptr, ng, n, idx, ptr0, ptr1;

  /* Check if response matrix calculation is on */

  if ((long)RDB[DATA_RMTX_CALC] == NO)
    return;

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    return;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to weight window structure */

  rmx = (long)RDB[DATA_PTR_RMX0];
  CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

  /* Check particle type */

  if ((long)RDB[part + PARTICLE_TYPE] != (long)RDB[rmx + RMX_PARTICLE_TYPE])
    return;

  /* Get number of energy groups */

  ng = (long)RDB[rmx + RMX_NG];
  CheckValue(FUNCTION_NAME, "ng", "", ng, 1, 100000);

  /* Get pointer to few-group structure */

  if ((ptr = (long)RDB[rmx + RMX_PTR_EGRID]) > VALID_PTR)
    {
      /* Get group index */

      if ((n = GridSearch(ptr, E)) < 0)
        Error(rmx, "Particle energy %1.2E beyond energy boundaries", E);
    }
  else
    n = 0;

  /* Check value */

  CheckValue(FUNCTION_NAME, "n", "", n, 0, ng - 1);

  /* Pointer to mesh */

  if ((msh = (long)RDB[rmx + RMX_PTR_MESH]) < VALID_PTR)
    Die(FUNCTION_NAME, "No mesh");

  /* Reset event index */

  idx = 0.0;

  /* Find mesh cell */

  if ((loc0 = MeshPtr(msh, x, y, z)) > VALID_PTR)
    {
      /* Pointer to structure */

      loc0 = (long)RDB[loc0];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Store pointers and variables */

      WDB[part + PARTICLE_ICM_PTR_ICM] = (double)loc0;
      WDB[part + PARTICLE_ICM_IDX] = -1.0;
      WDB[part + PARTICLE_ICM_G] = (double)n;

      /* Score total source */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr + n, wgt, id);

      /* Check mode */

      if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_WWG)
        {
          /* WWG modes */

          ptr0 = (long)RDB[loc0 + RMX_CELL_MC_WWG_SRC0];
          CheckPointer(FUNCTION_NAME, "(ptr0)", RES2_ARRAY, ptr0);

          ptr1 = (long)RDB[loc0 + RMX_CELL_MC_WWG_SRC1];
          CheckPointer(FUNCTION_NAME, "(ptr1)", RES2_ARRAY, ptr1);

          idx = StoreRMXEvent(rmx, ptr0 + n, ptr1 + n, wgt, idx, id);

          /* Put event index */

          WDB[part + PARTICLE_ICM_EVENT_IDX] = (double)idx;
        }
    }
  else
    Error(rmx, "Source point is outside mesh");
}

/*****************************************************************************/
