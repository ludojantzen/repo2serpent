/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorermxresp.c                                 */
/*                                                                           */
/* Created:       2016/04/24 (JLe)                                           */
/* Last modified: 2020/03/04 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Scores response term for response matrix solver              */
/*                                                                           */
/* Comments: - Ei toimi pintadetektorien kanssa jos pinta yhtyy meshin       */
/*             rajoihin.                                                     */
/*                                                                           */
/*           - Detektoriin flagi että käytetäänkö tässä vai ei.              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreRMXResp:"

/*****************************************************************************/

void ScoreRMXResp(long part, long det0, double res, long id)
{
  long rmx, loc0, det, ptr, i, ng, n, idx;
  double vol;

  /* Check if response matrix calculation is on */

  if ((long)RDB[DATA_RMTX_CALC] == NO)
    return;

  /* Check if mfp iteration */

  if ((long)RDB[DATA_RMTX_MFP_CALC] == YES)
    return;

  /* Check if source convergence acceleraion is on */

  if ((long)RDB[DATA_RMX_CONVG_ACC] == YES)
    return;

  /* Check active cycle and corrector step */

  if ((RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]) ||
      ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
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

  /* Get pointer */

  loc0 = (long)RDB[part + PARTICLE_ICM_PTR_ICM];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Check if TLE mode */

  if ((long)RDB[rmx + RMX_USE_TLE] == YES)
    {
      /* Detector pointer is set to zero for TLE */

      if (det0 < VALID_PTR)
        {
          /* Use the first detector */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);
        }
      else
        {
          /* Do not score other contributions */

          return;
        }
    }
  else
    {
      /* Do not score TLE */

      if (det0 < VALID_PTR)
        return;

      /* Loop over detectors and find match */

      det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
      while (det > VALID_PTR)
        {
          /* Compare pointer */

          if ((long)RDB[det + RMX_DET_PTR_DET] == det0)
            break;

          /* Next item */

          det = NextItem(det);
        }

      /* Check */

      if (det < VALID_PTR)
        return;
    }

  /* Get detector volume */

  if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_GVR)
    vol = RDB[loc0 + RMX_CELL_RVOL];
  else if ((vol = RDB[det0 + DET_VOL]) == 0.0)
    vol = 1.0;

  /* Absolute values are needed for surface current detector responses */

  res = fabs(res)/vol;

  /* Check wwg mode */

  if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_WWG)
    {
      /* Get event index */

      idx = (long)RDB[part + PARTICLE_ICM_EVENT_IDX];
      CheckValue(FUNCTION_NAME, "idx", "", idx, 1, MAX_RMX_BUFF);

      /* Collect buffered events */

      StoreRMXEvent(rmx, -1, -1, res, idx, id);

      /* Exit subroutine */

      return;
    }

  /* Get energy group */

  n = (long)RDB[part + PARTICLE_ICM_G];
  CheckValue(FUNCTION_NAME, "n", "", n, 0, ng - 1);

  /* Get index */

  if ((i = (long)RDB[part + PARTICLE_ICM_IDX]) == -1)
    {
      /* Particle came from source, score direct contributions */

      ptr = (long)RDB[det + RMX_DET_MC_RES_SRCC];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr + n, res, id);

      /* Score total contribution */

      ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr + n, res, id);

      /* Add to number of scores */

      ptr = (long)RDB[det + RMX_DET_MC_RES_SCORE_N];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, 1.0, id);
    }
  else if (i > -1)
    {
      /* Particle came from boundary */

      ptr = (long)RDB[det + RMX_DET_MC_RES_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr + i*ng + n, res, id);

      /* Score total contribution */

      ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr + n, res, id);

      /* Add to number of scores */

      ptr = (long)RDB[det + RMX_DET_MC_RES_SCORE_N];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, 1.0, id);
    }
}

/*****************************************************************************/
