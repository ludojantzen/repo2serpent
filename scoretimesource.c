/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoretimesource.c                              */
/*                                                                           */
/* Created:       2018/09/28 (VVa)                                           */
/* Last modified: 2018/09/28 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores source rates from different components in time-       */
/*              dependent calculations.                                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreTimeSource:"

/*****************************************************************************/

void ScoreTimeSource(long part, long mat, long mt,
                     double x0, double y0, double z0, double E0, double t0,
                     double val, double wgt, long id)
{
  long det, loc0, idx, rbin, ptr, type;

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Check particle type and super-imposed */

      if (((long)RDB[det + DET_PARTICLE] != type) ||
          ((long)RDB[det + DET_PTR_SBINS] > VALID_PTR))
        {
          /* Next detector */

          det = NextItem(det);

          /* Cycle loop */

          continue;
        }

      /* Reset response index */

      rbin = 0;

      /* Get pointer to response functions */

      loc0 = (long)RDB[det + DET_PTR_RBINS];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Loop over responses */

      while (loc0 > VALID_PTR)
        {
          /* Check mt */

          if ((long)RDB[loc0 + DET_RBIN_MT] == mt)
            {
              /* Get bin index */

              if ((idx = DetBin(det, mat, part, x0, y0, z0, E0, t0, id)) < 0)
                break;

              /* Get pointer to statistics */

              ptr = (long)RDB[det + DET_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Score */

              AddBuf(val, wgt, ptr, id, -1, idx, rbin);

              /* Set flags */

              SetDetFlags(det, part);
            }

          /* Update response index */

          rbin++;

          /* Next response */

          loc0 = NextItem(loc0);
        }

      /* Next detector */

      det = NextItem(det);
    }
}

/*****************************************************************************/
