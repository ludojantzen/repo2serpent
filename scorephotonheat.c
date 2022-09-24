/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorephotonheat.c                              */
/*                                                                           */
/* Created:       2016/06/21 (VVa)                                           */
/* Last modified: 2019/01/17 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores photon energy deposition                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScorePhotonHeat:"

/*****************************************************************************/

void ScorePhotonHeat(long part, long mat, double dE0,
                     double x0, double y0, double z0, double E0, double t0,
                     double wgt, long id)
{
  long det, loc0, idx, rbin, ptr, mt, norm;
  double val;

  /* Check particle type */

  if ((long)RDB[part + PARTICLE_TYPE] != PARTICLE_TYPE_GAMMA)
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Set scored value */

   val = dE0*MEV;

  /* Score for normalization and multi-physics interface */

  if ((long)RDB[DATA_EDEP_MODE] == EDEP_MODE_NEUTRON_PHOTON)
    {
      /* Correction */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
        {
          val = val*RDB[part + PARTICLE_EDEP_RENORM];
        }

      /* Score for multi-physics interface */

      ScoreInterfacePower(val, wgt, x0, y0, z0, 0.0, part, id);

      /* Pointer to total energy deposition */

      ptr = (long)RDB[RES_TOT_FISSE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Score */

      AddBuf1D(val, wgt, ptr, id, 0);

      /* Check material burn flag */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        AddBuf1D(val, wgt, ptr, id, 1);
      else
        AddBuf1D(val, wgt, ptr, id, 2);

      /* Normalization to energy deposition in single material */

      if ((norm = (long)RDB[mat + MATERIAL_PTR_NORM]) > VALID_PTR)
        {
          ptr = (long)RDB[norm + NORM_PTR_FISSE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(val, wgt, ptr, id, 0);
        }
    }

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Reset response index */

      rbin = 0;

      /* Get pointer to response functions */

      loc0 = (long)RDB[det + DET_PTR_RBINS];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Check particle type and super-imposed */

      if ((((long)RDB[det + DET_PARTICLE] != PARTICLE_TYPE_GAMMA) &&
           ((long)RDB[loc0 + DET_RBIN_MT] != MT_MACRO_HEATTOT)) ||
          ((long)RDB[det + DET_PTR_SBINS] > VALID_PTR))
        {
          /* Next detector */

          det = NextItem(det);

          /* Cycle loop */

          continue;
        }

      /* Loop over responses */

      while (loc0 > VALID_PTR)
        {
          /* Check mt */

          mt = (long)RDB[loc0 + DET_RBIN_MT];

          if ((mt == MT_MACRO_HEATPHOTANA) || (mt == MT_MACRO_HEATTOT
              && (long)RDB[DATA_EDEP_MODE] == EDEP_MODE_NEUTRON_PHOTON))
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
