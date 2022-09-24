/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorescaresp.c                                 */
/*                                                                           */
/* Created:       2016/04/24 (JLe)                                           */
/* Last modified: 2017/10/08 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Scores response term for source convergence acceleration     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreSCAResp:"

/*****************************************************************************/

void ScoreSCAResp(long part, double res, double x, double y, double z, long id)
{
  long rmx, loc0, det, ptr, i, n, msh, nx, ny, nz;
 
  /* Check if source convergence acceleraion is on */

  if ((long)RDB[DATA_RMX_CONVG_ACC] == NO)
    return;

  /* Check if response matrix calculation is on */

  if ((long)RDB[DATA_RMTX_CALC] == NO)
    return;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Check particle type */

  if ((long)RDB[part + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
    return;

  /* Pointer to weight window structure */
  
  rmx = (long)RDB[DATA_PTR_RMX0];
  CheckPointer(FUNCTION_NAME, "(rmx)", DATA_ARRAY, rmx);

  /* Does not work with multiple structures */

  if (NextItem(rmx) > VALID_PTR)
    Die(FUNCTION_NAME, "Multiple structures");

  /* Get pointer */

  loc0 = (long)RDB[part + PARTICLE_ICM_PTR_ICM];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Pointer to sub-mesh */

  msh = (long)RDB[loc0 + RMX_CELL_PTR_SUBMESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get index */

  if ((n = MeshIndex(msh, x, y, z, 0.0)) < 0)
    {
      /* Print warning */

      Warn(FUNCTION_NAME, "Error in sub-index");

      /* Exit */

      return;
    }

  /* Local mesh size */

  nx = (long)RDB[rmx + RMX_CONVG_SUBMESH_NX];
  ny = (long)RDB[rmx + RMX_CONVG_SUBMESH_NY];
  nz = (long)RDB[rmx + RMX_CONVG_SUBMESH_NZ];

  /* Check */

  CheckValue(FUNCTION_NAME, "n", "", n, 0, nx*ny*nz - 1);

  /* Pointer to detector */

  det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
  CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

  /* Get index */

  if ((i = (long)RDB[part + PARTICLE_ICM_IDX]) == -1)
    {
      /* Particle came from source, score direct contributions */
      
      ptr = (long)RDB[det + RMX_DET_MC_RES_SRCC];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, res, id);

      /* Score total contribution */

      ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, res, id);

      /* Add to number of scores */

      ptr = (long)RDB[det + RMX_DET_MC_RES_SCORE_N];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, 1.0, id);

      /* Score form factor */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_SRC_FF];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr + n, res, id);
    }
  else if (i > -1)
    {
      /* Particle came from boundary */

      ptr = (long)RDB[det + RMX_DET_MC_RES_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr + i, res, id);

      /* Score total contribution */

      ptr = (long)RDB[det + RMX_DET_MC_RES_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, res, id);

      /* Add to number of scores */

      ptr = (long)RDB[det + RMX_DET_MC_RES_SCORE_N];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr, 1.0, id);

      /* Score form factor */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_SURF_FF];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      AddPrivateRes(ptr + i*nx*ny*nz + n, res, id);
    }
}

/*****************************************************************************/
