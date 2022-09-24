/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : srcdet.c                                       */
/*                                                                           */
/* Created:       2015/06/27 (JLe)                                           */
/* Last modified: 2016/08/01 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Scores source point distribution for detectors               */
/*                                                                           */
/* Comments: - Separated from coldet.c 15/06/27 (2.1.25)                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SrcDet:"

/*****************************************************************************/

void SrcDet(long part, long mat, double x0, double y0, double z0, double u0, 
            double v0, double w0, double E0, double t0, double wgt, long id)
{
  long det, loc0, idx, rbin, ptr, type;
  double val, u, v, w;

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Compare particle types and check super-imposed */
      
      if ((type != (long)RDB[det + DET_PARTICLE]) || 
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

          if ((long)RDB[loc0 + DET_RBIN_MT] == MT_SOURCE_RATE)
            {
              /* Get bin index */

              if ((idx = DetBin(det, mat, part, x0, y0, z0, E0, t0, id)) < 0)
                break;

              /* Set scored value */

              val = 1.0;

              /* Get direction vector (onkohan tää ihan turha?) */

              u = RDB[det + DET_DIRVEC_U];
              v = RDB[det + DET_DIRVEC_V];
              w = RDB[det + DET_DIRVEC_W];
              
              /* Check if non-zero */
              
              if ((u != 0.0) || (v != 0.0) || (w != 0.0))
                {
                  /* Calculate scalar product */
                  
                  val = val*(u*u0 + v*v0 + w*w0);
                  
                  /* Check negative */
                  
                  if (val < 0.0)
                    val = 0.0;
                }
              
              /* Get pointer to statistics */
              
              ptr = (long)RDB[det + DET_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              
              /* Score */

              AddBuf(val, wgt, ptr, id, -1, idx, rbin);
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
