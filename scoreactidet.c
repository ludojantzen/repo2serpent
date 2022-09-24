/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoreactidet.c                                 */
/*                                                                           */
/* Created:       2017/02/08 (JLe)                                           */
/* Last modified: 2017/02/08 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Scores transmutation cross sections for activation detectors */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreActiDet:"

/*****************************************************************************/

void ScoreActiDet(double flx, long mat0, long part, double x, double y, 
                  double z, double E, double wgt, double t, long id)
{
  long det, ptr, mat;

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Check if scored */

      if (DetBin(det, mat0, part, x, y, z, E, t, id) < 0)
        {
          /* Next detector */

          det = NextItem(det);

          /* Cycle loop */

          continue;
        }

      /* Loop over activation bins */

      ptr = (long)RDB[det + DET_PTR_ABINS];
      while (ptr > VALID_PTR)
        {
          /* Get material */

          mat = (long)RDB[ptr + DET_ABIN_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Score */

          ScoreTransmuXS(flx, mat, E, wgt, id);

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next */
      
      det = NextItem(det);
    } 
}

/*****************************************************************************/
