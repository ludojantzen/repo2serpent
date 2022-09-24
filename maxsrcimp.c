/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : maxsrcimp.c                                    */
/*                                                                           */
/* Created:       2017/01/28 (JLe)                                           */
/* Last modified: 2018/09/25 (Jle)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Finds maximum source importance by random sampling           */
/*                                                                           */
/* Comments: - NOTE: tää on turha ja pitäis poistaa?                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MaxSrcImp:"
#define MAX_SRC_RESAMPLE 1000000

/*****************************************************************************/

void MaxSrcImp(long type, double x,  double y,  double z,  double E, long id)
{
  long wwd, msh, i, n, ptr, erg, ne;
  double p; 

  /* Check calculation */

  if ((long)RDB[DATA_SRC_IMP_CALC] == NO)
    return;

  /* Get importance */

  p = WWImportance(type, x, y, z, 0.0, 0.0, 0.0, E, WWMESH_SRC);

  /* Pointer to weight window structure */
  
  if ((wwd = (long)RDB[DATA_PTR_WWD0]) < VALID_PTR)
    return;

  /* Loop over weight window structures */

  while (wwd > VALID_PTR)
    {
      /* Pointer to mesh */
      
      if ((msh = (long)RDB[wwd + WWD_PTR_MESH]) < VALID_PTR)
        {
          /* Pointer to next */

          wwd = NextItem(wwd);

          /* Cycle loop */

          continue;
        }

      /* Check particle type */

      if (((long)RDB[wwd + WWD_PARTICLE_TYPE] > 0) &&
          ((long)RDB[wwd + WWD_PARTICLE_TYPE] != type))
        {
          /* Pointer to next */

          wwd = NextItem(wwd);

          /* Cycle loop */

          continue;
        }
      
      /* Get mesh index */

      if ((i = MeshIndex(msh, x, y, z, -1.0)) > -1)
        {
          /* Pointer to energy array */
          
          if ((erg = (long)RDB[wwd + WWD_PTR_ERG]) < VALID_PTR)
            n = 0;
          else
            {
              /* Number of energy groups */

              ne = (long)RDB[wwd + WWD_NE];
              CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 10000);

              /* Get group */

              n = SearchArray(&RDB[erg], E, ne + 1);
            }

          /* Score */

          ptr = (long)RDB[wwd + WWD_PTR_SRC_IMP_DIS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(p, 1.0, ptr, id, -1, n, i);
        }

      /* Next */

      wwd = NextItem(wwd);
    }
}

/*****************************************************************************/
