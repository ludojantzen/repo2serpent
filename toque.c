/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : toque.c                                        */
/*                                                                           */
/* Created:       2011/03/09 (JLe)                                           */
/* Last modified: 2019/09/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Puts neutron / photon in que                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ToQue:"

/*****************************************************************************/

void ToQue(long ptr, long id)
{
  long lim;

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check if common que is used */

  if ((lim = (long)RDB[DATA_COMMON_QUE_LIM]) != 0)
    if ((long)RDB[ptr + PARTICLE_MULTIPLICITY] > lim - 1)
      {
        /* Put particle in common que instead */

        ToCommonQue(ptr);

        /* Exit */

        return;
      }

  /* Add item in list */

  AddItem(OMPPtr(DATA_PART_PTR_QUE, id), ptr);

  /* Sort list to transport neutrons before photons */

  if (((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES) &&
      ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES))
    {
      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_QUE, id)];
      SortList(ptr, PARTICLE_TYPE, SORT_MODE_ASCEND);
    }
}

/*****************************************************************************/
