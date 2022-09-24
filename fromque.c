/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fromque.c                                      */
/*                                                                           */
/* Created:       2011/03/09 (JLe)                                           */
/* Last modified: 2019/09/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Retrieves neutron / photon from que                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FromQue:"

/*****************************************************************************/

long FromQue(long id)
{
  long ptr, new;

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Get pointer */

  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_QUE, id)];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Try common que first */

  if (RDB[DATA_COMMON_QUE_LIM] != 0.0)
    if (ListSize(ptr) < 3)
      if ((new = FromCommonQue(id)) > VALID_PTR)
        return new;

  /* Get last particle */

  ptr = LastItem(ptr);
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check type */

  if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
    {
      /* Que is empty, return null */

      return -1;
    }

  /* Check multiplicity */

  if ((long)RDB[ptr + PARTICLE_MULTIPLICITY] > 0)
    {
      /* Subtract value */

      WDB[ptr + PARTICLE_MULTIPLICITY] =
        RDB[ptr + PARTICLE_MULTIPLICITY] - 1.0;

      /* Duplicate */

      ptr = DuplicateParticle(ptr, id);
      WDB[ptr + PARTICLE_MULTIPLICITY] = 0.0;
    }
  else
    {
      /* Remove particle from que */

      RemoveItem(ptr);
    }

  /* Return pointer */

  return ptr;
}

/*****************************************************************************/
