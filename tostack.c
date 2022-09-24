/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : tostack.c                                      */
/*                                                                           */
/* Created:       2011/03/09 (JLe)                                           */
/* Last modified: 2016/02/01 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Stores neutron / photon / precursor to stack                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ToStack:"

/*****************************************************************************/

void ToStack(long ptr, long id)
{
  long loc0, type;

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get particle type */

  type = (long)RDB[ptr + PARTICLE_TYPE];

  /* Avoid compiler warning */

  loc0 = -1;

  /* Get stack pointer */

  if ((long)RDB[DATA_STOP_AFTER_PLOT] == STOP_AFTER_PLOT_TRACKS)
    loc0 = OMPPtr(DATA_PART_PTR_TRK_BANK, id);
  else if (type == PARTICLE_TYPE_NEUTRON)
    loc0 = OMPPtr(DATA_PART_PTR_NSTACK, id);
  else if (type == PARTICLE_TYPE_GAMMA)
    loc0 = OMPPtr(DATA_PART_PTR_GSTACK, id);
  else if (type == PARTICLE_TYPE_PRECURSOR)
    loc0 = OMPPtr(DATA_PART_PTR_PSTACK, id);
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Put item to stack */

  AddItem(loc0, ptr);
}

/*****************************************************************************/
