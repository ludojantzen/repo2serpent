/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fromtrkbank.c                                  */
/*                                                                           */
/* Created:       2014/02/12 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Retrieves neutron / photon from track plotter bank           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FromTrkBank:"

/*****************************************************************************/

long FromTrkBank(long id)
{
  long ptr;

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Avoid compiler warning */

  ptr = -1;

  /* Get pointer to bank */

  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_TRK_BANK, id)];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get pointer to last item */

  ptr = LastItem(ptr);
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check type */
  
  if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
    return -1;

  /* Remove particle from bank */
  
  RemoveItem(ptr);
  
  /* Return pointer */
  
  return ptr;
}

/*****************************************************************************/
