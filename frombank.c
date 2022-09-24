/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : frombank.c                                     */
/*                                                                           */
/* Created:       2012/10/12 (JLe)                                           */
/* Last modified: 2012/10/12 (JLe)                                           */
/* Version:       2.1.9                                                      */
/*                                                                           */
/* Description: Retrieves particle from bank                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FromBank:"

/*****************************************************************************/

long FromBank(long id)
{
  long ptr;

  /* Check actual thread number */

  if (OMP_THREAD_NUM != 0)
    Die(FUNCTION_NAME, "Called from an OpenMP parallel loop");

  /* Get pointer */

  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get pointer to first item */

  ptr = FirstItem(ptr);
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check type */

  if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_DUMMY)
    Die(FUNCTION_NAME, "Error in list");

  /* Get pointer to next */

  if ((ptr = NextItem(ptr)) < VALID_PTR)
    {
      /* Bank is empty, return null */

      return -1;
    }

  /* Remove particle from bank */

  RemoveItem(ptr);

  /* Return pointer */ 
  
  return ptr;
}

/*****************************************************************************/
