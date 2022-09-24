/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fromstore.c                                    */
/*                                                                           */
/* Created:       2015/01/19 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Retrieves particle from store                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FromStore:"

/*****************************************************************************/

long FromStore(long id, long eoi)
{
  long ptr;

  /* Check actual thread number */

  if (OMP_THREAD_NUM != 0)
    Die(FUNCTION_NAME, "Called from an OpenMP parallel loop");

  /* Get pointer */

  if(eoi)
    ptr = (long)RDB[OMPPtr(DATA_PART_PTR_EOI_STORE, id)];
  else
    ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BOI_STORE, id)];

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

  /* Remove particle from store */

  RemoveItem(ptr);

  /* Return pointer */ 
  
  return ptr;
}

/*****************************************************************************/
