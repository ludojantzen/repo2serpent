/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : tostore.c                                      */
/*                                                                           */
/* Created:       2015/01/19 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Puts particle in store                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ToStore:"

/*****************************************************************************/

void ToStore(long ptr, long id, long eoi)
{
  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Add item in list */

  if(eoi)
    AddItem(OMPPtr(DATA_PART_PTR_EOI_STORE, id), ptr);
  else
    AddItem(OMPPtr(DATA_PART_PTR_BOI_STORE, id), ptr);
}

/*****************************************************************************/
