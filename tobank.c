/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : tobank.c                                       */
/*                                                                           */
/* Created:       2012/10/12 (JLe)                                           */
/* Last modified: 2012/10/12 (JLe)                                           */
/* Version:       2.1.9                                                      */
/*                                                                           */
/* Description: Puts particle in bank                                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ToBank:"

/*****************************************************************************/

void ToBank(long ptr, long id)
{
  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Add item in list */

  AddItem(OMPPtr(DATA_PART_PTR_BANK, id), ptr);
}

/*****************************************************************************/
