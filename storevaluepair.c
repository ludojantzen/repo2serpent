/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : storevaluepair.c                               */
/*                                                                           */
/* Created:       2011/01/05 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Stores a pair of values for later use, taking threading into */
/*              account                                                      */
/*                                                                           */
/* Comments: - Used for cross sections, energy grid indexes, etc.            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StoreValuePair:"

/*****************************************************************************/

void StoreValuePair(long ptr, double f0, double f1, long id)
{
  
#ifdef DEBUG

  /* Check pointer (ei voi käyttää VALID_PTR:ää) */
  
  if ((ptr < 1) || (ptr > (long)RDB[DATA_ALLOC_MAIN_SIZE] - 1))
    Die(FUNCTION_NAME, "Pointer error");

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "PRIVA array not ready for access");

#endif

  /* Get pointer to data */

  ptr = (long)RDB[ptr];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  /* Store pair */

  PutPrivateData(ptr++, f0, id);
  PutPrivateData(ptr++, f1, id);

  /***************************************************************************/
}

/*****************************************************************************/
