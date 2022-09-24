/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addvaluepair.c                                 */
/*                                                                           */
/* Created:       2011/07/30 (JLe)                                           */
/* Last modified: 2011/11/30 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Adds avalue in a pair of values for later use, taking        */
/*              threading into account                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddValuePair:"

/*****************************************************************************/

void AddValuePair(long ptr, double f0, double f1, long id)
{
  /* Check pointer (ei voi käyttää VALID_PTR:ää) */
  
#ifdef DEBUG
  
  if ((ptr < 1) || (ptr > (long)RDB[DATA_ALLOC_MAIN_SIZE] - 1))
    Die(FUNCTION_NAME, "Pointer error");

#endif

  /* Get pointer to data */

  ptr = (long)RDB[ptr];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  
  /* Check value */

  if (f1 == -INFTY)
    Die(FUNCTION_NAME, "Trying to store -INFTY");

  /* Put new key */

  PutPrivateData(ptr++, f0, id);

  /* Calculate new value */

  f1 = f1 + GetPrivateData(ptr, id);

  /* Put value */

  PutPrivateData(ptr, f1, id);

  /***************************************************************************/
}

/*****************************************************************************/
