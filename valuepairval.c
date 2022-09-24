/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : valuepairval.c                                 */
/*                                                                           */
/* Created:       2015/06/26 (JLe)                                           */
/* Last modified: 2015/06/26 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Returns stored value                                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ValuePairVal:"

/*****************************************************************************/

double ValuePairVal(long ptr, long id)
{
  /* Check pointer (ei voi käyttää VALID_PTR:ää) */
  
#ifdef DEBUG

  if ((ptr < 1) || (ptr > (long)RDB[DATA_ALLOC_MAIN_SIZE] - 1))
    Die(FUNCTION_NAME, "Pointer error");

#endif

  /* Get pointer to data */

  ptr = (long)RDB[ptr];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  /* Return value */

  return GetPrivateData(ptr + 1, id);

  /***************************************************************************/
}

/*****************************************************************************/
