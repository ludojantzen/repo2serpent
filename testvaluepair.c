/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testvaluepair.c                                */
/*                                                                           */
/* Created:       2011/01/05 (JLe)                                           */
/* Last modified: 2011/11/30 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns stored value if key matches or -INFTY if not         */
/*                                                                           */
/* Comments: - Used for cross sections, energy grid indexes, etc.            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TestValuePair:"

/*****************************************************************************/

double TestValuePair(long ptr, double f0, long id)
{
  /* Check pointer (ei voi käyttää VALID_PTR:ää) */
  
#ifdef DEBUG

  if ((ptr < 1) || (ptr > (long)RDB[DATA_ALLOC_MAIN_SIZE] - 1))
    Die(FUNCTION_NAME, "Pointer error");

#endif

  /* Get pointer to data */

  ptr = (long)RDB[ptr];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  /* Compare key */

  if (f0 == GetPrivateData(ptr++, id))
    return GetPrivateData(ptr, id);
  else
    return -INFTY;

  /***************************************************************************/
}

/*****************************************************************************/
