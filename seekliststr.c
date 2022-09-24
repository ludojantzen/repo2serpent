/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : seekliststr.c                                  */
/*                                                                           */
/* Created:       2010/09/20 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: - Seeks item with matching string                            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SeekListStr:"

/*****************************************************************************/

long SeekListStr(long ptr, long param, char *str)
{
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Move to beginning */

  ptr = FirstItem(ptr);

  /* Loop over list */
  
  while (ptr > VALID_PTR)
    {
      /* Compare */

      if (!strcmp(GetText(ptr + param), str))
	return ptr;
	      
      /* Next */
      
      ptr = NextItem(ptr);
    }

  /* Not found, return null */

  return -1;
}

/*****************************************************************************/
