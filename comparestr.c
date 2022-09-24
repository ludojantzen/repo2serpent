/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : comparestr.c                                   */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2012/09/22 (JLe)                                           */
/* Version:       2.1.9                                                      */
/*                                                                           */
/* Description: Compares two strings in DATA structure                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CompareStr:"

/*****************************************************************************/

long CompareStr(long ptr1, long ptr2)
{
  /* Check pointers (cannot use VALID_PTR here, because the string may */
  /* be stored in the fixed data block) */

  if (ptr1 < 1)
    Die(FUNCTION_NAME, "Invalid pointer 1");
  else if (ptr2 < 1)
    Die(FUNCTION_NAME, "Invalid pointer 2");

  /* Compare */
  
  return !strcmp(GetText(ptr1), GetText(ptr2));
}

/*****************************************************************************/
