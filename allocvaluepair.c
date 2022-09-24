/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : allocvaluepair.c                               */
/*                                                                           */
/* Created:       2011/01/05 (JLe)                                           */
/* Last modified: 2019/11/10 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Allocates memory for pair of values that can be used to      */
/*              store temporary data                                         */
/*                                                                           */
/* Comments: - Used for cross sections, energy grid indexes, etc.            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AllocValuePair:"

/*****************************************************************************/

void AllocValuePair(long loc0)
{
  long ptr, loc1;

  /* Check pointer (ei voi käyttää VALID_PTR:ää) */

  if ((loc0 < 1) || (loc0 > (long)RDB[DATA_ALLOC_MAIN_SIZE] - 1))
    Die(FUNCTION_NAME, "Pointer error");

  /* Check that memory is not already allocated */

  if ((long)RDB[loc0] > VALID_PTR)
    Die(FUNCTION_NAME, "Memory already allocated");

  /* Allocate memory from private array */

  ptr = AllocPrivateData(2, PRIVA_ARRAY);

  /* Put pointer */

  WDB[loc0] = (double)ptr;

  /* Add pointer in list */

  loc1 = NewItem(DATA_PTR_VP0, VP_BLOCK_SIZE);
  WDB[loc1 + VP_PTR] = (double)ptr;

  /***************************************************************************/
}

/*****************************************************************************/
