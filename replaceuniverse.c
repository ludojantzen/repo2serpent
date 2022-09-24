/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : replaceuniverse.c                              */
/*                                                                           */
/* Created:       2012/07/14 (JLe)                                           */
/* Last modified: 2012/07/14 (JLe)                                           */
/* Version:       2.1.8                                                      */
/*                                                                           */
/* Description: Replaces univere 1 with universe 2                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SwapUniverses:"

/*****************************************************************************/

void ReplaceItem(long uni1, long uni2)
{
  long n;
  double tmp;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(uni1)", DATA_ARRAY, uni1);
  CheckPointer(FUNCTION_NAME, "(uni2)", DATA_ARRAY, uni2);

  /* Loop over structure */
  
  for (n = LIST_DATA_SIZE; n < UNIVERSE_BLOCK_SIZE; n++)
    {
      /* Swap contents */

      tmp = RDB[uni1 + n];
      WDB[uni1 + n] = RDB[uni2 + n];
      WDB[uni2 + n] = tmp;
    }
}

/*****************************************************************************/
