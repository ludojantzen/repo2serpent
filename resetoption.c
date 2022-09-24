/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resetoption.c                                  */
/*                                                                           */
/* Created:       2010/11/22 (JLe)                                           */
/* Last modified: 2017/05/04 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Resets option (flag) in parameter                            */
/*                                                                           */
/* Comments: - From Serpent 1.1.12                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ResetOption:"

/*****************************************************************************/

void ResetOption(long ptr, long opt)
{
  /* Check pointer (voi olla myös < VALID_PTR) */

  if ((ptr < 0) || (ptr > (long)RDB[DATA_ALLOC_MAIN_SIZE] - 1))
    Die(FUNCTION_NAME, "Pointer error");

  /* Set option */

  WDB[ptr] = (double)((long)RDB[ptr] & ~opt);
}

/*****************************************************************************/
