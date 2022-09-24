/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setoption.c                                    */
/*                                                                           */
/* Created:       2010/11/19 (JLe)                                           */
/* Last modified: 2015/10/02 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Sets option (flag) in parameter                              */
/*                                                                           */
/* Comments: - From Serpent 1.1.12                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetOption:"

/*****************************************************************************/

void SetOption(long ptr, long opt)
{
  /* Check pointer (voi olla myös < VALID_PTR) */

  if ((ptr < 0) || (ptr > (long)RDB[DATA_ALLOC_MAIN_SIZE] - 1))
    Die(FUNCTION_NAME, "Pointer error");

  /* Set option */

  WDB[ptr] = (double)((long)RDB[ptr] | opt);
}

/*****************************************************************************/
