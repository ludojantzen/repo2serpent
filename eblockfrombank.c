/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : eblockfrombank.c                               */
/*                                                                           */
/* Created:       2017/06/02 (VVa)                                           */
/* Last modified: 2017/06/02 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Retrieves a new event block structure from bank              */
/*                                                                           */
/* Comments: - Structure should be cleared when returned                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "EBlockFromBank:"

/*****************************************************************************/

long EBlockFromBank(long part, long id)
{
  long ptr, next;

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Get pointer to bank */

  ptr = (long)RDB[OMPPtr(DATA_PTR_EBLOCK_BANK, id)];

  if (ptr  < VALID_PTR)
    Die(FUNCTION_NAME, "Event score matrix bank is empty. "
        "Increase the bank size with \"sens opt direct <sz>\".");

  /* Pointer to next */

  next = NextItem(ptr);

  /* Set pointer */

  WDB[OMPPtr(DATA_PTR_EBLOCK_BANK, id)] = (double)next;

  /* Attach pointer to particle */

  WDB[ptr + LIFO_LIST_PTR_NEXT] = RDB[part + PARTICLE_PTR_SENS_EBLOCK];
  WDB[part + PARTICLE_PTR_SENS_EBLOCK] = (double)ptr;

  /* Set counter */

  WDB[ptr + SENS_EBLOCK_HIS_COUNT] = 1.0;

  /* Return pointer */

  return ptr;
}

/*****************************************************************************/
