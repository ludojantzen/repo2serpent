/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : eventfrombank.c                                */
/*                                                                           */
/* Created:       2014/09/29 (JLe)                                           */
/* Last modified: 2017/06/02 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Retrieves a new event structure from bank                    */
/*                                                                           */
/* Comments: - Structure should be cleared when returned                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "EventFromBank:"

/*****************************************************************************/

long EventFromBank(long part, long id)
{
  long ptr, next;

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Get pointer to bank */

  ptr = (long)RDB[OMPPtr(DATA_PTR_EVENT_BANK, id)];

  if (ptr  < VALID_PTR)
    Die(FUNCTION_NAME, "Event bank is empty");

  /* Pointer to next */

  next = NextItem(ptr);

  /* Set pointer */

  WDB[OMPPtr(DATA_PTR_EVENT_BANK, id)] = (double)next;

  /* Attach pointer to particle */

  WDB[ptr + LIFO_LIST_PTR_NEXT] = RDB[part + PARTICLE_PTR_EVENTS];
  WDB[part + PARTICLE_PTR_EVENTS] = (double)ptr;

  /* Set counter */

  WDB[ptr + EVENT_HIS_COUNT] = 1.0;

  /* Return pointer */

  return ptr;
}

/*****************************************************************************/
