/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : eventtobank.c                                  */
/*                                                                           */
/* Created:       2014/09/29 (JLe)                                           */
/* Last modified: 2017/06/02 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Returns event to bank                                        */
/*                                                                           */
/* Comments: - NOTE: Tätä pitää kutsua OpenMP barrierin takaa                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "EventToBank:"

/*****************************************************************************/

void EventToBank(long ptr, long id)
{
  long next;

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Reset data */

  memset(&WDB[ptr + LIFO_LIST_DATA_SIZE], 0.0,
         (EVENT_BLOCK_SIZE - LIFO_LIST_DATA_SIZE)*sizeof(double));

  /* Pointer to next */

  next = (long)RDB[OMPPtr(DATA_PTR_EVENT_BANK, id)];

  /* Put pointers */

  WDB[OMPPtr(DATA_PTR_EVENT_BANK, id)] = (double)ptr;
  WDB[ptr + LIFO_LIST_PTR_NEXT] = (double)next;
}

/*****************************************************************************/
