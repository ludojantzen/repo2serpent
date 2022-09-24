/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : eblocktobank.c                                 */
/*                                                                           */
/* Created:       2017/06/02 (VVa)                                           */
/* Last modified: 2018/06/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns event block to bank                                  */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "EBlockToBank:"

/*****************************************************************************/

void EBlockToBank(long ptr, long id)
{
  long next, loc0, ptr2, maxi;

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Reset history count */

  WDB[ptr + SENS_EBLOCK_HIS_COUNT] = 0.0;

  /* Get pointer to sensitivity block */

  loc0 = (long)RDB[DATA_PTR_SENS0];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Reset data */

  maxi = (long)RDB[loc0 + SENS_MAX_LABEL];

  ptr2 = (long)RDB[ptr + SENS_EBLOCK_PTR_DATA];
  memset(&WDB[ptr2], 0.0, maxi*sizeof(double));

  /* Pointer to next */

  next = (long)RDB[OMPPtr(DATA_PTR_EBLOCK_BANK, id)];

  /* Put pointers */

  WDB[OMPPtr(DATA_PTR_EBLOCK_BANK, id)] = (double)ptr;
  WDB[ptr + LIFO_LIST_PTR_NEXT] = (double)next;
}

/*****************************************************************************/
