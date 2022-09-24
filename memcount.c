/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : memcount.c                                     */
/*                                                                           */
/* Created:       2012/05/29 (JLe)                                           */
/* Last modified: 2012/05/29 (JLe)                                           */
/* Version:       2.1.6                                                      */
/*                                                                           */
/* Description: Updates a counter on the number of bytes allocated since     */
/*              the last call.                                               */
/*                                                                           */
/* Comments: - Used for better understanding where allocated memory is       */
/*             going.                                                        */
/*           - Should be called from main() only to avoid mixing the         */
/*             counter.                                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MemCount:"

/*****************************************************************************/

long MemCount()
{
  long sz;

  /* Calculate total number of bytes */

  CalculateBytes();

  /* Calculate allocated bytes since last call */

  sz = (long)RDB[DATA_TOTAL_BYTES] - (long)RDB[DATA_BYTE_COUNT];

  /* Update byte count */

  WDB[DATA_BYTE_COUNT] = RDB[DATA_TOTAL_BYTES];

  /* Return block size */
  
  return sz;  
}

/*****************************************************************************/
