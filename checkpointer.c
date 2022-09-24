/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkpointer.c                                 */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2019/12/04 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Checks that pointer is within allowed limits                 */
/*                                                                           */
/* Comments: - Adopted from Serpent 1.1.12                                   */
/*                                                                           */
/*           - Alarajalla oletetaan minimiksi DATA arrayn fixed block size.  */
/*             Vastaava muisti pit‰‰ varata muita blokkeja varattaessa.      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckPointer:"

/*****************************************************************************/

#ifdef DEBUG

void CheckPointer (char *fname, char *comment, long type, long ptr)
{
  long max, min;

  /* Avoid compiler warning */

  min = -1;
  max = -1;

  /* Get minimum and maximum */

  if (abs(type) == DATA_ARRAY)
    {
      min = VALID_PTR;
      max = (long)RDB[DATA_ALLOC_MAIN_SIZE] - 1;
    }
  else if (abs(type) == RES1_ARRAY)
    {
      min = VALID_PTR;
      max = (long)RDB[DATA_ALLOC_RES1_SIZE] - STAT_BLOCK_SIZE;
    }
  else if (abs(type) == ACE_ARRAY)
    {
      min = 0;
      max = (long)RDB[DATA_ALLOC_ACE_SIZE] - 1;
    }
  else if (abs(type) == PRIVA_ARRAY)
    {
      min = VALID_PTR;
      max = (long)RDB[DATA_ALLOC_PRIVA_SIZE] - 1;
    }
  else if (abs(type) == RES2_ARRAY)
    {
      min = VALID_PTR;
      max = (long)RDB[DATA_ALLOC_RES2_SIZE] - 1;
    }
  else if (abs(type) == RES3_ARRAY)
    {
      min = 0;
      max = (long)RDB[DATA_ALLOC_RES3_SIZE] - 1;
    }
  else if (abs(type) == BUF_ARRAY)
    {
      min = VALID_PTR;
      max = (long)RDB[DATA_ALLOC_BUF_SIZE] - 1;
    }
  else
    Die(fname, "Invalid data type %ld %s", type, comment);

  /* Check pointer */

  if (ptr < 0)
    Die(fname, "Null pointer %s", comment);
  else if ((type < 0) && (ptr != max))
    Die(fname, "Pointer error: value %ld not equal to maximum %s", ptr,
        comment);
  else if (ptr < min)
    Die(fname, "Pointer error: value %ld below minimum = %ld %s", ptr, min,
        comment);
  else if (ptr > max)
    Die(fname, "Pointer error: value %ld above maximum = %ld %s %ld", ptr, max,
        comment, type);
}

#endif

/*****************************************************************************/
