/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : storesensevent.c                               */
/*                                                                           */
/* Created:       2017/03/22 (VVa)                                           */
/* Last modified: 2017/03/22 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Creates an event for particle for use in sensitivity calcs.  */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StoreSensEvent:"

/*****************************************************************************/

void StoreSensEvent(long part, long label, double E, double val, long id)
{
  long ptr;
  double hitmiss;

  /* Check if storing events for Sens */

  if (!((long)RDB[DATA_EVENT_RECORD_FLAGS] & RECORD_EVENT_SENS))
    return;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  if ((long)RDB[DATA_SENS_SCORE_TYPE] == SENS_SCORE_TYPE_EVENT)
    {
      /* New event from bank */

      ptr = EventFromBank(part, id);

      /* Put type */

      WDB[ptr + EVENT_LABEL] = (double)label;

      /* Put energy */

      WDB[ptr + EVENT_E] = E;

      /* Put value */

      WDB[ptr + EVENT_VAL] = val;
    }
  else if ((long)RDB[DATA_SENS_SCORE_TYPE] == SENS_SCORE_TYPE_DIRECT)
    {
      /* Add score directly to score matrix */

      /* Turn label into a nonnegative one */

      if (label < 0)
        {
          label *= -1;
          hitmiss = -1;
        }
      else
        hitmiss = 1;

      /* Get pointer to particle event block */

      ptr = (long)RDB[part + PARTICLE_PTR_SENS_EBLOCK];

      /* Link a new block if it does not exist */

      if (ptr < VALID_PTR)
        ptr = EBlockFromBank(part, id);

      /* Get pointer to score array */

      ptr = (long)RDB[ptr + SENS_EBLOCK_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr_2)", DATA_ARRAY, ptr);

      /* Put value to correct place */

      WDB[ptr + label] += hitmiss*val;
    }
}
