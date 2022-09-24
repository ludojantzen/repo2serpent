/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processevents.c                                */
/*                                                                           */
/* Created:       2011/10/01 (JLe)                                           */
/* Last modified: 2018/01/26 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Removes unused events from killed particles, allocates more  */
/*              memory for bank, etc.                                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessEvents:"

/*****************************************************************************/

void ProcessEvents()
{
  long ptr, part, loc0, prev, n, np, id, bankslow, cutoff;
  long id0, id1, min, max, loc1, sz, next, tot, ave;
#ifdef DEBUG
  long sz0, count;
#endif

  /* Avoid compiler warning without debug mode */

  part = 0;
  id = 0;

  /* Check if events are set */

  if ((long)RDB[DATA_EVENT_RECORD_FLAGS] == 0)
    return;

#ifdef DEBUG

  /***************************************************************************/

  /***** Check history counts ************************************************/

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */

  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Loop over particles */

  while (part > VALID_PTR)
    {
      /* Reset previous */

      n = 1;

      /* Loop over events */

      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
        {
          /* Check */

          if ((long)RDB[ptr + EVENT_HIS_COUNT] < n)
            Die(FUNCTION_NAME, "Something wrong here");
          else
            n = (long)RDB[ptr + EVENT_HIS_COUNT];

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next particle */

      part = NextItem(part);
    }

  /***************************************************************************/

  /***** Confirm count *******************************************************/

  /* Reset counter */

  count = 0;

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */

  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Loop over particles */

  while (part > VALID_PTR)
    {
      /* Loop over events */

      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
        {
          /* Check if event was already counted */

          if ((long)RDB[ptr + EVENT_HIS_COUNT] > 0)
            {
              /* Add to count */

              count++;

              /* Mark as counted */

              WDB[ptr + EVENT_HIS_COUNT] =
                -RDB[ptr + EVENT_HIS_COUNT];
            }
          else if ((long)RDB[ptr + EVENT_HIS_COUNT] == 0)
            Die(FUNCTION_NAME, "WTF?");

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next particle */

      part = NextItem(part);
    }

  /* Count banks for all threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    count = count + LIFOListSize(OMPPtr(DATA_PTR_EVENT_BANK, id), -1);

  /* Check that all events are accounted for */

  if (count != (long)RDB[DATA_EVENT_BANK_SZ])
    Die(FUNCTION_NAME, "Mismatch in number of events");

  /***************************************************************************/

  /***** Set positive values back to counters ********************************/

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */

  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Loop over particles */

  while (part > VALID_PTR)
    {
      /* Loop over events */

      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
        {
          /* Check if negative and switch to positive */

          if ((long)RDB[ptr + EVENT_HIS_COUNT] < 0)
            WDB[ptr + EVENT_HIS_COUNT] = -RDB[ptr + EVENT_HIS_COUNT];

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next particle */

      part = NextItem(part);
    }

#endif

  /***************************************************************************/

  /***** Remove old events ***************************************************/

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */

  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Reset id */

  id = 0;

  /* Loop over particles */

  while (part > VALID_PTR)
    {
      /* Reset count */

      n = 0;

      /* Loop over events and update history counters */

      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
        {
          /* Check if fission and update counter */

          if ((long)RDB[ptr + EVENT_TYPE] == EVENT_TYPE_FISS)
            n++;

          /* Check with limit */

          if (n > (long)RDB[DATA_EVENT_MAX_GEN] - 1)
            {
              /* Update history count */

              WDB[ptr + EVENT_HIS_COUNT] = RDB[ptr + EVENT_HIS_COUNT] - 1.0;
            }

          /* Check history counter */

          if ((long)RDB[ptr + EVENT_HIS_COUNT] == 0)
            {
              /* Remove item */

              loc0 = ptr;

              /* Pointer to next */

              ptr = NextItem(loc0);

              /* Back to bank */

              EventToBank(loc0, id++);

              /* Check next id */

              if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
                id = 0;

              /* Mark as removed */

              WDB[loc0 + EVENT_HIS_COUNT] = -2E+6;
            }
          else if ((long)RDB[ptr + EVENT_HIS_COUNT] > 0)
            {
              /* Pointer to next */

              ptr = NextItem(ptr);
            }
          else if ((long)RDB[ptr + EVENT_HIS_COUNT] == -2E+6)
            break;
          else
            Die(FUNCTION_NAME, "WTF, negative? %ld %ld",
                (long)RDB[ptr + EVENT_HIS_COUNT],
                (long)RDB[ptr + EVENT_TYPE]);

        }

      /* Next particle */

      part = NextItem(part);
    }

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */

  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Loop over particles */

  while (part > VALID_PTR)
    {
       /* Loop over events and detach pointers */

      prev = -1;

      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
        {
          /* Check if event was removed */

          if ((long)RDB[ptr + EVENT_HIS_COUNT] == -2E+6)
            {
              /* Detach previous */

              if (prev < VALID_PTR)
                WDB[part + PARTICLE_PTR_EVENTS] = NULLPTR;
              else
                WDB[prev + LIFO_LIST_PTR_NEXT] = NULLPTR;

              /* Break loop */

              break;
            }

          /* Copy pointer */

          prev = ptr;

          /* Pointer to next */

          ptr = NextItem(ptr);
        }

#ifdef DEBUG
      /* Check history counts */

      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
        {
          if ((long)RDB[ptr + EVENT_HIS_COUNT] < 1)
            Die(FUNCTION_NAME, "Error in history count");

          /* Pointer to next */

          ptr = NextItem(ptr);
        }

#endif

      /* Next particle */

      part = NextItem(part);
    }

  /***************************************************************************/

  /***** Allocate memory for more events if needed ***************************/

  bankslow = 0;

  /* Calculate initial event bank size per thread */

  sz = (long)(RDB[DATA_EVENT_BANK_SZ]/RDB[DATA_OMP_MAX_THREADS]);

  /* Calculate threshold for low banks */

  cutoff = 0.2*sz;

  /* Pointer to event bank */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      if ((ptr = (long)RDB[OMPPtr(DATA_PTR_EVENT_BANK,id)]) > VALID_PTR)
        {
          /* Check size */

          if (LIFOListSize(OMPPtr(DATA_PTR_EVENT_BANK,id), cutoff) < cutoff)
            {
              bankslow = 1;
              break;
            }
        }
    }

  /***************************************************************************/

  /* If one or more banks are low, first redistribute them */

  if (bankslow)
    {
      /***** Redistribute events *******************************************/

      loc0 = DATA_PTR_EVENT_BANK;

      /* Check pointer */

      if ((long)RDB[loc0] < VALID_PTR)
        Die(FUNCTION_NAME, "Pointer error");

#ifdef DEBUG

      sz0 = 0;

      /* Count banks first */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {
          /* Add bank size */

          ptr = (long)RDB[OMPPtr(loc0, id)];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          sz0 = sz0 + LIFOListSize(ptr, -1);
        }
#endif


      /* Loop */

      while (1 != 2)
        {
          /* Find bank with lowest and highest number of particles */

          id0 = -1;
          id1 = -1;

          min = 100000000000;
          max = -1;

          tot = 0;

          for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
            {
              /* Get list size */

              loc1 = (long)RDB[OMPPtr(loc0, id)];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
              sz = LIFOListSize(loc1, -1);

              /* Add to total */
              tot += sz;

              /* Compare to minimum */

              if (sz < min)
                {
                  min = sz;
                  id0 = id;
                }

              /* Compare to maximum */

              if (sz > max)
                {
                  max = sz;
                  id1 = id;
                }
            }

          /* Get average number of events per thread */

          ave = (long)(tot/RDB[DATA_OMP_MAX_THREADS]);

          /* Check condition */

          if (max < 2*min)
            break;
          else
            {
              /* Count the number of events to move as the minimum required */
              /* to get one of the banks to average */

              sz = max - ave;
              if (ave - min < sz)
                sz = ave - min;

              /* Move events to min bank */

              for (n = 0; n < sz; n++)
                {
                  /* Get event from bank */

                  ptr = (long)RDB[OMPPtr(loc0, id1)];
                  next = NextItem(ptr);
                  WDB[OMPPtr(loc0, id1)] = (double)next;

                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Put event to other bank */

                  next = (long)RDB[OMPPtr(loc0, id0)];
                  WDB[OMPPtr(loc0, id0)] = (double)ptr;
                  WDB[ptr + LIFO_LIST_PTR_NEXT] = (double)next;
                }
            }

        }

      /* Check if the banks are still low */

      if (min > cutoff)
        bankslow = 0;

#ifdef DEBUG
      /* Check that we did not lose any events */
      tot = 0;
      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {
          /* Add bank size */

          ptr = (long)RDB[OMPPtr(loc0, id)];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          tot = tot + LIFOListSize(ptr, -1);
        }

      if (tot != sz0)
        Die(FUNCTION_NAME, "Lost %ld events", tot - sz0);
#endif
    }

  /* If one or more banks are still low, allocate more particles */

  if (bankslow)
    {
      /* Pointer to event bank */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {
          if ((ptr = (long)RDB[OMPPtr(DATA_PTR_EVENT_BANK,id)]) > VALID_PTR)
            {
              /* Check size */

              if (LIFOListSize(OMPPtr(DATA_PTR_EVENT_BANK,id), cutoff) < cutoff)
                {
                  /* Print */

                  Note(0, "Adjusting event bank size for thread %ld...", id);

                  /* Set expansion size */

                  np = cutoff;

                  /* Allow memory allocation */

                  Mem(MEM_ALLOW);

                  /* Allocate memory for events */

                  for (n = 0; n < np; n++)
                    NewLIFOItem(OMPPtr(DATA_PTR_EVENT_BANK, id), EVENT_BLOCK_SIZE);

                  /* Deny memory allocation */

                  Mem(MEM_DENY);

                  /* Put new size */

                  WDB[DATA_EVENT_BANK_SZ] = RDB[DATA_EVENT_BANK_SZ] + (double)np;
                }
            }
        }
    }

  /*****************************************************************************/
}
