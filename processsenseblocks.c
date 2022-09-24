/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processsenseblocks.c                           */
/*                                                                           */
/* Created:       2017/06/02 (VVa)                                           */
/* Last modified: 2018/06/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Removes unused event blocks from killed particles            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessSensEBlocks:"

/*****************************************************************************/

void ProcessSensEBlocks()
{
  long ptr, part, loc0, prev, sens, n, np, id, maxi, bankslow, cutoff;
  long id0, id1, min, max, loc1, sz, next, maxt, tot, ave;
#ifdef DEBUG
  long sz0, count;
#endif

  /* Avoid compiler warning without debug mode */

  part = 0;
  id = 0;

  /* Check if recording events for sensitivity calculations */

  if (!((long)RDB[DATA_EVENT_RECORD_FLAGS] & RECORD_EVENT_SENS))
    return;

  /* Check if recording events using event blocks */

  if ((long)RDB[DATA_SENS_SCORE_TYPE] != SENS_SCORE_TYPE_DIRECT)
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

      /* Get the first (0-gen) event block */

      ptr = (long)RDB[part + PARTICLE_PTR_SENS_EBLOCK];

      while (ptr > VALID_PTR)
        {
          /* Check */

          if ((long)RDB[ptr + SENS_EBLOCK_HIS_COUNT] < n)
            Die(FUNCTION_NAME, "Something wrong here");
          else
            n = (long)RDB[ptr + SENS_EBLOCK_HIS_COUNT];

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
      /* Get the first (0-gen) event block */

      ptr = (long)RDB[part + PARTICLE_PTR_SENS_EBLOCK];

      while (ptr > VALID_PTR)
        {
          /* Check if event block was already counted */

          if ((long)RDB[ptr + SENS_EBLOCK_HIS_COUNT] > 0)
            {
              /* Add to count */

              count++;

              /* Mark as counted */

              WDB[ptr + SENS_EBLOCK_HIS_COUNT] =
                -RDB[ptr + SENS_EBLOCK_HIS_COUNT];
            }
          else if ((long)RDB[ptr + SENS_EBLOCK_HIS_COUNT] == 0)
            Die(FUNCTION_NAME, "WTF?");

          /* Next block */

          ptr = NextItem(ptr);
        }

      /* Next particle */

      part = NextItem(part);
    }

  /* Count banks for all threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    count = count + LIFOListSize(OMPPtr(DATA_PTR_EBLOCK_BANK, id), -1);

  /* Check that all events are accounted for */

  if (count != (long)RDB[DATA_EBLOCK_BANK_SZ])
    Die(FUNCTION_NAME, "Mismatch in number of event blocks");

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

      /* Get the first (0-gen) event block */

      ptr = (long)RDB[part + PARTICLE_PTR_SENS_EBLOCK];

      while (ptr > VALID_PTR)
        {
          /* Check if negative and switch to positive */

          if ((long)RDB[ptr + SENS_EBLOCK_HIS_COUNT] < 0)
            WDB[ptr + SENS_EBLOCK_HIS_COUNT] = -RDB[ptr + SENS_EBLOCK_HIS_COUNT];

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

  /* Get maximum number of threads */

  maxt = (long)RDB[DATA_OMP_MAX_THREADS];

  /* Reset id */

  id = OMP_THREAD_NUM;

  /* Loop over particles */

  while (part > VALID_PTR)
    {
        /* Reset count */

        n = 0;

        /* Get the event block */

        ptr = (long)RDB[part + PARTICLE_PTR_SENS_EBLOCK];

        /* Loop overe all event blocks */

        while (ptr > VALID_PTR)
          {
            /* Check with limit */

            if (n > (long)RDB[DATA_EVENT_MAX_GEN] - 1)
              {
                /* Update history count (should be removed) */

                WDB[ptr + SENS_EBLOCK_HIS_COUNT]--;
              }

            /* Check history counter */

            if ((long)RDB[ptr + SENS_EBLOCK_HIS_COUNT] > 0)
              {
                /* Pointer to next */

                ptr = NextItem(ptr);

              }
            else if ((long)RDB[ptr + SENS_EBLOCK_HIS_COUNT] == 0)
              {
                /* This was last particle that linked to this (old) */
                /* event block */

                /* Remove item */

                loc0 = ptr;

                /* Pointer to next */

                ptr = NextItem(loc0);

                /* Back to bank */

                EBlockToBank(loc0, id++);

                /* Start id from zero if gone through all ids */

                if (id > maxt - 1)
                  id = 0;

                /* Mark as removed */

                WDB[loc0 + SENS_EBLOCK_HIS_COUNT] = -2E+6;
              }
            else if ((long)RDB[ptr + SENS_EBLOCK_HIS_COUNT] == -2E+6)
              Die(FUNCTION_NAME, "This should not happen");
            else
              Die(FUNCTION_NAME, "WTF, negative? %ld %ld",
                  (long)RDB[ptr + SENS_EBLOCK_HIS_COUNT], ptr);

            /* Update counter */

            n++;
          }

        /* Next particle */

        part = NextItem(part);
    }

  /* Start parallel timer */

  StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private(part, id, prev, ptr)
#endif
  {
    /* Get Open MP thread id */

    id = OMP_THREAD_NUM;

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
        /* Check if particle is handled by some other thread */

        if (part % maxt != id)
          {
            /* Get next particle */

            part = NextItem(part);

            /* Cycle loop */

            continue;
          }

        /* Loop over events and detach pointers */

        prev = -1;

        ptr = (long)RDB[part + PARTICLE_PTR_SENS_EBLOCK];

        while (ptr > VALID_PTR)
          {
            /* Check if event was removed */

            if ((long)RDB[ptr + SENS_EBLOCK_HIS_COUNT] == -2E+6)
              {
                /* Detach previous */

                if (prev < VALID_PTR)
                  WDB[part + PARTICLE_PTR_SENS_EBLOCK] = NULLPTR;
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

        ptr = (long)RDB[part + PARTICLE_PTR_SENS_EBLOCK];
        while (ptr > VALID_PTR)
          {
            if ((long)RDB[ptr + SENS_EBLOCK_HIS_COUNT] < 1)
              Die(FUNCTION_NAME, "Error in history count");

            /* Pointer to next */

            ptr = NextItem(ptr);
          }
#endif

        /* Next particle */

        part = NextItem(part);
      }
  }

  /* Stop parallel timer */

  StopTimer(TIMER_OMP_PARA);

  /***************************************************************************/

  /***** Allocate memory for more events if needed ***************************/

  bankslow = 0;

  /* Calculate initial event bank size per thread */

  sz = (long)(RDB[DATA_EBLOCK_BANK_SZ]/RDB[DATA_OMP_MAX_THREADS]);

  /* Calculate threshold for low banks */

  cutoff = 0.2*sz;

  /* Pointer to event bank */

  loc0 = DATA_PTR_EBLOCK_BANK;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {

      if ((long)RDB[OMPPtr(loc0,id)] > VALID_PTR)
        {
          /* Check size */

          if (LIFOListSize(OMPPtr(loc0,id), cutoff) < cutoff)
            {
              bankslow = 1;
              break;
            }
        }
    }

  /* If some of the banks are low, first redistribute them */

  if (bankslow)
    {
      /***** Redistribute events *******************************************/

      loc0 = DATA_PTR_EBLOCK_BANK;

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

          /* Get average number of event blocks per thread */

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

  /* If some of the banks are still low, allocate more event blocks */

  if (bankslow)
    {
      /* Pointer to event bank */

      loc0 = DATA_PTR_EBLOCK_BANK;

      /* Get pointer to sensitivity block */

      sens = (long)RDB[DATA_PTR_SENS0];
      CheckPointer(FUNCTION_NAME, "(sens)", DATA_ARRAY, sens);

      /* Maximum event label */

      maxi = (long)RDB[sens + SENS_MAX_LABEL];

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
        {

          if ((long)RDB[OMPPtr(loc0,id)] > VALID_PTR)
            {
              /* Check size */

              if (LIFOListSize(OMPPtr(loc0,id), cutoff) < cutoff)
                {
                  /* Print */

                  Note(0, "Adjusting event score matrix bank size for thread %ld...", id);

                  /* Set expansion size */

                  np = cutoff;

                  /* Allow memory allocation */

                  Mem(MEM_ALLOW);

                  /* Allocate memory for event blocks */

                  for (n = 0; n < np; n++)
                    {
                      NewLIFOItem(OMPPtr(loc0, id), SENS_EBLOCK_BLOCK_SIZE);

                      /* Get pointer to newly created block */

                      loc1 = (long)RDB[OMPPtr(loc0,id)];

                      /* Allocate memory for scores */

                      ptr = ReallocMem(DATA_ARRAY, maxi);

                      /* Store pointer */

                      WDB[loc1 + SENS_EBLOCK_PTR_DATA] = (double)ptr;
                    }

                  /* Deny memory allocation */

                  Mem(MEM_DENY);

                  /* Put new size */

                  WDB[DATA_EBLOCK_BANK_SZ] = RDB[DATA_EBLOCK_BANK_SZ] + (double)np;
                }
            }
        }
    }

/***************************************************************************/

}
