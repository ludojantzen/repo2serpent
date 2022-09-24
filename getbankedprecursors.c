/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getbankedprecursors.c                          */
/*                                                                           */
/* Created:       2015/09/15 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Gets precursors from bank to PSOURCE, but leaves neutrons    */
/*              and gammas to bank                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetBankedPrecursors:"

/*****************************************************************************/

/* Name could also be precursor population control tms. */

void GetBankedPrecursors()
{
  long ptr, next, id, loc0;

  /***************************************************************************/

  if ((long)RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_NONE)
    return;

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

#ifdef DNPRINT
  fprintf(outp, "getbankedprecursors.c -->\n");
#endif

  /***** Move banked precursors to source ************************************/

  /* Loop over threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {

      /* Get pointer */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to first item */

      ptr = FirstItem(ptr);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Check type */

      if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_DUMMY)
        Die(FUNCTION_NAME, "Error in list");

      /* Get pointer to first item after dummy*/

      next = NextItem(ptr);

      /* Loop over bank */

      while (next > VALID_PTR)
        {
          ptr = next;

          next = NextItem(ptr);

          if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_PRECURSOR)
            {
              /* Remove precursor from bank */

              RemoveItem(ptr);

              /* Put precursor to precursor source (will be sorted later) */

              AddItem(DATA_PART_PTR_PSOURCE, ptr);

            }
        }

    }

#ifdef DNPRINT
  fprintf(outp, "<-- getbankedprecursors.c\n\n");
#endif

  /***************************************************************************/
}

/*****************************************************************************/
