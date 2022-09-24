/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : countdynsrc.c                                  */
/*                                                                           */
/* Created:       2015/05/15 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Counts total live weight reaching time interval boundary     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CountDynSrc:"

/*****************************************************************************/

void CountDynSrc()
{
  long loc0, ptr, id;
  double Wlive;

  /***************************************************************************/

  /***** Count banked neutrons and their weight ******************************/

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

#ifdef DNPRINT
  fprintf(outp, "countdynsrc.c -->\n");
#endif

  /* Reset total weight and source size */

  Wlive = 0.0;

  /* Loop over threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get particles from bank */

      /* Get pointer to bank */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to dummy */

      ptr = FirstItem(ptr);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to first item after dummy */

      while ((ptr = NextItem(ptr)) > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
            Die(FUNCTION_NAME, "Invalid particle type");

          /* Add to total weight */

          Wlive = Wlive + RDB[ptr + PARTICLE_WGT];
        }
    }

  /* Store live neutron weight */

  WDB[loc0 + PRECDET_W_LIVE] = Wlive;

#ifdef DNPRINT
  fprintf(outp, "<-- countdynsrc.c\n\n");
#endif

}
