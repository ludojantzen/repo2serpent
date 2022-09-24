/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dtmajorant.c                                   */
/*                                                                           */
/* Created:       2012/10/09 (JLe)                                           */
/* Last modified: 2018/06/21 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns delta-tracking majorant for neutrons and photons     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DTMajorant:"

/*****************************************************************************/

double DTMajorant(long type, double E, long id)
{
  long rea, rea1, ptr, nuc, mat;
  double maj, xs, f;

  /* Avoid compiler warning */

  rea = -1;

  /* Check type and get reaction pointer */

  if (type == PARTICLE_TYPE_NEUTRON)
    rea = (long)RDB[DATA_PTR_MAJORANT];
  else if (type == PARTICLE_TYPE_GAMMA)
    rea = (long)RDB[DATA_PTR_PHOTON_MAJORANT];
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Check Pointer (aikaisemmin on oletettu että st-moodissa   */
  /* majoranttia ei ole laskettu, joten toi voi mennä nulliksi */
  /* siitä syystä) */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get majorant cross section */

  maj = MajorantXS(rea, E, id);

  /* Add alpha cross section */

  if (type == PARTICLE_TYPE_NEUTRON)
    maj = maj + AlphaXS(E);

  /* Add extra cross sections */

  ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
  while (ptr > VALID_PTR)
    {
      /* Get factor */

      f = RDB[ptr + MAJORANT_EXTRA_FRAC];
      CheckValue(FUNCTION_NAME, "f", "", f, 0.0, INFTY);

      /* Check for zero */

      if (f == 0.0)
        {
          /* Pointer to next */

          ptr = NextItem(ptr);

          /* Cycle loop */

          continue;
        }

      /* Avoid compiler warning */

      xs = -1.0;

      /* Check type and get cross section */

      if ((nuc = (long)RDB[ptr + MAJORANT_EXTRA_PTR_NUC]) > VALID_PTR)
        {
          /* Pointer to reaction */

          rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Get pointer to separate tmp majorant if not using MGXS */
          /* Multi-group majorant mode stores TMP_MAJORANT directly */
          /* to NUCLIDE_PTR_TOTXS */

          if (RDB[DATA_OPTI_MG_MODE] == (double)NO)
            if ((rea1 = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT]) > VALID_PTR)
              rea = rea1;

          if (RDB[DATA_OPTI_MG_MODE] == (double)YES)
            {
              /* Multi-group majorant mode */

              /* Get pointer to MGXS */

              xs = MGXS(rea, E, -1);
            }
          else
            {
              /* Continuous energy mode */

              /* Get microscopic cross section */

              xs = MicroXS(rea, E, id);
            }
        }
      else if ((mat = (long)RDB[ptr + MAJORANT_EXTRA_PTR_MAT]) > VALID_PTR)
        {
          /* Pointer to reaction */

          rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Get cross section */

          xs = MacroXS(rea, E, id);
        }
      else
        Die(FUNCTION_NAME, "Null extra pointer");

      /* Add to majorant */

      maj = maj + f*xs;

      /* Pointer to next */

      ptr = NextItem(ptr);
    }

  /* Return cross section */

  return maj;
}

/*****************************************************************************/
