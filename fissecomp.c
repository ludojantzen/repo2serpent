/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fissecomp.c                                    */
/*                                                                           */
/* Created:       2018/09/28 (RTu)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates energy released in fission based on ENDF MT458    */
/*              data.                                                        */
/* Comments: Return value is in MeV:s. Nubar at "E = 0" is the constant term */
/*           in the polynomial or the first tabulated value of total nubar.  */
/*           Based on the ENDF manual it is not clear at which energy nubar  */
/*           at "E = 0" should be evaluated.                                 */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FissEComp:"

/*****************************************************************************/

double FissEComp(long fisse_ptr, long nu_ptr, long comp, double E, long id)
{
  long loc0, NPLY, ptr, type, i, loc, NR, scheme, j, NP;
  double E0, E1, edep0, edep, nu0, nu, x, val0, val1;

  /* Check pointer to fission energy deposition data */

  CheckPointer(FUNCTION_NAME, "(fisse_ptr)", DATA_ARRAY, fisse_ptr);

  /* Check that valid component is requested */

  CheckValue(FUNCTION_NAME, "comp", "", comp, 0, FISSE_COMP_BLOCK_SIZE - 1);

  /* Pointer to components */

  loc0 = (long)RDB[fisse_ptr + FISSE_DATA_COMP];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Avoid compiler warning */

  nu = 0.0;
  nu0 = 0.0;

  /* Initialize edep */

  edep = 0.0;

  /* Process different formats */

  if ((NPLY = (long)RDB[fisse_ptr + FISSE_DATA_NPLY]) == 0)
    {
      /* Get type */

      if ((long)RDB[fisse_ptr + FISSE_DATA_LFC] == 0)
        type = FISSE_TYPE_SHER_BECK;
      else
        {
          /* Get pointer to array with component types */

          ptr = (long)RDB[fisse_ptr + FISSE_DATA_COMP_TYPES];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get type from array */

          type = (long)RDB[ptr + comp];
        }

      if (type == FISSE_TYPE_SHER_BECK)
        {
          /* Get value at E = 0 */

          edep0 = RDB[loc0 + comp];

          /* Get nubar data if needed */

          if ((comp == FISSE_COMP_ENP) || (comp == FISSE_COMP_ET))
            {
              /* Nubar at E */

              CheckPointer(FUNCTION_NAME, "(nu_ptr)", DATA_ARRAY, nu_ptr);
              nu = Nubar(nu_ptr, E, id);
              CheckValue(FUNCTION_NAME, "nu", "", nu, 0.0, 10.0);

              /* Nubar at E = 0 */

              nu0 = RDB[fisse_ptr + FISSE_DATA_NUBAR0];
              CheckValue(FUNCTION_NAME, "nu0", "", nu0, 0.0, 10.0);
            }

          /* Calculate component */

          switch (comp)
            {
            case FISSE_COMP_ENP:
              edep = edep0 + 1.307*E - 8.07*(nu - nu0);
              break;
            case FISSE_COMP_EGD:
              edep = edep0 - 0.075*E;
              break;
            case FISSE_COMP_EB:
              edep = edep0 - 0.075*E;
              break;
            case FISSE_COMP_ENU:
              edep = edep0 - 0.1*E;
              break;
            case FISSE_COMP_ET:
              edep = edep0 + 1.057*E - 8.07*(nu - nu0);
              break;
            default:
              edep = edep0;
            }
        }
      else if (type == FISSE_TYPE_TAB)
        {
          /* Get pointer to tab data */

          loc = (long)RDB[loc0 + comp];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, loc);

          /* Get pointer to incident energy array */

          ptr = (long)RDB[loc + FISSE_TAB_ERG];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get number of points */

          NP = (long)RDB[loc + FISSE_TAB_NP];

          /* Find index from grid */

          i = SearchArray(&RDB[ptr], E, NP);

          /* Check if energy point was found */

          if (i < 0)
            {
              /* Check co-incident boundary points */

              if (fabs(RDB[ptr] - E) < 1E-9*RDB[ptr])
                {
                  i = 0;
                }
              else if (fabs(RDB[ptr + NP - 1] - E) < 1E-9*RDB[ptr + NP - 1])
                {
                  i = NP - 2;
                }
              else
                {
                  Warn(FUNCTION_NAME, "Energy point was not found on the grid.\n"
                       "comp = %ld E = %11.5E Emin = %11.5E Emax = %11.5E",
                       comp, E, RDB[ptr], RDB[ptr + NP - 1]);
                  return 0.0;
                }
            }

          /* Get energies */

          E0 = RDB[ptr + i];
          E1 = RDB[ptr + i + 1];

          /* Get pointer to values */

          ptr = (long)RDB[loc + FISSE_TAB_EDEP];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get values */

          val0 = RDB[ptr + i];
          val1 = RDB[ptr + i + 1];

          /* Get number of interpolation ranges */

          NR = (long)RDB[loc + FISSE_TAB_NR];

          /* Get interpolation scheme */

          /* Only one range so same scheme is used for all energies */

          if (NR == 1)
            scheme = (long)RDB[loc + FISSE_TAB_INT];
          else
            {
              /* Multiple ranges */

              /* Get pointer to NBT */

              ptr = (long)RDB[loc + FISSE_TAB_NBT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Find range */

              for (j = 0; j < NR; j++)
                if ((long)RDB[ptr + j] > (i + 1))
                  break;

              if (j == NR)
                {
                  Warn(FUNCTION_NAME, "Interpolation scheme not found");
                  return 0.0;
                }

              /* Get pointer to scheme data */

              ptr = (long)RDB[loc + FISSE_TAB_INT];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Get scheme */

              scheme = (long)RDB[ptr + j];
            }

          /* Interpolate */

          edep = ENDFInterp(scheme, E, E0, E1, val0, val1);
        }
      else
        Die(FUNCTION_NAME, "Invalid type %ld", type);
    }
  else
    {
      /* Polynomial format */

      x = 1.0;
      edep = 0.0;

      /* Get pointer to coefficient data */

      ptr = (long)RDB[loc0 + comp];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      for (i = 0; i < NPLY + 1; i++)
        {
          edep += x*RDB[ptr + i];
          x = x*E;
        }
    }

  /* Check value in MeV:s */

  CheckValue(FUNCTION_NAME, "edep", "", edep, 0.0, 300.0);

  /* Return value */

  return edep;
}

/*****************************************************************************/
