/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : feifc.c                                        */
/*                                                                           */
/* Created:       2018/01/29 (BWe)                                           */
/* Last modified: 2018/02/23 (BWe)                                           */
/* Version:       2.1.31                                                      */
/*                                                                            */
/* Description: FE-defined material factors for multi-physics interface       */
/*                                                                            */
/* Comments: - This subroutine is intended to be used for user-defined        */
/*             density and temperature distributions with the multi-physics   */
/*             interface (types 31 & 32).                                     */
/*                                                                            */
/*             The function arguments are:                                    */
/*                                                                            */
/*             params  : Vector of FE parameters                              */
/*             coefs   : Vector of the FE coefficients                        */
/*             mat     : Pointer to the material where the collision occurs   */
/*             f       : The density factor, i.e. the ratio of local density  */
/*                       to the value given in the main input file (variable  */
/*                       adens). The factor must be between 0.0 and 1.0.      */
/*             T       : Local temperature. The value must be between Tmin    */
/*                       and Tmax, or zero, if the value is not used          */
/*             x, y, z : Coordinates where the collision occurs relative to   */
/*                       the root universe (usually universe "0")             */
/*             type    : specifies temperature or density                     */
/*                                                                            */
/*             Additional variables:                                          */
/*                                                                            */
/*             adens   : Atomic density of the interface material, as given   */
/*                       in the main input file                               */
/*             Tmin    : Minimum allowed temperature                          */
/*             Tmin    : Maximum allowed temperature                          */
/*                                                                            */
/******************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FEIFC:"

/******************************************************************************/

void FEIFC(const double *const params, const double *const coefs,
           long mat, double *f, double *T, double x, double y, double z,
                 long type, long threadID)
{
  double adens;
  double value;
#ifdef DEBUG
  double Tmin, Tmax;
#endif /* DEBUG */

  /* Evaluate the functional expansion */

  value = ExpandFE(params, coefs, x, y, z, threadID);

  if (value == INFTY)
    return;

  /* Reset state variables */

  *f = 1.0;
  *T = 0.0;

  switch (type)
  {
    case IFC_TYPE_FET_DENSITY:
    {
      adens = RDB[mat + MATERIAL_ADENS];

      value /= adens;

      if (value < 0.0)
        Warn(FUNCTION_NAME, "Density factor \"%ld\" too small; using 1.0.", *f);
      else if (value > 1.0)
        Warn(FUNCTION_NAME, "Density factor \"%ld\" too large; using 1.0.", *f);
      else
        *f = value;

      break;
    }

    case IFC_TYPE_FET_TEMP:
    {
#ifdef DEBUG
      Tmin = RDB[mat + MATERIAL_TMS_TMIN];
      Tmax = RDB[mat + MATERIAL_TMS_TMAX];

      if ((Tmin < INFTY) && (value < Tmin))
      {
        Warn(FUNCTION_NAME, "Temperature \"%lf\" too small; using %lf.", value, Tmin);

        *T = Tmin;
      }
      else if((Tmax > -INFTY) && (value > Tmax))
      {
        Warn(FUNCTION_NAME, "Temperature \"%lf\" too large; using %lf.", value, Tmax);

        *T = Tmax;
      }
      else
#endif /* DEBUG */
        *T = value;

      break;
    }

#ifdef DEBUG
    default:
      Die(FUNCTION_NAME, "Unsupported FE multiphysics interface type %ld", type);
#endif /* DEBUG */
  }
}
