/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reamulti.c                                     */
/*                                                                           */
/* Created:       2015/11/01 (JLe)                                           */
/* Last modified: 2018/10/01 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns microscopic cross section multiplier for reaction    */
/*              rea for calculating the total of macroscopic reaction mt     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReaMulti:"

/*****************************************************************************/

double ReaMulti(long rea, long mt, double E, long id)
{
  long ptr, nuc;
  double f;

  /* Check reaction pointer */
              
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
              
  /* Set multiplier */
  
  if (mt == MT_MACRO_FISSE)
    {
      /* JLe: muutettiin 9.3.2016 / 2.1.25 */
      /*
      f = RDB[rea + REACTION_Q]*RDB[DATA_NORM_U235_FISSE]/U235_FISSQ; 
      */

      /* Pointer to nuclide */
      
      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Get fission energy */
      
      f = FissE(rea, E, id);
      CheckValue(FUNCTION_NAME, "f", "", f/MEV, 0.0, 300.0);
    }
  else if (mt == MT_MACRO_NSF)
    {
      /* Get pointer to total nubar data */
      
      ptr = (long)RDB[rea + REACTION_PTR_TNUBAR];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Get multiplier */
      
      f = Nubar(ptr, E, id);
      
      /* Subtract delayed nubar data if delayed neutron */
      /* emission is off */
      
      if (((long)RDB[DATA_USE_DELNU] == NO) && 
          ((ptr = (long)RDB[rea + REACTION_PTR_DNUBAR]) > VALID_PTR))
        f = f - Nubar(ptr, E, id);
    }
  else if (mt == MT_MACRO_INLPRODXS)
    f = RDB[rea + REACTION_WGT_F] - 1.0;
  else
    f = RDB[rea + REACTION_WGT_F];

  /* Check value (may be zero if energy is out of range) */

  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, INFTY);

  /* Return value */

  return f;
}

/*****************************************************************************/
