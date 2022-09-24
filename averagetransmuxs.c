/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : averagetransmuxs.c                             */
/*                                                                           */
/* Created:       2011/06/08 (AIs)                                           */
/* Last modified: 2017/02/08 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: calculates substep averaged one-group cross-sections and     */
/*              flux, which are stored to REACTION_PTR_TRANSMUXS and         */
/*              MATERIAL_BURN_FLUX_SSA (these are read by MakeBurnMatrix)    */
/*                                                                           */
/* Comments: Lin.ext. and quad.int. commonly produce necative xs for         */ 
/*           reactions with very high treshold energy and hence large        */
/*           statistical variations, but since this implies neclicible flux  */
/*           above the treshold, there is no problem. BOS value is used when */
/*           extrapolation/interpolation yields negative value.              */
/*                                                                           */
/* TODO?: -Should we set the negat. xs to zero instead?                      */
/*        -More sophisticated check to see if non-neclicible reactions get   */
/*         negative values (implying error in the code/input) ?              */
/*                                                                           */
/*                                                                           */
/* input: mat = pointer to material for which to store the values            */
/*        t1 = beginning time of substep (t=0 at BOS)                        */
/*        t2 = ending time of substep                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AverageTransmuXS:"

/*****************************************************************************/

void AverageTransmuXS(long mat, double t1, double t2, long id)
{
  long dep, rea, ptr;
  double wPS1, wBOS, wEOS, flx, xs;

  /* Check decay only mode */

  if ((long)RDB[DATA_BURN_DECAY_CALC] == YES)
    return;

  /* average (integrate) general 2nd order polynomial from t1 to t2      */
  /* This works with lower order fits too, the higher order coefficient, */
  /* or rather the weights for them just are 0 */
  
  wPS1 = RDB[DATA_BURN_FIT_C2W1]*(t2*t2*t2-t1*t1*t1)/(t2-t1)/3
    + RDB[DATA_BURN_FIT_C1W1]*(t2*t2-t1*t1)/(t2-t1)/2
    + RDB[DATA_BURN_FIT_C0W1];
  
  wBOS = RDB[DATA_BURN_FIT_C2W2]*(t2*t2*t2-t1*t1*t1)/(t2-t1)/3
    + RDB[DATA_BURN_FIT_C1W2]*(t2*t2-t1*t1)/(t2-t1)/2
    + RDB[DATA_BURN_FIT_C0W2];
  
  wEOS = RDB[DATA_BURN_FIT_C2W3]*(t2*t2*t2-t1*t1*t1)/(t2-t1)/3
    + RDB[DATA_BURN_FIT_C1W3]*(t2*t2-t1*t1)/(t2-t1)/2
    + RDB[DATA_BURN_FIT_C0W3];
  
  CheckValue(FUNCTION_NAME, "sum of weights - 1", "", wPS1 + wBOS + wEOS - 1,
             -1E-12, 1E+12);

  /***************************************************************************/

  /*** flux ******************************************************************/
  
  WDB[mat + MATERIAL_BURN_FLUX_SSA] = wPS1*RDB[mat + MATERIAL_BURN_FLUX_PS1]
    + wBOS*RDB[mat + MATERIAL_BURN_FLUX_BOS] 
    + wEOS*RDB[mat + MATERIAL_BURN_FLUX_EOS];
  
  if (RDB[mat + MATERIAL_BURN_FLUX_SSA] < 0.0)
    WDB[mat + MATERIAL_BURN_FLUX_SSA] = RDB[mat + MATERIAL_BURN_FLUX_BOS];

  /* Override for activation detectors */

  if ((ptr = (long)RDB[mat + MATERIAL_PTR_ACT_DET]) > VALID_PTR)
    if ((flx = RDB[ptr + DET_ABIN_FLUX]) >= 0.0)
      WDB[mat + MATERIAL_BURN_FLUX_SSA] = flx;

  CheckValue(FUNCTION_NAME, "Substep averaged flux", "", 
             RDB[mat + MATERIAL_BURN_FLUX_SSA] , 0.0, INFTY);

  /***************************************************************************/
  
  /*** Transmutation list ****************************************************/
  
  /* Pointer to transmutation list */

  dep = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST]; 

  CheckPointer(FUNCTION_NAME, "(dep)", DATA_ARRAY, dep);

  while(dep > VALID_PTR)
    {
      /* Get the reaction pointer */

      rea = (long)RDB[dep + DEP_TRA_PTR_REA];

      /* weighted (substep average) cross-section */

      xs = wPS1*RDB[dep + DEP_TRA_PS1] + wBOS*RDB[dep + DEP_TRA_BOS]
        + wEOS*RDB[dep + DEP_TRA_EOS];
      
      if (xs < 0.0)
        xs = RDB[dep + DEP_TRA_BOS];

      StoreValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, xs, id);
      
      /* Next reaction */

      dep = NextItem(dep); 
    }
  
  /***************************************************************************/

  /*** Fission list **********************************************************/
  
  /* Pointer to fission list */

  dep = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST]; 
  
  while(dep > VALID_PTR)
    {
      /* Get the reaction pointer */
      
      rea = (long)RDB[dep + DEP_TRA_PTR_REA];
      
      /* weighted (substep average) cross-section */
      
      xs = wPS1*RDB[dep + DEP_TRA_PS1] + wBOS*RDB[dep + DEP_TRA_BOS]
        + wEOS*RDB[dep + DEP_TRA_EOS]; 
      
      if (xs < 0.0)
        xs = RDB[dep + DEP_TRA_BOS];

      StoreValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, xs, id);

      /* Next reaction */
      
      dep = NextItem(dep); 
    }

  /***************************************************************************/
}

/*****************************************************************************/
