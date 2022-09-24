/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : relaxtransmuxs.c                               */
/*                                                                           */
/* Created:       2014/07/23 (VVa)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Relaxes one-group cross-sections _inside_ a transport        */
/*              solution, when using coupled calculation                     */
/*              NB: Between step relaxation using Stochastic implicit euler  */
/*              is done in StoreTransmuXS                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RelaxTransmuXS:"

/*****************************************************************************/

void RelaxTransmuXS(long mat, long id)
{
  long dep, ptr, rea;
  double flx, xs, old_ave_flux, new_ave_flux, old_ave_xs, alpha, old_ave_rr;
  double rr, new_ave_rr;

  /* Check burnup mode */

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO) 
    return;

  /* avoid compiler warnings */

  old_ave_flux =-INFTY; 
  new_ave_flux =-INFTY; 
  old_ave_xs = -INFTY;
  old_ave_rr = -INFTY; 
  new_ave_rr = -INFTY;

  /* Get relaxation factor alpha from memory */
  
  alpha = RDB[DATA_SOL_REL_ALPHA];

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Get flux */

  ptr = (long)RDB[mat + MATERIAL_PTR_BURN_FLUX];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  flx = Mean(ptr, 0);

  /* Divide by volume (may be zero if reprocessing is applied) */

  if (RDB[mat + MATERIAL_VOLUME] > 0.0)
    flx = flx/RDB[mat + MATERIAL_VOLUME];
  else if (flx > 0.0)
    Die(FUNCTION_NAME, "Zero volume (%s)", GetText(mat + MATERIAL_PTR_NAME));

  /* Check value */

  CheckValue(FUNCTION_NAME, "flx", "", flx, 0.0, INFTY);

  /* Get old relaxed flux flux */

  if(RDB[DATA_SOL_REL_ITER] > 0)
    old_ave_flux = RDB[mat + MATERIAL_BURN_FLUX_REL];
  else
    old_ave_flux = 0.0;

  /* Check value */

  CheckValue(FUNCTION_NAME, "old_ave_flux", "", old_ave_flux, 0.0, INFTY);

  /* If this is the first iteration use flx as an initial guess for future  */
  /* Otherwise relax */

  if (old_ave_flux == 0)
    new_ave_flux = flx;
  else
    new_ave_flux = old_ave_flux - alpha*3.0/4.0*(old_ave_flux - flx);  

  /* Check value */
  
  CheckValue(FUNCTION_NAME, "new_ave_flux", "", new_ave_flux, 0.0, INFTY);

  /* Store relaxed flux */

  WDB[mat + MATERIAL_BURN_FLUX_REL] = new_ave_flux;
  
  /**************************************************************************/

  /***** Transmutation list *************************************************/

  /* Pointer to transmutation list */
  
  dep = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];
  CheckPointer(FUNCTION_NAME, "(dep)", DATA_ARRAY, dep);

  /* Loop over list */

  while (dep > VALID_PTR)
    {
      /* Get the reaction and its XS */

      rea = (long)RDB[dep + DEP_TRA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      xs = TestValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, id);

      /* Calculate reaction rate */

      rr = xs*flx;

      /* Get old relaxed cross section */

      if(RDB[DATA_SOL_REL_ITER] > 0)
        old_ave_xs = RDB[dep + DEP_TRA_REL];
      else
        old_ave_xs = 0.0;

      /* Calculate old relaxed reaction rate */

      old_ave_rr = old_ave_xs*old_ave_flux;
      CheckValue(FUNCTION_NAME, "old_ave_rr", "", old_ave_rr, 0.0, INFTY);

      /* Relax reaction rate */

      if (old_ave_xs == 0)
        new_ave_rr = rr;
      else
        new_ave_rr = old_ave_rr - alpha*3.0/4.0*(old_ave_rr - rr);  

      CheckValue(FUNCTION_NAME, "new_ave_rr", "", new_ave_rr, 0.0, INFTY);

      /* Store relaxed xs */

      if (new_ave_flux > 0.0)
        WDB[dep + DEP_TRA_REL] = new_ave_rr/new_ave_flux;
      else
        WDB[dep + DEP_TRA_REL] = old_ave_xs;

      /* Next reaction */

      dep = NextItem(dep);
    }

  /**************************************************************************/

  /***** Fission list *******************************************************/

  /* Pointer to fission list (pointteri voi olla null jos materiaali */
  /* on ei-fissiili) */

  dep = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];
  
  /* Loop over list */

  while (dep > VALID_PTR)
    {
      /* Get the reaction and its XS */

      rea = (long)RDB[dep + DEP_TRA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      xs = TestValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, id);

      /* Calculate reaction rate */

      rr = xs*flx;

      /* Get old relaxed cross section */

      if(RDB[DATA_SOL_REL_ITER] > 0)
        old_ave_xs = RDB[dep + DEP_TRA_REL];
      else
        old_ave_xs = 0.0;

      /* Calculate old relaxed reaction rate */

      old_ave_rr = old_ave_xs*old_ave_flux;
      CheckValue(FUNCTION_NAME, "old_ave_rr", "", old_ave_rr, 0.0, INFTY);

      /* Relax reaction rate */

      if (old_ave_xs == 0)
        new_ave_rr = rr;
      else
        new_ave_rr = old_ave_rr - alpha*3.0/4.0*(old_ave_rr - rr);  

      CheckValue(FUNCTION_NAME, "new_ave_rr", "", new_ave_rr, 0.0, INFTY);

      /* Store relaxed xs */

      if (new_ave_flux > 0.0)
        WDB[dep + DEP_TRA_REL] = new_ave_rr/new_ave_flux;
      else
        WDB[dep + DEP_TRA_REL] = old_ave_xs;

      /* Next reaction */

      dep = NextItem(dep);
    }

  /**************************************************************************/

}

/*****************************************************************************/
