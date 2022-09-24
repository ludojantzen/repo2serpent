/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : storetransmuxs.c                               */
/*                                                                           */
/* Created:       2011/06/08 (AIs)                                           */
/* Last modified: 2017/03/21 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Stores the one-group cross-sections and flux written by      */
/*              ClaculateTransmuXS to REACTION_PTR_TRANSMUXS to              */
/*              DEP_TRA_PS1/BOS/EOS and moves BOS to PS1 when needed         */
/*                                                                           */
/* Comments: This should maybe be part of CalculateTransmuXS                 */
/*           parempi että pidetään erillisenä aliohjelmana niin on selkeämpi */
/*           (JLe)                                                           */
/*                                                                           */
/*           Tätä ei kutsuta ollenkaan decay-stepeille versiosta 2.1.15      */
/*           eteenpäin (if-lauseke burnaterials.c:ssä) JLE 30.7.2013         */
/*                                                                           */
/* input: mat = pointer to material for which to store the values            */
/*        step = #step in this intervall                                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StoreTransmuXS:"

/*****************************************************************************/

void StoreTransmuXS(long mat, long step, long type, long id)
{
  long dep, ptr, rea, ciburn;
  double flx, xs, old_ave_flux, new_ave_flux, old_ave_xs, new_ave_xs, n;

  /* avoid compiler warnings */

  old_ave_flux = -INFTY; 
  new_ave_flux = -INFTY; 
  old_ave_xs = -INFTY; 
  n = -INFTY;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Get flag indicating, if we are currently in a CI-burn loop */

  ciburn = (long)RDB[DATA_BURN_CI_FLAG];

  /* Check burnup mode and burn flag */

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO) ||
      (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)))
    return;

  /* Get flux */

  if ((type != DEP_STEP_DEC_STEP) && (type != DEP_STEP_DEC_TOT))
    {
      ptr = (long)RDB[mat + MATERIAL_PTR_BURN_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      flx = Mean(ptr, 0);
    }
  else
    flx = 0.0;

  /* Divide by volume (may be zero if reprocessing is applied) */

  if ((ptr = (long)RDB[mat + MATERIAL_PTR_ACT_DET]) > VALID_PTR)
    flx = flx/RDB[ptr + DET_ABIN_VOL];
  else if (RDB[mat + MATERIAL_VOLUME] > 0.0)
    flx = flx/RDB[mat + MATERIAL_VOLUME];
  else if (flx > 0.0)
    Die(FUNCTION_NAME, "Zero volume (%s)", GetText(mat + MATERIAL_PTR_NAME));

  /* Use relaxed flux if available (already divided by volume) */

  if(RDB[mat + MATERIAL_BURN_FLUX_REL] != 0)
    flx = RDB[mat + MATERIAL_BURN_FLUX_REL];

  /* Check value */

  CheckValue(FUNCTION_NAME, "flx", "", flx, 0.0, INFTY);

  /* store value */

  if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    {
      /* Predictor */

      if (step == 0)
        WDB[mat + MATERIAL_BURN_FLUX_PS1] = -INFTY; 
      else
        WDB[mat + MATERIAL_BURN_FLUX_PS1] = RDB[mat + MATERIAL_BURN_FLUX_BOS];
      
      WDB[mat + MATERIAL_BURN_FLUX_BOS] = flx;
      WDB[mat + MATERIAL_BURN_FLUX_EOS] = -INFTY; 
    }
  else
    {
      /* Corrector */

      old_ave_flux = RDB[mat + MATERIAL_BURN_FLUX_AVE];
      CheckValue(FUNCTION_NAME, "old_ave_flux", "", old_ave_flux, 0.0, INFTY);

      /* Get number of corrector iteration for averaging */
      /* If not iterating, this will not average         */

      n = RDB[DATA_BURN_CI_I] + 1.0;

      /* If called from BurnMaterialsCI do not use averaged values */

      if (ciburn == YES)
        n = 1.0;

      /* First iteration/no iteration */
      /* -> no averaging              */

      new_ave_flux = (n-1.0)/n*old_ave_flux + (1.0/n)*flx;
      CheckValue(FUNCTION_NAME, "new_ave_flux", "", new_ave_flux, 0.0, INFTY);

      /* Store to EOS (can be zero) */

      WDB[mat + MATERIAL_BURN_FLUX_EOS] = new_ave_flux;

      /* Store to weighted average flux (for CI)    */
      /* Don't store if called from BurnMaterialsCI */

      if (ciburn == NO)
        WDB[mat + MATERIAL_BURN_FLUX_AVE] = new_ave_flux;
    }

  /**************************************************************************/

  /***** Transmutation list *************************************************/

  /* Pointer to transmutation list */
  
  dep = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];
  CheckPointer(FUNCTION_NAME, "(dep)", DATA_ARRAY, dep);

  /* Loop over list */

  while (dep > VALID_PTR)
    {
      /* Get the reaction and its XS */

      if ((type != DEP_STEP_DEC_STEP) && (type != DEP_STEP_DEC_TOT))
        {
          rea = (long)RDB[dep + DEP_TRA_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          xs = TestValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, id);

          /* use relaxed xs if available */

          if(RDB[dep + DEP_TRA_REL] != 0)
            xs = RDB[dep + DEP_TRA_REL];
        }
      else
        xs = 0.0;

      /* Check value */

      CheckValue(FUNCTION_NAME, "xs", "", xs, 0.0, INFTY);

      if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
        {
          /* Predictor */

          if (step == 0)
            WDB[dep + DEP_TRA_PS1] = -INFTY; 
          else
            WDB[dep + DEP_TRA_PS1] = RDB[dep + DEP_TRA_BOS];

          WDB[dep + DEP_TRA_BOS] = xs;
          WDB[dep + DEP_TRA_EOS] = -INFTY; 

          /* Store average one-group microxs for CI stopping criterion */
          /* to AV0. only used on first step */

          if (step == 0)
            WDB[dep + DEP_TRA_AV0] = xs;          
        }
      else
        {
          /* Corrector */

          /* Get old weighted average from memory   */
          /* zero if first iteration / no iteration */

          old_ave_xs =  WDB[dep + DEP_TRA_AVE];
          CheckValue(FUNCTION_NAME, "old_ave_xs", "", old_ave_xs, 0.0, INFTY);
      
          /* Get number of corrector iteration for averaging */
          /* If not iterating, this will not average         */

          n = RDB[DATA_BURN_CI_I] + 1.0;

          /* If called from BurnMaterialsCI do not use averaged values */

          if (ciburn == YES)
            n = 1.0;

          /* Average the xs (conserving reaction rates)        */
          /* In normal burnup calculation this will just store */
          /* flx*xs/flx = xs it just looks a bit cryptic       */

          /* NB: new_ave_flux may be zero if no neutrons have  */
          /* hit this depletion zone                           */

          if (new_ave_flux > 0.0)
            new_ave_xs =
              ((n-1.0)*old_ave_flux*old_ave_xs + flx*xs) / (new_ave_flux*n);
          else
            new_ave_xs = 0.0;

          CheckValue(FUNCTION_NAME, "new_ave_xs", "", new_ave_xs, 0.0, INFTY);

          /* Store to EOS if there was some contribution from this iteration */
          /* i.e. flux was non-zero this on iteration                        */

          if (flx > 0.0)
            WDB[dep + DEP_TRA_EOS] = new_ave_xs;
          else
            WDB[dep + DEP_TRA_EOS] = old_ave_xs;

          /* Store to weighted average*/

          if (ciburn == NO)
            {
              /* Store to weighted average if flux was nonzero */

              if (flx > 0.0)
                WDB[dep + DEP_TRA_AVE] = new_ave_xs;

              /* Store average one-group microxs for CI stopping criterion */

              /* Calculate new non-flux-weighted average of cross section */
              /* TODO: Pitäisikö tän olla kuitenkin tuo reaktiotaajuuden  */
              /* säilyttävä? VAST: tätä käytetään relaxoitujen nuklidi-   */
              /* tiheyksien kanssa */

              if (flx > 0.0)
                WDB[dep + DEP_TRA_AV1] = (n - 1.0)/n*RDB[dep + DEP_TRA_AV1] + 
                  1.0/n*xs;
            }
        }

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

      if ((type != DEP_STEP_DEC_STEP) && (type != DEP_STEP_DEC_TOT))
        {
          rea = (long)RDB[dep + DEP_TRA_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          xs = TestValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, id);

          /* use relaxed xs if available */

          if (RDB[dep + DEP_TRA_REL] != 0)
            xs = RDB[dep + DEP_TRA_REL];
        }
      else
        xs = 0.0;

      /* Check value */

      CheckValue(FUNCTION_NAME, "xs (fiss)", "", xs, 0.0, INFTY);

      if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
        {
          /* Predictor */

          if (step == 0)
            WDB[dep + DEP_TRA_PS1] = -INFTY;
          else
            WDB[dep + DEP_TRA_PS1] = RDB[dep + DEP_TRA_BOS];

          WDB[dep + DEP_TRA_BOS] = xs;
          WDB[dep + DEP_TRA_EOS] = -INFTY;

          /* Store average one-group microxs for CI stopping criterion */
          /* Only used on first step */

          if (step == 0)
            WDB[dep + DEP_TRA_AV0] = xs;
        }
      else
        {
          /* Get old weighted average from memory   */
          /* zero if first iteration / no iteration */

          old_ave_xs =  WDB[dep + DEP_TRA_AVE];
          CheckValue(FUNCTION_NAME, "old_ave_xs (fiss)", "", 
                     old_ave_xs, 0.0, INFTY);

          /* Get number of corrector iteration for averaging */
          /* If not iterating, this will not average         */

          n = RDB[DATA_BURN_CI_I] + 1.0;

          /* If called from BurnMaterialsCI do not use averaged values */

          if (ciburn == YES)
            n = 1.0;

          /* Average the xs (conserving reaction rates)        */
          /* In normal burnup calculation this will just store */
          /* flx*xs/flx = xs it just looks a bit cryptic       */

          /* NB: new_ave_flux may be zero if no neutrons have  */
          /* hit this depletion zone                           */

          if (new_ave_flux > 0.0)
            new_ave_xs =
              ((n - 1.0)*old_ave_flux*old_ave_xs + flx*xs) / (new_ave_flux*n);
          else
            new_ave_xs = 0.0;
          
          CheckValue(FUNCTION_NAME, "old_ave_xs (fiss)", "", 
                     old_ave_xs, 0.0, INFTY);

          /* Store to EOS */
          
          if (flx > 0.0)
            WDB[dep + DEP_TRA_EOS] = new_ave_xs;
          else
            WDB[dep + DEP_TRA_EOS] = old_ave_xs;

          /* Store to weighted average*/

          if (ciburn == NO)
            {
              /* Store to weighted average if flux was nonzero */

              if (flx > 0.0)
                WDB[dep + DEP_TRA_AVE] = new_ave_xs;

              /* Calculate new non-flux-weighted average of cross section */
              /* TODO: Pitäisikö tän olla kuitenkin tuo reaktiotaajuuden  */
              /* säilyttävä? VAST: tätä käytetään relaxoitujen nuklidi-   */
              /* tiheyksien kanssa */

              if (flx > 0.0)
                WDB[dep + DEP_TRA_AV1] = (n - 1.0)/n*RDB[dep + DEP_TRA_AV1] 
                  + 1.0/n*xs;
            }
        }
      
      /* Next reaction */

      dep = NextItem(dep);
    }

  /**************************************************************************/
}

/*****************************************************************************/
