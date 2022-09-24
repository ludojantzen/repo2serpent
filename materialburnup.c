/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : materialburnup.c                               */
/*                                                                           */
/* Created:       2011/07/07 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Calculates material-wise burnup for time step                */
/*                                                                           */
/* Comments: -This has been designed to be as consistent as possible with    */
/*            the calculation of the total burnup (in setdepstepsize.c) as   */
/*            possible without material-wise fission rate tallies. (AIs)     */
/*           -Note that both these, and the total burnup updates, are only   */
/*            as accurate as the rest of the calculation. To get the exact   */
/*            burnups, one would have to "tally" the exact number of fissions*/
/*            in the depletion calculations by using virtual nuclides (AIs)  */
/*           -Because of the above, material burnups do not sum up exactly   */
/*            to the total burnup. Not even with only one material region.   */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/* input:                                                                    */
/*    mat:  "pointer" to the material block                                  */
/*    Nbos: beginning of step (not substep) material composition array       */
/*    Neos: End of step material composition array                           */
/*    t1,t2:Beginning and end times of this substep with t=0 at BOS          */
/*    ss:   current substep (ss=0 being the first)                           */
/*    id:   Corrector iteration in dufek's scheme                            */
/*                                                                           */
/* NOTE: Nbos and Neos only need to be correct on the first substep of the   */
/*       predictor and corrector, respectively.                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MaterialBurnup:"

/*****************************************************************************/

void MaterialBurnup(long mat, double *Nbos, double *Neos, double t1, 
                    double t2, long ss, long id)
{
  long dep, nuc, rea, i;
  double pow, flx, adens, mdens, xs, Q, c0, c1, c2;

  /* Check material pointer */
  
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Reset power */

  pow = 0.0;

  /* Check values */
  
  flx = RDB[mat + MATERIAL_BURN_FLUX_BOS];
  mdens = RDB[mat + MATERIAL_INI_FISS_MDENS];
  
  if ((mdens <= 0.0) || (flx <= 0.0))
    return;

  /***************************************************************************/

  /***** calculate BOS power and store prev step power ***********************/

  if ((ss == 0) && ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP))
    {      
      /* Loop over fission reactions */

      /* Unknown: (tota listaa ei oo jos on se totfiss-juttu?) */

      dep = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];
      while (dep > VALID_PTR)
        {
          /* BOS flux */

          flx = RDB[mat + MATERIAL_BURN_FLUX_BOS];
          CheckValue(FUNCTION_NAME, "flx1", "", flx, 0.0, INFTY);

          /* Pointer to reaction */
          
          rea = (long)RDB[dep + DEP_TRA_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          
          /* Pointer to nuclide data */
          
          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
          
          /* Get index to composition vector */
          
          if ((i = (long)TestValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX, 
                                       (double)mat, id)) < 0)
            Die(FUNCTION_NAME, "i < 0");
          
          /* Get atomic density from composition vector */
          
          adens = Nbos[i];
          CheckValue(FUNCTION_NAME, "adens1", "", adens, 0.0, INFTY);

          /* Get deposited energy (JLe: muutettu 9.3.2016 / 2.1.26) */
          /*
          Q = RDB[rea + REACTION_Q]*RDB[DATA_NORM_U235_FISSE]/U235_FISSQ;
          */
          
          Q = RDB[nuc + NUCLIDE_FISSE];
          CheckValue(FUNCTION_NAME, "Q1", "", Q, 0.0, INFTY);

          /* Get BOS cross-section */
          
          xs = RDB[dep + DEP_TRA_BOS];
          CheckValue(FUNCTION_NAME, "xs1", "", xs, 0.0, INFTY);
          
          /* Add power from this reaction */
          
          pow = pow + adens*flx*xs*Q; /* W/cm^3 */
          CheckValue(FUNCTION_NAME, "pow1", "", pow, 0.0, INFTY);

          /* Next reaction */
          
          dep = NextItem(dep);
        }
      
      /* Store previous BOS value as PS value */
      
      WDB[mat + MATERIAL_BURN_POW_PS1] = RDB[mat + MATERIAL_BURN_POW_BOS];
      
      /* Store the new BOS value */
  
      WDB[mat + MATERIAL_BURN_POW_BOS] = pow;
    }
  
  /***************************************************************************/
  
  /*** Calculate EOS power ***************************************************/
  
  if ((ss == 0) && ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
    { 
      /* Loop over fission reactions */
      
      dep = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];
      while (dep > VALID_PTR)
        {
          /* EOS flux */
          
          flx = RDB[mat + MATERIAL_BURN_FLUX_EOS];
          CheckValue(FUNCTION_NAME, "flx2", "", flx, 0.0, INFTY);

          /* Pointer to reaction */
          
          rea = (long)RDB[dep + DEP_TRA_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          
          /* Pointer to nuclide data */
          
          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
          
          /* Get index to composition vector */
          
          if ((i = (long)TestValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX, 
                                       (double)mat, id)) < 0)
            Die(FUNCTION_NAME, "i < 0");
          
          /* Get atomic density from composition vector */
          
          adens = Neos[i];
          CheckValue(FUNCTION_NAME, "adens2", "", adens, 0.0, INFTY);

          /* Get deposited energy (JLe: muutettu 9.3.2016 / 2.1.26) */
          /*
          Q = RDB[rea + REACTION_Q]*RDB[DATA_NORM_U235_FISSE]/U235_FISSQ;
          */

          Q = RDB[nuc + NUCLIDE_FISSE];
          CheckValue(FUNCTION_NAME, "Q2", "", Q, 0.0, INFTY);

          /* Get EOS cross-section */
          
          xs = RDB[dep + DEP_TRA_EOS];
          CheckValue(FUNCTION_NAME, "xs2", "", xs, 0.0, INFTY);
          
          /* Add power from this reaction */
          
          pow = pow + adens*flx*xs*Q; /* W/cm^3 */
          CheckValue(FUNCTION_NAME, "pow2", "", pow, 0.0, INFTY);

          /* Next reaction */
          
          dep = NextItem(dep);
        }

      /* Store the new EOS value */
      
      WDB[mat + MATERIAL_BURN_POW_EOS] = pow;
    }

  /***************************************************************************/

  /* burnup is not updated on the predictor if there is a corrector */

  if (!((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE) &&
      ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP))
    return;

  /* Burnup is only incremented on the first corrector iteration */

  if ((RDB[DATA_BURN_SIE] == (double)YES) && (DATA_BURN_CI_I > 0))
    return;

  /***************************************************************************/

  /****** Update burnup ******************************************************/
  
  /* Polynomial coefficients for material power (see depletionpolyfit.c) */

  c2 = RDB[DATA_BURN_FIT_C2W1]*RDB[mat + MATERIAL_BURN_POW_PS1]
     + RDB[DATA_BURN_FIT_C2W2]*RDB[mat + MATERIAL_BURN_POW_BOS]
     + RDB[DATA_BURN_FIT_C2W3]*RDB[mat + MATERIAL_BURN_POW_EOS];
  
  c1 = RDB[DATA_BURN_FIT_C1W1]*RDB[mat + MATERIAL_BURN_POW_PS1]
     + RDB[DATA_BURN_FIT_C1W2]*RDB[mat + MATERIAL_BURN_POW_BOS]
     + RDB[DATA_BURN_FIT_C1W3]*RDB[mat + MATERIAL_BURN_POW_EOS];
  
  c0 = RDB[DATA_BURN_FIT_C0W1]*RDB[mat + MATERIAL_BURN_POW_PS1]
     + RDB[DATA_BURN_FIT_C0W2]*RDB[mat + MATERIAL_BURN_POW_BOS]
     + RDB[DATA_BURN_FIT_C0W3]*RDB[mat + MATERIAL_BURN_POW_EOS];

  CheckValue(FUNCTION_NAME, "c0", "", c0, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c1", "", c1, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c2", "", c2, -INFTY, INFTY);
  
  /* Integral of general 2nd order polynomial from t1 to t2 */
  /* This is the energy released in the material during the substep [J/cm3] */

  pow = c2/3*(t2*t2*t2 - t1*t1*t1) + c1/2*(t2*t2 - t1*t1) + c0*(t2 -t1);
  
  /* NOTE: in very noisy calculation or with very large changes in step length
     it might be possible to interpolate/extrapolate a negative value, it might
     be desirable to change this so thet power is set to zero or BOS value if
     this happens, so quick tests can be done (AIs) */

  CheckValue(FUNCTION_NAME,"pow","Fission energy in material",pow,0.0,INFTY);

  /* Increment the burnup [MWd/kgHM] */

  WDB[mat + MATERIAL_BURNUP] = RDB[mat + MATERIAL_BURNUP]
    + pow/mdens/86400/1.0E3;
  
  /***************************************************************************/
}

/*****************************************************************************/
