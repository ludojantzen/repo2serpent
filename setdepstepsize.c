/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setdepstepsize.c                               */
/*                                                                           */
/* Created:       2011/05/24 (JLe)                                           */
/* Last modified: 2017/06/07 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Sets the length of depletion step in seconds.                */
/*                                                                           */
/* Comments: - Teho/FMASS tarkastukset Calculate_t_from_bu:ssa pitäisi ehkä  */
/*             siirtää aikaisempaan (AIs 22.6.2011)                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetDepStepSize:"

double Calculate_bu_from_t(double);
double Calculate_t_from_bu(double);

/*****************************************************************************/

void SetDepStepSize(long dep, long step)
{
  long ptr, loc0, type;
  double t, pow, bu;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(dep)", DATA_ARRAY, dep);

  /* Check interval */

  CheckValue(FUNCTION_NAME, "step", "", step, 0, 
             (long)RDB[dep + DEP_HIS_N_STEPS] - 1);

  /* Pointer to steps */

  loc0 = (long)RDB[dep + DEP_HIS_PTR_STEPS];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get step type */

  type = (long)RDB[dep + DEP_HIS_STEP_TYPE];

  /* Avoid compiler warning */

  t = -1.0;
  bu = -1.0;

  /* Get power in burnable materials [W] */

  ptr = (long)RDB[RES_TOT_NEUTRON_POWER];
  pow = Mean(ptr, 1);

  /* Check value */

  CheckValue(FUNCTION_NAME, "pow", "", pow, 0.0, INFTY);
  
  /* store power */

  if ( (long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP )
    {
      /* Store the value from previous step */

      if (step == 0)
        {
          /* No prev. step value at first step */

          WDB[DATA_BURN_POW_PS1] = -INFTY;
        }
      else
        WDB[DATA_BURN_POW_PS1] = RDB[DATA_BURN_POW_BOS];

      /* Store BOS value */

      WDB[DATA_BURN_POW_BOS] = pow;

      /* No EOS value at predictor */

      WDB[DATA_BURN_POW_EOS] = -INFTY; 
    }
  else
    {
      /* Corrector step, store EOS value */

      WDB[DATA_BURN_POW_EOS] = pow; /* MARK (AIs) */
    } 

  /* Check step type */
  
  switch (type)
    {
    case DEP_STEP_DAY_STEP:
    case DEP_STEP_DEC_STEP:
    case DEP_STEP_ACT_STEP:
      {
        /* Step is given in days, convert to [s] and [MWd/kgHM] */

        t = RDB[loc0 + step]*60.0*60.0*24.0;
        bu = Calculate_bu_from_t(t);

        /* Break case */

        break; 
      }
    case DEP_STEP_DAY_TOT:
    case DEP_STEP_DEC_TOT:
    case DEP_STEP_ACT_TOT:
      {
        /* Step is given in days, convert to [s] and [MWd/kgHM] */

        if (RDB[loc0 + step]*60.0*60.0*24.0 <= RDB[DATA_BURN_CUM_BURNTIME])
          Error(dep, "Listed value %1.5E is below or equal to cumulative time",
                RDB[loc0 + step]);

        t = RDB[loc0 + step]*60.0*60.0*24.0 - RDB[DATA_BURN_CUM_BURNTIME];
        bu = Calculate_bu_from_t(t); 

        /* Break case */

        break;
      }
    case DEP_STEP_BU_STEP:
      {
        /* Step is given in burnup */

        bu = RDB[loc0 + step];
        t = Calculate_t_from_bu(bu);

        /* Break case */

        break;
      }
    case DEP_STEP_BU_TOT:
      {
        /* Step is given in cumulative burnup */
        
        if (RDB[loc0 + step] <=  RDB[DATA_BURN_CUM_BURNUP])
          Error(dep, 
                "Listed value %1.5E is below or equal to cumulative burnup",
                RDB[loc0 + step]);
        
        bu = RDB[loc0 + step] - RDB[DATA_BURN_CUM_BURNUP];
        t = Calculate_t_from_bu(bu);

        /* Break case */

        break;
      }
    default:
      Die(FUNCTION_NAME, "Invalid step type");
    }

  /* Check values */

  CheckValue(FUNCTION_NAME, "t", "", t, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "bu", "", bu, 0.0, INFTY);
  
  /* Set step size (this is read by burnmaterials) */

  WDB[DATA_BURN_TIME_INTERVAL] = t;
  WDB[DATA_BURN_BURNUP_INTERVAL] = bu;

  /* These are needed by depletionpolyfit.c on next step / corrector.      */
  /* The PS1 length will not be set if PC method is converted to CE due to */
  /* zero flux, but it wont be needed with CE so this is not a problem.    */

  if((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    WDB[DATA_BURN_PRED_LENGTH] = RDB[DATA_BURN_TIME_INTERVAL];
  
  if ( (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) && 
        ((long)RDB[DATA_BURN_CI_LAST] == YES)) ||
      ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE))
    WDB[DATA_BURN_PS1_LENGTH] = RDB[DATA_BURN_TIME_INTERVAL];

}

/*****************************************************************************/

/*****************************************************************************/
/* Calculate the mean power and resulting burnup when step length is known.  */
/* In: step length [s], Out: step length [MWd/kgHM]  (AIs)                   */
/*****************************************************************************/

double Calculate_bu_from_t(double t)
{
  double c2, c1, c0, pMean, bu;

  /* Check time */

  CheckValue(FUNCTION_NAME, "t2", "", t, ZERO, INFTY);
  
  /* Form the actual polynomial coefficients for total power */
  
  c2 =  RDB[DATA_BURN_FIT_C2W1]*RDB[DATA_BURN_POW_PS1]
       +RDB[DATA_BURN_FIT_C2W2]*RDB[DATA_BURN_POW_BOS]
       +RDB[DATA_BURN_FIT_C2W3]*RDB[DATA_BURN_POW_EOS];

  c1 =  RDB[DATA_BURN_FIT_C1W1]*RDB[DATA_BURN_POW_PS1]
       +RDB[DATA_BURN_FIT_C1W2]*RDB[DATA_BURN_POW_BOS]
       +RDB[DATA_BURN_FIT_C1W3]*RDB[DATA_BURN_POW_EOS];

  c0 =  RDB[DATA_BURN_FIT_C0W1]*RDB[DATA_BURN_POW_PS1]
       +RDB[DATA_BURN_FIT_C0W2]*RDB[DATA_BURN_POW_BOS]
       +RDB[DATA_BURN_FIT_C0W3]*RDB[DATA_BURN_POW_EOS];

  CheckValue(FUNCTION_NAME, "c0", "", c0, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c1", "", c1, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c2", "", c2, -INFTY, INFTY);

  /* Integral of general 2nd order polynomial from 0 (BOS) to t (EOS) */
  /* divided by the step length t */

  pMean = c2/3*t*t + c1/2*t + c0;  /* [W] */
  CheckValue(FUNCTION_NAME, "pMean", "", pMean, 0.0, INFTY);

  /* The burnup corresponding to pMean */

  if (RDB[DATA_INI_BURN_FMASS] > 0.0)
    bu = t/86400*pMean/1.0E6/RDB[DATA_INI_BURN_FMASS]; /* [MWd/kgHM] */
  else
    bu = 0.0; /* Burnup is undefined if there is no fissile material */

  CheckValue(FUNCTION_NAME, "bu2", "", bu, 0.0, INFTY);

  /*
  fprintf(outp,"c0=%.4e, c1=%.4e c2=%.4e\n",c0,c1,c2);
  fprintf(outp,"powers: ps1=%.4e, bos=%.4e, eos=%.4e\n",RDB[DATA_BURN_POW_PS1],
          RDB[DATA_BURN_POW_BOS],RDB[DATA_BURN_POW_EOS]);
  fprintf(outp,"   C0W1=%.4e, C0W2=%.4e, C0W3=%.4e\n",RDB[DATA_BURN_FIT_C0W1],
          RDB[DATA_BURN_FIT_C0W2],RDB[DATA_BURN_FIT_C0W3]);
  fprintf(outp,"   C1W1=%.4e, C1W2=%.4e, C1W3=%.4e\n",RDB[DATA_BURN_FIT_C1W1],
          RDB[DATA_BURN_FIT_C1W2],RDB[DATA_BURN_FIT_C1W3]);
  fprintf(outp,"   C2W1=%.4e, C2W2=%.4e, C2W3=%.4e\n",RDB[DATA_BURN_FIT_C2W1],
          RDB[DATA_BURN_FIT_C2W2],RDB[DATA_BURN_FIT_C2W3]);
  fprintf(outp,"fmass=%.4e, bu=%.4e\n",RDB[DATA_INI_BURN_FMASS],bu); 
  */

  return bu;
}

/******************************************************************************/
/* Calculate the  mean power and step length in seconds when step length is   */
/* known in units of burnup. In: step length [MWd/kgHM], Out: step length [s] */
/*                                                                            */
/* TODO: This has worked without problems, but there must be a better way.    */
/******************************************************************************/

double Calculate_t_from_bu(double bu)
{
  double c2, c1, c0, pMean, t, pMeanOld, tOld;
  
  long i;

  long MAXITER = 1000;  /* Max number of iterations, should only need a few*/
  double RTOL = 1.0e-6; /* Relative convergence limit for powr and step length*/

  /* Check that burnup is well enough defined */

  if (RDB[DATA_INI_BURN_FMASS] <= 0.0)
    Error(0, "No fissile burnable material, burnup step length undefined");

  if (RDB[DATA_BURN_POW_BOS] <= 0.0)
    Die(FUNCTION_NAME,
        "Zero power, step length in burnup ill defined.\n");

  /* Form the polynomial coefficients for total power */
  
  c2 =  RDB[DATA_BURN_FIT_C2W1]*RDB[DATA_BURN_POW_PS1]
       +RDB[DATA_BURN_FIT_C2W2]*RDB[DATA_BURN_POW_BOS]
       +RDB[DATA_BURN_FIT_C2W3]*RDB[DATA_BURN_POW_EOS];

  c1 =  RDB[DATA_BURN_FIT_C1W1]*RDB[DATA_BURN_POW_PS1]
       +RDB[DATA_BURN_FIT_C1W2]*RDB[DATA_BURN_POW_BOS]
       +RDB[DATA_BURN_FIT_C1W3]*RDB[DATA_BURN_POW_EOS];

  c0 =  RDB[DATA_BURN_FIT_C0W1]*RDB[DATA_BURN_POW_PS1]
       +RDB[DATA_BURN_FIT_C0W2]*RDB[DATA_BURN_POW_BOS]
       +RDB[DATA_BURN_FIT_C0W3]*RDB[DATA_BURN_POW_EOS];
  

  /*
  fprintf(outp,"c0=%.4e, c1=%.4e c2=%.4e\n",c0,c1,c2);
  fprintf(outp,"powers: ps1=%.4e, bos=%.4e, eos=%.4e\n",RDB[DATA_BURN_POW_PS1],
          RDB[DATA_BURN_POW_BOS],RDB[DATA_BURN_POW_EOS]);
  fprintf(outp,"   C0W1=%.4e, C0W2=%.4e, C0W3=%.4e\n",RDB[DATA_BURN_FIT_C0W1],
          RDB[DATA_BURN_FIT_C0W2],RDB[DATA_BURN_FIT_C0W3]);
  fprintf(outp,"   C1W1=%.4e, C1W2=%.4e, C1W3=%.4e\n",RDB[DATA_BURN_FIT_C1W1],
          RDB[DATA_BURN_FIT_C1W2],RDB[DATA_BURN_FIT_C1W3]);
  fprintf(outp,"   C2W1=%.4e, C2W2=%.4e, C2W3=%.4e\n",RDB[DATA_BURN_FIT_C2W1],
          RDB[DATA_BURN_FIT_C2W2],RDB[DATA_BURN_FIT_C2W3]);
  fprintf(outp,"fmass=%.4e, bu=%.4e\n",RDB[DATA_INI_BURN_FMASS],bu); 
  */

  /***************************************************************************/
  /* neither step length nor total power is known, so they need to be solved */
  /* together. This leads to the pair:                                       */
  /*    pMean = c2/3*t*t + c1/2*t + c0                                       */
  /*    t = bu/pMean                                                         */
  /* These yield a third order equations which could be solved analytically, */
  /* but the solution is quite nasty, so we iterate instead.                 */
  /*                                                                         */
  /***************************************************************************/

  /* initial quess is to use BOS power for whole step */

  pMean = RDB[DATA_BURN_POW_BOS];            /* [W] */
  t = 86400E+6*bu*RDB[DATA_INI_BURN_FMASS]/pMean; /* [s] */
  
  /* Simple fixed point iterateion */
  
  for (i=0; i< MAXITER; i++)
    {
      /* store old values */
      pMeanOld = pMean;
      tOld    = t;
      
      /* update values */
      pMean = c2/3*t*t + c1/2*t + c0;             /* [W] */
      t = 86400E+6*bu*RDB[DATA_INI_BURN_FMASS]/pMean; /* [s] */

      /* check for convergence */
      if (fabs(pMeanOld/pMean-1)<RTOL && fabs(tOld/t-1)<RTOL)
        break;
    }
  
  /* As far as I have seen the iteration always converges (and fast),
     but who knows. */

  if (i >= MAXITER)
      Die(FUNCTION_NAME, "Step length iteration did not converge."); 

  return t;

}
