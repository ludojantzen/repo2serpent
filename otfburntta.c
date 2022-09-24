/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : otfburntta.c                                   */
/*                                                                           */
/* Created:       2018/04/04 (JLe)                                           */
/* Last modified: 2018/04/04 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: TTA solver for on-the-fly burnup routine                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "OTFBurnTTA:"

void OTFBurnTTALoop (long iso, double t, long n, double *l,
                     double B, double Pin, double adens);

double OTFBurnSimpleTTAChain(long, double, double *, double *);

#define MAX_TTA_CHAIN 500

/*****************************************************************************/

void OTFBurnTTA(long mat)
{
  long loc0, dep, type, ptr, n;
  double t, l[MAX_TTA_CHAIN], adens;

  /* Check depletion step type */

  if (((long)RDB[DATA_BURN_CORR_TYPE] != CORR_TYPE_NONE) &&
      ((long)RDB[DATA_BURN_STEP_PC] != CORRECTOR_STEP))
    Die(FUNCTION_NAME, "Predictor step");

  /* Avoid compiler warning */
      
  t = -1.0;
  
  /* Check step type */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    {
      /* Corrector step, use value calculated for predictor */

      t = 0.5*RDB[DATA_BURN_TIME_INTERVAL];
    }
  else
    {
      /* Get pointer to depletion interval */
      
      dep = (long)RDB[DATA_BURN_PTR_CURRENT_DEP];
      CheckPointer(FUNCTION_NAME, "(dep)", DATA_ARRAY, dep);
  
      /* Get type */

      type = (long)RDB[dep + DEP_HIS_STEP_TYPE];
      
      /* Pointer to steps */
      
      ptr = (long)RDB[dep + DEP_HIS_PTR_STEPS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Get current step */
      
      n = (long)RDB[DATA_BURN_PTR_CURRENT_STEP];
      
      /* Check if last */
      
      if ((NextItem(dep) < VALID_PTR) && 
          (n == (long)RDB[dep + DEP_HIS_N_STEPS]))
        {
          /* Loop over data */
          
          loc0 = (long)RDB[mat + MATERIAL_PTR_OTF_BURN];
          while (loc0 > VALID_PTR)
            {
              /* Check material pointer */
              
              if ((long)RDB[loc0 + OTF_BURN_PTR_MAT] != mat)
                break;
              
              /* Set atomic density BOS value */
              
              WDB[loc0 + OTF_BURN_ADENS] = RDB[loc0 + OTF_BURN_ADENS0];
              
              /* Pointer to next */
              
              loc0 = NextItem(loc0);
            }
          
          /* Exit subroutine */
          
          return;
        }
      
      /* Check type and get time interval */
      
      if ((type == DEP_STEP_DAY_STEP) || (type == DEP_STEP_ACT_STEP))
        t = RDB[ptr + n]*60.0*60.0*24.0;
      else if ((type == DEP_STEP_DAY_TOT) || (type == DEP_STEP_ACT_TOT))
        t = RDB[ptr + n]*60.0*60.0*24.0 - RDB[DATA_BURN_CUM_BURNTIME];
      else
        Die(FUNCTION_NAME, "Burnup units not supported");

      /* Take mid-point */
      
      t = 0.5*t;
    }

  /* Loop over data */

  loc0 = (long)RDB[mat + MATERIAL_PTR_OTF_BURN];
  while (loc0 > VALID_PTR)
    {
      /* Check material pointer */

      if ((long)RDB[loc0 + OTF_BURN_PTR_MAT] != mat)
        break;

      /* Get initial atomic density */
      
      adens = RDB[loc0 + OTF_BURN_ADENS0];

      /* Check density and deplete */

      if (adens > 0.0)
        OTFBurnTTALoop(loc0, t, 0, l, 1.0, 1.0, adens);

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }
}
  
/*****************************************************************************/

void OTFBurnTTALoop (long loc0, double t, long n, double *l, 
                     double B, double Pin, double adens)
{
  long  loc1, nuc, mat, tgt, ptr, yld, rea, i;
  double lambda, f, br, Pout, norm, vol;

  /* Chain cut-off by remaining fraction */

  if (Pin < 1E-6)
    return;

  /* Get normalization factor */

  norm = NormCoef(PARTICLE_TYPE_NEUTRON);
  CheckValue(FUNCTION_NAME, "norm", "", norm, 0.0, INFTY);
 
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Pointer to material */

  mat = (long)RDB[loc0 + OTF_BURN_PTR_MAT];
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Get volume */

  if ((vol = RDB[mat + MATERIAL_VOLUME]) == 0.0)
    Die(FUNCTION_NAME, "Volume is zero");

  /* Adjust normalization factor */
  
  norm = BARN*norm/vol;

  /* Tarkastetaan ketjun pituus */

  if (n > MAX_TTA_CHAIN - 1)
    Die(FUNCTION_NAME, "Maximum length of TTA chain exceeded");

  /* Total transmutation cross section */ 

  ptr = (long)RDB[loc0 + OTF_BURN_PTR_TOT_XS];
  lambda = norm*BufVal(ptr, 0);
  CheckValue(FUNCTION_NAME, "lambda", "", lambda, 0.0, INFTY);

  /* Pointer to nuclide */

  nuc = (long)RDB[loc0 + OTF_BURN_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Add decay constant */

  lambda = lambda + RDB[nuc + NUCLIDE_LAMBDA];

  /* Check if stable */
  
  if (lambda == 0.0)
    {
      /* Add to atomic density */

      WDB[loc0 + OTF_BURN_ADENS] = RDB[loc0 + OTF_BURN_ADENS] + Pin*adens;

      /* Exit subroutine */

      return;
    }

  /* Check for closed loop */
  
  for(i = 0; i < n; i++)
    if (fabs((lambda - l[i]) / l[i]) < 1E-8)
      {
        /* Check remaining fraction */

        if (Pin > 1E-4)
          Warn(FUNCTION_NAME, "Truncation error from closed loop: %E", Pin);
        
        /* Exit subroutine */
        
        return;
      }

  /* Put coeficients */

  l[n] = lambda;
  n++;

  /* Solve chain */

  f = OTFBurnSimpleTTAChain(n, t, l, &Pout);
  Pout = Pout*B;

  /* Check value for failed numerical precision */

  if ((f > 1.0) || (f < 0.0) || (Pout < 0.0) || (Pout > Pin*1.001))
    {    
      Warn(FUNCTION_NAME, "TTA failed %E %E %E", f, Pout, Pin);
      return;
    }

  /* Add to atomic density */

  WDB[loc0 + OTF_BURN_ADENS] = RDB[loc0 + OTF_BURN_ADENS] + B*f*adens;

  /* Loop over transmutation list */

  loc1 = (long)RDB[loc0 + OTF_BURN_PTR_REA];
  while (loc1 > VALID_PTR)
    {
      /* Pointer to reaction data */

      rea = (long)RDB[loc1 + OTF_BURN_REA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      
      /* Pointer to target */

      if ((tgt = (long)RDB[loc1 + OTF_BURN_REA_PTR_TGT]) > VALID_PTR)
        {
          /*******************************************************************/

          /***** Transmutation reaction **************************************/

          /* Avoid compiler warning */

          br = -1.0;

          /* Calculate branching ratio */

          if ((ptr = (long)RDB[loc1 + OTF_BURN_REA_PTR_XS]) > VALID_PTR)
            br = norm*BufVal(ptr, 0)/lambda;
          else if ((f = RDB[nuc + NUCLIDE_LAMBDA]) > 0.0)
            br = f*RDB[rea + REACTION_BR]/lambda;
          else
            Die(FUNCTION_NAME, "WTF?");

          /* Check */

          CheckValue(FUNCTION_NAME, "br", "", br, 0.0, 1.0);

          /* Step forward */
                  
          OTFBurnTTALoop(tgt, t, n, l, B*br, Pout*br, adens);

          /*******************************************************************/
        }
      else if ((yld = (long)RDB[loc1 + OTF_BURN_REA_PTR_FISSY]) > VALID_PTR)
        {
          /*******************************************************************/

          /***** Fission *****************************************************/

          /* Avoid compiler warning */

          br = -1.0;

          /* Calculate branching ratio */

          if ((ptr = (long)RDB[loc1 + OTF_BURN_REA_PTR_XS]) > VALID_PTR)
            br = norm*BufVal(ptr, 0)/lambda;
          else if ((f = RDB[nuc + NUCLIDE_LAMBDA]) > 0.0)
            br = f*RDB[rea + REACTION_BR]/lambda;
          else
            Die(FUNCTION_NAME, "WTF?");

          /* Check */

          CheckValue(FUNCTION_NAME, "br", "", br, 0.0, 1.0);

          /* Loop over products */

          while (yld > VALID_PTR)
            {
              /* Pointer to target */

              tgt = (long)RDB[yld + OTF_BURN_FISSY_PTR_TGT];
              CheckPointer(FUNCTION_NAME, "(tgt)", DATA_ARRAY, tgt);

              /* Pointer to yield data */

              ptr = (long)RDB[yld + OTF_BURN_FISSY_PTR_YLD];
              CheckPointer(FUNCTION_NAME, "(tgt)", DATA_ARRAY, tgt);

              /* Get fraction */

              f = RDB[ptr + FY_INDEPENDENT_FRAC];
              CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 0.1);

              /* Step forward */

              OTFBurnTTALoop(tgt, t, n, l, B*br*f, Pout*br*f, adens);

              /* Next */

              yld = NextItem(yld);
            }

          /*******************************************************************/
        }

      /* Next reaction */

      loc1 = NextItem(loc1);
    }

  /* Exit subroutine */
  
  return;
}

/*****************************************************************************/

/***** Solution to not accounting for loops **********************************/

double OTFBurnSimpleTTAChain(long n, double t, double *l, double *Pout)
{
  long i, j;
  double a, b, sum;
  double sumP;

  /* Calculate coefficient for denominator */
  
  b = 1.0;
  for (i = 0; i < n-1; i++)
    b = b*l[i];

  /* Calculate sum */

  sum  = 0.0;
  sumP = 0.0;

  for (i = 0; i < n; i++)
    {
      /* Check lambda */

      CheckValue(FUNCTION_NAME, "l", "", l[i], ZERO, INFTY);

      /* Denominator of alpha */
      
      a = 1.0;
      
      for (j = 0; j < n; j++)
        if (i != j)
          a = a/(l[j] - l[i]);

      sum  = sum  + a*exp(-l[i]*t);
      sumP = sumP + a*exp(-l[i]*t)/l[i];
    }

  /* Return values */
  
  if (l[n - 1] == 0.0)
    *Pout = 0.0;
  else
    *Pout = 1 - l[n - 1]*b*sumP;

  return b*sum;
}

/*****************************************************************************/
