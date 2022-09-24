/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : depletionpolyfit.c                             */
/*                                                                           */
/* Created:       2011/06/07 (AIs)                                           */
/* Last modified: 2018/05/04 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Makes a (low order) polynomial approximation of the time     */
/*              development of cross-sections ect. in burnup calculations.   */
/*              Extrapolation on the predictor, interpolation on corrector.  */
/*                                                                           */
/* Comments: -This could probably be merged with SetDepStepSize to reduce    */
/*            the number of source files.                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DepletionPolyFit:"

/*****************************************************************************/

void DepletionPolyFit(long dep, long step)
{
  double t1,t2,t3, c0w1,c0w2,c0w3, c1w1,c1w2,c1w3, c2w1,c2w2,c2w3;
  long order;

  /* Avoid compiler warnings */

  c2w1 = 0.0;
  c2w2 = 0.0;
  c2w3 = 0.0;
  c1w1 = 0.0;
  c1w2 = 0.0;
  c1w3 = 0.0;
  c0w1 = 0.0;
  c0w2 = 0.0;
  c0w3 = 0.0;
  order = -666;

  /* Select "effective" type (i.e., due to lack of previous step data)*/

  if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    {
      if ((long)RDB[DATA_BURN_PRED_TYPE] == PRED_TYPE_CONSTANT)
        {
          /* Constant extrapolation*/

          order = 0;
        }
      else if (((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_DEC_STEP) ||
               ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_DEC_TOT) ||
               ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_ACT_STEP) ||
               ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_ACT_TOT))
        {
          /* Decay or activation step, CE is sufficient */

          order = 0;
        }
      else if ((Mean((long)RDB[RES_TOT_NEUTRON_FLUX], 0) <= 0.0))
        {
          /* Zero flux, use CE */

          order = 0;
        }
      else if ((long)RDB[DATA_BURN_PRED_TYPE] == PRED_TYPE_LINEAR)
        {
          /* Linear extrapolation */

          if (step == 0)
            {
              /* First step, no previous step data => use CE */

              order = 0;
            }
          else
            order = -1;
        }
      else
        Die(FUNCTION_NAME, "Unknown predictor type");
    }
  else if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    {
      /* Decay only steps are supposed to be done with simple Euler's method */

      if (((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_DEC_STEP) ||
          ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_DEC_TOT))
        Die(FUNCTION_NAME,"Corrector step on decay interval");

      if (((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_ACT_STEP) ||
          ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_ACT_TOT))
        Die(FUNCTION_NAME,"Corrector step on activation interval");

      if ((Mean( (long)RDB[RES_TOT_NEUTRON_FLUX], 0) <= 0.0))
        Die(FUNCTION_NAME,"Corrector step with zero flux");

      /* Check corrector type */

      if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_CONSTANT)
        {
          /* Constant backwards extrapolation */

          order = 100; /* change numbering if the method stays */
        }
      else if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_LINEAR)
        {
          /* Linear interpolation */

          order = 1;
        }
      else if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_QUADRATIC)
        {
          if (step < 2)
            {
              /* First and second step use linear instead */

              order = 1;
            }
          else
            {
              /* quadratic interpolation */

              order = 2;
            }
        }
      else if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE)
        Die(FUNCTION_NAME,"Corrector step of type NONE");
      else
        Die(FUNCTION_NAME, "Unknown corrector type");
    }
  else
    Die(FUNCTION_NAME, "Invalid PC mode");

  /***************************************************************************/
  /* NEW (AIs)                                                               */
  /* In inner iteration mode the corrector if first iterated with backwards  */
  /* extrapolation of the EOS values, only final round uses the selected     */
  /* scheme                                                                  */
  /*                                                                         */
  /***************************************************************************/

  if (((long)RDB[DATA_BURN_CI_TYPE] == CI_TYPE_INNER) &&
      ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
      ((long)RDB[DATA_BURN_CI_LAST] == NO))
    order = 100;

  /***************************************************************************/
  /* With Stochastic Implicit Euler constant backwards extrapolation is used */
  /***************************************************************************/

  if ((long)RDB[DATA_BURN_SIE] == YES)
    {
      if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
        order = 100;
      else
        order = 0;
    }

  /***************************************************************************/
  /*                                                                         */
  /* Set times corresponding to the PS1 (t1), BOS (t2) and EOS (t3)          */
  /* The EOS refers to the end of preceding predictor step, and it is only   */
  /* used at corrector steps                                                 */
  /* t=0 is set to BOS, i.e., t2=0 by definition/convention                  */
  /*                                                                         */
  /***************************************************************************/

  /* Time corresponding to BOS of previous step */

  if (step > 0)
    {
      t1 = -RDB[DATA_BURN_PS1_LENGTH];
      CheckValue(FUNCTION_NAME, "t1", "", t1, -INFTY, ZERO);
    }
  else
    t1 = -INFTY;

  t2 = 0.0; /* not actually used */

  /* Avoid compiler warning */

  t2 = t2 + 0.0;

  /* If on corrector, get predictor length (EOS time in the fit) it differs */
  /* from the corrector step length in some normalization modes */

  if ( (long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP )
    {
      t3 = RDB[DATA_BURN_PRED_LENGTH];
      CheckValue(FUNCTION_NAME,"t3","",t3,ZERO,INFTY);
    }
  else
    t3 = -INFTY;

  /***************************************************************************/
  /*                                                                         */
  /* Fit polynomial through the PS1, BOS and EOS (from predictor) points     */
  /*                                                                         */
  /* polynomial fit is defined by its coefficients which depend on the       */
  /* PS1 (n=1), BOS (n=2) and EOS (n=3) values. c<order>w<n> are the weights */
  /* of these values in the zeroth, first and second order coefficients.     */
  /*                                                                         */
  /* I.e., a=c2w1*y1+c2w2*y2+c2w3*y3 would be the second order term of the   */
  /* polynomial fitted to points (t1,y1) (t2,y2) and (t3,y3)                 */
  /* y1, y2 and y3 are the PS1, BOS and EOS values of the quality being      */
  /* interpolated/extrapolated (e.g., some cross-section)                    */
  /*                                                                         */
  /***************************************************************************/

  if (((long)RDB[DATA_OTF_BURN_MODE] == YES) &&
      ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
    {
      /* Corrector step for On-the-fly mode */
      printf("%s %E\n", FUNCTION_NAME, t3);
      c2w1 = 0.0;
      c2w2 = 0.0;
      c2w3 = 0.0;
      c1w1 = 0.0;
      c1w2 = 0.0;
      c1w3 = 2.0/t3;
      c0w1 = 0.0;
      c0w2 = 0.0;
      c0w3 = 0.0;
    }
  else if (order == 0)
    {
      /* Constant extrapolation */

      c2w1 = 0.0;
      c2w2 = 0.0;
      c2w3 = 0.0;
      c1w1 = 0.0;
      c1w2 = 0.0;
      c1w3 = 0.0;
      c0w1 = 0.0;
      c0w2 = 1.0;
      c0w3 = 0.0;
    }
  else if (order == -1)
    {
      /* Linear extrapolation */

      c2w1 = 0.0;
      c2w2 = 0.0;
      c2w3 = 0.0;
      c1w1 = 1.0/t1;
      c1w2 = -1.0/t1;
      c1w3 = 0.0;
      c0w1 = 0.0;
      c0w2 = 1.0;
      c0w3 = 0.0;
    }
  else if (order == 100)
    {
      /* Constant backwards extrapolation */

      c2w1 = 0.0;
      c2w2 = 0.0;
      c2w3 = 0.0;
      c1w1 = 0.0;
      c1w2 = 0.0;
      c1w3 = 0.0;
      c0w1 = 0.0;
      c0w2 = 0.0;
      c0w3 = 1.0;

    }
  else if (order == 1)
    {
      /* Linear interpolataion */

      c2w1 = 0.0;
      c2w2 = 0.0;
      c2w3 = 0.0;
      c1w1 = 0.0;
      c1w2 = -1.0/t3;
      c1w3 = 1.0/t3;
      c0w1 = 0.0;
      c0w2 = 1.0;
      c0w3 = 0.0;
    }
  else if (order==2)
    {
      /* Quadratic interpolation, second order term */

      c2w1 = (   - t3/t1 )/(-t1*t3 + t3*t3);
      c2w2 = (-1 + t3/t1 )/(-t1*t3 + t3*t3);
      c2w3 = ( 1         )/(-t1*t3 + t3*t3);

      /* First order term */

      c1w1 =       - c2w1*t3;
      c1w2 = -1/t3 - c2w2*t3;
      c1w3 =  1/t3 - c2w3*t3;

      /* Constant term */

      c0w1 = 0.0;
      c0w2 = 1.0;
      c0w3 = 0.0;
    }
  else
    Die(FUNCTION_NAME, "Invalid effective order");

  /* Check values */

  CheckValue(FUNCTION_NAME, "c0w1", "", c0w1, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c0w2", "", c0w2, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c0w3", "", c0w3, -INFTY, INFTY);

  CheckValue(FUNCTION_NAME, "c1w1", "", c1w1, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c1w2", "", c1w2, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c1w3", "", c1w3, -INFTY, INFTY);

  CheckValue(FUNCTION_NAME, "c2w1", "", c2w1, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c2w2", "", c2w2, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "c2w3", "", c2w3, -INFTY, INFTY);

  /* Save the coefficients */

  WDB[DATA_BURN_FIT_C2W1] = c2w1;
  WDB[DATA_BURN_FIT_C2W2] = c2w2;
  WDB[DATA_BURN_FIT_C2W3] = c2w3;

  WDB[DATA_BURN_FIT_C1W1] = c1w1;
  WDB[DATA_BURN_FIT_C1W2] = c1w2;
  WDB[DATA_BURN_FIT_C1W3] = c1w3;

  WDB[DATA_BURN_FIT_C0W1] = c0w1;
  WDB[DATA_BURN_FIT_C0W2] = c0w2;
  WDB[DATA_BURN_FIT_C0W3] = c0w3;

  /* store the type of fit */

  if (order == 0)
    WDB[DATA_BURN_FIT_TYPE] = (double)PRED_TYPE_CONSTANT;
  else if (order == -1)
    WDB[DATA_BURN_FIT_TYPE] = (double)PRED_TYPE_LINEAR;
  else if (order == 100)
    WDB[DATA_BURN_FIT_TYPE] = (double)CORR_TYPE_CONSTANT;
  else if (order == 1)
    WDB[DATA_BURN_FIT_TYPE] = (double)CORR_TYPE_LINEAR;
  else if (order == 2)
    WDB[DATA_BURN_FIT_TYPE] = (double)CORR_TYPE_QUADRATIC;
  else
    Die(FUNCTION_NAME, "wtf?");
}

/*****************************************************************************/
