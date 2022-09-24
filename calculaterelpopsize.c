/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculaterelpopsize.c                          */
/*                                                                           */
/* Created:       2014/10/07 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Calculates next iter population size for coupled calculation */
/*              with solution relaxation                                     */
/*                                                                           */
/* Comments: - TODO: Toi batchikoon kasvattaminen ei toimi microgroupXS:n    */
/*             kanssa kun sen batchikokoa pitää myös kasvattaa. Ei ehkä ole  */
/*             hankala homma.                                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateRelPopSize:"

#define MAX_ITER 100

#define nobatch

/*****************************************************************************/

void CalculateRelPopSize()
{
  long loc0, nb;
  double bleft;
  double popleft;

#ifndef nobatch
  long loc1, ptr, ptr1, nz, nr, i, j, k, l, na, nt, type, idx, ptr2, iter;
  double alpha, ncur, ntot, pn, prod, upper, lower, vn1, x0, fx, dfx, sn1;

  /* Avoid compiler error */

  upper = 0.0;
  lower = 0.0;
#endif

  /* Check that some interfaces are defined */

  if ((loc0 = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  /************************************************************/
  /* Calculate pn for adaptive population size                */
  /* See Dufek & Gudowski NSE vol 152, pp. 274-183 (2006)     */


  /* Loop over interfaces */

#ifndef nobatch

  /* Get relaxation alpha from memory */

  alpha = RDB[DATA_SOL_REL_ALPHA];

  while(loc0 > VALID_PTR)
    {
      /* Get interface type */

      type = (long)RDB[loc0 + IFC_TYPE];

      /* Calculation of pn currently only works with fuel behavior */
      /* interface                                                 */

      if((type == IFC_TYPE_FUEP) || (type == IFC_TYPE_FPIP))
        {

          /* Get pointer to first pin */

          loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];

          /* Reset numerator and denominator*/

          upper = 0.0;
          lower = 0.0;

          /* Loop over pins to calculate pn*/

          while(loc1 > VALID_PTR)
            {

              /* Get pointer to output limits*/

              ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];
              CheckPointer(FUNCTION_NAME, "(limptr)", DATA_ARRAY, ptr);

              /* get number of bins*/

              nz = (long)RDB[ptr + FUEP_NZ];

              na = (long)RDB[ptr + FUEP_NA];

              nr = (long)RDB[ptr + FUEP_NR];

              nt = (long)RDB[ptr + FUEP_NT];

              /* Get pointer to iteration power */

              ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER];
              CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);

              /* Get pointer to relaxed power */

              ptr1 = (long)RDB[loc1 + IFC_FUEP_PTR_POWER_REL];
              CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr1);

              /* Get pointer to gradient of relaxed power */

              ptr2 = (long)RDB[loc1 + IFC_FUEP_PTR_POWER_GRD];
              CheckPointer(FUNCTION_NAME, "(Pptr2)", DATA_ARRAY, ptr2);

              /* g(n+1) = ptr - ptr1 g(n) = same from prev. relaxation */
              /* See Eq. 25 for calculation of gradient */

              for(i=0; i < nz; i++)
                {
                  for(j=0; j < na; j++)
                    {
                      for(k=0; k < nr; k++)
                        {
                          for(l=0; l < nt; l++)
                            {

                              /* Calculate bin index */

                              idx = i + j*nz + k*nz*na + l*nz*na*nr;

                              /* Calculate product of gradients of this iteration */
                              /* and prev. iteration */

                              prod = (Mean(ptr,i,j,k,l) - RDB[ptr1 + idx])*
                                RDB[ptr2 + idx];

                              /* Store new gradient */

                              WDB[ptr2 + idx] = Mean(ptr,i,j,k,l) - RDB[ptr1 + idx];

                              /* Add product to numerator or denominator */
                              /* based on sign */

                              if(prod >= 0)
                                {
                                  upper = upper + prod;
                                  /*lower = lower + 0.0;*/
                                }
                              else
                                {
                                  /*upper = upper + 0.0*/
                                  lower = lower - prod;
                                }

                            }

                        }

                    }

                }

              /* Next pin */

              loc1 = NextItem(loc1);
            }

        }
      else if(type == IFC_TYPE_TET_MESH)
        {
          Die(FUNCTION_NAME, "Smart population size increse not implemented for tet mesh");
        }

      /* Next interface */

      loc0 = NextItem(loc0);

    }

  /* Handle zero denominator */

  if(lower == 0)
    {
      if(upper == 0)
        pn = 0;
      else
        pn = 1e3;
    }
  /* Calculate scalar pn */
  else
    pn = log(upper/lower);

  /* Decide step size control based on pn */

  if(pn > 0)
    vn1 = -0.3;
  else if(pn < 0)
    vn1 = -0.65;
  else
    vn1 = -0.5;

  /* Get current and total population sizes */

  ncur = RDB[DATA_SOL_REL_NCUR];
  ntot = RDB[DATA_SOL_REL_NTOT];

  /* Solve Sn+1 with Newtons method */

  /* Initial guess is to double the current population */

  x0 = ntot + 2*ncur;

  /* Initial value for the function */

  fx = (x0 - ntot)/x0 - alpha*pow((x0/ntot),vn1);

  /* Initial value for the derivative */

  dfx = ntot/(x0*x0) - alpha*vn1/pow(ntot,vn1)*pow(x0,vn1-1.0);

  /* Reset iteration counter */

  iter = 0;

  /* Iteration loop */

  while(fabs(fx) > 1e-4)
    {

      if(iter > MAX_ITER)
        {
          Warn(FUNCTION_NAME, "Newton's iteration reached MAX_ITER %ld", (long)MAX_ITER);
          break;
        }

      /* Calculate new guess for the root */

      x0 = x0 - fx/dfx;

      /* Update function and derivative values */

      fx = (x0 - ntot)/x0 - alpha*pow((x0/ntot),vn1);
      dfx = ntot/(x0*x0) - alpha*vn1/pow(ntot,vn1)*pow(x0,vn1-1.0);

      /* Increment iteration number */

      iter++;
    }

  /*
  fprintf(outp, "iter = %ld, fx = %E, pn = %E (%E/%E), vn1 = %E, Sn = %ld, Sn+1 = %ld, sn+1 = %ld, s1 = %ld, alpha = %E\n", iter, fx, pn,upper,lower, vn1, (long)ntot, (long)x0, (long)(x0-ntot), (long)(alpha*alpha*ntot), alpha);
  fprintf(outp, "%.6E %ld %.6E %% Res\n",pn, (long)ntot, alpha);
  */

  /* Calculate number of neutrons on next iteration */

  sn1 = (long)(x0-ntot);

#endif

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Calculate the population that will be still simulated */

      popleft = RDB[DATA_SOL_REL_MAX_POP] - RDB[DATA_SOL_REL_NTOT];

      /* Calculate the number of batches that the population will fill */

      bleft = popleft/RDB[DATA_CRIT_POP];

      /* If we cannot simulate any more batches */

      if (bleft < 1)
        {

          /* Set next iteration size to zero */

          WDB[DATA_SOL_REL_NCUR] = 0.0;

          WDB[DATA_CRIT_CYCLES] = (double)0;

          /* Disable iteration flag */

          WDB[DATA_ITERATE] = (double)NO;

          return;

        }
    }

#ifndef nobatch
/* Get next sample size */

  if(sn1 < popleft)
    ncur = sn1;
  else
    ncur = popleft;

  /* Round the batch number */

  nb = round(ncur/RDB[DATA_CRIT_POP]);

#else

  /* Non-increasing batch size */

  nb = (long)RDB[DATA_CRIT_CYCLES];

#endif

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Store number of neutrons to be simulated in next iteration */

      WDB[DATA_SOL_REL_NCUR] = ((double)nb)*RDB[DATA_CRIT_POP];

      /* Store number of criticality cycles */

      WDB[DATA_CRIT_CYCLES] = (double)nb;
    }
  else
    {
      /* Dynamic simulation mode */

      WDB[DATA_SOL_REL_NCUR] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

    }

  /*
  fprintf(outp, "New batch interval %ld\n", nb);
  */

  return;

  /***************************************************************************/
}

/*****************************************************************************/
