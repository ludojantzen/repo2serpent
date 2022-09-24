/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reinitrng.c                                    */
/*                                                                           */
/* Created:       2011/03/03 (TVi)                                           */
/* Last modified: 2015/12/04 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: This routine is used to skip n*2^STRIDE steps forward in the */
/*              random number sequence beginning from parent seed.           */ 
/*                                                                           */
/* Comments: Algorithm taken from MCNP5. Basically calculates                */
/*           g^(n*2^STRIDE)*parentseed+(g^(n*2^STRIDE)-1)/(g-1)*12345mod2^64 */
/*           efficiently.                                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReInitRNG:"

/*****************************************************************************/

unsigned long ReInitRNG(long n0)
{
  unsigned long gen, g, inc, c, gp, n;

  /* Copy to unsigned variable */

  n = n0;
  
  /* Adjust index in MPI mode */
  
  if ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO)
    {
      /* Check maximum (should not be a problem, because the input */
      /* value is related to particle history index, etc.) */

      if (n > (ULONG_MAX - mpiid)/mpitasks)
        Die(FUNCTION_NAME, "Maximum rng sequence reached");
      else
        n = n0*mpitasks + mpiid;
    }

  /* Re-initialize RNG */

  n = n << STRIDE;

  gen = 1;
  g = 2862933555777941757;
  inc = 0;
  c = 12345;

  while(n > 0)
    {
      if((n%2) == 1)
        {
          gen *= g;
          inc = inc*g + c;
        }
      
      gp = g + 1;
      g *= g;
      c *= gp;
      
      n = n>>1;
    }  

  return parent_seed*gen + inc;
}

/*****************************************************************************/
