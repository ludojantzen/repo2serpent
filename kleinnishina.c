/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : kleinnishina.c                                 */
/*                                                                           */
/* Created:       2011/04/15 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Samples photon momentum transfer from the Klein-Nishina      */
/*              formula                                                      */
/*                                                                           */
/* Comments: Uses method described in Appendix 3A of reference:              */
/*                                                                           */
/*           I. Lux and L. Koblinger. "Monte Carlo particle transport        */
/*           methods: neutron and photon calculations" CRC Press 1991.       */
/*                                                                           */
/*           (toi jälkimmäinen rutiini on oikeasti suoraan MCNP5:sta)        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "KleinNishina:"

/*****************************************************************************/

void KleinNishina(double Ein, double *Eout, double *mu, long id)
{
  double alpha0, alpha, rnd1, rnd2, rnd3, x, y, a, b, t;

  /* Photon energy relative to neutron rest mass */

  alpha0 = Ein/E_RESTMASS;

  /* Compare to limit 1 + sqrt(3) */

  if (alpha0 < 2.73205080756888)
    {
      /***** Use Kahn's rejection method (Fig. 3A.1, pp. 68) *****************/

      /* Rejection loop */

      while (1 != 2)
	{
	  /* Sample three random numbers */

	  rnd1 = RandF(id);
	  rnd2 = RandF(id);
	  rnd3 = RandF(id);

	  /* Branch */

	  if (rnd1 > (1.0 + 2.0*alpha0)/(9.0 + 2.0*alpha0))
	    {
	      /* Calculate x */

	      x = (1.0 + 2.0*alpha0)/(1.0 + 2.0*rnd2*alpha0);

	      /* Calculate coefficient */
	      
	      y = (1.0/alpha0 - x/alpha0 + 1.0);

	      /* Rejection test */

	      if (2.0*rnd3 <= y*y + 1.0/x)
		break;
	    }
	  else
	    {
	      /* Calculate x */

	      x = 1.0 + 2.0*rnd2*alpha0;

	      /* Rejection test */

	      if (rnd3 <= 4.0*(1.0/x - 1.0/x/x))
		break;
	    }
	}

      /* Calculate alpha */

      alpha = alpha0/x;

      /***********************************************************************/
    }
  else
    {
      /***** Use Koblinger's method (pp.66) **********************************/

      /* Calculate two coefficients */

      a = 1.0/alpha0;
      b = 1.0/(1.0 + 2.0*alpha0);

      /* Sample distribution function */

      t = 4.0*a + 0.5*(1.0 - b*b) - (1.0 - 2.0*(1.0 + alpha0)*a*a)*log(b);
      t = RandF(id)*t;

      /* Pick interval */
      
      if (t <= 2.0*a)
	alpha = 1.0/(a + 2.0*RandF(id));
      else if (t <= 4.0*a)
	alpha = alpha0*(1.0 + RandF(id)*(b - 1.0));
      else if(t <= 4.0*a + 0.5*(1.0-b*b))
	alpha = alpha0*sqrt(1.0 - RandF(id)*(1.0 - b*b));
      else
	alpha = alpha0*pow(b, RandF(id));
    }

  /* Calculate new energy */

  *Eout = alpha*E_RESTMASS;
  
  /* Calculate direction cosine */
  
  *mu = 1.0 + 1.0/alpha0 - 1.0/alpha;
}

/*****************************************************************************/
