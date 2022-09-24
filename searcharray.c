/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : searcharray.c                                  */
/*                                                                           */
/* Created:       2010/12/09 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Finds interval from array of doubles                         */
/*                                                                           */
/* Comments: - Based on binary search                                        */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "SearchArray:"

/*****************************************************************************/

long SearchArray(const double *dat, double val, long N)
{
  long lo, hi, n;

  /* Check that values are in ascending order */

#ifdef DEBUG
  /*
  for (n = 0; n < N - 1; n++)
    if (dat[n] > dat[n + 1])
      Die(FUNCTION_NAME, "Array is not sorted");
  */
#endif

  /* Check boundaries */
      
  if ((val < dat[0]) || (val >= dat[N - 1]))
    return -1;
  else if (val == dat[0])
    return 0;

  /* Initial boundaries */

  lo = 0;
  hi = N - 1;
      
  /* Search loop */

  while(1 != 2)
    {
      /* Check number of remaining points */

      if (hi - lo < 20)
	{
	  /***** Loop to interval ********************************************/

	  /* Check mid-point */

	  if (val < (dat[hi - 1] + dat[lo])/2.0)
	    {
	      /* Loop to interval from beginning */

	      for (n = lo; n < hi + 1; n++)
		{
		  if (dat[n] == val)
		    return n;
		  else if (dat[n] > val)
		    return n - 1;
		}

	      /* Error */

	      Die(FUNCTION_NAME, "Interval not found in array");
	    }
	  else
	    {
	      /* Loop to interval from end */

	      for (n = hi; n > lo - 1; n--)
		if (dat[n] <= val)
		  return n;

	      /* Error */

	      Die(FUNCTION_NAME, "Interval not found in array");
	    }
      
	  /*******************************************************************/
	}
      else
	{
	  /***** Split interval **********************************************/
	  
	  /* New guess */
	  
	  n = (long)(((double)hi + (double)lo)/2.0);
	  	  
	  /* New boundaries */
	  
	  if (val < dat[n])
	    {
	      /* Update upper boundary */
	      
	      hi = n;
	      
	      /* Check hit */

	      if ((val >= dat[hi - 1]) && (val < dat[hi]))
		return hi - 1;
	      else if (val == dat[hi])
		return hi;
	    }
	  else
	    {
	      /* Update lower boundary */
	      
	      lo = n;
	      
	      /* Check hit */
	      
	      if ((val >= dat[lo]) && (val < dat[lo + 1]))
		return lo;
	      else if (val == dat[lo + 1])
		return lo + 1;
	    }

	  /*******************************************************************/
	}
    }
    
  /* Shouldn't be here */

  Die(FUNCTION_NAME, "WTF?");

  /* Return index */
  
  return lo;
}

/*****************************************************************************/
