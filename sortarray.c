/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sortarray.c                                    */
/*                                                                           */
/* Created:       2010/09/12 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Sorts an array of doubles                                    */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "SortArray:"

/*****************************************************************************/

void SortArray(double *vec, long sz) 
{
  long n, i, j;
  double tmp;

  /* Perform insertion sort */
  
  for (i = 1; i < sz; i++)
    {
      j = i;
      
      tmp = vec[j];
      
      while (j > 0 && vec[j - 1] > tmp)
	{
	  vec[j] = vec[j - 1];
	  j--;
	}
      
      vec[j] = tmp;
    }  
  
  /* Avoid unused variable warning message in compilation */

  n = 0;

#ifdef DEBUG

  /* Check that array is sorted */

  for (n = 1; n < sz; n++)
    if (vec[n - 1] > vec[n])
      Die(FUNCTION_NAME, "Sorting failed");
  
#endif
}

/*****************************************************************************/
