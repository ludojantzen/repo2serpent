/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : finixptrfromfpe.c                              */
/*                                                                           */
/* Created:       2015/09/03 (VVa)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Gets pointer to FINIX block based on IFC_FUEP_BLOCK          */
/*                                                                           */
/* Comments:  -                                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UpdateFinixIFC:"

/*****************************************************************************/

long FinixPtrFromFpe(long fpe)
{
  long fib, found, nu, ptr, i;
  
  /* Check that some FINIX blocks are defined */
  
  if ((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return NULLPTR;

  /* Find correct pin */

  while (fib > VALID_PTR)
    {
      /* Reset found flag */

      found = 0;

      /* Pointer to number of rod segments */

      nu = (long)RDB[fpe + IFC_FUEP_N_UNI];

      /* Loop over rod segments */

      ptr = (long)RDB[fpe + IFC_FUEP_PTR_UNI_LIST];

      for (i=0; i < nu; i++)
        {
          /* Compare finix universe name and IFC_FUEP segment names */

          if (CompareStr(ptr + i, fib + FINIX_PTR_UNI_NAME))
            found = 1;
        }

      /* Break if found, otherwise test next FINIX block */

      if (found==1)
        break;
      else
        fib = NextItem(fib);         

    }

  /* Check if found */

  if (fib < VALID_PTR)
    return NULLPTR;
  else
    return fib;
}
