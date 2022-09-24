/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : newstat.c                                      */
/*                                                                           */
/* Created:       2010/10/22 (JLe)                                           */
/* Last modified: 2011/11/30 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Creates data structures for scores                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NewStat:"

/*****************************************************************************/

long NewStat(char *name, long dim, ...)
{
  long n, loc0, ptr, ntot, bins;
  va_list argp;
  va_start (argp, dim);

  /* Avoid compiler warning */

  loc0 = -1;

  /* Check type */

  if (dim > 0)
    {
      /***********************************************************************/

      /***** Regular statistical variable ************************************/

      /* New item */

      loc0 = NewItem(DATA_PTR_SCORE0, SCORE_BLOCK_SIZE);

      /* Put name */
  
      WDB[loc0 + SCORE_PTR_NAME] = (double)PutText(name);
  
      /* Allocate memory for bin sizes */

      ptr = ReallocMem(DATA_ARRAY, dim);

      /* Put pointer */

      WDB[loc0 + SCORE_PTR_NMAX] = (double)ptr;
      
      /* Put dimension */
      
      WDB[loc0 + SCORE_DIM] = (double)dim;
      
      /* Reset size */
      
      ntot = 1;
      
      /* Loop over dimensions */
      
      for (n = 0; n < dim; n++)
	{
	  /* Get number of bins */
	  
	  bins = va_arg(argp, int);
	  
	  /* Store value */
	  
	  WDB[ptr++] = (double)bins;
	  
	  /* Update total number of values */
	  
	  ntot = ntot*bins;
	}
      
      /* Put stat size */
      
      WDB[loc0 + SCORE_STAT_SIZE] = (double)ntot;
      
      /* Allocate memory for data */

      ptr = ReallocMem(RES1_ARRAY, ntot);
      
      /* Put pointer */
      
      WDB[loc0 + SCORE_PTR_DATA] = (double)ptr;
      
      /* Allocate memory for buffer */
      
      ptr = AllocPrivateData(ntot*BUF_BLOCK_SIZE, BUF_ARRAY);
      
      /* Put pointer */
      
      WDB[loc0 + SCORE_PTR_BUF] = (double)ptr;
      
      /***********************************************************************/
    }
  else if (dim < 0)
    {
      /***********************************************************************/

      /***** Distributed array or mesh ***************************************/

      /* New item */

      loc0 = NewItem(DATA_PTR_SCORE0, SCORE_BLOCK_SIZE);

      /* Put name */
  
      WDB[loc0 + SCORE_PTR_NAME] = (double)PutText(name);

      /* Reset dimension to identify type */

      WDB[loc0 + SCORE_DIM] = -1.0;

      /* Get number of values */
	  
      ntot = va_arg(argp, int);

      /* Allocate memory for data */
      
      ptr = AllocPrivateData(ntot, RES2_ARRAY);
      
      /* Put pointer */
      
      WDB[loc0 + SCORE_PTR_DATA] = (double)ptr;

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Dimension is zero");
  
  /* Return pointer */

  return loc0;
}

/*****************************************************************************/
