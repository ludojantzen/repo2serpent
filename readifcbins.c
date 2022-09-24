/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcbins.c                                  */
/*                                                                           */
/* Created:       2014/07/07 (VVa)                                           */
/* Last modified: 2016/02/01 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Reads user defined bin edges for multi-physics interfaces    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCBins:"

/*****************************************************************************/

long ReadIFCBins(long np, FILE *fp, double min, long update)
{
  long ptr0 = -1, ptr, n;
  double point;

  if (!update)
    ptr0 = ReallocMem(DATA_ARRAY, np+1);

  ptr = ptr0;

  for (n = 0; n < np + 1; n++)
    {
      if(fscanf(fp, "%lf", &point) == EOF)
        Die(FUNCTION_NAME, "Could not read output meshing");

      CheckValue(FUNCTION_NAME, "ifcmesh", "", point, min, INFTY);

      if(!update)
        WDB[ptr++] = point;

      min = point;
    }
          
  return ptr0;

}
