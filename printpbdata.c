/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printpbdata.c                                  */
/*                                                                           */
/* Created:       2011/11/20 (JLe)                                           */
/* Last modified: 2011/11/20 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Prints power distributions for explicit geometry             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintPBData:"

/*****************************************************************************/

void PrintPBData()
{
  long pbd, pbl, ptr, uni, idx;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Break if corrector step */

  if ((long)RDB[DATA_BURN_PRED_STEP] == CORRECTOR_STEP)
    return;

  /* Loop over geometries */
  
  pbd = (long)RDB[DATA_PTR_PB0];  
  while (pbd > VALID_PTR)
    {
      /* Check if results are requested */
      
      if ((long)RDB[pbd + PBED_CALC_RESULTS] == YES)
	{
	  /* Pointer to results */

	  ptr = (long)RDB[pbd + PBED_PTR_POW];

	  /* Open file for writing */

	  sprintf(tmpstr, "%s.out", GetText(pbd + PBED_PTR_FNAME));
	  
	  if ((fp = fopen(tmpstr, "w")) == NULL) 
	    Die(FUNCTION_NAME, "Unable to open file for writing");
 
	  /* Loop over pebbles */
	  
	  pbl = (long)RDB[pbd + PBED_PTR_PEBBLES];
	  while (pbl > 0)
	    {
	      /* Print coordinates */
	      
	      fprintf(fp, "%12.5E %12.5E %12.5E ", RDB[pbl + PEBBLE_X0],
		      RDB[pbl + PEBBLE_Y0], RDB[pbl + PEBBLE_Z0]);
	      
	      /* Print radius */
	      
	      fprintf(fp, "%12.5E ", RDB[pbl + PEBBLE_RAD]);

	      /* Print universe */
	      
	      uni = (long)RDB[pbl + PEBBLE_PTR_UNIV];
	      fprintf(fp, "%6s ", GetText(uni + UNIVERSE_PTR_NAME));

	      /* Pebble index */

	      idx = (long)RDB[pbl + PEBBLE_IDX];
	      
	      /* Print fission power */

	      fprintf(fp, "%12.5E %7.5f\n", Mean(ptr, idx), RelErr(ptr, idx));

	      /* Next pebble */
	      
	      pbl = NextItem(pbl);
	    }

	  /* Close file */

	  fclose(fp);
	}
      
      /* Next geometry */
      
      pbd = NextItem(pbd);
    }  
}

/*****************************************************************************/
