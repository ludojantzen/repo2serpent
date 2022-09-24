/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : intersectionlist.c                             */
/*                                                                           */
/* Created:       2010/10/18 (JLe)                                           */
/* Last modified: 2011/11/11 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Generates surface intersection list for cell.                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IntersectionList:"

/*****************************************************************************/

void IntersectionList(long cell, long *infix, long nparam)
{
  long i, n, ptr, ptp;

  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

  /* Check number of parameters */

  if (nparam < 1)
    Die(FUNCTION_NAME, "No parameters");

  /* Reset side */

  i = 1;
  
  /* Loop over operators */

  for (n = 0; n < nparam; n++)
    {
      /* Check for complement */
      
      if (infix[n] == SURF_OP_NOT)
	i = -1;
      else if (infix[n] == SURF_OP_AND)
	{
	  /* Do nothing */
	}
      else if (infix[n] > 0)
	{
	  /* Create new item in intersection list */
	  
	  ptr = NewItem(cell + CELL_PTR_SURF_INSC, CELL_INSC_BLOCK_SIZE);
	  
	  /* Put side */
	  
	  WDB[ptr + CELL_INSC_SIDE] = (double)i;
	  
	  /* Put surface pointer */
	  
	  WDB[ptr + CELL_INSC_PTR_SURF] = (double)infix[n];
	  
	  /* Allocate memory for counter from private array */
	  
	  ptp = AllocPrivateData(1, PRIVA_ARRAY);
	  
	  /* Put pointer */

	  WDB[ptr + CELL_INSC_PTR_OUT_COUNT] = (double)ptp;

	  /* Reset side */
	  
	  i = 1;
	}
      else
	Die(FUNCTION_NAME, "Invalid operator %ld", infix[n]);
    }
}

/*****************************************************************************/
