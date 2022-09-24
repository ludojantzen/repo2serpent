/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nestvolumes.c                                  */
/*                                                                           */
/* Created:       2011/07/03 (JLe)                                           */
/* Last modified: 2012/02/01 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Calculates nest volumes                                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NestVolumes:"

/*****************************************************************************/

void NestVolumes()
{
  long nst, reg, cell, in, out;
  double vol;

  /* Loop over nests */

  nst = (long)RDB[DATA_PTR_NST0];
  while (nst > VALID_PTR)
    {
      /* Get pointer to regions */
	
      reg = (long)RDB[nst + NEST_PTR_REGIONS];
      CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);
      
      /* Loop over regions */ 
	
      while (reg > VALID_PTR)
	{
	  /* Get pointer to cell */
	  
	  cell = (long)RDB[reg + NEST_REG_PTR_CELL];
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
	  
	  /* Pointers to inside and outside surfaces */
	  
	  in = (long)RDB[reg + NEST_REG_PTR_SURF_IN];
	  out = (long)RDB[reg + NEST_REG_PTR_SURF_OUT];
	  
	  /* outside last surface or infinite single-region cell */
	  
	  if (in < VALID_PTR)
	    vol = INFTY;
	  else
	    {
	      /* Volume of inside surface */
	      
	      vol = SurfaceVol(in);
	      
	      /* Minus volume of outside region */
	      
	      if (out > VALID_PTR)
		vol = vol - SurfaceVol(out);
	    }
	  
	  /* Put volume */
	  
	  if (vol > 0.0)
	    WDB[cell + CELL_VOLUME] = vol;
	  else
	    WDB[cell + CELL_VOLUME] = 0.0;

	  /* Next region */

	  reg = NextItem(reg);
	}

      /* Next nest */

      nst = NextItem(nst);
    }      
}

/*****************************************************************************/
