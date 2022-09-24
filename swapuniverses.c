/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : swapuniverses.c                                */
/*                                                                           */
/* Created:       2012/07/09 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Swaps the contents of two universes                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SwapUniverses:"

/*****************************************************************************/

void SwapUniverses(long uni1, long uni2)
{
  long cell, nst, reg, lat, pbl, pbd, ptr, n;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(uni1)", DATA_ARRAY, uni1);
  CheckPointer(FUNCTION_NAME, "(uni2)", DATA_ARRAY, uni2);

  /* Loop over cells */

  cell = (long)RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /* Check fill pointer and change */
            
      if ((long)RDB[cell + CELL_PTR_FILL] == uni1)
        WDB[cell + CELL_PTR_FILL] = (double)uni2;
      else if ((long)RDB[cell + CELL_PTR_FILL] == uni2)
        WDB[cell + CELL_PTR_FILL] = (double)uni1;
      
      /* Next cell */
      
      cell = NextItem(cell);
    }
  
  /* Loop over nests */

  nst = (long)RDB[DATA_PTR_NST0];
  while (nst > VALID_PTR)
    {
      /* Loop over regions */ 

      reg = (long)RDB[nst + NEST_PTR_REGIONS];              
      while (reg > VALID_PTR)
        {
          /* Check fill pointer and change */

          if ((long)RDB[reg + NEST_REG_PTR_FILL] == uni1)
            WDB[reg + NEST_REG_PTR_FILL] = (double)uni2;
          else if ((long)RDB[reg + NEST_REG_PTR_FILL] == uni2)
            WDB[reg + NEST_REG_PTR_FILL] = (double)uni1;
          
          /* Next region */
          
          reg = NextItem(reg);
        }
      
      /* Next nest */

      nst = NextItem(nst);
    }

  /* Loop over lattices */

  lat = (long)RDB[DATA_PTR_L0];
  while (lat > VALID_PTR)
    {
      /* Check type */
      
      if ((long)RDB[lat + LAT_TYPE] == LAT_TYPE_CLU)
        {
          /* Loop over rings */ 
          
          reg = (long)RDB[lat + LAT_PTR_FILL];              
          while (reg > VALID_PTR)
            {
              /* Loop over items */ 
              
              ptr = (long)RDB[reg + RING_PTR_FILL];
              while ((long)RDB[ptr] > VALID_PTR)
                {
                  /* Check fill pointer and change */

                  if ((long)RDB[ptr] == uni1)
                    WDB[ptr] = (double)uni2;
                  else if ((long)RDB[ptr] == uni2)
                    WDB[ptr] = (double)uni1;
                  
                  /* Next */

                  ptr++;
                }

              /* Next region */
              
              reg = NextItem(reg);
            }
        }
      else
        {
          /* Pointer to items */ 
            
          ptr = (long)RDB[lat + LAT_PTR_FILL];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Loop over items and call recursively (tuolla voi olla */
          /* se NULLPTR välissä) */
          
          for (n = 0; n < (long)RDB[lat + LAT_NTOT]; n++)
            {
              /* Check fill pointer and change */
              
              if ((long)RDB[ptr] == uni1)
                WDB[ptr] = (double)uni2;
              else if ((long)RDB[ptr] == uni2)
                WDB[ptr] = (double)uni1;
              
              /* Next */
              
              ptr++;
            }
        }
      
      /* Next lattice */
      
      lat = NextItem(lat);
    }

  /* Loop over pebble bed geometries */

  pbd = (long)RDB[DATA_PTR_PB0];
  while (pbd > VALID_PTR)
    {
      /* Check background universe pointer and change */
            
      if ((long)RDB[pbd + PBED_PTR_BG_UNIV] == uni1)
        WDB[pbd + PBED_PTR_BG_UNIV] = (double)uni2;
      else if ((long)RDB[pbd + PBED_PTR_BG_UNIV] == uni2)
        WDB[pbd + PBED_PTR_BG_UNIV] = (double)uni1;

      /* Loop over pebbles */
          
      pbl = (long)RDB[pbd + PBED_PTR_PEBBLES];
      while (pbl > VALID_PTR)
        {
          /* Check pointer and change */
          
          if ((long)RDB[pbl + PEBBLE_PTR_UNIV] == uni1)
            WDB[pbl + PEBBLE_PTR_UNIV] = (double)uni2;
          else if ((long)RDB[pbl + PEBBLE_PTR_UNIV] == uni2)
            WDB[pbl + PEBBLE_PTR_UNIV] = (double)uni1;
          
          /* Next pebble */
          
          pbl = NextItem(pbl);
        }

      /* Next geometry */

      pbd = NextItem(pbd);
    }

  /***************************************************************************/

}

/*****************************************************************************/
