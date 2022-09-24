/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processcpd.c                                   */
/*                                                                           */
/* Created:       2011/05/15 (JLe)                                           */
/* Last modified: 2018/06/09 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Process full-core power distributions                        */
/*                                                                           */
/* Comments: - Tää ei nyt välttämättä toimi kaikissa mahdollisissa           */
/*             geometrioissa                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessCPD:"

/*****************************************************************************/

void ProcessCPD()
{
  long l, l1, l2, lat, uni, n1, n2, n3, ptr, nd;

  /* Check if distribution is given */

  if ((nd = (long)RDB[DATA_CORE_PDE_DEPTH]) < 0.0)
    return;

  /* Find core lattice (lowest level) */

  l1 = MAX_GEOMETRY_LEVELS;

  /* Loop overl lattices */
  
  lat = (long)RDB[DATA_PTR_L0];
  while (lat > VALID_PTR)
    {
      /* Pointer to universe */
      
      uni = (long)RDB[lat + LAT_PTR_UNI];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Get level */
      
      if ((long)RDB[uni + UNIVERSE_LEVEL] < l1)
        l1 = (long)RDB[uni + UNIVERSE_LEVEL];
      
      /* Next lattice */
      
      lat = NextItem(lat);
    }

  /* Check level */

  if (l1 == MAX_GEOMETRY_LEVELS)
    Error(0, "No lattices for core power distribution");

  /* Find pin lattice (highest level) */

  l2 = l1;

  /* Loop overl lattices */
  
  lat = (long)RDB[DATA_PTR_L0];
  while (lat > VALID_PTR)
    {
      /* Pointer to universe */
      
      uni = (long)RDB[lat + LAT_PTR_UNI];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Get level */
      
      if ((long)RDB[uni + UNIVERSE_LEVEL] > l2)
        l2 = (long)RDB[uni + UNIVERSE_LEVEL];
      
      /* Next lattice */
      
      lat = NextItem(lat);
    }

  /* Override with user-given values */

  if ((l = (long)RDB[DATA_CORE_PDE_L1]) > -1)
    {
      if ((l < l1) || (l > l2))
        Error(0, "Invalid level %ld in core power distribution", l);
      else
        l1 = l;
    }

  if ((l = (long)RDB[DATA_CORE_PDE_L2]) > -1)
    {
      if ((l < l1) || (l > l2))
        Error(0, "Invalid level %ld in core power distribution", l);
      else
        l2 = l;
    }

  /* Check pointer */

  if (l2 == l1)
    l2 = -1;

  /* Reset counters */

  n1 = 0;
  n2 = 0;

  /* Read lattices */

  lat = (long)RDB[DATA_PTR_L0];
  while (lat > VALID_PTR)
    {
      /* Pointer to universe */
      
      uni = (long)RDB[lat + LAT_PTR_UNI];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Get level */

      if ((long)RDB[uni + UNIVERSE_LEVEL] == l1)
        {
          /* Add core lattice */

          ptr = NewItem(DATA_CORE_PDE_PTR_CORE, CPD_BLOCK_SIZE);

          /* Put pointer */

          WDB[ptr + CPD_PTR_LAT] = (double)lat;

          /* Compare to size */

          if ((long)RDB[lat + LAT_NTOT] > n1)
            n1 = (long)RDB[lat + LAT_NTOT];
        }
      else if (((long)RDB[uni + UNIVERSE_LEVEL] == l2) && (nd > 1))
        {
          /* Add assembly lattice */

          ptr = NewItem(DATA_CORE_PDE_PTR_ASS, CPD_BLOCK_SIZE);

          /* Put pointer */

          WDB[ptr + CPD_PTR_LAT] = (double)lat;

          /* Compare to size */

          if ((long)RDB[lat + LAT_NTOT] > n2)
            n2 = (long)RDB[lat + LAT_NTOT];
        }
      
      /* Next lattice */
      
      lat = NextItem(lat);
    }

  /* Check core level */

  if (n1 < 1)
    Die(FUNCTION_NAME, "No core level bins");

  /* Get number of axial regions */

  n3 = (long)RDB[DATA_CORE_PDE_NZ];

  /* Allocate memory for core level results */

  ptr = NewStat("CPD_ASS", 1, n1);
  WDB[DATA_CORE_PDE_PTR_RES0] = (double)ptr;  

  /* Allocate memory for pin level results */

  if (n2 > 0)
    {  
      ptr = NewStat("CPD_PIN", 2, n1, n2);
      WDB[DATA_CORE_PDE_PTR_RES1] = (double)ptr;  
    }
  
  /* Allocate memory for region level results */

  if (n3 > 0)
    {
      /* Check pin level */

      if (n2 < 1)
        {
          ptr = NewStat("CPD_REG", 3, n1, 1, n3);
          WDB[DATA_CORE_PDE_PTR_RES2] = (double)ptr;  
        }
      else
        {
          ptr = NewStat("CPD_REG", 3, n1, n2, n3);
          WDB[DATA_CORE_PDE_PTR_RES2] = (double)ptr;  
        }
    }    

  /* Put sizes */
  
  WDB[DATA_CORE_PDE_N0] = (double)n1;
  WDB[DATA_CORE_PDE_N1] = (double)n2;
  WDB[DATA_CORE_PDE_N2] = (double)n3;
}

/*****************************************************************************/
