/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkduplicates.c                              */
/*                                                                           */
/* Created:       2010/10/06 (JLe)                                           */
/* Last modified: 2018/09/08 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Checks duplicate input definitions                           */
/*                                                                           */
/* Comments: - This routine is used for checking only, calling it should not */
/*             change anything.                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckDuplicates:"
#define MAX_LIST_SIZE 5000

/*****************************************************************************/

void CheckDuplicates()
{
  long loc0, loc1, loc2, ptr, nm, nc, nn;

  /* Reset counters */

  nm = 0;
  nc = 0;
  nn = 0;

  /* This subroutine is a major bottleneck if the number of materials or */
  /* cells is large --> skip. */

  loc0 = (long)RDB[DATA_PTR_M0];
  if ((loc0 > VALID_PTR) && ListSize(loc0) > MAX_LIST_SIZE)
    {
      /* Print note */

      Note(0, "Large input file, duplicate definitions not checked");

      /* Exit */

      return;
    }

  loc0 = (long)RDB[DATA_PTR_C0];
  if ((loc0 > VALID_PTR) && (ListSize(loc0) > MAX_LIST_SIZE))
    {
      /* Print note */

      Note(0, "Large input file, duplicate definitions not checked");

      /* Exit */

      return;
    }

  fprintf(outp, "Checking duplicate input definitions...\n");  

  /***************************************************************************/

  /***** Materials ***********************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_M0];
  while (loc0 > VALID_PTR)
    {
      /* Add counter */

      nm++;

      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + MATERIAL_PTR_NAME, ptr + MATERIAL_PTR_NAME))
            Error(loc0, "Duplicate definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
      
      /* Loop over thermal scattering data and check duplicate ZA's */

      loc1 = (long)RDB[loc0 + MATERIAL_PTR_SAB];
      while (loc1 > VALID_PTR)
        {
          /* Second loop */

          ptr = NextItem(loc1);
          
          while (ptr > VALID_PTR)
            {
              /* Compare */
              
              if (RDB[loc1 + THERM_ZA] == RDB[ptr + THERM_ZA])
                Error(loc0, "S(a,b) ZA %ld used with multiple libraries",
                      (long)RDB[loc1 + THERM_ZA]);
              
              /* Next */
              
              ptr = NextItem(ptr);
            }

          /* Next */
          
          loc1 = NextItem(loc1);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Thermal scattering data *********************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_T0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + THERM_PTR_ALIAS, ptr + THERM_PTR_ALIAS))
            Error(loc0, "Duplicate definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Surfaces ************************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_S0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + SURFACE_PTR_NAME, ptr + SURFACE_PTR_NAME))
            Error(loc0, "Duplicate definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Cells ***************************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_C0];
  while (loc0 > VALID_PTR)
    {
      /* Add counter */

      nc++;

      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + CELL_PTR_NAME, ptr + CELL_PTR_NAME))
            Error(loc0, "Duplicate definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Energy grids ********************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_ENE0];
  while (loc0 > VALID_PTR)
    {
      /* Add counter */

      nc++;

      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + ENE_PTR_NAME, ptr + ENE_PTR_NAME))
            Error(loc0, "Duplicate definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Nests ***************************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_NST0];
  while (loc0 > VALID_PTR)
    {
      /* Add counter */

      nn++;

      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + NEST_PTR_NAME, ptr + NEST_PTR_NAME))
            Error(loc0, 
                  "Duplicate universe definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Lattices ************************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_L0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + LAT_PTR_NAME, ptr + LAT_PTR_NAME))
            Error(loc0, 
                  "Duplicate universe definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Pebble-bed geometries************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_PB0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + PBED_PTR_NAME, ptr + PBED_PTR_NAME))
            Error(loc0, 
                  "Duplicate universe definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Pbeds vs. nests *****************************************************/

  /* First loop */

  loc0 = (long)RDB[DATA_PTR_PB0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = (long)RDB[DATA_PTR_NST0];
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + PBED_PTR_NAME, ptr + NEST_PTR_NAME))
            Error(loc0, 
                  "Duplicate universe definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Pbeds vs. lattices **************************************************/

  /* First loop */

  loc0 = (long)RDB[DATA_PTR_PB0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = (long)RDB[DATA_PTR_L0];
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + PBED_PTR_NAME, ptr + LAT_PTR_NAME))
            Error(loc0, 
                  "Duplicate universe definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Lattices vs. nests **************************************************/

  /* First loop */

  loc0 = (long)RDB[DATA_PTR_L0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = (long)RDB[DATA_PTR_NST0];
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + LAT_PTR_NAME, ptr + NEST_PTR_NAME))
            Error(loc0, 
                  "Duplicate universe definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Cell universes vs. nests, lattices and pbeds ************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_C0];
  while (loc0 > VALID_PTR)
    {
      /* Nests */

      ptr = (long)RDB[DATA_PTR_NST0];
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + CELL_PTR_UNI, ptr + NEST_PTR_NAME))
            Error(loc0, 
                  "Duplicate universe definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Lattices */

      ptr = (long)RDB[DATA_PTR_L0];
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + CELL_PTR_UNI, ptr + LAT_PTR_NAME))
            Error(loc0, 
                  "Duplicate universe definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Pbeds */

      ptr = (long)RDB[DATA_PTR_PB0];
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + CELL_PTR_UNI, ptr + PBED_PTR_NAME))
            Error(loc0, 
                  "Duplicate universe definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Sources *************************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_SRC0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + SRC_PTR_NAME, ptr + SRC_PTR_NAME))
            Error(loc0, "Duplicate definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }
  
  /***************************************************************************/

  /***** Detectors ***********************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_DET0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + SRC_PTR_NAME, ptr + SRC_PTR_NAME))
            Error(loc0, "Duplicate definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }
  
  /***************************************************************************/

  /***** Material divisors ***************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_DIV0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + DIV_PTR_MAT, ptr + DIV_PTR_MAT))
            Error(loc0, "Duplicate definition on line %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }
  
  /***************************************************************************/

  /***** Variance reduction **************************************************/
  
  /* First loop */

  loc0 = (long)RDB[DATA_PTR_RMX0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + RMX_PTR_NAME, ptr + RMX_PTR_NAME))
            Error(loc0, "Duplicate definition on %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }

  /* First loop */

  loc0 = (long)RDB[DATA_PTR_WWD0];
  while (loc0 > VALID_PTR)
    {
      /* Second loop */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */
          
          if (CompareStr(loc0 + WWD_PTR_NAME, ptr + WWD_PTR_NAME))
            Error(loc0, "Duplicate definition on %ld in file \"%s\"",
                  (long)RDB[ptr + PARAM_LINE], 
                  GetText(ptr + PARAM_PTR_FNAME));

          /* Next */

          ptr = NextItem(ptr);
        }
 
      /* Next */

      loc0 = NextItem(loc0);
    }
  
  /***************************************************************************/

  /***** Identical burnable materials in cells and nests *********************/

   /* First loop */

  loc0 = (long)RDB[DATA_PTR_C0];
  while (loc0 > VALID_PTR)
    {
      /* Check material pointer */

      if ((long)RDB[loc0 + CELL_PTR_MAT] > VALID_PTR)
        {
          /* Nests */
          
          loc1 = (long)RDB[DATA_PTR_NST0];
          while (loc1 > VALID_PTR)
            {
              /* Loop over regions */

              ptr = (long)RDB[loc1 + NEST_PTR_REGIONS];
              while (ptr > VALID_PTR)
                {
                  /* Material name is in cell pointer */

                  if ((long)RDB[ptr + NEST_REG_PTR_CELL] > VALID_PTR)
                    {
                      /* Compare */

                      if (CompareStr(loc0 + CELL_PTR_MAT, 
                                     ptr + NEST_REG_PTR_CELL))
                        {
                          /* Find material */

                          loc2 = (long)RDB[DATA_PTR_M0];
                          while (loc2 > VALID_PTR)
                            {
                              /* Compare name and check burn-flag */
                              
                              if ((CompareStr(loc0 + CELL_PTR_MAT, 
                                              loc2 + MATERIAL_PTR_NAME)) &&
                                  (long)RDB[loc2 + MATERIAL_OPTIONS] &
                                  OPT_BURN_MAT)
                                Error(loc1, "Burnable material %s is not allowed both in pin structure and cell %s", GetText(loc2 + MATERIAL_PTR_NAME),
                                      GetText(loc0 + CELL_PTR_NAME));
                                

                              /* Next material */

                              loc2 = NextItem(loc2);
                            }
                        }
                    }

                  /* Next region */
                  
                  ptr = NextItem(ptr);
                }
              
              /* Next nest */
              
              loc1 = NextItem(loc1);
            }
        }
 
      /* Next cell */

      loc0 = NextItem(loc0);
    }
  
  /***************************************************************************/
  
  /***** Check counters ******************************************************/

  /* Check all counters */

  if ((nc == 0) && (nn == 0) && (nm == 0))
    Error(0, "Not a Serpent input file");

  /* Check geometry counters */

  if ((nc == 0) && (nn == 0))
    Error(0, "No geometry definition");

  /* Check material counter */

  if (nm == 0)
    Error(0, "No material definitions");

  /* Everything OK. */

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
