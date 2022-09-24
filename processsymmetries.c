/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processsymmetries.c                            */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2018/08/10 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Links symmetries to universes, etc.                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessSymmetries:"

/*****************************************************************************/

void ProcessSymmetries()
{
  long sym, uni, n;

  /* Loop over symmetries */

  sym = (long)RDB[DATA_PTR_SYM0];
  while (sym > VALID_PTR)
    {
      /* Find universe */

      uni = (long)RDB[DATA_PTR_U0];
      while (uni > VALID_PTR)
        {
          /* Compare names */

          if (CompareStr(uni + UNIVERSE_PTR_NAME, sym + SYMMETRY_PTR_UNI))
            break;

          /* Next universe */

          uni = NextItem(uni);
        }

      /* Check if found */

      if (uni < VALID_PTR)
        {
          /* Print warning */

          Note(sym, "Universe %s not defined", 
                GetText(sym + SYMMETRY_PTR_UNI));

          /* Skip symmetry */

          sym = NextItem(sym);

          /* Cycle loop */

          continue;
        }
      else if ((long)RDB[uni + UNIVERSE_PTR_SYM] > VALID_PTR)
        Error(0, "Universe %s already associated with a symmetry",
              GetText(uni + UNIVERSE_PTR_NAME));      

      /* Check coordinate transformation */
      
      if ((long)RDB[sym + SYMMETRY_COORD_TRANS] == YES)
        {
          /* Check that universe is root */

          if (uni != (long)RDB[DATA_PTR_ROOT_UNIVERSE])
            Error(sym, 
                  "Coordinate transformations only allowed on root universe");
        }

      /* Set pointers */
        
      WDB[sym + SYMMETRY_PTR_UNI] = (double)uni;
      WDB[uni + UNIVERSE_PTR_SYM] = (double)sym;

      /* Check Serpent 1-type symmetries */

      if ((n = (long)RDB[sym + SYMMETRY_SYM]) != 0)
        {
          /* Put axis and boundary condition */

          WDB[sym + SYMMETRY_AXIS] = (double)PutText("z");
          WDB[sym + SYMMETRY_BC] = (double)PutText("refl");
          
          /* Put angles */

          if (n == 4)
            {
              WDB[sym + SYMMETRY_THETA0] = 0.0;
              WDB[sym + SYMMETRY_ROT] = 90.0;
             }
        }

      /* Set axis */

      if (!strcasecmp(GetText(sym + SYMMETRY_AXIS), "x"))
        WDB[sym + SYMMETRY_AXIS] = 1.0;
      else if (!strcasecmp(GetText(sym + SYMMETRY_AXIS), "1"))
        WDB[sym + SYMMETRY_AXIS] = 1.0;
      else if (!strcasecmp(GetText(sym + SYMMETRY_AXIS), "y"))
        WDB[sym + SYMMETRY_AXIS] = 2.0;
      else if (!strcasecmp(GetText(sym + SYMMETRY_AXIS), "2"))
        WDB[sym + SYMMETRY_AXIS] = 2.0;
      else if (!strcasecmp(GetText(sym + SYMMETRY_AXIS), "z"))
        WDB[sym + SYMMETRY_AXIS] = 3.0;
      else if (!strcasecmp(GetText(sym + SYMMETRY_AXIS), "3"))
        WDB[sym + SYMMETRY_AXIS] = 3.0;
      else
        Error(0, "Invalid axis mode \"%s\"", GetText(sym + SYMMETRY_AXIS));

      /* Set boundary condition */
      
      if (!strcasecmp(GetText(sym + SYMMETRY_BC), "refl"))
        WDB[sym + SYMMETRY_BC] = (double)BC_REFLECTIVE;
      else if (!strcasecmp(GetText(sym + SYMMETRY_BC), "2"))
        WDB[sym + SYMMETRY_BC] = (double)BC_REFLECTIVE;
      else if (!strcasecmp(GetText(sym + SYMMETRY_BC), "peri"))
        WDB[sym + SYMMETRY_BC] = (double)BC_PERIODIC;
      else if (!strcasecmp(GetText(sym + SYMMETRY_BC), "3"))
        WDB[sym + SYMMETRY_BC] = (double)BC_PERIODIC;
      else
        Error(0, "Invalid boundary condition \"%s\"", 
              GetText(sym + SYMMETRY_BC));
      
      /* Convert angles to rad */

      WDB[sym + SYMMETRY_THETA0] = RDB[sym + SYMMETRY_THETA0]*PI/180.0;
      WDB[sym + SYMMETRY_ROT] = RDB[sym + SYMMETRY_ROT]*PI/180.0;

      /* Check origin (doesn't work if not zeros) */

      if ((RDB[sym + SYMMETRY_X0] != 0.0) || (RDB[sym + SYMMETRY_Y0]))
        Die(FUNCTION_NAME, "Not centered at origin");
      
      /* TODO: Tarkista sopivuus reunaehtojen kanssa */

      /* Next symmetry */

      sym = NextItem(sym);
    }
}

/*****************************************************************************/
