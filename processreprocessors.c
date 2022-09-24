/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processreprocessors.c                          */
/*                                                                           */
/* Created:       2012/07/14 (JLe)                                           */
/* Last modified: 2015/04/10 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: - Processing stuff for reprocessors                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessReprocessors:"

/*****************************************************************************/

void ProcessReprocessors()
{
  long dep, n, rep, loc0, ptr, uni, max, mat1, mat2, iso1, iso2;
  
  /***************************************************************************/

  /***** Link reprocessors to depletion histories ****************************/

  /* Loop over depletion histories */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Check pointer */

      if ((long)RDB[dep + DEP_HIS_PTR_REPROC] > VALID_PTR)
        {
          /* Find reprocessor */

          rep = (long)RDB[DATA_PTR_REP0];
          while (rep > VALID_PTR)
            {
              /* Compare */

              if (CompareStr(dep + DEP_HIS_PTR_REPROC, rep + REPROC_PTR_NAME))
                break;

              /* Next */

              rep = NextItem(rep);
            }
          
          /* Check pointer */

          if (rep < VALID_PTR)
            Error(dep, "Reprocessor %s not defined", 
                  GetText(dep + DEP_HIS_PTR_REPROC));

          /* Set pointer */

          WDB[dep + DEP_HIS_PTR_REPROC] = (double)rep;

          /* Set used-flag */

          SetOption(rep + REPROC_OPTIONS, OPT_USED);
        }

      /* Next */

      dep = NextItem(dep);
    }

  /* Remove unused */

  rep = (long)RDB[DATA_PTR_REP0];
  RemoveFlaggedItems(rep, REPROC_OPTIONS, OPT_USED, NO);

  /***************************************************************************/

  /***** Reprocessed universes ***********************************************/

  /* Remember maximum level */

  max = (long)RDB[DATA_GEOM_LEVELS];

  /* Loop over reprocessors */

  rep = (long)RDB[DATA_PTR_REP0];
  while (rep > VALID_PTR)
    {
      /***********************************************************************/

      /***** Swaps ***********************************************************/
      
      /* Loop over list */

      loc0 = (long)RDB[rep + REPROC_PTR_SWAP_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Loop over universes */

          for (n = 0; n < 2; n++)
            {
              /* Get pointer */
              
              if (n == 0)
                ptr = loc0 + REPROC_SWAP_PTR_UNI1;
              else
                ptr = loc0 + REPROC_SWAP_PTR_UNI2;
             
              /* Check pointer */

              if ((long)RDB[ptr] < VALID_PTR)
                Die(FUNCTION_NAME, "Universe pointer not set");
 
              /* Loop over existing */
                  
              uni = (long)RDB[DATA_PTR_U0];
              while (uni > VALID_PTR)
                {
                  /* Compare names */
                  
                  if (!strcmp(GetText(ptr), GetText(uni + UNIVERSE_PTR_NAME)))
                    break;
                      
                  /* Next universe */
                  
                  uni = NextItem(uni);
                }
              
              /* Check if found */
                  
              if (uni < VALID_PTR)
                {
                  /* Create new universe */
                  
                  uni = CreateUniverse(rep, GetText(ptr), 0);
                  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
                }

              /* Put pointer */

              WDB[ptr] = (double)uni;
            }
          
          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /***** Replacements ****************************************************/
      
      loc0 = (long)RDB[rep + REPROC_PTR_RPL_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Loop over universes */

          for (n = 0; n < 2; n++)
            {
              /* Get pointer */
              
              if (n == 0)
                ptr = loc0 + REPROC_RPL_PTR_UNI1;
              else
                ptr = loc0 + REPROC_RPL_PTR_UNI2;
              
              /* Check pointer */
              
              if ((long)RDB[ptr] < VALID_PTR)
                Die(FUNCTION_NAME, "Universe pointer not set");          
              
              /* Loop over existing */
              
              uni = (long)RDB[DATA_PTR_U0];
              while (uni > VALID_PTR)
                {
                  /* Compare names */
                  
                  if (!strcmp(GetText(ptr), GetText(uni + UNIVERSE_PTR_NAME)))
                    break;
                  
                  /* Next universe */
                  
                  uni = NextItem(uni);
                }
              
              /* Check if found */
              
              if (uni < VALID_PTR)
                {
                  /* Create new universe */
                  
                  uni = CreateUniverse(rep, GetText(ptr), 0);
                  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
                }

              /* Put pointer */

              WDB[ptr] = (double)uni;
            }
          
          /* Next */
          
          loc0 = NextItem(loc0);
        }
      
      /***********************************************************************/

      /* Next reprocessor */

      rep = NextItem(rep);
    }
  
  /* Check levels */

  if ((long)RDB[DATA_GEOM_LEVELS] > max)
    Die(FUNCTION_NAME, "Change in maximum levels");

  /***************************************************************************/

  /***** Removed materials ***************************************************/

  /* Loop over reprocessors */

  rep = (long)RDB[DATA_PTR_REP0];
  while (rep > VALID_PTR)
    {
      /* Loop over list */

      loc0 = (long)RDB[rep + REPROC_PTR_REM_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Find first material */
              
          if ((long)RDB[loc0 + REPROC_REM_PTR_MAT1] < VALID_PTR)
            Die(FUNCTION_NAME, "Material pointer not set");          
              
          /* Loop over existing */
          
          mat1 = (long)RDB[DATA_PTR_M0];
          while (mat1 > VALID_PTR)
            {
              /* Compare names */
              
              if (!strcmp(GetText(loc0 + REPROC_REM_PTR_MAT1), 
                          GetText(mat1 + MATERIAL_PTR_NAME)))
                break;
              
              /* Next material */
                  
              mat1 = NextItem(mat1);
            }

          /* Check */

          if (mat1 < VALID_PTR)
            Error(rep, "Material %s does not exist", 
                  GetText(loc0 + REPROC_REM_PTR_MAT1));
          else
            WDB[loc0 + REPROC_REM_PTR_MAT1] = (double)mat1;

          /* Find second material */
              
          if ((long)RDB[loc0 + REPROC_REM_PTR_MAT2] < VALID_PTR)
            Die(FUNCTION_NAME, "Material pointer not set");          
              
          /* Loop over existing */
          
          mat2 = (long)RDB[DATA_PTR_M0];
          while (mat2 > VALID_PTR)
            {
              /* Compare names */
              
              if (!strcmp(GetText(loc0 + REPROC_REM_PTR_MAT2), 
                          GetText(mat2 + MATERIAL_PTR_NAME)))
                break;
              
              /* Next material */
                  
              mat2 = NextItem(mat2);
            }

          /* Check */

          if (mat2 < VALID_PTR)
            {
              /* Create new */
              
              mat2 = NewItem(DATA_PTR_M0, MATERIAL_BLOCK_SIZE);
            }

          /* Put pointer */

          WDB[loc0 + REPROC_REM_PTR_MAT2] = (double)mat2;

          /* Merge */
          
          iso1 = (long)RDB[mat1 + MATERIAL_PTR_COMP];
          while (iso1 > VALID_PTR)
            {
              /* Find existing */

              iso2 = (long)RDB[mat2 + MATERIAL_PTR_COMP];
              if (SeekListStr(iso2, COMPOSITION_PTR_NUCLIDE, 
                              GetText(iso1 + COMPOSITION_PTR_NUCLIDE)) 
                  < VALID_PTR)
                {
                  /* Create new */
                  
                  iso2 = NewItem(mat2 + MATERIAL_PTR_COMP, 
                                 COMPOSITION_BLOCK_SIZE);

                  /* Put name */

                  WDB[iso2 + COMPOSITION_PTR_NUCLIDE] =
                    RDB[iso1 + COMPOSITION_PTR_NUCLIDE];
                }
              
              /* Next */

              iso1 = NextItem(iso1);
            }

          /* Set used-flags */

          SetOption(mat1 + MATERIAL_OPTIONS, OPT_USED);
          SetOption(mat2 + MATERIAL_OPTIONS, OPT_USED);

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Next reprocessor */

      rep = NextItem(rep);
    }

  /***************************************************************************/
}

/*****************************************************************************/
