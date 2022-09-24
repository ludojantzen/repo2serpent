/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processmsr.c                                   */
/*                                                                           */
/* Created:       2012/07/14 (JLe)                                           */
/* Last modified: 2017/03/21 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: - Processing stuff for MSR calculations                      */
/*                                                                           */
/* Comments: - Tätä pitää kutsua ennen kuin materiaalit on jaettu palama-    */
/*             alueisiin ja nimet on korvattu pointtereilla.                 */
/*                                                                           */
/*           - Tässä on se potentiaalinen ongelma että muodoissa 92235 ja    */
/*             U-235 annettuja nuklidinimiä ei tunnisteta samaksi.           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessMSR:"

void LoopFlow(long, long, long *);

/*****************************************************************************/

void ProcessMSR()
{
  long n, m, rep, loc0, loc1, ptr, mat1, mat2, iso1, iso2, ZAI;
  char str[MAX_STR];
  
  /***************************************************************************/

  /***** Process mass flow lists *********************************************/

  /* Loop over definitions */

  loc0 = (long)RDB[DATA_PTR_MFLOW0];
  while (loc0 > VALID_PTR)
    {
      /* Check duplicates */

      ptr = NextItem(loc0);
      while (ptr > VALID_PTR)
        {
          /* Compare */

          if (CompareStr(loc0 + MFLOW_PTR_NAME, ptr + MFLOW_PTR_NAME))
            Error(loc0, "Duplicate definition of mass flow %s",
                  GetText(loc0 + MFLOW_PTR_NAME));

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Loop over list */

      ptr = (long)RDB[loc0 + MFLOW_PTR_DATA];
      while (ptr > VALID_PTR)
        {
          /* Pointer to text */

          sprintf(str, "%s", GetText(ptr + MFLOW_LIST_ZAI));

          /* Remove dashes */

          n = 0;
          m = 0;
          
          while (str[n] != '\0')
            {
              if (str[n] != '-')
                str[m++] = str[n];
              
              n++;
            }
      
          str[m] = '\0';

          /* Get ZA */

          if (!strcasecmp(str, "all"))
            ZAI = -1;
          else if ((ZAI = IsotoZAI(str)) < 0)
            ZAI = atoi(str);

          /* Check value */

          if ((ZAI != -1) && (ZAI < 1))
            Error(loc0, "Error in mass flow list: \"%s\" is not a valid entry",
                  GetText(ptr + MFLOW_LIST_ZAI));

          /* Put ZAI */

          WDB[ptr + MFLOW_LIST_ZAI] = (double)ZAI;

          /* Pointer to next */

          ptr = NextItem(ptr);
        }

      /* Next */
      
      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Continuous reprocessing *********************************************/

  /* Loop over reprocessors */

  rep = (long)RDB[DATA_PTR_REP0];
  while (rep > VALID_PTR)
    {
      /* Loop over entries */

      loc0 = (long)RDB[rep + REPROC_PTR_CON_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Find first material */
              
          if ((long)RDB[loc0 + REPROC_CON_PTR_MAT1] < VALID_PTR)
            Die(FUNCTION_NAME, "Material pointer not set");          
              
          /* Loop over materials */
          
          mat1 = (long)RDB[DATA_PTR_M0];
          while (mat1 > VALID_PTR)
            {
              /* Check material */
              
              if (CompareStr(loc0 + REPROC_CON_PTR_MAT1,
                             mat1 + MATERIAL_PTR_NAME))
                break;

              /* Check parent */

              if ((ptr = (long)RDB[mat1 + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
                if (CompareStr(loc0 + REPROC_CON_PTR_MAT1, 
                               ptr + MATERIAL_PTR_NAME))
                  break;
              
              /* Next material */
                  
              mat1 = NextItem(mat1);
            }

          /* Check */

          if (mat1 < VALID_PTR)
            Error(rep, "Material %s does not exist", 
                  GetText(loc0 + REPROC_CON_PTR_MAT1));
          else
            WDB[loc0 + REPROC_CON_PTR_MAT1] = (double)mat1;

          /* Set used- and physical-flags */

          SetOption(mat1 + MATERIAL_OPTIONS, OPT_USED);
          SetOption(mat1 + MATERIAL_OPTIONS, OPT_PHYSICAL_MAT);

          /* Check burnable */

          if (!((long)RDB[mat1 + MATERIAL_OPTIONS] & OPT_BURN_MAT))
            Error(rep, "Material \"%s\" must be burnable", 
                  GetText(mat1 + MATERIAL_PTR_NAME));
      
          /* Find second material */
              
          if ((long)RDB[loc0 + REPROC_CON_PTR_MAT2] < VALID_PTR)
            Die(FUNCTION_NAME, "Material pointer not set");          
              
          /* Loop over materials */

          mat2 = (long)RDB[DATA_PTR_M0];
          while (mat2 > VALID_PTR)
            {
              /* Check material */
              
              if (CompareStr(loc0 + REPROC_CON_PTR_MAT2,
                             mat2 + MATERIAL_PTR_NAME))
                break;

              /* Check parent */

              if ((ptr = (long)RDB[mat2 + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
                if (CompareStr(loc0 + REPROC_CON_PTR_MAT2, 
                               ptr + MATERIAL_PTR_NAME))
                  break;
              
              /* Next material */
                  
              mat2 = NextItem(mat2);
            }

          /* Check */

          if (mat2 < VALID_PTR)
            Error(rep, "Material %s does not exist", 
                  GetText(loc0 + REPROC_CON_PTR_MAT2));
          else
            WDB[loc0 + REPROC_CON_PTR_MAT2] = (double)mat2;

          /* Set used- and physical-flags */

          SetOption(mat2 + MATERIAL_OPTIONS, OPT_USED);
          SetOption(mat2 + MATERIAL_OPTIONS, OPT_PHYSICAL_MAT);

          /* Check that materials are not the same */
          
          if (mat1 == mat2)
            Error(rep, "Duplicate material \"%s\" in continuous reprocessing",
                  GetText(mat1 + MATERIAL_PTR_NAME));

          /* Check burnable */

          if (!((long)RDB[mat2 + MATERIAL_OPTIONS] & OPT_BURN_MAT))
            Error(rep, "Material \"%s\" must be burnable", 
                  GetText(mat2 + MATERIAL_PTR_NAME));

          /* Find mass flow */
              
          if ((long)RDB[loc0 + REPROC_CON_PTR_MFLOW] < VALID_PTR)
            Die(FUNCTION_NAME, "Mass flow pointer not set");          
              
          /* Loop over mass flows */
          
          loc1 = (long)RDB[DATA_PTR_MFLOW0];
          while (loc1 > VALID_PTR)
            {
              /* Compare names */
              
              if (!strcmp(GetText(loc0 + REPROC_CON_PTR_MFLOW), 
                          GetText(loc1 + MFLOW_PTR_NAME)))
                break;
              
              /* Next material */
                  
              loc1 = NextItem(loc1);
            }

          /* Check */

          if (loc1 < VALID_PTR)
            Error(rep, "Mass flow %s does not exist", 
                  GetText(loc0 + REPROC_CON_PTR_MFLOW));
          else
            WDB[loc0 + REPROC_CON_PTR_MFLOW] = (double)loc1;

          /* Set used-flag */

          SetOption(loc1 + MFLOW_OPTIONS, OPT_USED);

          /* Create outflow list */

          ptr = NewItem(mat1 + MATERIAL_PTR_OUTFLOW, REPROC_CON_BLOCK_SIZE);
          for (n = LIST_DATA_SIZE; n < REPROC_CON_BLOCK_SIZE; n++)
            WDB[ptr + n] = RDB[loc0 + n];

          /* Pot pointer to reprocessor */

          WDB[ptr + REPROC_CON_PTR_REP] = (double)rep;

          /* Create inflow list */

          ptr = NewItem(mat2 + MATERIAL_PTR_INFLOW, REPROC_CON_BLOCK_SIZE);
          for (n = LIST_DATA_SIZE; n < REPROC_CON_BLOCK_SIZE; n++)
            WDB[ptr + n] = RDB[loc0 + n];

          /* Pot pointer to reprocessor */

          WDB[ptr + REPROC_CON_PTR_REP] = (double)rep;
          
          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Next reprocessor */

      rep = NextItem(rep);
    }

  /* Remove mass flows */

  loc0 = (long)RDB[DATA_PTR_MFLOW0];
  RemoveFlaggedItems(loc0, MFLOW_OPTIONS, OPT_USED, NO);

  /* Loop over mass flows */

  mat1 = (long)RDB[DATA_PTR_M0];
  while (mat1 > VALID_PTR)
    {
      /* Reset index */

      n = -1;

      /* Loop over flow */

      LoopFlow(mat1, mat1, &n);

      /* Set maximum number of materials */

      if (n > 0)
        WDB[mat1 + MATERIAL_FLOW_N] = (double)n;

      /* Next */

      mat1 = NextItem(mat1);
    }

  /***************************************************************************/

  /***** Add missing nuclides ************************************************/

  /* Loop over materials */

  mat1 = (long)RDB[DATA_PTR_M0];
  while (mat1 > VALID_PTR)
    {
      /* Skip materials not involved in chains */

      if ((long)RDB[mat1 + MATERIAL_FLOW_PTR_FIRST] < VALID_PTR)
        {
          /* Skip */
          
          mat1 = NextItem(mat1);
          
          /* Cycle loop */
          
          continue;
        }

      /* Loop over in same chain */

      mat2 = (long)RDB[DATA_PTR_M0];
      while (mat2 > VALID_PTR)
        {
          /* Check pointer */
          
          if ((long)RDB[mat1 + MATERIAL_FLOW_PTR_FIRST] !=
              (long)RDB[mat2 + MATERIAL_FLOW_PTR_FIRST])
            {
              /* Skip */
          
              mat2 = NextItem(mat2);
              
              /* Cycle loop */

              continue;
            }
          
          /* Loop over first composition */

          iso1 = (long)RDB[mat1 + MATERIAL_PTR_COMP];
          while (iso1 > VALID_PTR)
            {
              /* Loop over second composition and find match */
              
              iso2 = (long)RDB[mat2 + MATERIAL_PTR_COMP];
              while (iso2 > VALID_PTR)
                {
                  /* Compare */

                  if (CompareStr(iso1 + COMPOSITION_PTR_NUCLIDE,
                                 iso2 + COMPOSITION_PTR_NUCLIDE))
                    break;

                  /* Next */
                  
                  iso2 = NextItem(iso2);
                }

              /* Check if found */

              if (iso2 < VALID_PTR)
                {
                  /* Add new */

                  iso2 = NewItem(mat2 + MATERIAL_PTR_COMP, 
                                 COMPOSITION_BLOCK_SIZE);

                  /* Copy nuclide name */

                  WDB[iso2 + COMPOSITION_PTR_NUCLIDE] = 
                    RDB[iso1 + COMPOSITION_PTR_NUCLIDE];

                  /* Density is zero */

                  WDB[iso2 + COMPOSITION_ADENS] = 0.0;
                }

              /* Next */

              iso1 = NextItem(iso1);
            }          

          /* Next material */
          
          mat2 = NextItem(mat2);
        }

      /* Next material */

      mat1 = NextItem(mat1);
    }

  /***************************************************************************/

  /***** Print flow chains ***************************************************/

#ifdef DEBUG

  if ((long)RDB[DATA_PTR_MFLOW0] > VALID_PTR)
    {
      fprintf(outp, "MSR mass flows:\n\n");

      mat1 = (long)RDB[DATA_PTR_M0];
      while (mat1 > VALID_PTR)
        {
          if ((mat2 = (long)RDB[mat1 + MATERIAL_FLOW_PTR_FIRST]) > VALID_PTR)
            fprintf(outp, "%s : %ld %s\n", GetText(mat2 + MATERIAL_PTR_NAME),
                   (long)RDB[mat1 + MATERIAL_FLOW_IDX],
                   GetText(mat1 + MATERIAL_PTR_NAME));
          
          /* Next material */
          
          mat1 = NextItem(mat1);
        }
      fprintf(outp, "\n");
    }

#endif

  /***************************************************************************/
}

/*****************************************************************************/

/***** Loop over material flows in continuous reprocessing *******************/

void LoopFlow(long mat, long mat0, long *idx)
{
  long loc0, ptr;

  /* Check material pointer */
  
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
 
  /* Check if index was already assigned */

  if ((long)RDB[mat + MATERIAL_FLOW_IDX] > 0)
    return;
  
  /* Check index */

  if (*idx < 0)
    {
      /* Not a recursive call, check material inflow */

      if ((long)RDB[mat + MATERIAL_PTR_INFLOW] > VALID_PTR)
        return;

      /* Check outflow */

      if ((long)RDB[mat + MATERIAL_PTR_OUTFLOW] < VALID_PTR)
        return;

      /* Reset index */

      *idx = 0;
    }

  /* Update index */

  *idx = *idx + 1;

  /* Set index and pointer to first */

  WDB[mat + MATERIAL_FLOW_IDX] = (double)(*idx);
  WDB[mat + MATERIAL_FLOW_PTR_FIRST] = (double)(mat0);

  /* Loop over outflows */

  loc0 = (long)RDB[mat + MATERIAL_PTR_OUTFLOW];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to material */

      ptr = (long)RDB[loc0 + REPROC_CON_PTR_MAT2];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Call recursively */

      LoopFlow(ptr, mat0, idx);

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
