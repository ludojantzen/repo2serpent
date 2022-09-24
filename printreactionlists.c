/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printreactionclists.c                          */
/*                                                                           */
/* Created:       2013/01/23 (JLe)                                           */
/* Last modified: 2019/03/01 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Prints reaction lists for with thresholds and cutoffs        */
/*                                                                           */
/* Comments: - Used for debugging                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintReactionLists:"

/*****************************************************************************/

void PrintReactionLists()
{

#ifdef DEBUG

  long mat, nuc, rea, iso, n, mode, loc0, loc1, idx, imax;
  char tmpstr[MAX_STR];
  FILE *fp;

  return;

 /* Check mpi task */

  if (mpiid > 0)
    return;

  fprintf(outp, "Printing reaction lists for debugging...\n");
  
  /* Open file for writing */
  
  sprintf(tmpstr, "%s_rls%ld.m", GetText(DATA_PTR_INPUT_FNAME),
          (long)RDB[DATA_BURN_STEP]);
  
  if ((fp = fopen(tmpstr, "w")) == NULL)
    Warn(FUNCTION_NAME, "Unable to open file for writing");

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get pointer to composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Get maximum index */

      imax = ListSize(iso) - 1;

      /* Loop over lists */

      for (n = 0; n < 15; n++)
        {
          /* Get mode */
          
          if (n == 0)
            mode = MATERIAL_PTR_TOT_REA_LIST;
          else if (n == 1)
            mode = MATERIAL_PTR_ELA_REA_LIST;
          else if (n == 2)
            mode = MATERIAL_PTR_ABS_REA_LIST;
          else if (n == 3)
            mode = MATERIAL_PTR_FISS_REA_LIST;
          else if (n == 4)
            mode = MATERIAL_PTR_HEATT_REA_LIST;
          else if (n == 5)
            mode = MATERIAL_PTR_PHOTP_REA_LIST;
          else if (n == 6)
            mode = MATERIAL_PTR_INLP_REA_LIST;
          else if (n == 7)
            mode = MATERIAL_PTR_PHOT_TOT_LIST;
          else if (n == 8)
            mode = MATERIAL_PTR_PHOT_HEAT_LIST;
          else if (n == 9)
            mode = MATERIAL_PTR_TOT_URES_LIST;
          else if (n == 10)
            mode = MATERIAL_PTR_ABS_URES_LIST;
          else if (n == 11)
            mode = MATERIAL_PTR_ELA_URES_LIST;
          else if (n == 12)
            mode = MATERIAL_PTR_FISS_URES_LIST;
          else if (n == 13)
            mode = MATERIAL_PTR_HEAT_URES_LIST;
          else if (n == 14)
            mode = MATERIAL_PTR_TMP_MAJORANT_LIST;
          else
            Die(FUNCTION_NAME, "Overflow");

          /* Get pointer to list */

          if ((loc0 = (long)RDB[mat + mode]) < VALID_PTR)
            continue;
          
          if (mode == MATERIAL_PTR_TOT_REA_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_TOT_REA_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_ELA_REA_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_ELA_REA_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_ABS_REA_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_ABS_REA_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_FISS_REA_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_FISS_REA_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_HEATT_REA_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_HEAT_REA_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_PHOTP_REA_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_PHOTP_REA_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_INLP_REA_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_INLP_REA_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_PHOT_TOT_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_PHOT_TOT_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_PHOT_HEAT_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_PHOT_HEAT_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_TOT_URES_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_TOT_URES_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_ABS_URES_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_ABS_URES_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_ELA_URES_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_ELA_URES_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_FISS_URES_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_FISS_URES_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_HEAT_URES_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_HEAT_URES_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));
          else if (mode == MATERIAL_PTR_TMP_MAJORANT_LIST)
            fprintf(fp, "%s: MATERIAL_PTR_TMP_MAJORANT_LIST\n", 
                    GetText(mat + MATERIAL_PTR_NAME));

          /* Loop over list*/

          loc1 = (long)RDB[loc0 + RLS_PTR_REA0];
          while (loc1 > VALID_PTR)
            {
              /* Get composition index */

              idx = (long)RDB[loc1 + RLS_DATA_COMP_IDX];

              /* Check */

              if ((idx < 0) || (idx > imax))
                Die(FUNCTION_NAME, "Error in composition index");

              /* Get pointer to composition */

              iso = ListPtr(iso, idx);
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

              /* Get pointer to nuclide */

              nuc = (long)RDB[loc1 + RLS_DATA_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Pointer to reaction */

              rea = (long)RDB[loc1 + RLS_DATA_PTR_REA];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Print */

              fprintf(fp, "%15s %4ld : %1.5E / %1.5E : [%1.5E %1.5E] ",
                      GetText(nuc + NUCLIDE_PTR_NAME),
                      (long)RDB[rea + REACTION_MT],
                      RDB[iso + COMPOSITION_ADENS],
                      RDB[loc1 + RLS_DATA_MAX_ADENS],
                      RDB[loc1 + RLS_DATA_EMIN],
                      RDB[loc1 + RLS_DATA_EMAX]);

              if ((long)RDB[loc1 + RLS_DATA_CUT] == YES)
                fprintf(fp, "(cut)\n");
              else
                fprintf(fp, "\n");

              /* Next */
              
              loc1 = NextItem(loc1);
            }
        }

      /* Loop over transmutation reaction list */
      
      loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];
      
      if (loc0 > VALID_PTR)
        fprintf(fp, "%s: MATERIAL_PTR_DEP_TRA_LIST\n", 
                GetText(mat + MATERIAL_PTR_NAME));

      while(loc0 > VALID_PTR)
        {
          /* Get pointer to reaction */

          rea = (long)RDB[loc0 + DEP_TRA_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
              
          /* Get pointer to nuclide */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Print */
          
          fprintf(fp, "%15s %4ld : %1.5E : %s\n",
                  GetText(nuc + NUCLIDE_PTR_NAME),
                  (long)RDB[rea + REACTION_MT],
                  RDB[loc0 + DEP_TRA_E0],
                  ReactionMT((long)RDB[rea + REACTION_MT], NO));

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Loop over fission reaction list */
      
      loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];
      
      if (loc0 > VALID_PTR)
        fprintf(fp, "%s: MATERIAL_PTR_DEP_FISS_LIST\n", 
                GetText(mat + MATERIAL_PTR_NAME));

      while(loc0 > VALID_PTR)
        {
          /* Get pointer to reaction */

          rea = (long)RDB[loc0 + DEP_TRA_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
              
          /* Get pointer to nuclide */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Print */
          
          fprintf(fp, "%15s %4ld : %1.5E : %s\n",
                  GetText(nuc + NUCLIDE_PTR_NAME),
                  (long)RDB[rea + REACTION_MT],
                  RDB[loc0 + DEP_TRA_E0],
                  ReactionMT((long)RDB[rea + REACTION_MT], NO));

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Loop over fission neutron production list */
      
      loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_NSF_LIST];
      
      if (loc0 > VALID_PTR)
        fprintf(fp, "%s: MATERIAL_PTR_DEP_NSF_LIST\n", 
                GetText(mat + MATERIAL_PTR_NAME));

      while(loc0 > VALID_PTR)
        {
          /* Get pointer to reaction */

          rea = (long)RDB[loc0 + DEP_TRA_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
              
          /* Get pointer to nuclide */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Print */
          
          fprintf(fp, "%15s %4ld : %1.5E : %s\n",
                  GetText(nuc + NUCLIDE_PTR_NAME),
                  (long)RDB[rea + REACTION_MT],
                  RDB[loc0 + DEP_TRA_E0],
                  ReactionMT((long)RDB[rea + REACTION_MT], NO));

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Newline */

      fprintf(fp, "\n");

      /* Next material */
      
      mat = NextItem(mat);
    }          

  /* Close file */

  fclose(fp);

  /* Exit OK */

  fprintf(outp, "OK.\n\n");

#endif
}

/*****************************************************************************/
