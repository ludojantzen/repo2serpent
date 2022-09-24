/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectmicrodepdata.c                          */
/*                                                                           */
/* Created:       2016/12/30 (JLe)                                           */
/* Last modified: 2016/12/30 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Collects cross section for micro-depletion output            */
/*                                                                           */
/* Comments: - Koko systeemi muutettiin ennen jakelua                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectMicroDepData:"

/*****************************************************************************/

void CollectMicroDepData()
{
#ifdef mmmmmmmmmmmmmmmmmmmmmmm

  long loc0, loc1, mat, iso, nuc, rea, id, key, mat0;
  double br, flx, adens, rr, xs, vol;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Check micro-depletion data */

  if ((long)RDB[DATA_PTR_MDEP0] < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

 /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) 
    return;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Open file for writing */

  sprintf(tmpstr, "%s_mdx%ld.m", GetText(DATA_PTR_INPUT_FNAME),
          (long)RDB[DATA_BURN_STEP]);

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  fprintf(fp, "\n --- Material-wise data: -----------------\n\n");

  /***************************************************************************/
  
  /***** Reset and prepare data **********************************************/

  /* Loop over micro depletion data */
  
  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      /* Sort list by key */
      
      loc1 = (long)RDB[loc0 + MDEP_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
      SortList(loc1, MDEP_REA_KEY, SORT_MODE_ASCEND);
      
      /* Reset data */

      loc1 = (long)RDB[loc0 + MDEP_PTR_REA];
      while (loc1 > VALID_PTR)
        {
          /* Reset reaction rate and divisor */

          WDB[loc1 + MDEP_REA_DIV] = 0.0;
          WDB[loc1 + MDEP_REA_RR] = 0.0;

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Collect reaction rates **********************************************/

  /* Loop over micro depletion data */

  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over materials */
      
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if material list is given */
          
          if ((loc1 = (long)RDB[loc0 + MDEP_PTR_MAT]) > VALID_PTR)
            {
              /* Find material in list */

              while ((mat0 = (long)RDB[loc1++]) > VALID_PTR)
                {
                  /* Compare material */
                  
                  if (mat == mat0)
                    break;
                  
                  /* Compare parent */
                  
                  if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] == mat0)
                    break;
                }
              
              /* Check if found */

              if (mat0 < VALID_PTR)
                {
                  /* Pointer to next */
                  
                  mat = NextItem(mat);
                  
                  /* Cycle loop */
                  
                  continue;
                }
            }

          /* Get flux */
          
          flx = Truncate(RDB[mat + MATERIAL_BURN_FLUX_SSA], 6);
          CheckValue(FUNCTION_NAME, "flx", "", flx, 0.0, INFTY);      
          
          /* Check flux */
          
          if (flx == 0.0)
            {
              /* Pointer to next */
              
              mat = NextItem(mat);
              
              /* Cycle loop */
              
              continue;
            }
          
          /* Get volume */
          
          vol = RDB[mat + MATERIAL_VOLUME];
          CheckValue(FUNCTION_NAME, "flx", "", flx, ZERO, INFTY);
                    
          fprintf(fp, "material %s\n", GetText(mat + MATERIAL_PTR_NAME));
          
          /* Loop over composition */
          
          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Get atomic density */
              
              adens = RDB[iso + COMPOSITION_ADENS];
              
              /* Pointer to nuclide data */
              
              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /***************************************************************/

              /***** Collect reaction data ***********************************/

              /* Loop over reactions */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  /* Check type and target (fission is not included) */

                  if (((long)RDB[rea + REACTION_TYPE] != 
                       REACTION_TYPE_PARTIAL) ||
                      ((long)RDB[rea + REACTION_PTR_TGT] < VALID_PTR))
                    {
                      /* Next reaction */
                      
                      rea = NextItem(rea);
                      
                      /* Cycle loop */
                      
                      continue;
                    }
                  
                  /* Get branching ratio */
                  
                  br = RDB[rea + REACTION_BR];
                  
                  /* Check pointer to xs data */
                  
                  if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] < VALID_PTR)
                    Die(FUNCTION_NAME, "Pointer err1or");
                  
                  /* Sum data over threads */
                  
                  rr = 0.0;
                  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
                    if ((xs = TestValuePair(rea + REACTION_PTR_TRANSMUXS, 
                                            (double)mat, id)) > 0.0)
                      rr = rr + flx*xs;
                  
                  fprintf(fp, "%s %E %E %ld\n", GetText(nuc + NUCLIDE_PTR_NAME),
                          adens, rr/flx, (long)RDB[rea + REACTION_MT]);
                  
                  /* Search key */
                  
                  key = 1000*((long)RDB[nuc + NUCLIDE_ZAI]) +
                    (long)RDB[rea + REACTION_MT];

                  /* Check */
                  
                  if (rr > 0.0)
                    {
                      /* Find reaction */
                      
                      loc1 = (long)RDB[loc0 + MDEP_PTR_REA];

                      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
                      loc1 = SeekList(loc1, MDEP_REA_KEY, (double)key, 
                                      SORT_MODE_ASCEND);

                      /* Check */
                      
                      if (loc1 > VALID_PTR)
                        {
                          /* Store data */

                          WDB[loc1 + MDEP_REA_RR] = RDB[loc1 + MDEP_REA_RR] +
                            rr*adens*vol;
                          WDB[loc1 + MDEP_REA_DIV] = RDB[loc1 + MDEP_REA_DIV] +
                            flx*adens*vol;
                        }
                    }

                  /* Next reaction */
              
                  rea = NextItem(rea);
                }

              /***************************************************************/
              
              /***** Collect total fission data ******************************/

              /* Check pointer to total fission */
              
              if ((rea = (long)RDB[nuc + NUCLIDE_PTR_TOTFISS_REA]) > VALID_PTR)
                {
                  /* Sum data over threads */
                  
                  rr = 0.0;
                  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
                    if ((xs = TestValuePair(rea + REACTION_PTR_TRANSMUXS, 
                                            (double)mat, id)) > 0.0)
                      rr = rr + flx*xs;
                  
                  fprintf(fp, "%s %E %E 18\n", GetText(nuc + NUCLIDE_PTR_NAME),
                          adens, rr/flx);
                  
                  /* Search key */
                  
                  key = 1000*((long)RDB[nuc + NUCLIDE_ZAI]) + 18;
                  
                  /* Check */
                  
                  if (rr > 0.0)
                    {
                      /* Find reaction */
                      
                      loc1 = (long)RDB[loc0 + MDEP_PTR_REA];
                      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
                      loc1 = SeekList(loc1, MDEP_REA_KEY, (double)key, 
                                      SORT_MODE_ASCEND);
                      
                      /* Check */
                      
                      if (loc1 > VALID_PTR)
                        {
                          /* Store data */
                          
                          WDB[loc1 + MDEP_REA_RR] = 
                            RDB[loc1 + MDEP_REA_RR] + rr*adens*vol;
                          WDB[loc1 + MDEP_REA_DIV] = 
                            RDB[loc1 + MDEP_REA_DIV] + flx*adens*vol;
                        }
                    }
                }
              
              /***************************************************************/
              
              /* Next */
              
              iso = NextItem(iso);
            }
          
          /* Next */
          
          mat = NextItem(mat);
        }

      /* Next */
      
      loc0 = NextItem(loc0);
    }

  fprintf(fp, "\n --- Final output: -----------------------\n\n");

  /* Loop over micro depletion data */
  
  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      loc1 = (long)RDB[loc0 + MDEP_PTR_GCU];
      loc1 = (long)RDB[loc1 + GCU_PTR_UNIV];
      fprintf(fp, "Universe: %s\n", GetText(loc1 + UNIVERSE_PTR_NAME));
      
      if ((loc1 = (long)RDB[loc0 + MDEP_PTR_MAT]) > VALID_PTR)
        {
          fprintf(fp, "materials:");
          while ((mat = (long)RDB[loc1++]) > VALID_PTR)
            fprintf(fp, " %s", GetText(mat + MATERIAL_PTR_NAME));
          
          fprintf(fp, "\n");
        }

      /* Sort list by index */
      
      loc1 = (long)RDB[loc0 + MDEP_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
      SortList(loc1, MDEP_REA_IDX, SORT_MODE_ASCEND);
      
      /* Print data */

      loc1 = (long)RDB[loc0 + MDEP_PTR_REA];
      while (loc1 > VALID_PTR)
        {
          /* Calculate average */

          if (RDB[loc1 + MDEP_REA_DIV] > 0.0)
            WDB[loc1 + MDEP_REA_RR] = 
              RDB[loc1 + MDEP_REA_RR]/RDB[loc1 + MDEP_REA_DIV];
          else if (RDB[loc1 + MDEP_REA_RR] > 0.0)
            Die(FUNCTION_NAME, "WTF?");

          fprintf(fp, "%6ld %3ld %E\n", (long)RDB[loc1 + MDEP_REA_ZAI],
                  (long)RDB[loc1 + MDEP_REA_MT], RDB[loc1 + MDEP_REA_RR]);

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Close file and exit */

  fclose(fp);

#endif
}

/*****************************************************************************/
