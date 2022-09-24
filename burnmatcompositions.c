/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : burnmatcompositions.c                          */
/*                                                                           */
/* Created:       2012/05/14 (JLe)                                           */
/* Last modified: 2020/03/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Puts nuclide compositions to burnable materials              */
/*                                                                           */
/* Comments: - This was separated from processburnmat.c 14.5.2012 (2.1.6)    */
/*                                                                           */
/*           - Toi minimienergia riippuu nyt spectrum collapse menetelmästä. */
/*             Jos menetelmä on käytössä, energia on ures-alueen minimi, jos */
/*             ei, se on reaktion kynnysenergia.                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BurnMatCompositions:"

/*****************************************************************************/

void BurnMatCompositions()
{
  long mat, iso, nuc, last, iso1, iso2, nuc1, nuc2, rea, ptr, n, nmax;
  double tmp, mem;

  /* Check burnup mode */

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
    return;

  /* Check burnup step */

  if (((long)RDB[DATA_BURN_STEP] > 0) &&
      ((long)RDB[DATA_PTR_HISV0]) < VALID_PTR)
    Die(FUNCTION_NAME, "Should not be here");

  /* Print */

  fprintf(outp, "Setting nuclide compositions for burnable materials...\n");

  /* Reset maximum */

  nmax = -1.0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /***********************************************************************/

      /***** Check flags, etc. ***********************************************/

      /* Check burn flag and divisor mode */

      if ((!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)) ||
          ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR))
        {
          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check composition */

      if ((long)RDB[mat + MATERIAL_PTR_COMP] < VALID_PTR)
        Die(FUNCTION_NAME, "Material %s has no nuclides",
            GetText(mat + MATERIAL_PTR_NAME));

      /* Check mixture */

      if ((long)RDB[mat + MATERIAL_PTR_MIX] > VALID_PTR)
        Die(FUNCTION_NAME, "Burnable mixture %s?",
            GetText(mat + MATERIAL_PTR_NAME));

      /* Check div type */

      if (((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_S1) ||
          ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW))
        Die(FUNCTION_NAME, "Divided material");

      /* Calculate bytes */

      CalculateBytes();

      /* Get memory size */

      mem = RDB[DATA_TOTAL_BYTES];

      /***********************************************************************/

      /***** Add nuclides from decay and transmutation chains ****************/

      /* Reset used-flags on nuclides (these are used to identify */
      /* processed chains in AddChains() subroutine). */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Reset used- and exists-flags */

          ResetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);
          ResetOption(nuc + NUCLIDE_OPTIONS, OPT_EXISTS);

          /* Next nuclide */

          nuc = NextItem(nuc);
        }

      /* Set exists-flags for nuclides initial composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Set flag */

          SetOption(nuc + NUCLIDE_OPTIONS, OPT_EXISTS);

          /* Next */

          iso = NextItem(iso);
        }

      /* Get pointer to composition and last initial nuclide */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      last = LastItem(iso);
      CheckPointer(FUNCTION_NAME, "(last)", DATA_ARRAY, last);

      /* Loop over initial composition and add chains */

      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Add chain to material */

          AddChains(mat, nuc, 0);

          /* Pointer to next or break if last */

          if (iso == last)
            break;
          else
            iso = NextItem(iso);
        }

      /***********************************************************************/

      /***** Sort composition by ZAI and add nuclide for lost data ***********/

      /* Add nuclide for lost data (must be at the beginning to get */
      /* matrix shape right) */

      iso = NewItem(mat + MATERIAL_PTR_COMP, COMPOSITION_BLOCK_SIZE);
      WDB[iso + COMPOSITION_PTR_NUCLIDE] = RDB[DATA_PTR_NUCLIDE_LOST];
      WDB[iso + COMPOSITION_ADENS] = 0.0;

      /* Swap nuclide ZAI and composition atomic density (nuclide ZAI */
      /* is used as temporary storage space) */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Swap data */

          tmp = RDB[nuc + NUCLIDE_ZAI];
          WDB[nuc + NUCLIDE_ZAI] = RDB[iso + COMPOSITION_ADENS];
          WDB[iso + COMPOSITION_ADENS] = tmp;

          /* Next */

          iso = NextItem(iso);
        }

      /* Sort list by ZAI stored in composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      SortList(iso, COMPOSITION_ADENS, SORT_MODE_ASCEND);

      /* Swap nuclide ZAI and composition atomic density */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Swap data */

          tmp = RDB[nuc + NUCLIDE_ZAI];
          WDB[nuc + NUCLIDE_ZAI] = RDB[iso + COMPOSITION_ADENS];
          WDB[iso + COMPOSITION_ADENS] = tmp;

          /* Next */

          iso = NextItem(iso);
        }

      /***********************************************************************/

      /***** Check and remove duplicate ZAI's ********************************/

      /* Loop until all duplicates are removed */

      do
        {
          /* Reset count */

          n = 0;

          /* Loop over composition */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Copy pointers */

              iso1 = iso;
              iso2 = NextItem(iso);

              /* Check if next exists */

              if (iso2 > VALID_PTR)
                {
                  /* Get nuclide pointers */

                  nuc1 = (long)RDB[iso1 + COMPOSITION_PTR_NUCLIDE];
                  nuc2 = (long)RDB[iso2 + COMPOSITION_PTR_NUCLIDE];

                  CheckPointer(FUNCTION_NAME, "(nuc1)", DATA_ARRAY, nuc1);
                  CheckPointer(FUNCTION_NAME, "(nuc2)", DATA_ARRAY, nuc2);

                  /* Check oneis burnable absorber */

                  if (((long)RDB[nuc1 + NUCLIDE_TYPE_FLAGS] &
                       NUCLIDE_FLAG_BA) ||
                      ((long)RDB[nuc2 + NUCLIDE_TYPE_FLAGS] &
                       NUCLIDE_FLAG_BA))
                    {
                      /* Pointer to next */

                      iso = iso2;

                      /* Cycle loop */

                      continue;
                    }

                  /* Compare ZAI's */

                  if ((long)RDB[nuc1 + NUCLIDE_ZAI] >
                      (long)RDB[nuc2 + NUCLIDE_ZAI])
                    Die(FUNCTION_NAME, "Composition not sorted");

                  /* Toi lämpötilan testaus lisättiin 4.2.2014 / 2.1.17 */

                  else if (((long)RDB[nuc1 + NUCLIDE_ZAI] ==
                            (long)RDB[nuc2 + NUCLIDE_ZAI]) &&
                           (RDB[nuc1 + NUCLIDE_TEMP] ==
                            RDB[nuc2 + NUCLIDE_TEMP]))
                    {
                      /*******************************************************/

                      /***** Tää on uutta 4.2.2016 / 2.1.25 (JLe) ************/

                      /* Check material S(a,b) pointer */

                      if ((long)RDB[mat + MATERIAL_PTR_SAB] > VALID_PTR)
                        {
                          /* Check nuclide S(a,b) flags */

                          if (((long)RDB[nuc1 + NUCLIDE_TYPE_FLAGS] &
                               NUCLIDE_FLAG_SAB_DATA) &&
                              ((long)RDB[nuc2 + NUCLIDE_TYPE_FLAGS] &
                               NUCLIDE_FLAG_SAB_DATA))
                            {
                              /* Both set, remove the one with zero */
                              /* density */

                              if (RDB[iso1 + COMPOSITION_ADENS] > 0.0)
                                {
                                  /* Get new next item since removing current */

                                  ptr = NextItem(iso2);

                                  /* Remove current next item */

                                  RemoveItem(iso2);

                                  /* Deal with the next next item next */

                                  iso2 = ptr;
                                }
                              else
                                RemoveItem(iso1);
                            }
                          else if ((long)RDB[nuc1 + NUCLIDE_TYPE_FLAGS] &
                                   NUCLIDE_FLAG_SAB_DATA)
                            {
                              /* Get new next item since removing current */

                              ptr = NextItem(iso2);

                              /* Remove current next item */

                              RemoveItem(iso2);

                              /* Deal with the next next item next */

                              iso2 = ptr;
                            }
                          else if ((long)RDB[nuc2 + NUCLIDE_TYPE_FLAGS] &
                                   NUCLIDE_FLAG_SAB_DATA)
                            RemoveItem(iso1);
                          else
                            Die(FUNCTION_NAME, "No S(a,b) flags? (%s %s) 1",
                                GetText(mat + MATERIAL_PTR_NAME),
                                GetText(nuc + NUCLIDE_PTR_NAME));

                          n++;
                        }
                      else
                        {
                          /* Remove the one with S(a,b) flag */

                          if ((long)RDB[nuc1 + NUCLIDE_TYPE_FLAGS] &
                                   NUCLIDE_FLAG_SAB_DATA)
                            RemoveItem(iso1);
                          else if ((long)RDB[nuc2 + NUCLIDE_TYPE_FLAGS] &
                                   NUCLIDE_FLAG_SAB_DATA)
                            {
                              /* Get new next item since removing current */

                              ptr = NextItem(iso2);

                              /* Remove current next item */

                              RemoveItem(iso2);

                              /* Deal with the next next item next */

                              iso2 = ptr;
                            }
                          else
                            Die(FUNCTION_NAME, "No S(a,b) flags? (%s %s) 2",
                                GetText(mat + MATERIAL_PTR_NAME),
                                GetText(nuc + NUCLIDE_PTR_NAME));

                          n++;
                        }


                      /*******************************************************/
#ifdef WANHAA
                      /* Check S(a,b) flag on first nuclide */

                      if ((long)RDB[nuc1 + NUCLIDE_TYPE_FLAGS] &
                          NUCLIDE_FLAG_SAB_DATA)
                        {
                          /* Check density for second nuclide */

                          if (RDB[iso2 + COMPOSITION_ADENS] > 0.0)
                            Die(FUNCTION_NAME, "adens2 > 0");

                          /* Remove second nuclide */

                          RemoveItem(iso2);

                          /* Add count */

                          n++;
                        }

                      /* Check S(a,b) flag on second nuclide */

                      else if ((long)RDB[nuc2 + NUCLIDE_TYPE_FLAGS] &
                               NUCLIDE_FLAG_SAB_DATA)
                        {
                          /* Check density for first nuclide */

                          if (RDB[iso1 + COMPOSITION_ADENS] > 0.0)
                            Die(FUNCTION_NAME, "adens1 > 0");

                          /* Remove first nuclide */

                          Warn(FUNCTION_NAME, "removing %s",
                               GetText(nuc1 + NUCLIDE_PTR_NAME));

                          RemoveItem(iso1);

                          /* Add count */

                          n++;
                        }

                      /* Something wrong here */

                      else
                        Warn(FUNCTION_NAME, "Duplicates in %s: %s %s",
                             GetText(mat + MATERIAL_PTR_NAME),
                             GetText(nuc1 + NUCLIDE_PTR_NAME),
                             GetText(nuc2 + NUCLIDE_PTR_NAME));

                      /* Break loop */

                      break;

#endif
                    }
                }

              /* Next */

              iso = iso2;
            }
        }
      while (n > 0);

      /***********************************************************************/

      /***** Create transmutation and fission lists **************************/

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Check type and target ZAI */

              if (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL) &&
                  ((long)RDB[rea + REACTION_TGT_ZAI] > 0))
                {
                  /* Add to transmutation list */

                  ptr  = NewItem(mat + MATERIAL_PTR_DEP_TRA_LIST,
                                 DEP_TRA_BLOCK_SIZE);

                  /* Put reaction pointer */

                  WDB[ptr + DEP_TRA_PTR_REA] = (double)rea;

                  /* Put minimum energy */

                  if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == NO)
                    WDB[ptr + DEP_TRA_E0] = RDB[rea + REACTION_EMIN];
                }
              else if (((long)RDB[rea + REACTION_TYPE] ==
                        REACTION_TYPE_TRA_BRANCH) &&
                       ((long)RDB[rea + REACTION_PTR_ISO_BRA] > VALID_PTR))
                {
                  /* Add branching to isomeric state to transmutation list */

                  ptr  = NewItem(mat + MATERIAL_PTR_DEP_TRA_LIST,
                                 DEP_TRA_BLOCK_SIZE);

                  /* Put reaction pointer */

                  WDB[ptr + DEP_TRA_PTR_REA] = (double)rea;

                  /* Put minimum energy */

                  if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == NO)
                    WDB[ptr + DEP_TRA_E0] = RDB[rea + REACTION_EMIN];
                }
              else if (((long)RDB[rea + REACTION_TYPE]
                        != REACTION_TYPE_DECAY) &&
                       ((long)RDB[rea + REACTION_TGT_ZAI] < 0))
                {
                  /* Add to fission list */

                  ptr  = NewItem(mat + MATERIAL_PTR_DEP_FISS_LIST,
                                 DEP_TRA_BLOCK_SIZE);

                  /* Put reaction pointer */

                  WDB[ptr + DEP_TRA_PTR_REA] = (double)rea;

                  /* Put minimum energy */

                  if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == NO)
                    {
                      if (RDB[rea + REACTION_FISSY_IE0]
                          < RDB[rea + REACTION_EMIN])
                        WDB[ptr + DEP_TRA_E0] = RDB[rea + REACTION_EMIN];
                      else
                        WDB[ptr + DEP_TRA_E0] = RDB[rea + REACTION_FISSY_IE0];
                    }

                  /* Add to fission neutron production list */

                  if (((long)RDB[DATA_BURN_CALC_NSF] == YES) &&
                      (RDB[rea + REACTION_FISSY_IE0] == -INFTY))
                    {
                      ptr  = NewItem(mat + MATERIAL_PTR_DEP_NSF_LIST,
                                     DEP_TRA_BLOCK_SIZE);

                      /* Put reaction pointer */

                      WDB[ptr + DEP_TRA_PTR_REA] = (double)rea;

                      /* Put minimum energy */

                      WDB[ptr + DEP_TRA_E0] = RDB[rea + REACTION_EMIN];
                    }

                }

              /* Next reaction */

              rea = NextItem(rea);
            }

          /* Next nuclide */

          iso = NextItem(iso);
        }

      /* Close lists and sort by energy */

      if ((ptr = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST]) > VALID_PTR)
        {
          CloseList(ptr);

          if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == NO)
            SortList(ptr, DEP_TRA_E0, SORT_MODE_ASCEND);
        }

      if ((ptr = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST]) > VALID_PTR)
        {
          CloseList(ptr);

          if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == NO)
            SortList(ptr, DEP_TRA_E0, SORT_MODE_ASCEND);
        }

      if ((ptr = (long)RDB[mat + MATERIAL_PTR_DEP_NSF_LIST]) > VALID_PTR)
        {
          CloseList(ptr);

          if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == NO)
            SortList(ptr, DEP_TRA_E0, SORT_MODE_ASCEND);
        }

      /***********************************************************************/

      /* Calculate bytes */

      CalculateBytes();

      /* Update memory size */

      WDB[mat + MATERIAL_MEM_SIZE] = RDB[mat + MATERIAL_MEM_SIZE] +
        RDB[DATA_TOTAL_BYTES] - mem;

      /* Get compisition size and compare to maximum */

      ptr = (long)RDB[mat + MATERIAL_PTR_COMP];
      if (ListSize(ptr) > nmax)
        nmax = ListSize(ptr);

      /* Next material */

      mat = NextItem(mat);
    }

  /* Allocate memory for work arrays */

  if (nmax < 1)
    Error(0, "Burnup calculation without burnable materials in geometry");
  else
    {
      /*
        WorkArray(DATA_PTR_WORK_COMP1, PRIVA_ARRAY, nmax, 0);
        WorkArray(DATA_PTR_WORK_COMP2, PRIVA_ARRAY, nmax, 0);
      */
    }

  /***************************************************************************/

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
