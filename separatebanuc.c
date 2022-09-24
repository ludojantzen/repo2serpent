/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : separatebanuc.c                                */
/*                                                                           */
/* Created:       2019/11/13 (JLe)                                           */
/* Last modified: 2019/11/15 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Separates nuclides marked as burnable absorbers into two     */
/*              different species.                                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SeparateBANuc:"

/*****************************************************************************/

void SeparateBANuc()
{
  long loc0, mat, iso0, nuc0, iso1, nuc1, rea, n, ZAI0;

  /* Check if list is provided */

  if ((loc0 = (long)RDB[DATA_PTR_BA_CHAIN]) < VALID_PTR)
    return;

  fprintf(outp, "Separating burnable absorbers from chains:\n\n");

  /* Loop over chain */

  while ((ZAI0 = (long)RDB[loc0++]) > 0)
    {
      /* Print */

      fprintf(outp, "- %s\n", ZAItoIso(ZAI0, 1));

      /***********************************************************************/

      /***** Handle compositions *********************************************/

      /* Reset nuclide pointer */

      nuc1 = -1;

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check if domain decomposition is in use */

          if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
              ((long)RDB[mat + MATERIAL_MPI_ID] > -1) &&
              ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
            {
              /* Pointer to next */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Check burn flag */

          if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
            {
              /* Next material */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Loop over composition */

          iso0 = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso0 > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc0 = (long)RDB[iso0 + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc0)", DATA_ARRAY, nuc0);

              /* Check for non-zero initial density */

              if (RDB[iso0 + COMPOSITION_ADENS] != 0.0)
                {
                  /* Compare ZAI and type */

                  if ((long)RDB[nuc0 + NUCLIDE_ZAI] == ZAI0)
                    {
                      /* Check for fission yield data */

                      if (((long)RDB[nuc0 + NUCLIDE_PTR_NFY_DATA] > VALID_PTR)
                          ||
                          ((long)RDB[nuc0 + NUCLIDE_PTR_SFY_DATA] > VALID_PTR))
                        Die(FUNCTION_NAME, "Fission yields not processed");

                      /* Duplicate composition and reset density */

                      iso1 = DuplicateItem(iso0);
                      WDB[iso0 + COMPOSITION_ADENS] = 0.0;

                      /* Check if duplicate nuclide exists */

                      if (nuc1 < VALID_PTR)
                        {
                          /* Duplicate */

                          nuc1 = DuplicateItem(nuc0);

                          /* Reset reaction pointer */

                          WDB[nuc1 + NUCLIDE_PTR_REA] = 0.0;

                          /* Re-read ACE file */

                          ReadACEFile(nuc1);

                          /* Add branching */

                          AddBranching(nuc1);

                          /* Set mt 4 type to special (see addnuclide.c) */

                          rea = (long)RDB[nuc1 + NUCLIDE_PTR_REA];
                          while (rea > VALID_PTR)
                            {
                              /* Check mt */

                              if (((long)RDB[rea + REACTION_MT] == 4) &&
                                  ((long)RDB[rea + REACTION_TYPE] ==
                                   REACTION_TYPE_PARTIAL))
                                WDB[rea + REACTION_TYPE] =
                                  (double)REACTION_TYPE_SPECIAL;

                              /* Next */

                              rea = NextItem(rea);
                            }

                          /* Set BA flag */

                          SetOption(nuc1 + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_BA);
                        }

                      /* Convert A, ZA and ZAI */

                      WDB[nuc1 + NUCLIDE_ZAI] = RDB[nuc1 + NUCLIDE_ZAI]
                        + 1000.0;
                      WDB[nuc1 + NUCLIDE_A] = RDB[nuc1 + NUCLIDE_A] + 100.0;
                      WDB[nuc1 + NUCLIDE_ZA] = RDB[nuc1 + NUCLIDE_ZA] + 100.0;

                      /* Put pointer */

                      WDB[iso1 + COMPOSITION_PTR_NUCLIDE] = (double)nuc1;

                      /* Pointer to next */

                      iso0 = iso1;
                    }
                }

              /* Next */

              iso0 = NextItem(iso0);
            }

          /* Pointer to next */

          mat = NextItem(mat);
        }

      /* Check if found */

      if (nuc1 < VALID_PTR)
        Die(FUNCTION_NAME, "not found");

      /* Re-index nuclides */

      n = 1;

      nuc0 = (long)RDB[DATA_PTR_NUC0];
      while (nuc0 > VALID_PTR)
        {
          /* Put new index */

          WDB[nuc0 + NUCLIDE_IDX] = (double)(n++);

          /* Next */

          nuc0 = NextItem(nuc0);
        }

      /***********************************************************************/
    }

  /****************************************************************************/

  /* Newline */

  fprintf(outp, "\n");

  /****************************************************************************/
}

/*****************************************************************************/
