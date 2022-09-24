/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processiternucs.c                              */
/*                                                                           */
/* Created:       2017/02/28 (VVa)                                           */
/* Last modified: 2018/06/05 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Adds data for the calculation of critical concentrations for */
/*              user chosen nuclides                                         */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIterNucs:"

/*****************************************************************************/

void ProcessIterNucs()
{
  long mat, iso, nuc, ptr, inuc0, inuc, ZAI, ZAIlist;
  long maxmat, i, ncompfound, loc0, N;
  long *compositions, *materials;
  double adens, maxadens;

  /* Check mode */

  if ((long)RDB[DATA_ITER_MODE] != ITER_MODE_NUCLIDE)
    return;

  fprintf(outp, "Processing nuclides to be iterated...\n\n");

  /* Get total number of materials (maximum number of materials */
  /* included in the iteration) */
  /* TODO: Materials divided for depletion */

  maxmat = 0;

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      maxmat++;

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Create material list if not given ***********************************/

  /* Get pointer to first (dummy) block*/

  inuc0 = (long)RDB[DATA_PTR_ITER_NUC0];

  if ((ptr = (long)RDB[inuc0 + ITER_NUCLIDE_PTR_MATERIAL_LIST]) > VALID_PTR)
    {
      /* Material list given, swap names to pointers */

      /* Loop over material list */
      while ((long)RDB[ptr] > VALID_PTR)
        {

          /* Find this material from global materials */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              if (CompareStr(ptr, mat + MATERIAL_PTR_NAME))
                break;

              mat = NextItem(mat);
            }

          /* Check if found */

          if (mat < VALID_PTR)
            Error(inuc0, "Unable to find material \"%s\" linked to iterable"
                  " nuclide", GetText(ptr));

          /* Swap material name pointer for material pointer */

          WDB[ptr] = (double)mat;

          /* Check that material is physical */

          if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT) == 0)
            Error(inuc0, "Material %s manually flagged for density iteration not "
                  "physical. Not included in calculation?",
                  GetText(mat + MATERIAL_PTR_NAME));

          /* Check that material is included in majorant */

          if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_INCLUDE_MAJORANT) == 0)
            Error(inuc0, "Material %s manually flagged for density iteration not "
                  "included in majorant. Non-physical?",
                  GetText(mat + MATERIAL_PTR_NAME));

            /* Next material */

          ptr++;
        }

    }
  else
    {
      /* Material list not given */
      /* Find all materials that contain the ZAIs to be flagged */

      /* Reset number of found materials */

      N = 0;

      /* Allocate temporary list for collecting materials */

      materials = (long *)Mem(MEM_ALLOC, maxmat, sizeof(long));

      /* Loop over the material list */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {

          /* Skip materials that are not included in majorant or not physical */

          if ((((long)RDB[mat + MATERIAL_OPTIONS] & OPT_INCLUDE_MAJORANT) == 0) ||
              (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT) == 0))
            {
              /* Get next material */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Loop over composition list to try and find this ZAI */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {

              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Compare ZAI to all included ZAIs */

              /* Get pointer to ZAI list */

              ZAIlist = (long)RDB[inuc0 + ITER_NUCLIDE_ZAI];

              /* Smallest acceptable ZAI is 10010 */
              while ((ZAI = (long)RDB[ZAIlist]) > 10000)
                {
                  if ((long)RDB[nuc + NUCLIDE_ZAI] == ZAI)
                    break;

                  /* Compare next ZAI */

                  ZAIlist++;
                }

              /* Check if we found a match */

              if (ZAI > 10000)
                break;

              iso = NextItem(iso);
            }

          /* Check if we found a match */

          if (iso > VALID_PTR)
            {
              fprintf(outp, "Including material %s\n",
                      GetText(mat + MATERIAL_PTR_NAME));

              /* Found a match, include this material */

              materials[N] = mat;
              N++;
            }

          mat = NextItem(mat);

        }

      /***********************/
      /* Store material list */
      /***********************/

      loc0 = ReallocMem(DATA_ARRAY, N + 1);

      /* Set pointer */

      WDB[inuc0 + ITER_NUCLIDE_PTR_MATERIAL_LIST] = (double)loc0;

      /* Read materials */

      for (i = 0; i < N; i++)
        WDB[loc0++] = materials[i];

      /* Put NULLPTR to terminate list */

      WDB[loc0] = NULLPTR;

      /* Free temporary array of materials */

      Mem(MEM_FREE, materials);

      fprintf(outp, "\n");
    }

  /***************************************************************************/

  /***** Get unique nuclides for critical concentration calculation **********/

  /* Get pointer to first (dummy) block*/

  inuc0 = (long)RDB[DATA_PTR_ITER_NUC0];

  /* Reset number of iterable nuclides found */

  N = 0;

  /* Get pointer to ZAI list */

  ZAIlist = (long)RDB[inuc0 + ITER_NUCLIDE_ZAI];

  /* Smallest acceptable ZAI is 10010 */

  while ((ZAI = (long)RDB[ZAIlist]) > 10000)
    {

      /* Loop over material list */

      ptr = (long)RDB[inuc0 + ITER_NUCLIDE_PTR_MATERIAL_LIST];
      while ((mat = (long)RDB[ptr]) > VALID_PTR)
        {

          /* Loop over composition list to try and find this ZAI */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {

              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Compare ZAI */

              if ((long)RDB[nuc + NUCLIDE_ZAI] == ZAI)
                {
                  /* Found one, only process nuclides with non-zero atomic densities */
                  /* mix-card use sometimes leaves zero atomic density nuclides */
                  /* into the composition */

                  if (RDB[iso + COMPOSITION_ADENS] != 0.0)
                    {
                      /* Flag material */

                      WDB[mat + MATERIAL_PTR_ITER_ISO_LIST] = (double)YES;

                      /* Find this nuclide from inuc list */
                      /* Start from the one after dummy */

                      inuc = NextItem(inuc0);
                      while (inuc > VALID_PTR)
                        {

                          if ((long)RDB[inuc + ITER_NUCLIDE_PTR_NUCLIDE] == nuc)
                            break;

                          inuc = NextItem(inuc);
                        }

                      /* If not found, create a new block */

                      if (inuc < VALID_PTR)
                        {
                          inuc = NewItem(DATA_PTR_ITER_NUC0,
                                         ITER_NUCLIDE_BLOCK_SIZE);

                          /* Increment number of found iterable nuclides */

                          N++;

                          /* Copy parameter information */

                          WDB[inuc + PARAM_PTR_NAME] = RDB[inuc0 + PARAM_PTR_NAME];
                          WDB[inuc + PARAM_PTR_FNAME] = RDB[inuc0 + PARAM_PTR_FNAME];
                          WDB[inuc + PARAM_LINE] = RDB[inuc0 + PARAM_LINE];

                          /* Put ZAI */

                          WDB[inuc + ITER_NUCLIDE_ZAI] = ZAI;

                          /* Link nuclide */

                          WDB[inuc + ITER_NUCLIDE_PTR_NUCLIDE] = nuc;

                          /* Flag nuclide */

                          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_CRIT_ITER);

                          /* Link material list */

                          WDB[inuc + ITER_NUCLIDE_PTR_MATERIAL_LIST] =
                            RDB[inuc0 + ITER_NUCLIDE_PTR_MATERIAL_LIST];
                        }
                    }
                }

              iso = NextItem(iso);
            }

          /* Next material name */

          ptr++;
        }

      /* Next ZAI */

      ZAIlist++;
    }

  /* Check if we did not find any iterable nuclides */

  if (N == 0)
    Error(inuc0, "Could not find iterable nuclides from materials");

  /* Now we can remove the dummy block */

  RemoveItem(inuc0);

  /***************************************************************************/
  /***** Create nuclide wise composition lists and check atomic density ******/

  /* Get pointer to first iterable nuclide block */

  inuc = (long)RDB[DATA_PTR_ITER_NUC0];

  while (inuc > VALID_PTR)
    {
      /* Allocate memory for a temporary list of nuclide-wise */
      /* composition list for this nuclide */

      compositions = (long *)Mem(MEM_ALLOC, maxmat, sizeof(long));

      /* Reset maximum adens (used in dtmajorant) */

      maxadens = 0.0;

      /* reset number of found compositions for this nuclide */

      ncompfound = 0;

      /* Loop over material list */
      ptr = (long)RDB[inuc0 + ITER_NUCLIDE_PTR_MATERIAL_LIST];
      while ((mat = (long)RDB[ptr]) > VALID_PTR)
        {
          /* Loop over composition list to try and find this ZAI */

          N = 0.0;
          /* Get pointer to nuclide block */

          nuc = (long)RDB[inuc + ITER_NUCLIDE_PTR_NUCLIDE];

          /* Reset atomic density for this nuclide in this material */

          adens = 0.0;

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];

          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Compare nuclide pointer */

              if (nuc == (long)RDB[inuc + ITER_NUCLIDE_PTR_NUCLIDE])
                {
                  /* This material contains a nuclide that is iterated */

                  fprintf(outp, "Nuclide %s in material %s.\n",
                          GetText(nuc + NUCLIDE_PTR_NAME),
                          GetText(mat + MATERIAL_PTR_NAME));

                  /* Add to total density of this nuclide in this material */

                  adens = adens + RDB[iso + COMPOSITION_ADENS];

                  /* Check zero atomic density for iterable nuclide */

                  /***********************************************************/
                  /* Do something to material that contains iterable nuclide */
                  /* and maybe also to the nuclide                           */
                  /***********************************************************/

                  /* Store direct pointer to composition */

                  compositions[ncompfound] = iso;
                  ncompfound++;

                  /* Increment number of nuclides found for this material */

                  N++;

                }

              /* Next */

              iso = NextItem(iso);
            }

          /* If iterable nuclide was not found from the material */

          if (adens == 0.0)
            {
              /* Get pointer to nuclide block */

              nuc = (long)RDB[inuc + ITER_NUCLIDE_PTR_NUCLIDE];

              /* Write an error */

              Warn(FUNCTION_NAME, "Zero atomic density for iterable nuclide \"%s\" "
                    "in material \"%s\". If the material is a mixture, everything "
                    "might still be OK.", GetText(nuc + NUCLIDE_PTR_NAME),
                    GetText(mat + MATERIAL_PTR_NAME));
            }

          /* Check if atomic density at this material is the largest */

          if (adens > maxadens)
            maxadens = adens;

          /* Next material */

          ptr++;
        }

      /********************************************************/
      /* Store composition pointers for this iterable nuclide */
      /********************************************************/

      loc0 = ReallocMem(DATA_ARRAY, ncompfound + 1);

      /* Set pointer */

      WDB[inuc + ITER_NUCLIDE_PTR_COMPOSITION_LIST] = (double)loc0;

      /* Read materials */

      for (i = 0; i < ncompfound; i++)
        WDB[loc0++] = compositions[i];

      /* Put NULLPTR to terminate list */

      WDB[loc0] = NULLPTR;

      /* Free temporary array of materials */

      Mem(MEM_FREE, compositions);

      /**********************************************************/
      /* Store maximum atomic density for this iterable nuclide */
      /**********************************************************/

      WDB[inuc + ITER_NUCLIDE_MAX_ADENS] = maxadens;

      inuc = NextItem(inuc);
    }

  /*****************************************************************/
  /* Loop over materials to create material-wise composition lists */
  /* for iterable nuclides (needed for quick access in EquiXS)     */
  /*****************************************************************/

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Skip materials that are not included in majorant or not physical */

      if ((((long)RDB[mat + MATERIAL_OPTIONS] & OPT_INCLUDE_MAJORANT) == 0) ||
          (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT) == 0) ||
          ((long)RDB[mat + MATERIAL_PTR_ITER_ISO_LIST] == 0))
        {
          /* Get next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /*************************/
      /* Loop over composition */
      /*************************/

      /* Loop over composition list to count iter nuclides */
      N = 0;

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide data */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_CRIT_ITER)
            N++;

          /* Next composition block */

          iso = NextItem(iso);
        }

      /*****************************/
      /* Read and process nuclides */
      /*****************************/

      if (N > 0)
        {
          loc0 = ReallocMem(DATA_ARRAY, N + 1);

          /* Set pointer */

          WDB[mat + MATERIAL_PTR_ITER_ISO_LIST] = (double)loc0;

          /* Read nuclides  */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide data */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Store composition pointer if part of the iteration */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_CRIT_ITER)
                WDB[loc0++] = iso;

              /* Next composition block */

              iso = NextItem(iso);
            }

          /* Put NULLPTR to terminate list */

          WDB[loc0] = NULLPTR;
        }
      else if (RDB[mat + MATERIAL_PTR_ITER_ISO_LIST] == (double)YES)
        Die(FUNCTION_NAME, "Something wrong when linking iter nuclides");

      mat = NextItem(mat);
    }

  /**********************************************/
  /* Link iterable nuclides to extra dtmajorant */
  /**********************************************/

  inuc = (long)RDB[DATA_PTR_ITER_NUC0];

  while (inuc > VALID_PTR)
    {
      if ((nuc = (long)RDB[inuc + ITER_NUCLIDE_PTR_NUCLIDE])
          > VALID_PTR)
        {
          /* Get nuclide maximum atomic density */

          adens = RDB[inuc + ITER_NUCLIDE_MAX_ADENS];

          /* Find nuclide in majorant extra list */

          ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
          while (ptr > VALID_PTR)
            {
              /* Check nuclide pointer */

              if ((long)RDB[ptr + MAJORANT_EXTRA_PTR_NUC] == nuc)
                break;

              /* Next */

              ptr = NextItem(ptr);
            }

          /* Check if found */

          if (ptr < VALID_PTR)
            {
              /* New entry */

              ptr = NewItem(DATA_MAJORANT_PTR_EXTRA_XS, MAJORANT_EXTRA_BLOCK_SIZE);

              /* Put pointer and type */

              WDB[ptr + MAJORANT_EXTRA_PTR_NUC] = (double)nuc;
              WDB[ptr + MAJORANT_EXTRA_TYPE] =
                (double)MAJORANT_EXTRA_NUCLIDE_ITER;
              WDB[ptr + MAJORANT_EXTRA_FRAC] = adens;
            }

        }

      inuc = NextItem(inuc);
    }

  /* Print newline */

  fprintf(outp, "\n");
}

/*****************************************************************************/
