/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdatainterfaces.c                        */
/*                                                                           */
/* Created:       2018/06/04 (VVa)                                           */
/* Last modified: 2018/06/05 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Prepares general data-fields that can be updated             */
/*                                                                           */
/* Comments:  -                                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessDataInterfaces:"

/*****************************************************************************/

void ProcessDataInterfaces()
{
  long ptr, msh, mat, loc0, arr, iso, nuc;
  long n0, n1, n2, n, ZAI;

  if ((loc0 = (long)RDB[DATA_PTR_DATAIFC0]) < VALID_PTR)
    return;

  fprintf(outp, "Processing general data-interfaces...\n\n");

  /* Define variables to appease the compiler */

  nuc = 0;

  while (loc0 > VALID_PTR)
    {

      fprintf(outp, "Data-interface with data file \"%s\":\n",
              GetText(loc0 + DATAIFC_PTR_INPUT_FNAME));

      /******************************/
      /* Find and link the datamesh */
      /******************************/

      /* Get pointer to mesh name */

      ptr = (long)(loc0 + DATAIFC_PTR_MESH);

      /* Find this material from global materials */

      msh = (long)RDB[DATA_PTR_MESH0];
      if ((msh = SeekListStr(msh, MESH_PTR_NAME, GetText(ptr)))
          < VALID_PTR)
        Die(FUNCTION_NAME, "Unable to find mesh \"%s\" linked to "
            "data interface.", GetText(ptr));

      /* Get mesh size */

      n0 = (long)RDB[msh + MESH_N0];
      n1 = (long)RDB[msh + MESH_N1];
      n2 = (long)RDB[msh + MESH_N2];

      /* Print some output for the data mesh */

      fprintf(outp, " - Linked to data mesh %s with %ld x %ld x %ld bins (%ld bins total)\n",
              GetText(loc0 + DATAIFC_PTR_MESH), n0, n1, n2, n0*n1*n2);

      /* Link the mesh to interface */

      WDB[loc0 + DATAIFC_PTR_MESH] = (double)msh;

      /* Store data-memory size to interface (for MPI transfer) */

      WDB[loc0 + DATAIFC_DATA_MEM_SIZE] = (double)(n0*n1*n2);

      /* Allocate space for interface data */

      WDB[loc0 + DATAIFC_PTR_DATA] = (double)ReallocMem(DATA_ARRAY, n0*n1*n2);

      /**************************************************/
      /* Find and link the materials for this datafield */
      /**************************************************/

      if ((ptr = (long)RDB[loc0 + DATAIFC_PTR_MAT_ARR]) > VALID_PTR)
        {
          /* Loop over material list */

          while ((long)RDB[ptr] > VALID_PTR)
            {
              /* Find this material from global materials */

              mat = (long)RDB[DATA_PTR_M0];
              if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, GetText(ptr)))
                  < VALID_PTR)
                Die(FUNCTION_NAME, "Unable to find material \"%s\" linked to "
                    "a data interface."
                    , GetText(ptr));

              /* Swap material name pointer for material pointer */

              WDB[ptr] = (double)mat;

              /* Print something for this material */

              fprintf(outp, " - Provides data for material %s\n",
                      GetText(mat + MATERIAL_PTR_NAME));

              /* Add to number of data interfaces for the material  */
              /* Linking will be done after all data interfaces are */
              /* processed */

              WDB[mat + MATERIAL_PTR_DATAIFC_ARR]++;

              /* Next material */

              ptr++;
            }
        }

      /********************************************/
      /* Data type specific processing in the end */
      /********************************************/

      if ((long)RDB[loc0 + DATAIFC_DATA_TYPE] == DATA_TYPE_IFC_ADENS)
        {
          /************************************************/
          /* Find and link the nuclide for this datafield */
          /************************************************/

          ZAI = (long)RDB[loc0 + DATAIFC_PTR_NUCLIDE];

          /* Loop over material array */

          if ((arr = (long)RDB[loc0 + DATAIFC_PTR_MAT_ARR]) < VALID_PTR)
            Error(loc0, "Need to link at least one material to dataifc "
                  "providing adens data");

          while ((mat = (long)RDB[arr]) > VALID_PTR)
            {

              /* Loop over composition list to try and find this ZAI */

              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

              while (iso > VALID_PTR)
                {

                  /* Pointer to nuclide data */

                  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                  /* Compare material ZAI to ZAI to be linked */

                  if ((long)RDB[nuc + NUCLIDE_ZAI] == ZAI)
                    break;

                  /* Next nuclide from composition */

                  iso = NextItem(iso);
                }

              /* Check if nuclide was found */

              if (iso < VALID_PTR)
                Error(loc0, "Could not find ZAI %ld from material %s linked to a data interface",
                      ZAI, GetText(mat + MATERIAL_PTR_NAME));

              /* This should be available */

              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Print out information for the linked nuclide */

              fprintf(outp, " - Provides atomic density data for %s in material %s\n",
                      GetText(nuc + NUCLIDE_PTR_NAME), GetText(mat + MATERIAL_PTR_NAME));

              /* Flag material */

              SetOption(mat + MATERIAL_OPTIONS, OPT_EXT_ADENS_MAT);

              /* Link nuclide */

              WDB[loc0 + DATAIFC_PTR_NUCLIDE] = (double)nuc;

              /* Flag nuclide */

              SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_EXT_ADENS);

              /* Check critical adens iteration flag from nuclide */

              if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_CRIT_ITER)
                Error(loc0, "Nuclide %s linked to receive densities from a "
                      "data interface is also included in a critical atomic density "
                      "iteration", GetText(nuc + NUCLIDE_PTR_NAME));

              /* Next material from array */

              arr++;
            }
        }

      /* Next data interface */

      loc0 = NextItem(loc0);
    }

  /*************************************/
  /* Link data interfaces to materials */
  /*************************************/

  /* Loop over all materials */

  mat = (long)RDB[DATA_PTR_M0];

  while (mat > VALID_PTR)
    {
      /* If the material was not linked to any data interface, skip it here */

      if ((n = (long)RDB[mat + MATERIAL_PTR_DATAIFC_ARR]) == 0)
        {
          mat = NextItem(mat);

          continue;
        }

      /* Reserve memory for the NULLPTR terminated array */

      arr = (double)ReallocMem(DATA_ARRAY, n + 1);

      /* Store the array to the material */

      WDB[mat + MATERIAL_PTR_DATAIFC_ARR] = (double)arr;

      /* Terminate the array */

      WDB[arr + n] = NULLPTR;

      /* Decrement n */

      n--;

      /* Loop over the data interfaces to find and link the n interfaces */

      loc0 = (long)RDB[DATA_PTR_DATAIFC0];

      while (loc0 > VALID_PTR)
        {
          /* Skip data interfaces with no material list */

          if ((ptr = (long)RDB[loc0 + DATAIFC_PTR_MAT_ARR]) < VALID_PTR)
            {
              loc0 = NextItem(loc0);

              continue;
            }

          /* Loop over material array */

          while ((long)RDB[ptr] > VALID_PTR)
            {
              if ((long)RDB[ptr] == mat)
                {
                  /* Found a match, link it */

                  WDB[arr + n] = (double)loc0;

                  /* Decrement n */

                  n--;
                }

              /* Next material from array */

              ptr++;
            }

          /* Next data interface */

          loc0 = NextItem(loc0);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /**********************************************/
  /* Link external nuclides to extra dtmajorant */
  /* This only does the initial linking         */
  /**********************************************/

  loc0 = (long)RDB[DATA_PTR_DATAIFC0];

  while (loc0 > VALID_PTR)
    {
      if ((long)RDB[loc0 + DATAIFC_DATA_TYPE] == DATA_TYPE_IFC_ADENS)
        {
          if ((nuc = (long)RDB[loc0 + DATAIFC_PTR_NUCLIDE]) > VALID_PTR)
            {
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
                    (double)MAJORANT_EXTRA_NUCLIDE_EXT;
                  WDB[ptr + MAJORANT_EXTRA_FRAC] = 0.0;
                }
            }
        }

      loc0 = NextItem(loc0);
    }

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
