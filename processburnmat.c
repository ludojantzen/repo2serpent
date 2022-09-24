/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processburnmat.c                               */
/*                                                                           */
/* Created:       2011/01/25 (JLe)                                           */
/* Last modified: 2020/03/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Additional processing for materials used in burnup           */
/*              calculation                                                  */
/*                                                                           */
/* Comments: - Mass flow lists are reconstructed with nuclide pointers       */
/*             to simplify and speed-up MakeBurnMatrixMSR().                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessBurnMat:"

/*****************************************************************************/

void ProcessBurnMat()
{
  long mat, ptr, loc0, loc1, loc2, ne, mat0, mat1, sz, ZAI, iso0, iso1;
  long nuc0, nuc1;
  double mem;

  /* Check burnup mode */

 if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
   return;

  /* Check burnup step */

 if (((long)RDB[DATA_BURN_STEP] > 0) &&
     ((long)RDB[DATA_PTR_HISV0]) < VALID_PTR)
   Die(FUNCTION_NAME, "Should not be here");

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Reset burn-flags for non-physical materials (this is needed for */
      /* detector and source materials) */

      if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        ResetOption(mat + MATERIAL_OPTIONS, OPT_BURN_MAT);

      /* Next material */

      mat = NextItem(mat);
    }

  /* Print */

  fprintf(outp, "Allocating memory for burnup calculation...\n");

  /***************************************************************************/

  /***** Allocate memory from RES2 block to speed up calculation *************/

  /* Calculate memory size */

  sz = 0;

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

      /* Check burn flag and pointer to parent */

      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
          ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR))
        {
          /* Add sum of spectrum */

          sz = sz + 1;

          /* Add transmutation cross sections */

          if ((ptr = (long)RDB[mat0 + MATERIAL_PTR_DEP_TRA_LIST]) > VALID_PTR)
            sz = sz + ListSize(ptr);

          /* Add fission cross sections */

          if ((ptr = (long)RDB[mat0 + MATERIAL_PTR_DEP_FISS_LIST]) > VALID_PTR)
            sz = sz + ListSize(ptr);

          /* Add fission neutron production cross sections */

          if ((ptr = (long)RDB[mat0 + MATERIAL_PTR_DEP_NSF_LIST]) > VALID_PTR)
            sz = sz + ListSize(ptr);

          /* Add flux spectrum */

          if (((long)RDB[DATA_BURN_DECAY_CALC] == NO) &&
              ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES))
            {
              /* Pointer to unionized energy grid */

              ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Add number of points */

              sz = sz + (long)RDB[ptr + ENERGY_GRID_NE];
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Allocate memory */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    PreallocMem(sz, RES3_ARRAY);
  else
    PreallocMem(sz, RES2_ARRAY);

  /***************************************************************************/

  /***** Process material-wise data ******************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /***********************************************************************/

      /***** Check flags, etc. ***********************************************/

      /* Allocate memory for burn flux (NOTE: This must be here because of */
      /* domain decomposition) */

      ptr = NewStat("BURN_FLUX", 1, 1);
      WDB[mat + MATERIAL_PTR_BURN_FLUX] = (double)ptr;

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

      /* Check composition */

      if ((long)RDB[mat + MATERIAL_PTR_COMP] < VALID_PTR)
        Die(FUNCTION_NAME, "Material %s has no nuclides",
            GetText(mat + MATERIAL_PTR_NAME));

      /* Check mixture */

      if ((long)RDB[mat + MATERIAL_PTR_MIX] > VALID_PTR)
        Die(FUNCTION_NAME, "Burnable mixture %s?",
            GetText(mat + MATERIAL_PTR_NAME));

      /* Check div type */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
        Die(FUNCTION_NAME, "Parent material");

      /* Put output flag if not parent or divided */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NONE)
        WDB[mat + MATERIAL_BURN_PRINT_OUTPUT] = (double)YES;

      /* Calculate bytes */

      CalculateBytes();

      /* Get memory size */

      mem = RDB[DATA_TOTAL_BYTES];

      /***********************************************************************/

      /***** Copy transmutation and fission lists from parent ****************/

      /* Get pointer to parent material */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
        {
          /* Copy doppler mode */

          WDB[mat + MATERIAL_TMS_MODE] = RDB[mat0 + MATERIAL_TMS_MODE];

          /* Reset pointer */

          loc1 = -1;

          /* Loop over transmutation list */

          loc0 = (long)RDB[mat0 + MATERIAL_PTR_DEP_TRA_LIST];
          while (loc0 > VALID_PTR)
            {
              /* Allocate memory */

              loc1  = NewItem(mat + MATERIAL_PTR_DEP_TRA_LIST,
                              DEP_TRA_BLOCK_SIZE);

              /* Copy values */

              WDB[loc1 + DEP_TRA_PTR_REA] = RDB[loc0 + DEP_TRA_PTR_REA];
              WDB[loc1 + DEP_TRA_E0] = RDB[loc0 + DEP_TRA_E0];

              /* Next */

              loc0 = NextItem(loc0);
            }

          /* Close list */

          if (loc1 > VALID_PTR)
            CloseList(loc1);

          /* Reset pointer */

          loc1 = -1;

          /* Loop over fission list */

          loc0 = (long)RDB[mat0 + MATERIAL_PTR_DEP_FISS_LIST];
          while (loc0 > VALID_PTR)
            {
              /* Allocate memory */

              loc1  = NewItem(mat + MATERIAL_PTR_DEP_FISS_LIST,
                              DEP_TRA_BLOCK_SIZE);

              /* Copy values */

              WDB[loc1 + DEP_TRA_PTR_REA] = RDB[loc0 + DEP_TRA_PTR_REA];
              WDB[loc1 + DEP_TRA_E0] = RDB[loc0 + DEP_TRA_E0];

              /* Next */

              loc0 = NextItem(loc0);
            }

          /* Close list */

          if (loc1 > VALID_PTR)
            CloseList(loc1);

          /* Reset pointer */

          loc1 = -1;

          /* Loop over fission neutron production list */

          loc0 = (long)RDB[mat0 + MATERIAL_PTR_DEP_NSF_LIST];
          while (loc0 > VALID_PTR)
            {
              /* Allocate memory */

              loc1  = NewItem(mat + MATERIAL_PTR_DEP_NSF_LIST,
                              DEP_TRA_BLOCK_SIZE);

              /* Copy values */

              WDB[loc1 + DEP_TRA_PTR_REA] = RDB[loc0 + DEP_TRA_PTR_REA];
              WDB[loc1 + DEP_TRA_E0] = RDB[loc0 + DEP_TRA_E0];

              /* Next */

              loc0 = NextItem(loc0);
            }

          /* Close list */

          if (loc1 > VALID_PTR)
            CloseList(loc1);
        }

      /************************************************************************/

      /***** Allocate memory for results **************************************/

      /* Transmutation cross sections */

      loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Allocate memory */

          if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
            ptr = AllocPrivateData(1, RES3_ARRAY);
          else
            ptr = AllocPrivateData(1, RES2_ARRAY);

          WDB[loc0 + DEP_TRA_PTR_RESU] = (double)ptr;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Fission cross sections */

      loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Allocate memory */

          if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
            ptr = AllocPrivateData(1, RES3_ARRAY);
          else
            ptr = AllocPrivateData(1, RES2_ARRAY);

          WDB[loc0 + DEP_TRA_PTR_RESU] = (double)ptr;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Fission neutron production cross sections */

      loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_NSF_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Allocate memory */

          if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
            ptr = AllocPrivateData(1, RES3_ARRAY);
          else
            ptr = AllocPrivateData(1, RES2_ARRAY);

          WDB[loc0 + DEP_TRA_PTR_RESU] = (double)ptr;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Allocate memory for flux spectrum */

      if (((long)RDB[DATA_BURN_DECAY_CALC] == NO) &&
          ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES))
        {
          /* Pointer to unionized grid */

          if ((ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID]) < VALID_PTR)
            Die(FUNCTION_NAME, "Energy grid not unionized");

          /* Get number of points */

          ne = (long)RDB[ptr + ENERGY_GRID_NE];

          /* Allocate memory */

          if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
            ptr = AllocPrivateData(ne, RES3_ARRAY);
          else
            ptr = AllocPrivateData(ne, RES2_ARRAY);

          WDB[mat + MATERIAL_PTR_FLUX_SPEC] = (double)ptr;
        }

      /* Allocate memory for sum of spectrum */

      if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
        ptr = AllocPrivateData(1, RES3_ARRAY);
      else
        ptr = AllocPrivateData(1, RES2_ARRAY);

      WDB[mat + MATERIAL_PTR_FLUX_SPEC_SUM] = (double)ptr;

      /************************************************************************/

      /* Calculate bytes */

      CalculateBytes();

      /* Update memory size */

      WDB[mat + MATERIAL_MEM_SIZE] = RDB[mat + MATERIAL_MEM_SIZE] +
        RDB[DATA_TOTAL_BYTES] - mem;

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Process flow lists **************************************************/

  /* Loop materials */

  mat0 = (long)RDB[DATA_PTR_M0];
  while (mat0 > VALID_PTR)
    {
      /* Loop over outflow */

      loc0 = (long)RDB[mat0 + MATERIAL_PTR_OUTFLOW];
      while (loc0 > VALID_PTR)
        {
          /* Check if domain decomposition is in use */

          if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
            Die(FUNCTION_NAME, "No way this is going to work with DD!");

          /* Check pointer */

          if ((long)RDB[loc0 + REPROC_CON_PTR_MAT1] != mat0)
            Die(FUNCTION_NAME, "Pointer error");

          /* Get second material pointer */

          mat1 = (long)RDB[loc0 + REPROC_CON_PTR_MAT2];
          CheckPointer(FUNCTION_NAME, "(mat1)", DATA_ARRAY, mat1);

          /* Pointer to flow */

          loc1 = (long)RDB[loc0 + REPROC_CON_PTR_MFLOW];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Reset pointer */

          WDB[loc0 + REPROC_CON_PTR_MFLOW] = 0.0;

          /* Loop over data */

          loc2 = (long)RDB[loc1 + MFLOW_PTR_DATA];
          while (loc2 > VALID_PTR)
            {
              /* Get ZAI */

              ZAI = (long)RDB[loc2 + MFLOW_LIST_ZAI];

              /* Pointer to data */

              iso0 = (long)RDB[mat0 + MATERIAL_PTR_COMP];
              CheckPointer(FUNCTION_NAME, "(iso0)", DATA_ARRAY, iso0);

              /* Skip lost and loop over source composition */

              iso0 = NextItem(iso0);
              while (iso0 > VALID_PTR)
                {
                  nuc0 = (long)RDB[iso0 + COMPOSITION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc0)", DATA_ARRAY, nuc0);

                  /* Avoid compiler warning */

                  nuc1 = -1;

                  /* Compare */

                  if ((ZAI == -1) || (ZAI == (long)RDB[nuc0 + NUCLIDE_ZAI]) ||
                      ((ZAI < 1000) && (ZAI == (long)RDB[nuc0 + NUCLIDE_Z])))
                    {
                      /* Find target nuclide */

                      iso1 = (long)RDB[mat1 + MATERIAL_PTR_COMP];
                      while (iso1 > VALID_PTR)
                        {
                          nuc1 = (long)RDB[iso1 + COMPOSITION_PTR_NUCLIDE];
                          CheckPointer(FUNCTION_NAME, "(nuc1)", DATA_ARRAY,
                                       nuc1);

                          /* Compare */

                          if ((long)RDB[nuc1 + NUCLIDE_ZAI] ==
                              (long)RDB[nuc0 + NUCLIDE_ZAI])
                            break;

                          /* Next */

                          iso1 = NextItem(iso1);
                        }

                      /* Check */

                      if (iso1 < VALID_PTR)
                        Die(FUNCTION_NAME, "Nuclide %ld not found in %s",
                            ZAI, GetText(mat1 + MATERIAL_PTR_NAME));

                      /* Create new item */

                      ptr = NewItem(loc0 + REPROC_CON_PTR_MFLOW,
                                    MFLOW_LIST_BLOCK_SIZE);

                      /* Copy data */

                      memcpy(&WDB[ptr + LIST_DATA_SIZE],
                             &RDB[loc2 + LIST_DATA_SIZE],
                             (MFLOW_LIST_BLOCK_SIZE - LIST_DATA_SIZE)*
                             sizeof(double));

                      /* Put pointer to original (used for calculating */
                      /* total atomic density) */

                      WDB[ptr + MFLOW_LIST_PTR_ORIG] = (double)loc2;

                      /* Check pointers one more time */

                      CheckPointer(FUNCTION_NAME, "(iso0)", DATA_ARRAY, iso0);
                      CheckPointer(FUNCTION_NAME, "(iso1)", DATA_ARRAY, iso1);

                      /* Put nuclide pointers */

                      WDB[ptr + MFLOW_LIST_PTR_ISO0] = (double)iso0;
                      WDB[ptr + MFLOW_LIST_PTR_ISO1] = (double)iso1;
                    }

                  /* Next */

                  iso0 = NextItem(iso0);
                }

              /* Next */

              loc2 = NextItem(loc2);
            }

          /* Pointer to next */

          loc0 = NextItem(loc0);
        }

      /* Pointer to next */

      mat0 = NextItem(mat0);
    }

  /***************************************************************************/

  /* Calculate data size */

  CalculateBytes();

  /* Exit OK */

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
