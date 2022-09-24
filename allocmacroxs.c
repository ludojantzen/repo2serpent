/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : allocmacroxs.c                                 */
/*                                                                           */
/* Created:       2011/06/21 (JLe)                                           */
/* Last modified: 2018/03/27 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Allocates memory for reaction lists and material totals      */
/*                                                                           */
/* Comments: - Lists are filled in processreactionlists.c and cross section  */
/*             calculated in materialtotas.c                                 */
/*                                                                           */
/*           - List generation option is checked in newrealist.c             */
/*                                                                           */
/*           - Rutiinia muutettu radikaalisti 2.11.2011 (2.0.37)             */
/*                                                                           */
/*           - Tää ei nyt luo ures-listaa majorantille, eli oletetaan että   */
/*             Serpent 1:n käyttämää päivitysmoodia ei implementoida         */
/*                                                                           */
/*           - Selvitä jossain välissä että voisiko noille fissioreaktioille */
/*             käyttää yhtä ja samaa listaa (ja onko siitä hyötyä)           */
/*                                                                           */
/*           - Light element production added 11.8.2017 / 2.1.30 / JLe       */
/*                                                                           */
/*           - TODO:  + Muistipaikan pointteri talteen MPI-jakoa varten      */
/*                    + Muuttujien nimet "mode" ja "lst" on epäloogiset      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AllocMacroXS:"

/*****************************************************************************/

void AllocMacroXS()
{
  long mat, mat0, sz, mode, loc0, loc1, n, lst, rea, ptr, erg, ne, nr;
  double Emin, Emax, mem;

  /* Check decay only mode */

  if ((long)RDB[DATA_BURN_DECAY_CALC] == YES)
    return;

  fprintf(outp, "Allocating memory for macroscopic cross section data...\n");

  /***************************************************************************/

  /***** Estimate and pre-allocate memory ************************************/

  /* Check pointer to unionized grid */

  if ((ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID]) > VALID_PTR)
    {
      /* Get number of energy points */
          
      ne = (long)RDB[ptr + ENERGY_GRID_NE];

      /* Memory allocated for majorant */

      sz = ne;

      /* Memory allocated for material-wise cross sections */

      if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES)
        {
          /* Loop over materials */

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

              /* Exclude divisors */
              
              if ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
                {
                  /* Check CE TMS mode and fissile flag */
                  
                  if ((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_CE)
                    sz = sz + ne;
                  if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT)
                    sz = sz + 6*ne;
                  else
                    sz = sz + 4*ne;
                }

              /* Next material */
              
              mat = NextItem(mat);
            }
        }

      /* Pre-allocate memory */

      PreallocMem(sz, DATA_ARRAY);
    }

  /***************************************************************************/

  /***** Create reaction lists for materials *********************************/

  /* NOTE: Tää luo sen RLS-perusrakenteen kaikille materiaaleille, mutta  */
  /* jättää jaettujen materiaalien RLS_DATA -rakenteet pois. Ne linkataan */
  /* parent-materiaaleista mukaan seuraavassa silmukassa. */

  /* Loop over materials */

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

      /* Calculate data size */
  
      CalculateBytes();

      /* Get memory size */
          
      mem = RDB[DATA_TOTAL_BYTES];

      /* Check that composition exists */

      if ((long)RDB[mat + MATERIAL_PTR_COMP] < VALID_PTR)
        Die(FUNCTION_NAME, "Composition not defined");
      
      /* Avoid compiler warning */

      mode = -1;

      /* Loop over lists */

      for (n = 0; n < 20; n++)
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
            mode = MATERIAL_PTR_PROTP_REA_LIST;
          else if (n == 7)
            mode = MATERIAL_PTR_DEUTP_REA_LIST;
          else if (n == 8)
            mode = MATERIAL_PTR_TRITP_REA_LIST;
          else if (n == 9)
            mode = MATERIAL_PTR_HE3P_REA_LIST;
          else if (n == 10)
            mode = MATERIAL_PTR_HE4P_REA_LIST;
          else if (n == 11)
            mode = MATERIAL_PTR_INLP_REA_LIST;
          else if (n == 12)
            mode = MATERIAL_PTR_PHOT_TOT_LIST;
          else if (n == 13)
            mode = MATERIAL_PTR_PHOT_HEAT_LIST;
          else if (n == 14)
            mode = MATERIAL_PTR_TOT_URES_LIST;
          else if (n == 15)
            mode = MATERIAL_PTR_ABS_URES_LIST;
          else if (n == 16)
            mode = MATERIAL_PTR_ELA_URES_LIST;
          else if (n == 17)
            mode = MATERIAL_PTR_FISS_URES_LIST;
          else if (n == 18)
            mode = MATERIAL_PTR_HEAT_URES_LIST;
          else if (n == 19)
            mode = MATERIAL_PTR_TMP_MAJORANT_LIST;
          else
            Die(FUNCTION_NAME, "Overflow");

          /* Create list */
          
          NewReaList(mat, mode);
        }

      /* Calculate data size */
  
      CalculateBytes();

      /* Update memory size */
      
      WDB[mat + MATERIAL_MEM_SIZE] = RDB[mat + MATERIAL_MEM_SIZE] 
        + RDB[DATA_TOTAL_BYTES] - mem;

      /* Next material */
      
      mat = NextItem(mat);
    }

  /****************************************************************************/

  /***** Reaction lists for divided materials *********************************/

  /* Loop over materials */

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

      /* Pointer to parent */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        {
          /* Skip material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check that composition exists */

      if ((long)RDB[mat + MATERIAL_PTR_COMP] < VALID_PTR)
        Die(FUNCTION_NAME, "Composition not defined");
      
      /* Avoid compiler warning */

      mode = -1;

      /* Loop over lists */

      for (n = 0; n < 20; n++)
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
            mode = MATERIAL_PTR_PROTP_REA_LIST;
          else if (n == 7)
            mode = MATERIAL_PTR_DEUTP_REA_LIST;
          else if (n == 8)
            mode = MATERIAL_PTR_TRITP_REA_LIST;
          else if (n == 9)
            mode = MATERIAL_PTR_HE3P_REA_LIST;
          else if (n == 10)
            mode = MATERIAL_PTR_HE4P_REA_LIST;
          else if (n == 11)
            mode = MATERIAL_PTR_INLP_REA_LIST;
          else if (n == 12)
            mode = MATERIAL_PTR_PHOT_TOT_LIST;
          else if (n == 13)
            mode = MATERIAL_PTR_PHOT_HEAT_LIST;
          else if (n == 14)
            mode = MATERIAL_PTR_TOT_URES_LIST;
          else if (n == 15)
            mode = MATERIAL_PTR_ABS_URES_LIST;
          else if (n == 16)
            mode = MATERIAL_PTR_ELA_URES_LIST;
          else if (n == 17)
            mode = MATERIAL_PTR_FISS_URES_LIST;
          else if (n == 18)
            mode = MATERIAL_PTR_HEAT_URES_LIST;
          else if (n == 19)
            mode = MATERIAL_PTR_TMP_MAJORANT_LIST;
          else
            Die(FUNCTION_NAME, "Overflow");

          /* Check if list was erased */

          if ((loc0 = (long)RDB[mat0 + mode]) < VALID_PTR)
            WDB[mat + mode] = NULLPTR;
          else
            {
              /* Get pointer to list */

              loc1 = (long)RDB[mat + mode];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
              
              /* Copy data pointer */

              WDB[loc1 + RLS_PTR_REA0] = RDB[loc0 + RLS_PTR_REA0];
            }
        }
      
      /* Next material */
      
      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Allocate memory for reaction structures *****************************/

  /* Loop over materials */

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

      /* Check divisor type */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
        {
          /* Skip material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Calculate data size */
  
      CalculateBytes();

      /* Get memory size */
          
      mem = RDB[DATA_TOTAL_BYTES];
      
      /* Loop over reaction modes */
      
      for (n = 0; n < 17; n++)
        {
          /* Avoid compiler warning */

          mode = -1;
          lst = -1;

          /* Get mode and list pointer */
          
          if (n == 0)
            {
              mode = MATERIAL_PTR_TOTXS;
              lst = MATERIAL_PTR_TOT_REA_LIST;
            }
          else if (n == 1)
            {
              mode = MATERIAL_PTR_ABSXS;
              lst = MATERIAL_PTR_ABS_REA_LIST;
            }
          else if (n == 2)
            {
              mode = MATERIAL_PTR_ELAXS;
              lst = MATERIAL_PTR_ELA_REA_LIST;
            }
          else if (n == 3)
            {
              mode = MATERIAL_PTR_FISSXS;
              lst = MATERIAL_PTR_FISS_REA_LIST;
            }
          else if (n == 4)
            {
              mode = MATERIAL_PTR_HEATTXS;
              lst = MATERIAL_PTR_HEATT_REA_LIST;
            }
          else if (n == 5)
            {
              mode = MATERIAL_PTR_PHOTPXS;
              lst = MATERIAL_PTR_PHOTP_REA_LIST;
            }
          else if (n == 6)
            {
              mode = MATERIAL_PTR_PROTPXS;
              lst = MATERIAL_PTR_PROTP_REA_LIST;
            }
          else if (n == 7)
            {
              mode = MATERIAL_PTR_DEUTPXS;
              lst = MATERIAL_PTR_DEUTP_REA_LIST;
            }
          else if (n == 8)
            {
              mode = MATERIAL_PTR_TRITPXS;
              lst = MATERIAL_PTR_TRITP_REA_LIST;
            }
          else if (n == 9)
            {
              mode = MATERIAL_PTR_HE3PXS;
              lst = MATERIAL_PTR_HE3P_REA_LIST;
            }
          else if (n == 10)
            {
              mode = MATERIAL_PTR_HE4PXS;
              lst = MATERIAL_PTR_HE4P_REA_LIST;
            }
          else if (n == 11)
            {
              mode = MATERIAL_PTR_INLPXS;
              lst = MATERIAL_PTR_INLP_REA_LIST;
            }
          else if (n == 12)
            {
              mode = MATERIAL_PTR_FISSE;
              lst = MATERIAL_PTR_FISS_REA_LIST;
            }
          else if (n == 13)
            {
              mode = MATERIAL_PTR_NSF;
              lst = MATERIAL_PTR_FISS_REA_LIST;
            }
          else if (n == 14)
            {
              mode = MATERIAL_PTR_TOTPHOTXS;
              lst = MATERIAL_PTR_PHOT_TOT_LIST;
            }
          else if (n == 15)
            {
              mode = MATERIAL_PTR_HEATPHOTXS;
              lst = MATERIAL_PTR_PHOT_HEAT_LIST;
            }
          else if (n == 16)
            {
              mode = MATERIAL_PTR_TMP_MAJORANTXS;
              lst = MATERIAL_PTR_TMP_MAJORANT_LIST;
            }
          else
            Die (FUNCTION_NAME, "Invalid reaction mode");
          
          /* Check that reaction pointer doesn't exist */

          if ((long)RDB[mat + mode] > VALID_PTR)
            Die(FUNCTION_NAME, "Reaction already exists");

          /* Cycle loop if nuclide has no reactions of this type */

          if ((long)RDB[mat + lst] < VALID_PTR)
            continue;

          /* Allocate memory for block */

          rea = NewItem(mat + mode, REACTION_BLOCK_SIZE);
          
          /* Put type */
              
          WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_SUM;
              
          /* Put mt */
              
          if (mode == MATERIAL_PTR_TOTXS)
            WDB[rea + REACTION_MT] = MT_MACRO_TOTXS;
          else if (mode == MATERIAL_PTR_ABSXS)
            WDB[rea + REACTION_MT] = MT_MACRO_ABSXS;
          else if (mode == MATERIAL_PTR_ELAXS)
            WDB[rea + REACTION_MT] = MT_MACRO_ELAXS;
          else if (mode == MATERIAL_PTR_FISSXS)
            WDB[rea + REACTION_MT] = MT_MACRO_FISSXS;
          else if (mode == MATERIAL_PTR_HEATTXS)
            WDB[rea + REACTION_MT] = MT_MACRO_HEATXS;
          else if (mode == MATERIAL_PTR_PHOTPXS)
            WDB[rea + REACTION_MT] = MT_MACRO_PHOTXS;
          else if (mode == MATERIAL_PTR_PROTPXS)
            WDB[rea + REACTION_MT] = MT_MACRO_PROTPXS;
          else if (mode == MATERIAL_PTR_DEUTPXS)
            WDB[rea + REACTION_MT] = MT_MACRO_DEUTPXS;
          else if (mode == MATERIAL_PTR_TRITPXS)
            WDB[rea + REACTION_MT] = MT_MACRO_TRITPXS;
          else if (mode == MATERIAL_PTR_HE3PXS)
            WDB[rea + REACTION_MT] = MT_MACRO_HE3PXS;
          else if (mode == MATERIAL_PTR_HE4PXS)
            WDB[rea + REACTION_MT] = MT_MACRO_HE4PXS;
          else if (mode == MATERIAL_PTR_INLPXS)
            WDB[rea + REACTION_MT] = MT_MACRO_INLPRODXS;
          else if (mode == MATERIAL_PTR_FISSE)
            WDB[rea + REACTION_MT] = MT_MACRO_FISSE;
          else if (mode == MATERIAL_PTR_NSF)
            WDB[rea + REACTION_MT] = MT_MACRO_NSF;
          else if (mode == MATERIAL_PTR_TOTPHOTXS)
            WDB[rea + REACTION_MT] = MT_MACRO_TOTPHOTXS;
          else if (mode == MATERIAL_PTR_HEATPHOTXS)
            WDB[rea + REACTION_MT] = MT_MACRO_HEATPHOTXS;
          else if (mode == MATERIAL_PTR_TMP_MAJORANTXS)
            WDB[rea + REACTION_MT] = MT_MACRO_TMP_MAJORANTXS;
          else
            Die (FUNCTION_NAME, "Invalid reaction mode");

          /* Put mode */

          WDB[rea + REACTION_MODE] = (double)lst;
              
          /* Put material pointer */
          
          WDB[rea + REACTION_PTR_MAT] = (double)mat;
          
          /* Put pointer to partial list */
          
          WDB[rea + REACTION_PTR_PARTIAL_LIST] = RDB[mat + lst];
          
          /* Reset minimum and maximum energy */
          
          WDB[rea + REACTION_EMIN] = INFTY;
          WDB[rea + REACTION_EMAX] = -INFTY;
              
          /* Reset ures energy boundaries */
          
          WDB[rea + REACTION_URES_EMIN] = INFTY;
          WDB[rea + REACTION_URES_EMAX] = -INFTY;
          
          /* Reset pointer to energy grid */
          
          WDB[rea + REACTION_PTR_EGRID] = NULLPTR;
          
          /* Reset pointer to xs data */

          WDB[rea + REACTION_PTR_XS] = NULLPTR;
          
          /* Reset first point and number of points */
          
          WDB[rea + REACTION_XS_I0] = -1.0;
          WDB[rea + REACTION_XS_NE] = -1.0;
          
          /* Allocate memory for previous value */
              
          AllocValuePair(rea + REACTION_PTR_PREV_XS);
        }
      
      /* Calculate number of bytes */

      CalculateBytes();

      /* Store beginning of data block */

      WDB[mat + MATERIAL_PTR_DATA_BLOCK] = RDB[DATA_ALLOC_MAIN_SIZE];
    
      /***********************************************************************/
      
      /***** Allocate memory for multi-group total cross sections ************/
      
      if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
        {
          /* Get pointer to total cross section */

          rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Get number of points */

          ne = (long)RDB[DATA_COARSE_MG_NE];
          CheckValue(FUNCTION_NAME, "np", "", ne, 10, 50000);

          /* Allocate memory for data */

          ptr = ReallocMem(DATA_ARRAY, ne);

          /* Put pointer */

          WDB[rea + REACTION_PTR_MGXS] = (double)ptr;
        }

      /* Check if macroscopic cross sections are reconstructed */

      if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == NO)
        {
          /* Calculate number of bytes */

          CalculateBytes();
          
          /* Update memory size */
          
          WDB[mat + MATERIAL_MEM_SIZE] = RDB[mat + MATERIAL_MEM_SIZE] 
            + RDB[DATA_TOTAL_BYTES] - mem;

          /* Store data block size */

          WDB[mat + MATERIAL_DATA_BLOCK_SIZE] = 
            RDB[DATA_ALLOC_MAIN_SIZE] - WDB[mat + MATERIAL_PTR_DATA_BLOCK];

          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /***********************************************************************/
      
      /***** Allocate memory for pre-calculated neutron cross sections *******/

      /* Get pointer to unionized neutron energy grid */
  
      if ((erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID]) > VALID_PTR)
        {
          /* Get number of energy points */
          
          ne = (long)RDB[erg + ENERGY_GRID_NE];
          
          /* Get minimum and maximum energy */
          
          Emin = RDB[erg + ENERGY_GRID_EMIN];
          Emax = RDB[erg + ENERGY_GRID_EMAX];
        }
      else
        {
          /* Avoid compiler warning */
          
          ne = -1;
          Emin = INFTY;
          Emax = -INFTY;
        }

      /* Set number of modes */

      if ((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_CE)
        nr = 1;
      else
        nr = 15;

      /* Loop over reaction modes */
      
      for (n = 0; n < nr; n++)
        {
          /* Avoid compiler warning */

          mode = -1;
          lst = -1;

          /* Get mode and list pointer */

          if (n == 0)
            {
              mode = MATERIAL_PTR_TMP_MAJORANTXS;
              lst = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANT_LIST];
            }
          else if (n == 1)
            {
              mode = MATERIAL_PTR_TOTXS;
              lst = (long)RDB[mat + MATERIAL_PTR_TOT_REA_LIST];
            }
          else if (n == 2)
            {
              mode = MATERIAL_PTR_ABSXS;
              lst = (long)RDB[mat + MATERIAL_PTR_ABS_REA_LIST];
            }
          else if (n == 3)
            {
              mode = MATERIAL_PTR_ELAXS;
              lst = (long)RDB[mat + MATERIAL_PTR_ELA_REA_LIST];
            }
          else if (n == 4)
            {
              mode = MATERIAL_PTR_FISSXS;
              lst = (long)RDB[mat + MATERIAL_PTR_FISS_REA_LIST];
            }
          else if (n == 5)
            {
              mode = MATERIAL_PTR_HEATTXS;
              lst = (long)RDB[mat + MATERIAL_PTR_HEATT_REA_LIST];
            }
          else if (n == 6)
            {
              mode = MATERIAL_PTR_PHOTPXS;
              lst = (long)RDB[mat + MATERIAL_PTR_PHOTP_REA_LIST];
            }
          else if (n == 7)
            {
              mode = MATERIAL_PTR_PROTPXS;
              lst = (long)RDB[mat + MATERIAL_PTR_PROTP_REA_LIST];
            }
          else if (n == 8)
            {
              mode = MATERIAL_PTR_DEUTPXS;
              lst = (long)RDB[mat + MATERIAL_PTR_DEUTP_REA_LIST];
            }
          else if (n == 9)
            {
              mode = MATERIAL_PTR_TRITPXS;
              lst = (long)RDB[mat + MATERIAL_PTR_TRITP_REA_LIST];
            }
          else if (n == 10)
            {
              mode = MATERIAL_PTR_HE3PXS;
              lst = (long)RDB[mat + MATERIAL_PTR_HE3P_REA_LIST];
            }
          else if (n == 11)
            {
              mode = MATERIAL_PTR_HE4PXS;
              lst = (long)RDB[mat + MATERIAL_PTR_HE4P_REA_LIST];
            }
          else if (n == 12)
            {
              mode = MATERIAL_PTR_INLPXS;
              lst = (long)RDB[mat + MATERIAL_PTR_INLP_REA_LIST];
            }
          else if (n == 13)
            {
              mode = MATERIAL_PTR_FISSE;
              lst = (long)RDB[mat + MATERIAL_PTR_FISS_REA_LIST];
            }
          else if (n == 14)
            {
              mode = MATERIAL_PTR_NSF;
              lst = (long)RDB[mat + MATERIAL_PTR_FISS_REA_LIST];
            }
          else
            Die (FUNCTION_NAME, "Invalid reaction mode");

          /* Cycle loop if list is not defined */
              
          if (lst < VALID_PTR)
            continue;

          /* Check unionized grid pointer */

          if (erg < VALID_PTR)
            Die(FUNCTION_NAME, "Reconstruction without unionized grid (1)");

          /* Pointer to reaction data */

          rea = (long)RDB[mat + mode];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Put minimum and maximum energy */
          
          WDB[rea + REACTION_EMIN] = Emin;
          WDB[rea + REACTION_EMAX] = Emax;
              
          /* Put pointer to energy grid */
          
          WDB[rea + REACTION_PTR_EGRID] = (double)erg;
          
          /* Put first point and number of points */
          
          WDB[rea + REACTION_XS_I0] = 0.0;
          WDB[rea + REACTION_XS_NE] = (double)ne;
          
          /* Allocate memory */

          ptr = ReallocMem(DATA_ARRAY, ne);
          WDB[rea + REACTION_PTR_XS] = (double)ptr;
        }
      
      /***********************************************************************/

      /***** Allocate memory for pre-calculated photon cross sections ********/

      /* Get pointer to unionized neutron energy grid */
  
      if ((erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_PGRID]) > VALID_PTR)
        {
          /* Get number of energy points */
          
          ne = (long)RDB[erg + ENERGY_GRID_NE];
          
          /* Get minimum and maximum energy */
          
          Emin = RDB[erg + ENERGY_GRID_EMIN];
          Emax = RDB[erg + ENERGY_GRID_EMAX];
        }
      else
        {
          /* Avoid compiler warning */
          
          ne = -1;
          Emin = INFTY;
          Emax = -INFTY;
        }
      
      /* Loop over reaction modes */
      
      for (n = 0; n < 2; n++)
        {
          /* Avoid compiler warning */

          mode = -1;
          lst = -1;
          
          /* Get mode */
          
          if (n == 0)
            {
              mode = MATERIAL_PTR_TOTPHOTXS;
              lst = (long)RDB[mat + MATERIAL_PTR_PHOT_TOT_LIST];
            }
          else if (n == 1)
            {
              mode = MATERIAL_PTR_HEATPHOTXS;
              lst = (long)RDB[mat + MATERIAL_PTR_PHOT_HEAT_LIST];
            }
          else
            Die (FUNCTION_NAME, "Invalid reaction mode");
      
          /* Cycle loop if list is not defined */
              
          if (lst < VALID_PTR)
            continue;

          /* Check unionized grid pointer */

          if (erg < VALID_PTR)
            Die(FUNCTION_NAME, "Reconstruction without unionized grid (2)");
          
          /* Pointer to reaction data */

          rea = (long)RDB[mat + mode];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Put minimum and maximum energy */
          
          WDB[rea + REACTION_EMIN] = Emin;
          WDB[rea + REACTION_EMAX] = Emax;
              
          /* Put pointer to energy grid */
          
          WDB[rea + REACTION_PTR_EGRID] = (double)erg;
          
          /* Put first point and number of points */
          
          WDB[rea + REACTION_XS_I0] = 0.0;
          WDB[rea + REACTION_XS_NE] = (double)ne;
          
          /* Allocate memory */
          
          ptr = ReallocMem(DATA_ARRAY, ne);
          WDB[rea + REACTION_PTR_XS] = (double)ptr;
        }
      
      /***********************************************************************/

      /* Calculate number of bytes */

      CalculateBytes();

      /* Update memory size */
      
      WDB[mat + MATERIAL_MEM_SIZE] = RDB[mat + MATERIAL_MEM_SIZE] 
        + RDB[DATA_TOTAL_BYTES] - mem;

      /* Store data block size */

      WDB[mat + MATERIAL_DATA_BLOCK_SIZE] = 
        RDB[DATA_ALLOC_MAIN_SIZE] - WDB[mat + MATERIAL_PTR_DATA_BLOCK];
      
      /* Next material */
      
      mat = NextItem(mat);
    }

  /* Update counter */

  WDB[DATA_TOT_MAT_BYTES] = RDB[DATA_TOT_MAT_BYTES] + (double)MemCount();
      
  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
  
  /***** Allocate memory for DT neutron majorant *****************************/

  /* Check DT flag */

  if ((long)RDB[DATA_OPT_USE_DT] == YES)
    {        
      /* Check number of nuclides */
      
      if ((long)RDB[DATA_N_TRANSPORT_NUCLIDES] > 0)
        {
          /* Get pointer to unionized neutron energy grid */
  
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
          CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

          /* Create reaction block */

          rea = NewItem(DATA_PTR_MAJORANT, REACTION_BLOCK_SIZE);
        
          /* Put mt */
          
          WDB[rea + REACTION_MT] = (double)MT_MACRO_MAJORANT;
          
          /* Allocate memory for previous value */
          
          AllocValuePair(rea + REACTION_PTR_PREV_XS);
                    
          /* Copy energy grid pointer */
          
          WDB[rea + REACTION_PTR_EGRID] = (double)erg;
          
          /* Get number of energy points */
          
          ne = (long)RDB[erg + ENERGY_GRID_NE];
          
          /* Put number of energy points */
          
          WDB[rea + REACTION_XS_NE] = (double)ne;
          
          /* Put minimum and maximum energy */
          
          WDB[rea + REACTION_EMIN] = RDB[erg + ENERGY_GRID_EMIN];
          WDB[rea + REACTION_EMAX] = RDB[erg + ENERGY_GRID_EMAX];
          
          /* Allocate memory for data */
          
          ptr = ReallocMem(DATA_ARRAY, ne);
          WDB[rea + REACTION_PTR_MAJORANT_XS] = (double)ptr;
          
          /* Allocate memory for coarse multi-group majorants */
          
          if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
            {
              /* Get number of groups */
              
              ne = (long)RDB[DATA_COARSE_MG_NE];
              CheckValue(FUNCTION_NAME, "np", "", ne, 10, 50000);
              
              /* Allocate memory for data */
              
              ptr = ReallocMem(DATA_ARRAY, ne);
              
              /* Put pointer */
              
              WDB[rea + REACTION_PTR_MGXS] = (double)ptr;
            }
        }
    }

  /***************************************************************************/

  /***** Allocate memory for DT photon majorant ******************************/

  /* Check DT flag */

  if ((long)RDB[DATA_OPT_USE_DT] == YES)
    {        
      /* Check number of nuclides */
      
      if ((long)RDB[DATA_N_PHOTON_NUCLIDES] > 0)
        {
          /* Get pointer to unionized neutron energy grid */
      
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_PGRID];
          CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

          /* Create reaction block */

          rea = NewItem(DATA_PTR_PHOTON_MAJORANT, REACTION_BLOCK_SIZE);
        
          /* Put mt */
          
          WDB[rea + REACTION_MT] = (double)MT_MACRO_MAJORANT;
          
          /* Allocate memory for previous value */
          
          AllocValuePair(rea + REACTION_PTR_PREV_XS);
          
          /* Copy energy grid pointer */
          
          WDB[rea + REACTION_PTR_EGRID] = (double)erg;
          
          /* Get number of energy points */
          
          ne = (long)RDB[erg + ENERGY_GRID_NE];
          
          /* Put number of energy points */
          
          WDB[rea + REACTION_XS_NE] = (double)ne;
          
          /* Put minimum and maximum energy */
          
          WDB[rea + REACTION_EMIN] = RDB[erg + ENERGY_GRID_EMIN];
          WDB[rea + REACTION_EMAX] = RDB[erg + ENERGY_GRID_EMAX];
          
          /* Allocate memory for data */
          
          ptr = ReallocMem(DATA_ARRAY, ne);
          WDB[rea + REACTION_PTR_MAJORANT_XS] = (double)ptr;
        }
    }

  /* Update counter (count majorants in xs) */
  
  WDB[DATA_TOT_XS_BYTES] = RDB[DATA_TOT_XS_BYTES] + (double)MemCount();

  /***************************************************************************/
}

/*****************************************************************************/
