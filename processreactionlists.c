/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processreactionlists.c                         */
/*                                                                           */
/* Created:       2012/05/28 (JLe)                                           */
/* Last modified: 2019/03/01 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Prepares reaction lists for transport calculation            */
/*                                                                           */
/* Comments: - Re-written from scratch 28.5.2012 (2.1.6)                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessReactionLists:"

/* Use local function to simplify OpenMP implementation */

void ProcessReactionLists0(long);

/*****************************************************************************/

void ProcessReactionLists()
{
  long mat, mat0, nuc, rea, iso, imax, n, mode, loc0, loc1, idx, nused, nacti;
  long ptr;
  double f, Emin, Emax;

  /***************************************************************************/

  /***** Calculate maximum densities *****************************************/

  fprintf(outp, "Calculating maximum densities...\n");

  /* Put densities for undivided mateials and reset values for divided */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check division */

      if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
        {
          /* Skip material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check that material is handled in all domains */

      if (((long)RDB[DATA_DD_DECOMPOSE] == YES) && 
          ((long)RDB[mat + MATERIAL_MPI_ID] > -1))
        Die(FUNCTION_NAME, "Not included in domain");          

      /* Get pointer to composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Get maximum index */

      imax = ListSize(iso) - 1;

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

          /* Get pointer to list */

          if ((loc0 = (long)RDB[mat + mode]) < VALID_PTR)
            continue;
          
          /* Loop over list*/

          loc1 = (long)RDB[loc0 + RLS_PTR_REA0];
          while (loc1 > VALID_PTR)
            {
              /* Get composition index */

              idx = (long)RDB[loc1 + RLS_DATA_COMP_IDX];

              /* Check */

              if ((idx < 0) || (idx > imax))
                Die(FUNCTION_NAME, "Error in composition index");

              /* Get pointer to data */

              iso = ListPtr(iso, idx);
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

              /* Check material div type and set density */

              if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
                WDB[loc1 + RLS_DATA_MAX_ADENS] = 0.0;
              else
                WDB[loc1 + RLS_DATA_MAX_ADENS] = RDB[iso + COMPOSITION_ADENS];

              /* Next */

              loc1 = NextItem(loc1);
            }
        }          

      /* Next material */
      
      mat = NextItem(mat);
    }          

  /* Densities for divided */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check division */

      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
        {
          /* Skip material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

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

      /* Get pointer to composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Get maximum index */

      imax = ListSize(iso) - 1;

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

          /* Get pointer to list */

          if ((loc0 = (long)RDB[mat0 + mode]) < VALID_PTR)
            continue;

          /* Loop over list */

          loc1 = (long)RDB[loc0 + RLS_PTR_REA0];
          while (loc1 > VALID_PTR)
            {
              /* Get composition index */

              idx = (long)RDB[loc1 + RLS_DATA_COMP_IDX];

              /* Check */

              if ((idx < 0) || (idx > imax))
                Die(FUNCTION_NAME, "Error in composition index");

              /* Get pointer to data */

              iso = ListPtr(iso, idx);
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
              
              /* Compare density */

              if (RDB[iso + COMPOSITION_ADENS] > RDB[loc1 + RLS_DATA_MAX_ADENS])
                WDB[loc1 + RLS_DATA_MAX_ADENS] = RDB[iso + COMPOSITION_ADENS];

              /* Next */

              loc1 = NextItem(loc1);
            }
        }          

      /* Next material */
      
      mat = NextItem(mat);
    }          

  fprintf(outp, "OK.\n\n");

  /* Share data with other domains */

  DistributeDDMaxAdens();

  /***************************************************************************/

  /***** Density cut-off *****************************************************/

  fprintf(outp, "Performing density cut-off...\n");

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check division */

      if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
        {
          /* Skip material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check that material is handled in all domains */

      if (((long)RDB[DATA_DD_DECOMPOSE] == YES) && 
          ((long)RDB[mat + MATERIAL_MPI_ID] > -1))
        Die(FUNCTION_NAME, "Not included in domain");

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

          /* Get pointer to list */

          if ((loc0 = (long)RDB[mat + mode]) < VALID_PTR)
            continue;

          /* Loop over data and set cut-off flags */

          loc1 = (long)RDB[loc0 + RLS_PTR_REA0];
          while (loc1 > VALID_PTR)
            {
              /* Reset flag */

              WDB[loc1 + RLS_DATA_CUT] = (double)NO;

              /* Pointer to reaction */

              rea = (long)RDB[loc1 + RLS_DATA_PTR_REA];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Pointer to nuclide */

              nuc = (long)RDB[loc1 + RLS_DATA_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Density cut-off */

              if (RDB[loc1 + RLS_DATA_MAX_ADENS]*RDB[nuc + NUCLIDE_MAX_TOTXS] 
                  < RDB[mat + MATERIAL_ADENS]*RDB[DATA_MIN_TOTXS])
                WDB[loc1 + RLS_DATA_CUT] = (double)YES;

              /* Remove cut-off flags for equilibrium poisons */

              if (((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == YES) &&
                  ((long)RDB[nuc + NUCLIDE_ZAI] == 541350))
                WDB[loc1 + RLS_DATA_CUT] = (double)NO;

              if (((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == YES) &&
                  ((long)RDB[nuc + NUCLIDE_ZAI] == 621490))
                WDB[loc1 + RLS_DATA_CUT] = (double)NO;

              /* Pointer to OTF burn ZAI list */

              if ((ptr = (long)RDB[DATA_OTF_BURN_NUC_LIST]) > VALID_PTR)
                {
                  /* Loop over ZAI list to find match */

                  while (RDB[ptr] > 0.0)
                    {              
                      /* Compare */
              
                      if (RDB[nuc + NUCLIDE_ZAI] == RDB[ptr])
                        {
                          /* Reset cut flag */
                          
                          WDB[loc1 + RLS_DATA_CUT] = (double)NO;

                          /* Break loop */

                          break;
                        }
                      else
                        ptr++;
                    }
                }

              /* Next */
              
              loc1 = NextItem(loc1);
            }
        }          
      
      /* Next material */
      
      mat = NextItem(mat);
    }

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
  
  /***** Ures dilution cut-off ***********************************************/
  
  if (((long)RDB[DATA_USE_URES] == YES) && ((long)RDB[DATA_URES_USED] > 0))
    {
      fprintf(outp,
              "Performing ures dilution cut-off with %1.1E threshold...\n",
              RDB[DATA_URES_DILU_CUT]);
      
      /* Set sampling flags on nuclides */
      
      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Check if ures data is used */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED)
            {
              /* Sanity check on energy boundaries */

              if (RDB[nuc + NUCLIDE_URES_EMIN] < ZERO)
                Die(FUNCTION_NAME, "Error in ures Emin");
              if (RDB[nuc + NUCLIDE_URES_EMAX] < RDB[nuc + NUCLIDE_URES_EMIN])
                Die(FUNCTION_NAME, "Error in ures Emax");
              
              /* Get maximum atomic fraction */
              
              f = RDB[nuc + NUCLIDE_MAX_AFRAC];
              CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
              
              /* Test cut-off and set sampling flag on nuclide */
              
              if (f < RDB[DATA_URES_DILU_CUT])
                WDB[nuc + NUCLIDE_URES_SAMPLING] = (double)NO;            
              else
                WDB[nuc + NUCLIDE_URES_SAMPLING] = (double)YES;
            }

          /* Next nuclide */
          
          nuc = NextItem(nuc);
        }
      
      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Reset ures boundaries for material */

          WDB[mat + MATERIAL_URES_EMIN] = INFTY;
          WDB[mat + MATERIAL_URES_EMAX] = -INFTY;

          /* Check division */

          if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
            {
              /* Skip material */

              mat = NextItem(mat);
              
              /* Cycle loop */
              
              continue;
            }

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

          /* Avoid compiler warning */

          mode = -1;

          /* Loop over lists */

          for (n = 0; n < 5; n++)
            {
              /* Get mode */
              
              if (n == 0)
                mode = MATERIAL_PTR_TOT_URES_LIST;
              else if (n == 1)
                mode = MATERIAL_PTR_ABS_URES_LIST;
              else if (n == 2)
                mode = MATERIAL_PTR_ELA_URES_LIST;
              else if (n == 3)
                mode = MATERIAL_PTR_FISS_URES_LIST;
              else if (n == 4)
                mode = MATERIAL_PTR_HEAT_URES_LIST;
              else
                Die(FUNCTION_NAME, "Overflow");
              
              /* Get pointer to list */
              
              if ((loc0 = (long)RDB[mat + mode]) < VALID_PTR)
                continue;

              /* Loop over data and do cut-offs */

              loc1 = (long)RDB[loc0 + RLS_PTR_REA0];
              while (loc1 > VALID_PTR)
                {
                  /* Check if already cut by density */

                  if ((long)RDB[loc1 + RLS_DATA_CUT] == NO)
                    {
                      /* Pointer to nuclide */
                      
                      nuc = (long)RDB[loc1 + RLS_DATA_PTR_NUCLIDE];
                      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
                      
                      /* Check sampling flag and do cut-off */

                      if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO)
                        WDB[loc1 + RLS_DATA_CUT] = YES;
                      else
                        {
                          /* Compare energy boundaries to ures limits */
          
                          if (RDB[nuc + NUCLIDE_URES_EMIN] <
                              RDB[mat + MATERIAL_URES_EMIN])
                            WDB[mat + MATERIAL_URES_EMIN] =
                              RDB[nuc + NUCLIDE_URES_EMIN];

                          if (RDB[nuc + NUCLIDE_URES_EMAX] >
                              RDB[mat + MATERIAL_URES_EMAX])
                            WDB[mat + MATERIAL_URES_EMAX] =
                              RDB[nuc + NUCLIDE_URES_EMAX];
                        }
                    }
                  
                  /* Next */
                  
                  loc1 = NextItem(loc1);
                }
            }
          
          /* Next material */
          
          mat = NextItem(mat);
        }

      /* Set material- and reaction-wise ures boundaries for divided */
      /* materials */

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

          /* Check division */
          
          if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
            {
              /* Copy material-wise boundaries from parent */

              WDB[mat + MATERIAL_URES_EMIN] = RDB[mat0 + MATERIAL_URES_EMIN];
              WDB[mat + MATERIAL_URES_EMAX] = RDB[mat0 + MATERIAL_URES_EMAX];
            }

          /* Avoid compiler warning */

          rea = -1;

          /* Loop over lists */

          for (n = 0; n < 7; n++)
            {
              /* Get mode (HUOM! indeksi testataan myös tuolla alempana) */
              
              if (n == 0)
                rea = MATERIAL_PTR_TOTXS;
              else if (n == 1)
                rea = MATERIAL_PTR_ABSXS;
              else if (n == 2)
                rea = MATERIAL_PTR_ELAXS;
              else if (n == 3)
                rea = MATERIAL_PTR_FISSXS;
              else if (n == 4)
                rea = MATERIAL_PTR_HEATTXS;
              else if (n == 5)
                rea = MATERIAL_PTR_FISSE;
              else if (n == 6)
                rea = MATERIAL_PTR_NSF;
              else
                Die(FUNCTION_NAME, "Overflow");
              
              /* Get pointer to reaction data */
              
              if ((rea = (long)RDB[mat + rea]) < VALID_PTR)
                continue;

              /* Reset boundaries */

              Emin = INFTY;
              Emax = -INFTY;
              
              /* Check for fission channels */

              if ((n == 3) || (n == 5) || (n == 6))
                {
                  /* Find boundaries for fissile isotopes (NOTE: tämä    */
                  /* siksi että materiaalissa voi olla ures-isotooppeja  */
                  /* ja fissiilejä isotooppeja joilla ei ole ures-dataa. */

                  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
                  while (iso > VALID_PTR)
                    {
                      /* Pointer to nuclide data */

                      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                      /* Check sampling flag and fission pointer */
                      
                      if (((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES) &&
                          ((long)RDB[nuc + NUCLIDE_PTR_FISSXS] > VALID_PTR))
                        {
                          /* Compare boundaries */

                          if (RDB[nuc + NUCLIDE_URES_EMIN] < Emin)
                            Emin = RDB[nuc + NUCLIDE_URES_EMIN];
                          
                          if (RDB[nuc + NUCLIDE_URES_EMAX > Emax])
                            Emax = RDB[nuc + NUCLIDE_URES_EMAX];
                        }

                      /* Next */

                      iso = NextItem(iso);
                    }
                }
              else
                {
                  /* Use material-wise values */

                  Emin = RDB[mat + MATERIAL_URES_EMIN];
                  Emax = RDB[mat + MATERIAL_URES_EMAX];
                }

              /* Put boundaries */

              WDB[rea + REACTION_URES_EMIN] = Emin;
              WDB[rea + REACTION_URES_EMAX] = Emax;
            }
          
          /* Next material */
          
          mat = NextItem(mat);
        }

      /* Reset global ures energy boundaries */

      WDB[DATA_URES_EMIN] = INFTY;
      WDB[DATA_URES_EMAX] = -INFTY;

      /* Count number of nuclides and set boundaries */
  
      nused = 0;
      nacti = 0;
  
      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Check sampling flag */
          
          if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES)
            {
              /* Compare to global boundaries */

              if (RDB[nuc + NUCLIDE_URES_EMIN] < RDB[DATA_URES_EMIN])
                WDB[DATA_URES_EMIN] = RDB[nuc + NUCLIDE_URES_EMIN];
                  
              if (RDB[nuc + NUCLIDE_URES_EMAX] > RDB[DATA_URES_EMAX])
                WDB[DATA_URES_EMAX] = RDB[nuc + NUCLIDE_URES_EMAX];
 
              nacti++;
            }

          /* Check if data is used */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED)
            nused++;
          
          /* Next nuclide */

          nuc = NextItem(nuc);
        }
      
      /* Print */

      if (nacti == 1)
        fprintf(outp, "\n - Ures ptable data used with One nuclide:\n\n");
      else if (nacti > 0)
        fprintf(outp, "\n - Ures ptable data used with %ld nuclides:\n\n", 
                nacti);
  
      /* Loop over nuclides */
      
      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Check flag */
          
          if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES)
            fprintf(outp, "%12s (from %1.2E to %1.2E MeV)\n", 
                    GetText(nuc + NUCLIDE_PTR_NAME), 
                    RDB[nuc + NUCLIDE_URES_EMIN], 
                    RDB[nuc + NUCLIDE_URES_EMAX]);
      
          /* Next nuclide */
          
          nuc = NextItem(nuc);
        }
      
      /* Print number of excluded nuclides */

      if (nused - nacti == 1)
        fprintf(outp, 
                "\n - Ures ptable data excluded in one nuclide.\n");
      else if ((nacti == 0) && (nused > 0))
        fprintf(outp, "\n - All ures ptables excluded by cut-offs.\n");
      else if (nused - nacti > 0)
        fprintf(outp,
                "\n - Ures ptable data excluded in %ld nuclides.\n",
                nused - nacti);
      else
        fprintf(outp, "\n");

      /* Print global boundaries */

      if (nused > 0)
        fprintf(outp, 
                " - Sampling used from %1.2E to %1.2E MeV\n\n",
                RDB[DATA_URES_EMIN], RDB[DATA_URES_EMAX]);
      else
        fprintf(outp, "\n");

      /* Set count */
      
      WDB[DATA_URES_USED] = (double)nacti;

      fprintf(outp, "OK.\n\n");
    }

  /***************************************************************************/

  /***** Ures-energy boundaries in depletion lists ***************************/

  /* Check spectrum-collapse method */
  
  if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
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
          
          /* Transmutation list */

          loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];
          while (loc0 > VALID_PTR)
            {
              /* Pointer to reaction */

              rea = (long)RDB[loc0 + DEP_TRA_PTR_REA];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Pointer to nuclide */

              nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check sampling flag and set energy */

              if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES)
                WDB[loc0 + DEP_TRA_E0] = RDB[rea + REACTION_URES_EMIN];
              else
                WDB[loc0 + DEP_TRA_E0] = INFTY;

              /* Next */

              loc0 = NextItem(loc0);
            }

          /* Fission list */

          loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];
          while (loc0 > VALID_PTR)
            {
              /* Pointer to reaction */

              rea = (long)RDB[loc0 + DEP_TRA_PTR_REA];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Pointer to nuclide */

              nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check sampling flag and set energy */

              if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES)
                WDB[loc0 + DEP_TRA_E0] = RDB[rea + REACTION_URES_EMIN];
              else
                WDB[loc0 + DEP_TRA_E0] = INFTY;

              /* Next */

              loc0 = NextItem(loc0);
            }

          /* Fission neutron production list */

          loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_NSF_LIST];
          while (loc0 > VALID_PTR)
            {
              /* Pointer to reaction */

              rea = (long)RDB[loc0 + DEP_TRA_PTR_REA];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Pointer to nuclide */

              nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check sampling flag and set energy */

              if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES)
                WDB[loc0 + DEP_TRA_E0] = RDB[rea + REACTION_URES_EMIN];
              else
                WDB[loc0 + DEP_TRA_E0] = INFTY;

              /* Next */

              loc0 = NextItem(loc0);
            }

          /* Next material */

          mat = NextItem(mat);
        }
    }

  /***************************************************************************/

  /***** Sort lists **********************************************************/

  fprintf(outp, "Sorting material-wise reaction lists:\n\n");

  /* Reset thread numbers */
  
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check thread number */
      
      WDB[mat + MATERIAL_OMP_ID] = -1.0;
      
      /* Next material */
      
      mat = NextItem(mat);
    }          

  /* Start parallel timer */

  StartTimer(TIMER_OMP_PARA);
    
#ifdef OPEN_MP
#pragma omp parallel private (mat)
#endif
  {
    /* Print */

    PrintProgress(0, 0);
    
    /* Loop over materials */
    
    mat = (long)RDB[DATA_PTR_M0];
    while (mat > VALID_PTR)
      {
        /* Test parallel id's (reaction lists are not included in the */
        /* MPI distributed data block, so the sameprocessing is done  */
        /* for all tasks. */

        if (MyParallelMat(mat, NO) == YES)
          {
            /* Process */
            
            ProcessReactionLists0(mat);
            
            /* Print */
            
            PrintProgress(mat, 1);
          }
        
        /* Next material */
        
        mat = NextItem(mat);
      }
  }
  
  /* Stop parallel timer */

  StopTimer(TIMER_OMP_PARA);

  /* Print */

  PrintProgress(0, 100);

  /***************************************************************************/

  /***** Reconfigure pointers for divided materials ***************************/

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

          /* Check pointer */

          if ((loc0 = (long)RDB[mat0 + mode]) < VALID_PTR)
            continue;
        
          /* Get pointer to list */

          loc1 = (long)RDB[mat + mode];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
              
          /* Copy data pointer */
          
          WDB[loc1 + RLS_PTR_REA0] = RDB[loc0 + RLS_PTR_REA0];
        }
      
      /* Next material */
      
      mat = NextItem(mat);
    }

  /***************************************************************************/
}

/*****************************************************************************/

void ProcessReactionLists0(long mat)
{
  long n, mode, loc0, loc1, nmax;

  /* Avoid compiler warning */

  mode = -1;

  /* Sort only transmutation and fission lists for divided materials */

  if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
    nmax = 3;
  else
    nmax = 23;
  
  /* Loop over lists */
  
  for (n = 0; n < nmax; n++)
    {
      /* Get mode */

      if (n == 0)
        mode = MATERIAL_PTR_DEP_TRA_LIST;
      else if (n == 1)
        mode = MATERIAL_PTR_DEP_FISS_LIST;
      else if (n == 2)
        mode = MATERIAL_PTR_DEP_NSF_LIST;
      else if (n == 3)
        mode = MATERIAL_PTR_TOT_REA_LIST;
      else if (n == 4)
        mode = MATERIAL_PTR_ELA_REA_LIST;
      else if (n == 5)
        mode = MATERIAL_PTR_ABS_REA_LIST;
      else if (n == 6)
        mode = MATERIAL_PTR_FISS_REA_LIST;
      else if (n == 7)
        mode = MATERIAL_PTR_HEATT_REA_LIST;
      else if (n == 8)
        mode = MATERIAL_PTR_PHOTP_REA_LIST;
      else if (n == 9)
        mode = MATERIAL_PTR_PROTP_REA_LIST;
      else if (n == 10)
        mode = MATERIAL_PTR_DEUTP_REA_LIST;
      else if (n == 11)
        mode = MATERIAL_PTR_TRITP_REA_LIST;
      else if (n == 12)
        mode = MATERIAL_PTR_HE3P_REA_LIST;
      else if (n == 13)
        mode = MATERIAL_PTR_HE4P_REA_LIST;
      else if (n == 14)
        mode = MATERIAL_PTR_INLP_REA_LIST;
      else if (n == 15)
        mode = MATERIAL_PTR_PHOT_TOT_LIST;
      else if (n == 16)
        mode = MATERIAL_PTR_PHOT_HEAT_LIST;
      else if (n == 17)
        mode = MATERIAL_PTR_TOT_URES_LIST;
      else if (n == 18)
        mode = MATERIAL_PTR_ABS_URES_LIST;
      else if (n == 19)
        mode = MATERIAL_PTR_ELA_URES_LIST;
      else if (n == 20)
        mode = MATERIAL_PTR_FISS_URES_LIST;
      else if (n == 21)
        mode = MATERIAL_PTR_HEAT_URES_LIST;
      else if (n == 22)
        mode = MATERIAL_PTR_TMP_MAJORANT_LIST;
      else
        Die(FUNCTION_NAME, "Overflow");
      
      /* Get pointer to list */
      
      if ((loc0 = (long)RDB[mat + mode]) < VALID_PTR)
        continue;
      
      /* Fission and transmutation lists */

      if ((mode == MATERIAL_PTR_DEP_TRA_LIST) || 
          (mode == MATERIAL_PTR_DEP_FISS_LIST) ||
          (mode == MATERIAL_PTR_DEP_NSF_LIST))
        {
          /* Sort */
          
          SortList(loc0, DEP_TRA_E0, SORT_MODE_ASCEND);
        }
      else
        {
          /* Get pointer to data */
      
          loc1 = (long)RDB[loc0 + RLS_PTR_REA0];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
          
          /* Sort total list by maximum density and others by energy */
          
          if (mode == MATERIAL_PTR_TOT_REA_LIST)
            SortList(loc1, RLS_DATA_MAX_ADENS, SORT_MODE_DESCEND);
          else if (mode == MATERIAL_PTR_TMP_MAJORANT_LIST)
            SortList(loc1, RLS_DATA_MAX_ADENS, SORT_MODE_DESCEND);
          else
            SortList(loc1, RLS_DATA_EMIN, SORT_MODE_ASCEND);
      
          /* Sort everything by cut-off flag */
          
          SortList(loc1, RLS_DATA_CUT, SORT_MODE_ASCEND);
        }
    }
}

/*****************************************************************************/
