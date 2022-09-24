/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processotfburn.c                               */
/*                                                                           */
/* Created:       2018/03/25 (JLe)                                           */
/* Last modified: 2018/04/12 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Creates data structures used for the on-the-fly burnup       */
/*              solver.                                                      */
/*                                                                           */
/* Comments: - Tässä on varmaan ongelma tuon NewStat():in kanssa silloin     */
/*             kun domain decomposition on käytössä.                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessOTFBurn:"

void LoopOTFBurnChain(long loc0, long l);

/*****************************************************************************/

void ProcessOTFBurn()
{
  long mat, iso, nuc, loc0, loc1, ptr, n, m;
  char *str;

  /* Check burnup mode */

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
    return;

  /* Check OTF burn mode */

  if ((long)RDB[DATA_OTF_BURN_MODE] == NO)
    return;

  fprintf(outp, "Processing data for on-the-fly burnup calculation...\n");

  /***************************************************************************/
  
  /***** Process nuclide list ************************************************/

  /* Pointer to list */

  if ((ptr = (long)RDB[DATA_OTF_BURN_NUC_LIST]) < VALID_PTR)
    Die(FUNCTION_NAME, "Nuclide list is empty");

  /* Loop over list */

  while ((long)RDB[ptr] > VALID_PTR)
    {
      /* Pointer to text */

      str = &ASCII[(long)RDB[ptr]];

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

      /* Get ZAI */

      if ((n = IsotoZAI(str)) < 0)
        n = atoi(str);
      
      /* Put ZAI in list */

      if (n > 0)
        WDB[ptr] = n;
      else        
        Error(0, "Error in OTF burn nuclide list: \"%s\" is not a valid entry",
              str);
      
      /* Next */

      ptr++;
    }

  /***************************************************************************/
  
  /***** Process data ********************************************************/

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
      
      /* Check burn flag and divisor type */

      if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) ||
          ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT))
        {
          /* Pointer to next */
          
          mat = NextItem(mat);

          /* Cycle loop */
          
          continue;
        }

      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Pointer to ZAI list */

          ptr = (long)RDB[DATA_OTF_BURN_NUC_LIST];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over ZAI list to find match */

          while (RDB[ptr] > 0.0)
            {              
              /* Compare */
              
              if (RDB[nuc + NUCLIDE_ZAI] == RDB[ptr])
                break;
              else
                ptr++;
            }
          
          /* Check if found */

          if (RDB[ptr] > 0.0)
            {
              /* Find nuclide in majorant extra list */

              loc0 = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
              while (loc0 > VALID_PTR)
                {
                  /* Check nuclide pointer */
                  
                  if ((long)RDB[loc0 + MAJORANT_EXTRA_PTR_NUC] == nuc)
                    break;
                  
                  /* Next */
                  
                  loc0 = NextItem(loc0);
                }

              /* Check if found */
      
              if ((loc0 < VALID_PTR) && 
                  ((long)RDB[nuc + NUCLIDE_PTR_TOTXS] > VALID_PTR))
                {
                  /* New entry */
                  
                  loc0 = NewItem(DATA_MAJORANT_PTR_EXTRA_XS, 
                                 MAJORANT_EXTRA_BLOCK_SIZE);
                  
                  /* Put pointer and type */
                  
                  WDB[loc0 + MAJORANT_EXTRA_PTR_NUC] = (double)nuc;
                  WDB[loc0 + MAJORANT_EXTRA_TYPE] =
                    (double)MAJORANT_EXTRA_OTF_BURN;
                }

              /* Create structure */

              loc1 = NewItem(DATA_PTR_OTF_BURN0, OTF_BURN_BLOCK_SIZE);

              /* Put pointers */

              WDB[loc1 + OTF_BURN_PTR_COMP] = (double)iso;
              WDB[loc1 + OTF_BURN_PTR_NUCLIDE] = (double)nuc;
              WDB[loc1 + OTF_BURN_PTR_MAT] = (double)mat;
              WDB[loc1 + OTF_BURN_PTR_MAJ] = (double)loc0;

              /* Reset processed flag */

              WDB[loc1 + OTF_BURN_PROCESSED] = (double)NO;

              /* Check if material pointer has been set */

              if ((long)RDB[mat + MATERIAL_PTR_OTF_BURN] < VALID_PTR)
                WDB[mat + MATERIAL_PTR_OTF_BURN] = (double)loc1;
            }
          
          /* Next */

          iso = NextItem(iso);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Process chains ******************************************************/

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
      
      /* Check burn flag and divisor type */

      if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) ||
          ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT))
        {
          /* Pointer to next */
          
          mat = NextItem(mat);

          /* Cycle loop */
          
          continue;
        }

      /* Loop over data */

      loc0 = (long)RDB[mat + MATERIAL_PTR_OTF_BURN];
      while (loc0 > VALID_PTR)
        {
          /* Check material pointer */

          if ((long)RDB[loc0 + OTF_BURN_PTR_MAT] != mat)
            break;

          /* Pointer to nuclide */

          nuc = (long)RDB[loc0 + OTF_BURN_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Loop over chains */

          LoopOTFBurnChain(loc0, 0);

          /* Next */

          loc0 = NextItem(loc0);
        }











      /* Loop over data */

      loc0 = (long)RDB[mat + MATERIAL_PTR_OTF_BURN];
      while (loc0 > VALID_PTR)
        {
          /* Check material pointer */

          if ((long)RDB[loc0 + OTF_BURN_PTR_MAT] != mat)
            break;

          /* Pointer to reactions */

          loc1 = (long)RDB[loc0 + OTF_BURN_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Close list */

          CloseList(loc1);

          /* Sort */
          
          SortList(loc1, OTF_BURN_REA_MT, SORT_MODE_ASCEND);

          /* Pointer to nuclide */
          
          nuc = (long)RDB[loc0 + OTF_BURN_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/

void LoopOTFBurnChain(long loc0, long l)
{
  long mat, nuc, ptr, rea, tgt, yld, loc1, loc2, mt;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Check length */

  if (l > 100)
    Die(FUNCTION_NAME, "Infinite loop");

  /* Check if already processed */

  if ((long)RDB[loc0 + OTF_BURN_PROCESSED] == YES)
    return;
  else
    WDB[loc0 + OTF_BURN_PROCESSED] = (double)YES;

  /* Pointer to material */

  mat = (long)RDB[loc0 + OTF_BURN_PTR_MAT];
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Pointer to nuclide */

  nuc = (long)RDB[loc0 + OTF_BURN_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Allocate memory for total transmutation xs */
  
  ptr = NewStat("OTF_BURN_XS", 1, 1);
  WDB[loc0 + OTF_BURN_PTR_TOT_XS] = (double)ptr;

  /*
  printf("%s ", GetText(mat + MATERIAL_PTR_NAME));
  for (n = 0; n < l; n++)
    printf("  ");
  printf("%s\n", GetText(nuc + NUCLIDE_PTR_NAME));
  */
  /* Loop over reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Get mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Skip branch reactions */

      if (mt > 20000)
        {
          /* Pointer to next reaction */

          rea = NextItem(rea);

          /* Cycle loop */

          continue;
        }

      /* Check target pointer */

      if ((tgt = (long)RDB[rea + REACTION_PTR_TGT]) > VALID_PTR)
        {
          /* Create new reaction */
          
          loc1 = NewItem(loc0 + OTF_BURN_PTR_REA, OTF_BURN_REA_BLOCK_SIZE);
          
          /* Put pointer and mt */
          
          WDB[loc1 + OTF_BURN_REA_PTR_REA] = (double)rea;
          WDB[loc1 + OTF_BURN_REA_MT] = (double)mt;

          /* Find target */
          
          ptr = (long)RDB[mat + MATERIAL_PTR_OTF_BURN];
          while (ptr > VALID_PTR)
            {
              /* Compare pointer */
              
              if (tgt == (long)RDB[ptr + OTF_BURN_PTR_NUCLIDE])
                break;
              
              /* Next */
              
              ptr = NextItem(ptr);
            }
          
          /* Check pointer */
          
          if (ptr > VALID_PTR)
            {
              /* Put target */
              
              WDB[loc1 + OTF_BURN_REA_PTR_TGT] = (double)ptr;

              /* Allocate memory for results */

              if (mt < 10000)
                {
                  loc2 = NewStat("OTF_BURN_XS", 1, 1);
                  WDB[loc1 + OTF_BURN_REA_PTR_XS] = (double)loc2;
                }
              else
                WDB[loc1 + OTF_BURN_REA_PTR_XS] = NULLPTR;

              /*
              printf("%s ", GetText(mat + MATERIAL_PTR_NAME));
              for (n = 0; n < l; n++)
                printf("  ");
              printf("%5ld --> %s\n", (long)RDB[rea + REACTION_MT],
                     GetText(tgt + NUCLIDE_PTR_NAME));
              */
              /* Move forward in chain */
              
              LoopOTFBurnChain(ptr, l + 1);
            }
          else
            {
              /* No target, put pointers */
              
              WDB[loc1 + OTF_BURN_REA_PTR_TGT] = NULLPTR;
              WDB[loc1 + OTF_BURN_REA_PTR_XS] = NULLPTR;
            }
        }
      else if (((mt > 17) && (mt < 22)) || (mt == 38) || (mt == 10006))
        {
          /* Create new reaction */
          
          loc1 = NewItem(loc0 + OTF_BURN_PTR_REA, OTF_BURN_REA_BLOCK_SIZE);
          
          /* Put pointer and mt */
          
          WDB[loc1 + OTF_BURN_REA_PTR_REA] = (double)rea;
          WDB[loc1 + OTF_BURN_REA_MT] = (double)mt;

          /* Allocate memory for results */

          if (mt < 10000)
            {
              loc2 = NewStat("OTF_BURN_XS", 1, 1);
              WDB[loc1 + OTF_BURN_REA_PTR_XS] = (double)loc2;
            }
          else
            WDB[loc1 + OTF_BURN_REA_PTR_XS] = NULLPTR;
          
          /* Check Pointer to yield data */

          if ((yld = (long)RDB[rea +  REACTION_PTR_FISSY]) > VALID_PTR)
            {
              /* Pointer to distribution */

              yld = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];
              CheckPointer(FUNCTION_NAME, "(yld)", DATA_ARRAY, yld);

              /* Loop over distribution */

              while (yld > VALID_PTR)
                {
                  /* Pointer to target */

                  tgt = (long)RDB[yld + FY_PTR_TGT];
                  CheckPointer(FUNCTION_NAME, "(tgt)", DATA_ARRAY, tgt);

                  /* Find target */
          
                  ptr = (long)RDB[mat + MATERIAL_PTR_OTF_BURN];
                  while (ptr > VALID_PTR)
                    {
                      /* Compare pointer */
                      
                      if (tgt == (long)RDB[ptr + OTF_BURN_PTR_NUCLIDE])
                        break;
                      
                      /* Next */
                      
                      ptr = NextItem(ptr);
                    }

                  /* Check pointer */
                  
                  if (ptr > VALID_PTR)
                    {
                      /*
                      printf("%s ", GetText(mat + MATERIAL_PTR_NAME));
                      for (n = 0; n < l; n++)
                        printf("  ");
                      printf("%5ld --> %s\n", (long)RDB[rea + REACTION_MT],
                             GetText(tgt + NUCLIDE_PTR_NAME));
                      */
                      /* Create structure */

                      loc2 = NewItem(loc1 + OTF_BURN_REA_PTR_FISSY,
                                     OTF_BURN_FISSY_BLOCK_SIZE);

                      /* Put pointers */

                      WDB[loc2 + OTF_BURN_FISSY_PTR_YLD] = (double)yld;
                      WDB[loc2 + OTF_BURN_FISSY_PTR_TGT] = (double)ptr;

                      /* Move forward in chain */
                      
                      LoopOTFBurnChain(ptr, l + 1);
                    }

                  /* Next */

                  yld = NextItem(yld);
                }
            }
        }

      /* Next Reaction */

      rea = NextItem(rea);
    }
}

/*****************************************************************************/
