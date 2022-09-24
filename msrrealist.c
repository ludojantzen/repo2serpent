/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : msrrealist.c                                   */
/*                                                                           */
/* Created:       2015/03/28 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Stores transmutation reactions into a separate structure     */
/*                                                                           */
/* Comments: - Needed because the original burnup routine is not compatible  */
/*             with MSR continuous-reprocessing                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MSRReaList:"

/*****************************************************************************/

double **MSRReaList(double **R, long mat, long id)
{
  long mat0, mat1, nm, sz, iso, nuc, rea, n, i, type;

  /* Get pointer to first material */
  
  mat0 = (long)RDB[mat + MATERIAL_FLOW_PTR_FIRST];
  CheckPointer(FUNCTION_NAME, "(mat0)", DATA_ARRAY, mat0);
  
  /* Get number of materials in chain */
  
  nm = (long)RDB[mat0 + MATERIAL_FLOW_N];
  CheckValue(FUNCTION_NAME, "nm", "", nm, 1, 20);

  /***************************************************************************/

  /***** Allocate memory for structure ***************************************/

  /* Check pointer */

  if (R == NULL)
    {
      /* Reset number of reactions */

      sz = 0;

      /* Loop over composition */
      
      iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */
          
          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              
              /* Check pointer  (pointteri on PRIVA-blokkiin) */
              
              if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
                sz++;
              
              /* Next reaction */
              
              rea = NextItem(rea);
            }

          /* Check for total fission */

          if ((long)RDB[nuc + NUCLIDE_PTR_NFY_DATA] < VALID_PTR)
            if ((rea = (long)RDB[nuc + NUCLIDE_PTR_TOTFISS_REA]) > VALID_PTR)
              if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
                sz++;
          
          /* Next */
          
          iso = NextItem(iso);
        }

      /* Allocate memory */
          
      R = (double **)Mem(MEM_ALLOC, 4, sizeof(double *));        
          
      for(n = 0; n < 4; n++)
        R[n] = (double *)Mem(MEM_ALLOC, nm*sz, sizeof(double));
    }

  /***************************************************************************/
  
  /***** Read data ***********************************************************/

  /* Reset count */
  
  n = 0;

  /* Loop over composition */
  
  iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Pointer to nuclide */
      
      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
      
      /* Loop over materials */

      for (i = 0; i < nm; i++)
        {
          /* Find corresponding material */

          mat1 = (long)RDB[DATA_PTR_M0];
          while (mat1 > VALID_PTR)
            {
              /* Check index and that material belongs to same chain */
              
              if ((long)RDB[mat1 + MATERIAL_FLOW_PTR_FIRST] == mat0)
                if ((i % nm) == (long)RDB[mat1 + MATERIAL_FLOW_IDX] - 1)
                  break;
                 
              /* Pointer to next */

              mat1 = NextItem(mat1);
            }

          /* Check pointer */
          
          CheckPointer(FUNCTION_NAME, "(mat0)", DATA_ARRAY, mat0);
          
          /* Loop over reactions */
      
          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {          
              /* Get type */

              type = (long)RDB[rea + REACTION_TYPE];
              
              /* Check fission yield and reaction type */

              if ((long)RDB[rea + REACTION_PTR_FISSY] > VALID_PTR)
                if (type != REACTION_TYPE_DECAY)
                  {
                    /* Check material pointer */
                    
                    if (mat1 == mat)
                      {
                        /* Store material, nuclide and reaction pointers */
                        
                        R[0][n] = (double)mat;
                        R[1][n] = (double)nuc;
                        R[2][n] = (double)rea;
                        
                        /* Store rate */
                        
                        R[3][n] = TestValuePair(rea + REACTION_PTR_TRANSMUXS, 
                                                (double)mat, id);
                      }
                    
                    /* Update index */
                    
                    n++;
                  }

              /* Check target, type and pointer */
              /* (pointteri on PRIVA-blokkiin)  */

              if ((long)RDB[rea + REACTION_PTR_TGT] > VALID_PTR)
                if ((type == REACTION_TYPE_PARTIAL) || 
                    (type == REACTION_TYPE_TRA_BRANCH))
                  if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
                    {
                      /* Check material pointer */
                      
                      if (mat1 == mat)
                        {
                          /* Store material, nuclide and reaction pointers */
                          
                          R[0][n] = (double)mat;
                          R[1][n] = (double)nuc;
                          R[2][n] = (double)rea;
                          
                          /* Store rate */
                          
                          R[3][n] = TestValuePair(rea + REACTION_PTR_TRANSMUXS, 
                                                  (double)mat, id);
                        }
                      
                      /* Update index */
                      
                      n++;
                    }
              
              /* Next reaction */
          
              rea = NextItem(rea);
            }
      
          /* Check for total fission */

          if ((long)RDB[nuc + NUCLIDE_PTR_NFY_DATA] < VALID_PTR)
            if ((rea = (long)RDB[nuc + NUCLIDE_PTR_TOTFISS_REA]) > VALID_PTR)
              if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
                {
                  /* Check material pointer */

                  if (mat1 == mat)
                    {
                      /* Store material, nuclide and reaction pointers */

                      R[0][n] = (double)mat;
                      R[1][n] = (double)nuc;
                      R[2][n] = (double)rea;

                      /* Store rate */

                      R[3][n] = TestValuePair(rea + REACTION_PTR_TRANSMUXS, 
                                              (double)mat, id);
                    }
                  
                  /* Update index */
                  
                  n++;
                }
        }
      
      /* Next */
      
      iso = NextItem(iso);
    }

  /***************************************************************************/
  
  /***** Reset data **********************************************************/

  /* Reset count */
  
  n = 0;

  /* Loop over composition */
  
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
          /* Check pointer  (pointteri on PRIVA-blokkiin) */
          
          if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
            StoreValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, 0.0, id);
          
          /* Next reaction */
          
          rea = NextItem(rea);
        }
      
      /* Check for total fission */
      
      if ((long)RDB[nuc + NUCLIDE_PTR_NFY_DATA] < VALID_PTR)
        if ((rea = (long)RDB[nuc + NUCLIDE_PTR_TOTFISS_REA]) > VALID_PTR)
          if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
            StoreValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, 0.0, id);
      
      /* Next */
      
      iso = NextItem(iso);
    }

  /* Return table */

  return R;
  
  /***************************************************************************/
}

/*****************************************************************************/
