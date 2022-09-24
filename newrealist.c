/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : newrealist.c                                   */
/*                                                                           */
/* Created:       2011/11/02 (JLe)                                           */
/* Last modified: 2017/08/11 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Creates and allocates memory for reaction lists              */
/*                                                                           */
/* Comments: - This routine is called only once, and it loops over all       */
/*             possible reactions in the list.                               */
/*                                                                           */
/*           - Completely changed 28.5.2012 (2.1.6)                          */
/*                                                                           */
/*           - RLS_DATA_PTR_NUCLIDE ei S(a,b) vaikutusaloille ole sama kuin  */
/*             RLS_DATA_PTR_REA --> REACTION_PTR_NUCLIDE                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NewReaList:"

/*****************************************************************************/

void NewReaList(long mat, long mode)
{
  long loc0, loc1, iso, nuc, rea, idx, ures, mt, ty, mul, ptr;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Check if reaction list already exists */

  if ((long)RDB[mat + mode] > VALID_PTR)
    Die(FUNCTION_NAME, "Reaction list already exists");

  /* Allocate memory for structure */

  loc0 = NewItem(mat + mode, RLS_BLOCK_SIZE);

  /* Put list pointer in material */

  WDB[mat + mode] = (double)loc0;

  /* Put material pointer and mode */

  WDB[loc0 + RLS_PTR_MAT] = (double)mat;
  WDB[loc0 + RLS_REA_MODE] = (double)mode;

  /* Allocate memory for next pointer */
 
  AllocValuePair(loc0 + RLS_PTR_NEXT);
 
  /* Check division */

  if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
    return;

  /* Reset composition index (index is updated at the beginning */
  /* if loop, because recycle can occur mid-way) */

  idx = -1;

  /* Loop over composition */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Update composition index */
      
      idx++;
      
      /* Pointer to nuclide data */
      
      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      
      /* Get ures-flag */
      
      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED)
        ures = YES;
      else
        ures = NO;
      
      /* Skip all but neutron transport and gamma nuclides */

      if (((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_TRANSPORT) &&
          ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON))
        {
          /* Next nuclide in composition */
          
          iso = NextItem(iso);
          
          /* Cycle loop */
          
          continue;
        }

      /* Reset reaction pointer */

      rea = -1;
      
      /* Check mode */

      if (mode == MATERIAL_PTR_TOT_REA_LIST)
        {
          /* Neutron total cross section */
          
          if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
            {
              /* Get pointer to reaction data */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            }
        }
      else if (mode == MATERIAL_PTR_TOT_URES_LIST)
        {
          /* Neutron total for ures list */
          
          if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) &&
              (ures == YES))
            {
              /* Get pointer to reaction data */
                
              rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            }
        }
      else if (mode == MATERIAL_PTR_ELA_URES_LIST)
        {
          /* Neutron elastic for ures list */
          
          if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) &&
              (ures == YES))
            {
              /* Get pointer to reaction data */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            }
        }
      else if (mode == MATERIAL_PTR_ABS_URES_LIST)
        {
          /* Neutron capture for ures list (NOTE: reaktiota isomeeriseen */
          /* tilaan ei tarvitse linkittää koska molemmissa on pointteri  */
          /* samaan factoriin) */

          if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) &&
              (ures == YES))
            {
              /* Get pointer to reaction data */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            }
        }
      else if (mode == MATERIAL_PTR_FISS_URES_LIST)
        {
          /* Neutron fission for ures list */
          
          if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) &&
              (ures == YES))
            {
              /* Get pointer to reaction data (no pointer check) */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];
            }
        }
      else if (mode == MATERIAL_PTR_HEAT_URES_LIST)
        {
          /* Neutron heating for ures list */
          
          if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) &&
              (ures == YES))
            {
              /* Get pointer to reaction data (no pointer check) */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS];
            }
        }
      else if (mode == MATERIAL_PTR_PHOT_TOT_LIST)
        {
          /* Photon total cross section (get pointer to photon nuclide */
          /* if type is different) */
          
          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON)   
            nuc = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_DATA];
          
          /* Check pointer */
          
          if (nuc > VALID_PTR)
            {
              /* Get pointer to reaction data */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            }
        }
      else if (mode == MATERIAL_PTR_PHOT_HEAT_LIST)
        {
          /* Photon heating cross section (get pointer to photon nuclide */
          /* if type is different) */
          
          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_PHOTON)   
            nuc = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_DATA];
          
          /* Check pointer */
          
          if (nuc > VALID_PTR)
            {
              /* Get pointer to reaction data */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_HEATPRODXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            }
        }
      else if (mode == MATERIAL_PTR_TMP_MAJORANT_LIST)
        {
          /* Neutron CE temperature majorant cross section */
          
          if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) &&
              ((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_CE))
            {
              /* Get pointer to reaction data */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
              
              /* Get pointer to majorant (tota ei välttämättä ole) */

              rea = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT];
            }
        }
    
      /* Check reaction pointer */

      if (rea > VALID_PTR)
        {
          /* Allocate memory for structure */

          loc1 = NewItem(loc0 + RLS_PTR_REA0, RLS_DATA_BLOCK_SIZE);
          
          /* Put pointer, composition index, and energy boundaries */
          
          WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
          WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
          WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
          WDB[loc1 + RLS_DATA_EMIN] = -INFTY;
          WDB[loc1 + RLS_DATA_EMAX] = INFTY;

          /* Allocate memory for reaction counter */
          
          ptr = AllocPrivateData(1, PRIVA_ARRAY);
          WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
        }
      else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
        {
          /* Other neutron reactions, Loop over reactions */
          
          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Get mt, ty and multiplication */
              
              mt = (long)RDB[rea + REACTION_MT];
              ty = (long)RDB[rea + REACTION_TY];
              mul = (long)RDB[rea + REACTION_WGT_F];

              /* Check reaction type */

              if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
                {
                  /* Absorption */
                  
                  if (mode == MATERIAL_PTR_ABS_REA_LIST)  
                    if (ty == 0)
                      {
                        /* Allocate memory for structure */

                        loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                       RLS_DATA_BLOCK_SIZE);
          
                        /* Put data */
          
                        WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                        WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                        WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                        WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                        WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

                        /* Allocate memory for reaction counter */
                        
                        ptr = AllocPrivateData(1, PRIVA_ARRAY);
                        WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                      }
                  
                  /* Elastic */
                  
                  if (mode == MATERIAL_PTR_ELA_REA_LIST)
                    if ((mt == 2) || (mt == 1002) || (mt == 1004)
                        || (mt == 2002) || (mt == 2004))
                      {
                        /* Allocate memory for structure */
                        
                        loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                       RLS_DATA_BLOCK_SIZE);
                        
                        /* Put data */
                        
                        WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                        WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                        WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                        WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                        WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

                        /* Allocate memory for reaction counter */
                        
                        ptr = AllocPrivateData(1, PRIVA_ARRAY);
                        WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                      }
                  
                  /* Fission */
                  
                  if (mode == MATERIAL_PTR_FISS_REA_LIST)
                    if (((mt > 17) && (mt < 22)) || (mt == 38))
                      {
                        /* Allocate memory for structure */
                        
                        loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                       RLS_DATA_BLOCK_SIZE);
                        
                        /* Put data */
                        
                        WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                        WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                        WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                        WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                        WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

                        /* Allocate memory for reaction counter */
                        
                        ptr = AllocPrivateData(1, PRIVA_ARRAY);
                        WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                      }
                  
                  /* Inelastic scattering production */
                  
                  if (mode == MATERIAL_PTR_INLP_REA_LIST)
                    if (mul > 1)
                      {
                        /* Allocate memory for structure */
                        
                        loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                       RLS_DATA_BLOCK_SIZE);
                        
                        /* Put data */
                        
                        WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                        WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                        WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                        WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                        WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

                        /* Allocate memory for reaction counter */
                        
                        ptr = AllocPrivateData(1, PRIVA_ARRAY);
                        WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                      }
                }
              else if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_SPECIAL)
                {
                  /* Heating */
                  
                  if (mode == MATERIAL_PTR_HEATT_REA_LIST)
                    if (mt == 301)
                      {
                        /* Allocate memory for structure */

                        loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                       RLS_DATA_BLOCK_SIZE);
                        
                        /* Put data */
                        
                        WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                        WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                        WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                        WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                        WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

                        /* Allocate memory for reaction counter */
                        
                        ptr = AllocPrivateData(1, PRIVA_ARRAY);
                        WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                      }

                  /* Photon production */
                  
                  if (mode == MATERIAL_PTR_PHOTP_REA_LIST)
                    if (mt == 202)
                      {
                        /* Allocate memory for structure */

                        loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                       RLS_DATA_BLOCK_SIZE);
                        
                        /* Put data */
                        
                        WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                        WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                        WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                        WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                        WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

                        /* Allocate memory for reaction counter */
                        
                        ptr = AllocPrivateData(1, PRIVA_ARRAY);
                        WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                      }

                  /* Proton production */

                  if ((long)RDB[DATA_INCLUDE_PROT_PROD_XS] == YES)
                    if (mode == MATERIAL_PTR_PROTP_REA_LIST)
                      if (mt == 203)
                        {
                          /* Allocate memory for structure */
                          
                          loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                         RLS_DATA_BLOCK_SIZE);
                          
                          /* Put data */
                          
                          WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                          WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                          WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                          WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                          WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];
                          
                          /* Allocate memory for reaction counter */
                          
                          ptr = AllocPrivateData(1, PRIVA_ARRAY);
                          WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                        }

                  /* Deuterium production */

                  if ((long)RDB[DATA_INCLUDE_DEUT_PROD_XS] == YES)
                    if (mode == MATERIAL_PTR_DEUTP_REA_LIST)
                      if (mt == 204)
                        {
                          /* Allocate memory for structure */
                          
                          loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                         RLS_DATA_BLOCK_SIZE);
                          
                          /* Put data */
                          
                          WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                          WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                          WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                          WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                          WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];
                          
                          /* Allocate memory for reaction counter */
                          
                          ptr = AllocPrivateData(1, PRIVA_ARRAY);
                          WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                        }

                  /* Tritium production */

                  if ((long)RDB[DATA_INCLUDE_TRIT_PROD_XS] == YES)
                    if (mode == MATERIAL_PTR_TRITP_REA_LIST)
                      if (mt == 205)
                        {
                          /* Allocate memory for structure */
                          
                          loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                         RLS_DATA_BLOCK_SIZE);
                          
                          /* Put data */
                          
                          WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                          WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                          WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                          WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                          WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];
                          
                          /* Allocate memory for reaction counter */
                          
                          ptr = AllocPrivateData(1, PRIVA_ARRAY);
                          WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                        }

                  /* He-3 production */

                  if ((long)RDB[DATA_INCLUDE_HE3_PROD_XS] == YES)
                    if (mode == MATERIAL_PTR_HE3P_REA_LIST)
                      if (mt == 206)
                        {
                          /* Allocate memory for structure */
                          
                          loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                         RLS_DATA_BLOCK_SIZE);
                          
                          /* Put data */
                          
                          WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                          WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                          WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                          WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                          WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];
                          
                          /* Allocate memory for reaction counter */
                          
                          ptr = AllocPrivateData(1, PRIVA_ARRAY);
                          WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                        }

                  /* He-4 production */

                  if ((long)RDB[DATA_INCLUDE_HE4_PROD_XS] == YES)
                    if (mode == MATERIAL_PTR_HE4P_REA_LIST)
                      if (mt == 207)
                        {
                          /* Allocate memory for structure */
                          
                          loc1 = NewItem(loc0 + RLS_PTR_REA0, 
                                         RLS_DATA_BLOCK_SIZE);
                          
                          /* Put data */
                          
                          WDB[loc1 + RLS_DATA_PTR_NUCLIDE] = (double)nuc;
                          WDB[loc1 + RLS_DATA_PTR_REA] = (double)rea;
                          WDB[loc1 + RLS_DATA_COMP_IDX] = (double)idx;
                          WDB[loc1 + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
                          WDB[loc1 + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];
                          
                          /* Allocate memory for reaction counter */
                          
                          ptr = AllocPrivateData(1, PRIVA_ARRAY);
                          WDB[loc1 + RLS_DATA_PTR_COUNT] = (double)ptr;
                        }
                } 
              
              /* Next reaction */
              
              rea = NextItem(rea);
            }
        }
      
      /* Next nuclide in composition */
      
      iso = NextItem(iso);
    }

  /* Check list pointer */

  if ((loc1 = (long)RDB[loc0 + RLS_PTR_REA0]) > VALID_PTR)
    {
      /* Close list */

      CloseList(loc1);
    }
  else
    {
      /* No reactions, set null pointer to erase list */

      WDB[mat + mode] = NULLPTR;
    }
}

/*****************************************************************************/
