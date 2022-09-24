/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoretransmuxs.c                               */
/*                                                                           */
/* Created:       2018/04/02 (JLe)                                           */
/* Last modified: 2018/05/15 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores transmutation cross sections and fission product      */
/*              yield weights for on-the-fly burnup calculation              */
/*                                                                           */
/* Comments: - Tää on erotettu omaksi aliohjelmaksi jotta ei rikota          */
/*             scoretransmuxs.c:tä.                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreOTFBurn:"

/*****************************************************************************/

void ScoreOTFBurn(double flx, long mat, double E, double wgt, long id)
{
  long rea0, rea, loc0, loc1, ptr, RFS;
  double val, E0, E1, E2, f, T, Er;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return;

  /* Pointer to OTF burn data */

  if ((loc0 = (long)RDB[mat + MATERIAL_PTR_OTF_BURN]) < VALID_PTR)
    return;

  /* Check burnup mode */

 if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
   return;

 /* Check if in coefficient calculation */

 if ((long)RDB[DATA_COEF_CALC_IDX] > 0)
   return;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return;

  /* Check material burn flag */

  if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
    return;  

  /* Check predictor-corrector calculation */

  if (((long)RDB[DATA_BURN_CORR_TYPE] != CORR_TYPE_NONE) &&
      ((long)RDB[DATA_BURN_STEP_PC] != CORRECTOR_STEP))
    return;

  /* Check parent type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
    Die(FUNCTION_NAME, "Parent type");

  /* Check MPI id when domain decomposition is in use */

  if (((long)RDB[DATA_DD_DECOMPOSE] == YES) && 
      ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
    Die(FUNCTION_NAME, "Error in MPI id");

  /* Check flux */

  CheckValue(FUNCTION_NAME, "flx", "", flx, ZERO, INFTY);

  /***************************************************************************/

  /***** Score transmutation rates *******************************************/

  /* Loop over nuclides */

  while (loc0 > VALID_PTR)
    {
      /* Check material pointer */

      if ((long)RDB[loc0 + OTF_BURN_PTR_MAT] != mat)
        break;

      /* Loop over transmutation reactions */

      loc1 = (long)RDB[loc0 + OTF_BURN_PTR_REA];
      while (loc1 > VALID_PTR)
        {
          /* Pointer to reaction data */
      
          rea0 = (long)RDB[loc1 + OTF_BURN_REA_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

          /* Get pointer to parent data if branching to isomeric state */

          if ((rea = (long)RDB[rea0 + REACTION_PTR_BRANCH_PARENT]) < VALID_PTR)
            rea = rea0;

          /* Check energy */

          if (E < RDB[rea + REACTION_EMIN])
            {
              /* Pointer to next */

              loc1 = NextItem(loc1);

              /* Cycle loop */

              continue;
            }

          /* Check mt */

          if (((long)RDB[rea + REACTION_MT] < 16) && 
              ((long)RDB[rea + REACTION_MT] != 4))
            Die(FUNCTION_NAME, "Error in mt");

          /* Check type */
          
          if ((long)RDB[rea + REACTION_MT] > 10000)
            break;

          /* Get cross section */
          
          if ((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_NONE)
            val = MicroXS(rea, E, id);
          else if ((T = GetTemp(mat, id)) < ZERO)
            val = MicroXS(rea, E, id);
          else            
            val = DopMicroXS(mat, rea, E, &Er, T, id);
          
          /* Check for fission */

          if ((long)RDB[loc1 + OTF_BURN_REA_PTR_FISSY] < VALID_PTR)
            {
              /***************************************************************/

              /***** Non-fission transmutation *******************************/

              /* Check if reaction is associated with isomeric branching */
          
              if ((long)RDB[rea + REACTION_PTR_ISO_BRA] > VALID_PTR)
                {
                  /* Get fraction */
                  
                  RFS = (long)RDB[rea0 + REACTION_RFS];
                  f = BranchFrac(rea0, RFS, E, id);
                }
              else
                f = 1.0;

              /* Score total */
              
              ptr = (long)RDB[loc0 + OTF_BURN_PTR_TOT_XS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf1D(f*val*flx, wgt, ptr, id, 0);

              /* Score reaction rate */

              if ((ptr = (long)RDB[loc1 + OTF_BURN_REA_PTR_XS]) > VALID_PTR)
                AddBuf1D(f*val*flx, wgt, ptr, id, 0);

              /***************************************************************/
            }
          else
            {
              /***************************************************************/

              /***** Fission *************************************************/

              /* Get interpolation energies */

              E0 = RDB[rea0 + REACTION_FISSY_IE0];
              E1 = RDB[rea0 + REACTION_FISSY_IE1];
              E2 = RDB[rea0 + REACTION_FISSY_IE2];

              /* Check energy */

              if (E < E0)
                {
                  /* Pointer to next */
                  
                  loc1 = NextItem(loc1);
                  
                  /* Cycle loop */
                  
                  continue;
                }

              /* Check values */
              
              CheckValue(FUNCTION_NAME, "E1" ,"", E1, E0, E2);
              CheckValue(FUNCTION_NAME, "E" ,"", E, E0, INFTY);
              
              /* Check interval */
              
              if (E < E2)
                {
                  /* Interpolation factor */
                  
                  if (E < E1)
                    {
                      /* Below interpolation energy */
                      
                      if (E0 < 0.0)
                        f = 1.0;
                      else
                        f = (E - E0)/(E1 - E0);
                    }
                  else
                    {
                      /* Above interpolation energy */
                      
                      if (E2 > 1000.0)
                        f = 1.0;
                      else
                        f = (E2 - E)/(E2 - E1);
                    }
                  
                  /* Check factor */
                  
                  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);

                  /* Score total */
                  
                  ptr = (long)RDB[loc0 + OTF_BURN_PTR_TOT_XS];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  AddBuf1D(f*val*flx, wgt, ptr, id, 0);
                  
                  /* Score reaction rate */
              
                  if ((ptr = (long)RDB[loc1 + OTF_BURN_REA_PTR_XS]) > 
                      VALID_PTR)
                    AddBuf1D(f*val*flx, wgt, ptr, id, 0);
                }
            }
          
          /* Next reaction */

          loc1 = NextItem(loc1);
        }
      
      /* Next nuclide */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/
}

/*****************************************************************************/
