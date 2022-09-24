/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorescattering.c                              */
/*                                                                           */
/* Created:       2011/03/03 (JLe)                                           */
/* Last modified: 2018/05/24 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores parameters related to scattering reactions            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreScattering:"

/*****************************************************************************/

void ScoreScattering(long mat, long rea, double mu, double E0, double E1, 
                     double wgt1, double wgt2, long id)
{
  long gcu, ptr, loc0, uni, ntot, ng0, ng1, mt, ncol, nmu, iso, nuc;
  double P0, P1, P2, P3, P4, P5, P6, P7;

  /* Check weight (mt 5 coi asettaa wgt2 < wgt1) */

  CheckValue(FUNCTION_NAME, "wgt1", "", wgt1, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "wgt2", "", wgt2, ZERO, INFTY);

  /* Check mu (toi mu lasketaan suuntavektoreiden skalaaritulona, ja se   */
  /* voi numeriikasta johtuen poiketa ykkösestä, vaikka se tarkastetaankin */
  /* kaikissa sämpläysrutiineissa. */

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.000001, 1.000001);
  
  /***************************************************************************/

  /***** Analog reaction rates ***********************************************/

  /* Check modes */

  if ((long)RDB[DATA_OPTI_IMPLICIT_RR] == NO)
    {
      /* Total reaction rate */

      ptr = (long)RDB[RES_TOT_NEUTRON_RR];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, wgt1, ptr, id, 0);

      /* Get mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Score elastic or inelastic scattering */

      if ((mt == 2) || (mt == 1002) || (mt == 1004))
        {
          ptr = (long)RDB[RES_TOT_ELARATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, wgt1, ptr, id, 0);
        }
      else if (wgt2 > wgt1)
        {
          ptr = (long)RDB[RES_TOT_INLPRODRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(wgt2/wgt1 - 1.0, wgt1, ptr, id, 0);
        }

      /* Exit subroutine */

      return;
    }

  /***************************************************************************/

  /***** Check if data is required *******************************************/

  /* Check active cycle and corrector step */

  if ((RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]) ||
      (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
       ((long)RDB[DATA_B1_BURNUP_CORR] == NO)))
    return;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Check for multiple levels */

  if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
    {
      /* Single level, get pointer */

      if ((gcu = (long)TestValuePair(DATA_GCU_PTR_UNI, (double)ncol, id)) 
          < VALID_PTR)
        return;
    }
  else
    {
      /* Multiple levels, get pointer to list */

      gcu = (long)RDB[DATA_PTR_GCU0];
      CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);
    }

  /* Loop over universes */
  
  while (gcu > VALID_PTR)
    {
      /* Check multi-level mode */
          
      if ((long)RDB[DATA_MULTI_LEVEL_GCU] == YES)
        {
          /* Pointer to universe */
          
          uni = (long)RDB[gcu + GCU_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
          
          /* Check collision */
          
          if (TestValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, id) < 0.0)
            {
              /* Next universe */
              
              gcu = NextItem(gcu);
              
              /* Cycle loop */
              
              continue;
            }
        }
      
      /* Legendre polynomial coefficients */
      
      P0 = 1.0;
      P1 = mu;
      P2 = 0.5*(3.0*mu*mu - 1.0);
      P3 = 0.5*(5.0*mu*mu*mu - 3.0*mu);
      P4 = (35.0*mu*mu*mu*mu - 30.0*mu*mu + 3.0)/8.0;
      P5 = (63.0*mu*mu*mu*mu*mu - 70.0*mu*mu*mu + 15*mu)/8.0;
      P6 = (231.0*mu*mu*mu*mu*mu*mu - 315.0*mu*mu*mu*mu + 105*mu*mu 
            - 5.0)/16.0;
      P7 = (429.0*mu*mu*mu*mu*mu*mu*mu - 693.0*mu*mu*mu*mu*mu + 315*mu*mu*mu 
            - 35.0*mu)/16.0;
      
      /***********************************************************************/
      
      /***** MORA data *******************************************************/
      
      /* Check pointer */
      
      if ((loc0 = (long)RDB[gcu + GCU_PTR_MORA]) > VALID_PTR)
        {
          /* Get pointer to energy grid */
          
          ptr = (long)RDB[loc0 + MORA_PTR_EG];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Get group indexes (E0 --> E1) */
          
          ng0 = GridSearch(ptr, E0);
          ng1 = GridSearch(ptr, E1);
          
          /* Get mu bin */

          if (mu == 1.0)
            nmu = (long)RDB[loc0 + MORA_N_COS] - 1;
          else
            nmu = (long)(0.5*(mu + 1.0)*RDB[loc0 + MORA_N_COS]);
          
          /* Score scattering rate */
          
          if ((ng0 > -1) && (ng1 > -1))
            {
              ptr = (long)RDB[loc0 + MORA_PTR_SCATTP];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(1.0, wgt1, ptr, id, -1, ng1, ng0, nmu);
              
              ptr = (long)RDB[loc0 + MORA_PTR_SCATTW];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              AddBuf(wgt2/wgt1, wgt1, ptr, id, -1, ng1, ng0, nmu);
            }
        }

      /***********************************************************************/

      /***** Micro-group data ************************************************/
      
      /* Get pointer to microgroup energy grid */
      
      ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);        
      
      /* Number of groups */
      
      ntot = (long)RDB[ptr + ENERGY_GRID_NE] - 1;
      
      /* Get micro-group indexes (E0 --> E1) */
      
      ng0 = GridSearch(ptr, E0);
      ng1 = GridSearch(ptr, E1);
      
      /* Check incident group */
      
      if ((ng0 > -1) && (ng1 > -1))
        {
          /* Invert */

          ng0 = ntot - ng0 - 1;
          CheckValue(FUNCTION_NAME, "ng0", "", ng0, 0, ntot - 1);
          
          ng1 = ntot - ng1 - 1;
          CheckValue(FUNCTION_NAME, "ng1", "", ng1, 0, ntot - 1);
          
          /* Scattering matrixes */
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATT0];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt1*P0, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP0];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt2*P0, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATT1];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt1*P1, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP1];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt2*P1, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATT2];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt1*P2, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP2];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt2*P2, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATT3];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt1*P3, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP3];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt2*P3, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATT4];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt1*P4, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP4];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt2*P4, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATT5];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt1*P5, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP5];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt2*P5, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATT6];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt1*P6, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP6];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt2*P6, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATT7];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt1*P7, id);
          
          ptr = (long)RDB[gcu + GCU_MICRO_SCATTP7];
          CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
          AddPrivateRes(ptr + ng0*ntot + ng1, wgt2*P7, id);

          /* Transport correction (NOTE: Täällä tosta vaikutusalasta */
          /* vähennetään P1-komponentti niiltä nuklideilta joilla    */
          /* korjausta ei tehdä. Kokonaisvaikutusala lisätään sitten */
          /* scoregc.c:ssä. */

          if ((long)RDB[DATA_PTR_TRC0] > VALID_PTR)
            {
              /* Check if material is associated with transport correction */

              if ((ptr = (long)RDB[mat + MATERIAL_PTR_TRANSP_CORR])
                  < VALID_PTR)
                {
                  /* No correction for material, score */

                    ptr = (long)RDB[gcu + GCU_MICRO_TRC];
                    CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                    AddPrivateRes(ptr + ng0, -wgt1*P1, id);
                }
              else if (E0 < RDB[ptr + TRANSP_CORR_EMIN])
                {
                  /* Below limit, score */

                    ptr = (long)RDB[gcu + GCU_MICRO_TRC];
                    CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                    AddPrivateRes(ptr + ng0, -wgt1*P1, id);
                }
              else if ((loc0 = (long)RDB[ptr + TRANSP_CORR_PTR_ISO]) 
                       > VALID_PTR)
                {
                  /* Get nuclide */

                  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                  /* Loop over list */
                  
                  while ((iso = (long)RDB[loc0++]) > VALID_PTR)
                    {
                      /* Pointer to nuclide data */
                      
                      ptr = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                      CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);
                      
                      /* Compare and break loop */
                      
                      if ((long)RDB[nuc + NUCLIDE_ZAI] != 
                          (long)RDB[ptr + NUCLIDE_ZAI])
                        break;
                    }

                  /* Check if found */

                  if (iso < VALID_PTR)                  
                    {
                      /* No correction for nuclide, score */
                      
                      ptr = (long)RDB[gcu + GCU_MICRO_TRC];
                      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
                      AddPrivateRes(ptr + ng0, -wgt1*P1, id);
                    }
                }
            }
        }
      
      /***********************************************************************/

      /* Next universe */

      if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
        break;
      else
        gcu = NextItem(gcu);
    }
}

/*****************************************************************************/
