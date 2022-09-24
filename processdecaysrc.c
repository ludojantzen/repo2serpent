/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdecaysrc.c                              */
/*                                                                           */
/* Created:       2015/05/23 (JLe)                                           */
/* Last modified: 2018/10/29 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes decay source into a sampling list                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessDecaySrc:"

/*****************************************************************************/

void ProcessDecaySrc()
{  
  long mat, iso, nuc, ptr, src, rad, type, idx;
  double vol, adens, lambda, totp, totn, prev, I;

  fprintf(outp, "Processing decay source...\n");
  
  /* Check total source rate */
  
  if ((RDB[DATA_TOT_PHOTON_DEC_SRC_RATE] == 0.0) && 
      (RDB[DATA_TOT_NEUTRON_DEC_SRC_RATE] == 0.0))
    Error(0, "Total decay source rate is zero in decay source mode");

  /* Reset nuclide index */

  idx = -1;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get volume */

      vol = RDB[mat + MATERIAL_VOLUME];

      /* Check total photon source rate (if material has neutron emission */
      /* reactions, it also should have photon emission reactions). */

      if ((vol == 0.0) || ((RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE]/
                            RDB[DATA_TOT_PHOTON_DEC_SRC_RATE]) < 1E-19))
        {
          /* Next material */
          
          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }
      
      /* Check that pointers are not defined (NOTE: Tää on lähinnä sitä */
      /* varten että jos tätä rutiinia joskus kutsutaan silmukasta niin */
      /* muistia ei varata turhaan. */

      if ((long)RDB[mat + MATERIAL_PTR_PHOTON_DECAY_SRC] > VALID_PTR)
        Die(FUNCTION_NAME, "Pointer to decay source already exists");

      if ((long)RDB[mat + MATERIAL_PTR_NEUTRON_DECAY_SRC] > VALID_PTR)
        Die(FUNCTION_NAME, "Pointer to decay source already exists");

      /* Avoid compiler warning */

      src = -1;

      /* Reset totals */

      totp = 0.0;
      totn = 0.0;

      /* Pointer to composition (in photon transport mode the original */
      /* nuclide composition is in another list) */

      if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) < VALID_PTR)
        iso = (long)RDB[mat + MATERIAL_PTR_COMP];      

      /* Loop over composition */

      while (iso > VALID_PTR)
        {
          /* Get atomic density */

          adens = RDB[iso + COMPOSITION_ADENS]*1E+24;

          /* Get pointer to nuclide data */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Get decay constant */

          lambda = RDB[nuc + NUCLIDE_LAMBDA];

          /* Check radiations */

          if ((long)RDB[nuc + NUCLIDE_PTR_RADIATIONS] > VALID_PTR)
            idx++;

          /* Loop over radiations */

          rad = (long)RDB[nuc + NUCLIDE_PTR_RADIATIONS];
          while (rad > VALID_PTR)
            {
              /* Get type */

              type = (long)RDB[rad + NUCLIDE_RAD_TYPE];

              /* Get total intensity */

              I = RDB[rad + NUCLIDE_RAD_SPEC_I];

              /* Check type */

              if (type == PARTICLE_TYPE_GAMMA)
                {
                  /***********************************************************/

                  /***** Photon source ***************************************/
          
                  /* Check intensity */

                  if (I*lambda*adens*vol/
                      RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE] > 1E-18)
                    {
                      /* Create entry */
                      
                      src = NewItem(mat + MATERIAL_PTR_PHOTON_DECAY_SRC, 
                                    SRC_DECCAY_BLOCK_SIZE);

                      /* Put index */

                      WDB[src + SRC_DECCAY_IDX] = (double)idx;

                      /* Put nuclide pointer */
                      
                      WDB[src + SRC_DECCAY_PTR_NUCLIDE] = (double)nuc;
                      
                      /* Put emission rate */
                      
                      WDB[src + SRC_DECCAY_I] = I*lambda*adens*vol;
              
                      /* Reset weight */
                      
                      WDB[src + SRC_DECCAY_WGT] = 1.0;
                      
                      /* Add to total */
                      
                      totp = totp + RDB[src + SRC_DECCAY_I];
                      
                      /* Put pointer to emission spectra */
              
                      ptr = (long)RDB[rad + NUCLIDE_RAD_PTR_SPEC];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      WDB[src + SRC_DECCAY_PTR_SPEC] = (double)ptr;
                    }
                  
                  /***********************************************************/
                }
              else if (type == PARTICLE_TYPE_NEUTRON)
                {
                  /***********************************************************/
                  
                  /***** Neutron source **************************************/
                  
                  /* Check intensity */
                  
                  if (I*lambda*adens*vol/
                      RDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE] > 1E-18)
                    {
                      /* Create entry */
                      
                      src = NewItem(mat + MATERIAL_PTR_NEUTRON_DECAY_SRC, 
                                    SRC_DECCAY_BLOCK_SIZE);
                      
                      /* Put nuclide pointer */
                      
                      WDB[src + SRC_DECCAY_PTR_NUCLIDE] = (double)nuc;
                      
                      /* Put emission rate */
                      
                      WDB[src + SRC_DECCAY_I] = I*lambda*adens*vol;
                      
                      /* Reset weight */
                      
                      WDB[src + SRC_DECCAY_WGT] = 1.0;
                      
                      /* Add to total */
                      
                      totn = totn + RDB[src + SRC_DECCAY_I];
                      
                      /* Put pointer to emission spectra */
                      
                      ptr = (long)RDB[rad + NUCLIDE_RAD_PTR_SPEC];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      WDB[src + SRC_DECCAY_PTR_SPEC] = (double)ptr;
                    }

                  /***********************************************************/
                }
              else if ((type != PARTICLE_TYPE_ELECTRON) &&
                       (type != PARTICLE_TYPE_POSITRON) &&
                       (type != PARTICLE_TYPE_ALPHA))
                Die(FUNCTION_NAME, "Invalid radiation type");

              /* Next radiation */

              rad = NextItem(rad);
            }
          
          /* Next nuclide */
          
          iso = NextItem(iso);
        }

      /***********************************************************************/
     
      /***** Check if photon source was created ******************************/

      src = (long)RDB[mat + MATERIAL_PTR_PHOTON_DECAY_SRC];
      if (src > VALID_PTR)
        {
          /* Close list */

          CloseList(src);

          /* Sort */

          SortList(src, SRC_DECCAY_I, SORT_MODE_DESCEND);

          /* Check total */

          if (totp == 0.0)
            Die(FUNCTION_NAME, "WTF?");
          
          /* Calculate cumulative probabilities */

          prev = 0.0;
          
          src = (long)RDB[mat + MATERIAL_PTR_PHOTON_DECAY_SRC];
          while (src > VALID_PTR)
            {
              /* Calculate cumulative probability */
              
              WDB[src + SRC_DECCAY_CUM_P] = 
                (prev + RDB[src + SRC_DECCAY_I])/totp;
              
              /* Update previous */
              
              prev = prev + RDB[src + SRC_DECCAY_I];

              /* Next */
              
              src = NextItem(src);
            }
        }

      /***********************************************************************/
      
      /***** Check if neutron source was created *****************************/
      
      src = (long)RDB[mat + MATERIAL_PTR_NEUTRON_DECAY_SRC];
      if (src > VALID_PTR)
        {
          /* Close list */

          CloseList(src);

          /* Sort */

          SortList(src, SRC_DECCAY_I, SORT_MODE_DESCEND);

          /* Check total */

          if (totn == 0.0)
            Die(FUNCTION_NAME, "WTF?");
          
          /* Calculate cumulative probabilities */

          prev = 0.0;
          
          src = (long)RDB[mat + MATERIAL_PTR_NEUTRON_DECAY_SRC];
          while (src > VALID_PTR)
            {
              /* Calculate cumulative probability */
              
              WDB[src + SRC_DECCAY_CUM_P] = 
                (prev + RDB[src + SRC_DECCAY_I])/totn;
              
              /* Update previous */
              
              prev = prev + RDB[src + SRC_DECCAY_I];
              
              /* Next */
              
              src = NextItem(src);
            }
        }

      /***********************************************************************/
      
      /* Allocate memory for sampled stats */
      
      if ((totp > 0.0) || (totn > 0.0))
        {
          ptr = NewStat("SAMPLED_DECAY_SRC", 1, 2);
          WDB[mat + MATERIAL_SAMPLED_DECAY_SRC] = (double)ptr;  
        }

      /* Next material */

      mat = NextItem(mat);
    } 

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
