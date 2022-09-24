/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : pulsedet.c                                     */
/*                                                                           */
/* Created:       2015/06/27 (JLe)                                           */
/* Last modified: 2015/08/01 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Analog pulse-height detector for photons                     */
/*                                                                           */
/* Comments: - Called with part = -1 after all histories in a batch are      */
/*             completed.                                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PulseDet:"

/*****************************************************************************/

void PulseDet(long part, long mat, double dE, double x0, double y0, double z0, 
              double wgt0, long id)
{
  long det, loc0, loc1, idx, ptr, pts, hix, i;
  double Etot, wgt;

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Check type */
      
      if (((long)RDB[det + DET_PARTICLE] != PARTICLE_TYPE_GAMMA) || 
          ((long)RDB[det + DET_PTR_SBINS] > VALID_PTR))
        {
          /* Next detector */
      
          det = NextItem(det);

          /* Cycle loop */

          continue;
        }

      /* Get pointer to response functions */

      loc0 = (long)RDB[det + DET_PTR_RBINS];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Check mt (no multiple reaction bins allowed) */

      if ((long)RDB[loc0 + DET_RBIN_MT] != MT_PHOTON_PULSE_HEIGHT)
        {
          /* Next detector */

          det = NextItem(det);

          /* Cycle loop */

          continue;
        }
          
      /* Pointer to pulse data */

      loc1 = (long)RDB[loc0 + DET_RBIN_PTR_PULSE_DATA];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
      
      /* Get pointer to statistics */
              
      pts = (long)RDB[det + DET_PTR_STAT];
      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

      /* Check call type */
      
      if (part == -1)
        {
          /*******************************************************************/

          /***** Called after last batch *************************************/

          /* Loop over OpenMP threads */
                  
          for (i = 0; i < (long)RDB[DATA_OMP_MAX_THREADS]; i++)
            {
              /* Retrieve stored energy and weight */
              
              Etot = GetPrivateData(loc1 + DET_PULSE_EDEP, i);
              wgt = GetPrivateData(loc1 + DET_PULSE_WGT, i);
              
              /* Pointer to energy grid */
              
              if ((ptr = (long)RDB[det + DET_PTR_EGRID]) < VALID_PTR)
                idx = 0;
              else
                {
                  /* Grid search */
                  
                  idx = GridSearch(ptr, Etot);
                }
              
              /* Check index */
              
              if (idx > -1)
                {
                  /* Score pulse */
                  
                  AddBuf(1.0, wgt, pts, i, -1, idx, 0);
                }
              
              /* Reset pulse data */
                  
              PutPrivateData(loc1 + DET_PULSE_PHOTON_IDX, -1, i);
              PutPrivateData(loc1 + DET_PULSE_EDEP, 0.0, i);
              PutPrivateData(loc1 + DET_PULSE_WGT, -1.0, i);
            }
          
          /***************************************************************/
        }
      else
        {
          /*******************************************************************/

          /***** Called from collision point *********************************/

          /* Check particle pointer */

          CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
          
          /* Check bin (this is just to check that the collision is */
          /* inside the detector, the index is not relevant here).  */
          
          if (DetBin(det, mat, part, x0, y0, z0, -INFTY, 0.0, id) < 0)
            {
              /* Next detector */

              det = NextItem(det);

              /* Cycle loop */
              
              continue;
            }
          
          /* Get history index */

          hix = (long)RDB[part + PARTICLE_HISTORY_IDX] + 1;
          
          /* Check index */
          
          if (hix == (long)GetPrivateData(loc1 + DET_PULSE_PHOTON_IDX, id))
            {
              /* Still same history, add to temporary storage */
              
              AddPrivateData(loc1 + DET_PULSE_EDEP, dE, id);

              /* Compare weight */

              if (wgt0 != GetPrivateData(loc1 + DET_PULSE_WGT, id))
                Die(FUNCTION_NAME, "Mismatch in weight");
            }
          else
            {
              /* New history, retrieve stored energy and weight */
              
              Etot = GetPrivateData(loc1 + DET_PULSE_EDEP, id);
              wgt = GetPrivateData(loc1 + DET_PULSE_WGT, id);

              /* Pointer to energy grid */
              
              if ((ptr = (long)RDB[det + DET_PTR_EGRID]) < VALID_PTR)
                idx = 0;
              else
                {
                  /* Grid search */
                  
                  idx = GridSearch(ptr, Etot);
                }
              
              /* Check index */
              
              if (idx > -1)
                {
                  /* Score pulse */
                  
                  AddBuf(1.0, wgt, pts, id, -1, idx, 0);
                }
              
              /* Store new data */

              PutPrivateData(loc1 + DET_PULSE_PHOTON_IDX, hix, id);
              PutPrivateData(loc1 + DET_PULSE_EDEP, dE, id);
              PutPrivateData(loc1 + DET_PULSE_WGT, wgt0, id);
            }

          /*******************************************************************/
        }
      
      /* Next detector */
      
      det = NextItem(det);
    }
}

/*****************************************************************************/
