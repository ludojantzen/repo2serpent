/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdephis.c                                */
/*                                                                           */
/* Created:       2011/06/14 (JLe)                                           */
/* Last modified: 2017/03/14 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Links normalization, sets  pcc modes, etc. for depletion     */
/*              histories                                                    */
/*                                                                           */
/* Comments: TODO: jos normalisaatio on asetettu niin, että neutronivuo on   */
/*                 nolla, askel pitäisi tulkita hajoamisaskeleeksi           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessDepHis:"

/*****************************************************************************/

void ProcessDepHis()
{
  long dep, norm, n, ptr, ntot, bra;

  /* Check burnup mode */

 if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
   return;

 /* Check pointer to depletion histories */

 if ((long)RDB[DATA_BURN_PTR_DEP] < VALID_PTR)
   return;

  /* Reset counter */

  n = 0;

  /* Loop over depletion intervals */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Add to number of burnup intervals */
      
      if (((long)RDB[dep + DEP_HIS_STEP_TYPE] != DEP_STEP_DEC_TOT) &&
          ((long)RDB[dep + DEP_HIS_STEP_TYPE] != DEP_STEP_DEC_STEP) &&
          ((long)RDB[dep + DEP_HIS_STEP_TYPE] != DEP_STEP_ACT_TOT) &&
          ((long)RDB[dep + DEP_HIS_STEP_TYPE] != DEP_STEP_ACT_STEP))
        n++;
      
      /* Check activation intervals */

      if (((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_ACT_TOT) ||
          ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_ACT_STEP))
        {
          /* Must be preceded by at least one burn interval */

          if (n == 0)
            Error(dep, 
                  "Activation interval must be preceded by burn interval");

          /* Check previous type */

          ptr = PrevItem(dep);
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          if (((long)RDB[ptr + DEP_HIS_STEP_TYPE] == DEP_STEP_DEC_TOT) ||
              ((long)RDB[ptr + DEP_HIS_STEP_TYPE] == DEP_STEP_DEC_STEP))
            Error(dep, "Decay interval cannot precede activation interval");
      
          /* Warn about burnup */
          
          Note(dep, 
               "Burnup is not calculated correctly in activation intervals");
        }

      /* Remember previous */

      ptr = dep;

      /* Next */
      
      dep = NextItem(dep);
    }

  /* Avoid compiler warning */

  norm = -1;

  /* Set decay only mode and check normalization (toi restart file */
  /* testataan että saataisiin tehtyä restartit hajoamislaskuina)  */

  if ((n == 0) && ((long)RDB[DATA_READ_RESTART_FILE] == NO))
    WDB[DATA_BURN_DECAY_CALC] = (double)YES;
  else if (((norm = (long)RDB[DATA_PTR_NORM]) < VALID_PTR) &&
           ((long)RDB[DATA_USE_DECAY_SRC] == NO))
    Error(0, "Normalization must be defined in burnup mode");

  /* Tää on viritys jolla saadaan matlaboutput.c tulostamaan jos */
  /* restartataan decay-laskuun */

  if ((n == 0) && ((long)RDB[DATA_READ_RESTART_FILE] == YES))
    WDB[DATA_CYCLE_IDX] = RDB[DATA_CRIT_SKIP] + 1.0;

  /* Set pointers to normalization */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Set pointer to normalization */

      WDB[dep + DEP_HIS_PTR_NORM] = (double)norm;

      /* Next normalization */

      if (norm > VALID_PTR)
        if (NextItem(norm) > VALID_PTR)
          norm = NextItem(norm);
      
      /* Next */

      dep = NextItem(dep);
    }

  /* Check steps */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Get pointer to steps */

      ptr = (long)RDB[dep + DEP_HIS_PTR_STEPS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get total number of steps */
      
      ntot = (long)RDB[dep + DEP_HIS_N_STEPS];

      /* Check step type */

      if (((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_BU_TOT) ||
          ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_DAY_TOT) ||
          ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_DEC_TOT) ||
          ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_ACT_TOT))
        {
          /* Loop over steps and check cumulative */
          
          for (n = 1; n < ntot; n++)
            if (RDB[ptr + n] <= RDB[ptr + n - 1])
              Error(dep, "Steps not in ascending order after %1.5E", 
                    RDB[ptr + n - 1]);
        }

      /* Check values */
      
      for (n = 0; n < ntot; n++)
        if (RDB[ptr + n] <= 0.0)
          Error(dep, "Invalid step size %E", RDB[ptr + n]);

      /* Next */

      dep = NextItem(dep);
    }

  /* Link branches */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Check pointer to branch */

      if ((long)RDB[dep + DEP_HIS_PTR_BRANCH] > VALID_PTR)
        {
          /* Find branch */

          bra = (long)RDB[DATA_PTR_BRA0];
          while (bra > VALID_PTR)
            {
              /* Compare */

              if (CompareStr(dep + DEP_HIS_PTR_BRANCH, bra + DEP_BRA_PTR_NAME))
                break;
              
              /* Next branch */

              bra = NextItem(bra);
            }

          /* Check pointer */

          if (bra < VALID_PTR)
            Error(dep, "Branch %s is not defined", 
                  GetText(dep + DEP_HIS_PTR_BRANCH));
          else
            WDB[dep + DEP_HIS_PTR_BRANCH] = (double)bra;
        }
      
      /* Next */

      dep = NextItem(dep);
    }
}

/*****************************************************************************/
