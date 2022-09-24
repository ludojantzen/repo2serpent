/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writecimomfluxes.c                             */
/*                                                                           */
/* Created:       2014/07/23 (VVa)                                           */
/* Last modified: 2014/07/23 (VVa)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Writes material-wise one group burn fluxes to momfluxes.m    */
/*              in Corrector Iteration (SIE) calculations                    */
/*                                                                           */
/* Comments: - Used for debugging etc.                                       */
/*           - Printed fluxes are NOT SIE averaged but they are relaxed      */
/*             inside the current transport cycle if available               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteCIMomFluxes:"

/*****************************************************************************/

void WriteCIMomFluxes(long step)
{
#ifdef DEBUG

  char tmpstr[MAX_STR];
  long mat, pc, ptr;
  double flx;
  FILE *fp;

  if (!((RDB[DATA_BURN_SIE] == YES) || (RDB[DATA_BURN_CI_MAXI] > 1.0)))
    return;

  /* Print file name */

  sprintf(tmpstr, "momfluxes.m");

  /* Open file for writing or appending */

  if(RDB[DATA_BURN_STEP_TOT] < 2)
    {
    if ((fp = fopen(tmpstr, "w")) == NULL) 
      Die(FUNCTION_NAME, "Unable to open file for writing");
    }
  else
    {
    if ((fp = fopen(tmpstr, "a")) == NULL) 
      Die(FUNCTION_NAME, "Unable to open file for writing");

    }

  /* Get predictor/corrector mode */

  pc = (long)RDB[DATA_BURN_STEP_PC];

  /* Loop over materials to write burn fluxes */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {

      /* Skip non-burnable materials */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        {

          /* Get pointer to burn-flux */

          ptr = (long)RDB[mat + MATERIAL_PTR_BURN_FLUX];

          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get value of burn flux */
          
          flx = Mean(ptr, 0);

          /* Use relaxed flux if available (already divided by volume) */
          
          if(RDB[mat + MATERIAL_BURN_FLUX_REL] != 0)
            {
              flx = RDB[mat + MATERIAL_BURN_FLUX_REL];
              flx = flx*RDB[mat + MATERIAL_VOLUME];
            }

          /* Write flux to disk */

          if (((pc == PREDICTOR_STEP) && (step == 0)) || (pc == CORRECTOR_STEP))
            fprintf(fp,"%s_flux(%ld)=%.12e;\n",GetText(mat+MATERIAL_PTR_NAME),
                    (long)RDB[DATA_BURN_STEP_TOT], flx);
        
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Close file */

  fclose(fp);

#endif
}
