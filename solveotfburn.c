/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : solveotfburn.c                                 */
/*                                                                           */
/* Created:       2018/03/25 (JLe)                                           */
/* Last modified: 2018/03/25 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Solves on-the-fly burnup problem.                            */
/*                                                                           */
/* Comments: - Tällä hetkellä pelkkä testirutiini joka katsoo että muut      */
/*             osat toimii.                                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SolveOTFBurn:"

/*****************************************************************************/

void SolveOTFBurn()
{
  long loc0, ptr, iso, mat, id, nuc;
  double old, val;

  /* Check OTF burn mode */

  if ((long)RDB[DATA_OTF_BURN_MODE] == NO)
    return;

  /* Check predictor-corrector calculation */

  if (((long)RDB[DATA_BURN_CORR_TYPE] != CORR_TYPE_NONE) &&
      ((long)RDB[DATA_BURN_STEP_PC] != CORRECTOR_STEP))
    return;

  /* Reset maximum concentrations */

  ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
  while (ptr > VALID_PTR)
    {
      /* Check type */
      
      if ((long)RDB[ptr + MAJORANT_EXTRA_TYPE] == MAJORANT_EXTRA_OTF_BURN)
        WDB[ptr + MAJORANT_EXTRA_FRAC] = 0.0;
      
      /* Pointer to next */

      ptr = NextItem(ptr);
    }

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_OTF_BURN0];
  while (loc0 > VALID_PTR)
    {
      /* Reset concentration */

      WDB[loc0 + OTF_BURN_ADENS] = 0.0;

      /* Next material */

      loc0 = NextItem(loc0);
    }

  /* Loop materials and set OpenMP id */

  id = 0;

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check pointer to data */
      
      if ((long)RDB[mat + MATERIAL_PTR_OTF_BURN] > VALID_PTR)
        {
          /* Put OpenMP id */

          WDB[mat + MATERIAL_OMP_ID] = (double)(id++);

          /* Check id */
          
          if (id == (long)RDB[DATA_OMP_MAX_THREADS])
            id = 0;
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Start timer */

  StartTimer(TIMER_OMP_PARA);
  
  /* Loop materials  */

#ifdef OPEN_MP
#pragma omp parallel private (mat)
#endif
  {
    mat = (long)RDB[DATA_PTR_M0];
    while (mat > VALID_PTR)
      {
        /* Get pointer to data and solve TTA */
        
        if ((long)RDB[mat + MATERIAL_PTR_OTF_BURN] > VALID_PTR)
          if ((long)RDB[mat + MATERIAL_OMP_ID] == OMP_THREAD_NUM)
            OTFBurnTTA(mat);
        
        /* Next material */
        
        mat = NextItem(mat);
      }
  }

  /* Stop timer */

  StopTimer(TIMER_OMP_PARA);

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_OTF_BURN0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to composition */

      iso = (long)RDB[loc0 + OTF_BURN_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Replace concentration */

      old = RDB[iso + COMPOSITION_ADENS];
      WDB[iso + COMPOSITION_ADENS] = RDB[loc0 + OTF_BURN_ADENS];
      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      /* Pointer to material */

      mat = (long)RDB[loc0 + OTF_BURN_PTR_MAT];
      CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
      
      /* Adjust material density */

      val = RDB[iso + COMPOSITION_ADENS];
      WDB[mat + MATERIAL_ADENS] = RDB[mat + MATERIAL_ADENS] - old + val;

      /* Pointer to majorant */

      if ((ptr = (long)RDB[loc0 + OTF_BURN_PTR_MAJ]) > VALID_PTR)
        {
          /* Compare fraction */
          
          if (val > RDB[ptr + MAJORANT_EXTRA_FRAC])
            WDB[ptr + MAJORANT_EXTRA_FRAC] = val;
        }

      /* Next */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
