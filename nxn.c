/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nxn.c                                          */
/*                                                                           */
/* Created:       2011/03/07 (JLe)                                           */
/* Last modified: 2017/03/21 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Handles (n,xn) reactions                                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Nxn:"

/*****************************************************************************/

void Nxn(long rea, long part, double *E0, double x, double y, double z, 
         double *u0, double *v0, double *w0, double wgt1, double *wgt2, 
         double t, double *dE, long id)
{
  long n, new, nuc, ptr;
  double E, u, v, w, mul;

  /* Get multiplication */

  mul = RDB[rea + REACTION_WGT_F];
  CheckValue(FUNCTION_NAME, "mul", "", mul, 2.0, 4.0);

  /* Check mode */

  if ((long)RDB[DATA_OPT_IMPL_NXN] == YES)
    {
      /***********************************************************************/

      /***** Implicit (nxn) **************************************************/

      /* Perform inelastic scattering on incident neutron */

      InelasticScattering(rea, E0, u0, v0, w0, id);

      /* Multiply weight */

      *wgt2 = wgt1*mul;

      /* Calculate change in energy for deposition */

      *dE = *dE - *E0*mul;

      /* Score particle balance */

      ptr = (long)RDB[RES_N_BALA_SRC];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf(*wgt2 - wgt1, 1.0, ptr, id, -1, BALA_N_SRC_NXN, 1);

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Analog (n,xn) ***************************************************/
     
      /* This mode is disabled because it doesn't work with scattering */
      /* matrixes */

      if ((long)RDB[DATA_OPTI_GC_CALC] == YES)
        Die(FUNCTION_NAME, "Analog (n,xn) does not work with scattering mtx");

      /* Loop over product nuclides */

      for (n = 0; n < (long)mul - 1; n++)
        {
          /* Duplicate incident neutron */
          
          new = DuplicateParticle(part, id);
          
          /* Copy energy and angular variables */

          E = *E0;
          u = *u0;
          v = *v0;
          w = *w0;

          /* Perform inelastic scattering on neutron */

          InelasticScattering(rea, &E, &u, &v, &w, id);
          
          /* Adjust minimum and maximum */

          if (E < 1.0000001*RDB[DATA_NEUTRON_EMIN])
            E = 1.000001*RDB[DATA_NEUTRON_EMIN];
          else if (E > 0.999999*RDB[DATA_NEUTRON_EMAX])
            E = 0.999999*RDB[DATA_NEUTRON_EMAX];

          /* Put variables */

          WDB[new + PARTICLE_X] = x;
          WDB[new + PARTICLE_Y] = y;
          WDB[new + PARTICLE_Z] = z;
          
          WDB[new + PARTICLE_U] = u;
          WDB[new + PARTICLE_V] = v;
          WDB[new + PARTICLE_W] = w;
          
          WDB[new + PARTICLE_E] = E;
          WDB[new + PARTICLE_WGT] = wgt1;
          WDB[new + PARTICLE_T] = t;
          
          /* Put neutron in que */

          ToQue(new, id);

          /* Calculate change in energy for deposition */

          *dE = *dE - E;

          /* Score particle balance */

          ptr = (long)RDB[RES_N_BALA_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(1.0, 1.0, ptr, id, -1, BALA_N_SRC_NXN, 0);
          AddBuf(wgt1, 1.0, ptr, id, -1, BALA_N_SRC_NXN, 1);
        }

      /* Perform inelastic scattering on incident neutron */

      InelasticScattering(rea, E0, u0, v0, w0, id);
 
      /* Calculate change in energy for deposition */

      *dE = *dE - *E0;

      /* Weight is preserved */

      *wgt2 = wgt1;

      /***********************************************************************/
    }
  
  /* Get nuclide pointer */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
  
  /* Check ZAI */
  
  if (((long)RDB[nuc + NUCLIDE_ZAI] == 922330) || 
      ((long)RDB[nuc + NUCLIDE_ZAI] == 922350) || 
      ((long)RDB[nuc + NUCLIDE_ZAI] == 942390) || 
      ((long)RDB[nuc + NUCLIDE_ZAI] == 942410))
    {
      /* Score analog fissile loss rate */
      
      ptr = (long)RDB[RES_ANA_CONV_RATIO];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, wgt1, ptr, id, 2);
    }

  /* Check ZAI */
  
  if ((((long)RDB[nuc + NUCLIDE_ZAI] == 922340) || 
       ((long)RDB[nuc + NUCLIDE_ZAI] == 922360) || 
       ((long)RDB[nuc + NUCLIDE_ZAI] == 942400) || 
       ((long)RDB[nuc + NUCLIDE_ZAI] == 942420)) && 
      ((long)RDB[rea + REACTION_MT] == 16))
    {
      /* Score analog fissile production rate */
      
      ptr = (long)RDB[RES_ANA_CONV_RATIO];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, wgt1, ptr, id, 1);
    }
}

/*****************************************************************************/
