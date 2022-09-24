/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : complexrea.c                                   */
/*                                                                           */
/* Created:       2014/10/16 (JLe)                                           */
/* Last modified: 2017/02/23 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Handles reaction MT 5, which includes several channels and   */
/*              neutron multiplication                                       */
/*                                                                           */
/* Comments: - Weight multiplication handled implicitly                      */
/*           - Contribution to energy deposition may be incorrect            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ComplexRea:"

/*****************************************************************************/

void ComplexRea(long rea, long part, double *E0, double x, double y, double z,
                double *u0, double *v0, double *w0, double wgt1, double *wgt2, 
                double t, double *dE, long id)
{
  long ptr, ne, erg, i, n, new, mul;
  double Emin, Emax, f, r, E, u, v, w;

  /* Get pointer to data */

  ptr = (long)RDB[rea + REACTION_PTR_MULT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of energies */

  ne = (long)RDB[ptr++];
  CheckValue(FUNCTION_NAME, "(ne)", "", ne, 2, 10000);

  /* Get pointer to energy grid */
  
  erg = (long)RDB[ptr++];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check number of energy points */

  if (ne != (long)RDB[erg + ENERGY_GRID_NE])
    Die(FUNCTION_NAME, "Mismatch in number of energy groups");

  /* Get minimum and maximum energy */

  Emin = RDB[erg + ENERGY_GRID_EMIN];
  Emax = RDB[erg + ENERGY_GRID_EMAX];

  if ((r = GridFactor(erg, *E0, id)) < 0.0)
   {
      /* Avoid compiler warning */

      i = 0;
      r = 0.0;

      /* Check if energy is above or below limits */

      if (*E0 > Emax)
        {
          i = ne - 2;
          r = 1.0;
        }
      else if (*E0 < Emin)
        {
          i = 0;
          r = 0.0;
        }
      else
        Die(FUNCTION_NAME, "wtf?!?");
    }
  else
    {
      /* Get bin index and adjust factor */

      i = (long)r;
      r = r - (double)i;
    }

  /* Check values */

  CheckValue(FUNCTION_NAME, "i", "", i, 0, ne - 2);
  CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 1.0);  

  /* Get multiplication */

  f = RDB[ptr + i]*(1.0 - r) + RDB[ptr + i + 1]*r;
  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 25.0);

  /* Analog or implicit treatment */

  if ((long)RDB[DATA_OPT_IMPL_NXN] == YES)
    {
      /***********************************************************************/

      /***** Implicit method *************************************************/
      
      /* Adjust weight */
      
      if ((*wgt2 = f*wgt1) > 0.0)
        {
          /* Perform inelastic scattering on incident neutron */
          
          InelasticScattering(rea, E0, u0, v0, w0, id);
          
          /* Calculate change in energy for deposition */
          
          *dE = *dE - *E0*f;

          /* Score particle balance */

          ptr = (long)RDB[RES_N_BALA_SRC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(*wgt2 - wgt1, 1.0, ptr, id, -1, BALA_N_SRC_NXN, 1);
        }
      
      /***********************************************************************/
    }
  else if (f > 0.0)
    {
      /***********************************************************************/

      /***** Analog method ***************************************************/

      /* This mode is disabled because it doesn't work with scattering */
      /* matrixes */
      
      if ((long)RDB[DATA_OPTI_GC_CALC] == YES)
        Die(FUNCTION_NAME, "Analog (n,xn) does not work with scattering mtx");

      /* Get integer part */

      mul = (long)f;
      
      /* Sample extra neutron */

      if (RandF(id) < f - (double)mul)
        mul++;

      /* Loop over product nuclides */

      for (n = 0; n < mul - 1; n++)
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
 
     /* Weight is preserved */

      *wgt2 = wgt1;

      /***********************************************************************/
    }
}

/*****************************************************************************/
