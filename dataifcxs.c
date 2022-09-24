/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dataifcxs.c                                    */
/*                                                                           */
/* Created:       2018/06/05 (VVa)                                           */
/* Last modified: 2018/06/21 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns sum of macroscopic cross sections for nuclides       */
/*              coming from interface for a certain material.                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DataIFCXS:"

/*****************************************************************************/

double DataIFCXS(long mat, double E, long mt, long ng, long id)
{
  long nuc, rea, arr, ifc, mul, ty, reamt;
  double xs, adens, sum;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Reset cross section */

  sum = 0.0;

  /* Avoid compiler warning */

  rea = -1;

  /* Get pointer to material data ifc array */

  if ((arr = (long)RDB[mat + MATERIAL_PTR_DATAIFC_ARR]) > VALID_PTR)
    {

      while ((ifc = (long)RDB[arr]) > VALID_PTR)
        {
          /* Skip data interfaces that do not provide atomic density */

          if ((long)RDB[ifc + DATAIFC_DATA_TYPE] != DATA_TYPE_IFC_ADENS)
            {
              /* Next interface */

              arr++;

              continue;
            }

          /* Pointer to nuclide */

          nuc = (long)RDB[ifc + DATAIFC_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Get atomic density */

          adens = DataIFCAdens(mat, nuc, id);

          /* Get pointer to reaction */

          if ((mt == MT_MACRO_INLPRODXS) || (mt == MT_MACRO_FISSXS))
            {
              /* Needs to be calculated from partials */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  /* Get mt, ty and multiplication */

                  reamt = (long)RDB[rea + REACTION_MT];
                  ty = (long)RDB[rea + REACTION_TY];
                  mul = (long)RDB[rea + REACTION_WGT_F];

                  /* Check reaction type */

                  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
                    {
                      /* Fission */

                      if (mt == MT_MACRO_FISSXS)
                        if (((reamt > 17) && (reamt < 22)) || (reamt == 38))
                          {
                            if (RDB[DATA_OPTI_MG_MODE] == (double)YES)
                              {
                                /* Multi-group majorant mode */

                                /* Get pointer to MGXS */

                                xs = MGXS(rea, E, ng);

                              }
                            else
                              {
                                /* Continuous energy mode */

                                /* Get microscopic cross section */

                                xs = MicroXS(rea, E, id);
                              }

                            /* Add to scattering production cross section */

                            sum += xs*adens;
                          }

                      /* Inelastic scattering production */

                      if (mt == MT_MACRO_INLPRODXS)
                        if (mul >= 2)
                          {
                            /* Get cross section */

                            if (RDB[DATA_OPTI_MG_MODE] == (double)YES)
                              {
                                /* Multi-group majorant mode */

                                /* Get pointer to MGXS */

                                xs = MGXS(rea, E, ng);

                              }
                            else
                              {
                                /* Continuous energy mode */

                                /* Get microscopic cross section */

                                xs = MicroXS(rea, E, id);
                              }

                            /* Add to scattering production cross section */

                            sum += xs*(mul - 1.0)*adens;

                          }
                    }

                  /* Next reaction */

                  rea = NextItem(rea);
                }

            }
          else
            {
              if (mt == MT_MACRO_TOTXS)
                rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
              else if (mt == MT_MACRO_ABSXS)
                rea = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];
              else if (mt == MT_MACRO_ELAXS)
                rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
              else if (mt == MT_MACRO_TMP_MAJORANTXS)
                {
                  /* Get pointer to reaction data */

                  rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
                  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

                  /* Get pointer to separate majorant if not using MGXS */
                  /* Multi-group majorant mode stores TMP_MAJORANT directly */
                  /* to NUCLIDE_PTR_TOTXS */

                  if (RDB[DATA_OPTI_MG_MODE] == (double)NO)
                    if ((rea = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT]) < VALID_PTR)
                      Die(FUNCTION_NAME, "Temperature majorant did not exist");
                }
              else
                Die(FUNCTION_NAME, "Invalid reaction mode %ld", mt);

              /* Check pointer */

              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              if (RDB[DATA_OPTI_MG_MODE] == (double)YES)
                {
                  /* Multi-group majorant mode */

                  /* Get pointer to MGXS */

                  xs = MGXS(rea, E, ng);

                }
              else
                {
                  /* Continuous energy mode */

                  /* Get microscopic cross section */

                  xs = MicroXS(rea, E, id);
                }

              /* Add to total */

              sum = sum + xs*adens;
            }

          /* Increment pointer */

          arr++;
        }
    }

  return sum;
}

/*****************************************************************************/
