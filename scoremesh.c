/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoremesh.c                                    */
/*                                                                           */
/* Created:       2011/03/25 (JLe)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores distributions for mesh plot                           */
/*                                                                           */
/* Comments: - Collision point ja collision weight ei tarkista hiukkasen     */
/*             tyyppiä                                                       */
/*                                                                           */
/*           - dE = 0  --> kutsuttu score.c:stä                              */
/*             dE > 0  --> kutsuttu collision.c:stä                          */
/*             dE = -1 --> kutsuttu samplesrcpoint.c:stä                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreMesh:"

/*****************************************************************************/

void ScoreMesh(long part, long mat, double flx0, double dE, double x0,
               double y0, double z0, double E, double t, double wgt, double g,
               long id)
{
  long mpl, ax, mode, ptr, rea, det, type, evn;
  double x, y, z, xs, H, f, T, majorant, V, flx, wgt0;

  /* Mesh plots are not scored in void */

  if (mat < VALID_PTR)
    return;

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Loop over plots */

  mpl = (long)RDB[DATA_PTR_MPL0];
  while(mpl > VALID_PTR)
    {
      /* Get axis */

      ax = (long)RDB[mpl + MPL_AX];

      /* Check is and set coordinates */

      if (ax == 1)
        {
          /* distribution in yz-plane */

          x = y0;
          y = z0;
          z = x0;
        }
      else if (ax == 2)
        {
          /* distribution in xz-plane */

          x = x0;
          y = z0;
          z = y0;
        }
      else
        {
          /* distribution in xy-plane */

          x = x0;
          y = y0;
          z = z0;
        }

      /* Get volume for cylindrical type */

      if (ax == 4)
        {
          /* Pointer to mesh */

          ptr = (long)RDB[mpl + MPL_PTR_VAL1];
          CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

          /* Get volume */

          V = MeshCellVol(ptr, x, y, z);
        }
      else
        V = 1.0;

      /* Get mode */

      mode = (long)RDB[mpl + MPL_TYPE];

      /* Check mode */

      if (mode == MPL_TYPE_FLUXPOW)
        {
          /*******************************************************************/

          /**** Thermal flux / fission rate mesh *****************************/

          /* Check particle type */

          if (type == PARTICLE_TYPE_NEUTRON)
            {
              /* Check material pointer */

              CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

              /* Fission rate */

              if ((rea = (long)RDB[mat + MATERIAL_PTR_FISSXS]) > VALID_PTR)
                {
                  /* Get total fission cross section */

                  xs = g*MacroXS(rea, E, id);

                  /* Pointer to data */

                  ptr = (long)RDB[mpl + MPL_PTR_VAL1];
                  CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

                  /* Add value to mesh */

                  AddMesh(ptr, xs*wgt*flx0/V, x, y, z, id);
                }

              /* Thermal flux */

              else if (E < 0.625E-6)
                {
                  /* Pointer to data */

                  ptr = (long)RDB[mpl + MPL_PTR_VAL2];
                  CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);

                  /* Add value to mesh */

                  AddMesh(ptr, wgt*flx0/V, x, y, z, id);
                }
            }

          /*******************************************************************/
        }
      else if (mode == MPL_TYPE_FLUXTEMP)
        {
          /*******************************************************************/

          /**** Thermal flux / temperature mesh ******************************/

          /* Check particle type */

          if (type == PARTICLE_TYPE_NEUTRON)
            {
              /* Get temperature */

              if ((T = GetTemp(mat, id)) > 0.0)
                {
                  /* Pointer to data */

                  ptr = (long)RDB[mpl + MPL_PTR_VAL1];
                  CheckPointer(FUNCTION_NAME, "(ptr3)", DATA_ARRAY, ptr);

                  /* Add value to mesh */

                  AddMesh(ptr, T/V, x, y, z, id);

                  /* Pointer to weight distribution */

                  ptr = (long)RDB[mpl + MPL_PTR_DIV1];
                  CheckPointer(FUNCTION_NAME, "(ptr4)", DATA_ARRAY, ptr);

                  /* Add value to mesh */

                  AddMesh(ptr, 1.0/V, x, y, z, id);
                }

              /* Thermal flux */

              else if (E < 0.625E-6)
                {
                  /* Pointer to data */

                  ptr = (long)RDB[mpl + MPL_PTR_VAL2];
                  CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);

                  /* Add value to mesh */

                  AddMesh(ptr, wgt*flx0/V, x, y, z, id);

                  /* Pointer to weight distribution */

                  ptr = (long)RDB[mpl + MPL_PTR_DIV2];
                  CheckPointer(FUNCTION_NAME, "(ptr6)", DATA_ARRAY, ptr);

                  /* Add value to mesh */

                  AddMesh(ptr, 1.0/V, x, y, z, id);
                }
            }

          /*******************************************************************/
        }
      else if (mode == MPL_TYPE_DT_NEFF)
        {
          /*******************************************************************/

          /**** Delta-tracking efficiency for neutrons ***********************/

          /* Check particle type */

          if (type == PARTICLE_TYPE_NEUTRON)
            {
              /* Check material pointer */

              CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

              /* Total cross section */

              rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(re)", DATA_ARRAY, rea);
              xs = g*MacroXS(rea, E, id);

              /* Majorant cross section */

              majorant = DTMajorant(type, E, id);

              /* Check limit */

              if (majorant > 0.0)
                if (xs/majorant > RDB[DATA_DT_NTHRESH])
                  {
                    /* Pointer to data */

                    ptr = (long)RDB[mpl + MPL_PTR_VAL1];
                    CheckPointer(FUNCTION_NAME, "(ptr9)", DATA_ARRAY, ptr);

                    /* Add value to mesh */

                    AddMesh(ptr, xs*flx0, x, y, z, id);

                    /* Pointer to data */

                    ptr = (long)RDB[mpl + MPL_PTR_DIV1];
                    CheckPointer(FUNCTION_NAME, "(ptr10)", DATA_ARRAY, ptr);

                    /* Add value to mesh */

                    AddMesh(ptr, majorant*flx0, x, y, z, id);
                  }
            }

          /*******************************************************************/
        }
      else if (mode == MPL_TYPE_DT_GEFF)
        {
          /*******************************************************************/

          /**** Delta-tracking efficiency for photons ************************/

          /* Check particle type */

          if (type == PARTICLE_TYPE_GAMMA)
            {
              /* Check material pointer */

              CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

              /* Total cross section */

              rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
              xs = g*MacroXS(rea, E, id);

              /* Majorant cross section */

              majorant = DTMajorant(type, E, id);

              /* Check limit */

              if (majorant > 0.0)
                if (xs/majorant > RDB[DATA_DT_PTHRESH])
                  {
                    /* Pointer to data */

                    ptr = (long)RDB[mpl + MPL_PTR_VAL1];
                    CheckPointer(FUNCTION_NAME, "(ptr9)", DATA_ARRAY, ptr);

                    /* Add value to mesh */

                    AddMesh(ptr, xs*flx0/V, x, y, z, id);

                    /* Pointer to data */

                    ptr = (long)RDB[mpl + MPL_PTR_DIV1];
                    CheckPointer(FUNCTION_NAME, "(ptr10)", DATA_ARRAY, ptr);

                    /* Add value to mesh */

                    AddMesh(ptr, majorant*flx0/V, x, y, z, id);
                  }
            }

          /*******************************************************************/
        }
      else if (mode == MPL_TYPE_GAMMAHEAT)
        {
          /*******************************************************************/

          /***** Photon energy deposition ************************************/

          /* Check particle type */

          if (type == PARTICLE_TYPE_GAMMA)
            {
              /* Pointer to heating cross section */

              if ((rea = (long)RDB[mat + MATERIAL_PTR_HEATPHOTXS]) > VALID_PTR)
                {
                  /* Get heating value */

                  H = g*PhotonMacroXS(rea, E, id);

                  /* Pointer to data */

                  ptr = (long)RDB[mpl + MPL_PTR_VAL1];
                  CheckPointer(FUNCTION_NAME, "(ptr13)", DATA_ARRAY, ptr);

                  /* Add value to mesh */

                  AddMesh(ptr, wgt*H*flx0/V, x, y, z, id);
                }
            }

          /*******************************************************************/
        }
      else if (mode == MPL_TYPE_COLPT)
        {
          /*******************************************************************/

          /***** Number of collisions ****************************************/

          /* Pointer to data */

          ptr = (long)RDB[mpl + MPL_PTR_VAL1];
          CheckPointer(FUNCTION_NAME, "(ptr14)", DATA_ARRAY, ptr);

          /* Add value to mesh */

          AddMesh(ptr, 1.0/V, x, y, z, id);

          /*******************************************************************/
        }
      else if (mode == MPL_TYPE_COLWGT)
        {
          /*******************************************************************/

          /***** Collided weight *********************************************/

          /* Pointer to data */

          ptr = (long)RDB[mpl + MPL_PTR_VAL1];
          CheckPointer(FUNCTION_NAME, "(ptr15)", DATA_ARRAY, ptr);

          /* Add value to mesh */

          AddMesh(ptr, wgt/V, x, y, z, id);

          /*******************************************************************/
        }
      else if ((mode == MPL_TYPE_DET) || (mode == MPL_TYPE_DET_IMP))
        {
          /*******************************************************************/

          /***** Detector or detector importance *****************************/

          /* Pointer to detector */

          det = (long)RDB[mpl + MPL_PTR_DET];
          CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

          /* Check particle type */

          if (type != (long)RDB[det + DET_PARTICLE])
            {
              /* Next plot */

              mpl = NextItem(mpl);

              /* Cycle loop */

              continue;
            }

          /* Get pointer to response function */

          ptr = (long)RDB[det + DET_PTR_RBINS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Reset flux */

          flx = flx0;

          /* Get value of response function */

          if ((long)RDB[ptr + DET_RBIN_MT] == MT_MACRO_RECOILE)
            {
              /* Analog estimator of scattering recoil energy deposition */
              /* (NOTE: only collision.c passes non-zero value). */

              f = dE*MEV;
              flx = 1.0;
            }
          else if ((long)RDB[ptr + DET_RBIN_MT] == MT_SOURCE_RATE)
            {
              /* Source rate (only samplesrcpoint.c passes g = -1) */

              if (g == -1.0)
                {
                  f = 1.0;
                  flx = 1.0;
                }
              else
                f = 0.0;
            }
          else if (DetBin(det, mat, part, x0, y0, z0, E, t, id) > -1)
            f = DetResponse(det, ptr, part, mat, E, g, id);
          else
            f = 0.0;

          /* Add value to mesh */

          if (mode == MPL_TYPE_DET)
            {
              /* Pointer to data */

              ptr = (long)RDB[mpl + MPL_PTR_VAL1];
              CheckPointer(FUNCTION_NAME, "(ptr16)", DATA_ARRAY, ptr);

              /* Check responce function */

              if (f != 0.0)
                {
                  /* Check importance weighting */

                  if ((long)RDB[det + DET_TYPE] == DETECTOR_TYPE_IMP_WGT)
                    AddMesh(ptr, flx*f/V, x, y, z, id);
                  else
                    AddMesh(ptr, wgt*flx*f/V, x, y, z, id);
                }
            }
          else
            {
              /* Pointer to flux distribution */

              ptr = (long)RDB[mpl + MPL_PTR_DIV1];
              CheckPointer(FUNCTION_NAME, "(ptr17)", DATA_ARRAY, ptr);

              /* Score flux */

              AddMesh(ptr, wgt*flx/V, x, y, z, id);

              /* Check response function */

              if (f != 0.0)
                {
                  /* Pointer to data */

                  ptr = (long)RDB[mpl + MPL_PTR_VAL1];
                  CheckPointer(FUNCTION_NAME, "(ptr18)", DATA_ARRAY, ptr);

                  /* Loop over events */

                  evn = (long)RDB[part + PARTICLE_PTR_EVENTS];
                  while (evn > VALID_PTR)
                    {
                      /* Get weight and flux */

                      wgt0 = RDB[evn + EVENT_WGT];
                      flx0 = RDB[evn + EVENT_FLX];

                      /* Check */

                      if (flx0*wgt0 < ZERO)
                        {
                          /* Pointer to next */

                          evn = NextItem(evn);

                          /* Cycle loop */

                          continue;
                        }

                      /* Get coordinates of previous collision point */

                      if (ax == 1)
                        {
                          /* distribution in yz-plane */

                          x = RDB[evn + EVENT_Y];
                          y = RDB[evn + EVENT_Z];
                          z = RDB[evn + EVENT_X];
                        }
                      else if (ax == 2)
                        {
                          /* distribution in xz-plane */

                          x = RDB[evn + EVENT_X];
                          y = RDB[evn + EVENT_Z];
                          z = RDB[evn + EVENT_Y];
                        }
                      else
                        {
                          /* distribution in xy-plane */

                          x = RDB[evn + EVENT_X];
                          y = RDB[evn + EVENT_Y];
                          z = RDB[evn + EVENT_Z];
                        }

                      /* Score responce function */

                      AddMesh(ptr, wgt*flx*f*wgt0*flx0/V, x, y, z, id);

                      /* Next event */

                      evn = NextItem(evn);
                    }
                }
            }

          /*******************************************************************/
        }
      else if (mode == MPL_TYPE_DENSITY)
        Die(FUNCTION_NAME, "Not supported anymore");
      else
        Die(FUNCTION_NAME, "Invalid mesh type");

      /* Next plot */

      mpl = NextItem(mpl);
    }
}

/*****************************************************************************/
