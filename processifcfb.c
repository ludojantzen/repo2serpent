/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifcfb.c                                 */
/*                                                                           */
/* Created:       2015/02/02 (VVa)                                           */
/* Last modified: 2018/01/25 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Processes fuel behavior multi-physics interfaces             */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCFB:"

/*****************************************************************************/

void ProcessIFCFB(long loc0, long update)
{
  long loc1, loc2, mat, ptr, ptr0, uni, nst, reg, cell, axi, surf, ang, m;
  long nr, i, n;
  double f, r1, r0, T, T0, phi;

  /********************************************************************/

  /* Check if this is a 2D calculation and adjust the axial limits */

  if ((long)RDB[DATA_GEOM_DIM] == 2)
    {

      /* Pointer to fuel rod interface definition */

      loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop over all pins */

      while (loc1 > VALID_PTR)
        {

          /* Check that the user has put in a single */
          /* axial zone */

          if ((long)RDB[loc1 + IFC_FUEP_N_AX] > 1)
            Die(FUNCTION_NAME,
                "Multiple axial segments in 2D geometry?");

          /* Get pointer to the output limit structure */

          loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Only one axial output bin */

          WDB[loc2 + FUEP_NZ] = (double)1;

          /* Put axial output limits */

          WDB[loc2 + FUEP_ZMIN] = -INFTY;
          WDB[loc2 + FUEP_ZMAX] = INFTY;

          /* Limits for axial flux binning */

          loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FLIM];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Only one axial output bin */

          WDB[loc2 + FUEP_NZ] = (double)1;

          /* Put axial output limits */

          WDB[loc2 + FUEP_ZMIN] = -INFTY;
          WDB[loc2 + FUEP_ZMAX] = INFTY;

          /* Put limits for the axial input zone */

          axi = (long)RDB[loc1 + IFC_FUEP_PTR_AX];
          CheckPointer(FUNCTION_NAME, "axi", DATA_ARRAY, axi);

          /* Put lower limit */

          WDB[axi + IFC_FUEP_AX_ZMIN] = -INFTY;

          /* Put upper limit */

          WDB[axi + IFC_FUEP_AX_ZMAX] = INFTY;

          /* Next interface rod */

          loc1 = NextItem(loc1);
        }
    }


  /* This part is not executed when updating the interface */
  /* It is only for setting the TMS limits that currently cannot be */
  /* modified during the simulation*/

  if (!update)
    {

      fprintf(outp, "\nFuel behavior interface \"%s\":\n\n",
              GetText(loc0 + IFC_PTR_INPUT_FNAME));

      /* Pointer to structure */

      loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop */

      while (loc1 > VALID_PTR)
        {

          /* Loop over rod segments to link interface to universes */

          for (m = 0; m < (long)RDB[loc1 + IFC_FUEP_N_UNI]; m++)
            {
              /* Get pointer to ifc nests */

              ptr = (long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

              /* Find pin universe */

              uni = (long)RDB[DATA_PTR_U0];
              while (uni > VALID_PTR)
                {
                  /* Compare */

                  if (CompareStr(uni + UNIVERSE_PTR_NAME,
                                 ptr + m))
                    break;

                  /* Next */

                  uni = NextItem(uni);
                }

              /* Check pointer */

              if (uni < VALID_PTR)
                Error(loc0, "Universe %s does not exist",
                      GetText((long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST] + m));

              /* Check universe type (tää testaa vaan että on nesti, ei */
              /* sitä onko pinnan tyyppinä sylinteri) */

              if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_NEST)
                Error(loc0, "Universe %s is not pin type",
                      GetText(uni + UNIVERSE_PTR_NAME));

              /* Check that universe is not associated with another ifc */

              if ((long)RDB[uni + UNIVERSE_PTR_IFC_FUEP] > VALID_PTR)
                Error(loc0, "Multiple interfaces for universe %s",
                      GetText(uni + UNIVERSE_PTR_NAME));

              /* Put pointers */

              WDB[loc1 + IFC_FUEP_PTR_UNI] = (double)uni;
              WDB[uni + UNIVERSE_PTR_IFC_FUEP] = (double)loc1;

              /* Get pointer to nest structure */

              nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
              CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

              axi = (long)RDB[loc1 + IFC_FUEP_PTR_AX];

              /* Loop over axial zones to set TMS Min/Max temperatures for nest materials */
              /* Pilkoin näitä tarkistuksia ja prosessointia vähän useampaan eri looppiin */
              /* Selkeyden vuoksi 9-Jul-13 (VVa) */

              while (axi > VALID_PTR)
                {
                  /* Pointer to angular zones */

                  ang = (long)RDB[axi + IFC_FUEP_AX_PTR_ANG];
                  CheckPointer(FUNCTION_NAME, "(ang)", DATA_ARRAY, ang);

                  /* Loop over angular zones */

                  while (ang > VALID_PTR)
                    {
                      /* Get pointer to nest regions */

                      reg = (long)RDB[nst + NEST_PTR_REGIONS];

                      while (reg > VALID_PTR)
                        {
                          r0 = 0.0;
                          r1 = 0.0;

                          /* Pointer to outer surface*/

                          surf = (long)RDB[reg + NEST_REG_PTR_SURF_IN];

                          /* Do not process the outermost (non-bounded) region */
                          /* which is usually coolant or such */

                          if (surf < VALID_PTR)
                            {
                              reg = NextItem(reg);
                              continue;
                            }

                          /* Outer surface exists, get parameters */

                          ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];

                          /* Outer radius */

                          r1 = RDB[ptr + 2];

                          /* Pointer to inner surface */

                          surf = (long)RDB[reg + NEST_REG_PTR_SURF_OUT];

                          /* Get inner radius if available*/

                          if (surf > VALID_PTR)
                            {
                              ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
                              r0 = RDB[ptr + 2];
                            }

                          /* Pointer to cell */

                          cell = (long)RDB[reg + NEST_REG_PTR_CELL];
                          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

                          /* Pointer to material */

                          if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                            {

                              loc2 = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];

                              /* Do not check the last node as it is an artificial */
                              /* node added to allow coordinate transformations    */
                              /* even if cladding creeps inwards */

                              for (i = 0; i < (long)RDB[ang + IFC_FUEP_ANG_N_RAD] - 1; i++)
                                {
                                  /* if this node is not in this nest region, */
                                  /* skip NOTE: There might not be nodes in   */
                                  /* every region */

                                  if((RDB[loc2 + i] < r0) || (RDB[loc2 + i] > r1))
                                    continue;

                                  /* Get pointer to EOI temperature array */

                                  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TMP];
                                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                                  /* Get pointer to BOI temperature array */

                                  ptr0 = (long)RDB[ang + IFC_FUEP_ANG_PTR_TMP_BOI];
                                  CheckPointer(FUNCTION_NAME, "(ptr0)", DATA_ARRAY, ptr0);

                                  /* Get larger temperature */

                                  if ((T = RDB[ptr + i]) < RDB[ptr0 + i])
                                    T = RDB[ptr0 + i];

                                  if (T > 0.0)
                                    {

                                      /* Adjust material temperatures */

                                      if (T > RDB[mat + MATERIAL_TMS_TMAX])
                                        WDB[mat + MATERIAL_TMS_TMAX] = T;

                                      if (T < RDB[mat + MATERIAL_TMS_TMIN])
                                        WDB[mat + MATERIAL_TMS_TMIN] = T;

                                    }
                                  /* End of radial for loop */
                                }

                              /* If there is no node in this region, the maximum */
                              /* and minimum temperatures are at the inner/outer */
                              /* surfaces */

                              for (n = 0; n < 2; n++)
                                {
                                  /**********************************************/
                                  /* Check temperature at outer radius (n == 0) */
                                  /* or inner radius (n == 1)                   */
                                  /**********************************************/

                                  loc2 = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];

                                  /* Get node index, surface will be between nodes */
                                  /* i and i+1 */

                                  if (n == 0)
                                    i = SearchArray(&RDB[loc2], r1,
                                                    (long)RDB[ang + IFC_FUEP_ANG_N_RAD]);
                                  else
                                    i = SearchArray(&RDB[loc2], r1,
                                                    (long)RDB[ang + IFC_FUEP_ANG_N_RAD]);

                                  if (i > -1)
                                    {

                                      /* Get pointer to EOI temperature array */

                                      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TMP];
                                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                                      /* Get pointer to BOI temperature array */

                                      ptr0 = (long)RDB[ang + IFC_FUEP_ANG_PTR_TMP_BOI];
                                      CheckPointer(FUNCTION_NAME, "(ptr0)", DATA_ARRAY, ptr0);

                                      /* Get larger temperature */

                                      if ((T = RDB[ptr + i + 1]) < RDB[ptr0 + i + 1])
                                        T = RDB[ptr0 + i + 1];

                                      /* Interpolate temperature if IFC_TYPE_FPIP */

                                      if ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FPIP)
                                        {
                                          /* Both nodes at same location? */

                                          if (RDB[loc2 + i + 1] - RDB[loc2 + i] == 0)
                                            {
                                              /* Get larger temperature at inner radius */

                                              if ((T = RDB[ptr + i]) < RDB[ptr0 + i])
                                                T = RDB[ptr0 + i];
                                            }
                                          else
                                            {
                                              /* Interpolate EOI temperatures */

                                              T = RDB[ptr + i] +
                                                (RDB[ptr + i + 1] - RDB[ptr + i])*
                                                (r1 - RDB[loc2 + i])/
                                                (RDB[loc2 + i+1] - RDB[loc2 + i]);

                                              /* Interpolate BOI temperatures */

                                              T0 = RDB[ptr0 + i] +
                                                (RDB[ptr0 + i + 1] - RDB[ptr0 + i])*
                                                (r1 - RDB[loc2 + i])/
                                                (RDB[loc2 + i+1] - RDB[loc2 + i]);

                                              /* Get larger temperature */

                                              if (T < T0)
                                                T = T0;
                                            }
                                        }

                                      /* Check resulting temperature */

                                      if (T < 0)
                                        Die(FUNCTION_NAME,
                                            "Negative temperature when interpolating between: %E %E or %E %E",
                                            RDB[ptr + i + 1], RDB[ptr + i],
                                            RDB[ptr0 + i + 1], RDB[ptr0 + i]);

                                      /* Adjust TMS temperatures */

                                      if (T > RDB[mat + MATERIAL_TMS_TMAX])
                                        WDB[mat + MATERIAL_TMS_TMAX] = T;

                                      if (T < RDB[mat + MATERIAL_TMS_TMIN])
                                        WDB[mat + MATERIAL_TMS_TMIN] = T;

                                    }
                                }
                            }

                          /* Next region */

                          reg = NextItem(reg);
                        }

                      /* Next angular zone */

                      ang = NextItem(ang);
                    }

                  /* Next axial zone */

                  axi = NextItem(axi);
                }

            }
          /* Next pin */

          loc1 = NextItem(loc1);
        }

      /***********************************************/
      /* Put final limits and set TMS mode on or off */
      /***********************************************/

      /* Pointer to structure */

      loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop */

      while (loc1 > VALID_PTR)
        {

          for (m = 0; m < (long)RDB[loc1 + IFC_FUEP_N_UNI]; m++)
            {
              /* Get pointer to ifc nests */

              ptr = (long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

              /* Find pin universe */

              uni = (long)RDB[DATA_PTR_U0];
              while (uni > VALID_PTR)
                {
                  /* Compare */

                  if (CompareStr(uni + UNIVERSE_PTR_NAME,
                                 ptr + m))
                    break;

                  /* Next */

                  uni = NextItem(uni);
                }

              /* Check pointer */

              if (uni < VALID_PTR)
                Error(loc0, "Universe %s does not exist",
                      GetText((long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST] + m));

              /* Get pointer to nest structure */

              nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
              CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

              /* Get pointer to nest regions */

              reg = (long)RDB[nst + NEST_PTR_REGIONS];
              CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

              while (reg > VALID_PTR)
                {

                  /* Pointer to cell */

                  cell = (long)RDB[reg + NEST_REG_PTR_CELL];
                  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

                  /* Pointer to material */

                  if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                    {

                      /* Put flag */

                      if (RDB[mat + MATERIAL_TMS_TMIN] <=
                           RDB[mat + MATERIAL_TMS_TMAX] )
                          WDB[mat + MATERIAL_USE_IFC] = (double)YES;

                      /* Set on-the-fly Doppler-broadening mode */

                      if ((RDB[mat + MATERIAL_TMS_TMIN] <
                           RDB[mat + MATERIAL_TMS_TMAX] ))
                        {

                          /* Check that material doppler temperature is not already */
                          /* set by tmp-card */

                          if (RDB[mat + MATERIAL_DOPPLER_TEMP] >= 0)
                            Error(mat, "Material temperature set by tmp-card but a temperature distribution is also given by interface %s", GetText(loc0 + IFC_PTR_INPUT_FNAME));

                          fprintf(outp, " - Material %s TMS limits: %.2f K and %.2f K\n",
                                  GetText(mat + MATERIAL_PTR_NAME),
                                  RDB[mat + MATERIAL_TMS_TMIN], RDB[mat + MATERIAL_TMS_TMAX]);

                          /* Set TMS-mode on */

                          WDB[mat + MATERIAL_TMS_MODE] = (double)YES;

                          /* Set Doppler-preprocessor off */

                          WDB[mat + MATERIAL_DOPPLER_TEMP] = -1.0;

                        }
                      else if (RDB[mat + MATERIAL_TMS_TMIN] == RDB[mat + MATERIAL_TMS_TMAX])
                        {

                          /* Check that material doppler temperature is not already */
                          /* set by tmp-card */

                          if (RDB[mat + MATERIAL_DOPPLER_TEMP] >= 0)
                            Error(mat, "Material temperature set by tmp-card but a temperature distribution is also given by interface %s", GetText(loc0 + IFC_PTR_INPUT_FNAME));

                          WDB[mat + MATERIAL_DOPPLER_TEMP] = RDB[mat + MATERIAL_TMS_TMIN];

                          /* Set Doppler preprocessor on */

                          WDB[DATA_USE_DOPPLER_PREPROCESSOR] = (double)YES;


                        }
                    }

                  reg = NextItem(reg);
                }
            }

          /* Next pin */

          loc1 = NextItem(loc1);
        }

    }

  /***************************************************************************/

  /***** Additional processing for fuep type *********************************/

  /* Loop over pins */

  loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
  while (loc1 > VALID_PTR)
    {

      /* Loop over axial zones */

      axi = (long)RDB[loc1 + IFC_FUEP_PTR_AX];

      while (axi > VALID_PTR)
        {
          /* Loop over angular zones */

          ang = (long)RDB[axi + IFC_FUEP_AX_PTR_ANG];

          while (ang > VALID_PTR)
            {

              /* Calculate the parameters for the limiting planes of */
              /* the angular segment */

              phi = RDB[ang + IFC_FUEP_ANG_AMIN];

              if (cos(phi)==0)
                WDB[ang + IFC_FUEP_ANG_CMIN] = -sin(phi)*INFTY;
              else
                WDB[ang + IFC_FUEP_ANG_CMIN] = -sin(phi)/cos(phi);

              phi = RDB[ang + IFC_FUEP_ANG_AMAX];

              if (cos(phi)==0)
                WDB[ang + IFC_FUEP_ANG_CMAX] = -sin(phi)*INFTY;
              else
                WDB[ang + IFC_FUEP_ANG_CMAX] = -sin(phi)/cos(phi);

              /* Get number of radial nodes */

              nr = (long)RDB[ang + IFC_FUEP_ANG_N_RAD];

              /* Loop over radial zones to set square radius and DF */

              for (i = 0; i < nr; i++)
                {

                  /* Calculate square radius */

                  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];
                  WDB[ptr + i] = RDB[ptr + i]*RDB[ptr + i];

                  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];
                  WDB[ptr + i] = RDB[ptr + i]*RDB[ptr + i];

                  /* Calculate density factor */

                  if (i > 0)
                    {

                      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];
                      f = RDB[ptr + i - 1] - RDB[ptr + i];

                      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];
                      f = f/(RDB[ptr + i - 1] - RDB[ptr + i]);

                    }
                  else /* First zone */
                    {

                      f = 1.0;

                      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];

                      if (RDB[ptr + i] != 0.0)
                        {
                          f = RDB[ptr + i];

                          ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];
                          f = f/RDB[ptr + i];

                        }
                    }

                  /* Check value */

                  if (f > 1.0)
                    {
                      /*Warn(FUNCTION_NAME, "Density factor larger than 1.0, setting 1.0, outradius of zone %E",
                        sqrt(RDB[ptr + i]));*/
                      f = 1.0;
                    }

                  /* Put density factor */

                  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_DF];

                  WDB[ptr + i] = f;

                  /* Next radial zone */

                }

              /* In time dependent coupled calculation, we'll create separate lists for BOI */
              /* To allow interpolation of temperature/density data inside the interval     */

              if ((RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN) && (update == 0))
                {

                  /* Allocate memory for temperatures */

                  ptr0 = ReallocMem(DATA_ARRAY, nr);

                  /* Put temperature list pointer to memory */

                  WDB[ang + IFC_FUEP_ANG_PTR_TMP_BOI] = (double)ptr0;

                  /* Get pointer to EOI temperature array */

                  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TMP];

                  /* Copy temperature data from EOI to BOI */

                  memcpy(&WDB[ptr0], &RDB[ptr], nr*sizeof(double));

                  /**************************************/

                  /* Allocate memory for density factor */

                  ptr0 = ReallocMem(DATA_ARRAY, nr);

                  /* Put temperature list pointer to memory */

                  WDB[ang + IFC_FUEP_ANG_PTR_DF_BOI] = (double)ptr0;

                  /* Get pointer to EOI density factor array */

                  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_DF];

                  /* Copy density factor data from EOI to BOI */

                  memcpy(&WDB[ptr0], &RDB[ptr], nr*sizeof(double));

                  /****************************************/

                  /* Allocate memory for hot radius array */

                  ptr0 = ReallocMem(DATA_ARRAY, nr);

                  /* Put hot radius array pointer to memory */

                  WDB[ang + IFC_FUEP_ANG_PTR_HOT_R2_BOI] = (double)ptr0;

                  /* Get pointer to EOI hot radius factor array */

                  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];

                  /* Copy hot square radius data from EOI to BOI */

                  memcpy(&WDB[ptr0], &RDB[ptr], nr*sizeof(double));

                }

              /* Next angular zone */

              ang = NextItem(ang);
            }

          /* Next axial zone */

          axi = NextItem(axi);

        }
      /* Next pin */

      loc1 = NextItem(loc1);
    }

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
