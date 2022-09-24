/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : samplereaction.c                               */
/*                                                                           */
/* Created:       2011/01/04 (JLe)                                           */
/* Last modified: 2019/11/07 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Samples reaction after collision                             */
/*                                                                           */
/* Comments: - Nuclide and reaction lists must be sorted according to sample */
/*             counter.                                                      */
/*                                                                           */
/*           - Tää on kirjoitettu käytännössä kokonaan uusiksi versioon      */
/*             2.1.13 (2/2013).                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleReaction:"

/*****************************************************************************/

long SampleReaction(long mat, long type, double E, double wgt, long id)
{
  long lst, rea, rls, ptr, i, nuc, ng, TMS, mt;
  double totxs, adens, xs, absxs, f, g, Er, Emin, Emax, T, E0;

  /* Avoid compiler warning */

  g = 1.0;
  T = -1.0;
  E0 = INFTY;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Check energy */

  CheckValue(FUNCTION_NAME, "E", "", E, ZERO, INFTY);

  /* Get TMS mode */

  TMS = (long)RDB[mat + MATERIAL_TMS_MODE];

  /***************************************************************************/

  /***** Sample target nuclide ***********************************************/

  /* Check mode */

  if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
    {
      /***********************************************************************/

      /***** Multi-group mode ************************************************/

      /* Pointer to energy groups */

      ptr = (long)RDB[DATA_COARSE_MG_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Find group index */

      if ((ng = GridSearch(ptr, E)) < 0)
        return -1;

      /* Check particle type */

      if (type == PARTICLE_TYPE_NEUTRON)
        {
          /* Get (macroscopic) total neutron cross section */

          rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          totxs = MGXS(rea, E, ng);

          /* Get mt */

          mt = (long)RDB[rea + REACTION_MT];
          CheckValue(FUNCTION_NAME, "mt", "", mt, -30, -1);

          /* Add poison cross section (miksi tossa on toi MACRO_ABSXS?) */

          if (((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] == YES) ||
              ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] == YES))
            if ((mt == MT_MACRO_TOTXS) || (mt == MT_MACRO_ABSXS) ||
                (mt == MT_MACRO_TMP_MAJORANTXS))
              totxs = totxs + PoisonXS(mat, E, mt, id);

          if ((long)RDB[DATA_ITER_MODE] == ITER_MODE_NUCLIDE)
            if ((mt == MT_MACRO_TOTXS) || (mt == MT_MACRO_ABSXS) ||
                (mt == MT_MACRO_TMP_MAJORANTXS))
              totxs = totxs + IterNucXS(mat, E, mt, ng, id);

          /* Add cross sections from data interface nuclides */

          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_EXT_ADENS_MAT)
            totxs = totxs + DataIFCXS(mat, E, mt, ng, id);

          /* Add cross sections for on-the-fly burnup mode */

          if ((long)RDB[DATA_OTF_BURN_MODE] == YES)
            if ((mt == MT_MACRO_TOTXS) || (mt == MT_MACRO_ABSXS) ||
                (mt == MT_MACRO_TMP_MAJORANTXS))
              totxs = totxs + OTFBurnXS(mat, E, mt, id);

          /* Check implicit capture mode */

          if (RDB[DATA_OPT_IMPL_CAPT] == YES)
            Die(FUNCTION_NAME, "Implicit capture not working in this mode");
        }
      else
        {
          /* Get total gamma cross section */

          rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          totxs = MGXS(rea, E, ng);
        }

      /* Check cross section */

      CheckValue(FUNCTION_NAME, "totxs", "", totxs, 0.0, INFTY);

      /* Value may be zero if energy is outside the boundaries */

      if (totxs == 0.0)
        return -2;

      /* Sample fraction of total cross section */

      totxs = RandF(id)*totxs;

      /* Get pointer to partial list */

      lst = (long)RDB[rea + REACTION_PTR_PARTIAL_LIST];
      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

      /* Avoid compiler warning */

      xs = -1.0;
      rls = -1;

      /* Reset reaction pointer (rewind list) */

      rea = -1;

      /* Loop over reactions */

      while ((rls = NextReaction(lst, &rea, &adens, &Emin, &Emax, id))
             > VALID_PTR)
        {
          /* Check reaction pointer and energy */

          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          CheckValue(FUNCTION_NAME, "E", "", E, Emin, Emax);

          /* Get microscopic cross section */

          xs = MGXS(rea, E, ng);

          /* Adjust total and check */

          if ((totxs = totxs - adens*xs) < 0.0)
            break;
        }

      /* Check pointers */

      if ((rea < VALID_PTR) || (rls < VALID_PTR))
        {
          /* Print warning */

          Warn(FUNCTION_NAME, "Failed to sample target nuclide (mat %s)",
               GetText(mat + MATERIAL_PTR_NAME));

          /* Score failure */

          ptr = (long)RDB[RES_REA_SAMPLING_FAIL];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, 1.0, ptr, id, type - 1);

          /* Reject sample */

          return -3;
        }

      /* Pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Microscopic total cross section used for sampling nuclide */

      totxs = xs;
      CheckValue(FUNCTION_NAME, "totxs", "", totxs, ZERO, INFTY);

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Continuous-energy mode ******************************************/

      /* Check particle type */

      if (type == PARTICLE_TYPE_NEUTRON)
        {
          /* Check TMS mode and get pointer to reaction channel */

          if (TMS == TMS_MODE_MG)
            Die(FUNCTION_NAME, "Invalid TMS mode");
          else if (TMS == TMS_MODE_CE)
            {
              /* Use continous-energy majorant */

              rea = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            }
          else
            {
              /* Use total neutron cross section */

              rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
            }

          /* Get total cross section */

          totxs = MacroXS(rea, E, id);

          /* Check implicit capture mode and get absorption xs*/

          if (RDB[DATA_OPT_IMPL_CAPT] == YES)
            {
              /* Check TMS */

              if (TMS != TMS_MODE_NONE)
                Die(FUNCTION_NAME, "TMS not allowed with implicit capture");

              /* Get total absorption cross section (may be zero for He) */

              if ((ptr = (long)RDB[mat + MATERIAL_PTR_ABSXS]) > VALID_PTR)
                absxs = MacroXS(ptr, E, id);
              else
                absxs = 0.0;

              /* Calculate total cross section reduced by absorption */

              totxs = totxs - absxs;
            }
        }
      else
        {
          /* Get total photon cross section */

          rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          totxs = PhotonMacroXS(rea, E, id);
        }

      /* Check cross section */

      CheckValue(FUNCTION_NAME, "totxs", "", totxs, 0.0, INFTY);

      /* Value may be zero if energy is outside the boundaries */

      if (totxs == 0.0)
        return -4;

      /* Sample fraction of total cross section */

      totxs = RandF(id)*totxs;

      /* Avoid compiler warning */

      xs = -1.0;
      rls = -1;

      /* Get pointer to partial list */

      lst = (long)RDB[rea + REACTION_PTR_PARTIAL_LIST];
      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

      /* Reset reaction pointer (rewind list) */

      rea = -1;

      /* Loop over reactions */

      while ((rls = NextReaction(lst, &rea, &adens, &Emin, &Emax, id))
             > VALID_PTR)
        {
          /* Check reaction pointer and energy */

          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          CheckValue(FUNCTION_NAME, "E", "", E, Emin, Emax);

          /* Check particle type */

          if (type == PARTICLE_TYPE_NEUTRON)
            {
              /* Get total neutron cross section */

              if (TMS == TMS_MODE_CE)
                xs = MicroMajorantXS(rea, E, id);
              else
                xs = MicroXS(rea, E, id);

              /* Check implicit capture mode and get absorption xs*/

              if (RDB[DATA_OPT_IMPL_CAPT] == YES)
                {
                  /* Pointer to nuclide */

                  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                  /* Get total absorption cross section */

                  ptr = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  absxs = MicroXS(ptr, E, id);

                  /* Calculate total cross section reduced by absorption */

                  xs = xs - absxs;
                }
            }
          else
            {
              /* Get total gamma cross section */

              xs = PhotonMicroXS(rea, E, id);
            }

          /* Adjust total and check */

          if ((totxs = totxs - adens*xs) < 0.0)
            break;
        }

      /* Check reaction pointer */

      if ((rea < VALID_PTR) || (rls < VALID_PTR))
        {
          /* Print warning */

          Warn(FUNCTION_NAME, "Failed to sample target nuclide (%s %ld %ld)",
               GetText(mat + MATERIAL_PTR_NAME), rea, rls);

          /* Score failure */

          ptr = (long)RDB[RES_REA_SAMPLING_FAIL];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, 1.0, ptr, id, type - 1);

          /* Return null */

          return -5;
        }

      /* Pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Microscopic total cross section used for sampling nuclide */

      totxs = xs;
      CheckValue(FUNCTION_NAME, "totxs", "", totxs, ZERO, INFTY);

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Handle rejections ***************************************************/

  /* Check nuclide, reaction and list pointers */

  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  CheckPointer(FUNCTION_NAME, "(rls)", DATA_ARRAY, rls);

  /* Check total cross section */

  CheckValue(FUNCTION_NAME, "totxs", "", totxs, ZERO, INFTY);

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      /* Check mode */

      if (TMS != TMS_MODE_NONE)
        {
          /* Get material temperature for temperature rejection */

          T = GetTemp(mat, id);
          CheckValue(FUNCTION_NAME, "T", "", T, ZERO, INFTY);

          /* Get total cross section */

          rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          xs = DopMicroXS(mat, rea, E, &Er, T, id);
          CheckValue(FUNCTION_NAME, "xs", "", xs, 0, INFTY);

          /* Low-energy correction factor */

          g = PotCorr(nuc, E, T*KELVIN);
          CheckValue(FUNCTION_NAME, "g", "", g, ZERO, INFTY);

          /* Score total number of TMS samples for efficiency */

          ptr = (long)RDB[RES_TMS_SAMPLING_EFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, 1.0, ptr, id, 1);

          /* Score total number of TMS samples for failure fraction */

          ptr = (long)RDB[RES_TMS_FAIL_STAT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, 1.0, ptr, id, 0);

          /* Calculate rejection probability */

          f = xs/totxs;

          /* Check majorant failure */

          if (f > 1.000001)
            {
              /* Score */

              AddBuf1D(1.0, 1.0, ptr, id, 1);

              /* Print warning */

#ifdef DEBUG

              Warn(FUNCTION_NAME,
                   "TMS majorant exceeded %s E: %E Er: %E T %f f: %f\n",
                   GetText(nuc + NUCLIDE_PTR_NAME), E, Er, T, f);
#endif
            }

          /* Rejection sampling */

          if (RandF(id) > f)
            return -6;

          /* Score total number of TMS sampling efficiency */

          ptr = (long)RDB[RES_TMS_SAMPLING_EFF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf1D(1.0, 1.0, ptr, id, 0);

          /* Put energy and cross section for reaction sampling */

          E0 = E;
          E = Er;
          totxs = xs/g;
        }
      else if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
        {
          /* Normal multi-group rejection, get total cross section */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          xs = MicroXS(ptr, E, id);

          /* Check cross sections */

          CheckValue(FUNCTION_NAME, "xs", "", xs, 0.0, INFTY);
          CheckValue(FUNCTION_NAME, "totxs", "", totxs, xs, INFTY);

          /* Calculate rejection probability */

          f = xs/totxs;
          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);

          /* Rejection sampling */

          if (RandF(id) > f)
            return -7;

          /* Put cross section for reaction sampling */

          totxs = xs;
        }
    }

  /* Add to nuclide counter */

  ptr = (long)RDB[rls + RLS_DATA_PTR_COUNT];
  CheckPointer(FUNCTION_NAME, "(ptr)",PRIVA_ARRAY, ptr);
  AddPrivateData(ptr, 1.0, id);

  /***************************************************************************/

  /***** Sample reaction *****************************************************/

  /* Check nuclide pointer and total cross section */

  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
  CheckValue(FUNCTION_NAME, "totxs", "", totxs, ZERO, INFTY);

  /* Sample fraction of microscopic total cross section */

  totxs = RandF(id)*totxs;

  /* Pointer to partial reaction list */

  lst = (long)RDB[nuc + NUCLIDE_PTR_SAMPLE_REA_LIST];
  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

  /* Loop over partials */

  i = 0;
  while ((rls = ListPtr(lst, i++)) > VALID_PTR)
    {
      /* Pointer to reaction data */

      rea = (long)RDB[rls + RLS_DATA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get cross section */

      if (type == PARTICLE_TYPE_NEUTRON)
        {
          /* On-the-fly S(a,b) reactions require a special treatment */

          if (((long)RDB[rea + REACTION_MT] == 2004) ||
              ((long)RDB[rea + REACTION_MT] == 2002))
            {

              /* Kontribuutio pitää jakaa g:llä koska myös sämpläyksessä */
              /* käytettävä kokonaisvaikutusala jaettiin sillä */

              xs = OTFSabXS(rea, E0, T, id)/g;
            }
          else if ((E0 < RDB[nuc + NUCLIDE_SAB_EMAX]) &&
                   ((long)RDB[rea + REACTION_MT] == 2))
            {
              /* Skip free atom elastic reaction if OTF S(a,b) treatment */
              /* is in use and energy is not between Emaxlow and Emax,   */
              /* the region in which only 1/2 S(a,b) reaction is active  */

              if((RDB[nuc + NUCLIDE_SAB_EMAXLOW] == 0.0) ||
                 (E0 < RDB[nuc + NUCLIDE_SAB_EMAXLOW]))
                continue;

              /* If energy is between Emaxlow and Emax, */
              /* xs = elastic - active S(a,b) */

              else
                xs = MicroXS(rea, E, id) - OTFSabXS(rea, E0, T, id)/g;
            }
          else
            xs = MicroXS(rea, E, id);
        }
      else
        xs = PhotonMicroXS(rea, E, id);

      /* Adjust total and check */

      if ((totxs = totxs - xs) < 0.0)
        {
          /* Add to counter */

          ptr = (long)RDB[rls + RLS_DATA_PTR_COUNT];
          CheckPointer(FUNCTION_NAME, "(ptr)",PRIVA_ARRAY, ptr);
          AddPrivateData(ptr, 1.0, id);

          /* Score analog reaction rate estimator */

          if ((ptr = (long)RDB[rea + REACTION_PTR_ANA_RATE]) > VALID_PTR)
            AddBuf1D(1.0, wgt, ptr, id, 0);

          /* Return reaction pointer */

          return rea;
        }
    }

  /* Check remaining total xs */

  if (totxs > 1E-9)
    {
      /* Print warning */

      Warn(FUNCTION_NAME,
          "Reaction sampling failed miserably (%s in %s E = %E totxs = %E)",
          GetText(nuc + NUCLIDE_PTR_NAME), GetText(mat + MATERIAL_PTR_NAME), E,
          totxs);

      /* Check reaction list sum */

      CheckReaListSum(mat, PARTICLE_TYPE_NEUTRON, E, YES, id);
    }
  else
    {
      /* Check reaction list sum */

      CheckReaListSum(mat, PARTICLE_TYPE_NEUTRON, E, NO, id);

      /* Print warning */

      Warn(FUNCTION_NAME, "Reaction sampling failed (%s in %s)",
           GetText(nuc + NUCLIDE_PTR_NAME), GetText(mat + MATERIAL_PTR_NAME));
    }

  /* Score failure */

  ptr = (long)RDB[RES_REA_SAMPLING_FAIL];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(1.0, 1.0, ptr, id, type - 1);

  /* Reject sample */

  return -8;

  /***************************************************************************/
}
