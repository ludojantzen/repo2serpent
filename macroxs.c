/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : macroxs.c                                      */
/*                                                                           */
/* Created:       2011/01/02 (JLe)                                           */
/* Last modified: 2019/11/10 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Interpolates macroscopic cross section                       */
/*                                                                           */
/* Comments: - Toimiiko toi poison optio ollenkaan ilman pre-generated       */
/*             moodia?                                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MacroXS:"

/*****************************************************************************/

double MacroXS(long rea0, double E, long id)
{
  long i, ptr, rea, erg, ne, mat, nuc, ncol, mt;
  double xs0, xs1, xs, adens, f, mult, Emin, Emax, Er, T;

  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

  /* Get mt */

  mt = (long)RDB[rea0 + REACTION_MT];
  CheckValue(FUNCTION_NAME, "mt", "", mt, -57, -1);

  /* Get pointer to material */

  mat = (long)RDB[rea0 + REACTION_PTR_MAT];
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Get pointer to data */

  if ((ptr = (long)RDB[rea0 + REACTION_PTR_XS]) > VALID_PTR)
    {
      /***********************************************************************/

      /***** Interpolate pre-calculated data *********************************/

      /* Test existing data (NOTE: added 10.11.2019 / 2.1.32 / JLe) */

      if ((xs = TestValuePair(rea0 + REACTION_PTR_PREV_XS, E, id)) > -INFTY)
        return xs;

#ifdef DEBUG

      /* Sanity check for TMS */

      if (((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_MG) ||
          (((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_CE) &&
           (mt != MT_MACRO_TMP_MAJORANTXS)))
        Die(FUNCTION_NAME, "Pre-calculated data in %s %ld",
            GetText(mat + MATERIAL_PTR_NAME), mt);

#endif

      /* Test existing data (ei voi käyttää TMS:n kanssa) */
      /* ei voi myöskään käyttää data interfacen atomitiheyksien kanssa */

      if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_EXT_ADENS_MAT))
        if ((xs = TestValuePair(rea0 + REACTION_PTR_PREV_XS, E, id)) > -INFTY)
          return xs;

      /* Get pointer to energy grid */

      erg = (long)RDB[rea0 + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get interpolation factor */

      if ((f = GridFactor(erg, E, id)) < 0)
        xs = 0.0;
      else
        {
          /* Check interpolation factor */

          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, MAX_EGRID_NE);

          /* Separate integer and decimal parts of interpolation factor */

          i = (long)f;
          f = f - (double)i;

          /* Get number of points */

          ne = (long)RDB[rea0 + REACTION_XS_NE];
          CheckValue(FUNCTION_NAME, "ne", "", ne, 2, MAX_EGRID_NE);

          /* Check boundaries */

          if ((i < 0) || (i > ne - 1))
            xs = 0.0;
          else
            {
              /* Get tabulated cross sections */

              xs0 = RDB[ptr + i];
              xs1 = RDB[ptr + i + 1];

              if (mt != MT_MACRO_TMP_MAJORANTXS)
                {
                  /* Interpolate in normal case */

                  if (i == ne - 1)
                    xs = (1.0 - f)*xs0;
                  else
                    xs = f*(xs1 - xs0) + xs0;
                }
              else
                {
                  /* TMS-tapauksessa majoranttia ei interpoloida */
                  /* (histogrammimajorantti) */

                  xs = xs0;
                }
                  }
        }

      /* Add poison cross section */

      if (((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] == YES) ||
          ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] == YES))
        if ((mt == MT_MACRO_TOTXS) || (mt == MT_MACRO_ABSXS))
          xs = xs + PoisonXS(mat, E, mt, id);

      if ((long)RDB[DATA_ITER_MODE] == ITER_MODE_NUCLIDE)
        if ((mt == MT_MACRO_TOTXS) || (mt == MT_MACRO_ABSXS) ||
            (mt == MT_MACRO_TMP_MAJORANTXS))
          xs = xs + IterNucXS(mat, E, mt, -1, id);

      /* Add cross sections from data interface nuclides */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_EXT_ADENS_MAT)
        xs = xs + DataIFCXS(mat, E, mt, -1, id);

      /* Add cross sections for on-the-fly burnup mode */

      if ((long)RDB[DATA_OTF_BURN_MODE] == YES)
        if ((mt == MT_MACRO_TOTXS) ||
            (mt == MT_MACRO_ELAXS) ||
            (mt == MT_MACRO_ABSXS) ||
            (mt == MT_MACRO_FISSXS) ||
            (mt == MT_MACRO_FISSE) ||
            (mt == MT_MACRO_NSF) ||
            (mt == MT_MACRO_INLPRODXS) ||
            (mt == MT_MACRO_TMP_MAJORANTXS))
          xs = xs + OTFBurnXS(mat, E, mt, id);

      /* Perform ures correction */

      xs = MacroUresCorr(rea0, xs, E, id);

      /* Store value */

      StoreValuePair(rea0 + REACTION_PTR_PREV_XS, E, xs, id);

      /* Return value */

      return xs;

      /***********************************************************************/
    }
  else if (((long)RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE) &&
           (mt != MT_MACRO_TMP_MAJORANTXS))
    {
      /***********************************************************************/

      /***** TMS-moodi ******************************************************/

#ifdef DEBUG

      /* Sanity check for TMS */

      if ((RDB[mat + MATERIAL_TMS_TMIN] == 0.0) ||
          (RDB[mat + MATERIAL_TMS_TMIN] > RDB[mat + MATERIAL_TMS_TMAX]))
        Die(FUNCTION_NAME, "Error in temperature in %s",
            GetText(mat + MATERIAL_PTR_NAME));

#endif

      /* Get collision number */

      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      ncol = (long)GetPrivateData(ptr, id);

      /* Test existing data (HUOM! ncol, koska lämpötila voi olla eri) */

      if ((xs = TestValuePair(rea0 + REACTION_PTR_PREV_XS, (double)ncol, id))
          > -INFTY)
        return xs;

      /* Reset cross section */

      xs = 0.0;

      /* Get material temperature for on-the-fly temperature treatment */

      if ((T = GetTemp(mat, id)) > 0.0)
        {
          /* Get pointer to partial list */

          ptr = (long)RDB[rea0 + REACTION_PTR_PARTIAL_LIST];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Reset reaction pointer (rewind list) */

          rea = -1;

          /* Loop over reactions */

          while (NextReaction(ptr, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
            {
              /* Check reaction pointer */

              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* get multiplier */

              mult = ReaMulti(rea, mt, E, id);

              /* Pointer to nuclide */

              nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Add to cross section */

              xs = xs + mult*adens*DopMicroXS(mat, rea, E, &Er, T, id);

              /* Check energy cut-off */

              if (E < Emin)
                break;
            }

          /* Store cross section */

          StoreValuePair(rea0 + REACTION_PTR_PREV_XS, (double)ncol, xs, id);

          /* Return interpolated value */

          return xs;
        }
      /***********************************************************************/
    }
  else if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_EXT_ADENS_MAT)
    {
      /***********************************************************************/

      /************ Atomic densities given through interface *****************/
      /************ and summing up from partials             *****************/

      /* Get collision number */

      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      ncol = (long)GetPrivateData(ptr, id);

      /***** Calculate sum of partials *************************************/

      /* Test existing data (need to use ncol as atomic density can vary)  */

      if ((xs = TestValuePair(rea0 + REACTION_PTR_PREV_XS, (double)ncol, id))
          > -INFTY)
        return xs;

      /* Reset cross section */

      xs = 0.0;

      /* Get pointer to partial list */

      ptr = (long)RDB[rea0 + REACTION_PTR_PARTIAL_LIST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Reset reaction pointer (rewind list) */

      rea = -1;

      /* Loop over reactions */

      while (NextReaction(ptr, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
        {
          /* Check reaction pointer */

          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Get multiplier */

          mult = ReaMulti(rea, mt, E, id);

          /* Add to cross section */

          if (mt != MT_MACRO_TMP_MAJORANTXS)
            xs = xs + mult*adens*MicroXS(rea, E, id);

          /* In case of majorantxs, use MicroMajorantXS */

          else
            xs = xs + mult*adens*MicroMajorantXS(rea, E, id);

          /* Check energy cut-off */

          if (E < Emin)
            break;
        }

      /* Store cross section */

      StoreValuePair(rea0 + REACTION_PTR_PREV_XS, (double)ncol, xs, id);

      /* Return value */

      return xs;

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Calculate sum of partials *******************************************/

  /* Test existing data  */

  if ((xs = TestValuePair(rea0 + REACTION_PTR_PREV_XS, E, id)) > -INFTY)
      return xs;

  /* Reset cross section */

  xs = 0.0;

  /* Get pointer to partial list */

  ptr = (long)RDB[rea0 + REACTION_PTR_PARTIAL_LIST];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Reset reaction pointer (rewind list) */

  rea = -1;

  /* Loop over reactions */

  while (NextReaction(ptr, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
    {
      /* Check reaction pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get multiplier */

      mult = ReaMulti(rea, mt, E, id);

      /* Add to cross section */

      if (mt != MT_MACRO_TMP_MAJORANTXS)
        xs = xs + mult*adens*MicroXS(rea, E, id);

      /* In case of majorantxs, use MicroMajorantXS */

      else
        xs = xs + mult*adens*MicroMajorantXS(rea, E, id);

      /* Check energy cut-off */

      if (E < Emin)
        break;
    }

  /* Store cross section */

  StoreValuePair(rea0 + REACTION_PTR_PREV_XS, E, xs, id);

  /* Return value */

  return xs;

  /****************************************************************************/
}

/******************************************************************************/
