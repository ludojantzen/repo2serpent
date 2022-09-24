/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkrealistsum.c                              */
/*                                                                           */
/* Created:       2011/01/05 (JLe)                                           */
/* Last modified: 2019/11/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Check that reaction list sum matches total                   */
/*                                                                           */
/* Comments: - This is used to debug SampleReaction()                        */
/*           - Kutsu SampleMu:hun on väliaikainen ratkaisu                   */
/*           - Tästä tulee joskus errori kun E == 20 ja moodi on 0 tai 1     */
/*             (majorantti on nolla mutta reaktioiden summa ei)              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckReaListSum:"

/*****************************************************************************/

void CheckReaListSum(long mat, long type, double E, long kill, long id)
{
  long ptr, loc0, rea, n, nuc, lst1, rls1, iso, mt;
  double totxs, absxs, adens, sum, sum2, xs, Emin, Emax;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Check energy */

  CheckValue(FUNCTION_NAME, "E", "", E, ZERO, INFTY);

  if (kill == YES)
    {
      fprintf(errp, "\n\n%s Fatal error encountered...\n\n", FUNCTION_NAME);
      fprintf(errp, "mat = %s E = %E\n", GetText(mat + MATERIAL_PTR_NAME), E);
    }

  /***************************************************************************/

  /***** Nuclide totals ******************************************************/

  /* Loop over composition */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Skip lost */

      if (nuc == (long)RDB[DATA_PTR_NUCLIDE_LOST])
        {
          /* Pointer to next */

          iso = NextItem(iso);

          /* Cycle loop */

          continue;
        }

      /* Break at decay data */

      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
        break;

      /* Reset sum */

      sum = 0.0;
      sum2 = 0.0;

      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Skip all but partial reactions */

          if ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_PARTIAL)
            {
              /* Pointer to next */

              rea = NextItem(rea);

              /* Cycle loop */

              continue;
            }

          /* Get cross section and add to sum */

          xs = MicroXS(rea, E, id);
          sum = sum + xs;

          /* Get reaction mt */

          mt = (long)RDB[rea + REACTION_MT];

          /* Add total capture */

          if ((mt > 101) && (mt < 200))
            sum2 = sum2 + xs;

          /* Print */

          if (kill == YES)
            {
              fprintf(errp, "nuc = %10s mt = %4ld ",
                      GetText(nuc + NUCLIDE_PTR_NAME), mt);

              if ((E > RDB[rea + REACTION_URES_EMIN]) &&
                  (E < RDB[rea + REACTION_URES_EMAX]))
                fprintf(errp, "(u) ");
              else
                fprintf(errp, "    ");

              fprintf(errp, "Emin = %E ", RDB[rea + REACTION_EMIN]);
              fprintf(errp, "xs = %E sum = %E", xs, sum);
              /*
              if ((urs = (long)RDB[rea + REACTION_PTR_URES]) > VALID_PTR)
                fprintf(errp, "  %E %E",
                        TestValuePair(rea + REACTION_PTR_PREV_XS, E, id),
                        UresDiluMicroXS(rea, E, id));
              */
              fprintf(errp, "\n");
            }

          /* Next */

          rea = NextItem(rea);
        }

      /* Get total absorption cross section */

      rea = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      absxs = MicroXS(rea, E, id);

      /* Print */

      if (kill == YES)
        {
          fprintf(errp, "nuc = %10s mt = %4ld ",
                  GetText(nuc + NUCLIDE_PTR_NAME),
                  (long)RDB[rea + REACTION_MT]);
          if ((E > RDB[rea + REACTION_URES_EMIN]) &&
              (E < RDB[rea + REACTION_URES_EMAX]))
            fprintf(errp, "(u) ");
          else
            fprintf(errp, "    ");

          fprintf(errp, "Emin = %E ", RDB[rea + REACTION_EMIN]);
          fprintf(errp, "xs = %E sum = %E", absxs, sum2);
          /*
          if ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR)
            fprintf(errp, "  %E %E %E",
                    TestValuePair(rea + REACTION_PTR_PREV_XS, E, id),
                    UresDiluMicroXS(rea, E, id),
                    UresFactor(rea, E, id)
                    );
          */
          fprintf(errp, "\n");
        }

      /* Get total xs */

      rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      totxs = MicroXS(rea, E, id);

      /* Print */

      if (kill == YES)
        {
          fprintf(errp, "nuc = %10s mt = %4ld ",
                  GetText(nuc + NUCLIDE_PTR_NAME),
                  (long)RDB[rea + REACTION_MT]);
          if ((E > RDB[rea + REACTION_URES_EMIN]) &&
              (E < RDB[rea + REACTION_URES_EMAX]))
            fprintf(errp, "(u) ");
          else
            fprintf(errp, "    ");

          fprintf(errp, "Emin = %E ", RDB[rea + REACTION_EMIN]);
          fprintf(errp, "xs = %E sum = %E", totxs, sum);
          /*
          if ((urs = (long)RDB[rea + REACTION_PTR_URES]) > VALID_PTR)
            fprintf(errp, "  %E %E",
                    TestValuePair(rea + REACTION_PTR_PREV_XS, E, id),
                    UresDiluMicroXS(rea, E, id));
          */
          fprintf(errp, "\n");
        }

      /* Print or compare */

      if (kill == YES)
        fprintf(errp, "diff = %E %E\n\n", totxs - sum, absxs - sum2);
      else if ((fabs(sum - totxs) > 1E-10) && (fabs(sum/totxs - 1.0) > 1E-10))
        CheckReaListSum(mat, type, E, YES, id);

      /* Next */

      iso = NextItem(iso);
    }

  /***************************************************************************/

  /***** Check material total ************************************************/

  /* Get total absorption cross section */

  ptr = (long)RDB[mat + MATERIAL_PTR_ABSXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  absxs = MacroXS(ptr, E, id);

  if (kill == YES)
    fprintf(errp, "absxs    = %E (ures [%E %E])\n", absxs,
            RDB[ptr + REACTION_URES_EMIN], RDB[ptr + REACTION_URES_EMAX]);

  /* Get total cross section */

  ptr = (long)RDB[mat + MATERIAL_PTR_TOTXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  totxs = MacroXS(ptr, E, id);

  /* Check implicit capture mode */

  if (RDB[DATA_OPT_IMPL_CAPT] == YES)
    {
      /* Subtract absorption */

      totxs = totxs - absxs;
    }

  if (kill == YES)
    fprintf(errp, "totxs    = %E (ures [%E %E])\n", totxs,
            RDB[ptr + REACTION_URES_EMIN], RDB[ptr + REACTION_URES_EMAX]);

  /* Check type and get reaction pointer */

  if (type == PARTICLE_TYPE_NEUTRON)
    ptr = (long)RDB[DATA_PTR_MAJORANT];
  else if (type == PARTICLE_TYPE_GAMMA)
    ptr = (long)RDB[DATA_PTR_PHOTON_MAJORANT];
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Get majorant */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  xs = MajorantXS(ptr, E, id);

  if (kill == YES)
    fprintf(errp, "majorant = %E (ures [%E %E])\n\n", xs,
            RDB[ptr + REACTION_URES_EMIN], RDB[ptr + REACTION_URES_EMAX]);

  /* Reset sum */

  sum = 0.0;

  if ((loc0 = (long)RDB[mat + MATERIAL_PTR_TOT_REA_LIST]) > VALID_PTR)
    {
      if (kill == YES)
        fprintf(errp, "Calculating sum over nuclides and reactions:\n\n");

      /* Pointer to material total (toimii vaan neutroneille) */

      ptr = (long)RDB[mat + MATERIAL_PTR_TOTXS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to partial list */

      lst1 = (long)RDB[ptr + REACTION_PTR_PARTIAL_LIST];
      CheckPointer(FUNCTION_NAME, "(lst1)", DATA_ARRAY, lst1);

      /* Reset reaction pointer */

      rea = -1;

      /* Loop over reactions */

      while ((rls1 = NextReaction(lst1, &rea, &adens, &Emin, &Emax, id)) >
             VALID_PTR)
        {
          /* Check reaction pointer */

          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Pointer to nuclide */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Pointer to partial reaction list */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_SAMPLE_REA_LIST];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over partials */

          while (ptr > VALID_PTR)
            {
              /* Pointer to reaction data */

              rea = (long)RDB[ptr + RLS_DATA_PTR_REA];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Get cross section */

              xs = MicroXS(rea, E, id);

              /* Sample mu */

              if (xs > 0.0)
                for (n = 0; n < 2; n++)
                  SampleMu(rea, -1, E, NULL, NULL, id);

              /* Check energy */

              if (E < RDB[rea + REACTION_EMIN])
                {
                  if (xs > 0.0)
                    Die(FUNCTION_NAME,
                        "nonzero xs %E below boundary (%s %ld %E %E)",
                        xs, GetText(nuc + NUCLIDE_PTR_NAME),
                        (long)RDB[rea + REACTION_MT], E,
                        RDB[rea + REACTION_EMIN]);
                }
              else
                {
                  /* Add partial cross section */

                  sum = sum + adens*xs;

                  if (kill == YES)
                    {
                      fprintf(errp, "nuc = %10s adens = %E mt = %4ld ",
                              GetText(nuc + NUCLIDE_PTR_NAME), adens,
                              (long)RDB[rea + REACTION_MT]);

                      if ((E > RDB[rea + REACTION_URES_EMIN]) &&
                          (E < RDB[rea + REACTION_URES_EMAX]))
                        fprintf(errp, "(u) ");
                      else
                        fprintf(errp, "    ");

                      fprintf(errp, "Emin = %E ", RDB[rea + REACTION_EMIN]);
                      fprintf(errp, "xs = %E sum = %E\n", adens*xs, sum);
                    }
                }

              /* Next */

              ptr = NextItem(ptr);
            }
        }
    }
  else
    Die(FUNCTION_NAME, "Material %s has no partial or total list",
        GetText(mat + MATERIAL_PTR_NAME));

  /***************************************************************************/

  if (kill == YES)
    {
      fprintf(errp, "\nOptimization mode = %ld\n",
              (long)RDB[DATA_OPTI_MODE]);
      fprintf(errp, "Reconstruction tolerance = %E\n", RDB[DATA_ERG_TOL]);
      fprintf(errp, "totxs    = %E  sum = %E diff = %E\n", totxs, sum,
              sum/totxs - 1.0);

      /* Check type and get reaction pointer */

      if (type == PARTICLE_TYPE_NEUTRON)
        rea = (long)RDB[DATA_PTR_MAJORANT];
      else if (type == PARTICLE_TYPE_GAMMA)
        rea = (long)RDB[DATA_PTR_PHOTON_MAJORANT];
      else
        Die(FUNCTION_NAME, "Invalid particle type");

      if (rea > VALID_PTR)
        fprintf(errp, "majorant = %E diff = %E\n",
                MajorantXS(rea, E, id), totxs/MajorantXS(rea, E, id) - 1.0);

      Die(FUNCTION_NAME,
          "Error in total or majorant, or reaction sampling failed");
    }

  /* Tähän kosahtaa helposti opt = 3 -moodissa */

  /* Compare sum and total */

  if ((fabs(sum - totxs) > 1E-10) && (fabs(sum/totxs - 1.0) > 1E-10))
    CheckReaListSum(mat, type, E, YES, id);

  /* Check type and get reaction pointer */

  if (type == PARTICLE_TYPE_NEUTRON)
    rea = (long)RDB[DATA_PTR_MAJORANT];
  else if (type == PARTICLE_TYPE_GAMMA)
    rea = (long)RDB[DATA_PTR_PHOTON_MAJORANT];
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Get majorant cross section and check */

  if (rea > VALID_PTR)
    if (totxs/MajorantXS(rea, E, id) - 1.0 > 1E-6)
      CheckReaListSum(mat, type, E, YES, id);
}

/*****************************************************************************/
