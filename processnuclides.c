/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processnuclides.c                              */
/*                                                                           */
/* Created:       2010/09/10 (JLe)                                           */
/* Last modified: 2020/04/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Selects which nuclides are included in the calculation,      */
/*              ace abd decat data, sets up decay / transmutation paths,     */
/*              constructs fission product lists, etc.                       */
/*                                                                           */
/* Comments: - TestDosFile pitää jotenkin ujuttaa opendatafileen             */
/*           - Termisen sironnan datasta luetaan pelkkä nimi (nuklidin vois  */
/*             lukea osuuteen?)                                              */
/*           - Toi OPT_BURN -juttu sotkee (optiota ei kopioida, ym)          */
/*           - External burnup moodi ei varmaan toimi                        */
/*                                                                           */
/*           - Muuta toi DATA_PTR_FP_LIB_ID_LIST paremmaksi                  */
/*                                                                           */
/*           - Noi aktinidi ja fp-listat pitää tehdä pelkille fissiileille   */
/*             materiaaleille, sillä muuten esim. eri lämpötilassa olevat    */
/*             absorbaattorit saa samat koostumukset                         */
/*                                                                           */
/*           - Täältä pitää siivota käyttämättömiä pointtereita, ym.         */
/*                                                                           */
/*           - Koko aliohjelman saisi yksinkertasemmaksi jos käsittelisi     */
/*             ensin alusta loppuun poltettavat materiaalit, ja sen jälkeen  */
/*             vasta lisäisi nuklidit muille. Tällöin ei tulisi esim. sitä   */
/*             ongelmaa että jos rakennemateriaalit on samassa lämpötilassa  */
/*             niiden vaikutusalalliset nuklidit muuttaa myös palavien       */
/*             materiaalien transmutaatioketjujen luomista.                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessNuclides:"

/*****************************************************************************/

void ProcessNuclides()
{
  long mat, src, iso, ace, nuc, rea, loc0, ptr, ZAI, n, TMS, warnflag;
  static char lib[MAX_STR], name[MAX_STR];
  double T, Tmin, Tmax;
  unsigned long seed;

  /* Get material pointer */

  if ((mat = (long)RDB[DATA_PTR_M0]) < VALID_PTR)
    Error(0, "No material definitions");

  /* Decompose elements into isotopes */

  DecomposeElements();

  /* Sort burnable materials at the beginning (this is needed to */
  /* get all processing done). (tää pois?) */

  SortList(mat, MATERIAL_BURN_SORT_FLAG, SORT_MODE_DESCEND);

  /* Read data from ACE directory files (pitää lukea myös silloin jos */
  /* decay-laskussa koostumukset on annettu vaikutusalamuodossa). */

  ReadDirectoryFile();

  /* Override nuclide id's for coefficient calculations */

  OverrideIDs();

  /* Init random number sequence */

  if (((long)RDB[DATA_BURN_RANDOMIZE_DEC] == YES) ||
      ((long)RDB[DATA_BURN_RANDOMIZE_FY] == YES))
    {
      seed = ReInitRNG(0);
      SEED[0] = seed;
    }

  /* Read decay data from ENDF format file */

  ReadDecayFile();

  /* Check that nuclear data exists */

  if ((long)RDB[DATA_PTR_ACE0] < 1)
    Error(0, "Missing file path for transport and/or decay data");

  /* Add stable nuclides at the end */

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES) ||
      (((long)RDB[DATA_USE_DECAY_SRC] == YES) &&
       ((long)RDB[DATA_READ_RESTART_FILE] == YES)))
    AddStableNuclides();

  /* Read fission yield data */

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES) ||
      ((long)RDB[DATA_OPTI_POISON_CALC] > -1) ||
      ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] == YES) ||
      ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] == YES))
    ReadFissionYields();

  /* Read isomeric branching ratios */

  ReadBRAFile();

  /* Set default branching ratios */

  DefaultBraData();

  /***************************************************************************/

  /***** Combine identical nuclides in composition ***************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check burn flags and neutron transport mode */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == NO)
          Error(mat,
                "Burnable materials allowed only in neutron transport mode");

      /* Check division */

      if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR)
        {
          /* Loop over composition */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Second loop */

              ptr = NextItem(iso);
              while (ptr > VALID_PTR)
                {
                  /* Compare names */

                  if (CompareStr(iso + COMPOSITION_PTR_NUCLIDE,
                                 ptr + COMPOSITION_PTR_NUCLIDE))
                    {
                      /* Combine densities */

                      WDB[iso + COMPOSITION_ADENS] =
                        RDB[iso + COMPOSITION_ADENS] +
                        RDB[ptr + COMPOSITION_ADENS];

                      /* Copy pointer */

                      n = ptr;

                      /* Pointer to next */

                      ptr = NextItem(ptr);

                      /* Remove duplicate */

                      RemoveItem(n);
                    }
                  else
                    {
                      /* Next */

                      ptr = NextItem(ptr);
                    }
                }

              /* Next */

              iso = NextItem(iso);
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Link nuclides to ACE data *******************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check division */

      if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR)
        {
          /* Loop over composition */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Get name */

              sprintf(name, "%s", GetText(iso + COMPOSITION_PTR_NUCLIDE));

              /* Try converting to ZAI (this will take care of decay nuclides */
              /* entered using symbolic names) */

              if ((ZAI = IsotoZAI(name)) > 0)
                sprintf(name, "%ld", ZAI);

              /* Loop over ace data (ei voi käyttää VALID_PTR) */

              ace = (long)RDB[DATA_PTR_ACE0];
              while (ace > 0)
                {
                  /* Compare */

                  WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_ALIAS];
                  if (!strcmp(GetText(DATA_DUMMY), name))
                    {
                      /* Put name */

                      WDB[iso + COMPOSITION_PTR_NUCLIDE] =
                        ACE[ace + ACE_PTR_NAME];

                      /* Break loop */

                      break;
                    }

                  /* next */

                  ace = (long)ACE[ace + ACE_PTR_NEXT];
                }

              /* Check if found */

              if (ace < 0)
                {
                  /* Check if elemental decomposition is on */

                  if ((long)RDB[DATA_ELEM_DECOMP] == YES)
                    {
                      /* Print warning */

                      Note(mat, "Nuclide %s not found in data libraries",
                           GetText(iso + COMPOSITION_PTR_NUCLIDE));

                      /* Copy pointer */

                      ptr = iso;
                      iso = NextItem(iso);

                      /* Remove */

                      RemoveItem(ptr);
                    }
                  else
                    Error(mat, "Nuclide %s not found in data libraries",
                          GetText(iso + COMPOSITION_PTR_NUCLIDE));
                }

              /* Check type */

              if ((long)ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_SAB)
                Error(mat, "S(a,b) library %s in composition list",
                      GetText(iso + COMPOSITION_PTR_NUCLIDE));

              /* Next */

              iso = NextItem(iso);
            }
        }

      /* Next */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /*** Read data to SAB-structure ********************************************/

  /* Loop over therm definitions */

  ptr = (long)RDB[DATA_PTR_T0];
  while (ptr > VALID_PTR)
    {
      /* Loop over S(a,b) data */

      loc0 = (long)RDB[ptr + THERM_PTR_SAB];
      while (loc0 > VALID_PTR)
        {
          /* Get name of first nuclide */

          sprintf(name, "%s", GetText(loc0 + SAB_PTR_NAME));

          /* Loop over ace data (ei voi käyttää VALID_PTR) */

          ace = (long)RDB[DATA_PTR_ACE0];
          while (ace > 0)
            {
              /* Compare */

              WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_ALIAS];
              if (!strcmp(GetText(DATA_DUMMY), name))
                {
                  /* Put name and temperature */

                  WDB[loc0 + SAB_PTR_NAME] = ACE[ace + ACE_PTR_NAME];
                  WDB[loc0 + SAB_T] = ACE[ace + ACE_TEMP];

                  /* Break loop */

                  break;
                }

              /* next */

              ace = (long)ACE[ace + ACE_PTR_NEXT];
            }

          /* Check if found */

          if (ace < 0)
            Error(ptr, "S(a,b) data %s not found in libraries",
                  GetText(loc0 + SAB_PTR_NAME));

          /* Check type */

          if ((long)ACE[ace + ACE_TYPE] != NUCLIDE_TYPE_SAB)
            Error(ptr, "%s is not an S(a,b) library",
                  GetText(loc0 + SAB_PTR_NAME));

          loc0=NextItem(loc0);
        }

      /* Next */

      ptr = NextItem(ptr);
    }

  /***************************************************************************/

  /***** Copy TMS flags from mixtures to materials ***************************/

  /* NOTE: Tää on nyt tehty siten että noi TMS parametrit kopioidaan  */
  /*       tässä mixturelta materiaaleille ja resetoidaan mixturessa. */
  /*       processmaterials.c:stä kutsuttu processmixture.c sitten    */
  /*       palauttaa ne parametrit materiaalille joka luodaan tuosta  */
  /*       mixturesta. (JLe 28.8.2015 / 2.1.25) */

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check for mixture and TMS flag */

      if ((loc0 = (long)RDB[mat + MATERIAL_PTR_MIX]) > VALID_PTR)
        if (RDB[mat + MATERIAL_TMS_MODE] == YES)
          {
            /* Loop over materials */

            while (loc0 > VALID_PTR)
              {
                /* Pointer to material */

                ptr = (long)RDB[loc0 + MIXTURE_PTR_MAT];
                CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                /* Check mixture of mixtures */

                if ((long)RDB[ptr + MATERIAL_PTR_MIX] > VALID_PTR)
                  Error(mat, "Mixtures of mixtures not allowed with TMS");

                /* Set TMS flag */

                WDB[ptr + MATERIAL_TMS_MODE] = (double)YES;

                /* Set temperatures */

                if (RDB[mat + MATERIAL_TMS_TMIN]
                    < RDB[ptr + MATERIAL_TMS_TMIN])
                  WDB[ptr + MATERIAL_TMS_TMIN] = RDB[mat + MATERIAL_TMS_TMIN];

                if (RDB[mat + MATERIAL_TMS_TMAX]
                    > RDB[ptr + MATERIAL_TMS_TMAX])
                  WDB[ptr + MATERIAL_TMS_TMAX] = RDB[mat + MATERIAL_TMS_TMAX];

                /* Next material */

                loc0 = NextItem(loc0);
              }

            /* Reset TMS flag and parameters for mixture */

            WDB[mat + MATERIAL_TMS_MODE] = (double)NO;
            WDB[mat + MATERIAL_TMS_TMIN] = INFTY;
            WDB[mat + MATERIAL_TMS_TMAX] = -INFTY;
          }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Process initial composition *****************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check division */

      if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR)
        {
          fprintf(outp, "Adding nuclides in material %s...\n\n",
                  GetText(mat + MATERIAL_PTR_NAME));

          /* Set temperature */

          T = RDB[mat + MATERIAL_DOPPLER_TEMP];

          /* Check TMS (tässä vielä NO tai YES, ton moodin voi asettaa */
          /* vasta SetOptimization:in jälkeen) */

          if (RDB[mat + MATERIAL_TMS_MODE] == YES)
            {
              /* Check that material is not a mixture */

              if ((long)RDB[mat + MATERIAL_PTR_MIX] > VALID_PTR)
                Error(mat, "Materials associated with TMS cannot be mixtures");

              /* Set flag */

              TMS = NUCLIDE_FLAG_TMS;

              /* Set mode to material */

              WDB[mat + MATERIAL_TMS_MODE] = RDB[DATA_TMS_MODE];
            }
          else
            TMS = 0;

          /* Reset counter */

          n = 0;

          /* Loop over composition */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Add nuclide */

              if ((nuc = AddNuclide(GetText(iso + COMPOSITION_PTR_NUCLIDE), -1,
                                    NULL, T, -1, TMS)) == 0)
                Die(FUNCTION_NAME, "Nuclide %s not found",
                    GetText(iso + COMPOSITION_PTR_NUCLIDE));
              else if (!((long)RDB[nuc + NUCLIDE_OPTIONS] & OPT_USED))
                n++;

              /* Set initial and used-flags */

              SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);
              SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

              /* Check burnup flag */

              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
                {
                  /* Set depletion flag for nuclide */

                  SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DEP);

                  /* Set daughter flag if not fixed composition */

                  if (((long)RDB[mat + MATERIAL_DEFAULT_PTR_LIB_ID] < VALID_PTR)
                      || ((long)RDB[DATA_BURN_DECAY_CALC] == YES))
                      SetOption(nuc + NUCLIDE_TYPE_FLAGS,
                              NUCLIDE_FLAG_NEW_DAUGHTERS);

                  /* Set default library ID and temp for decay nuclides */

                  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
                    {
                      /* Check that material ID is given */

                      if ((long)RDB[mat + MATERIAL_DEFAULT_PTR_LIB_ID]
                          < VALID_PTR)
                        Error(mat, "Default library ID must be given when decay nuclides are present in composition");
                      else
                        {
                          WDB[nuc + NUCLIDE_PTR_LIB_ID] =
                            RDB[mat + MATERIAL_DEFAULT_PTR_LIB_ID];
                          WDB[nuc + NUCLIDE_TEMP] =
                            RDB[mat + MATERIAL_DEFAULT_TMP];
                          WDB[nuc + NUCLIDE_XS_TEMP] =
                            RDB[mat + MATERIAL_DEFAULT_TMP];
                        }
                    }
                }

              /* Set pointer */

              WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)nuc;

              /* Next isotope */

              iso = NextItem(iso);
            }

          /* Print newline if nuclides were added */

          if (n > 0)
            fprintf(outp, "\n");
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Separate burnable absorber nuclides */

  SeparateBANuc();

  /***************************************************************************/

  /***** Add activation products *********************************************/

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES) &&
      ((long)RDB[DATA_BURN_DECAY_CALC] == NO))
    {
      /* Print */

      fprintf(outp, "Adding activation products...\n\n");

      /* Reset used flags */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Reset flag */

          ResetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

          /* Next */

          nuc = NextItem(nuc);
        }

      /* Loop over nuclides and add activation products */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Check initial flag */

          if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_INITIAL))
            break;

          /* Get TMS flag */

          TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

          /* Form paths */

          if ((long)RDB[nuc + NUCLIDE_ZAI] < 900000)
            if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP)
              FormTransmuPaths(nuc, 0, -1.0, -1.0, 1, TMS);

          /* Next */

          nuc = NextItem(nuc);
        }

      /* Print newline if nuclides were added */

      if (nuc > VALID_PTR)
        fprintf(outp, "\n");
    }

  /* Get library id's and temperatures for burnup calculation */

  GetBurnIDs();

  /***************************************************************************/

  /***** Add actinide data ***************************************************/

  /* Skip this in decay mode */

  if ((long)RDB[DATA_BURN_DECAY_CALC] == NO)
    {
      /* Combine actinides */

      CombineActinides();

      /* Loop over ID's */

      loc0 = (long)RDB[DATA_PTR_FP_LIB_ID_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Get library id, temperature and TMS flag */

          sprintf(lib, "%s", GetText(loc0 + FP_IDENT_PTR_ID));
          T = RDB[loc0 + FP_IDENT_TEMP];
          TMS = (long)RDB[loc0 + FP_IDENT_TMS];

          /* Get pointer to list */

          if ((ptr = (long)RDB[DATA_PTR_AC_ZAI_LIST]) > VALID_PTR)
            {
              /* Print */

              if (PrevItem(loc0) < VALID_PTR)
                fprintf(outp, "Adding actinide data...\n\n");

              /* Loop over list */

              while ((ZAI = (long)RDB[ptr++]) > 0)
                {
              /* Add nuclide */

                  if ((nuc = AddNuclide(NULL, ZAI, lib, T,
                                        NUCLIDE_TYPE_TRANSPORT, TMS))
                      > VALID_PTR)
                    {
                      /* Set activation product flag */

                      SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DEP);
                      SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_AP);
                      SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);
                    }
                }

              /* Print newline */

              if (NextItem(loc0) < VALID_PTR)
                fprintf(outp, "\n");
            }

          /* Next ID and temperature */

          loc0 = NextItem(loc0);
        }
    }

  /***************************************************************************/

  /***** Add fission yield data **********************************************/

  /* Skip this in decay mode */

  if ((long)RDB[DATA_BURN_DECAY_CALC] == NO)
    {
      /* Combine fission yields */

      CombineFissionYields();

      /* Loop over ID's */

      loc0 = (long)RDB[DATA_PTR_FP_LIB_ID_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Get library id, temperature TMS flag */

          sprintf(lib, "%s", GetText(loc0 + FP_IDENT_PTR_ID));
          T = RDB[loc0 + FP_IDENT_TEMP];
          TMS = (long)RDB[loc0 + FP_IDENT_TMS];

          /* Get pointer to list */

          if ((ptr = (long)RDB[DATA_PTR_FP_ZAI_LIST]) > VALID_PTR)
            {
              /* Print */

              if (PrevItem(loc0) < VALID_PTR)
                fprintf(outp, "Adding fission product data...\n\n");

              /* Loop over list */

              while ((ZAI = (long)RDB[ptr++]) > 0)
                {
                  /* Add nuclide */

                  if ((nuc = AddNuclide(NULL, ZAI, lib, T,
                                        NUCLIDE_TYPE_TRANSPORT, TMS))
                      > VALID_PTR)
                    {
                      /* Set activation product flag */

                      SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DEP);
                      SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_FP);
                      SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);
                    }
                }

              /* Print newline */

              if (NextItem(loc0) < VALID_PTR)
                fprintf(outp, "\n");
            }

          /* Next ID and temperature */

          loc0 = NextItem(loc0);
        }
    }

  /***************************************************************************/

  /***** Add remaining nuclides **********************************************/

  /* Check burnup mode */

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
    {
      /* Print */

      fprintf(outp, "Adding remaining data and forming paths...\n\n");

      /* Generate nuclide for lost data */

      nuc = NewItem(DATA_PTR_NUCLIDE_LOST, NUCLIDE_BLOCK_SIZE);

      /* Put name and alias and library ID */

      WDB[nuc + NUCLIDE_PTR_NAME] = (double)PutText("lost");
      WDB[nuc + NUCLIDE_PTR_LIB_ID] = RDB[nuc + NUCLIDE_PTR_NAME];

      /* Put ZAI */

      WDB[nuc + NUCLIDE_ZAI] = -1.0;

      /* Allocate memory for matrix index */

      AllocValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX);

      /* Put pointer */

      WDB[DATA_PTR_NUCLIDE_LOST] = (double)nuc;

      /* Reset used flags */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Reset flag */

          ResetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

          /* Next */

          nuc = NextItem(nuc);
        }

      /* Loop over initial nuclides and set paths */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Check initial flag */

          if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_INITIAL))
            break;

          /* Get TMS flag */

          TMS = ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS);

          /* Form paths */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DEP)
            FormTransmuPaths(nuc, 0, -1.0, -1.0, 2, TMS);

          /* Next */

          nuc = NextItem(nuc);
        }

      /* Print */

      fprintf(outp, "OK.\n\n");
    }

  /* Set path levels and used flags for initial nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check initial flag */

      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_INITIAL)
        {
          /* Set used flag */

          SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

          /* Set paht levels */

          SetPathLevels(nuc, 0);
        }

      /* Next */

      nuc = NextItem(nuc);
    }

  /* Remove unused nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  RemoveFlaggedItems(nuc, NUCLIDE_OPTIONS, OPT_USED, NO);

  /***************************************************************************/

  /***** Add nuclides in sources *********************************************/

  /* Loop over sources */

  n = 0;

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Check reaction pointer */

      if ((long)RDB[src + SRC_PTR_XSDATA] > VALID_PTR)
        {
          /* Check number and print */

          if (n++ == 0)
            fprintf(outp, "Adding nuclides for sources...\n\n");

          /* Add nuclide */

          if ((nuc = AddNuclide(GetText(src + SRC_PTR_XSDATA),-1, NULL,
                                -1.0, -1, NO)) == 0)
            Error(src, "Nuclide %s not found", GetText(src + SRC_PTR_XSDATA));

          /* Set source and used-flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_SRC);
          SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

          /* Set pointer */

          WDB[src + SRC_PTR_XSDATA] = (double)nuc;
        }

      /* Next source */

      src = NextItem(src);
    }

  /* Newline */

  if (n > 0)
    fprintf(outp, "\n");

  /***************************************************************************/

  /***** Set maximum temperatures for TMS ***********************************/

  /* Check TMS mode */

  if ((long)RDB[DATA_TMS_MODE] != TMS_MODE_NONE)
    {
      fprintf(outp, "Setting majorant temperatures for TMS:\n\n");

      /* Reset count */

      n = 0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check flag and divisor */

          if ((RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE) &&
              ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR))
            {
              /* Get temperatures */

              Tmin = RDB[mat + MATERIAL_TMS_TMIN];
              Tmax = RDB[mat + MATERIAL_TMS_TMAX];

              /* Check */

              if (Tmin > Tmax)
                Die(FUNCTION_NAME, "TMS minimum %E in %s exceeds maximum %E",
                    Tmin, GetText(mat + MATERIAL_PTR_NAME), Tmax);
              else if ((Tmin == 0.0) && (Tmax == 0.0))
                Die(FUNCTION_NAME, "Zero TMS temperatures");

              /* Loop over paths */

              if (Tmin < Tmax)
                fprintf(outp, "Material %s between %1.1fK and %1.1fK\n",
                        GetText(mat + MATERIAL_PTR_NAME), Tmin, Tmax);
              else
                fprintf(outp, "Material %s at %1.1fK\n",
                        GetText(mat + MATERIAL_PTR_NAME), Tmin);

              /* Add counter */

              n++;

              /* Reset used flags */

              nuc = (long)RDB[DATA_PTR_NUC0];
              while (nuc > VALID_PTR)
                {
                  /* Reset flag */

                  ResetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

                  /* Next */

                  nuc = NextItem(nuc);
                }

              /* Loop over composition */

              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              while (iso > VALID_PTR)
                {
                  /* Pointer to nuclide */

                  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                  /* Check TMS flag */

                  if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                        NUCLIDE_FLAG_TMS))
                    Die(FUNCTION_NAME, "TMS flag is not set");

                  /* Loop over chain */

                  FormTransmuPaths(nuc, 0, Tmin, Tmax, 3, NUCLIDE_FLAG_TMS);

                  /* Next */

                  iso = NextItem(iso);
                }
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check if there are no TMS nuclides and remove flag. NOTE: tän   */
      /* muuttamisen vaikutuksia ei ole testattu, voi aiheuttaa ongelmia */
      /* (5.2.2014 / 2.1.19) */

      if (n == 0)
        {
          Warn(FUNCTION_NAME,
               "TMS treatment not required for any material. TMS disabled.");
          WDB[DATA_TMS_MODE] = (double)TMS_MODE_NONE;
        }

      /* Final loop */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
        {
          /* Check flag */

          if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS))
            {
              /* Pointer to next */

              nuc = NextItem(nuc);

              /* Cycle loop */

              continue;
            }

          /* Check if maximum value was not set */

          if (RDB[nuc + NUCLIDE_TMS_MAX_TEMP] < 0.0)
            {
              Die(FUNCTION_NAME, "Shouldn't be here at all");

              /* Check that minimum value was not set either */

              if (RDB[nuc + NUCLIDE_TMS_MIN_TEMP] != INFTY)
                Die(FUNCTION_NAME, "Minimum TMS temperature was set");

              /* Reset temperatures */

              WDB[nuc + NUCLIDE_TMS_MAX_TEMP] = RDB[nuc + NUCLIDE_TEMP];
              WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = RDB[nuc + NUCLIDE_TEMP];
            }

          /* Sanity checks */

          if (RDB[nuc + NUCLIDE_TMS_MIN_TEMP] >
              RDB[nuc + NUCLIDE_TMS_MAX_TEMP])
            Die(FUNCTION_NAME, "TMS minimum exceeds maximum");
          else if (RDB[nuc + NUCLIDE_TMS_MIN_TEMP] < RDB[nuc + NUCLIDE_TEMP])
            Error(0, "TMS temperature %1.1fK below nuclide %s minimum",
                  RDB[nuc + NUCLIDE_TMS_MIN_TEMP],
                  GetText(nuc + NUCLIDE_PTR_NAME));
          else if (RDB[nuc + NUCLIDE_TEMP] < RDB[nuc + NUCLIDE_XS_TEMP])
            Die(FUNCTION_NAME, "Doppler temperature below XS temperature");

          /* Set temperature for Doppler-broadening */

          if (RDB[nuc + NUCLIDE_TMS_MIN_TEMP] > RDB[nuc + NUCLIDE_TEMP])
            {
              /* Check type */

              if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
                {
                  /* Set temperature (this is needed to trick later checks) */

                  if (RDB[nuc + NUCLIDE_TEMP] > 0.0)
                    WDB[nuc + NUCLIDE_TEMP] = RDB[nuc + NUCLIDE_TMS_MIN_TEMP];
                }
              else
                {
                  /* Check that nuclide is not already set for broadening */

                  if (RDB[nuc + NUCLIDE_TEMP] != RDB[nuc + NUCLIDE_XS_TEMP])
                    Die(FUNCTION_NAME, "Error in Doppler temperature");
                  else
                    WDB[nuc + NUCLIDE_TEMP] = RDB[nuc + NUCLIDE_TMS_MIN_TEMP];

                  /* Override if 0K data */

                  if (RDB[nuc + NUCLIDE_XS_TEMP] == 0.0)
                    {
                      WDB[nuc + NUCLIDE_TEMP] = 0.0;
                      WDB[nuc + NUCLIDE_TMS_MIN_TEMP] = 0.0;
                    }

                  /* Set preprocessor flag */

                  if (RDB[nuc + NUCLIDE_TEMP] != RDB[nuc + NUCLIDE_XS_TEMP])
                    WDB[DATA_USE_DOPPLER_PREPROCESSOR] = (double)YES;
                }
            }

          /* Next */

          nuc = NextItem(nuc);
        }

      /* Print newline */

      fprintf(outp, "\n");
    }

  /***************************************************************************/

  /***** Add zero-kelvin data for DBRC ***************************************/

  /* Check pointer to list. NOTE: Tää pitää tosiaan kutsua vasta täällä */
  /* ihan viimeiseksi tai muuten kaikkia reaktioita ei lueta mukaan jos */
  /* samaa nuklidia käytetään muuallakin. */

  if (((ptr = (long)RDB[DATA_PTR_DBRC]) > VALID_PTR) &&
      ((long)RDB[DATA_USE_DBRC] == YES))
    {
      /* Print */

      fprintf(outp, "Adding 0K data for DBRC...\n\n");

      /* Reset counter */

      /* Loop over DBRC nuclides */

      while ((long)RDB[ptr] > VALID_PTR)
        {
          /* Add nuclide */

          if ((nuc = AddNuclide(GetText(ptr), -1, NULL, -1,
                                NUCLIDE_TYPE_DBRC, NO)) == 0)
            Error(0, "Nuclide %s not found", GetText(ptr));

          /* Check temperature */

          if (RDB[nuc + NUCLIDE_TEMP] != 0.0)
            Error(0, "DBRC nuclide %s is not at 0K temperature",
                  GetText(nuc + NUCLIDE_PTR_NAME));

          /* Set used flag */

          SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

          /* Reset ures-flag */

          ResetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_URES_USED);

          /* Next nuclide in list */

          ptr++;
        }

      fprintf(outp, "\n");
    }

  /***************************************************************************/

  /***** Set fissile flags for materials *************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Set fissile material flag */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
            {
              /* Set flag */

              SetOption(mat + MATERIAL_OPTIONS, OPT_FISSILE_MAT);

              /* Break loop */

              break;
            }

          /* Next */

          iso = NextItem(iso);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Add data for equilibrium xenon calculation and poison production */

  ProcessPoisons();

  /* Doppler-broadening on cross sections (siirrettiin */
  /* tänne 21.8.2012 / 2.1.8) */

  DopplerBroad();

  /***************************************************************************/

  /***** Add S(a,b) data *****************************************************/

  /* Link thermal scattering data (must be called before the pointer */
  /* check to remove unused data) */

  LinkSabData();

  /* Check if data is defined */

  if ((ptr = (long)RDB[DATA_PTR_T0]) > VALID_PTR)
    {
      fprintf(outp, "Adding thermal scattering data...\n\n");

      /* Loop definitions */

      while (ptr > VALID_PTR)
        {
          /* Loop over S(a,b) data */

          loc0 = (long)RDB[ptr + THERM_PTR_SAB];

          while (loc0 > VALID_PTR)
            {
              /* Add nuclide (tähänkin suora pointteri?) */

              if ((nuc = AddNuclide(GetText(loc0 + SAB_PTR_NAME),
                                    -1, NULL, -1, -1, NO)) == 0)
                Die(FUNCTION_NAME, "S(a,b) library %s not found",
                    GetText(loc0 + SAB_PTR_NAME));

              /* Set pointer */

              WDB[loc0 + SAB_PTR_ISO] = (double)nuc;

              /* Set option */

              SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_INITIAL);

              /* Next S(a,b) data */

              loc0 = NextItem(loc0);
            }

          /* Check for makxsf temperature interpolation */

          if ((long)RDB[ptr + THERM_INTERP_MODE] == THERM_INTERP_MAKXSF)
            {
              /* Make a new interpolated THERM nuclide */

              if (( nuc = InterpolateSab(ptr) ) < VALID_PTR)
                Die(FUNCTION_NAME, "Interpolated S(a,b) nuclide not found");

              /* Update nuclide in first item */

              loc0 = (long)RDB[ptr + THERM_PTR_SAB];
              WDB[loc0 + SAB_PTR_ISO] = (double)nuc;

              /* Remove second item from list */

              loc0 = NextItem(loc0);
              RemoveItem(loc0);
            }

          /* Next therm definition */

          ptr = NextItem(ptr);
        }

      /* Add S(a,b) data into nuclides */

      AddSabData();

      /* Print newline */

      fprintf(outp, "\n");
    }

  /***************************************************************************/

  /***** Check for natural elements in burnable materials ********************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check burn flag */

      if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check A */

          if ((long)RDB[nuc + NUCLIDE_A] == 0)
            Note(mat, "Elemental XS (%s) in burnable material %s",
                 GetText(nuc + NUCLIDE_PTR_NAME),
                 GetText(mat + MATERIAL_PTR_NAME));

          /* Next */

          iso = NextItem(iso);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /****************************************************************************/

  /***** Remaining stuff *****************************************************/

  /* Set used flags */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Reset flag */

      SetOption(nuc + NUCLIDE_OPTIONS, OPT_USED);

      /* Next */

      nuc = NextItem(nuc);
    }

  /* Read additional photon data */

  ReadPhotonData();

  /* Close list */

  if ((nuc = (long)RDB[DATA_PTR_NUC0]) > VALID_PTR)
    CloseList(nuc);
  else
    Error(0, "No nuclide data read");

  /* Put number of nuclides */

  WDB[DATA_TOT_NUCLIDES] = (double)ListSize(nuc);

  /* Update reaction count */

  ReactionCount();

  /* Set global precursor group structure */

  SetPrecursorGroups();

  /* Set fission energies */

  SetFissE();

  /* Check and print */

  CheckNuclideData();

  /* Put indexes (sorting is done in CheckNuclideData(), lost is zero) */

  n = 1;
  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      WDB[nuc + NUCLIDE_IDX] = (double)(n++);
      WDB[nuc + NUCLIDE_INVENTORY_IDX] = -1.0;

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Adjust energy grid minimum and maximum */

  if (RDB[DATA_NEUTRON_EMIN] < RDB[DATA_NEUTRON_XS_EMIN])
    WDB[DATA_NEUTRON_EMIN] = RDB[DATA_NEUTRON_XS_EMIN];

  if (RDB[DATA_NEUTRON_EMAX] > RDB[DATA_NEUTRON_XS_EMAX])
    WDB[DATA_NEUTRON_EMAX] = RDB[DATA_NEUTRON_XS_EMAX];

  if (RDB[DATA_PHOTON_EMIN] < RDB[DATA_PHOTON_XS_EMIN])
    {
      /* Check if warning should be printed */

      warnflag = 0;
      if (fabs(RDB[DATA_PHOTON_EMIN]/RDB[DATA_PHOTON_XS_EMIN] - 1.0) > 1.0e-8)
        warnflag = 1;

      /* Tää haara antaa varoituksen aina jos fotonitransportmoodi  */
      /* ei oo päällä. Koko warnflag-hässäkästä olis muutenkin hyvä */
      /* päästä eroon (JLe/3.5.2018/2.1.31) */

      if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == NO)
        warnflag = 0;

      /* Print warning */

      if (warnflag)
        Warn(FUNCTION_NAME,
             "Minimum photon cross section energy %E MeV is\n"
             "above the energy grid minimum %E MeV.\n"
             "The energy grid minimum is set to %E MeV.\n"
             "Possible changes in energy cutoff cards (warned if any).",
             RDB[DATA_PHOTON_XS_EMIN], RDB[DATA_PHOTON_EMIN],
             RDB[DATA_PHOTON_XS_EMIN]);

      /* Update the minimum */

      WDB[DATA_PHOTON_EMIN] = RDB[DATA_PHOTON_XS_EMIN];

      /* Update the energy cutoff if needed */

      if (RDB[DATA_PHOTON_ECUT] < RDB[DATA_PHOTON_XS_EMIN])
        {
          /* Print warning */

          if (warnflag)
            Warn(FUNCTION_NAME, "Photon energy cutoff %E MeV is changed to %E.\n",
                 RDB[DATA_PHOTON_ECUT], RDB[DATA_PHOTON_XS_EMIN]);

          WDB[DATA_PHOTON_ECUT] = RDB[DATA_PHOTON_XS_EMIN];
        }
    }

  if (RDB[DATA_PHOTON_EMAX] > RDB[DATA_PHOTON_XS_EMAX])
    WDB[DATA_PHOTON_EMAX] = RDB[DATA_PHOTON_XS_EMAX];

  /* Link parent nuclide pointers if different (tämä siksi että S(a,b) */
  /* reaktioissa pointteri on sinne termisen sironnan nuklidiin, mikä  */
  /* sotkee esim. fotonituottoreaktiot */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Check pointer and set */

          if ((long)RDB[rea + REACTION_PTR_NUCLIDE] != nuc)
            WDB[rea + REACTION_PTR_PARENT_NUCLIDE] = (double)nuc;

          /* Next */

          rea = NextItem(rea);
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Read covariance data for nuclides that are included in the simulation */

  if ((ptr = (long)RDB[DATA_PTR_COVERXDATA_FNAME_LIST]) > VALID_PTR)
    while ((long)RDB[ptr] > VALID_PTR)
      {
        /* Read current file */

        ReadCOVERXFile(ptr);

        /* Next file */

        ptr++;
      }

  /***************************************************************************/
}

/*****************************************************************************/
