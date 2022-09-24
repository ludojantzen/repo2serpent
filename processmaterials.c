/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processmaterials.c                             */
/*                                                                           */
/* Created:       2010/12/28 (JLe)                                           */
/* Last modified: 2020/05/09 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Processes material compositions, sets nuclide and reaction   */
/*              lists.                                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessMaterials:"

/*****************************************************************************/

void ProcessMaterials()
{
  long mat, ptr, nmat, nbumat, i, TMS, iso, nuc;
  double T;

  /* Check burnup step */

  if (((long)RDB[DATA_BURN_STEP] > 0) &&
      ((long)RDB[DATA_PTR_HISV0]) < VALID_PTR)
    Die(FUNCTION_NAME, "Should not be here");

  /* Pointer to material list */

  if ((mat = (long)RDB[DATA_PTR_M0]) < VALID_PTR)
    Error(0, "No material definitions");

  /***************************************************************************/

  /***** Finalize material compositions **************************************/

  fprintf(outp, "Normalizing compositions and processing mixtures...\n");

  /* Calculate fractions */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      IsotopeFractions(mat);
      mat = NextItem(mat);
    }

  /* Process mixtures */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      ProcessMixture(mat, 0);
      mat = NextItem(mat);
    }

  fprintf(outp, "OK.\n\n");

  /* Add data for critical concentration calculation of user chosen nuclides */

  ProcessIterNucs();

  /* Print decomposed mixtures */

  PrintMixtures();

  /* Process transport corrections */

  ProcessTranspCorr();

  /* Replace isotopic with atomic data in photon transport calculation */
  /* (NOTE: tätä kutsutaan myöhemmin myös tuolla alempana). */

  ReplacePhotonData();

  /***************************************************************************/

  /***** Sort composition for better memory management ***********************/

  /* 13.3.2018 / 2.1.30 / JLe: Tämän sorttauksen idea on jo unohtunut, */
  /* mutta pidetään vielä mukana. */

  /* Check if macroscopic cross sections are pregenerated (if */
  /* number of materials is large, the calculation hangs here) */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES)
    {
      /* Use color index to remember initial order */

      i = 0;

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check divided and burn flags */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
            WDB[mat + MATERIAL_SORT_IDX] = 1E+12;
          else if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
            WDB[mat + MATERIAL_SORT_IDX] = 1E+11;
          else
            WDB[mat + MATERIAL_SORT_IDX] = (double)(i++);

          /* Next material */

          mat = NextItem(mat);
        }

      /* Use OpenMP thread number for sort */

      i = 0;

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check divided and burn flags */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
            WDB[mat + MATERIAL_OMP_ID] = 1E+12;
          else if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
            WDB[mat + MATERIAL_OMP_ID] = 1E+11;
          else
            {
              /* Set id */

              WDB[mat + MATERIAL_OMP_ID] = (double)(i++);

              /* Check id */

              if (i == (long)RDB[DATA_OMP_MAX_THREADS])
                i = 0;
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Sort */

      mat = (long)RDB[DATA_PTR_M0];
      SortList(mat, MATERIAL_OMP_ID, SORT_MODE_ASCEND);
    }

  /***************************************************************************/

  /***** Allocate memory and process *****************************************/

  /* Process compositions of burnable materials */

  BurnMatCompositions();

  /* Put composition to divided materials */

  PutCompositions();

  /* Process burnable materials */

  ProcessBurnMat();

  /* Calculate masses (to get the value on output) */

  CalculateMasses();

  /* Close composition lists */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if domain decomposition is in use */

      if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
          ((long)RDB[mat + MATERIAL_MPI_ID] > -1) &&
          ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Pointer to composition */

      ptr = (long)RDB[mat + MATERIAL_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      CloseList(ptr);

      /* Tämä on lista alkuperäiseen nuklidikoostumukseen joka korvataan */
      /* alkuainekohtaisella koostumuksella fotonitransportmoodissa jos  */
      /* inputti on annettu neutronidatana (JLe 3.6.2015 / 2.1.24). */

      if ((ptr = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) > VALID_PTR)
        CloseList(ptr);

      mat = NextItem(mat);
    }

  /* Check temperatures and TMS flags */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Skip if no neutron transport mode */

      if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == NO)
        break;

      /* Get temperature */

      T = RDB[mat + MATERIAL_DOPPLER_TEMP];

      /* Get TMS flag */

      if (RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE)
        TMS = NUCLIDE_FLAG_TMS;
      else
        TMS = 0;

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
              /* Next */

              iso = NextItem(iso);

              /* Cycle loop */

              continue;
            }

          /* Compare temperature */

          if ((T > 0) && (T != RDB[nuc + NUCLIDE_TEMP]))
            Die(FUNCTION_NAME, "Error in temperature: %s %s %E %E",
                GetText(mat + MATERIAL_PTR_NAME),
                GetText(nuc + NUCLIDE_PTR_NAME), T, RDB[nuc + NUCLIDE_TEMP]);

          /* Check TMS flag */

          if (TMS != ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS))
            Die(FUNCTION_NAME, "Error in TMS flag %s %s %ld %ld",
                GetText(mat + MATERIAL_PTR_NAME), GetText(nuc + NUCLIDE_PTR_NAME), TMS,
                ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS));

          /* Check TMS limits */

          if (TMS != NO)
            {
              /* Check minimum */

              if (RDB[nuc + NUCLIDE_TMS_MIN_TEMP] >
                  RDB[mat + MATERIAL_TMS_TMIN])
                Die(FUNCTION_NAME, "Error in TMS Tmin: %s %s %E %E",
                    GetText(mat + MATERIAL_PTR_NAME),
                    GetText(nuc + NUCLIDE_PTR_NAME),
                    RDB[nuc + NUCLIDE_TMS_MIN_TEMP],
                    RDB[mat + MATERIAL_TMS_TMIN]);

              /* Check maximum */

              if (RDB[nuc + NUCLIDE_TMS_MAX_TEMP] <
                  RDB[mat + MATERIAL_TMS_TMAX])
                Die(FUNCTION_NAME, "Error in TMS Tmax: %s %s %E %E",
                    GetText(mat + MATERIAL_PTR_NAME),
                    GetText(nuc + NUCLIDE_PTR_NAME),
                    RDB[nuc + NUCLIDE_TMS_MAX_TEMP],
                    RDB[mat + MATERIAL_TMS_TMAX]);
            }

          /* Next */

          iso = NextItem(iso);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Allocate memory for reaction lists and macroscopic cross sections */

  AllocMacroXS();

  /* Sort composition to get initial order */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES)
    {
      mat = (long)RDB[DATA_PTR_M0];
      SortList(mat, MATERIAL_SORT_IDX, SORT_MODE_ASCEND);
    }

  /* Re-read compositions from restart file (JLe: tämä pitää lukea uudestaan  */
  /* siitä syystä että alialuejako tehdään vasta PutCompositions():issa, jota */
  /* kutsutaan tuolla ylempänä). */

  ReadRestartFile(RESTART_OVERRIDE);

  /* Replace isotopic with atomic data in photon transport calculation */
  /* (JLe: Tätä joudutaan kutsumaan uudestaan että myös alialueiden    */
  /* koostumukset menee oikein). */

  ReplacePhotonData();

  /* Check if decay source is used or multi-group spectra requested */

  if (((long)RDB[DATA_USE_DECAY_SRC] == YES) ||
      ((long)RDB[DATA_GSPEC_PTR_EGRID] > VALID_PTR))
    {
      /* Check domain decomposition */

      if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
        Die(FUNCTION_NAME, "Domain decomposition in use");

      /* Calculate activities (Aktiivisuudet lasketaan uudestaan */
      /* transportcycle.c:n lopussa.) */

      CalculateActivities();

      /* Process decay source */

      ProcessDecaySrc();

      /* Print gamma source (spectra for testing) */

      PrintGammaSpectra();
    }

  /* Process photon data if photon transport mode */

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    {
      /* Check domain decomposition */

      if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
        Die(FUNCTION_NAME, "Domain decomposition in use");

      /* Process photon attenuation data */

      ProcessPhotonAtt();

      /* Process electron data */

      if ((long)RDB[DATA_PHOTON_USE_TTB] == YES)
        ProcessElectrons();

      /* Process atomic relaxation data */

      ProcessRelaxation();
    }

  /***************************************************************************/

  /***** Assign MPI numbers to materials *************************************/

  /* NOTE: tää pitää miettiä uudestaan niin että toi järjestys menee */
  /* jotenkin fiksusti */

  /* Set MPI id's */

  i = 0;

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check domain decomposition */

      if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
        {
          /* ID's should already have been set, check MSR mode */

          if ((long)RDB[mat + MATERIAL_FLOW_IDX] > 0)
            Error(mat,
                  "Domain decomposition not available with material flows");
        }
      else if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        {
          /* Set id (MSR-palamarutiinit pitää suorittaa yhdellä */
          /* MPI-taskilla) */

          if ((long)RDB[mat + MATERIAL_FLOW_IDX] == 0)
            WDB[mat + MATERIAL_MPI_ID] = (double)(i++);
          else
            WDB[mat + MATERIAL_MPI_ID] = 0.0;

          /* Check id */

          if (i == mpitasks)
            i = 0;
        }
      else
        WDB[mat + MATERIAL_MPI_ID] = -1.0;

      /* Reset OpenMP thread number */

      WDB[mat + MATERIAL_OMP_ID] = -1.0;

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Set pointers to first divided material ******************************/

  /* NOTE: Tuota käytetään vain calculatedtmajorants.c:ssä ettei */
  /* jouduta optimointimoodeissa 1 ja 2 luuppaamaan kaikkien yli */
  /* (atomitiheydet haetaan RLS_MAX_ADENS:eistä) */

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Pointer to parent */

      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
        {
          /* Get index */

          if ((i = (long)RDB[mat + MATERIAL_DIV_ZONE_IDX]) < 1)
            Die(FUNCTION_NAME, "Index not set");
          else if ((long)RDB[ptr + MATERIAL_DIV_PTR_FIRST] < VALID_PTR)
            {
              /* Check special conditions for domain decomposition */

              if (((long)RDB[DATA_DD_DECOMPOSE] == NO) ||
                  ((long)RDB[mat + MATERIAL_MPI_ID] < 0) ||
                  ((long)RDB[mat + MATERIAL_MPI_ID] == mpiid))
                WDB[ptr + MATERIAL_DIV_PTR_FIRST] = (double)mat;
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Count materials and set indexes for printout ************************/

  /* Reset counters */

  nmat = 0;
  nbumat = 0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Add total counter */

      nmat++;

      /* Put index */

      WDB[mat + MATERIAL_PROC_IDX] = (double)nmat;

      /* Check burn flag */

      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
          ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT))
        {
          /* Add burn counter */

          nbumat++;

          /* Put index */

          WDB[mat + MATERIAL_BURN_IDX] = (double)nmat;
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Put counters */

  WDB[DATA_N_MATERIALS] = (double)nmat;
  WDB[DATA_N_BURN_MATERIALS] = (double)nbumat;

  /***************************************************************************/

  /***** Calculate memory size and print summary *****************************/

  /* Calculated divided material total memory */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check pointer to parent */

      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
        WDB[ptr + MATERIAL_TOT_DIV_MEM_SIZE] =
          RDB[ptr + MATERIAL_TOT_DIV_MEM_SIZE] + RDB[mat + MATERIAL_MEM_SIZE];

      /* Next material */

      mat = NextItem(mat);
    }

  /* Set data structures for on-the-fly burnup solver */

  ProcessOTFBurn();

  /* Print material data */

  PrintMaterialData();

  /***************************************************************************/
}

/*****************************************************************************/
