/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printmaterialdata.c                            */
/*                                                                           */
/* Created:       2011/01/20 (JLe)                                           */
/* Last modified: 2020/06/09 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Prints material compositions to output                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintMaterialData:"

/*****************************************************************************/

void PrintMaterialData()
{
  long mat, iso, ptr, ptf, table, nuc, sab, mix, ndiv, nsub, ntot, mat0;
  double sum1, sum2, sum3, mem;
  char outfile[MAX_STR];
  FILE *fp;

  /*********************************************************************/

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Reset table counter (TODO: lue to DATA:sta) */

  table = 0;

  /* Set output file */

  sprintf(outfile, "%s.out", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file for writing */

  if ((fp = fopen(outfile, "a")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /***************************************************************************/

  /***** Material data *******************************************************/

  table++;

  fprintf(fp, "\n --- Table %2ld: Summary of material compositions: \n\n",
          table);

  ptr = (long)RDB[DATA_PTR_M0];
  fprintf(fp, " %ld materials included in calculation\n\n", ListSize(ptr));

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if material was produced by division */

      if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Get number of divisor regions */

      if ((ndiv = (long)RDB[mat + MATERIAL_DIV_N_ZONES]) < 1)
        ndiv = 1;

      if ((nsub = (long)RDB[mat + MATERIAL_DIV_N_SUB_ZONES]) < 1)
        nsub = 1;

      if ((ntot = (long)RDB[mat + MATERIAL_DIV_N_TOT_ZONES]) < 1)
        ntot = 1;

      /* Reset sums */

      sum1 = 0;
      sum2 = 0;
      sum3 = 0;

      /* Reset sab flag */

      sab = 0;

      /* Calculate sums */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Add to sums */

          sum1 = sum1 + RDB[iso + COMPOSITION_ADENS];

          if (RDB[mat + MATERIAL_ADENS] > 0.0)
            sum2 = sum2 +
              RDB[iso + COMPOSITION_ADENS]/RDB[mat + MATERIAL_ADENS];

          if (RDB[mat + MATERIAL_MDENS] > 0.0)
            sum3 = sum3 + RDB[nuc + NUCLIDE_AW]*RDB[iso + COMPOSITION_ADENS]
              /RDB[mat + MATERIAL_MDENS]/N_AVOGADRO;

          /* Check S(a,b) flag */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SAB_DATA)
            sab++;

          /* Next isotope */

          iso = NextItem(iso);
        }

      if (1 != 2)
        {
          fprintf(fp, "\nMaterial \"%s\":\n\n",
                  GetText(mat + MATERIAL_PTR_NAME));

          if (ntot > 1)
            {
              if ((ndiv == 1) || (nsub == 1))
                fprintf(fp, " - Material is divided into %ld depletion zones\n",
                      ntot);
              else
                fprintf(fp, " - Material is divided into %ld depletion zones (%ld cells and %ld sub-zones)\n",
                        ntot, ndiv, nsub);
            }

          if ((long)RDB[mat + MATERIAL_DT_MODE] == DT_MAT_BLOCK)
            fprintf(fp, " - Delta-tracking is never used inside material\n");
          else if ((long)RDB[mat + MATERIAL_DT_MODE] == DT_MAT_FORCE)
            fprintf(fp, " - Delta-tracking is always used inside material\n");

          /* Pointer to first divided */

          if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_FIRST]) < VALID_PTR)
            mat0 = mat;

          /* Use pointer to divided for these */

          if ((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            fprintf(fp, " - Material is burnable\n");
          else
            fprintf(fp, " - Material is not burnable\n");

          if ((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_INCLUDE_MAJORANT)
            fprintf(fp, " - Material is included in majorant\n");
          else
            fprintf(fp, " - Material is not included in majorant\n");

          if ((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT)
            fprintf(fp, " - Material is included in geometry\n");
          else
            fprintf(fp, " - Material is not included in geometry\n");

          if ((long)RDB[mat0 + MATERIAL_PTR_TRANSP_CORR] > VALID_PTR)
            fprintf(fp, " - Material is linked to transport correction data\n");

          fprintf(fp, " - Atom density %1.5E 1/barn*cm\n",
                  RDB[mat + MATERIAL_ADENS]);
          fprintf(fp, " - Mass density %1.5E g/cm3\n",
                  RDB[mat + MATERIAL_MDENS]);

          fprintf(fp, " - Volume %1.5E cm3\n", RDB[mat + MATERIAL_VOLUME]);

          fprintf(fp, " - Mass %1.5E g\n", RDB[mat + MATERIAL_MASS]);

          fprintf(fp, " - Photon emission rate %1.5E 1/s\n",
                  RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE]);

          fprintf(fp, " - Neutron emission rate %1.5E 1/s\n",
                  RDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE]);

          fprintf(fp, " - %ld nuclides in composition\n",
                  ListSize((long)RDB[mat + MATERIAL_PTR_COMP]));

          if (sab > 0)
            fprintf(fp, " - %ld nuclides associated with S(a,b) data\n", sab);
          else
            fprintf(fp, " - No nuclides associated with S(a,b) data\n");

          if ((mem = RDB[mat + MATERIAL_MEM_SIZE]/KILO) < KILO)
            fprintf(fp, " - %1.2f kb of memory allocated for data\n", mem);
          else
            fprintf(fp, " - %1.2f Mb of memory allocated for data\n", mem/KILO);

          if (ntot > 1)
            {
              if ((mem = RDB[mat + MATERIAL_TOT_DIV_MEM_SIZE]/KILO) < KILO)
                fprintf(fp, " - %1.2f kb of memory allocated for all zones\n",
                        mem);
              else
                fprintf(fp, " - %1.2f Mb of memory allocated for all zones\n",
                        mem/KILO);
            }

          /* Inflow list */

          if ((ptf = (long)RDB[mat + MATERIAL_PTR_INFLOW]) > VALID_PTR)
            {
              fprintf(fp, " - Inflow from materials:");

              /* Loop over list */

              while (ptf > VALID_PTR)
                {
                  /* Check pointer to self */

                  if (mat != (long)RDB[ptf + REPROC_CON_PTR_MAT2])
                    Die(FUNCTION_NAME, "No pointer to self");

                  /* Check other */

                  if ((ptr = (long)RDB[ptf + REPROC_CON_PTR_MAT1]) < VALID_PTR)
                    Die(FUNCTION_NAME, "Pointer error");

                  /* Print name */

                  fprintf(fp, " %s", GetText(ptr + MATERIAL_PTR_NAME));

                  /* Next */

                  ptf = NextItem(ptf);
                }

              fprintf(fp, "\n");
            }

          /* Outflow list */

          if ((ptf = (long)RDB[mat + MATERIAL_PTR_OUTFLOW]) > VALID_PTR)
            {
              fprintf(fp, " - Outflow to materials:");

              /* Loop over list */

              while (ptf > VALID_PTR)
                {
                  /* Check pointer to self */

                  if (mat != (long)RDB[ptf + REPROC_CON_PTR_MAT1])
                    Die(FUNCTION_NAME, "No pointer to self");

                  /* Check other */

                  if ((ptr = (long)RDB[ptf + REPROC_CON_PTR_MAT2]) < VALID_PTR)
                    Die(FUNCTION_NAME, "Pointer error");

                  /* Print name */

                  fprintf(fp, " %s", GetText(ptr + MATERIAL_PTR_NAME));

                  /* Next */

                  ptf = NextItem(ptf);
                }

              fprintf(fp, "\n");
            }

          /* Check mixture pointer */

          if ((mix = (long)RDB[mat + MATERIAL_PTR_MIX]) > VALID_PTR)
            {

              fprintf(fp, "\nMaterial is a mixture of %ld components:\n\n",
                      ListSize(mix));

              fprintf(fp, "-----------------------------------------\n");
              fprintf(fp, "Material             a. frac      m. frac\n");
              fprintf(fp, "-----------------------------------------\n");

              /* Loop over mixture */

              while (mix > VALID_PTR)
                {
                  /* Pointer to material */

                  ptr = (long)RDB[mix + MIXTURE_PTR_MAT];

                  fprintf(fp, "%-15s %12.5E %12.5E\n",
                          GetText(ptr + MATERIAL_PTR_NAME),
                          RDB[mix + MIXTURE_VFRAC],
                          RDB[mix + MIXTURE_MFRAC]);

                  /* Next material */

                  mix = NextItem(mix);
                }

              fprintf(fp, "-------------------------------------------------------------------\n");
            }

          fprintf(fp, "\nIsotopic composition (non-zero densities):\n\n");

          fprintf(fp, "-------------------------------------------------------------------\n");
          fprintf(fp, " Nuclide    a. weight   temp      a. dens      a. frac      m. frac\n");
          fprintf(fp, "-------------------------------------------------------------------\n");

          /* Loop over composition */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Check type */

              if (RDB[iso + COMPOSITION_ADENS] > 0.0)
                {
                  fprintf(fp, "%10s ", GetText(nuc + NUCLIDE_PTR_NAME));
                  fprintf(fp, "%10.5f ", RDB[nuc + NUCLIDE_AW]);

                  if (RDB[nuc + NUCLIDE_TEMP] < 0.0)
                    fprintf(fp, "   N/A ");
                  else
                    fprintf(fp, "%6.1f ", RDB[nuc + NUCLIDE_TEMP]);

                  fprintf(fp, "%12.5E ", RDB[iso + COMPOSITION_ADENS]);

                  if (RDB[mat + MATERIAL_ADENS] > 0.0)
                    fprintf(fp, "%12.5E ", RDB[iso + COMPOSITION_ADENS]
                            /RDB[mat + MATERIAL_ADENS]);
                  else
                    fprintf(fp, "%12.5E ", 0.0);

                  if (RDB[mat + MATERIAL_MDENS] > 0.0)
                    fprintf(fp, "%12.5E ", RDB[nuc + NUCLIDE_AW]
                            *RDB[iso + COMPOSITION_ADENS]
                            /RDB[mat + MATERIAL_MDENS]/N_AVOGADRO);
                  else
                    fprintf(fp, "%12.5E ", 0.0);

                  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                      NUCLIDE_FLAG_SAB_DATA)
                    fprintf(fp, "(s)");
                  else
                    fprintf(fp, "   ");

                  fprintf(fp, "\n");
                }

              /* Next isotope */

              iso = NextItem(iso);
            }

          fprintf(fp,"       sum                   %12.5E %12.5E %12.5E\n", sum1, sum2, sum3);


          fprintf(fp,"-------------------------------------------------------------------\n");

        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /* Close file */

  fclose(fp);
}

/*****************************************************************************/
