/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processprecdet.c                               */
/*                                                                           */
/* Created:       2015/05/07 (VVa)                                           */
/* Last modified: 2018/06/05 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sets up stats for implicit treatment of delayed neutron      */
/*              precursor groups                                             */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessPrecDet:"

/*****************************************************************************/


void ProcessPrecDet()
{
  long loc0, loc1, nuc, rea, lst, i, n, ptr, nprec;
  long nx, ny, nz, nt, tmplong;
  long *realist, *preclist;
  long reaarray, precarray, lamarray;
  char tmpstr[MAX_STR];
  double *lamlist, tmpdouble;
  FILE *fp;

  /* Get pointer to precursor detector or return */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

#ifdef DNPRINT
  fprintf(outp, "processprecdet.c-->\n");
#endif

  /***********************************************************************/
  /*********Calculate number of unique precursor groups*******************/
  /***********************************************************************/

  /* Reset number of unique precursors */

  nprec = 0;

  /* Reset temporary lists */

  lamlist  = (double *)Mem(MEM_ALLOC, 10000, sizeof(double));
  realist  = (long *)Mem(MEM_ALLOC, 10000, sizeof(long));
  preclist = (long *)Mem(MEM_ALLOC, 10000, sizeof(long));

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];

  while (nuc > VALID_PTR)
    {
      /* Skip nuclide if no ACE data */

      if ((long)RDB[nuc + NUCLIDE_PTR_ACE] < VALID_PTR)
        {
          /* Pointer to next nuclide */

          nuc = NextItem(nuc);

          /* Loop */

          continue;
        }

      /* Skip nuclide if not transport nuclide */

      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_TRANSPORT)
        {
          /* Pointer to next nuclide */

          nuc = NextItem(nuc);

          /* Loop */

          continue;

        }

      /* Pointer to reaction data */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Loop over reactions */

      while (rea > VALID_PTR)
        {

          /*
          printf("Nuclide %s, reaction %ld\n",
                 GetText(nuc + NUCLIDE_PTR_NAME),
                 (long)RDB[rea + REACTION_MT]);
          */

          /* Pointer to precursor group data */

          if ((lst = (long)RDB[rea + REACTION_PTR_PREC_LIST]) < VALID_PTR)
            {
              /* No precursor group list */

              /* Get next reaction */

              rea = NextItem(rea);

              /* Cycle loop */

              continue;
            }

          /* Loop over groups */

          n = 0;

          while ((loc1 = ListPtr(lst, n++)) > VALID_PTR)
            {

#ifdef DNPRINT
              fprintf(outp,
                "Nuclide %s, reaction %ld, prec ptr %ld, prec group %ld, prec lambda %E\n",
                     GetText(nuc + NUCLIDE_PTR_NAME),
                     (long)RDB[rea + REACTION_MT], loc1,
                     (long)RDB[loc1 + PREC_IDX], RDB[loc1 + PREC_LAMBDA]);
#endif
              /* Check if this precursor group is already in our list */
              /* Comparison based on lambda */

              for (i = 0; i < nprec; i++)
                {
                  /* Compare lambdas and break */

                  if (RDB[loc1 + PREC_LAMBDA] == lamlist[i])
                    break;

                }

              /* Check if this was found */

              if (i == nprec)
                {
                  /* Was not found, add to list */

                  lamlist[nprec]  = RDB[loc1 + PREC_LAMBDA];
                  realist[nprec]  = rea;
                  preclist[nprec] = loc1;

                  /* Increment number of precursor groups in list */

                  nprec++;

#ifdef DNPRINT
                  fprintf(outp,"Was not in list -> added\n");
#endif
                }
#ifdef DNPRINT
              else
                fprintf(outp,"Was already in list\n");
#endif
            }

          /* Next reaction */

          rea = NextItem(rea);
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Store number of precursor groups */

  WDB[loc0 + PRECDET_NG] = (double)nprec;

  /********************************************************/
  /* Store lambdas, precursor group pointers and reaction */
  /* pointers for each group                              */
  /********************************************************/

  /* Allocate memory for arrays */

  reaarray  = ReallocMem(DATA_ARRAY, nprec);
  precarray = ReallocMem(DATA_ARRAY, nprec);
  lamarray  = ReallocMem(DATA_ARRAY, nprec);

  /* Move data from temporary lists to precdet arrays */

  for (i = 0; i < nprec; i++)
    {
      WDB[reaarray + i]  = (double)realist[i];
      WDB[lamarray + i]  = (double)lamlist[i];
      WDB[precarray + i] = (double)preclist[i];
    }

  /* Store pointers to arrays */

  WDB[loc0 + PRECDET_PTR_REA_ARRAY]  = (double)reaarray;
  WDB[loc0 + PRECDET_PTR_PREC_ARRAY] = (double)precarray;
  WDB[loc0 + PRECDET_PTR_LAM_ARRAY]  = (double)lamarray;

  /* Free temporary lists */

  Mem(MEM_FREE, realist);
  Mem(MEM_FREE, preclist);
  Mem(MEM_FREE, lamlist);

  /**********************************************/
  /* Set the mesh size for precursor population */
  /**********************************************/

  /* Get mesh size either from "savesrc" card or */
  /* as defined in "dynsrc" file                 */

  if (RDB[loc0 + PRECDET_PTR_IN_FNAME] < VALID_PTR)
    {
      /* No initial precursor concentrations   */
      /* Use mesh size from savesrc parameters */
      /* Set at readinput.c */

      nx = (long)RDB[loc0 + PRECDET_N0];
      ny = (long)RDB[loc0 + PRECDET_N1];
      nz = (long)RDB[loc0 + PRECDET_N2];
    }
  else
    {
      /* Initial populations set by "set dynsrc" */
      /* Read meshing from <source>.main file */

      sprintf(tmpstr, "%s.main", GetText(loc0 + PRECDET_PTR_IN_FNAME));

      /* Try to open file */

      if ((fp = fopen(tmpstr, "r")) == NULL)
        Die(FUNCTION_NAME, "Dynamic source file %s cannot be opened", tmpstr);

      /* Loop over live neutron population and relerr */

      for (i = 0; i < 2; i++)
        tmplong = fscanf(fp, "%lf", &tmpdouble);

      /* Loop over time and group binning */

      for (i = 0; i < 2; i++)
        tmplong = fscanf(fp, "%ld", &tmplong);

      /* Get mesh binning */

      if (fscanf(fp, "%ld %ld %ld\n", &nx, &ny, &nz) != 3)
        Die(FUNCTION_NAME, "Could not read mesh size from file %s", tmpstr);

      fclose(fp);
    }

  /*********************/

  /* Get number of time intervals */

  nt = (long)RDB[DATA_DYN_NB];

  /* Store number of time bins */

  WDB[loc0 + PRECDET_NT] = (double)(nt + 1);

  /* Allocate memory for statistics */

  ptr = NewStat("PRECURSOR_DETECTOR", 3, nt + 1, nprec, nx*ny*nz);

  /* Store pointer to statistics */

  WDB[loc0 + PRECDET_PTR_STAT] = (double)ptr;

  /* Create mesh (now cartesian) */
  /* Allocate memory */

  loc1 = NewItem(loc0 + PRECDET_PTR_MESH, MESH_BLOCK_SIZE);

  /* Set mesh data */

  WDB[loc1 + MESH_N0] = (double)nx;
  WDB[loc1 + MESH_MIN0] = RDB[DATA_GEOM_MINX];
  WDB[loc1 + MESH_MAX0] = RDB[DATA_GEOM_MAXX];

  WDB[loc1 + MESH_N1] = (double)ny;
  WDB[loc1 + MESH_MIN1] = RDB[DATA_GEOM_MINY];
  WDB[loc1 + MESH_MAX1] = RDB[DATA_GEOM_MAXY];

  WDB[loc1 + MESH_N2] = (double)nz;
  WDB[loc1 + MESH_MIN2] = RDB[DATA_GEOM_MINZ];
  WDB[loc1 + MESH_MAX2] = RDB[DATA_GEOM_MAXZ];

  AllocValuePair(loc1 + MESH_PREV_COL_IDX);

  /* Put type */

  WDB[loc1 + MESH_TYPE] = (double)MESH_TYPE_CARTESIAN;

  /* Create list used for sampling mesh bin */

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_MESH)
    {

      for (i = 0; i < nprec*nx*ny*nz; i++)
        {
          ptr = NewItem(loc0 + PRECDET_PTR_MESH_LIST, VALUE_PAIR_BLOCK_SIZE);

          WDB[ptr + VALUE_PAIR_VAL1] = (double)i;
          WDB[ptr + VALUE_PAIR_VAL2] = 0.0;
        }

      /* Close list for sorting */

      CloseList(ptr);
    }
#ifdef DNPRINT
  fprintf(outp, "<-- processprecdet.c\n\n");
#endif
}
