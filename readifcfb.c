/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcfb.c                                    */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2018/11/07 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads fuel behavior multi-physics interfaces                 */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCFB:"

/*****************************************************************************/

void ReadIFCFB(long loc0, long update)
{
  long loc1, loc2, loc3, ptr, type, np, n, m, i, nr, nax, na, uni, l;
  long DFptr, Tptr, CRptr, HRptr;
  double zmin, zmax, T, r1, r2, amin, amax, old, new;
  double maxeps, maxdiff, L2abs, L2rel;
  char tmpstr[MAX_STR];
  FILE *fp, *fout;

  /* Reset maximum of convergence criterion */

  maxeps  = 0.0;
  maxdiff = 0.0;
  L2abs   = 0.0;
  L2rel   = 0.0;

  /* Avoid compiler warning */

  loc1 = -1;
  loc3 = -1;

  /* Open file for reading */

  if ((fp = fopen(GetText(loc0 + IFC_PTR_INPUT_FNAME), "r")) == NULL)
    Error(loc0, "Multi-physics interface file \"%s\" does not exist",
          GetText(loc0 + IFC_PTR_INPUT_FNAME));

  /* Read interface type */

  if (fscanf(fp, "%ld", &type) == EOF)
    Die(FUNCTION_NAME, "fscanf error");


  /* Interface to fuel performance codes */

  WDB[loc0 + IFC_TYPE] = (double)type;
  WDB[loc0 + IFC_PTR_MAT_ARR] = NULLPTR;
  WDB[loc0 + IFC_CALC_OUTPUT] = (double)YES;

  /* Reset output flag to avoid going to next loop */

  n = NO;

  /* Read output file name */

  if (fscanf(fp, "%s", tmpstr) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  WDB[loc0 + IFC_PTR_OUTPUT_FNAME] = (double)PutText(tmpstr);

  /*******************************************************************/

  /***** Interface to fuel performance codes *************************/

  /* Reset maximum relative difference */

  maxeps = 0;

  /* Read number of pins */

  if (fscanf(fp, "%ld", &np) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Loop */

  for (i = 0; i < np; i++)
    {
      /* Allocate memory for structure */
      if (!update)
        loc1 = NewItem(loc0 + IFC_PTR_FUEP, IFC_FUEP_LIST_BLOCK_SIZE);
      else if (i>0)
        loc1 = NextItem(loc1);
      else
        loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];

      /* Read universe */

      if (fscanf(fp, "%s", tmpstr) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      na = atol(tmpstr);

      if (na < 0)
        {
          /* Segmented rod */

          na = -na;
          uni = 0;
        }
      else if (na == 0)
        {
          /* Rod universe name is not a number */

          uni = -1;

          /* One axial segment */

          na = 1;
        }
      else
        {
          /* Rod universe is a number (na == uni) */

          uni = na;

          /* One axial segment */

          na = 1;
        }

      /* Store ifc type to rod block */

      WDB[loc1 + IFC_FUEP_TYPE] = RDB[loc0 + IFC_TYPE];

      /* Store number of pin segments */

      WDB[loc1 + IFC_FUEP_N_UNI] = (double)na;

      /* Allocate memory for pin universes */
      if (!update)
        {
          ptr = ReallocMem(DATA_ARRAY, na);
          WDB[loc1 + IFC_FUEP_PTR_UNI_LIST] = (double)ptr;
        }
      else
        ptr = (long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST];

      if(uni != 0)
        {
          WDB[ptr] = (double)PutText(tmpstr);
        }
      else
        {
          for(n = 0; n < na; n++)
            {
              if(fscanf(fp, "%s", tmpstr) == EOF)
                Die(FUNCTION_NAME, "fscanf error");

              WDB[ptr++] = (double)PutText(tmpstr);
            }
        }

      /* Read output meshing */

      ReadIFCFBLims(fp, loc1, update);

      /* Read number of axial zones */

      if (fscanf(fp, "%ld", &nax) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      /* Put value */

      WDB[loc1 + IFC_FUEP_N_AX] = (double)nax;

      /* Avoid compiler warning */

      loc2 = -1;

      /* Loop over axial zones */

      for (n = 0; n < nax; n++)
        {
          /* Allocate memory */

          if (!update)
            loc2 = NewItem(loc1 + IFC_FUEP_PTR_AX,
                           IFC_FUEP_AX_BLOCK_SIZE);
          else if (n > 0)
            loc2 = NextItem(loc2);
          else
            loc2 = (long)RDB[loc1 + IFC_FUEP_PTR_AX];

          /* Read minimum and maximum height */

          if (fscanf(fp, "%lf %lf", &zmin, &zmax) == EOF)
            Die(FUNCTION_NAME, "fscanf error");

          /* Put values */

          WDB[loc2 + IFC_FUEP_AX_ZMIN] = zmin;
          WDB[loc2 + IFC_FUEP_AX_ZMAX] = zmax;

          /* Read number of angular zones */

          if (fscanf(fp, "%ld", &na) == EOF)
            Die(FUNCTION_NAME, "fscanf error");

          /* Put value */

          WDB[loc2 + IFC_FUEP_AX_N_ANG] = (double)na;

          /* Loop over angular zones */

          for (m = 0; m < na; m++)
            {
              /* Allocate memory */

              if(!update)
                loc3 = NewItem(loc2 + IFC_FUEP_AX_PTR_ANG,
                               IFC_FUEP_ANG_BLOCK_SIZE);
              else if(m > 0)
                loc3 = NextItem(loc3);
              else
                loc3 = (long)RDB[loc2 + IFC_FUEP_AX_PTR_ANG];

              /* Read minimum and maximum angle */

              if (fscanf(fp, "%lf %lf", &amin, &amax) == EOF)
                Die(FUNCTION_NAME, "Was expecting angular input limits");

              /* Put values */

              WDB[loc3 + IFC_FUEP_ANG_AMIN] = amin*2*PI/360.0;
              WDB[loc3 + IFC_FUEP_ANG_AMAX] = amax*2*PI/360.0;

              /* Read number of radial points */

              if (fscanf(fp, "%ld", &nr) == EOF)
                Die(FUNCTION_NAME, "Was expecting the number of radial points");

              /* Put number */

              WDB[loc3 + IFC_FUEP_ANG_N_RAD] = (double)nr;

              if (fscanf(fp, "%lf %lf %lf", &r1, &r2, &T)
                  == EOF)
                Die(FUNCTION_NAME, "fscanf error");

              if(r1 != 0.0)
                nr=nr+1;

              /* Add point outside cladding */

              nr = nr + 1;

              /* Store total number of points */

              WDB[loc3 + IFC_FUEP_ANG_N_RAD] = (double)nr;

              /* Allocate arrays for the points */

              if (!update)
                {
                  /* Density factor */

                  DFptr = ReallocMem(DATA_ARRAY, nr);
                  WDB[loc3 + IFC_FUEP_ANG_PTR_DF]     = (double)DFptr;
                  WDB[loc3 + IFC_FUEP_ANG_PTR_DF_BOI] = (double)DFptr;

                  /* Temperature at the end of time step */

                  Tptr  = ReallocMem(DATA_ARRAY, nr);
                  WDB[loc3 + IFC_FUEP_ANG_PTR_TMP]     = (double)Tptr;
                  WDB[loc3 + IFC_FUEP_ANG_PTR_TMP_BOI] = (double)Tptr;

                  /* Cold radius */

                  CRptr = ReallocMem(DATA_ARRAY, nr);
                  WDB[loc3 + IFC_FUEP_ANG_PTR_COLD_R2] = (double)CRptr;

                  /* Hot radius */

                  HRptr = ReallocMem(DATA_ARRAY, nr);
                  WDB[loc3 + IFC_FUEP_ANG_PTR_HOT_R2]     = (double)HRptr;
                  WDB[loc3 + IFC_FUEP_ANG_PTR_HOT_R2_BOI] = (double)HRptr;
                }
              else
                {
                  /* Density factor */

                  DFptr = (long)RDB[loc3 + IFC_FUEP_ANG_PTR_DF];

                  /* Temperature */

                  Tptr = (long)RDB[loc3 + IFC_FUEP_ANG_PTR_TMP];

                  /* Cold radius */

                  CRptr = (long)RDB[loc3 + IFC_FUEP_ANG_PTR_COLD_R2];

                  /* Hot radius */

                  HRptr = (long)RDB[loc3 + IFC_FUEP_ANG_PTR_HOT_R2];

                }

              /* BOI values are copied in processinterface */

              if (r1 != 0.0)
                {
                  /* No centerline data given in interface */

                  if (update == YES)
                    {
                      /**********************************************************/
                      /* Convergence criterions based on momentary distribution */
                      /**********************************************************/

                      /* Get previous iteration CL temperature */

                      old = RDB[Tptr];

                      /* Get new iteration CL temperature */

                      new = T;

                      /* Calculate convergence criteria */

                      CalcConvCriteria(new, old,
                                       &maxdiff, &maxeps, &L2abs, &L2rel);

                    }

                  /* Write CL temperature as the same as first radial node */

                  WDB[Tptr]  = T;
                  WDB[CRptr] = 0.0;
                  WDB[HRptr] = 0.0;

                  /* Write first radial node data */

                  WDB[Tptr+1]  = T;
                  WDB[CRptr+1] = r1;
                  WDB[HRptr+1] = r2;

                  /* We have written data for two nodes now */

                  l = 2;
                }
              else
                {

                  /* Centerline data given in interface */

                  if (update == YES)
                    {
                      /* Get previous iteration CL temperature */

                      old = RDB[Tptr];

                      /* Get new iteration CL temperature */

                      new = T;

                      /* Calculate convergence criteria */

                      CalcConvCriteria(new, old,
                                       &maxdiff, &maxeps, &L2abs, &L2rel);
                    }

                  /* Write CL data based on interface file */

                  WDB[Tptr]  = T;
                  WDB[CRptr] = r1;
                  WDB[HRptr] = r2;

                  /* We have written data for one node now */

                  l = 1;

                }

              for (; l < nr - 1; l++)
                {

                  /* Read next coldradius, hotradius and temperature */

                  if (fscanf(fp, "%lf %lf %lf", &r1, &r2, &T)
                      == EOF)
                    Die(FUNCTION_NAME, "fscanf error");

                  /* Check that the radii are in ascending order */

                  if (r1 < RDB[CRptr + l - 1])
                    Error(loc0, "Radii must be in ascending order %f %f",
                          RDB[CRptr + l - 1], r1);

                  if (r2 < RDB[HRptr + l - 1])
                    Error(loc0, "Radii must be in ascending order %f %f",
                          RDB[HRptr + l - 1], r2);

                  if (update == YES)
                    {
                      /**********************************************************/
                      /* Convergence criterions based on momentary distribution */
                      /**********************************************************/

                      /* Get previous iteration temperature */

                      old = RDB[Tptr + l];

                      /* Get new iteration temperature */

                      new = T;

                      /* Calculate convergence criteria */

                      CalcConvCriteria(new, old,
                                       &maxdiff, &maxeps, &L2abs, &L2rel);

                    }

                  /* Put values for this radial node */

                  WDB[DFptr + l] = 1.0;
                  WDB[Tptr + l]  = T;
                  WDB[CRptr + l] = r1;
                  WDB[HRptr + l] = r2;

                  /* Next radial point */
                }

              /* Add static point outside cladding       */
              /* This allows to do coordinate transforms */
              /* Even if cladding creeps inwards         */

              WDB[DFptr + l] = 1.0;
              WDB[Tptr + l]  = T;
              WDB[CRptr + l] = r1*2.0;
              WDB[HRptr + l] = r1*2.0;

              /* Next angular region */
            }

          /* Next axial region */
        }

    }

  /* Set TMS on */
  if (!update)
    WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;

  if (update)
    {

      /* Open output file for convergence */

      if (WDB[DATA_SOL_REL_ITER] == (double)0)
        {
          /* On first iteration open a new file */

          sprintf(tmpstr, "%s_Tconv%ld.m", GetText(loc0 + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fout = fopen(tmpstr, "w");

          /* Reset idx */

          fprintf(fout,"\nidx = 1;\n\n");
        }
      else
        {
          /* On subsequent iterations append to old file */

          sprintf(tmpstr, "%s_Tconv%ld.m", GetText(loc0 + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fout = fopen(tmpstr, "a");

          /* Increment idx */

          fprintf(fout,"\nidx = idx + 1;\n\n");
        }

      /* Write out convergence criteria of momentary power distribution */

      fprintf(fout, "T_eps(idx) = %E;\n", maxeps);
      fprintf(fout, "T_delta(idx) = %E;\n", maxdiff);
      fprintf(fout, "T_L2_of_absolute(idx) = %E;\n", sqrt(L2abs));
      fprintf(fout, "T_L2_of_relative(idx) = %E;\n", sqrt(L2rel));

      /* Close output file */

      fclose(fout);
    }

  /*******************************************************************/

  /* Close file */

  fclose(fp);

}

/*****************************************************************************/
