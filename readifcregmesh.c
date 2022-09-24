/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcregmesh.c                               */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads regular mesh based multi-physics interfaces            */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCRegMesh:"

/*****************************************************************************/

void ReadIFCRegMesh(long ifc, long update)
{
  long loc1, msh, msh1, type, np, n, nmat, ptr, nmsh;
  long nx, ny, nz, i, j, arrsize, nr, tmplist, rholist;
  double zmin, zmax, dmax, Tmax, Tmin;
  double old, new, maxeps, maxdiff, L2abs, L2rel;
  double d, T, *lims, x0, y0, z0, pitchx, pitchy, pitchz;
  char tmpstr[MAX_STR];
  FILE *fp, *fout;

  /* Open file for reading */

  if ((fp = fopen(GetText(ifc + IFC_PTR_INPUT_FNAME), "r")) == NULL)
    Error(ifc, "Multi-physics interface file \"%s\" does not exist",
          GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Read interface type */

  if (fscanf(fp, "%ld", &type) == EOF)
    Die(FUNCTION_NAME, "Could not read interface type from %s",
        GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Read material name */

  if (fscanf(fp, "%s", tmpstr) == EOF)
    Die(FUNCTION_NAME, "Could not read material name from %s",
        GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Allocate memory for the material array and store material name */
  /* May be already given using setmat */

  if ((update == NO) && ((long)RDB[ifc + IFC_PTR_MAT_ARR] < VALID_PTR))
    {
      /* IFC files always contain only one material name */

      nmat = 1;

      /* Store number of materials */

      WDB[ifc + IFC_N_MAT] = (double)nmat;

      /* Allocate memory for material array */

      ptr = ReallocMem(DATA_ARRAY, nmat + 1);

      /* Store pointer to material array */

      WDB[ifc + IFC_PTR_MAT_ARR] = (double)ptr;

      /* Store material name */

      WDB[ptr++] = (double)PutText(tmpstr);

      /* Null terminate the array */

      WDB[ptr] = NULLPTR;
    }

  /* Read output flag */

  if (fscanf(fp, "%ld", &n) == EOF)
    Die(FUNCTION_NAME, "Could not read output flag from %s",
        GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Store output flag */

  WDB[ifc + IFC_CALC_OUTPUT] = (double)n;

  /* Read output file name and axial binning for output */

  if (n == IFC_OUTPUT_SAME_MESH)
    {
      WDB[ifc + IFC_CALC_OUTPUT] = (double)YES;

      /* Read just the file name */

      if (fscanf(fp, "%s", tmpstr) == EOF)
        Die(FUNCTION_NAME, "Could not read output parameters");

      /* Store output file name */

      if(!update)
        WDB[ifc + IFC_PTR_OUTPUT_FNAME] = (double)PutText(tmpstr);

      /* Store output type */

      WDB[ifc + IFC_OUTPUT_TYPE] = (double)IFC_OUTPUT_SAME_MESH;

      /* Loop over previous interfaces and check duplicate file names */

      loc1 = PrevItem(ifc);
      while (loc1 > VALID_PTR)
        {
          /* Compare file names */

          if ((long)RDB[loc1 + IFC_PTR_OUTPUT_FNAME] > VALID_PTR)
            if (CompareStr(ifc + IFC_PTR_OUTPUT_FNAME,
                           loc1 + IFC_PTR_OUTPUT_FNAME))
              Error(ifc,
                    "Duplicate output file name with distribution \"%s\"",
                    GetText(loc1 + IFC_PTR_INPUT_FNAME));

          /* Pointer to previous */

          loc1 = PrevItem(loc1);
        }

    }
  else if (n != NO)
    {
      /* Read file name and binning for other types */

      if (fscanf(fp, "%s %ld %lf %lf %ld", tmpstr, &nz, &zmin,
                 &zmax, &nr) == EOF)
        Die(FUNCTION_NAME, "Could not read output parameters");

      if(!update)
        WDB[ifc + IFC_PTR_OUTPUT_FNAME] = (double)PutText(tmpstr);

      /* Put number of axial zones */

      if (nz < 1)
        Error(ifc, "Error in number of axial zones");
      else
        WDB[ifc + IFC_NZ] = (double)nz;

      /* Put axial limits */

      if (zmin < zmax)
        {
          WDB[ifc + IFC_ZMIN] = zmin;
          WDB[ifc + IFC_ZMAX] = zmax;
        }
      else
        Error(ifc, "Error in axial boundaries");

      /* Put number of radial zones */

      if (nr < 1)
        Error(ifc, "Error in number of radial zones");
      else
        WDB[ifc + IFC_NR] = (double)nr;

      /* Loop over previous and check duplicate file names */

      loc1 = PrevItem(ifc);
      while (loc1 > VALID_PTR)
        {
          /* Compare file names */

          if ((long)RDB[loc1 + IFC_PTR_OUTPUT_FNAME] > VALID_PTR)
            if (CompareStr(ifc + IFC_PTR_OUTPUT_FNAME,
                           loc1 + IFC_PTR_OUTPUT_FNAME))
              Error(ifc,
                    "Duplicate output file name with distribution \"%s\"",
                    GetText(loc1 + IFC_PTR_INPUT_FNAME));

          /* Pointer to previous */

          loc1 = PrevItem(loc1);
        }
    }

  /* Reset limiting values */

  Tmax = -INFTY;
  Tmin = INFTY;
  dmax = 0.0;

  /* Read number of nested meshes */

  if (type == IFC_TYPE_REG_MESH_MULTILVL)
    {
      if (fscanf(fp, "%ld", &nmsh) == EOF)
        Die(FUNCTION_NAME, "Could not read number of nested multilevel"
            " meshes");
    }
  else
    nmsh = 1;

  /* Store number of nested meshes */

  WDB[ifc + IFC_N_MESH_LVL] = (double)nmsh;

  /* Reset number of cells */

  np = 1;

  /* Read all meshes */

  for (i = 0; i < nmsh; i++)
    {
      /* Get mesh type */

      if (fscanf(fp, "%ld", &n) == EOF)
        Die(FUNCTION_NAME, "Could not read mesh type");

      /* Avoid compiler warning */

      lims = NULL;
      x0 = 0;
      y0 = 0;
      z0 = 0;

      pitchx = 0;
      pitchy = 0;
      pitchz = 0;

      /* Check type */

      if (n == MESH_TYPE_CARTESIAN)
        {
          /* Cartesian, allocate memory for limits */

          lims = (double *)Mem(MEM_ALLOC, 6, sizeof(double));

          /* Read data */

          if (fscanf(fp, "%ld %lf %lf %ld %lf %lf %ld %lf %lf",
                     &nx, &lims[0], &lims[1], &ny, &lims[2], &lims[3],
                     &nz, &lims[4], &lims[5]) == EOF)
            Die(FUNCTION_NAME, "Could not read mesh data");

          x0 = 0.5*(lims[0] + lims[1]);
          y0 = 0.5*(lims[2] + lims[3]);
          z0 = 0.5*(lims[4] + lims[5]);

          pitchx = (lims[1] - lims[0])/(double)nx;
          pitchy = (lims[3] - lims[2])/(double)ny;
          pitchz = (lims[5] - lims[4])/(double)nz;

        }
      else if ((n == MESH_TYPE_HEXX) || (n == MESH_TYPE_HEXY))
        {
          /* Cartesian, allocate memory for limits */

          lims = (double *)Mem(MEM_ALLOC, 6, sizeof(double));

          /* Read data */

          if (fscanf(fp, "%lf %lf %lf %lf %lf %ld %ld %ld",
                     &lims[0], &lims[2], &lims[1], &lims[4], &lims[5],
                     &nx, &ny, &nz) == EOF)
            Die(FUNCTION_NAME, "Could not read mesh data");

          x0 = lims[0];
          y0 = lims[2];
          z0 = 0.5*(lims[4] + lims[5]);

          pitchx = lims[1];
          pitchy = lims[1];
          pitchz = (lims[5] - lims[4])/(double)nz;
        }
      else if (n == MESH_TYPE_ORTHOGONAL)
        {
          /* Orthogonal mesh */

          if (fscanf(fp, "%ld %ld %ld", &nx, &ny, &nz) == EOF)
            Die(FUNCTION_NAME, "fscanf error");

          /* Allocate memory */

          lims = (double *)Mem(MEM_ALLOC, nx + ny + nz + 3, sizeof(double));

          /* Read data (NOTE: noiden järjestys pitäis tarkistaa) */

          for (i = 0; i < nx + ny + nz + 3; i++)
            if (fscanf(fp, "%lf", &lims[i]) == EOF)
              Die(FUNCTION_NAME, "fscanf error");
        }
      else
        Error(ifc, "Invalid mesh type %ld", n);

      /* Calculate number of values */

      np = np*nx*ny*nz;

      /* Check */

      if (np < 1)
        Error(ifc, "Zero mesh size");

      if(!update)
        {
          /* Create blank mesh structure */

          msh = NewItem(ifc + IFC_PTR_SEARCH_MESH_LIST, MESH_BLOCK_SIZE);

          /* Create correct mesh structure */

          msh1 = CreateMesh(n, MESH_CONTENT_PTR, -1, nx, ny, nz, lims, -1);

          /* Copy content to the one in the list */

          memcpy(&WDB[msh + LIST_DATA_SIZE], &RDB[msh1 + LIST_DATA_SIZE],
                (MESH_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

          /* Create lattice for multi-level mesh */

          if (type == IFC_TYPE_REG_MESH_MULTILVL)
            {
              /* Create new lattice item */

              loc1 = NewItem(ifc + IFC_PTR_LATTICE_LIST, LAT_BLOCK_SIZE);

              /* Print name */

              sprintf(tmpstr, "Interface %s level %ld",
                      GetText(ifc + IFC_PTR_INPUT_FNAME), i + 1);

              /* Put name */

              WDB[loc1 + LAT_PTR_NAME] = (double)PutText(tmpstr);

              /* Put type */

              if (n == MESH_TYPE_CARTESIAN)
                WDB[loc1 + LAT_TYPE] = (double)LAT_TYPE_CUBOID;
              else if (n == MESH_TYPE_HEXX)
                WDB[loc1 + LAT_TYPE] = (double)LAT_TYPE_XPRISM;
              else if (n == MESH_TYPE_HEXY)
                WDB[loc1 + LAT_TYPE] = (double)LAT_TYPE_YPRISM;
              else
                Error(ifc, "Mesh type %ld not allowed in multi-level mesh "
                      "interfaces");

              /* Put lattice origin */

              WDB[loc1 + LAT_ORIG_X0] = x0;
              WDB[loc1 + LAT_ORIG_Y0] = y0;
              WDB[loc1 + LAT_ORIG_Z0] = z0;

              /* Reset pointer */

              ptr = -1;

              /* Check type */

              /* Lattice size */

              WDB[loc1 + LAT_NX] = (double)nx;
              WDB[loc1 + LAT_NY] = (double)ny;
              WDB[loc1 + LAT_NZ] = (double)nz;
              WDB[loc1 + LAT_NTOT] = (double)(nx*ny*nz);

              /* Lattice pitch */

              WDB[loc1 + LAT_PITCHX] = pitchx;
              WDB[loc1 + LAT_PITCHY] = pitchy;
              WDB[loc1 + LAT_PITCHZ] = pitchz;

              /* Allocate memory for data */

              ptr = ReallocMem(DATA_ARRAY, nx*ny*nz + 1);

              /* Set pointer */

              WDB[loc1 + LAT_PTR_FILL] = ptr;

              /* Put items */

              for (j = 0; j < nx*ny*nz; j++)
                WDB[ptr++] = (double)PutText("-1");

              /* Put null pointer */

              WDB[ptr] = NULLPTR;

              /***************************************************************/
            }
        }
      /* Free allocated memory */

      Mem(MEM_FREE, lims);
    }

  /* Store number of cells */

  WDB[ifc + IFC_NP] = (double)np;

  arrsize = np;

  /* Increase array size if needed for remapping */

  if (!update)
    {

      WDB[ifc + IFC_ARR_SIZE] = (double)arrsize;
    }
  else
    arrsize = (long)RDB[ifc + IFC_ARR_SIZE];

  /* Reset maximum of convergence criterion */

  maxeps  = 0.0;
  maxdiff = 0.0;
  L2abs   = 0.0;
  L2rel   = 0.0;

  /* Allocate or get pointer to T and rho lists */

  if (!update)
    {
      fprintf(outp, "Allocating list sizes of %ld\n", arrsize);
      /* Allocate memory for temperatures */

      tmplist = ReallocMem(DATA_ARRAY, arrsize);

      /* Put temperature list pointer to memory */

      WDB[ifc + IFC_PTR_TMP_LIST] = (double)tmplist;

      /* Allocate memory for densities */

      rholist = ReallocMem(DATA_ARRAY, arrsize);

      /* Put density list pointer to memory */

      WDB[ifc + IFC_PTR_DF_LIST] = (double)rholist;

    }
  else
    {
      /* Get temperature list pointer from memory */

      tmplist = (long)RDB[ifc + IFC_PTR_TMP_LIST];

      /* Get density list pointer from memory */

      rholist = (long)RDB[ifc + IFC_PTR_DF_LIST];
    }

  /* Loop over points and read data */

  n = 0;

  for (n = 0; n < np; n++)
    {
      /* Read values */

      if (fscanf(fp, "%lf %lf", &d, &T) == EOF)
        Die(FUNCTION_NAME, "fscanf error for value pair %ld", n);

      /* Compare to limits */

      if (fabs(d) > fabs(dmax))
        dmax = d;

      /* Translate negative temperatures to zero (no TMS) */
      if (T < 0)
        T = 0;

      if (T > Tmax)
        Tmax = T;

      if (T < Tmin)
        Tmin = T;


      /**********************************************************/
      /* Convergence criterions based on momentary distribution */
      /**********************************************************/

      if (update)
        {
          old = RDB[tmplist + n];
          new = T;

          CalcConvCriteria(new, old,
                           &maxdiff, &maxeps, &L2abs, &L2rel);
        }

      /* Put data */

      WDB[rholist + n] = d;
      WDB[tmplist + n] = T;

      /* Update index */
    }

  /* Put maximum density and temperature                  */
  /* For updates these are checked in processifcregmesh.c */

  WDB[ifc + IFC_MAX_DENSITY] = dmax;
  WDB[ifc + IFC_MAX_TEMP] = Tmax;
  WDB[ifc + IFC_MIN_TEMP] = Tmin;

  /* Write some values to the end of the arrays */

  for (n = np; n < arrsize; n++)
    {
      WDB[rholist + n] = RDB[ifc + IFC_MAX_DENSITY];
      WDB[tmplist + n] = RDB[ifc + IFC_MAX_TEMP];
    }

  /* Set TMS on */

  if (Tmax > 0.0)
    if(!update)
      WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;


  if (update)
    {

      /* Open output file for convergence */

      if (WDB[DATA_SOL_REL_ITER] == (double)0)
        {
          /* On first iteration open a new file */

          sprintf(tmpstr, "%s_Tconv%ld.m", GetText(ifc + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fout = fopen(tmpstr, "w");

          /* Reset idx */

          fprintf(fout,"\nidx = 1;\n\n");
        }
      else
        {
          /* On subsequent iterations append to old file */

          sprintf(tmpstr, "%s_Tconv%ld.m", GetText(ifc + IFC_PTR_INPUT_FNAME),
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
