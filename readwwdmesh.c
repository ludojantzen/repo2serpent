/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readwwdmesh.c                                  */
/*                                                                           */
/* Created:       2016/07/04 (JLe)                                           */
/* Last modified: 2019/09/26 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reads weight window mesh from file                           */
/*                                                                           */
/* Comments: - MCNP format does not yet support cylindrical mesh or fine     */
/*             structure                                                     */
/*                                                                           */
/*           - Toi MCNP meshin luku tehdään kahteen kertaan että neutronien  */
/*             ja fotonien meshit saadaan eri rakenteisiin --> mahdollistaa  */
/*             myöhemmin eri geometrian käytön muissa tyypeissä.             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadWWDMesh:"

/* Recursive read function */

long ReadWWMesh0(FILE *fp, long rmx, long wwd);

/*****************************************************************************/

void ReadWWDMesh(long wwd)
{
  long msh, rmx, loc0, loc1, ni, ptr, type, n, nfx, nfy, nfz, ncx, ncy, ncz;
  long i, j, k, np, m, ne[2], mty, nx, ny, nz, ng, i0, j0, k0, idx, nr, sz, nc;
  long nmax;
  double x0, y0, z0, f, r, *lims, max;
  char fname[MAX_STR], line[MAX_STR];
  FILE *fp;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "wwd", DATA_ARRAY, wwd);

  /* Check file name */

  if ((long)RDB[wwd + WWD_PTR_FNAME] < VALID_PTR)
    Die(FUNCTION_NAME, "File name not provided");

  /* Get file name */

  sprintf(fname, "%s", GetText(wwd + WWD_PTR_FNAME));

  /* Print */

  fprintf(outp, "Reading weight-window mesh from file \"%s\"...\n", fname);

  /* Check type */

  if ((type = (long)RDB[wwd + WWD_TYPE]) == WWD_MESH_TYPE_MCNP)
    {
      /***********************************************************************/

      /***** MCNP Weight window file (Appendix J in MCNP4C manual) ***********/

      /* Loop over particle types */

      for (ni = 0; ni < 2; ni++)
        {
          /* Open file for reading */

          if ((fp = fopen(fname, "r")) == NULL)
            Error(wwd, "Weight window data file \"%s\" does not exist", fname);

          /* Check format */

          TestDOSFile(fname);

          /* First line may or may not contain comments, so it must be read */
          /* into a string variable */

          if (fgets(line, 82, fp) == NULL)
            Die(FUNCTION_NAME, "fgets error");

          /* Read four first variables (second one is unused) */

          if (sscanf(line, "%ld %ld %ld %ld", &n, &m, &np, &nr) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");

          /* Check type (must be 1) */

          if (n != 1)
            Error(wwd, "Unsupported file type %ld", n);

          /* Check number of particle types */

          if (np > 2)
            Error(wwd, "Multiple particle types not supported");

          /* Check number of mesh parameters */

          if (nr != 10)
            Error(wwd, "Only rectangular mesh supported");

          /* Duplicate item */

          if (ni > 0)
            {
              /* Duplicate structure */

              ptr = DuplicateItem(wwd);
              wwd = ptr;

              /* Reset mesh pointers */

              WDB[wwd + WWD_PTR_MESH] = NULLPTR;
              WDB[wwd + WWD_PTR_MESH_DATA] = NULLPTR;

              /* Reset file name to avoid second call from ProcessVR() */

              WDB[wwd + WWD_PTR_FNAME] = NULLPTR;
            }

          /* Put particle type */

          if ((np > 1) && (ni == 0))
            WDB[wwd + WWD_PARTICLE_TYPE] = (double)PARTICLE_TYPE_NEUTRON;
          else if ((np > 1) && (ni == 1))
            WDB[wwd + WWD_PARTICLE_TYPE] = (double)PARTICLE_TYPE_GAMMA;

          /* Loop over particle types */

          for (m = 0; m < np; m++)
            {
              /* Read number of energies to temporary array */

              if (fscanf(fp, "%ld", &ne[m]) == EOF)
                Error(wwd, "Unexpected EOF in weight window data file");
              else if ((m > 0) && (ne[m] < 1))
                Error(wwd, "Invalid number of photon energy groups");
              else if (ne[m] < 0)
                Error(wwd, "Invalid number of neutron energy groups");

              CheckValue(FUNCTION_NAME, "ne", "", ne[m], 0, 10000);
            }

          /* Skip neutron data if number of groups is zero */

          if ((ni == 0) && (ne[ni] == 0))
            continue;

          /* Avoid compiler warning */

          nfx = 0;
          nfy = 0;
          nfz = 0;

          ncx = 0;
          ncy = 0;
          ncz = 0;

          /* Read mesh parameters (for some reason these are all in float) */

          if (fscanf(fp, "%lf", &f) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");
          else
            nfx = (long)f;

          if (fscanf(fp, "%lf", &f) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");
          else
            nfy = (long)f;

          if (fscanf(fp, "%lf", &f) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");
          else
            nfz = (long)f;

          if (fscanf(fp, "%lf %lf %lf", &x0, &y0, &z0) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");

          if (fscanf(fp, "%lf", &f) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");
          else
            ncx = (long)f;

          if (fscanf(fp, "%lf", &f) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");
          else
            ncy = (long)f;

          if (fscanf(fp, "%lf", &f) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");
          else
            ncz = (long)f;

          /* Read geometry type */

          if (fscanf(fp, "%lf", &f) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");
          else if (f != 1.0)
            Error(wwd, "Only reactangular mesh supported %lf", f);

          /* Check some values */

          if (nfx != ncx)
            Die(FUNCTION_NAME, "nfx = %ld ncx = %ld", nfx, ncx);
          else if (nfy != ncy)
            Die(FUNCTION_NAME, "nfy = %ld ncy = %ld", nfy, ncy);
          else if (nfz != ncz)
            Die(FUNCTION_NAME, "nfz = %ld ncz = %ld", nfz, ncz);

          CheckValue(FUNCTION_NAME, "ncx", "", ncx, 1, 10000);
          CheckValue(FUNCTION_NAME, "ncy", "", ncy, 1, 10000);
          CheckValue(FUNCTION_NAME, "ncz", "", ncz, 1, 10000);

          /* Allocate memory for mesh data */

          lims = (double *)Mem(MEM_ALLOC, ncx + ncy + ncz + 3,
                               sizeof(double));

          /*******************************************************************/

          /***** Read mesh boundaries ****************************************/

          /* Reset mesh index */

          idx = 0;

          /* Loop over x-boundaries */

          if (fscanf(fp, "%lf", &lims[idx++]) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");

          for (n = 0; n < ncx; n++)
            {
              /* Read parameters */

              if (fscanf(fp, "%lf %lf %lf", &f, &lims[idx++], &r) == EOF)
                Error(wwd, "Unexpected EOF in weight window data file");

              /* Check for fine mesh parameters */

              if ((f != 1.0) || (r != 1.0))
                Error(wwd, "Fine mesh structure not supported");
            }

          /* Loop over y-boundaries */

          if (fscanf(fp, "%lf", &lims[idx++]) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");

          for (n = 0; n < ncy; n++)
            {
              /* Read parameters */

              if (fscanf(fp, "%lf %lf %lf", &f, &lims[idx++], &r) == EOF)
                Error(wwd, "Unexpected EOF in weight window data file");

              /* Check for fine mesh parameters */

              if ((f != 1.0) || (r != 1.0))
                Error(wwd, "Fine mesh structure not supported");
            }

          /* Loop over z-boundaries */

          if (fscanf(fp, "%lf", &lims[idx++]) == EOF)
            Error(wwd, "Unexpected EOF in weight window data file");

          for (n = 0; n < ncz; n++)
            {
              /* Read parameters */

              if (fscanf(fp, "%lf %lf %lf", &f, &lims[idx++], &r) == EOF)
                Error(wwd, "Unexpected EOF in weight window data file");

              /* Check for fine mesh parameters */

              if ((f != 1.0) || (r != 1.0))
                Error(wwd, "Fine mesh structure not supported");
            }

          /* Check order */

          idx = 0;

          for (n = 0; n < ncx + 1; n++)
            {
              if (n > 0)
                if (lims[idx - 1] > lims[idx])
                  Error(wwd, "Error in weight window mesh bounds");

              idx++;
            }

          for (n = 0; n < ncy + 1; n++)
            {
              if (n > 0)
                if (lims[idx - 1] > lims[idx])
                  Error(wwd, "Error in weight window mesh bounds");

              idx++;
            }

          for (n = 0; n < ncz + 1; n++)
            {
              if (n > 0)
                if (lims[idx - 1] > lims[idx])
                  Error(wwd, "Error in weight window mesh bounds");

              idx++;
            }

          /*******************************************************************/

          /**** Skip first particle data *************************************/

          /* Loop over previous particle types */

          for (m = 0; m < ni; m++)
            {
              /* Energy array */

              for (n = 0; n < ne[m]; n++)
                if (fscanf(fp, "%lf", &f) == EOF)
                  Error(wwd, "Unexpected EOF in weight window data file");

              /* Data */

              for (n = 0; n < ne[m]; n++)
                for (k = 0; k < ncz; k++)
                  for (j = 0; j < ncy; j++)
                    for (i = 0; i < ncx; i++)
                      if (fscanf(fp, "%lf", &f) == EOF)
                        Error(wwd,
                              "Unexpected EOF in weight window data file");
            }

          /*******************************************************************/

          /***** Read energy group structure *********************************/

          /* Put number of energy groups */

          WDB[wwd + WWD_NE] = (double)ne[ni];

          /* Allocate memory */

          ptr = ReallocMem(DATA_ARRAY, ne[ni] + 1);
          WDB[wwd + WWD_PTR_ERG] = (double)ptr;

          /* Put minimum energy (use zero just to be safe) */

          WDB[ptr] = 0.0;

          /* Loop over values */

          for (n = 0; n < ne[ni]; n++)
            {
              /* Read value */

              if (fscanf(fp, "%lf", &WDB[ptr + n + 1]) == EOF)
                Error(wwd, "Unexpected EOF in weight window data file");

              /* Check */

              if (RDB[ptr + n + 1] <= RDB[ptr + n])
                Error(wwd, "Energy group boundaries not in ascending order");
            }

          /*******************************************************************/

          /***** Create mesh and read data ***********************************/

          /* Create mesh */

          msh = CreateMesh(MESH_TYPE_ORTHOGONAL, MESH_CONTENT_PTR, -1,
                           ncx, ncy, ncz, lims, -1);
          WDB[wwd + WWD_PTR_MESH] = (double)msh;

          /* Free memory */

          Mem(MEM_FREE, lims);

          /* Allocate memory */

          for (k = 0; k < ncz; k++)
            for (j = 0; j < ncy; j++)
              for (i = 0; i < ncx; i++)
                {
                  /* Create structure */
                  /*
                  Die(FUNCTION_NAME, "Tää tehtiin jo");
                  */
                  loc0 = NewItem(wwd + WWD_PTR_MESH_DATA, WWD_MESH_BLOCK_SIZE);

                  /* Put pointer */

                  ptr = ReadMeshPtr(msh, i, j, k);
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                  WDB[ptr] = (double)loc0;

                  /* Allocate memory for values importance vector */

                  ptr = ReallocMem(DATA_ARRAY, ne[ni]);
                  WDB[loc0 + WWD_MESH_PTR_IMP] = (double)ptr;
                }

          /* Read data */

          for (n = 0; n < ne[ni]; n++)
            for (k = 0; k < ncz; k++)
              for (j = 0; j < ncy; j++)
                for (i = 0; i < ncx; i++)
                  {
                    /* Get pointer */

                    ptr = ReadMeshPtr(msh, i, j, k);
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    ptr = (long)RDB[ptr];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    ptr = (long)RDB[ptr + WWD_MESH_PTR_IMP];
                    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                    /* Read value */

                    if (fscanf(fp, "%lf", &f) == EOF)
                      Error(wwd, "Unexpected EOF in weight window data file");

                    /* Check value and store */
                    /*
                    CheckValue(FUNCTION_NAME, "f", "", f, 1E-100, 1E+100);
                    WDB[ptr + n] = RDB[DATA_WWD_LOWER_BOUND]/f;
                    */

                    /* Muutettu 9.8.2017 / 2.1.30 / JLe */

                    if (f > 0.0)
                      WDB[ptr + n] = RDB[DATA_WWD_LOWER_BOUND]/f;
                    else if (f < 0.0)
                      WDB[ptr + n] = 1E-100;
                    else
                      WDB[ptr + n] = 0.0;
                  }

                /* Check for major differences */

          if (1 == 2)
            {
              max = -1.0;

              for (k = 0; k < ncz; k++)
                for (j = 0; j < ncy; j++)
                  for (i = 0; i < ncx; i++)
                    {
                      /* Get pointer */

                      ptr = ReadMeshPtr(msh, i, j, k);
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      ptr = (long)RDB[ptr];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      ptr = (long)RDB[ptr + WWD_MESH_PTR_IMP];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      for (n = 1; n < ne[ni]; n++)
                        if ((RDB[ptr + n] > 0.0) && (RDB[ptr + n - 1] > 0.0))
                          if (fabs(log10(RDB[ptr + n]/RDB[ptr + n - 1])) > max)
                            max = fabs(log10(RDB[ptr + n]/RDB[ptr + n - 1]));
                    }

              for (n = 0; n < ne[ni]; n++)
                for (k = 0; k < ncz; k++)
                  for (j = 0; j < ncy; j++)
                    for (i = 0; i < ncx; i++)
                      {
                        /* Get pointer */

                        ptr = ReadMeshPtr(msh, i, j, k);
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                        ptr = (long)RDB[ptr];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                        ptr = (long)RDB[ptr + WWD_MESH_PTR_IMP];
                        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                        /* Get value */

                        f = RDB[ptr + n];

                        /* Check with neighbours */

                        if (f > 0.0)
                          for (i0 = i - 1; i0 < i + 2; i0++)
                          for (j0 = j - 1; j0 < j + 2; j0++)
                          for (k0 = k - 1; k0 < k + 2; k0++)
                            if ((i0 > -1) && (j0 > -1) && (k0 > -1) &&
                                (i0 < ncx) && (j0 < ncy) && (k0 < ncz))
                              {
                                /* Get pointer */

                                ptr = ReadMeshPtr(msh, i0, j0, k0);
                                CheckPointer(FUNCTION_NAME, "(ptr)",
                                             DATA_ARRAY, ptr);

                                ptr = (long)RDB[ptr];
                                CheckPointer(FUNCTION_NAME, "(ptr)",
                                             DATA_ARRAY, ptr);

                                ptr = (long)RDB[ptr + WWD_MESH_PTR_IMP];
                                CheckPointer(FUNCTION_NAME, "(ptr)",
                                             DATA_ARRAY, ptr);

                                if (RDB[ptr + n] > 0.0)
                                  if (fabs(log10(RDB[ptr + n]/f)) > max)
                                    max = fabs(log10(RDB[ptr + n]/f));
                              }
                      }

              /* Check */

              if ((max > 4.0) && (np == 1))
                Note(wwd,
                     "Maximum increase in importance is very large (> 1E%1.0f)",
                     max);
              else if ((max > 4.0) && (np > 1) && (ni == 0))
                Note(wwd,
                     "Maximum increase in neutron importance is very large (> 1E%1.0f)",
                     max);
              else if ((max > 4.0) && (np > 1) && (ni == 1))
                Note(wwd,
                     "Maximum increase in photon importance is very large (> 1E%1.0f)",
                     max);
            }

          /*******************************************************************/

          /* Check that everything is read */

          if (ni == np)
            if (fscanf(fp, "%lf", &f) != EOF)
              Die(FUNCTION_NAME, "EOF not reached");

          /* Close file */

          fclose(fp);

          /* Break loop if only one particle type is given */

          if (np == 1)
            break;
        }

      /***********************************************************************/
    }
  else if (type == WWD_MESH_TYPE_SSS)
    {
      /***********************************************************************/

      /***** Serpent wwd-format **********************************************/

      /* Check if rmx data already exists */

      if ((rmx = (long)RDB[DATA_PTR_RMX0]) > VALID_PTR)
        {
          /* Check for multiple */

          if (NextItem(rmx) > VALID_PTR)
            Error(rmx, "Multiple rmx definitions not allowed with wwd files");
        }
      else
        {
          /* Create response matrix structure */

          rmx = NewItem(DATA_PTR_RMX0, RMX_BLOCK_SIZE);
        }

      /* Put pointer */

      WDB[wwd + WWD_PTR_RMX] = (double)rmx;

      /* Put flag */

      WDB[rmx + RMX_FROM_FILE] = (double)YES;

      /* Open file for reading */

      if ((fp = fopen(fname, "r")) == NULL)
        Error(wwd, "Weight window data file \"%s\" does not exist", fname);

      /**************************/

      /* Väliaikainen viritys joka sallii vanhojen mesh-tyyppien lukemisen */

      if (1 != 2)
        {
          /* Read particle type */

          if ((sz = fread(&n, sizeof(long), 1, fp)) == 0)
            Error(wwd, "Error in weight window data file");
        }
      else
        n = PARTICLE_TYPE_GAMMA;

      /***************************/

      /* Put type */

      if ((n != PARTICLE_TYPE_NEUTRON) && (n != PARTICLE_TYPE_GAMMA))
        Error(wwd, "Error in particle type");
      else
        WDB[wwd + WWD_PARTICLE_TYPE] = (double)n;

      /* Read number of energy groups */

      if ((sz = fread(&ng, sizeof(long), 1, fp)) == 0)
        Error(wwd, "Error in weight window data file");

      /* Check */

      if ((ng > 0) && (ng < 1001))
        {
          /* Allocate memory for group boundaries */

          ptr = ReallocMem(DATA_ARRAY, ng + 1);

          /* Read boundaries */

          if ((sz = fread(&WDB[ptr], sizeof(double), ng + 1, fp)) == 0)
            Error(wwd, "Error in weight window data file");

          /* Make into energy grid */

          ptr = MakeEnergyGrid(ng + 1, 0, 0, -1, &RDB[ptr], EG_INTERP_MODE_LIN);

          /* Put pointers */

          WDB[wwd + WWD_PTR_ERG] = (double)ptr;
          WDB[rmx + RMX_PTR_EGRID] = (double)ptr;
        }
      else if (ng > 1000)
        Error(wwd, "Invalid number of energy groups %ld\n", ng);
      else
        ng = 1;

      /* Put number of groups */

      WDB[wwd + WWD_NE] = (double)ng;
      WDB[rmx + RMX_NG] = (double)ng;

      /* Read number of cells */

      if ((sz = fread(&nc, sizeof(long), 1, fp)) == 0)
        Error(wwd, "Error in weight window data file");

      /* Check for existing mesh */

      if ((long)RDB[rmx + RMX_PTR_MESH_DATA] > VALID_PTR)
        Error(rmx, "Mesh type must be set to -1 when wwin/wf is used");

      /* Check for type 4 */

      if ((long)RDB[rmx + RMX_MODE] == RMX_MODE_WWG)
        Note(rmx, "Type %ld mesh generation may not work with wwin/wf",
             RMX_MODE_WWG);

      for (n = 0; n < nc; n++)
        {
          /* Create structure */

          loc0 = NewItem(rmx + RMX_PTR_MESH_DATA, RMX_CELL_BLOCK_SIZE);

          /* Read cell index */

          if ((sz = fread(&i, sizeof(long), 1, fp)) == 0)
            Error(wwd, "Error in weight window data file");

          WDB[loc0 + RMX_CELL_IDX] = (double)(i);

          /* Read mesh index */

          if ((sz = fread(&i, sizeof(long), 1, fp)) == 0)
            Error(wwd, "Error in weight window data file");

          WDB[loc0 + RMX_CELL_MESH_IDX] = (double)(i);

          /* Read volume */

          if ((sz = fread(&WDB[loc0 + RMX_CELL_RVOL], sizeof(double), 1, fp))
              == 0)
            Error(wwd, "Error in weight window data file");

          /* Read maximum number of neighbors */

          if ((sz = fread(&nmax, sizeof(long), 1, fp)) == 0)
            Error(wwd, "Error in weight window data file");

          WDB[loc0 + RMX_CELL_MTX_SIZE] = (double)nmax;

          /* Loop over bounds */

          for (m = 0; m < nmax; m++)
            {
              /* Create structure */

              loc1 = NewItem(loc0 + RMX_CELL_PTR_BOUNDS,
                            RMX_CELL_BOUND_BLOCK_SIZE);

              /* Read cell index */

              if ((sz = fread(&i, sizeof(long), 1, fp)) == 0)
                Error(wwd, "Error in weight window data file");

              WDB[loc1 + RMX_CELL_BOUND_PTR_CELL] = (double)(i);

              /* Read forward and adjoint indexes */

              if ((sz = fread(&i, sizeof(long), 1, fp)) == 0)
                Error(wwd, "Error in weight window data file");

              WDB[loc1 + RMX_CELL_BOUND_FWD_IDX] = (double)(i);

              if ((sz = fread(&i, sizeof(long), 1, fp)) == 0)
                Error(wwd, "Error in weight window data file");

              WDB[loc1 + RMX_CELL_BOUND_ADJ_IDX] = (double)(i);
            }

          /* Read source importances */

          ptr = ReallocMem(DATA_ARRAY, ng);
          WDB[loc0 + RMX_CELL_IMP_SRC] = (double)ptr;

          if ((sz = fread(&WDB[ptr], sizeof(double), ng, fp)) == 0)
            Error(wwd, "Error in weight window data file");

          ptr = ReallocMem(DATA_ARRAY, ng);
          WDB[loc0 + RMX_CELL_IMP_SRC_KEEP] = (double)ptr;

          if ((sz = fread(&WDB[ptr], sizeof(double), ng, fp)) == 0)
            Error(wwd, "Error in weight window data file");

          /* Read current importances */

          ptr = ReallocMem(DATA_ARRAY, ng);
          WDB[loc0 + RMX_CELL_IMP_CURR] = (double)ptr;

          if ((sz = fread(&WDB[ptr], sizeof(double), ng, fp)) == 0)
            Error(wwd, "Error in weight window data file");

          ptr = ReallocMem(DATA_ARRAY, ng);
          WDB[loc0 + RMX_CELL_IMP_CURR_KEEP] = (double)ptr;

          if ((sz = fread(&WDB[ptr], sizeof(double), ng, fp)) == 0)
            Error(wwd, "Error in weight window data file");

          /* Read partial solutions */

          ptr = ReallocMem(DATA_ARRAY, nmax*ng);
          WDB[loc0 + RMX_CELL_ADJ_SOL_IN_CURR] = (double)ptr;

          if ((sz = fread(&WDB[ptr], sizeof(double), nmax*ng, fp)) == 0)
            Error(wwd, "Error in weight window data file");
        }

      /* Close list */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
      CloseList(loc0);

      /* Sort by index */

      SortList(loc0, RMX_CELL_IDX, SORT_MODE_ASCEND);

      /* Put pointers */

      WDB[wwd + WWD_PTR_MESH_DATA] = (double)loc0;
      WDB[rmx + RMX_PTR_MESH_DATA] = (double)loc0;

      /* Loop over data */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Loop over bounds */

          loc1 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
          while (loc1 > VALID_PTR)
            {
              /* Get cell index */

              i = (long)RDB[loc1 + RMX_CELL_BOUND_PTR_CELL];
              CheckValue(FUNCTION_NAME, "i", "", i, 0, 1000000000);

              /* Pointer to data */

              ptr = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

              /* Find match */

              if ((ptr = SeekList(ptr, RMX_CELL_IDX, (double)i,
                                  SORT_MODE_ASCEND)) < VALID_PTR)
                Die(FUNCTION_NAME, "Mesh cell not found");

              /* Put pointer */

              WDB[loc1 + RMX_CELL_BOUND_PTR_CELL] = (double)ptr;

              /* Next */

              loc1 = NextItem(loc1);
            }
          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Read mesh type */

      if ((sz = fread(&mty, sizeof(long), 1, fp)) == 0)
        Error(wwd, "Error in weight window data file");

      /* Avoid compiler warning */

      lims = NULL;

      /* Check for adaptive mesh */

      if (mty == MESH_TYPE_ADAPTIVE)
        {
          /*******************************************************************/

          /***** Adaptive type ***********************************************/

          /* Create mesh recursively */

          msh = ReadWWMesh0(fp, rmx, wwd);

          /* Put mesh pointers */

          WDB[wwd + WWD_PTR_MESH] = (double)msh;
          WDB[rmx + RMX_PTR_MESH] = (double)msh;

          /*******************************************************************/
        }
      else
        {
          /*******************************************************************/

          /***** Conventional types ******************************************/

          /* Check type */

          if ((mty == MESH_TYPE_CARTESIAN) || (mty == MESH_TYPE_CYLINDRICAL) ||
              (mty == MESH_TYPE_SPHERICAL) || (mty == MESH_TYPE_HEXX) ||
              (mty == MESH_TYPE_HEXY))
            {
              /* Regular square, cylindrical and hex types, read sizes */

              if ((sz = fread(&nx, sizeof(long), 1, fp)) == 0)
                Error(wwd, "Error in weight window data file");

              if ((sz = fread(&ny, sizeof(long), 1, fp)) == 0)
                Error(wwd, "Error in weight window data file");

              if ((sz = fread(&nz, sizeof(long), 1, fp)) == 0)
                Error(wwd, "Error in weight window data file");

              /* Allocate memory for boundaries */

              lims = (double *)Mem(MEM_ALLOC, 6, sizeof(double));

              /* Read boundaries */

              if ((sz = fread(lims, sizeof(double), 6, fp)) == 0)
                Error(wwd, "Error in weight window data file");
            }
          else if ((mty == MESH_TYPE_ORTHOGONAL) ||
                   (mty == MESH_TYPE_ICYL))
            {
              /* Unevenly-distributed types , read sizes */

              if ((sz = fread(&nx, sizeof(long), 1, fp)) == 0)
                Error(wwd, "Error in weight window data file");

              if ((sz = fread(&ny, sizeof(long), 1, fp)) == 0)
                Error(wwd, "Error in weight window data file");

              if ((sz = fread(&nz, sizeof(long), 1, fp)) == 0)
                Error(wwd, "Error in weight window data file");

              /* Allocate memory */

              lims = (double *)Mem(MEM_ALLOC, nx + ny + nz + 3, sizeof(double));

              /* Read data */

              if ((sz = fread(lims, sizeof(double), nx + ny + nz + 3, fp)) == 0)
                Error(wwd, "Error in weight window data file");
            }
          else
            Error(wwd, "Mesh type %ld not supported", mty);

          /* Check sizes */

          if ((nx < 1) || (nx > 1000000))
            Error(wwd, "Invalid mesh size %ld", nx);
          else if ((ny < 1) || (ny > 1000000))
            Error(wwd, "Invalid mesh size %ld", ny);
          else if ((nz < 1) || (nz > 1000000))
            Error(wwd, "Invalid mesh size %ld", nz);

          /* Check mesh-specific */

          if (mty == MESH_TYPE_CARTESIAN)
            {
              if (lims[1] <= lims[0])
                Error(wwd, "Invalid limiting values %E and %E\n",
                      lims[0], lims[1]);
              else if (lims[3] <= lims[2])
                Error(wwd, "Invalid limiting values %E and %E\n",
                      lims[2], lims[3]);
              else if (lims[5] <= lims[4])
                Error(wwd, "Invalid limiting values %E and %E\n",
                      lims[4], lims[5]);
            }
          else if (mty == MESH_TYPE_SPHERICAL)
            {
              if (lims[0] < 0.0)
                Error(wwd, "Invalid radius %E\n", lims[0]);
              else if (lims[1] <= lims[0])
                Error(wwd, "Invalid limiting values %E and %E\n",
                      lims[0], lims[1]);
              else if (lims[3] <= lims[2])
                Error(wwd, "Invalid limiting values %E and %E\n",
                      lims[2], lims[3]);
              else if (lims[5] <= lims[4])
                Error(wwd, "Invalid limiting values %E and %E\n",
                      lims[4], lims[5]);
            }
          else if (mty == MESH_TYPE_CYLINDRICAL)
            {
              if (lims[0] < 0.0)
                Error(wwd, "Invalid radius %E\n", lims[0]);
              else if (lims[1] <= lims[0])
                Error(wwd, "Invalid limiting values %E and %E\n",
                      lims[0], lims[1]);
              else if (lims[3] <= lims[2])
                Error(wwd, "Invalid limiting values %E and %E\n",
                      lims[2], lims[3]);
              else if (lims[5] <= lims[4])
                Error(wwd, "Invalid limiting values %E and %E\n",
                      lims[4], lims[5]);
            }

          /* Create mesh */

          msh = CreateMesh(mty, MESH_CONTENT_PTR, -1, nx, ny, nz, lims, -1);

          /* Put pointers */

          WDB[wwd + WWD_PTR_MESH] = (double)msh;
          WDB[rmx + RMX_PTR_MESH] = (double)msh;

          /* Create structures */

          idx = 0;

          for (k = 0; k < nz; k++)
            for (j = 0; j < ny; j++)
              for (i = 0; i < nx; i++)
                {
                  /* Pointer to data */

                  loc0 = ReadMeshPtr(msh, i, j, k);
                  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                  /* Read index */

                  if ((sz = fread(&idx, sizeof(long), 1, fp)) == 0)
                    Error(wwd, "Error in weight window data file");

                  /* Pointer to data */

                  loc1 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
                  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

                  /* Find match */

                  if ((loc1 = SeekList(loc1, RMX_CELL_IDX, (double)idx,
                                       SORT_MODE_ASCEND)) < VALID_PTR)
                    Die(FUNCTION_NAME, "Mesh cell not found");

                  /* Put pointer */

                  WDB[loc0] = (double)loc1;

                  /* Put mesh pointer */

                  WDB[loc1 + RMX_CELL_PTR_MESH] = (double)msh;
                }

          /* Put data pointer */

          WDB[wwd + WWD_PTR_MESH_DATA] = WDB[rmx + RMX_PTR_MESH_DATA];

          /* Free temporary array */

          Mem(MEM_FREE, lims);

          /*******************************************************************/
        }

      /***********************************************************************/
    }
  else
    Error(wwd, "File type %ld not supported", type);

  /* Exit OK */

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/

/***** Recursive read function ***********************************************/

long ReadWWMesh0(FILE *fp, long rmx, long wwd)
{
  long sz, nx, ny, nz, i, j, k, idx, loc0, loc1, ptr, msh;
  double lims[6];

  /* read size */

  if ((sz = fread(&nx, sizeof(long), 1, fp)) == 0)
    Error(wwd, "Error in weight window data file");

  if ((sz = fread(&ny, sizeof(long), 1, fp)) == 0)
    Error(wwd, "Error in weight window data file");

  if ((sz = fread(&nz, sizeof(long), 1, fp)) == 0)
    Error(wwd, "Error in weight window data file");

  /* Read boundaries */

  if ((sz = fread(&lims, sizeof(double), 6, fp)) == 0)
    Error(wwd, "Error in weight window data file");

  /* Create mesh */

  msh = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_PTR, -1, nx, ny, nz,
                   lims, -1);

  /* Loop over mesh */

  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
        {
          /* Pointer to data */

          loc0 = ReadMeshPtr(msh, i, j, k);
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* Read index */

          if ((sz = fread(&idx, sizeof(long), 1, fp)) == 0)
            Error(wwd, "Error in weight window data file");

          /* Check content */

          if (idx > -1)
            {
              /* Pointer to data */

              loc1 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

              /* Find match */

              if ((loc1 = SeekList(loc1, RMX_CELL_IDX, (double)idx,
                                   SORT_MODE_ASCEND)) < VALID_PTR)
                Die(FUNCTION_NAME, "Mesh cell not found");

              /* Put pointer */

              WDB[loc0] = (double)loc1;

              /* Put mesh pointer */

              WDB[loc1 + RMX_CELL_PTR_MESH] = (double)msh;
            }
          else if (idx == -1)
            {
              /* Another mesh follows, call recursively */

              ptr = ReadWWMesh0(fp, rmx, wwd);
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Put pointer */

              WDB[loc0] = -(double)ptr;
            }
          else
            Die(FUNCTION_NAME, "Invalid index read");
        }

  /* Change type to adaptive */

  WDB[msh + MESH_TYPE] = (double)MESH_TYPE_ADAPTIVE;

  /* Return mesh pointer */

  return msh;
}

/*****************************************************************************/
