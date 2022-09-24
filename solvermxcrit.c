/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : solvermx.c                                     */
/*                                                                           */
/* Created:       2017/10/08 (JLe)                                           */
/* Last modified: 2019/04/04 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Solves criticality eigenvalue problem using the response     */
/*              matrix method.                                               */
/*                                                                           */
/* Comments: - Assumes single detector and energy group                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SolveRMXCrit:"

/*****************************************************************************/

void SolveRMXCrit(long rmx)
{
  long loc0, ptr, det, loc1, loc2, sz, i, j, k, it, imax, nmax, compl;
  long l, m, n, nx, ny, nz, msh, fail;
  double flim, sum, val, *Jin, *Jout, Stot, frac, avg;
  double div, *S, bala;
  const double *s, *r, *rs, *alpha;
  char fname[MAX_STR];
  FILE *fp;

  /* Check that mesh exists */

  if ((loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA]) < VALID_PTR)
    Die(FUNCTION_NAME, "no mesh data");

  /* Get number of mesh cells */

  sz = ListSize(loc0);
  CheckValue(FUNCTION_NAME, "sz", "", sz, 2, 1000000000);

  /***************************************************************************/

  /***** Iteration loop ******************************************************/

  /* Put initial k-eff */

  if ((long)RDB[DATA_RMX_CONVG_ITER_IDX] == 0)
    WDB[rmx + RMX_CONVG_PREV_KEFF] = 1.0;

  /* Set convergence limit */

  flim = RDB[rmx + RMX_CONV];

  /* Set maximum number of iterations */

  imax = (long)RDB[rmx + RMX_N_ITER];
  CheckValue(FUNCTION_NAME, "imax", "", imax, -1, 10000);

  /* Reset completed flag */

  compl = 0;

  /* Loop over mesh */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to detector */

      det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
      CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

      /* Reset previous data */

      WDB[det + RMX_DET_FWD_PREV] = 0.0;

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Eigevalue iterations ************************************************/

  for (l = 0; l < 1000; l++)
    {
      /***********************************************************************/

      /***** Renormalize source **********************************************/

      /* Reset total source */

      Stot = 0.0;

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get total source rate */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
          S = &RES2[ptr];

          /* Add to total source */

          Stot = Stot + *S;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Check */

      if (Stot == 0.0)
        Die(FUNCTION_NAME, "Zero source");

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get total source rate */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
          S = &RES2[ptr];

          /* Divide by total source */

          *S = *S/Stot;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /***** Source to forward currents **************************************/

      /* Reset current balance */

      bala = 0.0;

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get matrix size */

          nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

          /* Get source coefficient vector */

          ptr = (long)RDB[loc0 + RMX_CELL_COEF_FWD_SRCC];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          s = &RDB[ptr];

          /* Get total source rate */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
          S = &RES2[ptr];

          /* Get outward current vector */

          ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          Jout = &WDB[ptr];

          /* Reset */

          memset(Jout, 0.0, nmax*sizeof(double));

          /* Copy values */

          for (i = 0; i < nmax; i++)
            {
              /* Put source particles to outward currents */

              Jout[i] = Jout[i] + s[i]*(*S);

              /* Add to current balance */

              bala = bala + s[i]*(*S);
            }

          /* Get inward current vector */

          ptr = (long)RDB[loc0 + RMX_CELL_WRK_IN_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          Jin = &WDB[ptr];

          /* Reset */

          memset(Jin, 0.0, nmax*sizeof(double));

          /* Reset solution vector */

          ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          memset(&WDB[ptr], 0.0, nmax*sizeof(double));

          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /* Main loop */

      for (it = 0; it != imax; it++)
        {
          /*******************************************************************/

          /***** Distribute forward current to neighbours ********************/

#ifdef OPEN_MP
#pragma omp parallel private(loc0, loc1, loc2, ptr, i, j, k, nmax, Jin, Jout)
#endif
          {

#ifdef OPEN_MP
#pragma omp for
#endif
            /* Loop over mesh */

            for (k = 0; k < sz; k++)
              {
                /* Pointer to data */

                loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                loc0 = ListPtr(loc0, k);
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                /* Check empty flag */

                if ((long)RDB[loc0 + RMX_CELL_EMPTY] == YES)
                  continue;

                /* Get matrix size */

                nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
                CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

                /* Get outward current vector */

                ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                Jout = &WDB[ptr];

                /* Loop over neighbours */

                loc2 = (long)RDB[loc0 + RMX_CELL_PTR_BOUNDS];
                for (i = 0; i < nmax; i++)
                  {
                    /* Pointer to cell */

                    loc1 = (long)RDB[loc2 + RMX_CELL_BOUND_PTR_CELL];
                    CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                    /* Check pointers */

                    if (loc0 == loc1)
                      Die(FUNCTION_NAME, "Pointer to self");

                    /* Boundary index */

                    j = (long)RDB[loc2 + RMX_CELL_BOUND_FWD_IDX];

                    /* Get inward current vector */

                    ptr = (long)RDB[loc1 + RMX_CELL_WRK_IN_CURR];
                    CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                    Jin = &WDB[ptr];

                    /* Put data */

                    Jin[j] += Jout[i];

                    /* Next */

                    loc2 = NextItem(loc2);
                  }
              }
          }

          /*******************************************************************/

          /***** Move forward currents from inward to outward buffer *********/

#ifdef OPEN_MP
#pragma omp parallel private(loc0, loc1, ptr, i, j, k, nmax, Jin, Jout, alpha)
#endif
          {

#ifdef OPEN_MP
#pragma omp for
#endif
            /* Loop over mesh */

            for (k = 0; k < sz; k++)
              {
                /* Pointer to data */

                loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                loc0 = ListPtr(loc0, k);
                CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

                /* Check empty flag */

                if ((long)RDB[loc0 + RMX_CELL_EMPTY] == YES)
                  continue;

                /* Get matrix size */

                nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
                CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

                /* Get inward current vector */

                ptr = (long)RDB[loc0 + RMX_CELL_WRK_IN_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                Jin = &WDB[ptr];

                /* Add inward current to solution */

                ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

                for (i = 0; i < nmax; i++)
                  WDB[ptr + i] = RDB[ptr + i] + Jin[i];

                /* Get outward current vector */

                ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
                CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
                Jout = &WDB[ptr];

                /* Reset outward current */

                memset(Jout, 0.0, nmax*sizeof(double));

                /* Get coefficient matrix */

                ptr = (long)RDB[loc0 + RMX_CELL_COEF_FWD_TRANS_MTX];
                CheckPointer(FUNCTION_NAME, "(ptr4)", DATA_ARRAY, ptr);
                alpha = &RDB[ptr];

                /* Loop over matrix */

                  for (i = 0; i < nmax; i++)
                    if (Jin[i] != 0.0)
                      for (j = 0; j < nmax; j++)
                        Jout[j] = Jout[j] + alpha[i*nmax + j]*Jin[i];

                /* Reset inward current */

                memset(Jin, 0.0, nmax*sizeof(double));
              }
          }

          /*******************************************************************/

          /***** Check convergence *******************************************/

          /* Reset sum of inward currents and current fraction */

          sum = 0.0;
          frac = INFTY;

          /* Loop over mesh */

          loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          while (loc0 > VALID_PTR)
            {
              /* Get matrix size */

              nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
              CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

              /* Get outward current vector */

              ptr = (long)RDB[loc0 + RMX_CELL_WRK_OUT_CURR];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              Jout = &WDB[ptr];

              /* Add to sum */

              for (i = 0; i < nmax; i++)
                sum = sum + Jout[i];

              /* Next cell */

              loc0 = NextItem(loc0);
            }

          /* Calculate remaining fraction */

          if (bala == 0.0)
            Die(FUNCTION_NAME, "Error in current balance");
          else
            frac = sum/bala;

          /* Check convergence */

          if ((it > 10) && (frac < flim))
            break;

          /*******************************************************************/
        }

      /***********************************************************************/

      /***** Calculate source distribution for next cycle ********************/

      /* Loop over mesh */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Get matrix size */

          nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
          CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

          /* Get total source rate */

          ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
          CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
          S = &RES2[ptr];

          /* Pointer to detector */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

          /* Get direct source-to-response coefficient */

          ptr = (long)RDB[det + RMX_DET_COEF_DIR_SRCC_RES];
          CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);
          rs = &RDB[ptr];

          /* Put direct contribution */

          WDB[det + RMX_DET_FWD_CHK] = (*rs)*(*S);

          /* Get response vector */

          ptr = (long)RDB[det + RMX_DET_COEF_FWD_RES];
          CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);
          r = &RDB[ptr];

          /* Get inward current solution */

          ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
          Jin = &WDB[ptr];

          /* Add contribution */

            for (i = 0; i < nmax; i++)
              WDB[det + RMX_DET_FWD_CHK] = RDB[det + RMX_DET_FWD_CHK]
                + r[i]*Jin[i];

          /* Put value to new source (NOTE: olettaa että */
          /* energiaryhmiä on vain 1) */

          S[0] = RDB[det + RMX_DET_FWD_CHK];

          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /***** Check convergence ***********************************************/

      /* Calculate sums */

      val = 0.0;
      div = 0.0;

      /* Loop over data */

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Pointer to detector */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

          /* Check values */

          if (RDB[det + RMX_DET_FWD_CHK] > 0.0)
            {
              /* Add to sum */

              div = div + RDB[det + RMX_DET_FWD_CHK];
              val = val + 1.0;
            }

          /* Pointer to next */

          loc0 = NextItem(loc0);
        }

      /* Check */

      if (val == 0.0)
        Die(FUNCTION_NAME, "Error in sum");

      /* Count number of failed cells */

      fail = 0;

      loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
      while (loc0 > VALID_PTR)
        {
          /* Pointer to detector */

          det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
          CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

          /* Compare to previous */

          if (RDB[det + RMX_DET_FWD_PREV] > 0.0)
            if (fabs(val*RDB[det + RMX_DET_FWD_CHK]/div/
                     RDB[det + RMX_DET_FWD_PREV] - 1.0) >
                RDB[DATA_RMX_CONVG_OUT_TOL])
              fail++;

          /* Remember value */

          WDB[det + RMX_DET_FWD_PREV] = val*RDB[det + RMX_DET_FWD_CHK]/div;

          /* Divide by total to get probability */

          WDB[loc0 + RMX_CELL_SRC_PROB] = RDB[det + RMX_DET_FWD_CHK]/div;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /***********************************************************************/

      /***** Print distribution **********************************************/

      if ((1 == 2) && (mpiid == 0))
        {
          /* Sort list */

          loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
          SortList(loc0, RMX_CELL_IDX, SORT_MODE_ASCEND);

          /* Open file */

          sprintf(fname, "%s_rmx_src.m", GetText(DATA_PTR_INPUT_FNAME));

          if (((long)RDB[DATA_RMX_CONVG_ITER_IDX] == 0.0) && (l == 0))
            fp = fopen(fname, "w");
          else
            fp = fopen(fname, "a");

          /* Loop over data and print */

          fprintf(fp, "x%ld(:,%ld) = [\n", (long)RDB[DATA_RMX_CONVG_ITER_IDX],
                  l + 1);

          loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
          while (loc0 > VALID_PTR)
            {
              /* Get matrix size */

              nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
              CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

              /* Get total source rate */

              ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
              CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
              S = &RES2[ptr];

              /* Local mesh size */

              nx = (long)RDB[rmx + RMX_CONVG_SUBMESH_NX];
              ny = (long)RDB[rmx + RMX_CONVG_SUBMESH_NY];
              nz = (long)RDB[rmx + RMX_CONVG_SUBMESH_NZ];

              /* Get inward current solution */

              ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
              CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
              Jin = &WDB[ptr];

              m = 0;

              for (k = 0; k < nz; k++)
                for (j = 0; j < ny; j++)
                  for (i = 0; i < nx; i++)
                    {
                      ptr = (long)RDB[loc0 + RMX_CELL_COEF_SRC_FF];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      /* Divide by k-eff */

                      if (l == 0)
                        val = (*S)*RDB[ptr + m]/
                          RDB[rmx + RMX_CONVG_PREV_KEFF];
                      else
                        val = (*S)*RDB[ptr + m]/Stot;

                      ptr = (long)RDB[loc0 + RMX_CELL_COEF_SURF_FF];
                      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                      for (n = 0; n < nmax; n++)
                        if (Jin[n] != 0.0)
                          val = val + Jin[n]*RDB[ptr + n*nx*ny*nz + m];

                      /* Print */

                      fprintf(fp, "%1.5E\n", val);

                      /* Add to index */

                      m++;
                    }


              /*
              det = (long)RDB[loc0 + RMX_CELL_PTR_DET];
              fprintf(fp, "%1.5f\n", val*RDB[det + RMX_DET_FWD_CHK]/div);
              */
              loc0 = NextItem(loc0);
            }

          fprintf(fp, "];\n");
          fclose(fp);
        }

      /***********************************************************************/

      /* Check convergence */

      if ((l > 10) && (fail == 0))
        break;

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Calculate average ***************************************************/

  /* NOTE: Tätä tarvitaan jakauman normeeraamiseen UFS-moodissa */

  /* Reset variables */

  avg = 0.0;
  div = 0.0;

  /* Local mesh size */

  nx = (long)RDB[rmx + RMX_CONVG_SUBMESH_NX];
  ny = (long)RDB[rmx + RMX_CONVG_SUBMESH_NY];
  nz = (long)RDB[rmx + RMX_CONVG_SUBMESH_NZ];

  /* Loop over cells */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Get matrix size */

      nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Get total source rate */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
      S = &RES2[ptr];

      /* Pointer to mesh */

      msh = (long)RDB[loc0 + RMX_CELL_PTR_SUBMESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get inward current solution */

      ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
      Jin = &WDB[ptr];

      m = 0;

      /* Source contributions */

      ptr = (long)RDB[loc0 + RMX_CELL_COEF_SRC_FF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      rs = &WDB[ptr];

      /* Surface contributions */

      ptr = (long)RDB[loc0 + RMX_CELL_COEF_SURF_FF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      s = &WDB[ptr];

      /* Loop over distribution */

      for (k = 0; k < nz; k++)
        for (j = 0; j < ny; j++)
          for (i = 0; i < nx; i++)
            {
              /* Source */

              if (l == 0)
                val = (*S)*rs[m]/RDB[rmx + RMX_CONVG_PREV_KEFF];
              else
                val = (*S)*rs[m]/Stot;

              /* Surface contributions */

              for (n = 0; n < nmax; n++)
                if (Jin[n] != 0.0)
                  val = val + Jin[n]*s[n*nx*ny*nz + m];

              /* Add to sum and number of non-zeros */

              avg = avg + val;

              if (val > 0.0)
                div = div + 1.0;

              /* Add to index */

              m++;
            }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Check */

  if ((div == 0.0) || (avg == 0.0))
    Die(FUNCTION_NAME, "No non-zero zones");
  else
    avg = avg/div;

  /***************************************************************************/

  /***** Calculate local distributions ***************************************/

  /* Local mesh size */

  nx = (long)RDB[rmx + RMX_CONVG_SUBMESH_NX];
  ny = (long)RDB[rmx + RMX_CONVG_SUBMESH_NY];
  nz = (long)RDB[rmx + RMX_CONVG_SUBMESH_NZ];

  /* Loop over cells */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Get matrix size */

      nmax = (long)RDB[loc0 + RMX_CELL_MTX_SIZE];
      CheckValue(FUNCTION_NAME, "nmax", "", nmax, 1, 10000);

      /* Get total source rate */

      ptr = (long)RDB[loc0 + RMX_CELL_MC_SRCC_TOT];
      CheckPointer(FUNCTION_NAME, "ptr", RES2_ARRAY, ptr);
      S = &RES2[ptr];

      /* Pointer to mesh */

      msh = (long)RDB[loc0 + RMX_CELL_PTR_SUBMESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get inward current solution */

      ptr = (long)RDB[loc0 + RMX_CELL_FWD_SOL_IN_CURR];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
      Jin = &WDB[ptr];

      m = 0;

      /* Source contributions */

      ptr = (long)RDB[loc0 + RMX_CELL_COEF_SRC_FF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      rs = &WDB[ptr];

      /* Surface contributions */

      ptr = (long)RDB[loc0 + RMX_CELL_COEF_SURF_FF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      s = &WDB[ptr];

      /* Loop over distribution */

      for (k = 0; k < nz; k++)
        for (j = 0; j < ny; j++)
          for (i = 0; i < nx; i++)
            {
              /* Source */

              if (l == 0)
                val = (*S)*rs[m]/RDB[rmx + RMX_CONVG_PREV_KEFF];
              else
                val = (*S)*rs[m]/Stot;

              /* Surface contributions */

              for (n = 0; n < nmax; n++)
                if (Jin[n] != 0.0)
                  val = val + Jin[n]*s[n*nx*ny*nz + m];

              /* Put value */

              if ((long)RDB[DATA_UFS_MODE] == UFS_MODE_RMX)
                PutMeshIdx(msh, val/avg, i, j, k);
              else
                PutMeshIdx(msh, val, i, j, k);

              /* Add to index */

              m++;
            }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Create sampling lists ***********************************************/

  /* Loop over cells */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to sub-mesh */

      msh = (long)RDB[loc0 + RMX_CELL_PTR_SUBMESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get total */

      sum = MeshTot(msh);

      /* Get pointer to data */

      loc1 = (long)RDB[msh + MESH_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Get pointer to sample array */

      ptr = (long)RDB[loc0 + RMX_CELL_PTR_SUBMESH_SAMPLE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over data (sum is stored to index 0) */

      for (i = 0; i < nx*ny*nz; i++)
        {
          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Check values */

          CheckValue(FUNCTION_NAME, "p", "", RDB[loc1 + i + 1], 0.0, 1E+18);
          CheckValue(FUNCTION_NAME, "sum", "", sum, 0.0, 1E+18);

          /* Put index */

          WDB[ptr + RMX_SUBMESH_SAMPLE_IDX] = (double)i;

          /* Check sum and put probability */

          if (sum > 0.0)
            WDB[ptr + RMX_SUBMESH_SAMPLE_P] = RDB[loc1 + i + 1]/sum;
          else
            WDB[ptr + RMX_SUBMESH_SAMPLE_P] = 0.0;

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Sort list */

      ptr = (long)RDB[loc0 + RMX_CELL_PTR_SUBMESH_SAMPLE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      SortList(ptr, RMX_SUBMESH_SAMPLE_P, SORT_MODE_DESCEND);

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Sort list for sampling */

  loc0 = (long)RDB[rmx + RMX_PTR_MESH_DATA];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  SortList(loc0, RMX_CELL_SRC_PROB, SORT_MODE_DESCEND);

  /***************************************************************************/

  /***** Finalize ************************************************************/

  /* Remember previous k-eff */

  WDB[rmx + RMX_CONVG_PREV_KEFF] = Stot;

  /* Use value also in neutronics calculation */

  ptr = (long)RDB[DATA_PTR_CYCLE_EIG_KEFF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  for (n = 0; n < (long)RDB[DATA_N_POP_EIG]; n++)
    WDB[ptr + n] = Stot;

  /***************************************************************************/
}

/*****************************************************************************/
