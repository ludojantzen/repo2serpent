/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : homoflux.c                                     */
/*                                                                           */
/* Created:       2014/05/15 (JLe)                                           */
/* Last modified: 2019/12/18 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates homogeneous diffusion flux for ADF's and          */
/*              pin-power form functions.                                    */
/*                                                                           */
/* Comments: - Varmista että kulma-alueen koko on sama kuin dfpos.c:ssä.     */
/*                                                                           */
/*           - Ne negatiiviset dfsol():in arvot pitää selvittää              */
/*                                                                           */
/*           - BEAVRS testikeisseissä HZP_B150 antaa päättömiä tuloksia      */
/*             (sama dfsolver.m -skriptillä)                                 */
/*                                                                           */
/*           - Tota hexajuttua ei oo kunnolla testattu                       */
/*                                                                           */
/*           - PPW-juttu ei ehkä toimi jos ADF ei oo määritetty              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "HomoFlux:"

/* Local test function prototypes */

void CalcTotFlx(long surf, long nG, long nj, double **surfs_n,
                complex **T, complex **U, complex *c);

void PrintFlx(long gcu, long surf, long nG, long nj, double **surfs_n,
              complex **T, complex **U, complex *c);

/* Surfaces */

#define W 0
#define S 1
#define E 2
#define N 3

/* Corners */

#define NW 0
#define NE 1
#define SE 2
#define SW 3

/*****************************************************************************/

void HomoFlux(long gcu)
{
  long nG, adf, ptr, surf, ppw, type, ns, nc, np, n, m, i, j, k, Nr, ifCorn;
  long maps[6], mapc[6], nneg, ptr2, Nr2;
  double *vflx, *sflx, *cflx, *pflx, *D, *tot, *nsf, *chi, *scatt, *Js, *Jm;
  double *Jc, *L, **surfs_n, **surfs_r, *phi, *t, r[3], tn[3], div, sum;
  double vol, src, loss, keff, dc, nx, ny, *sflxsgn, sgns[6];
  const double *ds, *dm;
  complex z;
  complex **A,**T,**U,**F, *b, *c;

  /***************************************************************************/

  /***** Common stuff ********************************************************/

  /* Check pointer to group constant universe */

  CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);

  /* Get pointer to ADF structure */

  if ((adf = (long)RDB[gcu + GCU_PTR_ADF]) < VALID_PTR)
    return;

  /* Get number of surfaces and corners */

  ns = (long)RDB[adf + ADF_NSURF];
  CheckValue(FUNCTION_NAME, "ns", "", ns, 1, 6);

  nc = (long)RDB[adf + ADF_NCORN];
  CheckValue(FUNCTION_NAME, "nc", "", nc, 0, 6);

  /* Get number of energy groups */

  nG = (long)RDB[DATA_ERG_FG_NG];
  CheckValue(FUNCTION_NAME, "nG", "", nG, 1, 100000);

  /* Get areas of surface and mid-point segments (lengths in 2D calculation) */

  ptr = (long)RDB[adf + ADF_PTR_SURF_AREA];
  ds = &RDB[ptr];

  ptr = (long)RDB[adf + ADF_PTR_MID_AREA];
  dm = &RDB[ptr];

  /* Use single width for corner segments (check that all values are equal) */

  ptr = (long)RDB[adf + ADF_PTR_CORN_AREA];

  for (n = 1; n < nc; n++)
    {
      /* Check zero */

      if (RDB[ptr + n - 1] < ZERO)
        Die(FUNCTION_NAME, "Zero corner segment size");

      /* Check if differ */

      if (fabs(RDB[ptr + n]/RDB[ptr + n - 1] - 1.0) > 1E-6)
        Die(FUNCTION_NAME, "Mismatch in corner segment size");
    }

  dc = RDB[ptr];

  /* Get volume (cross-sectional area in 2D calculation) */

  vol = RDB[adf + ADF_VOL];
  CheckValue(FUNCTION_NAME, "vol", "", vol, ZERO, INFTY);

  /* Check pointer to pin-power distribution */

  if ((ppw = (long)RDB[gcu + GCU_PTR_PPW]) > VALID_PTR)
    {
      /* Number of pins */

      np = (long)RDB[ppw + PPW_NP];
      CheckValue(FUNCTION_NAME, "np", "", np, 1, 100000);
    }
  else
    np = 0;

  /* Check mode */

  if ((long)RDB[DATA_HOMOFLUX_SOLVER] == 1)
    {
      /* Include corners */

      ifCorn = 1;
    }
  else
    {
      /* No corners */

      ifCorn = 0;
    }

  /* Set number of poins for trapezoidal integration */

  Nr = (long)RDB[DATA_ADF_TRAPZ_PT];

  if (Nr % 2)
      Nr2 = Nr/2;
  else
      Nr2 = (Nr + 1)/2;

  /* Reset number of corners if not included */

  if (ifCorn == 0)
    nc = 0;

  /* Reset negative point counter */

  nneg = 0;

  /* Sign moments of discontinuity factors */

  sgns[0] = 0.0;
  sgns[1] = 0.0;
  sgns[2] = 0.0;
  sgns[3] = 0.0;
  sgns[4] = 0.0;
  sgns[5] = 0.0;

  /***************************************************************************/

  /***** Check for zero boundary currents ************************************/

  /* Reset number of non-zero currents */

  i = 0;

  /* Loop over energy groups */

  for (n = 0; n < nG; n++)
    {
      /* Get volume flux */

      ptr = (long)RDB[gcu + GCU_RES_FG_DF_HET_VOL_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);
      div = BufVal(ptr, n);

      /* Reset sum */

      sum = 0.0;

      /* Get pointer to net current over surfaces */

      ptr = (long)RDB[gcu + GCU_RES_FG_DF_SURF_NET_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

      /* Loop over boundaries and add to sum and count */

      for (m = 0; m < ns; m++)
        if (BufVal(ptr, m, n) != 0.0)
          {
            sum = sum + BufVal(ptr, m, n);
            i++;
          }

      /* Check */

      if (div == 0.0)
        {
          /* Zero flux, skip calculation */

          i = 0;

          /* Break outer loop */

          break;
        }
      else if (fabs(sum/div) < 1E-6)
        {
          /* Close to zero total net current, skip calculation */

          i = 0;

          /* Break outer loop */

          break;
        }
    }

  /* Get pointer to outer boundary surface */

  surf = (long)RDB[adf + ADF_PTR_SURF];
  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Get surface type */

  type = (long)RDB[surf + SURFACE_TYPE];

  /* Check */

  if ((type == SURF_CUBE) || (type == SURF_CUBOID) ||
      (type == SURF_HEXXPRISM) || (type == SURF_HEXYPRISM))
    {
      /* Print warning */

      Note(0, "Unsupported surface type, skipping diffusion flux solver");

      /* Skip calculation */

      i = 0;
    }

  /* Check count */

  if ((i == 0) || ((long)RDB[adf + ADF_ENFORCE_REF] == YES))
    {
      /* Allocate memory for volume flux */

      vflx = (double *)Mem(MEM_ALLOC, nG, sizeof(double));

      /* Read data */

      ptr = (long)RDB[gcu + GCU_RES_FG_DF_HET_VOL_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr);

      for (n = 0; n < nG; n++)
        vflx[n] = BufVal(ptr, n);

      /* Mark scoring buffer unreduced */

      WDB[DATA_BUF_REDUCED] = (double)NO;

      /* Store data */

      for (n = 0; n < nG; n++)
        {
          /* Copy value to surface fluxes */

          ptr = (long)RDB[gcu + GCU_RES_FG_DF_HOM_SURF_FLUX];
          CheckPointer(FUNCTION_NAME, "(ptr3)", DATA_ARRAY, ptr);

          for (m = 0; m < ns; m++)
            AddBuf(vflx[n], 1.0, ptr, 0, -1, m, n);

          /* Copy value to corner fluxes */

          if (nc > 0)
            {
              ptr = (long)RDB[gcu + GCU_RES_FG_DF_HOM_CORN_FLUX];
              CheckPointer(FUNCTION_NAME, "(ptr4)", DATA_ARRAY, ptr);

              for (m = 0; m < nc; m++)
                AddBuf(vflx[n], 1.0, ptr, 0, -1, m, n);
            }

          /* Homogeneous flux for pin-power form factors */

          if (np > 0)
            {
              /* Get pointer to flux distribution */

              ptr = (long)RDB[gcu + GCU_RES_FG_PPW_HOM_FLUX];
              CheckPointer(FUNCTION_NAME, "(ptr5)", DATA_ARRAY, ptr);

              /* Loop over pins and store value */

              for (m = 0; m < np; m++)
                AddBuf(vflx[n]/vol, 1.0, ptr, 0, -1, m, n);
            }
        }

      /* Reduce buffer */

      ReduceBuffer();

      /* Free memory */

      Mem(MEM_FREE, vflx);

      /* Exit subroutine */

      return;
    }

  /***************************************************************************/

  /***** Remap currents ******************************************************/

  /* NOTE: Tätä tarvitaan siitä syystä että ScoreDF() ym. indeksoi */
  /* pinnat ja kulmat eri tavalla kuin mitä täällä on käytetty. */

  /* Get pointer to outer boundary surface */

  surf = (long)RDB[adf + ADF_PTR_SURF];
  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Get surface type */

  type = (long)RDB[surf + SURFACE_TYPE];

  /* Remap currents */

  if ((type == SURF_SQC) || (type == SURF_RECT))
    {
      /* Square geometry: W-S-E-N --> E-N-W-S, NW-NE-SE-SW --> NE-NW-SW-SE */

      maps[W] = E;
      maps[S] = N;
      maps[E] = W;
      maps[N] = S;
      maps[4] = 0;
      maps[5] = 0;

      mapc[NW] = NE;
      mapc[NE] = NW;
      mapc[SE] = SW;
      mapc[SW] = SE;
      mapc[4] = 0;
      mapc[5] = 0;

      sgns[0] = -1.0;
      sgns[1] = 1.0;
      sgns[2] = 1.0;
      sgns[3] = -1.0;
    }
  else if ((type == SURF_CUBE) || (type == SURF_CUBOID))
    {
      /* Square geometry: W-S-E-N --> E-N-W-S, NW-NE-SE-SW --> NE-NW-SW-SE */

      maps[W] = E;
      maps[S] = N;
      maps[E] = W;
      maps[N] = S;
      maps[4] = 4;
      maps[5] = 5;

      mapc[NW] = NE;
      mapc[NE] = NW;
      mapc[SE] = SW;
      mapc[SW] = SE;
      mapc[4] = 0;
      mapc[5] = 0;

      sgns[0] = 1.0;
      sgns[1] = -1.0;
      sgns[2] = -1.0;
      sgns[3] = 1.0;
    }
  else if ((type == SURF_HEXXC) || (type == SURF_HEXYC))
    {
      /* Hex geometry */

      maps[0] = 1;
      maps[1] = 0;
      maps[2] = 5;
      maps[3] = 4;
      maps[4] = 3;
      maps[5] = 2;

      mapc[0] = 5;
      mapc[1] = 0;
      mapc[2] = 1;
      mapc[3] = 2;
      mapc[4] = 3;
      mapc[5] = 4;

      if (type == SURF_HEXXC)
        {
          /* The following comments are hopefully correct */

          /* Also the sgns are still to be confirmed */

          sgns[0] = 1.0; /* E  in Serpent, NE in homoflux */
          sgns[1] = 1.0; /* NE in Serpent, N  in homoflux */
          sgns[2] = 1.0; /* NW in Serpent, NW in homoflux */
          sgns[3] = -1.0; /* W  in Serpent, SE in homoflux */
          sgns[4] = -1.0; /* SW in Serpent, S  in homoflux */
          sgns[5] = -1.0; /* SE in Serpent, SW in homoflux */
        }
    }
  else
    {
      /* 1D geometries, etc. */

      maps[0] = 0;
      maps[1] = 1;
      maps[2] = 2;
      maps[3] = 3;
      maps[4] = 4;
      maps[5] = 5;

      mapc[0] = 0;
      mapc[1] = 1;
      mapc[2] = 2;
      mapc[3] = 3;
      mapc[4] = 4;
      mapc[5] = 5;
    }

  /***************************************************************************/

  /***** Allocate memory for temporary arrays ********************************/

  /* Homogenized group constants */

  D = (double *)Mem(MEM_ALLOC, nG, sizeof(double));
  tot = (double *)Mem(MEM_ALLOC, nG, sizeof(double));
  nsf = (double *)Mem(MEM_ALLOC, nG, sizeof(double));
  chi = (double *)Mem(MEM_ALLOC, nG, sizeof(double));
  scatt = (double *)Mem(MEM_ALLOC, nG*nG, sizeof(double));

  /* Surface currents and fluxes */

  Js = (double *)Mem(MEM_ALLOC, ns*nG, sizeof(double));
  Jm = (double *)Mem(MEM_ALLOC, ns*nG, sizeof(double));
  vflx = (double *)Mem(MEM_ALLOC, nG, sizeof(double));
  sflx = (double *)Mem(MEM_ALLOC, ns*nG, sizeof(double));
  sflxsgn = (double *)Mem(MEM_ALLOC, ns*nG, sizeof(double));

  if (nc == 0)
    {
      Jc = NULL;
      cflx = NULL;
    }
  else
    {
      Jc = (double *)Mem(MEM_ALLOC, nc*nG, sizeof(double));
      cflx = (double *)Mem(MEM_ALLOC, nc*nG, sizeof(double));
    }

  /* Pin-wise fluxes */

  if (np == 0)
    pflx = NULL;
  else
    pflx = (double *)Mem(MEM_ALLOC, np*nG, sizeof(double));

  /* Leakage rate */

  L = (double *)Mem(MEM_ALLOC, nG, sizeof(double));

  /* JLe: Tässä voisi kuvata vähän tarkemmin että mitä nämä muuttujat on, */
  /* samaan tyyliin kuin yllä, niin selitys löytyisi yhdestä paikasta. */

  A = (complex **)Mem(MEM_ALLOC, nG, sizeof(complex *));
  T = (complex **)Mem(MEM_ALLOC, nG, sizeof(complex *));
  U = (complex **)Mem(MEM_ALLOC, nG, sizeof(complex *));
  F = (complex **)Mem(MEM_ALLOC, nG, sizeof(complex *));

  for (i = 0; i < nG; i++)
    {
      A[i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex));
      T[i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex));
      U[i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex));
      F[i] = (complex *)Mem(MEM_ALLOC, Nr, sizeof(complex));
    }

  b = (complex *)Mem(MEM_ALLOC, nG*(ns + nc), sizeof(complex));
  c = (complex *)Mem(MEM_ALLOC, nG*(ns + nc), sizeof(complex));

  /* Surface vectors */

  surfs_n = (double **)Mem(MEM_ALLOC, ns + nc, sizeof(double *));
  surfs_r = (double **)Mem(MEM_ALLOC, ns + nc, sizeof(double *));

  for (i = 0; i < ns + nc; i++)
    {
      surfs_n[i] = (double *)Mem(MEM_ALLOC, 3, sizeof(double));
      surfs_r[i] = (double *)Mem(MEM_ALLOC, 3, sizeof(double));
    }

  /***************************************************************************/

  /***** Read data ***********************************************************/

  /* Loop over energy groups */

  for (n = 0; n < nG; n++)
    {
      /* Flux used for group constant generation (NOTE: This is not */
      /* necessarily the same as the flux inside the ADF region)    */

      ptr = (long)RDB[gcu + GCU_INF_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr6)", DATA_ARRAY, ptr);
      if ((div = BufVal(ptr, n)) <= 0.0)
        Die(FUNCTION_NAME, "Zero or negative flux in group %ld", n + 1);

      /* Diffusion coefficient */

      if ((long)RDB[DATA_HOMOFLUX_DIFFCOEF] == HOMOFLUX_DIFF_INF)
        ptr = (long)RDB[gcu + GCU_INF_TRANSPXS];
      else if ((long)RDB[DATA_HOMOFLUX_DIFFCOEF] == HOMOFLUX_DIFF_TRC)
        ptr = (long)RDB[gcu + GCU_TRC_TRANSPXS];
      else
        Die(FUNCTION_NAME, "Invalid option");

      CheckPointer(FUNCTION_NAME, "(ptr7)", DATA_ARRAY, ptr);

      D[n] =  BufVal(ptr, n)/div/3.0;
      CheckValue(FUNCTION_NAME, "D", "", D[n], ZERO, INFTY);

      /* Total cross section (Tässä otetaan huomioon (n,xn)-monistus) */

      ptr = (long)RDB[gcu + GCU_INF_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr8)", DATA_ARRAY, ptr);
      tot[n] = BufVal(ptr, n)/div;

      ptr = (long)RDB[gcu + GCU_INF_RABSXS];
      CheckPointer(FUNCTION_NAME, "(ptr9)", DATA_ARRAY, ptr);
      tot[n] = tot[n] - BufVal(ptr, n)/div;

      ptr = (long)RDB[gcu + GCU_INF_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr10)", DATA_ARRAY, ptr);
      tot[n] = tot[n] + BufVal(ptr, n)/div;

      CheckValue(FUNCTION_NAME, "tot", "", tot[n], ZERO, INFTY);

      /* Fission neutron production */

      ptr = (long)RDB[gcu + GCU_INF_NSF];
      CheckPointer(FUNCTION_NAME, "(ptr11)", DATA_ARRAY, ptr);
      nsf[n] =  BufVal(ptr, n)/div;
      CheckValue(FUNCTION_NAME, "nsf", "", nsf[n], 0.0, INFTY);

      /* Fission spectrum */

      ptr = (long)RDB[gcu + GCU_INF_CHIT];
      CheckPointer(FUNCTION_NAME, "(ptr12)", DATA_ARRAY, ptr);
      chi[n] =  BufVal(ptr, n);
      CheckValue(FUNCTION_NAME, "chi", "", chi[n], 0.0, INFTY);

      /* Scattering matrix (tossa ei oteta (n,xn)-monistusta huomioon) */

      ptr = (long)RDB[gcu + GCU_INF_S0];
      CheckPointer(FUNCTION_NAME, "(ptr13)", DATA_ARRAY, ptr);

      for (m = 0; m < nG; m++)
        {
          scatt[n*nG + m] = BufVal(ptr, n, m)/div;
          CheckValue(FUNCTION_NAME, "scatt", "", scatt[n*nG + m], 0.0, INFTY);
        }

      /* Net current over surfaces */

      ptr = (long)RDB[gcu + GCU_RES_FG_DF_SURF_NET_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr14)", DATA_ARRAY, ptr);

      for (m = 0; m < ns; m++)
        {
          /* Remap currents */

          i = maps[m];
          Js[i*nG + n] = BufVal(ptr, m, n);
          CheckValue(FUNCTION_NAME, "Js", "", Js[i*nG + n], -INFTY, INFTY);
        }

      /* Net current over surfaces mid-points */

      ptr = (long)RDB[gcu + GCU_RES_FG_DF_MID_NET_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr15)", DATA_ARRAY, ptr);

      for (m = 0; m < ns; m++)
        {
          /* Remap currents */

          i = maps[m];
          Jm[i*nG + n] = BufVal(ptr, m, n);
          CheckValue(FUNCTION_NAME, "Jm", "", Jm[i*nG + n], -INFTY, INFTY);
        }

      /* Net current over corners */

      if (nc > 0)
        {
          ptr = (long)RDB[gcu + GCU_RES_FG_DF_CORN_NET_CURR];
          CheckPointer(FUNCTION_NAME, "(ptr16)", DATA_ARRAY, ptr);

          for (m = 0; m < nc; m++)
            {
              /* Remap current */

              i = mapc[m];
              Jc[i*nG + n] = BufVal(ptr, m, n);
              CheckValue(FUNCTION_NAME, "Jc", "", Jc[i*nG + n], -INFTY, INFTY);
            }
        }

      /* Volume-averaged heterogeneous flux */

      ptr = (long)RDB[gcu + GCU_RES_FG_DF_HET_VOL_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr17)", DATA_ARRAY, ptr);
      vflx[n] = BufVal(ptr, n);
      CheckValue(FUNCTION_NAME, "vflx", "", vflx[n], ZERO, INFTY);
    }

  /* Re-normalize fission spectrum */

  sum = 0.0;
  for (n = 0; n < nG; n++)
    sum = sum + chi[n];

  if (sum > 0.0)
    for (n = 0; n < nG; n++)
      chi[n] = chi[n]/sum;

  /* Calculate leakage rate (pintavirta ulottuu kulma-alueelle) */

  for (n = 0; n < nG; n++)
    L[n] = 0.0;

  for (n = 0; n < nG; n++)
    for (m = 0; m < ns; m++)
      L[n] = L[n] - Js[m*nG + n];

  /* Calculate source term */

  src = 0.0;
  for (n = 0; n < nG; n++)
    src = src + vol*vflx[n]*nsf[n];

  /* Calculate loss term */

  loss = 0.0;
  for (n = 0; n < nG; n++)
    {
      loss = loss + vol*vflx[n]*tot[n] + L[n];

      for (m = 0; m < nG; m++)
        loss = loss - vol*vflx[n]*scatt[n*nG + m];
    }

  /* Avoid compiler warning */

  keff = -1.0;

  /* Calculate k-eff */

  if (loss == 0.0)
    Die(FUNCTION_NAME, "Zero loss rate");
  else if (src > 0.0)
    keff = src/loss;
  else
    keff = 0.0;

  /***************************************************************************/

  /***** Get surface vectors and parameters **********************************/

  /* Check geometry type */

  switch (type)
    {
    case SURF_PX:
    case SURF_PY:
    case SURF_PZ:
    case SURF_PLANE:
      {
        /* Set surface normal vector */

        surfs_n[0][0] = 1.0;
        surfs_n[0][1] = 0.0;
        surfs_n[0][2] = 0.0;

        /* Set origin for surface integration */

        surfs_r[0][0] = 0.0;
        surfs_r[0][1] = 0.0;
        surfs_r[0][2] = 0.0;

        /* No corners */

        ifCorn = 0;

        /* Break case */

        break;
      }
    case SURF_SQC:
    case SURF_RECT:
      {
        /* Set surface normal vectors (E-N-W-S) */

        surfs_n[0][0] = 1.0;
        surfs_n[0][1] = 0.0;
        surfs_n[0][2] = 0.0;

        surfs_n[1][0] = 0.0;
        surfs_n[1][1] = 1.0;
        surfs_n[1][2] = 0.0;

        surfs_n[2][0] = -1.0;
        surfs_n[2][1] = 0.0;
        surfs_n[2][2] = 0.0;

        surfs_n[3][0] = 0.0;
        surfs_n[3][1] = -1.0;
        surfs_n[3][2] = 0.0;

        /* Set corner vectors (NE-NW-SW-SE) */;

        if (nc > 0)
          {
            /* Calculate factors */

            nx = ds[1]/sqrt(ds[0]*ds[0] + ds[1]*ds[1]);
            CheckValue(FUNCTION_NAME, "nx", "", nx, 0.0, 1.0);

            ny = sqrt(1.0 - nx*nx);
            CheckValue(FUNCTION_NAME, "ny", "", ny, 0.0, 1.0);

            /* Set vectors */

            surfs_n[4][0] = nx;
            surfs_n[4][1] = ny;
            surfs_n[4][2] = 0.0;

            surfs_n[5][0] = -nx;
            surfs_n[5][1] = ny;
            surfs_n[5][2] = 0.0;

            surfs_n[6][0] = -nx;
            surfs_n[6][1] = -ny;
            surfs_n[6][2] = 0.0;

            surfs_n[7][0] = nx;
            surfs_n[7][1] = -ny;
            surfs_n[7][2] = 0.0;
          }

        /* Check that all values match for square prism */

        if (type == SURF_SQC)
          for (n = 1; n < 4; n++)
            if (ds[n] != ds[n - 1])
              Die(FUNCTION_NAME, "Mismatch in ds");

        /* Check that x- and y-values match for rectangular prism */

        if ((ds[0] != ds[2]) || (ds[1] != ds[3]))
          Die(FUNCTION_NAME, "Mismatch in ds");

        /* Set origin for surface integration (E-N-W-S) */

        surfs_r[0][0] = 0.5*ds[1];
        surfs_r[0][1] = -0.5*ds[0];
        surfs_r[0][2] = 0.0;

        surfs_r[1][0] = 0.5*ds[1];
        surfs_r[1][1] = 0.5*ds[0];
        surfs_r[1][2] = 0.0;

        surfs_r[2][0] = -0.5*ds[1];
        surfs_r[2][1] = 0.5*ds[0];
        surfs_r[2][2] = 0.0;

        surfs_r[3][0] = -0.5*ds[1];
        surfs_r[3][1] = -0.5*ds[0];
        surfs_r[3][2] = 0.0;

        /* Break case */

        break;
      }
    case SURF_HEXYC:
    case SURF_HEXXC:
      {
        /* Set surface normal vectors */

        surfs_n[0][0] = 0.5*SQRT3;
        surfs_n[0][1] = 0.5;
        surfs_n[0][2] = 0.0;

        surfs_n[1][0] = 0.0;
        surfs_n[1][1] = 1.0;
        surfs_n[1][2] = 0.0;

        surfs_n[2][0] = -0.5*SQRT3;
        surfs_n[2][1] = 0.5;
        surfs_n[2][2] = 0.0;

        surfs_n[3][0] = -0.5*SQRT3;
        surfs_n[3][1] = -0.5;
        surfs_n[3][2] = 0.0;

        surfs_n[4][0] = 0.0;
        surfs_n[4][1] = -1.0;
        surfs_n[4][2] = 0.0;

        surfs_n[5][0] = 0.5*SQRT3;
        surfs_n[5][1] = -0.5;
        surfs_n[5][2] = 0.0;

        /* Set corner vectors */

        if (nc > 0)
          {
            surfs_n[6][0] = 0.5;
            surfs_n[6][1] = 0.5*SQRT3;
            surfs_n[6][2] = 0.0;

            surfs_n[7][0] = -0.5;
            surfs_n[7][1] = 0.5*SQRT3;
            surfs_n[7][2] = 0.0;

            surfs_n[8][0] = -1.0;
            surfs_n[8][1] = 0.0;
            surfs_n[8][2] = 0.0;

            surfs_n[9][0] = -0.5;
            surfs_n[9][1] = -0.5*SQRT3;
            surfs_n[9][2] = 0.0;

            surfs_n[10][0] = 0.5;
            surfs_n[10][1] = -0.5*SQRT3;
            surfs_n[10][2] = 0.0;

            surfs_n[11][0] = 1.0;
            surfs_n[11][1] = 0.0;
            surfs_n[11][2] = 0.0;
          }

        /* Check that all values match */

        for (n = 1; n < 6; n++)
          if (ds[n] != ds[n - 1])
            Die(FUNCTION_NAME, "Mismatch in ds");

        /* Set origin for surface integration */

        surfs_r[0][0] = ds[0];
        surfs_r[0][1] = 0.0;
        surfs_r[1][0] = 0.5*ds[0];
        surfs_r[1][1] = 0.5*SQRT3*ds[0];
        surfs_r[2][0] = -0.5*ds[0];
        surfs_r[2][1] = 0.5*SQRT3*ds[0];
        surfs_r[3][0] = -ds[0];
        surfs_r[3][1] = 0.0;
        surfs_r[4][0] = -0.5*ds[0];
        surfs_r[4][1] = -0.5*SQRT3*ds[0];
        surfs_r[5][0] = 0.5*ds[0];
        surfs_r[5][1] = -0.5*SQRT3*ds[0];

        /* Break case */

        break;
      }
    case SURF_TRIAG:
      {
        Warn(FUNCTION_NAME,
             "Diffusion flux solver not supported for boundary type %ld",
             type);

        return;
      }
    default:
      {
        /* Not supported */

        Die(FUNCTION_NAME,
            "Diffusion flux solver not supported for boundary type %ld",
            type);
      }
    }

  /***************************************************************************/

  /***** Calculate coefficients for basis functions **************************/

  /* Check k-eff */

  CheckValue(FUNCTION_NAME, "keff", "", keff, 0.0, 5.0);

  /* Form matrix A */

  for (j = 0; j < nG; j++)
    for (i = 0; i < nG; i++)
      {
        /* Check diffusion coefficient and cross sections */

        CheckValue(FUNCTION_NAME, "D", "", D[i], ZERO, INFTY);
        CheckValue(FUNCTION_NAME, "scatt", "", scatt[i], 0.0, INFTY);
        CheckValue(FUNCTION_NAME, "tot", "", tot[i], ZERO, INFTY);
        CheckValue(FUNCTION_NAME, "nsf", "", D[i], 0.0, INFTY);
        CheckValue(FUNCTION_NAME, "chi", "", chi[i], 0.0, 1.0);


        /* Calculate A(i,j) */

        A[i][j].re = -1.0/D[i]*scatt[j*nG + i];

        if (keff > 0.0)
          A[i][j].re = A[i][j].re - 1.0/D[i]*1.0/keff*chi[i]*nsf[j];

        if (i == j)
          A[i][i].re = A[i][i].re + 1.0/D[i]*tot[i];
      }

  /* Compute Schur Factorization of A (only matrices T and U */
  /* needed from now on) */

  if (nG == 1)
    {
      T[0][0].re = 1.0;
      U[0][0].re = A[0][0].re;
    }
  else
    schurFactorization(nG, A, T, U);

  /* OBS! Testatessa oli eroa etumerkeissä, muista tämä jos tulee */
  /* outoja tuloksia? */

  /*
  for (i=0; i<nG; i++)
    for (j=0; j<nG; j++)
      printf("T(%ld,%ld) = (%f,%f)\n", i+1,j+1,T[i][j].re,T[i][j].im);

  for (i=0; i<nG; i++)
    for (j=0; j<nG; j++)
      printf("U(%ld,%ld) = (%f,%f)\n", i+1,j+1,U[i][j].re,U[i][j].im);
  */

  /* Solve coefficients c from boundary conditions. Compute RHS of */
  /* boundary condition system: b = -inv(D)*J (remapping doesn't   */
  /* change the order of segment lengths */

  k = 0;
  n = 0;
  for (j = 0; j < ns; j++)
    for (i = 0; i < nG; i++)
      b[k++].re = -1.0/D[i]*Js[n++]/ds[j];

  n = 0;
  for (j = 0; j < nc; j++)
    for (i = 0; i < nG; i++)
      b[k++].re = -1.0/D[i]*Jc[n++]/dc;

  /* Solve coefficients c from boundary conditions */

  DFSolver(ifCorn, nG, ns + nc, ds, dc, surfs_n, surfs_r, T, U, b, c);

  /***************************************************************************/

  /***** Compute homogeneous surface flux ************************************/

  /* Loop over surfaces */

  for (i = 0; i < ns; i++)
    {
      /* Create array of points for integration */

      t = MakeArray(0.0, ds[maps[i]], Nr, 1);

      /* Tangent vector (u,v,w) */

      tn[0] = -surfs_n[i][1];
      tn[1] = surfs_n[i][0];
      tn[2] = 0.0;

      /* Compute solution at Nr points along boundary surface */

      for (k = 0; k < Nr; k++)
        {
          /* Set coordinates (x,y,z) */

          r[0] = surfs_r[i][0] + t[k]*tn[0];
          r[1] = surfs_r[i][1] + t[k]*tn[1];
          r[2] = 0.0;
          /*
          printf("s%ld(%ld,:) = [%E %E];\n", i + 1, k + 1, r[0], r[1]);
          */
          /* Solution at point (returns flux for all energy groups) */

          phi = dfSol(nG, ns + nc, surfs_n, r, T, U, c);

          /* Store row-wise */

          for (n = 0; n < nG; n++)
            {
              /* Check value (small negative value possible) */
              /* JLe: Onko -0.1 pieni? Huomaa että toi on mielivaltaisesti */
              /* normeerattu */
              /*
              CheckValue(FUNCTION_NAME, "phi","", phi[n], -1E-1, INFTY);
              */

              /* Check negative */

              if (phi[n] < 0.0)
                nneg++;

              /* Store value */

              F[n][k].re = phi[n];
            }

          /* Free flux array */

          Mem(MEM_FREE, phi);
        }

      /* Loop over energy groups */

      for (n = 0; n < nG; n++)
        {
          /* Trapezoidal integration */

          z = trapz(Nr, t, F[n]);

          /* Store flux */

          sflx[i*nG + n] = z.re/ds[maps[i]];
        }

      /* Free points array */

      Mem(MEM_FREE, t);
    }

  /***************************************************************************/

  /***** Compute homogeneous corner flux *************************************/

  /* Create array of points for integration */

  t = MakeArray(0.0, 0.5*dc, Nr, 1);

  /* Loop over corners */

  for (i = 0; i < nc; i++)
    {
      /* Index to second tangent vector */

      if (i == nc - 1)
        j = 0;
      else
        j = i + 1;

      /* First tangent vector */

      tn[0] = -surfs_n[i][1];
      tn[1] = surfs_n[i][0];
      tn[2] = 0.0;

      /* Loop over points */

      for (k = 0; k < Nr; k++)
        {
          /* Set coordinates (x,y,z) */

          r[0] = surfs_r[i][0] + t[k]*tn[0];
          r[1] = surfs_r[i][1] + t[k]*tn[1];
          r[2] = 0.0;

          /* Solution at point (returns flux for all energy groups) */

          phi = dfSol(nG, ns + nc, surfs_n, r, T, U, c);

          /* Store row-wise */

          for (n = 0; n < nG; n++)
            {
              /* Check negative */

              if (phi[n] < 0.0)
                nneg++;

              /* Store value */

              F[n][k].re = phi[n];
            }

          /* Free flux array */

          Mem(MEM_FREE, phi);
        }

      /* Loop over energy groups */

      for (n = 0; n < nG; n++)
        {
          /* Trapezoidal integration */

          z = trapz(Nr, t, F[n]);

          /* Store flux */

          if (i > 0)
            cflx[(i - 1)*nG + n] = z.re/dc;
          else
            cflx[(nc - 1)*nG + n] = z.re/dc;
        }

      /* Second tangent vector */

      tn[0] = -surfs_n[j][1];
      tn[1] = surfs_n[j][0];
      tn[2] = 0.0;

      /* Loop over points */

      for (k = 0; k < Nr; k++)
        {
          /* Set coordinates (x,y,z) */

          r[0] = surfs_r[i][0] + t[k]*tn[0];
          r[1] = surfs_r[i][1] + t[k]*tn[1];
          r[2] = 0.0;

          /* Solution at point (returns flux for all energy groups) */

          phi = dfSol(nG, ns + nc, surfs_n, r, T, U, c);

          /* Store row-wise */

          for (n = 0; n < nG; n++)
            {
              /* Check negative */

              if (phi[n] < 0.0)
                nneg++;

              /* Store value */

              F[n][k].re = phi[n];
            }

          /* Free flux array */

          Mem(MEM_FREE, phi);
        }

      /* Loop over energy groups */

      for (n = 0; n < nG; n++)
        {
          /* Trapezoidal integration */

          z = trapz(Nr, t, F[n]);

          /* Store flux */

          if (i > 0)
            cflx[(i - 1)*nG + n] = cflx[(i - 1)*nG + n] + z.re/dc;
          else
            cflx[(nc - 1)*nG + n] = cflx[(nc - 1)*nG + n] + z.re/dc;
        }
    }

  /* Free points array */

  Mem(MEM_FREE, t);

  /***************************************************************************/

  /***** Compute homogeneous sign weighted surface flux **********************/

  /* Loop over surfaces */

  for (i = 0; i < ns; i++)
    {

      /* j == 0 is the first half and j == 1 is the second half, and sgns gives
         the multiplier such that if it is 1.0, the first half is plus and the
         second half is minus */

      for (j = 0; j < 2; j++)
        {
          /* Create array of points for integration */

          if (j == 0)
            {
              t = MakeArray(0.0, 0.5*ds[maps[i]], Nr2, 1);
            }
          else
            {
              t = MakeArray(0.5*ds[maps[i]], ds[maps[i]], Nr2, 1);
            }

          /* Tangent vector (u,v,w) */

          tn[0] = -surfs_n[i][1];
          tn[1] = surfs_n[i][0];
          tn[2] = 0.0;

          /* Compute solution at Nr2 points along boundary surface */

          for (k = 0; k < Nr2; k++)
          {
              /* Set coordinates (x,y,z) */

              r[0] = surfs_r[i][0] + t[k]*tn[0];
              r[1] = surfs_r[i][1] + t[k]*tn[1];
              r[2] = 0.0;

              /* Solution at point (returns flux for all energy groups) */

              phi = dfSol(nG, ns + nc, surfs_n, r, T, U, c);

              /* Store row-wise */

              for (n = 0; n < nG; n++)
              {
                  /* Store value */

                  F[n][k].re = phi[n];
              }

              /* Free flux array */

              Mem(MEM_FREE, phi);
          }

          /* Loop over energy groups */

          for (n = 0; n < nG; n++)
          {
              /* Trapezoidal integration */

              z = trapz(Nr2, t, F[n]);

              /* Store flux */

              if (j == 0)
              {
                  sflxsgn[i*nG + n] = sgns[i]*z.re/ds[maps[i]];
              }
              else
              {
                  sflxsgn[i*nG + n] -= sgns[i]*z.re/ds[maps[i]];
              }
          }

          /* Free points array */

          Mem(MEM_FREE, t);
        }
    }

  /***************************************************************************/

  /***** Compute homogeneous flux for pin-power form functions ***************/

  /* Loop over pins and energy groups */

  for (m = 0; m < np; m++)
    for (n = 0; n < nG; n++)
      {
        /* Get pointer to power distribution (used for weighting */
        /* factor for coordinates) */

        ptr = (long)RDB[gcu + GCU_RES_FG_PPW_POW];
        CheckPointer(FUNCTION_NAME, "(ptr18)", DATA_ARRAY, ptr);

        /* Get coordinates */

        ptr = (long)RDB[gcu + GCU_RES_FG_PPW_XYZ];
        CheckPointer(FUNCTION_NAME, "(ptr19)", DATA_ARRAY, ptr);
        r[0] = BufVal(ptr, 0, m, n);
        r[1] = BufVal(ptr, 1, m, n);
        r[2] = BufVal(ptr, 2, m, n);

        /* Value -1E+18 is used to indicate zero scores */

        if (r[0] > -1E+18)
          {
            /* Solution at point (returns flux for all energy groups) */

            phi = dfSol(nG, ns + nc, surfs_n, r, T, U, c);

            /* Check negative */

            if (phi[n] < 0.0)
              nneg++;

            /* Store homogeneous flux */

            pflx[m*nG + n] = phi[n];

            /* Free flux array */

            Mem(MEM_FREE, phi);
          }
      }

  /***************************************************************************/

  /***************************************************************************/

  /***** Print final values for debugging ************************************/

  /* Check total flux */
  /*
  CalcTotFlx(surf, nG, ns + nc, surfs_n, T, U, c);
  */
  /* Print flux */

#ifdef DEBUG

  PrintFlx(gcu, surf, nG, ns + nc, surfs_n, T, U, c);

#endif

  /***************************************************************************/

  /***** Store homogeneous surface, corner and pin flux **********************/

  /* Mark scoring buffer unreduced */

  WDB[DATA_BUF_REDUCED] = (double)NO;

  /* Surface flux */

  ptr = (long)RDB[gcu + GCU_RES_FG_DF_HOM_SURF_FLUX];
  CheckPointer(FUNCTION_NAME, "(ptr20)", DATA_ARRAY, ptr);

  /* Sign moments of discontinuity factors */

  ptr2 = (long)RDB[gcu + GCU_RES_FG_DF_SGN_HOM_SURF_FLUX];
  CheckPointer(FUNCTION_NAME, "(ptr2X)", DATA_ARRAY, ptr2);

  for (n = 0; n < nG; n++)
    for (m = 0; m < ns; m++)
      {
        /* NOTE: Tässä oli jotain häikkää, mutta ilmeisesti oikein */

        if (1 != 2)
          i = maps[m];
        else
          i = m;

        /* Store data */

        AddBuf(sflx[i*nG + n], 1.0, ptr, 0, -1, m, n);

        /* Sign moments of discontinuity factors */

        AddBuf(sflxsgn[i*nG + n], 1.0, ptr2, 0, -1, m, n);
      }

  /* Corner flux */

  if (nc > 0)
    {
      ptr = (long)RDB[gcu + GCU_RES_FG_DF_HOM_CORN_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr21)", DATA_ARRAY, ptr);

      for (n = 0; n < nG; n++)
        for (m = 0; m < nc; m++)
          {
            /* NOTE: Tässä oli jotain häikkää, mutta ilmeisesti oikein */

            if (1 != 2)
              i = mapc[m];
            else
              i = nc - m - 1;

            /* Store data */

            AddBuf(cflx[i*nG + n], 1.0, ptr, 0, -1, m, n);
          }
    }

  /* Pin-flux */

  if (np > 0)
    {
      ptr = (long)RDB[gcu + GCU_RES_FG_PPW_HOM_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr22)", DATA_ARRAY, ptr);

      for (n = 0; n < nG; n++)
        for (m = 0; m < np; m++)
          AddBuf(pflx[m*nG + n], 1.0, ptr, 0, -1, m, n);
    }

  /* Reduce buffer */

  ReduceBuffer();

  /***************************************************************************/

  /***** Free allocated memory ***********************************************/

  Mem(MEM_FREE, D);
  Mem(MEM_FREE, tot);
  Mem(MEM_FREE, nsf);
  Mem(MEM_FREE, chi);
  Mem(MEM_FREE, scatt);
  Mem(MEM_FREE, Js);
  Mem(MEM_FREE, Jm);
  Mem(MEM_FREE, vflx);
  Mem(MEM_FREE, sflx);
  Mem(MEM_FREE, sflxsgn);

  if (Jc != NULL)
    Mem(MEM_FREE, Jc);

  if (cflx != NULL)
    Mem(MEM_FREE, cflx);

  if (pflx != NULL)
    Mem(MEM_FREE, pflx);

  Mem(MEM_FREE, L);

  for (i = 0; i < nG; i++)
    {
      Mem(MEM_FREE, A[i]);
      Mem(MEM_FREE, T[i]);
      Mem(MEM_FREE, U[i]);
      Mem(MEM_FREE, F[i]);
    }

  Mem(MEM_FREE, A);
  Mem(MEM_FREE, T);
  Mem(MEM_FREE, U);
  Mem(MEM_FREE, F);

  for (i = 0; i < ns + nc; i++)
    {
      Mem(MEM_FREE, surfs_n[i]);
      Mem(MEM_FREE, surfs_r[i]);
    }

  Mem(MEM_FREE, surfs_n);
  Mem(MEM_FREE, surfs_r);

  Mem(MEM_FREE, b);
  Mem(MEM_FREE, c);

  /* Check number of negative values */

  if (nneg > 0)
    {
      /* Pointer to universe */

      ptr = (long)RDB[gcu + GCU_PTR_UNIV];
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Print warning */

      Note(0, "Negative values in homogeneous flux solution in universe %s",
           GetText(ptr + UNIVERSE_PTR_NAME));
    }

  /***************************************************************************/
}

/*****************************************************************************/

/***** Calculate mean integral flux for comparison ***************************/

void CalcTotFlx(long surf, long nG, long nj, double **surfs_n,
                complex **T, complex **U, complex *c)
{
  long type, ptr, np, n, i;
  double xmin, xmax, ymin, ymax, r[3], *f, *flx;

  /* Check surface pointer */

  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Get surface type */

  type = (long)RDB[surf + SURFACE_TYPE];

  /* Get pointer to parameters */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
  CheckPointer(FUNCTION_NAME, "(ptr23)", DATA_ARRAY, ptr);

  /* Avoid compiler warning */

  xmin = 0.0;
  xmax = 0.0;
  ymin = 0.0;
  ymax = 0.0;

  /* Check type */

  /* Check type */

  if (type == SURF_SQC)
    {
      /* Get bounding box */

      xmin = -RDB[ptr + 2] + RDB[ptr];
      xmax = RDB[ptr + 2] + RDB[ptr];
      ymin = -RDB[ptr + 2] + RDB[ptr + 1];
      ymax = RDB[ptr + 2] + RDB[ptr + 1];
    }
  else if (type == SURF_RECT)
    {
      /* Get bounding box */

      xmin = RDB[ptr];
      xmax = RDB[ptr + 1];
      ymin = RDB[ptr + 2];
      ymax = RDB[ptr + 3];
    }
  else if (type == SURF_HEXXC)
    {
      /* Get bounding box */

      xmin = -RDB[ptr + 2] + RDB[ptr];
      xmax = RDB[ptr + 2] + RDB[ptr];
      ymin = -2.0/SQRT3*RDB[ptr + 2] + RDB[ptr + 1];
      ymax = 2.0/SQRT3*RDB[ptr + 2] + RDB[ptr + 1];
    }
  else if (type == SURF_HEXYC)
    {
      /* Get bounding box */

      xmin = -2.0/SQRT3*RDB[ptr + 2] + RDB[ptr];
      xmax = 2.0/SQRT3*RDB[ptr + 2] + RDB[ptr];
      ymin = -RDB[ptr + 2] + RDB[ptr + 1];
      ymax = RDB[ptr + 2] + RDB[ptr + 1];
    }
  else
    Die(FUNCTION_NAME, "Invalid surface type");

  /* Set number of random points */

  np = 100000;

  /* Allocate memory for total flux */

  flx = (double *)Mem(MEM_ALLOC, nG, sizeof(double));

  /* Loop over points */

  for (n = 0; n < np; n++)
    {
      /* Sample coordinates inside surface */

      do
        {
          r[0] = RandF(0)*(xmax - xmin) + xmin;
          r[1] = RandF(0)*(ymax - ymin) + ymin;
          r[0] = 0.0;
        }
      while (TestSurface(surf, r[0], r[1], r[2], NO, 0) == NO);

      /* Calculate flux in point */

      f = dfSol(nG, nj, surfs_n, r, T, U, c);

      /* Add to sum */

      for (i = 0; i < nG; i++)
        flx[i] = flx[i] + f[i]/((double)np);

      /* Free array */

      Mem(MEM_FREE, f);
    }

  /* Print */

  printf("vflx (MC) = ");
  for (n = 0; n < nG; n++)
    printf("%1.5E ", flx[n]);
  printf("\n");

  /* Free array */

  Mem(MEM_FREE, flx);
}

/*****************************************************************************/

/***** Print flux for plotting ***********************************************/

void PrintFlx(long gcu, long surf, long nG, long nj, double **surfs_n,
              complex **T, complex **U, complex *c)
{
  long uni, type, ptr, nx, ny, n, m, i;
  double xmin, xmax, ymin, ymax, r[3], *f, x0, y0, z0;
  FILE *fp;
  char tmpstr[MAX_STR];

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);
  CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);

  /* Get pointer to universe */

  uni = (long)RDB[gcu + GCU_PTR_UNIV];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Get surface type */

  type = (long)RDB[surf + SURFACE_TYPE];

  /* Get pointer to parameters */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
  CheckPointer(FUNCTION_NAME, "(ptr24)", DATA_ARRAY, ptr);

  /* Avoid compiler warning */

  xmin = 0.0;
  xmax = 0.0;
  ymin = 0.0;
  ymax = 0.0;

  x0 = 0.0;
  y0 = 0.0;
  z0 = 0.0;

  /* Check type */

  if (type == SURF_SQC)
    {
      /* Get origin */

      x0 = RDB[ptr];
      y0 = RDB[ptr + 1];
      z0 = 0.0;

      /* Calculate bounding box */

      xmin = -RDB[ptr + 2] + 1E-6;
      xmax = RDB[ptr + 2] - 1E-6;
      ymin = -RDB[ptr + 2] + 1E-6;
      ymax = RDB[ptr + 2] - 1E-6;
    }
  else if (type == SURF_RECT)
    {
      /* Get origin */

      x0 = 0.5*(RDB[ptr] + RDB[ptr + 1]);
      y0 = 0.5*(RDB[ptr + 2] + RDB[ptr + 3]);
      z0 = 0.0;

      /* Calculate bounding box */

      xmin = -0.5*(RDB[ptr + 1] - RDB[ptr]) + 1E-6;
      xmax = 0.5*(RDB[ptr + 1] - RDB[ptr]) - 1E-6;
      ymin = -0.5*(RDB[ptr + 3] - RDB[ptr + 2]) + 1E-6;
      ymax = 0.5*(RDB[ptr + 3] - RDB[ptr + 2]) - 1E-6;
    }
  else if (type == SURF_HEXXC)
    {
      /* Get origin */

      x0 = RDB[ptr];
      y0 = RDB[ptr + 1];
      z0 = 0.0;

      /* Calculate bounding box */

      xmin = -RDB[ptr + 2];
      xmax = RDB[ptr + 2];
      ymin = -2.0/SQRT3*RDB[ptr + 2];
      ymax = 2.0/SQRT3*RDB[ptr + 2];
    }
  else if (type == SURF_HEXYC)
    {
      /* Get origin */

      x0 = RDB[ptr];
      y0 = RDB[ptr + 1];
      z0 = 0.0;

      /* Calculate bounding box */

      xmin = -2.0/SQRT3*RDB[ptr + 2];
      xmax = 2.0/SQRT3*RDB[ptr + 2];
      ymin = -RDB[ptr + 2];
      ymax = RDB[ptr + 2];
    }
  else
    Die(FUNCTION_NAME, "Invalid surface type");

  /* Set size */

  nx = 50;
  ny = 50;

  /* Open file for writing */

  sprintf(tmpstr,"U%s_hflx.m", GetText(uni + UNIVERSE_PTR_NAME));
  fp = fopen(tmpstr, "w");

  /* Print arrays */

  fprintf(fp, "x = [\n");
  for (n = 0; n < nx; n++)
    fprintf(fp, "%E\n", ((double)n)/((double)nx - 1.0)*(xmax - xmin)
            + xmin + x0);
  fprintf(fp, "];\n");

  fprintf(fp, "y = [\n");
  for (m = 0; m < ny; m++)
    fprintf(fp, "%E\n", ((double)m)/((double)ny - 1.0)*(ymax - ymin)
            + ymin + y0);
  fprintf(fp, "];\n");

  /* Init matrices */

  for (i = 0; i < nG; i++)
    fprintf(fp, "f%ld = zeros(%ld,%ld);\n", i + 1, nx, ny);

  /* Loop over points */

  for (n = 0; n < nx; n++)
    for (m = 0; m < ny; m++)
      {
        /* Set coordinates */

        r[0] = ((double)n)/((double)nx - 1.0)*(xmax - xmin) + xmin;
        r[1] = ((double)m)/((double)ny - 1.0)*(ymax - ymin) + ymin;
        r[2] = 0.0;

        /* Check if inside */

        if (TestSurface(surf, x0 + r[0], y0 + r[1], z0 + r[2], YES, 0) == NO)
          continue;

        /* Calculate flux in point */

        f = dfSol(nG, nj, surfs_n, r, T, U, c);

        /* Print values */

        for (i = 0; i < nG; i++)
          fprintf(fp, "f%ld(%ld, %ld) = %1.5E;\n", i + 1, n + 1, m + 1, f[i]);

        /* Free array */

        Mem(MEM_FREE, f);
      }

  /* Close file */

  fclose(fp);
}

/*****************************************************************************/
