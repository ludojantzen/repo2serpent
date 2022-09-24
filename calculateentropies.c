/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculateentropies.c                           */
/*                                                                           */
/* Created:       2011/11/19 (JLe)                                           */
/* Last modified: 2012/10/12 (JLe)                                           */
/* Version:       2.1.9                                                      */
/*                                                                           */
/* Description: Calculates fission source entropies                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateEntropies:"

/*****************************************************************************/

void CalculateEntropies()
{
  long nx, ny, nz, part, i, j, k, idx, ptr;
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, wgt, totw, totp;
  double *spt, *swg, entrp, entrw, sump, sumw;

  /* Check if entropies are calculated */

  if ((long)RDB[DATA_OPTI_ENTROPY_CALC] == NO)
    return;
  
  /* Get mesh size */

  nx = (long)RDB[DATA_ENTROPY_NX];
  ny = (long)RDB[DATA_ENTROPY_NY];
  nz = (long)RDB[DATA_ENTROPY_NZ];

  /* Check values */

  CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 500);
  CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 500);
  CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 500);

  /* Get boundaries */

  xmin = RDB[DATA_ENTROPY_XMIN];
  xmax = RDB[DATA_ENTROPY_XMAX];
  ymin = RDB[DATA_ENTROPY_YMIN];
  ymax = RDB[DATA_ENTROPY_YMAX];
  zmin = RDB[DATA_ENTROPY_ZMIN];
  zmax = RDB[DATA_ENTROPY_ZMAX];

  /* Check values */

  CheckValue(FUNCTION_NAME, "xmin", "", xmin, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "xmax", "", xmax, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "ymin", "", ymin, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "ymax", "", ymax, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "zmax", "", zmax, -INFTY, INFTY);

  /* Check order */

  if (xmin > xmax)
    Die(FUNCTION_NAME, "Error in boundaries");

  if (ymin > ymax)
    Die(FUNCTION_NAME, "Error in boundaries");

  if (zmin > zmax)
    Die(FUNCTION_NAME, "Error in boundaries");

  /* Allocate memory for temporary arrays  */

  spt = (double *)Mem(MEM_ALLOC, nx*ny*nz, sizeof(double));
  swg = (double *)Mem(MEM_ALLOC, nx*ny*nz, sizeof(double));

  /* Reset data */

  memset(spt, 0.0, nx*ny*nz*sizeof(double));
  memset(swg, 0.0, nx*ny*nz*sizeof(double));
  
  /* Reset totals */

  totw = 0.0;
  totp = 0.0;

  /* Pointer to first item after dummy */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  part = NextItem(part);

  /* Loop over source and collect points */

  while(part > VALID_PTR)
    {
      /* Get coordinates and weight */

      x = RDB[part + PARTICLE_X];
      y = RDB[part + PARTICLE_Y];
      z = RDB[part + PARTICLE_Z];

      wgt = RDB[part + PARTICLE_WGT];

      /* Get local normalised co-ordinates */
      
      if (xmin != xmax)
        x = (x - xmin)/(xmax - xmin);
      else
        x = 0.0;
      
      if (ymin != ymax)
        y = (y - ymin)/(ymax - ymin);
      else
        y = 0.0;

      if (zmin != zmax)
        z = (z - zmin)/(zmax - zmin);
      else
        z = 0.0;

      /* Check */

      if ((x >= 0.0) && (x < 1.0) && (y >= 0.0) && (y < 1.0) &&
          (z >= 0.0) && (z < 1.0))
        {
          /* Calculate indexes */
      
          i = (long)(x*nx);
          j = (long)(y*ny);
          k = (long)(z*nz);

          /* Calculate index */

          idx = i + nx*j + nx*ny*k;

          /* Add points and weight */

          spt[idx] = spt[idx] + 1.0;
          swg[idx] = swg[idx] + wgt;
          
          /* Add to totals */

          totp = totp + 1.0;
          totw = totw + wgt;
        }

      /* Next particle */

      part = NextItem(part);
    }

  /* X-component of entropy */
  
  entrp = 0.0;
  entrw = 0.0;

  for (i = 0; i < nx; i++)
    {
      /* Sum over values */
      
      sump = 0.0;
      sumw = 0.0;
      
      for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
          {
            /* Calculate index */
            
            idx = i + nx*j + nx*ny*k;
          
            /* Add to sums */

            sump = sump + spt[idx]/totp;
            sumw = sumw + swg[idx]/totw;
          }

      /* Add to entropies */

      if ((sump > 0.0) && (sump < 1.0))
        entrp = entrp - sump*log2(sump);

      if ((sumw > 0.0) && (sumw < 1.0))
        entrw = entrw - sumw*log2(sumw);
    }
  
  /* Divide by uniform source entropy */
  
  if (nx > 1)
    {
      entrp = entrp/log2(nx);
      entrw = entrw/log2(nx);
    }

  /* Add to statistics */

  ptr = (long)RDB[DATA_ENTROPY_PTR_SPT_STAT];  
  AddStat(entrp, ptr, 1);

  ptr = (long)RDB[DATA_ENTROPY_PTR_SWG_STAT];  
  AddStat(entrw, ptr, 1);

  /* Y-component of entropy */
  
  entrp = 0.0;
  entrw = 0.0;

  for (j = 0; j < ny; j++)
    {
      /* Sum over values */
      
      sump = 0.0;
      sumw = 0.0;

      for (i = 0; i < nx; i++)
        for (k = 0; k < nz; k++)
          {
            /* Calculate index */
            
            idx = i + nx*j + nx*ny*k;
          
            /* Add to sums */

            sump = sump + spt[idx]/totp;
            sumw = sumw + swg[idx]/totw;
          }

      /* Add to entropies */

      if ((sump > 0.0) && (sump < 1.0))
        entrp = entrp - sump*log2(sump);

      if ((sumw > 0.0) && (sumw < 1.0))
        entrw = entrw - sumw*log2(sumw);
    }
  
  /* Divide by uniform source entropy */
  
  if (ny > 1)
    {
      entrp = entrp/log2(ny);
      entrw = entrw/log2(ny);
    }

  /* Add to statistics */

  ptr = (long)RDB[DATA_ENTROPY_PTR_SPT_STAT];  
  AddStat(entrp, ptr, 2);

  ptr = (long)RDB[DATA_ENTROPY_PTR_SWG_STAT];  
  AddStat(entrw, ptr, 2);

  /* Z-component of entropy */
  
  entrp = 0.0;
  entrw = 0.0;

  for (k = 0; k < nz; k++)
    {
      /* Sum over values */
      
      sump = 0.0;
      sumw = 0.0;

      for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
          {
            /* Calculate index */
            
            idx = i + nx*j + nx*ny*k;
          
            /* Add to sums */

            sump = sump + spt[idx]/totp;
            sumw = sumw + swg[idx]/totw;
          }

      /* Add to entropies */

      if ((sump > 0.0) && (sump < 1.0))
        entrp = entrp - sump*log2(sump);

      if ((sumw > 0.0) && (sumw < 1.0))
        entrw = entrw - sumw*log2(sumw);
    }
  
  /* Divide by uniform source entropy */
  
  if (nz > 1)
    {
      entrp = entrp/log2(nz);
      entrw = entrw/log2(nz);
    }

  /* Add to statistics */

  ptr = (long)RDB[DATA_ENTROPY_PTR_SPT_STAT];  
  AddStat(entrp, ptr, 3);

  ptr = (long)RDB[DATA_ENTROPY_PTR_SWG_STAT];  
  AddStat(entrw, ptr, 3);

  /* Total entropy */
  
  entrp = 0.0;
  entrw = 0.0;

  /* Loop over values */
      
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++)
        {
          /* Calculate index */
          
          idx = i + nx*j + nx*ny*k;
          
          /* Get value */

            sump = spt[idx]/totp;
            sumw = swg[idx]/totw;
            
            /* Add to entropies */

            if ((sump > 0.0) && (sump < 1.0))
              entrp = entrp - sump*log2(sump);
            
            if ((sumw > 0.0) && (sumw < 1.0))
              entrw = entrw - sumw*log2(sumw);
        }
  
  /* Divide by uniform source entropy */
  
  if (nx*ny*nz > 1)
    {
      entrp = entrp/log2(nx*ny*nz);
      entrw = entrw/log2(nx*ny*nz);
    }

  /* Add to statistics */

  ptr = (long)RDB[DATA_ENTROPY_PTR_SPT_STAT];  
  AddStat(entrp, ptr, 0);

  ptr = (long)RDB[DATA_ENTROPY_PTR_SWG_STAT];  
  AddStat(entrw, ptr, 0);

  /* Free temporary arrays */

  Mem(MEM_FREE, spt);
  Mem(MEM_FREE, swg);
}

/*****************************************************************************/
