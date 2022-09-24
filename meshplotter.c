/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshplotter.c                                  */
/*                                                                           */
/* Created:       2011/03/25 (JLe)                                           */
/* Last modified: 2018/08/30 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Creates png-format mesh plots                                */
/*                                                                           */
/* Comments: - Noi minimit ja maksimit pitäis ottaa detektoridatasta että    */
/*             time binien kanssa tehdyt plotit normeerautuu oikein          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifndef NO_GFX_MODE
#include <gd.h>
#endif

#define FUNCTION_NAME "MeshPlotter:"

#define COLOURS 127

/*****************************************************************************/

void MeshPlotter()
{

#ifndef NO_GFX_MODE

  gdImagePtr im;
  long n, m, i, ptv1, ptd1, ptv2, ptd2, sym, type;
  long palette[COLOURS*2 + 2], R[COLOURS*2 + 2], G[COLOURS*2 + 2], 
    B[COLOURS*2 + 2];
  long sx, sy, mpl, i1, i2, ndis;
  double max1, min1, ave1, max2, min2, ave2;
  double tmpv1, tmpd1, tmpv2, tmpd2;
  double **val1, **val2;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) 
    return;

  /* Loop over meshes */

  mpl = (long)RDB[DATA_PTR_MPL0];
  while (mpl > VALID_PTR)
    {
      /* Reduce private results */

      ReducePrivateRes();

      /* Put filename */
 
      if ((long)RDB[DATA_RUN_VR_ITER] == YES)
        sprintf(tmpstr, "%s_vr%ld.png", GetText(mpl + MPL_PTR_FNAME),
                (long)RDB[DATA_VR_ITER_IDX]);
      else if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
        sprintf(tmpstr, "%s_bstep%ld.png", GetText(mpl + MPL_PTR_FNAME),
                (long)RDB[DATA_BURN_STEP]);
      else
        sprintf(tmpstr, "%s.png", GetText(mpl + MPL_PTR_FNAME));

      /* Open file for writing */

      if ((fp = fopen(tmpstr, "w")) == NULL)
        Die(FUNCTION_NAME, "Unable to open file for writing");

      /* Get pointers */
      
      ptv1 = (long)RDB[mpl + MPL_PTR_VAL1];
      ptd1 = (long)RDB[mpl + MPL_PTR_DIV1];
      ptv2 = (long)RDB[mpl + MPL_PTR_VAL2];
      ptd2 = (long)RDB[mpl + MPL_PTR_DIV2];

      /* Mesh size */

      sx = (long)RDB[mpl + MPL_NX];
      sy = (long)RDB[mpl + MPL_NY];

      /* Get type and number of distributions */

      type = (long)RDB[mpl + MPL_TYPE];
      ndis = (long)RDB[mpl + MPL_NDIST];

      /***********************************************************************/
      
      /***** Allocate memory for distribution matrices ***********************/

      if (ptv1 > VALID_PTR)
        {
          val1 = (double **)Mem(MEM_ALLOC, sx, sizeof(double *));
          
          for(n = 0; n < sx; n++)
            val1[n] = (double *)Mem(MEM_ALLOC, sy, sizeof(double));
        }
      else
        val1 = NULL;

      if (ptv2 > VALID_PTR)
        {
          val2 = (double **)Mem(MEM_ALLOC, sx, sizeof(double *));
          
          for(n = 0; n < sx; n++)
            val2[n] = (double *)Mem(MEM_ALLOC, sy, sizeof(double));
        }
      else
        val2 = NULL;

      /***********************************************************************/

      /***** Create distributions ********************************************/

      /* Check symmetry */

      if ((sym = (long)RDB[mpl + MPL_SYM]) == 0)
        {
          /***** No symmetry *************************************************/

          /* Loop over values */

          for (n = 0; n < sx; n++)
            for(m = 0; m < sy; m++)
              {
                /* Check pointer to first distribution */

                if (ptv1 > VALID_PTR)
                  {
                    /* Get value */
                    
                    if ((long)RDB[mpl + MPL_AX] < 4)
                      tmpv1 = ReadMesh(ptv1, n, m, 0);
                    else
                      tmpv1 = ReadMesh(ptv1, n, 0, sy - m - 1);

                    /* Check divider */

                    if (ptd1 > VALID_PTR)
                      if ((tmpd1 = ReadMesh(ptd1, n, m, 0)) != 0.0)
                        tmpv1 = tmpv1/tmpd1;

                    /* Write value in matrix */

                    val1[n][m] = val1[n][m] + tmpv1;
                  }

                /* Check pointer to second distribution */

                if (ptv2 > VALID_PTR)
                  {
                    /* Get value */
                    
                    tmpv2 = ReadMesh(ptv2, n, m, 0);
                
                    /* Check divider */

                    if (ptd2 > VALID_PTR)
                      if ((tmpd2 = ReadMesh(ptd2, n, m, 0)) != 0.0)
                        tmpv2 = tmpv2/tmpd2;

                    /* Write value in matrix */

                    val2[n][m] = val2[n][m] + tmpv2;
                  }

              }
          
          /*******************************************************************/
        }
      else if (sym == 2)
        {
          /***** 1/2 symmetry ************************************************/

          /* Check matrix size */

          if (sx != sy)
            Die(FUNCTION_NAME, "Symmetry error");

          /* Loop over values */

          for (n = 0; n < sx; n++)
            for(m = 0; m < sy; m++)
              {
                /* Check pointer to first distribution */

                if (ptv1 > VALID_PTR)
                  {
                    /* Get value */
                    
                    tmpv1 = ReadMesh(ptv1, n, m, 0);
                
                    /* Check divider */

                    if (ptd1 > VALID_PTR)
                      if ((tmpd1 = ReadMesh(ptd1, n, m, 0)) != 0.0)
                        tmpv1 = tmpv1/tmpd1;

                    /* Write value in matrix */

                    val1[n][m] = val1[n][m] + tmpv1;
                    val1[m][n] = val1[n][m] + tmpv1;
                  }

                /* Check pointer to second distribution */

                if (ptv2 > VALID_PTR)
                  {
                    /* Get value */
                    
                    tmpv2 = ReadMesh(ptv2, n, m, 0);
                
                    /* Check divider */

                    if (ptd2 > VALID_PTR)
                      if ((tmpd2 = ReadMesh(ptd2, n, m, 0)) != 0.0)
                        tmpv2 = tmpv2/tmpd2;

                    /* Write value in matrix */

                    val2[n][m] = val2[n][m] + tmpv2;
                    val2[m][n] = val2[n][m] + tmpv2;
                  }
              }
          
          /*******************************************************************/
        }

      else if (sym == 4)
        {
          /***** 1/4 symmetry ************************************************/

          /* Check matrix size */

          if (sx != sy)
            Die(FUNCTION_NAME, "Symmetry error");

          /* Loop over values */

          for (n = 0; n < sx; n++)
            for(m = 0; m < sy; m++)
              {
                /* Check pointer to first distribution */

                if (ptv1 > VALID_PTR)
                  {
                    /* Get value */
                    
                    tmpv1 = ReadMesh(ptv1, n, m, 0);
                
                    /* Check divider */

                    if (ptd1 > VALID_PTR)
                      if ((tmpd1 = ReadMesh(ptd1, n, m, 0)) != 0.0)
                        tmpv1 = tmpv1/tmpd1;

                    /* Write value in matrix */

                    val1[n][m] = val1[n][m] + tmpv1; 
                    val1[sx - n - 1][m] = val1[sx - n - 1][m] + tmpv1; 
                    val1[sx - n - 1][sy - m - 1] = 
                      val1[sx - n - 1][sy - m - 1] + tmpv1; 
                    val1[n][sy - m - 1] = val1[n][sy - m - 1] + tmpv1; 
                  }

                /* Check pointer to second distribution */

                if (ptv2 > VALID_PTR)
                  {
                    /* Get value */
                    
                    tmpv2 = ReadMesh(ptv2, n, m, 0);
                
                    /* Check divider */

                    if (ptd2 > VALID_PTR)
                      if ((tmpd2 = ReadMesh(ptd2, n, m, 0)) != 0.0)
                        tmpv2 = tmpv2/tmpd2;

                    /* Write value in matrix */

                    val2[n][m] = val2[n][m] + tmpv2; 
                    val2[sx - n - 1][m] = val2[sx - n - 1][m] + tmpv2; 
                    val2[sx - n - 1][sy - m - 1] = 
                      val2[sx - n - 1][sy - m - 1] + tmpv2; 
                    val2[n][sy - m - 1] = val2[n][sy - m - 1] + tmpv2; 
                  }
              }
          
          /*******************************************************************/
        }

      else if (sym == 8)
        {
          /***** 1/8 symmetry ************************************************/

          /* Check matrix size */

          if (sx != sy)
            Die(FUNCTION_NAME, "Symmetry error");

          /* Loop over values */

          for (n = 0; n < sx; n++)
            for(m = 0; m < sy; m++)
              {
                /* Check pointer to first distribution */

                if (ptv1 > VALID_PTR)
                  {
                    /* Get value */
                    
                    tmpv1 = ReadMesh(ptv1, n, m, 0);
                
                    /* Check divider */

                    if (ptd1 > VALID_PTR)
                      if ((tmpd1 = ReadMesh(ptd1, n, m, 0)) != 0.0)
                        tmpv1 = tmpv1/tmpd1;

                    /* Write value in matrix */

                    val1[n][m] = val1[n][m] + tmpv1; 
                    val1[sx-n-1][m] = val1[sx-n-1][m] + tmpv1; 
                    val1[sx-n-1][sy-m-1] = val1[sx-n-1][sy-m-1] + tmpv1; 
                    val1[n][sy-m-1] = val1[n][sy-m-1] + tmpv1; 
                    val1[m][n] = val1[m][n] + tmpv1; 
                    val1[sy-m-1][n] = val1[sy-m-1][n] + tmpv1; 
                    val1[sy-m-1][sx-n-1] = val1[sy-m-1][sx-n-1] + tmpv1; 
                    val1[m][sx-n-1] = val1[m][sx-n-1] + tmpv1; 
                  }

                /* Check pointer to second distribution */

                if (ptv2 > VALID_PTR)
                  {
                    /* Get value */
                    
                    tmpv2 = ReadMesh(ptv2, n, m, 0);
                
                    /* Check divider */

                    if (ptd2 > VALID_PTR)
                      if ((tmpd2 = ReadMesh(ptd2, n, m, 0)) != 0.0)
                        tmpv2 = tmpv2/tmpd2;

                    /* Write value in matrix */

                    val2[n][m] = val2[n][m] + tmpv2; 
                    val2[sx-n-1][m] = val2[sx-n-1][m] + tmpv2; 
                    val2[sx-n-1][sy-m-1] = val2[sx-n-1][sy-m-1] + tmpv2; 
                    val2[n][sy-m-1] = val2[n][sy-m-1] + tmpv2; 
                    val2[m][n] = val2[m][n] + tmpv2; 
                    val2[sy-m-1][n] = val2[sy-m-1][n] + tmpv2; 
                    val2[sy-m-1][sx-n-1] = val2[sy-m-1][sx-n-1] + tmpv2; 
                    val2[m][sx-n-1] = val2[m][sx-n-1] + tmpv2; 
                  }
              }
          
          /*******************************************************************/
        }
      else
        Die(FUNCTION_NAME, "Invalid symmetry option");
      
      /***********************************************************************/

      /***** Adjust and normalize distributions ******************************/

      if ((long)RDB[mpl + MPL_COLOR_SCALE] == COLOR_SCALE_LIN)
        {
          /*******************************************************************/
          
          /***** Linear color scheme *****************************************/
        
          /* Calculate average and get minimum */

          ave1 = 0.0;
          min1 = 0.0;
          
          ave2 = 0.0;
          min2 = 0.0;
          
          i1 = 0;
          i2 = 0;
          
          /* Loop over distributions */

          for (n = 0; n < sx; n++)
            for(m = 0; m < sy; m++)
              {
                /* Distribution 1 */
                
                if (val1 != NULL)
                  if (val1[n][m] > 0.0)
                    {
                      /* Add to average */
                      
                      ave1 = ave1 + val1[n][m];
                      
                      /* Compare to minimum */
                      
                      if (val1[n][m] < min1)
                        min1 = val1[n][m];
                      
                      /* Increase counter */
                      
                      i1++;
                    }
                
                /* Distribution 2 */
                
                if (val2 != NULL)
                  if (val2[n][m] > 0.0)
                    {
                      /* Add to average */
                      
                      ave2 = ave2 + val2[n][m];
                      
                      /* Compare to minimum */
                      
                      if (val2[n][m] < min2)
                        min2 = val2[n][m];
                      
                      /* Increase counter */
                      
                      i2++;
                    }
              }
          
          /* Divide average by total */
          
          if (i1 > 0)
            ave1 = ave1/((double)i1);
          else
            ave1 = 0.0;
          
          if (i2 > 0)
            ave2 = ave2/((double)i2);
          else
            ave2 = 0.0;
          
          /* Set maximum values */
          
          max1 = 1.00*ave1*2.0;
          max2 = 1.00*ave2*2.0;
          
          /* Some plots are special */
          
          if ((type == MPL_TYPE_DT_NEFF) || (type == MPL_TYPE_DT_GEFF) ||
              (type == MPL_TYPE_DENSITY))
            {
              min1 = 0.0;
              max1 = 1.0;
              min2 = 0.0;
              max2 = 1.0;
            }
          
          /* Check values */
          
          CheckValue(FUNCTION_NAME, "min1", "", min1, 0.0, INFTY);
          CheckValue(FUNCTION_NAME, "ave1", "", ave1, 0.0, INFTY);
          CheckValue(FUNCTION_NAME, "max1", "", max1, 0.0, INFTY);
          
          CheckValue(FUNCTION_NAME, "min2", "", min2, 0.0, INFTY);
          CheckValue(FUNCTION_NAME, "ave2", "", ave2, 0.0, INFTY);
          CheckValue(FUNCTION_NAME, "max2", "", max2, 0.0, INFTY);
          
          /* Check if last cycle */
          
          if (RDB[DATA_SIMULATION_COMPLETED] == 1.0)
            {
              /* Last cycle. Remember or store values */
              
              if (RDB[mpl + MPL_MIN1] < 0.0)
                WDB[mpl + MPL_MIN1] = min1;
              else
                min1 = RDB[mpl + MPL_MIN1];
              
              if (RDB[mpl + MPL_MAX1] < 0.0)
                WDB[mpl + MPL_MAX1] = max1;
              else
                max1 = RDB[mpl + MPL_MAX1];
              
              if (RDB[mpl + MPL_MIN2] < 0.0)
                WDB[mpl + MPL_MIN2] = min2;
              else
                min2 = RDB[mpl + MPL_MIN2];
              
              if (RDB[mpl + MPL_MAX2] < 0.0)
                WDB[mpl + MPL_MAX2] = max2;
              else
                max2 = RDB[mpl + MPL_MAX2];
            }
          
          /* Normalize */
          
          for (n = 0; n < sx; n++)
            for(m = 0; m < sy; m++)
              {
                /* Distribution 1 */
                
                if (val1 != NULL)
                  {
                    if (val1[n][m] < min1)
                      val1[n][m] = 0.0;
                    else if (val1[n][m] > max1)
                      val1[n][m] = 1.0;
                    else
                      val1[n][m] = (val1[n][m] - min1)/(max1 - min1);
                  }
                
                /* Distribution 2 */
                
                if (val2 != NULL)
                  {
                    if (val2[n][m] < min2)
                      val2[n][m] = 0.0;
                    else if (val2[n][m] > max2)
                      val2[n][m] = 1.0;
                    else
                      val2[n][m] = (val2[n][m] - min2)/(max2 - min2);
                  }
              }

          /*******************************************************************/
        }
      else
        {
          /*******************************************************************/

          /***** Log color scheme ********************************************/

          /* Reset minima and maxima */

          min1 = INFTY;
          max1 = -INFTY;
          
          /* Loop over distributions */

          for (n = 0; n < sx; n++)
            for(m = 0; m < sy; m++)
              {
                /* Check zero */
                
                if (val1[n][m] > 0.0)
                  {
                    /* Take log */

                    val1[n][m] = log(val1[n][m]);
                    
                    /* Compare to minimum and maximum */
                    
                    if (val1[n][m] < min1)
                      min1 = val1[n][m];
                    if (val1[n][m] > max1)
                      max1 = val1[n][m];
                  }
                else
                  val1[n][m] = INFTY;
              }

          /* Reconfigure zeros */

          for (n = 0; n < sx; n++)
            for(m = 0; m < sy; m++)
              if (val1[n][m] == INFTY)
                val1[n][m] = min1;
          
          /* Renormalize */

          for (n = 0; n < sx; n++)
            for(m = 0; m < sy; m++)
              val1[n][m] = (val1[n][m] - min1)/(max1 - min1) + ZERO;

          /*******************************************************************/
        }
      
      /***********************************************************************/
      
      /***** Create and draw image *******************************************/

      /* Create the image */
      
      if ((im = gdImageCreate(sx, sy)) == NULL)
        Die(FUNCTION_NAME, "Unable to create image for mesh plot");
      
      /* Generate palette */
      
      MakePalette(R, G, B, COLOURS*2 + 2, (long)RDB[mpl + MPL_COLMAP]);

      palette[0] = gdImageColorAllocate(im, 0, 0, 0);  
      
      for (n = 1; n < COLOURS*2 + 2; n++)
        palette[n] = gdImageColorAllocate(im, R[n], G[n], B[n]);

      /* Draw image */
      
      for (n = 0; n < sx; n++)
        for (m = 0; m < sy; m++)
          {
            /* Check number of distributions */

            if (ndis == 2)
              {
                /* Try first */
                
                tmpv1 = val1[n][m];
                
                if (tmpv1 > 0.0)
                  i = (long)((COLOURS - 1)*tmpv1) + COLOURS + 1;
                else
                  {
                    /* Zero value --> try second */
                    
                    tmpv2 = val2[n][m];
                    
                    if (tmpv2 > 0)
                      i = (long)((COLOURS - 1)*tmpv2) + 1;
                    else
                      i = 0;
                  }
              }
            else
              {
                /* Single distribution */
                
                tmpv1 = val1[n][m];
                
                if (tmpv1 > 0.0)
                  i = (long)((2*COLOURS - 1)*tmpv1 + 1);
                else
                  i = 0;
              }
            
            /* Set pixel */

            gdImageSetPixel(im, n, sy - m - 1, palette[i]);
          }
      
      /* Write image (png format) */
      
      gdImagePng(im, fp);

      /***********************************************************************/

      /***** Free memory etc. ************************************************/
      
      /* Free memory */
      
      gdImageDestroy(im);
      
      /* Close file */
      
      fclose(fp);

      /* Free distribution matrices */
      
      if (val1 != NULL)
        {
          for(n = 0; n < sx; n++)
            Mem(MEM_FREE, val1[n]);
          
          Mem(MEM_FREE, val1);
        }

      if (val2 != NULL)
        {
          for(n = 0; n < sx; n++)
            Mem(MEM_FREE, val2[n]);
          
          Mem(MEM_FREE, val2);
        }

      /* Next distribution */

      mpl = NextItem(mpl);

      /***********************************************************************/
    }

#endif
}

/*****************************************************************************/
