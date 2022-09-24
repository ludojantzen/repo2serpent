/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readsourcefile.c                               */
/*                                                                           */
/* Created:       2012/03/30 (JLe)                                           */
/* Last modified: 2018/08/10 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads source distribution from file                          */
/*                                                                           */
/* Comments:   -Added binary files in 2.1.25 as the limited precision in     */
/*              text files causes problems if the source points are close    */
/*              to a surface                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadSourceFile:"

/*****************************************************************************/

void ReadSourceFile(long src, double *x, double *y, double *z, double *u, 
                    double *v, double *w, double *E, double *wgt, double *t)
{
  long idx, sz, nb, ptr, pos, pos0, eof;
  double sum;
  FILE *fp;

  /* Check source pointer */

  CheckPointer(FUNCTION_NAME, "(src)", DATA_ARRAY, src);

  /* Set OpenMP barrier */

#ifdef OPEN_MP
#pragma omp critical
#endif
  {
    /* Get buffer size and number of points in buffer (pitää olla ton */
    /* omp criticalin sisällä) */

    sz = (long)RDB[src + SRC_READ_BUF_SZ];
    nb = (long)RDB[src + SRC_READ_BUF_PTS];
    
    CheckValue(FUNCTION_NAME, "nb", "", nb, 0, sz);

    /* Get index to next neutron in buffer */

    idx = (long)RDB[src + SRC_READ_BUF_IDX];

    /* Check if all particles fit in buffer */

    if ((nb > 0) && (nb < sz))
      {
        /*********************************************************************/

        /***** Adjust index without reading more points **********************/

        /* Check and reset index */

        if (idx > nb - 1)
          {
            /* Set index */

            idx = 0;

            /* Store value */

            WDB[src + SRC_READ_BUF_IDX] = (double)idx;
          }

        /*********************************************************************/
      }
    else if (idx > sz - 1)
      {
        /*********************************************************************/

        /***** Fill buffer with source points ********************************/

        /* Open file for reading */
      
        if ((fp = fopen(GetText(src + SRC_READ_PTR_FILE), "r")) == NULL)
          Error(src, "Unable to open source file \"%s\"", 
                GetText(src + SRC_READ_PTR_FILE));

        /* Get pointer to buffer */

        ptr = (long)RDB[src + SRC_READ_PTR_BUF];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
        
        /* Reset buffer index */
        
        idx = -1;

        /* Get starting position */
            
        pos0 = (long)RDB[src + SRC_READ_FILE_POS];

        /* Loop until buffer is full */
        
        do
          {
            /* Get position */
            
            pos = (long)RDB[src + SRC_READ_FILE_POS];

            /* Seek to position */
            
            fseek(fp, pos, SEEK_SET);
            
            /* Fill buffer with source points*/

            while (1 == 1)
              {
                /* Read next source point or break if cannot */

                if ((long)RDB[src + SRC_READ_BINARY] == YES)
                  {
                    /* Read BINARY file */

                    if (fread(x, sizeof(double), 1, fp) == 0)
                      break;

                    if (fread(y, sizeof(double), 1, fp) == 0)
                      break;

                    if (fread(z, sizeof(double), 1, fp) == 0)
                      break;

                    if (fread(u, sizeof(double), 1, fp) == 0)
                      break;

                    if (fread(v, sizeof(double), 1, fp) == 0)
                      break;

                    if (fread(w, sizeof(double), 1, fp) == 0)
                      break;

                    if (fread(E, sizeof(double), 1, fp) == 0)
                      break;

                    if (fread(wgt, sizeof(double), 1, fp) == 0)
                      break;

                    if (fread(t, sizeof(double), 1, fp) == 0)
                      break;
                  }
                else
                  {
                    /* Read ASCII file */

                    if ((eof = fscanf(fp, 
                                      "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                                      x, y, z, u, v, w, E, wgt, t)) == EOF)
                      break;
                  }
              
                /* Update index */
                
                idx++;
                
                /* Get position in file */
                
                pos = ftell(fp);

                /* Renormalization for type 2 */

                if ((long)RDB[src + SRC_READ_FILE_TYPE] == 
                    SRC_FILE_TYPE_S1_RENORM)
                  {
                    *wgt = 1.0;
                    *t = 0.0;
                  }
                else if ((long)RDB[src + SRC_READ_FILE_TYPE] == 
                    SRC_FILE_TYPE_WGT_RENORM)
                  {
                    *wgt = 1.0;
                  }
        
                /* Put data */
                
                WDB[ptr + SRC_BUF_X] = *x;
                WDB[ptr + SRC_BUF_Y] = *y;
                WDB[ptr + SRC_BUF_Z] = *z;
                WDB[ptr + SRC_BUF_U] = *u;
                WDB[ptr + SRC_BUF_V] = *v;
                WDB[ptr + SRC_BUF_W] = *w;
                WDB[ptr + SRC_BUF_E] = *E;
                WDB[ptr + SRC_BUF_WGT] = *wgt;
                WDB[ptr + SRC_BUF_T] = *t;
                
                /* Check if buffer is full */
                
                if (idx == sz - 1)
                  break;
              
                /* Update pointer */
                
                ptr = ptr + SRC_BUF_BLOCK_SIZE;                              
              }

            /* Check index */

            if (idx == -1)
              Error(src, "Unable to sample source point from file \"%s\"",
                    GetText(src + SRC_READ_PTR_FILE));

            /* Check for file errors */

            if (ferror(fp))
              Die(FUNCTION_NAME, "FERROR flag has been set");

            /* Put number of points read */

            nb = idx + 1;
            WDB[src + SRC_READ_BUF_PTS] = (double)nb;

            /* Check if all points were read */
            
            if ((pos0 == 0) && (feof(fp)))
              {
                /* Entire source fits in buffer */

                break;
              }
            else if (feof(fp))
              {
                /* Reached end of file, go back to beginning */
                
                WDB[src + SRC_READ_FILE_POS] = 0.0;
              }
            else
              {
                /* Remember position for next cycle */

                WDB[src + SRC_READ_FILE_POS] = (double)pos;
              }
          }
        while (idx < sz - 1);
        
        /* Reset buffer index */
        
        idx = 0;
        WDB[src + SRC_READ_BUF_IDX] = (double)idx;
        
        /* Close file */
        
        fclose(fp);
     
        /*********************************************************************/
      }
    
    /* Update buffer index */

    WDB[src + SRC_READ_BUF_IDX] += 1.0;
  }

  /***************************************************************************/

  /***** Get data from buffer ************************************************/

  /* Check index */

  if ((idx < 0) || (idx > nb - 1))
    Die(FUNCTION_NAME, "Error in buffer index");

  /* Get pointer to buffer */

  ptr = (long)RDB[src + SRC_READ_PTR_BUF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Get pointer to data */
  
  ptr = ptr + idx*SRC_BUF_BLOCK_SIZE;
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Get data */
  
  *x = RDB[ptr + SRC_BUF_X];
  *y = RDB[ptr + SRC_BUF_Y];
  *z = RDB[ptr + SRC_BUF_Z];
  *u = RDB[ptr + SRC_BUF_U];
  *v = RDB[ptr + SRC_BUF_V];
  *w = RDB[ptr + SRC_BUF_W];
  *E = RDB[ptr + SRC_BUF_E];
  *wgt = RDB[ptr + SRC_BUF_WGT];
  *t = RDB[ptr + SRC_BUF_T];

  /***************************************************************************/

  /***** Check values and find position in geometry **************************/

  /* Check energy */

  if ((*E < ZERO) || (*E > 1E+6))
    Error(src, "Invalid particle energy %E read from source file", *E);

  /* Check weight */

  if ((*wgt < ZERO) || (*wgt > 1E+6))
    Error(src, "Invalid particle weight read from source file");

  /* Normalize direction cosines */

  sum = (*u)*(*u) + (*v)*(*v) + (*w)*(*w);
  
  if (fabs(sum - 1.0) > 1E-4)
    Error(src, "Direction cosines read from source file are not normalized");
  else
    sum = sqrt(sum);

  *u = *u/sum;
  *v = *v/sum;
  *w = *w/sum;

  /***************************************************************************/
}

/*****************************************************************************/
