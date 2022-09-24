/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readplasmasrc.c                                */
/*                                                                           */
/* Created:       2014/03/04 (PSi)                                           */
/* Last modified: 2017/11/15 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Reads data for fusion plasma source                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadPlasmaSrc:"

/*****************************************************************************/

void ReadPlasmaSrc(long src)
{
#ifdef mmmmmmmmmmmmmmmmmmmmmmm
  long loc0, ptr_xb, ptr_yb, ptr_phi, loc1, ptr_rea_e, ptr_rea_p_e;
  long ptr_rea_rho, ptr_rea_p_rho, ptr_rea_p_xy;
  long loc2, N_iso, ptr_rea_phi, ptr_p_iso;
  long ptr_rea_xb, ptr_rea_yb, N_e, N_rz;
  long N_xyb, N_rho, N_rea_xyb, N_rea_phi, N_react, ptr_rea_p_phi;
  long ptr_rea_r, ptr_rea_z, ptr_rea_p_rz;
  long i ,n, m, l, co;
  double sum, max, X0, Y0, p, rea, T, phi, x, y, E, r, z, rho;
  double Ntot;
  FILE *fp;

  /* Check source pointer */
  
  CheckPointer(FUNCTION_NAME, "(src)", DATA_ARRAY, src);

  /* Check that file name is given in src card */

  if ((long)RDB[src + SRC_READ_PTR_FILE] < VALID_PTR)
    Die(FUNCTION_NAME, "File name not given");
  
  /* Open file for reading */
      
  if ((fp = fopen(GetText(src + SRC_READ_PTR_FILE), "r")) == NULL)
    Error(src, "Unable to open plasma source file \"%s\"", 
          GetText(src + SRC_READ_PTR_FILE));
  
  /* Create plasma source structure */
  
  /* General data (geometry) */

  loc0 = NewItem(src + SRC_PTR_PLASMA_SRC, SRC_PLASMA_BLOCK_SIZE); 
 
  /*coordinate system*/
  /*1: Rz or 2: RhoTheta*/
  
  if (fscanf(fp, "%ld", &co) == EOF)
    Error(src, "Missing co");
  
  WDB[loc0 + SRC_PLASMA_CO] = (double)co;
 
  /* Read general data X0, Y0, N_xyb, boundary x,y */
  /* (if RhoTheta system is used)*/
  
  printf("coordinate system ");
  printf("%ld\n", co); 

  switch (co)
    {
    case 1:
      {
        /* Total amount of neutrons*/  
        if (fscanf(fp, "%lf", &Ntot) == EOF)
          Error(src, "Missing Ntot");
        WDB[loc0 + SRC_PLASMA_NTOT] = (double)Ntot;  
        break;
      }
    case 2:
      {
        /* Total amount of neutrons*/          
        if (fscanf(fp, "%lf", &Ntot) == EOF)
          Error(src, "Missing Ntot");
        WDB[loc0 + SRC_PLASMA_NTOT] = (double)Ntot;  
        
        /*printf("total amount of n ");
          printf("%E\n", Ntot); */
        
        /*magnetic axis*/
        if (fscanf(fp, "%lf", &X0) == EOF)
          Error(src, "Missing X0");
        WDB[loc0 + SRC_PLASMA_X0] = (double)X0;
        
        /*printf("x0 ");
          printf("%E\n", X0); */
  
        if (fscanf(fp, "%lf", &Y0) == EOF)
          Error(src, "Missing Y0");
        WDB[loc0 + SRC_PLASMA_Y0] = (double)Y0;
        
        /*printf("y0 ");
          printf("%E\n", Y0); 
          exit(0); */
         
        /*size of boundary vector*/        

        if (fscanf(fp, "%ld", &N_xyb) == EOF)
          Error(src, "Missing number of values");

        if ((N_xyb < 1) || (N_xyb > 1000))
          Error(src, "Invalid size %ld for plasma boundary vector", N_xyb);
      
        ptr_xb = ReallocMem(DATA_ARRAY, N_xyb);
        ptr_yb = ReallocMem(DATA_ARRAY, N_xyb);
        
 
        WDB[loc0 + SRC_PLASMA_XYB] = (double)N_xyb;
        WDB[loc0 + SRC_PLASMA_PTR_XB] = (double)ptr_xb;
        WDB[loc0 + SRC_PLASMA_PTR_YB] = (double)ptr_yb;
        
        printf("N_xyb ");
        printf("%ld\n", N_xyb);
        /*exit(0); */ 
        /* boundary x */

        for (n = 0; n < N_xyb; n++)
          {
            if (fscanf(fp, "%lf", &x) == EOF)
              Error(src, "Unexpected end-of-file 1");

            /* ITER: ~15m kokoluokka riittää koordinaateille*/

            CheckValue(FUNCTION_NAME, "x", "", x, -15, 15);
            WDB[ptr_xb + n] = x;
          }

        /* printf("x boundary coordinate \n");
           for (n = 0; n < N_xyb; n++)
           printf("%ld %E\n", n, RDB[ptr_xb + n]); 
           exit(0); */

        /* boundary y */

        for (n = 0; n < N_xyb; n++)
          {
            if (fscanf(fp, "%lf", &y) == EOF)
              Error(src, "Unexpected end-of-file 2");
            

            CheckValue(FUNCTION_NAME, "y", "", y, -15, 15);
            WDB[ptr_yb + n] = y;
          }

        /* Break case */
        /*printf("y boundary coordinate \n");
          for (n = 0; n < N_xyb; n++)
          printf("%ld %E\n", n, RDB[ptr_yb + n]); */
        break;
      }
    }
  
  /* Read isotropic distributions for using isotropic theta */
  /* and phi distribution*/ 

  if (fscanf(fp, "%ld", &N_iso) == EOF) 
    Error(src, "Missing number of values");

  /*printf("N_xyb ");
    printf("%ld\n", N_xyb); */
   
  CheckValue(FUNCTION_NAME, "y", "", N_iso, 0, 15000);
  WDB[loc0 + SRC_PLASMA_N_ISO] = (double)N_iso;

  /*printf("N_iso ");
    printf("%ld\n", N_iso);*/
  
  ptr_phi = ReallocMem(DATA_ARRAY, N_iso); 
  WDB[loc0 + SRC_PLASMA_PTR_PHI] = (double)ptr_phi;

  ptr_p_iso = ReallocMem(DATA_ARRAY, N_iso + 1);
  WDB[loc0 + SRC_PLASMA_PTR_P_ISO] = (double)ptr_p_iso;
  
  for (n = 0; n < N_iso; n++)
    {
      if (fscanf(fp, "%lf", &phi) == EOF)
        Error(src, "Unexpected end-of-file 3");

      CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 2*PI);
      WDB[ptr_phi + n] = phi;
    }
        
  /*printf("phi \n");
    for (n = 0; n < N_iso; n++)
    printf("%ld %E\n", n, RDB[ptr_phi + n]); 
    exit(0); */
  
  for (n = 0; n < N_iso; n++)
    {
      if (fscanf(fp, "%lf", &p) == EOF)
        Error(src, "Unexpected end-of-file 4");
           
      CheckValue(FUNCTION_NAME, "p", "", p, 0.0, 1.0);
      WDB[ptr_p_iso + n + 1] = p;
    } 

  
  /* Calculate sum for normalization*/

  sum = 0.0;
  for (n = 0; n < N_iso + 1; n++)
    sum = sum + RDB[ptr_p_iso + n];
 
  CheckValue(FUNCTION_NAME, "sum", "", sum, ZERO, INFTY);
 
  /* Calculate normalized cumulative distribution */

  WDB[ptr_p_iso] = RDB[ptr_p_iso]/sum;

  for (n = 1; n < N_iso + 1; n++)
    WDB[ptr_p_iso + n] = (RDB[ptr_p_iso + n])/sum + RDB[ptr_p_iso + n - 1]; 
    
  /* Last value */

  WDB[ptr_p_iso + N_iso] = 1.0;
  
  /*printf("p_iso \n");
    for (n = 0; n < N_iso+1; n++)
    printf("%ld %E\n", n, RDB[ptr_p_iso + n]); 
    exit(0); */
        
 
  /* Read number of reactions*/
  
  if (fscanf(fp, "%ld", &N_react) == EOF)
    Error(src, "Missing number of values");

  if ((N_react < 1) || (N_react > 10))
    Error(src, "Invalid size %ld for reaction distribution", N_react);

  WDB[loc0 + SRC_PLASMA_REACTIONS] = (double)N_react;
  
  /*printf("N_react ");
    printf("%ld\n", N_react); 
    exit(0); */

  /* Reaction data*/
 
  for(m = 0; m < N_react; m++)          
    {
      /* Create reaction data structure */        
      printf("m = %ld %ld\n", m, N_react);
      loc1 = NewItem(loc0 + SRC_PLASMA_PTR_REACT, SRC_PLASMA_REA_BLOCK_SIZE);
          
      /* Reaction */

      if (fscanf(fp, "%lf", &rea) == EOF)    
        Error(src, "Unexpected end-of-file 5");

      CheckValue(FUNCTION_NAME, "rea", "", rea, 0.0, 1.0);
      WDB[loc1 + SRC_PLASMA_REA_PROB]=(double)rea; 
          
      printf("reac p ");
      printf("%E\n", rea); 
          
      /* Time */

      if (fscanf(fp, "%lf", &T) == EOF)
        Error(src, "Unexpected end-of-file 6");
          
      CheckValue(FUNCTION_NAME, "T", "", T, 0.0, 20.0);
      WDB[loc1 + SRC_PLASMA_REA_T]=(double)T; 
          
      /* Energy vector length*/

      if (fscanf(fp, "%ld", &N_e) == EOF) 
        Error(src, "Missing number of values");
      else if ((N_e < 1) || (N_e > 1000))
        Error(src, "Invalid size %ld for energy distribution", N_e);

      /*number of rz-grids*/

       if (fscanf(fp, "%ld", &N_rz) == EOF) 
             Error(src, "Missing number of values");
         
       if ((N_rz < 1) || (N_rz > 10000))
             Error(src, "Invalid size %ld for R vector", N_rz);


      CheckValue(FUNCTION_NAME, "N_e", "", T, 0, 100.0);
      WDB[loc1 + SRC_PLASMA_REA_N_E ] = (double)N_e;
        
      ptr_rea_e = ReallocMem(DATA_ARRAY, N_e);
      WDB[loc1 + SRC_PLASMA_REA_PTR_E ] = (double)ptr_rea_e;
        
      ptr_rea_p_e = ReallocMem(DATA_ARRAY, N_rz*(N_e + 1));
      WDB[loc1 + SRC_PLASMA_REA_PTR_P_E] = (double)ptr_rea_p_e;

      p = 0.0;

      /* Read energy vector values */
      printf("N_e = %ld\n", N_e);
      for (n = 0; n < N_e; n++)
        {
          if (fscanf(fp, "%lf", &E) == EOF)
            Error(src, "Unexpected end-of-file 7"); 

          CheckValue(FUNCTION_NAME, "E", "", E, 1.0, 20.0); 
          WDB[ptr_rea_e + n] = E;
        }

      printf("E ");
      printf("%E\n", E); 
        


      /* Read P_e values to every Rz grid from the input file */  

        for(i = 0; i<N_rz; i++){
              for (n = 1; n < N_e+1; n++){
                  if (fscanf(fp, "%lf", &p) == EOF)
                    Error(src, "Unexpected end-of-file 8"); 

                  CheckValue(FUNCTION_NAME, "p", "", p, 0.0, 1.0); 
                  WDB[ptr_rea_p_e + i*(N_e+1) + n] = p; 
                }

        /* Calculate cummulative distrFibutin*/ 

      sum = 0.0;
      WDB[ptr_rea_p_e + (N_e+1)*i] = sum;
      for (l = 1; l < N_e+1; l++){
        sum = sum + RDB[ptr_rea_p_e + i*(N_e+1) + l];
        WDB[ptr_rea_p_e + (N_e+1)*i + l] = sum;
        }

        CheckValue(FUNCTION_NAME, "sum", "", sum, ZERO, INFTY);


        /*normalization of cummulative distribution*/
        if(sum == 0){                                        
                for (l = 0; l < N_e+1; l++)
                        WDB[ptr_rea_p_e + (N_e+1)*i+l] = 0;
        }else{
                for (l = 0; l < N_e+1; l++)
                        WDB[ptr_rea_p_e + (N_e+1)*i+l] = RDB[ptr_rea_p_e + i*(N_e+1) + l]/sum;
        }

                
        }


        /* Geometrical distributions */   
      /* radial p */

      /* Pointer to geometrical distribution data */
    
      loc2 = NewItem(loc1 + SRC_PLASMA_REA_PTR_GEOM, 
                     SRC_PLASMA_REA_GEOM_BLOCK_SIZE); 

      switch (co) 
        {
        case 1:
          {
            /*****************************************************************/

            /***** Rz grid is used *******************************************/

            /* size of grid */

            if (fscanf(fp, "%ld", &N_rz) == EOF) 
              Error(src, "Missing number of values");
         
            if ((N_rz < 1) || (N_rz > 10000))
              Error(src, "Invalid size %ld for R vector", N_rz);

            WDB[loc2 + SRC_PLASMA_REA_N_RZ] = (double)N_rz;
                  
                printf("N_zr=%ld", N_rz);
            /* Pointer to R */  

            ptr_rea_r = ReallocMem(DATA_ARRAY, N_rz);           
            WDB[loc2 + SRC_PLASMA_REA_PTR_R] = (double)ptr_rea_r;
          
            /* Read R values*/

            for (n = 0; n < N_rz; n++)
              {
                if (fscanf(fp, "%lf", &r) == EOF)
                  Error(src, "Unexpected end-of-file 9"); 

                /* ITER: R: [-15.0 15.0] */

                CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 15.0); 
                WDB[ptr_rea_r + n] = r;
              }

            /* Pointer to Z */

            ptr_rea_z = ReallocMem(DATA_ARRAY, N_rz);
            WDB[loc2 + SRC_PLASMA_REA_PTR_Z] = (double)ptr_rea_z;
          
            /* Read Z values*/

            for (n = 0; n < N_rz; n++)
              {
                if (fscanf(fp, "%lf", &z) == EOF)
                  Error(src, "Unexpected end-of-file 10"); 

                /* ITER: Z: [-15.0 15.0] */

                CheckValue(FUNCTION_NAME, "z", "", z, -10.0, 15.0); 
                WDB[ptr_rea_z + n] = z;
              }
        
            /*P_RZ*/

            /* Pointer to P_Rz*/

            ptr_rea_p_rz = ReallocMem(DATA_ARRAY, N_rz + 1);
            WDB[loc2 + SRC_PLASMA_REA_PTR_P_RZ] = (double)ptr_rea_p_rz;    

            /* Read P_RZ from input file*/
            
            for (n = 0; n < N_rz; n++)
              {
                if (fscanf(fp, "%lf", &p) == EOF)
                  Error(src, "Unexpected end-of-file 11");

                CheckValue(FUNCTION_NAME, "p", "", p, 0.0, 1.0);
                WDB[ptr_rea_p_rz + n + 1] = p;   
              }

            /* Calculate maximum for normalization*/

            max = 0.0; 
            for (n = 0; n < N_rz + 1 ; n++)
              if (max < RDB[ptr_rea_p_rz + n])
                max = RDB[ptr_rea_p_rz + n];

            /* JLe 15/06/24: Tähän jotkut järkevät rajat */

            CheckValue(FUNCTION_NAME, "max", "", max, ZERO, INFTY);
            
            /* Calculate normalized cumulative distribution */

            for (n = 0; n < N_rz + 1; n++)
              WDB[ptr_rea_p_rz + n] = (RDB[ptr_rea_p_rz + n])/max;
                
            /* Last value forced to 1.0 */

            WDB[ptr_rea_p_rz + N_rz] = 1.0;

            /* Break case */

            break;
       
            /*****************************************************************/
          }
        case 2:
          {
            /*****************************************************************/

            /***** RhoTheta grid is used *************************************/

            /* size of rho grid */          

            if (fscanf(fp, "%ld", &N_rho) == EOF) 
              Error(src, "Missing number of values");
            printf("N_rho = %ld\n", N_rho);
            if ((N_rho < 1) || (N_rho > 1000))
              Error(src, "Invalid size %ld for rho distribution", N_rho);

            WDB[loc2 + SRC_PLASMA_REA_N_RHO] = (double)N_rho;
                
            /*printf("N_rho ");
              printf("%ld\n", N_rho); */
         
            /* Pointer to Rho and P_rho*/

            ptr_rea_rho = ReallocMem(DATA_ARRAY, N_rho);
            WDB[loc2 + SRC_PLASMA_REA_PTR_RHO] = (double)ptr_rea_rho;
            
            ptr_rea_p_rho = ReallocMem(DATA_ARRAY, N_rho + 1);
            WDB[loc2 + SRC_PLASMA_REA_PTR_P_RHO] = (double)ptr_rea_p_rho;
         
            for (n = 0; n < N_rho; n++)
              {
                if (fscanf(fp, "%lf", &rho) == EOF)
                  Error(src, "Unexpected end-of-file 12");

                /* rho: [0, 1] */

                CheckValue(FUNCTION_NAME, "rho", "", rho, 0.0, 1.01);
                WDB[ptr_rea_rho + n] = rho;
              }
                  
            /*printf("rho values \n");
              for (n = 0; n < N_rho; n++)
              printf("%ld %E\n", n, RDB[ptr_rea_rho + n]);   */
  
            for (n = 0; n < N_rho; n++)
              {
                if (fscanf(fp, "%lf", &p) == EOF)
                  Error(src, "Unexpected end-of-file 13");
        
                CheckValue(FUNCTION_NAME, "p", "", p, 0.0, 1.0);
                WDB[ptr_rea_p_rho + n + 1] = p;
              }
                  
            /*printf("p_rho values \n");
              for (n = 0; n < N_rho; n++)
              printf("%ld %E\n", n, RDB[ptr_rea_p_rho + n]); */
         
            /* Sum for normalization*/

            /* sum = 0.0;
               for (n = 0; n < N_rho + 1 ; n++)
               sum = sum + RDB[ptr_rea_p_rho + n];
            */
               /* Calculate normalized cumulative distribution */
            
            /*WDB[ptr_rea_p_rho] = RDB[ptr_rea_p_rho]/sum;
              for (n = 1; n < N_rho + 1; n++)
              WDB[ptr_rea_p_rho + n] = (RDB[ptr_rea_p_rho + n])/sum + 
              RDB[ptr_rea_p_rho + n - 1];*/

            /* Last value */

            WDB[ptr_rea_p_rho + N_rho]=1.0;
     
            /* size of boundary grid (poloidal xy)*/         

            if (fscanf(fp, "%ld", &N_rea_xyb) == EOF)
              Error(src, "Missing number of values");
         
            CheckValue(FUNCTION_NAME, "N_rea_xyb", "", N_rea_xyb, 0, 1000);
            printf("N_rea_xyb = %ld\n", N_rea_xyb);
            if (N_rea_xyb > 0)
              {         
                /* Anisotropic P_xy was given*/

                WDB[loc2 + SRC_PLASMA_REA_XYB] = (double)N_rea_xyb;

                /* Pointers to the boundary (x,y) vectors an P_xy*/

                ptr_rea_xb = ReallocMem(DATA_ARRAY, N_xyb);
                WDB[loc2 + SRC_PLASMA_REA_PTR_XB] = (double)ptr_rea_xb;

                ptr_rea_yb = ReallocMem(DATA_ARRAY, N_xyb);
                WDB[loc2 + SRC_PLASMA_REA_PTR_YB] = (double)ptr_rea_yb;

                ptr_rea_p_xy = ReallocMem(DATA_ARRAY, N_xyb + 1);
                WDB[loc2 + SRC_PLASMA_REA_PTR_P_XY] = (double)ptr_rea_p_xy;
                       
                /* Boundary x coordinates */

                for (n = 0; n < N_rea_xyb; n++)
                  {
                    if (fscanf(fp, "%lf", &x) == EOF)
                      Error(src, "Unexpected end-of-file 14");
                    
                    CheckValue(FUNCTION_NAME, "x", "", x, -10, 10);
                    WDB[ptr_rea_xb + n] = x;
                  }

                /* Boundary y coordinates */

                for (n = 0; n < N_rea_xyb; n++)
                  {
                    if (fscanf(fp, "%lf", &y) == EOF)
                      Error(src, "Unexpected end-of-file 15");
                    
                    CheckValue(FUNCTION_NAME, "y", "", y, -10, 10);
                    WDB[ptr_rea_yb + n] = y;
                  } 

                /* P_xy distribution */

                for (n = 0; n < N_rea_xyb+1; n++)
                  {          
                    if (fscanf(fp, "%lf", &p) == EOF)
                      Error(src, "Unexpected end-of-file 16");
                    
                    CheckValue(FUNCTION_NAME, "p", "", p, 0.0, 1.0);
                    WDB[ptr_rea_p_xy + n + 1] = p;
                  }

                /* Calculate sum for normalization*/
                    
                sum = 0.0;
                for (n = 0; n < N_xyb + 1; n++)
                  sum = sum + RDB[ptr_rea_p_xy + n];

                /* JLe 15/06/24: Tähän jotkut järkevät rajat */

                CheckValue(FUNCTION_NAME, "sum", "", sum, ZERO, INFTY);
                                            
                /* Calculate normalized cumulative distribution */

                WDB[ptr_rea_p_xy] = RDB[ptr_rea_p_xy]/sum;
                for (n = 1; n < N_rea_xyb + 1; n++)
                  WDB[ptr_rea_p_xy + n] = (RDB[ptr_rea_p_xy + n])/sum + 
                    RDB[ptr_rea_p_xy + n - 1]; 
                
                /* Last value */
                
                WDB[ptr_rea_p_xy + N_rea_xyb]=1.0;
              }
            
            /* Size (N_rea_xy) = 0 -> isotropic (general data), x, y, */
            /* P_xy are used*/  
          
            else
              {

                WDB[loc2 + SRC_PLASMA_REA_XYB] = (double)N_xyb;
                WDB[loc2 + SRC_PLASMA_REA_PTR_XB] = (double)ptr_xb; 
                WDB[loc2 + SRC_PLASMA_REA_PTR_YB] = (double)ptr_yb;
                WDB[loc2 + SRC_PLASMA_REA_PTR_P_XY] = (double)ptr_p_iso;
                 
                /* printf("p_xyb values \n");
                   for (n = 0; n < N_xyb; n++)
                   printf("%ld %E\n", n, RDB[ptr_p_iso + n]); 
                 
                   }*/
              }
                  
            /* Toroidal distribution P_phi*/

            if (fscanf(fp, "%ld", &N_rea_phi) == EOF) 
              Error(src, "Missing number of values");
            
            printf("N_rea_phi = %ld\n", N_rea_phi);
            
            if (N_rea_phi > 0)
              { 
                /* Anisotropic P_phi was given*/
                
                CheckValue(FUNCTION_NAME, "N_rea_z", "", N_rea_phi, 0, 1000);
                WDB[loc2 + SRC_PLASMA_REA_N_PHI] = (double)N_rea_phi;
                
                /* Pointers to phi and P_phi*/
                
                ptr_rea_phi = ReallocMem(DATA_ARRAY, N_rea_phi);
                WDB[loc2 + SRC_PLASMA_REA_PTR_PHI] = (double)ptr_rea_phi;
                
                ptr_rea_p_phi = ReallocMem(DATA_ARRAY, N_rea_phi + 1);
                WDB[loc2 + SRC_PLASMA_REA_PTR_P_PHI] = (double)ptr_rea_p_phi;
                
                /* Read phi*/ 
                
                for (n = 0; n < N_rea_phi; n++)
                  {
                    if (fscanf(fp, "%lf", &phi) == EOF)
                      Error(src, "Unexpected end-of-file 17");
                    
                    CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 2*PI);
                    WDB[ptr_rea_phi + n] = phi;
                  }
                
                /* Read P_phi */ 
                
                for (n = 0; n < N_rea_phi; n++)
                  {
                    if (fscanf(fp, "%lf", &p) == EOF)
                      Error(src, "Unexpected end-of-file 18");
                    
                    CheckValue(FUNCTION_NAME, "p", "", p, 0.0, 1.0);
                    WDB[ptr_rea_p_phi + n + 1] = p;
                  }
                
                /* Calculate sum for normalization*/
                
                sum = 0.0;
                for (n = 0; n < N_rea_phi + 1; n++)
                  sum = sum + RDB[ptr_rea_p_phi + n];
                
                /* JLe 15/06/24: Tähän jotkut järkevät rajat */
                
                CheckValue(FUNCTION_NAME, "sum", "", sum, ZERO, INFTY);
                
                /* Calculate normalized cumulative distribution */
                
                WDB[ptr_rea_p_phi] = RDB[ptr_rea_p_phi]/sum;
                for (n = 1; n < N_rea_phi + 1; n++)
                  WDB[ptr_rea_p_phi + n] = (RDB[ptr_rea_p_phi + n])/sum + 
                    RDB[ptr_rea_p_phi + n - 1]; 
                
                /* Last value*/
                
                WDB[ptr_rea_p_phi + N_rea_phi]=1.0;
              }
            else
              {
                /* size = 0 -> isotropic (general data) phi, P_phi are used*/
                
                WDB[loc2 + SRC_PLASMA_REA_N_PHI] = (double)N_iso; 
                WDB[loc2 + SRC_PLASMA_REA_PTR_PHI] = (double)ptr_phi;
                WDB[loc2 + SRC_PLASMA_REA_PTR_P_PHI] = (double)ptr_p_iso;
                
              }
          }
        }
    }

#endif
  /*****************************************************************************/
}
