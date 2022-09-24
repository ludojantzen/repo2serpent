/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sampleplasmasrc.c                              */
/*                                                                           */
/* Created:       2014/03/04 (PSi)                                           */
/* Last modified: 2017/11/15 (PSi)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Samples source point from plasma source                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SamplePlasmaSrc:"

/*****************************************************************************/

void SamplePlasmaSrc(long src, double *x, double *y, double *z, double *u,  
                     double *v, double *w, double *E, double *wgt, double *t,
                     long id)
{
#ifdef mmmmmmmmmmmmmmmmmmm

  long loc0, loc1, loc2, ptr_p_react, ptr_r, ptr_z, ptr_p_rz, ptr_p_xy; 
  long ptr_e, ptr_rho, ptr_xb, ptr_yb, ptr_phi;
  long ptr_p_e, ptr_p_rho, ptr_p_phi;
  long N_react, N_rho, N_rz, N_rea_xyb, N_e, N_phi;
  long i2, i3, i4, i5, co;
  double rnd1, rnd2, rnd3, rnd4, rnd5, rho, X0, Y0, phi, Rp, Zp, R, Z ,dens;

  /* Check source pointer */

  CheckPointer(FUNCTION_NAME, "(src)", DATA_ARRAY, src);

  /* Get pointer to plasma source structure */

  loc0 = (long)RDB[src + SRC_PTR_PLASMA_SRC];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* JLe 15/06/24: Voisiko nämä rivit poistaa kun loc1 ja loc2   */
  /* luetaan tuolla myöhemmin, eikä niitä käytetä tässä välissä? */
  
  
  loc2 = (long)RDB[loc1 + SRC_PLASMA_REA_PTR_GEOM];
  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);
  
  /* General data */
  
  /* Read sizes */

  N_react = (long)RDB[loc0 + SRC_PLASMA_REACTIONS];
  CheckValue(FUNCTION_NAME, "N_react", "", N_react, 1, 10);

  N_phi = (long)RDB[loc0 + SRC_PLASMA_N_ISO];
  CheckValue(FUNCTION_NAME, "N_phi", "", N_phi, 1, 1000); 
  
  /* JLe 15/06/24: Voisiko tämän poistaa kun sitä ei käytetä ollenkaan? */

  ptr_p_react = (long)RDB[loc0 + SRC_PLASMA_PTR_REACT];
  CheckPointer(FUNCTION_NAME, "(ptr_p_react)", DATA_ARRAY, ptr_p_react);

  /* Avoid compiler warning */

  X0 = 0.0;
  Y0 = 0.0;
  R = 0.0;
  Z = 0.0;
  rho = 0.0;

  ptr_xb = -1;
  ptr_yb = -1;
  
  
  /* Coordinate system*/

  co = (long)RDB[loc0 + SRC_PLASMA_CO];
  /*printf("coordinate system ");
  printf("%ld", co); */
 
  switch (co)
    {
    case 1:
      { 
        /* Break case */

        break;
      }
    case 2:
      {
        /* Read origin  */

        X0 = (double)RDB[loc0 + SRC_PLASMA_X0];
        Y0 = (double)RDB[loc0 + SRC_PLASMA_Y0];
        /*printf("X0 ");
    printf("%E\n", X0);
        printf("Y0 ");
        printf("%E\n", Y0);*/
        /* Break case */
    /*exit(0);*/
        break;
      }
    default:
      {
        /* Invalid mode */

        Die(FUNCTION_NAME, "co = %ld\n", co);
      }
    }
  
  /* Reset neutron coordinates */

  *x = 0.0;
  *y = 0.0;
  *z = 0.0;
  
  phi = 0.0;
  Rp = 0.0;
  Zp = 0.0;

  /* JLe: Se voisi ehkä selkeyttää tätä rakennetta jos nuo satunnaislukujen */
  /* arvonnat siirtäisi sinne missä niitä käytetään. */
  
  /* Random numbers */ 

  /* reaction */

  rnd1 = RandF(id);

  /* rho or Rz coordinate*/

  rnd2 = RandF(id);
 
  /* poloidal plane (xy) if it is needed */

  rnd3 = RandF(id);
 
  /* toroidal angle (phi) */

  rnd4 = RandF(id);

  /* energy */

  rnd5 = RandF(id);

 
  /* Pointer to 1. data structure */

  loc1 = (long)RDB[loc0 + SRC_PLASMA_PTR_REACT];
  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
 
  while (NextItem(loc1) > VALID_PTR)
    {
      if (rnd1 < RDB[loc1 + SRC_PLASMA_REA_PROB])
        break;
        
      loc1 = NextItem(loc1);
    }
  
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

  loc2 = (long)RDB[loc1 + SRC_PLASMA_REA_PTR_GEOM];
  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

  switch (co)
    {
    case 1:
      {
        N_rz = (long)RDB[loc2 + SRC_PLASMA_REA_N_RZ]; 
        CheckValue(FUNCTION_NAME, "N_rz", "", N_rz, 1, 10000);
        
        ptr_r = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_R];
        CheckPointer(FUNCTION_NAME, "(ptr_r)", DATA_ARRAY, ptr_r);
 
        ptr_z = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_Z];
        CheckPointer(FUNCTION_NAME, "(ptr_z)", DATA_ARRAY, ptr_z);
                
        ptr_p_rz = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_P_RZ];
        CheckPointer(FUNCTION_NAME, "(ptr_p_rz)", DATA_ARRAY, ptr_p_rz);

        /*checking if the density in the sampled grid is zero*/
        N_e = (long)RDB[loc1 + SRC_PLASMA_REA_N_E];
          CheckValue(FUNCTION_NAME, "N_e", "", N_e, 1, 1000); 
          ptr_p_e = (long)RDB[loc1 + SRC_PLASMA_REA_PTR_P_E];
          CheckPointer(FUNCTION_NAME, "(ptr_p_e)", DATA_ARRAY, ptr_p_e);

        dens = 0;
        while(dens == 0){
                /* Rz index (i2)*/   
                i2 = SearchArray(&RDB[ptr_p_rz], rnd2, N_rz + 1);
                CheckValue(FUNCTION_NAME, "i2", "", i2, 0, N_rz - 1);
                dens = RDB[ptr_p_e+ (N_e+1)*i2+N_e-1]; 
                rnd2 = RandF(id);
        }

        /* ITER suuruusluokka ~15m riittää*/

        R = RDB[ptr_r + i2];
        CheckValue(FUNCTION_NAME, "R", "", R, -15.00, 15.00);

        Z = RDB[ptr_z + i2];
        CheckValue(FUNCTION_NAME, "Z", "", Z, -15.00, 15.00);

        break;
      } 
    case 2:
      {
        /* JLe 15/12/26: Tässä ja seuraavassa loc1 --> loc2 */
        
        N_rho = (long)RDB[loc2 + SRC_PLASMA_REA_N_RHO];
        CheckValue(FUNCTION_NAME, "N_rho", "", N_rho, 1, 1000);

        N_rea_xyb = (long)RDB[loc2 + SRC_PLASMA_REA_XYB];
        CheckValue(FUNCTION_NAME, "N_rea_xyb", "", N_rea_xyb, 1, 1000);
                
        /* JLe 15/12/26: PTR_R --> PTR_RHO */
        
        ptr_rho = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_RHO];
        CheckPointer(FUNCTION_NAME, "(ptr_rho)", DATA_ARRAY, ptr_rho);
                
        ptr_p_rho = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_P_RHO];
        CheckPointer(FUNCTION_NAME, "(ptr_p_rho)", DATA_ARRAY, ptr_p_rho);        
                        
        /* rho index (i2)*/
        
        i2 = SearchArray(&RDB[ptr_p_rho], rnd2, N_rho + 1);
        CheckValue(FUNCTION_NAME, "i2", "", i2, 0, N_rho - 1);
        
        ptr_xb = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_XB];
        CheckPointer(FUNCTION_NAME, "(ptr_xb)", DATA_ARRAY, ptr_xb);

        ptr_yb = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_YB];
        CheckPointer(FUNCTION_NAME, "(ptr_yb)", DATA_ARRAY, ptr_yb);
        
        /*
        exit(0);        
        */
        
        ptr_p_xy = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_P_XY]; 
        CheckPointer(FUNCTION_NAME, "(ptr_p_xy)", DATA_ARRAY, ptr_p_xy);
 
        /* boundary index (i3)*/

        i3 = SearchArray(&RDB[ptr_p_xy], rnd3, N_rea_xyb + 1);
        CheckValue(FUNCTION_NAME, "i3", "", i3, 0, N_rea_xyb - 1);
        


        rho = RDB[ptr_rho + i2];
                
        *x = RDB[ptr_xb + i3];
        Rp = *x*rho + X0;
        *y = RDB[ptr_yb + i3];
        Zp = *y*rho + Y0;
        
        break;
      }
    }

  /* JLE 15/06/24: Näitä pointtereita ym. voisi ehkä siirrellä sinne */
  /* missä niitä käytetään ensimmäistä kertaa. */
 
  N_phi = (long)RDB[loc2 + SRC_PLASMA_REA_N_PHI];
  CheckValue(FUNCTION_NAME, "N_phi", "", N_phi, 1, 1000);

  N_e = (long)RDB[loc1 + SRC_PLASMA_REA_N_E];
  CheckValue(FUNCTION_NAME, "N_e", "", N_e, 1, 1000); 
  
  ptr_phi = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_PHI];
  CheckPointer(FUNCTION_NAME, "(ptr_phi)", DATA_ARRAY, ptr_phi);
  
  ptr_p_phi = (long)RDB[loc2 + SRC_PLASMA_REA_PTR_P_PHI];
  CheckPointer(FUNCTION_NAME, "(ptr_p_phi)", DATA_ARRAY, ptr_p_phi);
        
  ptr_p_e = (long)RDB[loc1 + SRC_PLASMA_REA_PTR_P_E];
  CheckPointer(FUNCTION_NAME, "(ptr_p_e)", DATA_ARRAY, ptr_p_e);

  ptr_e = (long)RDB[loc1 + SRC_PLASMA_REA_PTR_E];
  CheckPointer(FUNCTION_NAME, "(ptr_e)", DATA_ARRAY, ptr_e);

  
  i4 = -1;
                 
  /* toroidal angle */

  /* phi index (i4)*/
  i4 = SearchArray(&RDB[ptr_p_phi], rnd4, N_phi + 1);
  CheckValue(FUNCTION_NAME, "i4", "", i4, 0, N_phi - 1);

  phi = RDB[ptr_phi  + i4]; 
  
  /* phi = RandF(id)*2*PI; *//*
  printf("%E %E %E\n", R, phi, Z);
                             */
  switch(co)
    {
    case 1:
      {
        /* from RZPhi to Serpent xyz*/        

        *x =-R*100*cos(phi);
        *z = Z*100;
        *y = R*100*sin(phi);         

        break;
      }
    case 2:
      {
        /* from RhoThetaPhi to Serpent xyz*/

        *x = RDB[ptr_xb + i3];
        Rp = *x*rho + X0;
        *y = RDB[ptr_yb + i3];
        Zp = *y*rho + Y0; 
        *x = -Rp*100.0*sin(phi);
        *z = Zp*100.0;
        *y = Rp*100.0*cos(phi); 

        break;
      }
    }

  /* energy */

  /* JLe 15/06/24: Tää energia antaa kummallisia arvoja? Tuossa ekassa     */
  /* if-haarassa pitäis varmaan olla ptr_e, ei ptr_e + 1? Mutta silloinkin */
  /* ne arvotut arvot on 17.6 ja 3.27, eli varmaan reaktion Q-arvoja? */

  phi = RandF(id)*2.0*PI;

  Rp = sqrt(*x**x + *y**y);

  *x = Rp*cos(phi);
  *y = Rp*sin(phi);


  if (N_e == 1)
    *E = RDB[ptr_e];
  else
    {
        /*printf("%d/n", i2);*/
      i5 = SearchArray(&RDB[ptr_p_e+ (N_e+1)*i2], rnd5, N_e+1);
      CheckValue(FUNCTION_NAME, "i5", "", i5, 0, N_e-1);
      *E = RDB[ptr_e + i5];
        
        
        /*if((i2==12*50+25))
                printf("%lf ", RDB[ptr_e + i5]);
        printf("%d\t%d\n", i2%25, i2/25);
        printf(" rand5 = %f\n", rnd5);*/
    }

  /* Sample direction cosines isotropically */

  IsotropicDirection(u, v, w, id);
  /*
  printf("%E %E %E : %E\n", *x, *y, *z, *E);
  */
  /* Reset weight and time */

  *wgt = 1.0;
  *t = 0.0;

#endif
}


/*****************************************************************************/

        
