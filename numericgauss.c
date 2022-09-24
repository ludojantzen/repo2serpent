/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : numericgauss.c                                 */
/*                                                                           */
/* Created:       2011/05/04 (MPu)                                           */
/* Last modified: 2014/08/22 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Ratkaise (A - theta[i]*I)x = b Gaussin eliminaatiolla        */
/*              Ol: A on ccs-matriisi ja A:n symbolinen LU-hajotelma annettu */
/*                                                                           */
/* Comments: - Poistin t‰st‰ k‰ytt‰m‰ttˆm‰n muuttujan col, koska k‰‰nt‰j‰    */
/*             antaa siit‰ tarpeettoman varoituksen                          */
/*             (JLe 28.12.2012 / 2.0.12)                                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "NumericGauss:"

/*****************************************************************************/

long NumericGauss(struct ccsMatrix *cmat, complex *c, complex chv, complex *x)
{
  /* T‰m‰ funktio kirjoitettu kokonaan uusiksi crs->ccs varten */

  /* CRAM-kertoimet lis‰t‰‰n diagonaalille t‰ss‰ rutiinissa */    

  /* !! cmat = LU = A u F eli sis‰lt‰‰ paikat fill-in -alkioille !! */

  long  k,i,j,n,m,s,l;
  long *row, *first, *next, *colind; 
  complex z, z2, A_kj, A_jj, A_jl;
  complex *val, *array, *diag; 

  /*
  long *col;
  */ 

  /***************************************************************************/

  n   = cmat->n;
  m   = cmat->m;

  row     = cmat->rowind; 
  /*
  col     = cmat->colptr;
  */
  val     = cmat->values;  

  /* apumuuttuja */

  array = (complex *)Mem(MEM_ALLOC, n,sizeof(complex));

  /* tallenna diagonaalialkiot */ 

  diag  = (complex *)Mem(MEM_ALLOC, n,sizeof(complex));

  if (n != m)
    Die(FUNCTION_NAME, "Matrix not square\n");
  
  /* Tarkastetaan ett‰ arvot on reaalisia ja OK */

  for (i=0; i < cmat->nnz; i++)
    {
      CheckValue(FUNCTION_NAME, "val.re", "", val[i].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "val.im", "", val[i].im, 0.0, 0.0);
    }

  /***************************************************************************/

  /* Hae data riveitt‰in kenttiin cmat->rowptr ja cmat->next */
  /* rowptr antaa jokaisen rivin 1. alkion sijainnin */
  /* next[m] = -1 kun rivi loppuu */

  /* Miksi t‰t‰ kutsutaan t‰‰lt‰, kun riitt‰isi kai tehd‰ kerran... (16.11) */

  /* FindRowIndexes(cmat); */

  first  = cmat->rowptr; 
  next   = cmat->next;
  colind = cmat->colind;

  /*************************************************************************/

  /* Algoritmissa oletetaan, ett‰ matriisi on sortattu! */

  /* Ensimm‰iselle riville ei tehd‰ algoritmissa mit‰‰n => ota talteen 
     muuttujaan diag */


  if( row[0] != 0){ /* A_11 == 0 */
    diag[0].re = -chv.re; 
    diag[0].im = -chv.im; 
  }
  else 
    diag[0] = c_sub(val[0], chv); /* muuten val[0] = A_11*/

  /* Jos nollia diagonaalilla ... */
  for (k=1; k < n; k++){
    diag[k].re = -chv.re;
    diag[k].im = -chv.im;
  }
  

  /***************************************************************************/

  /* Algoritmi etenee riveitt‰in, aloita toiselta rivilt‰: */

  for (k=1; k < n; k++)
    { 

      /* Ota rivin k alkiot muuttujaan array ja lis‰‰ CRAM-kertoimet diagonaalille */

      m = first[k]; /* rivin k ensimm‰inen alkio */
      while(m > 0) /* rivin loppumisen merkkin‰ tulee vastaan arvo -1 */
        {
          j = colind[m]; 
          i = row[m]; 
          
          if (i != k)
            Die(FUNCTION_NAME, "Mismatch in row\n");

          
          if (j == k) /* diagonaalialkio*/
            { 
              /* CRAM-kertoimet lis‰t‰‰n t‰ss‰ */
              
              val[m].re = val[m].re - chv.re; 
              val[m].im = val[m].im - chv.im;  
              
              diag[j].re = val[m].re; 
              diag[j].im = val[m].im; 
            }  
                
          /* Ota talteen: array[j] = a_kj */
          array[j].re = val[m].re; 
          array[j].im = val[m].im;

          m = next[m]; /* seuraava t‰ll‰ rivill‰ */

        }

      /* Junaile nollat diagonaalin vas. puolelle rivill‰ k */
      m = first[k]; 
      while(m > 0)
        {
          j = colind[m]; 
          i = row[m]; 

          if (i != k)
            Die(FUNCTION_NAME, "Mismatch in row\n");

          if (j < k) /* nollasta eroava diagonaalin vas. puolella */
            {
              
              A_kj = array[j]; 
              A_jj = diag[j]; 

              /* pivot */

              z    = c_div(A_kj, A_jj); 


              /* p‰ivit‰ c */

              z2   = c_mul(z, c[j]); 
              c[k] = c_sub(c[k], z2);

              /* Lis‰‰ rivi j riville k s.e. A_kj menee nollaksi :  */

              /* A_kl = A_kl - A_kj / A_jj * A_jl */

              s = first[j]; 
              while(s > 0){
                l = colind[s]; 
                if (l > j){ /* t‰m‰ osa p‰ivittyy */
                  A_jl = val[s]; 
                  z2   = c_mul(A_jl, z); 
                  array[l] = c_sub(array[l], z2); /* A_kl */
                }
                
                s = next[s]; 
              }
            }
          m = next[m];           
        }

      /* siirr‰ p‰ivitettyneet arvot matriisiin */
      m=first[k]; 
      while(m > 0)
        {

          j = colind[m]; 

          val[m] = array[j]; 

          if (j == k)
            diag[j] = array[j]; 

          m = next[m]; 
        }
    }  /* k luuppi p‰‰ttyy t‰h‰n */


  /***************************************************************************/

  i = GaussianSubst(cmat, c, x, diag);

  /***************************************************************************/

  Mem(MEM_FREE, diag);
  Mem(MEM_FREE, array);

  return i;

  /***************************************************************************/
}

/*****************************************************************************/
