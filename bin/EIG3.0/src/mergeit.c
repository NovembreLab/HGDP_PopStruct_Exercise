#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>

#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  

#define WVERSION   "2001" 
#define MAXFL  50   
#define MAXSTR  512

// phasedmode added
// bugfix for pedcols
// bugfix for maxlinelength

char *trashdir = "/var/tmp" ;
int verbose = NO ;
int qtmode = NO ;
Indiv **indm1, **indm2 ;  
SNP **snpm1, **snpm2 ;
int nums1, nums2 ; 
int numi1, numi2 ; 

char  *snp1 = NULL ;
char  *ind1 = NULL ;
char  *geno1 = NULL ;

char  *snp2 = NULL ;
char  *ind2 = NULL ;
char  *geno2 = NULL ;

char *indoutfilename = NULL ;
char *snpoutfilename = NULL ;
char  *genooutfilename = NULL ;
char *badsnpname = NULL ;

int packout = -1 ;
int tersem  = YES ;
extern enum outputmodetype outputmode  ;
extern int checksizemode ;
char *omode = "packedancestrymap" ;
extern int packmode ;
int ogmode = NO ;
int docheck = YES ;
int strandcheck = YES ;
int phasedmode = NO ;

int xchrom = -1 ;
int lopos = -999999999 ; 
int hipos = 999999999 ;
int minchrom = 1 ; 
int maxchrom = 97 ;

/** 
 docheck  YES  allele flipping check
 hashcheck  YES  ... NO => snp names, Indiv names changed MUST retain order and number
 strandcheck (default YES) if NO then alleles are assumed on same strand
*/

char  unknowngender = 'U' ;

void setomode(enum outputmodetype *outmode, char *omode)  ;
void readcommands(int argc, char **argv) ;
void outfiles(char *snpname, char *indname, char *gname, SNP **snpm, 
  Indiv **indiv, int numsnps, int numind, int packem, int ogmode) ;
char compbase(char x) ;
int checkmatch(SNP *cupt1, SNP *cupt2) ;

int 
mergeit(SNP **snpm1, SNP **snpm2, Indiv ***pindm1, Indiv **indm2, 
 int nums1, int nums2, int numi1, int numi2)   ;


int main(int argc, char **argv)
{
  SNP **snpmarkers ;
  Indiv **indivmarkers ;
  int numsnps, numindivs ;
  unsigned char *packg1, *packg2 ;

  int **snppos ;
  int *snpindx ;
  int  lsnplist, lindlist, numeg ;
  int i,j; 
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;

  int ch1, ch2 ;
  int fmnum , lmnum ;
  int num, n1, n2 ;
  int nkill = 0 ;
  int t, k, x ;

  int nignore, numrisks = 1 ;

  char **genolist ;
  int numgenolist ;
  int maxmiss ; 

  tersem = YES ;     // no snp counts

  readcommands(argc, argv) ;

  setomode(&outputmode, omode) ;
  packmode = YES ;
  settersemode(tersem) ;

  nums1 = 
    getsnps(snp1, &snpm1, 0.0, NULL, &nignore, numrisks) ;

  putped(1) ;
  freeped() ;

  nums2 = 
    getsnps(snp2, &snpm2, 0.0, NULL, &nignore, numrisks) ;

  putped(2) ;
  freeped() ;

  for (x=0; x<nums1; ++x)  {  
   cupt1 = snpm1[x] ;
   cupt1 -> tagnumber = -1 ;
  }
  for (x=0; x<nums2; ++x)  {  
   cupt2 = snpm2[x] ;
   t = x %1000 ;   
// if (t==0) printf("zz %d %d\n", x, nums2) ;

   k = snpindex(snpm1, nums1, cupt2 -> ID) ;  
   if (k<0) { 
    cupt2 -> ignore = YES ;
    continue ;
   }
   cupt1 = snpm1[k] ;
   cupt1 -> tagnumber = x ;
   t = checkmatch(cupt1, cupt2) ;
   if (t==1) continue ;
   if (t==2) {  
    cupt2 -> isrfake = YES ;
    continue ;
   }
   if (t<0) {  
    cupt1  -> ignore = cupt2 -> ignore = YES ;
    continue ;
   }
   printf("allele funny: %s", cupt1 -> ID) ;
   printalleles(cupt1, stdout) ;
   printalleles(cupt2, stdout) ;
   printnl() ;
   cupt1  -> ignore = cupt2 -> ignore = YES ;
   continue ;
  }
  freesnpindex() ;
  numi1 = getindivs(ind1, &indm1) ;
  numi2 = getindivs(ind2, &indm2) ;

  for (x=0; x<numi2; ++x) {  
   k = indindex(indm1, numi1, indm2[x] -> ID) ;
// this code could be modified to allow duplicate individuals
   if (k>=0) fatalx("dup ind: %s\n", indm2[x] -> ID) ;  // fix later?  
  }

  setgenotypename(&geno1, ind1) ;
  getped(1) ;
  getgenos(geno1, snpm1, indm1, 
     nums1, numi1, nignore) ;

  packg1 = (unsigned char *) getpackgenos() ;
  clearpackgenos() ;

  setgenotypename(&geno2, ind2) ;
  getped(2) ;
  getgenos(geno2, snpm2, indm2, 
     nums2, numi2, nignore) ;

  packg2 = (unsigned char *) getpackgenos() ;
  numindivs = mergeit(snpm1, snpm2, &indm1, indm2, nums1, nums2, numi1, numi2) ;

  snpmarkers = snpm1 ; 
  numsnps = nums1 ;
  indivmarkers = indm1 ; 

  free(packg1) ;
  free(packg2) ;

  outfiles(snpoutfilename, indoutfilename, genooutfilename, 
   snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode) ;

  printf("##end of mergeit run\n") ;
  return 0 ;
}
int checkmatch(SNP *cupt1, SNP *cupt2) 
{

  char a1, a2, b1 , b2 ;

  if (docheck == NO) return 1 ;
  if (cupt1 -> alleles == NULL) return -1 ;
  if (cupt2 -> alleles == NULL) return -1 ;

  a1 = cupt1 -> alleles[0] ;
  a2 = cupt1 -> alleles[1] ;

  b1 = cupt2 -> alleles[0] ;
  b2 = cupt2 -> alleles[1] ;

  if ((a1 == 'X') && (a2 == 'X') && (strandcheck)) return -1 ;  // flipcheck impossible

  if (strandcheck) {

   if ((a1 == 'A') && (a2 == 'T')) return -1 ;
   if ((a1 == 'T') && (a2 == 'A')) return -1 ;
   if ((a1 == 'C') && (a2 == 'G')) return -1 ;
   if ((a1 == 'G') && (a2 == 'C')) return -1 ;

   if ((b1 == 'A') && (b2 == 'T')) return -1 ;
   if ((b1 == 'T') && (b2 == 'A')) return -1 ;
   if ((b1 == 'C') && (b2 == 'G')) return -1 ;
   if ((b1 == 'G') && (b2 == 'C')) return -1 ;
  }

  if ((a1 == b1) && (a2 == 'X')) {           
   cupt1 -> alleles[1] = b2 ;
   return 1 ;
  }

  if ((a1 == b2) && (a2 == 'X')) {           
   cupt1 -> alleles[1] = b1 ;
   return 2 ;
  }


  if ((a1 == b1) && (a2 == b2)) return 1 ;
  if ((a1 == b2) && (a2 == b1)) return 2 ;

  if ((a1 == b1) && (b2 =='X')) return 1 ;
  if ((a2 == b1) && (b2 =='X')) return 2 ;

  if (strandcheck == NO) return 0 ;

  b1 = compbase(b1) ;
  b2 = compbase(b2) ;

  if ((a1 == b1) && (a2 == 'X')) {           
   cupt1 -> alleles[1] = b2 ;
   return 1 ;
  }

  if ((a1 == b2) && (a2 == 'X')) {           
   cupt1 -> alleles[1] = b1 ;
   return 2 ;
  }

  if ((a1 == b1) && (a2 == b2)) return 1 ;
  if ((a1 == b2) && (a2 == b1)) return 2 ;

  if ((a1 == b1) && (b2 =='X')) return 1 ;
  if ((a2 == b1) && (b2 =='X')) return 2 ;

  return 0 ;

}

char compbase(char x) 
{
 if (x=='A') return 'T' ;
 if (x=='C') return 'G' ;
 if (x=='G') return 'C' ;
 if (x=='T') return 'A' ;

 return x ;

}

void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:vVf")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case 'f':
	phasedmode = YES ;                      
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring(ph, "geno1:", &geno1) ;
   getstring(ph, "snp1:", &snp1) ;
   getstring(ph, "ind1:", &ind1) ;

   getstring(ph, "geno2:", &geno2) ;
   getstring(ph, "snp2:", &snp2) ;
   getstring(ph, "ind2:", &ind2) ;

   getstring(ph, "indoutfilename:", &indoutfilename) ;
   getstring(ph, "indivoutname:", &indoutfilename) ;

   getstring(ph, "snpoutfilename:", &snpoutfilename) ;
   getstring(ph, "snpoutname:", &snpoutfilename) ;

   getstring(ph, "genooutfilename:", &genooutfilename) ; 
   getstring(ph, "genotypeoutname:", &genooutfilename) ; 

   getstring(ph, "outputformat:", &omode) ;

   getint(ph, "docheck:", &docheck) ;
   getint(ph, "hashcheck:", &hashcheck) ;
   getint(ph, "strandcheck:", &strandcheck) ;
   getint(ph, "phasedmode:", &phasedmode) ;

   writepars(ph) ;
   closepars(ph) ;

}
int 
mergeit(SNP **snpm1, SNP **snpm2, Indiv ***pindm1, Indiv **indm2, 
 int nums1, int nums2, int numi1, int numi2)   
{
   SNP *cupt1, *cupt2 ;
   int k, x, g, t  ;
   double y ;
   long rlen, packlen ; 
   static unsigned char *packg ;
   unsigned char *buff ;
   Indiv **indm1 ;
   static Indiv **indivmarkers ;
   int numindivs, numsnps ;

   indm1 = *pindm1 ;
   numindivs = numi1 + numi2 ;
   numsnps = nums1 ;
   ZALLOC(indivmarkers, numindivs, Indiv *) ;

   t = 0 ;
   for (x=0; x<numi1; ++x)  {  
    indivmarkers[t] = indm1[x]  ;
    ++t ;
   }
   for (x=0; x<numi2; ++x)  {  
    indivmarkers[t] = indm2[x]  ;
    ++t ;
   }
// we don't bother with a destructor here.   Sloppy code

   y = (double) (numindivs * 2) / (8 * (double) sizeof (char)) ;
   rlen = nnint(ceil(y)) ;
   rlen = MAX(rlen, 48)  ;
   packlen = numsnps*rlen ;
   ZALLOC(packg, packlen, unsigned char) ;
   cclear((unsigned char *) packg, 0XFF, packlen) ;
// wipe to invalid

   buff = packg ;
   for (k=0; k<nums1; k++) {  
    cupt1 = snpm1[k] ;
    x = cupt1 -> tagnumber ;
    if (x < 0 ) cupt1 -> ignore = YES ;
    if (cupt1 -> ignore) continue ;
    cupt2 = snpm2[x] ; 
    if (cupt2 -> isrfake) { 
     if (phasedmode == NO)  flipalleles(cupt2) ;
     if (phasedmode == YES)  flipalleles_phased(cupt2) ;
    }
    for (t=0; t<numi1; ++t) {  
      g = getgtypes(cupt1, t) ;
      if (g<0) continue ;
      wbuff((unsigned char *)buff, t, g) ;
    }
    for (t=0; t<numi2; ++t) {  
      g = getgtypes(cupt2, t) ;
      if (g<0) continue ;
      wbuff((unsigned char *)buff, numi1+t, g) ;
    }
    cupt1 -> ngtypes = numindivs ;
    cupt1 -> pbuff = (char *) buff ;
    buff += rlen ;
   }
   *pindm1 = indivmarkers ;
   return numindivs ;
}
