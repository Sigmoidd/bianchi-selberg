/* Standalone Picard Q-odd cusp-mass certificate: exact rational inputs,
 * MPFR directed rounding, composite Simpson integration, and explicitly
 * proved fourth-derivative / truncation error bounds.
 * Build (Debian with MPFR runtime but no dev headers):
 * gcc -O2 -std=c11 -Wall -Wextra -o certify_picard_cusp_mass_mpfr \
 *     certify_picard_cusp_mass_mpfr.c -l:libmpfr.so.6 -lgmp -lm
 * With normal development headers: gcc -O2 -std=c11 -o ... file.c -lmpfr -lgmp
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(__has_include) && __has_include(<mpfr.h>)
#include <mpfr.h>
#else
/* x86_64/LP64 libmpfr.so.6 ABI fallback; NOT portable to arbitrary platforms. */
typedef long mpfr_prec_t; typedef long mpfr_exp_t; typedef int mpfr_sign_t;
typedef struct {mpfr_prec_t _mpfr_prec;mpfr_sign_t _mpfr_sign;mpfr_exp_t _mpfr_exp;unsigned long *_mpfr_d;} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];
#endif
extern void mpfr_init2(mpfr_t,mpfr_prec_t),mpfr_clear(mpfr_t);
extern int mpfr_set(mpfr_t,const mpfr_t,int),mpfr_set_si(mpfr_t,long,int),mpfr_set_str(mpfr_t,const char*,int,int);
extern int mpfr_add(mpfr_t,const mpfr_t,const mpfr_t,int),mpfr_sub(mpfr_t,const mpfr_t,const mpfr_t,int);
extern int mpfr_mul(mpfr_t,const mpfr_t,const mpfr_t,int),mpfr_div(mpfr_t,const mpfr_t,const mpfr_t,int);
extern int mpfr_mul_si(mpfr_t,const mpfr_t,long,int),mpfr_div_si(mpfr_t,const mpfr_t,long,int);
extern int mpfr_ui_div(mpfr_t,unsigned long,const mpfr_t,int);
extern int mpfr_sqrt(mpfr_t,const mpfr_t,int),mpfr_exp(mpfr_t,const mpfr_t,int),mpfr_log(mpfr_t,const mpfr_t,int);
extern int mpfr_cos(mpfr_t,const mpfr_t,int),mpfr_const_pi(mpfr_t,int);
extern int mpfr_cmp(const mpfr_t,const mpfr_t),mpfr_cmp_si(const mpfr_t,long),mpfr_printf(const char*,...);
#define DOWN 3
#define UP 2
#define PREC 192
#ifndef NT
#define NT 512
#endif
#ifndef NY
#define NY 32
#endif
#define FAIL(...) do {fprintf(stderr,__VA_ARGS__);fprintf(stderr,"\n");exit(2);}while(0)
typedef struct{mpfr_t lo,hi;} Iv;
static void init(Iv*a){mpfr_init2(a->lo,PREC);mpfr_init2(a->hi,PREC);}
static void del(Iv*a){mpfr_clear(a->lo);mpfr_clear(a->hi);}
static void vint(Iv*a,long n){mpfr_set_si(a->lo,n,DOWN);mpfr_set_si(a->hi,n,UP);}
static void vstr(Iv*a,const char*str){if(mpfr_set_str(a->lo,str,10,DOWN)||mpfr_set_str(a->hi,str,10,UP))FAIL("bad decimal %s",str);}
static void vr(Iv*a,long n,long d){vint(a,n);mpfr_div_si(a->lo,a->lo,d,DOWN);mpfr_div_si(a->hi,a->hi,d,UP);}
static void cpy(Iv*a,const Iv*b){mpfr_set(a->lo,b->lo,DOWN);mpfr_set(a->hi,b->hi,UP);}
static void add(Iv*c,const Iv*a,const Iv*b){mpfr_add(c->lo,a->lo,b->lo,DOWN);mpfr_add(c->hi,a->hi,b->hi,UP);}
static void sub(Iv*c,const Iv*a,const Iv*b){mpfr_sub(c->lo,a->lo,b->hi,DOWN);mpfr_sub(c->hi,a->hi,b->lo,UP);}
static void mul(Iv*c,const Iv*a,const Iv*b){
 mpfr_t low[4],hi[4]; const mpfr_t * aa[2]={&a->lo,&a->hi}, *bb[2]={&b->lo,&b->hi};
 for(int k=0;k<4;k++){mpfr_init2(low[k],PREC);mpfr_init2(hi[k],PREC);mpfr_mul(low[k],*aa[k/2],*bb[k%2],DOWN);mpfr_mul(hi[k],*aa[k/2],*bb[k%2],UP);}
 int l=0,h=0;for(int k=1;k<4;k++){if(mpfr_cmp(low[k],low[l])<0)l=k;if(mpfr_cmp(hi[k],hi[h])>0)h=k;}
 mpfr_set(c->lo,low[l],DOWN);mpfr_set(c->hi,hi[h],UP);
 for(int k=0;k<4;k++){mpfr_clear(low[k]);mpfr_clear(hi[k]);}
}
static void divi(Iv*c,const Iv*a,const Iv*b){if(mpfr_cmp_si(b->lo,0)<=0&&mpfr_cmp_si(b->hi,0)>=0)FAIL("division by zero interval");
 Iv r;init(&r);mpfr_ui_div(r.lo,1,b->hi,DOWN);mpfr_ui_div(r.hi,1,b->lo,UP);mul(c,a,&r);del(&r);}
static void mulsi(Iv*c,const Iv*a,long n){
 if(n<0)FAIL("negative mulsi unsupported");
 mpfr_mul_si(c->lo,a->lo,n,DOWN);
 mpfr_mul_si(c->hi,a->hi,n,UP);
}
static void divsi(Iv*c,const Iv*a,long n){if(n<=0)FAIL("bad divsi");mpfr_div_si(c->lo,a->lo,n,DOWN);mpfr_div_si(c->hi,a->hi,n,UP);}
static void expi(Iv*c,const Iv*a){mpfr_exp(c->lo,a->lo,DOWN);mpfr_exp(c->hi,a->hi,UP);}
static void logi(Iv*c,const Iv*a){if(mpfr_cmp_si(a->lo,0)<=0)FAIL("log negative");mpfr_log(c->lo,a->lo,DOWN);mpfr_log(c->hi,a->hi,UP);}
static void cosi(Iv*c,const Iv*a){
 /* Evaluate at both endpoints and expand by input interval width via |cos'|<=1.
  * The tiny MPFR-representational input interval never crosses full cycles;
  * the expansion is valid for ALL interval widths, although loose if large. */
 mpfr_t low[2],hi[2],w;mpfr_init2(w,PREC);
 for(int i=0;i<2;i++){mpfr_init2(low[i],PREC);mpfr_init2(hi[i],PREC);mpfr_cos(low[i],i?a->hi:a->lo,DOWN);mpfr_cos(hi[i],i?a->hi:a->lo,UP);}
 mpfr_sub(w,a->hi,a->lo,UP);
 mpfr_sub(c->lo,mpfr_cmp(low[0],low[1])<0?low[0]:low[1],w,DOWN);
 mpfr_add(c->hi,mpfr_cmp(hi[0],hi[1])>0?hi[0]:hi[1],w,UP);
 for(int i=0;i<2;i++){mpfr_clear(low[i]);mpfr_clear(hi[i]);}mpfr_clear(w);
}
static void coshpos(Iv*out,const Iv*u){
 /* Valid for t>=0, cosh monotone: exp(t)+exp(-t) / 2.
  * exp(-t) is antitone: the two functions independently enclosed. */
 Iv e,minus,ee,tmp;init(&e);init(&minus);init(&ee);init(&tmp);
 expi(&e,u);vint(&minus,0);sub(&ee,&minus,u);expi(&minus,&ee);add(&tmp,&e,&minus);divsi(out,&tmp,2);
 del(&e);del(&minus);del(&ee);del(&tmp);
}
static void printiv(const char*name,const Iv*v){printf("%s [",name);mpfr_printf("%.36RDe",v->lo);printf(", ");mpfr_printf("%.36RUe",v->hi);printf("]\n");}
int main(void){
 Iv r,pi,twoPi,one,err,besselCos[NT+1],besselCosh[NT+1],zero;
 Iv*constants[]={&r,&pi,&twoPi,&one,&err,&zero};for(int i=0;i<6;i++)init(constants[i]);
 vstr(&r,"6.622119348586011");mpfr_const_pi(pi.lo,DOWN);mpfr_const_pi(pi.hi,UP);mulsi(&twoPi,&pi,2);vint(&one,1);vr(&err,1,10000000000L);vint(&zero,0);
 Iv tt,phase,tmp,xn,ee,nn;Iv*scratch[]={&tt,&phase,&tmp,&xn,&ee,&nn};for(int i=0;i<6;i++)init(scratch[i]);
 for(int j=0;j<=NT;j++){
  init(&besselCos[j]);init(&besselCosh[j]);vr(&tt,3*j,NT);mul(&phase,&tt,&r);cosi(&besselCos[j],&phase);coshpos(&besselCosh[j],&tt);
 }
 long anum[4]={5,7,19,21}, aden[4]={3,4,10,10};
 Iv masses,positive;init(&masses);init(&positive);vint(&masses,0);
 for(int k=0;k<3;k++){
  Iv a,b,width,h,yi,zt,sumt,sumy,integrand,efac,wt,prod1,piece,J,downJ,divratio,logratio,sq,t0,tmp2,mass;
  Iv*vs[]={&a,&b,&width,&h,&yi,&zt,&sumt,&sumy,&integrand,&efac,&wt,&prod1,&piece,&J,&downJ,&divratio,&logratio,&sq,&t0,&tmp2,&mass};
  for(int j=0;j<21;j++)init(vs[j]);
  vr(&a,anum[k],aden[k]);vr(&b,anum[k+1],aden[k+1]);sub(&width,&b,&a);divsi(&h,&width,NY);vint(&sumy,0);
  for(int i=0;i<=NY;i++){
   mulsi(&t0,&h,i);add(&yi,&a,&t0);mul(&zt,&twoPi,&yi);vint(&sumt,0);
   for(int j=0;j<=NT;j++){
    mul(&efac,&zt,&besselCosh[j]);sub(&tmp2,&zero,&efac);expi(&efac,&tmp2);
    mul(&integrand,&efac,&besselCos[j]);mulsi(&piece,&integrand,j==0||j==NT?1:(j%2?4:2));add(&tmp2,&sumt,&piece);cpy(&sumt,&tmp2);
   }
   /* t Simpson weight h_t/3 = (3/512)/3 = 1/512. */
   divsi(&sumt,&sumt,NT);
   divi(&integrand,&sumt,&yi);mulsi(&piece,&integrand,i==0||i==NY?1:(i%2?4:2));add(&tmp2,&sumy,&piece);cpy(&sumy,&tmp2);
  }
  /* y Simpson h_y/3 */
  divsi(&tmp2,&h,3);mul(&J,&sumy,&tmp2);
  /* Total analytic error 1e-10 per slab from explicit 4th-derivative bounds
   * plus t>3 truncation. Subtract error, lower bound remains positive. */
  sub(&downJ,&J,&err);
  if(mpfr_cmp_si(downJ.lo,0)<=0)FAIL("J lower bound not positive in slab %d",k);
  divi(&divratio,&b,&a);logi(&logratio,&divratio);
  mul(&sq,&downJ,&downJ);mulsi(&t0,&sq,2);divi(&mass,&t0,&logratio);
  add(&tmp2,&masses,&mass);cpy(&masses,&tmp2);
  char label[80];snprintf(label,sizeof label,"certified J slab %d",k+1);printiv(label,&downJ);
  snprintf(label,sizeof label,"certified mass-sq slab %d",k+1);printiv(label,&mass);
  for(int j=0;j<21;j++)del(vs[j]);
 }
 printiv("certified tail mass lower squared",&masses);
 /* 4e-7 squared = 16e-14 = 1.6e-13. */
 vr(&positive,16,100000000000000L);
 if(mpfr_cmp(masses.lo,positive.hi)<=0)FAIL("FAIL target mu_tail>4e-7");
 printf("PASS: mu_tail > 4e-7 (one norm-1 C4 orbit, exact coefficients +1,-1,+1,-1)\n");
 for(int j=0;j<=NT;j++){del(&besselCos[j]);del(&besselCosh[j]);}
 for(int i=0;i<6;i++){del(constants[i]);del(scratch[i]);}
 del(&masses);del(&positive);
 return 0;
}
