/* Standalone independently implemented MPFR-directed Picard floor-flux verifier.
   Uses ONLY MPFR and rational Gaussian-node brackets in gl_n*.txt; validates
   the Legendre roots and obtains rigorous weights from them. It does NOT use
   mpmath, python-flint, SciPy, or their interval evaluator at runtime.
   Build: gcc -O2 -std=c11 -o verify_floor_mpfr verify_floor_mpfr.c -l:libmpfr.so.6 -lgmp -lm
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(__has_include) && __has_include(<mpfr.h>)
#include <mpfr.h>
#else
/* x86_64/LP64 libmpfr.so.6 ABI fallback for this container, which lacks
   development headers. Use the official header for portable builds. */
typedef long mpfr_prec_t; typedef long mpfr_exp_t; typedef int mpfr_sign_t;
typedef struct {mpfr_prec_t _mpfr_prec; mpfr_sign_t _mpfr_sign; mpfr_exp_t _mpfr_exp; unsigned long *_mpfr_d;} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];
#endif
extern void mpfr_init2(mpfr_t,mpfr_prec_t); extern void mpfr_clear(mpfr_t);
extern int mpfr_set(mpfr_t,const mpfr_t,int); extern int mpfr_set_str(mpfr_t,const char*,int,int);
extern int mpfr_set_si(mpfr_t,long,int); extern int mpfr_set_ui(mpfr_t,unsigned long,int);
extern int mpfr_add(mpfr_t,const mpfr_t,const mpfr_t,int),mpfr_sub(mpfr_t,const mpfr_t,const mpfr_t,int);
extern int mpfr_mul(mpfr_t,const mpfr_t,const mpfr_t,int),mpfr_div(mpfr_t,const mpfr_t,const mpfr_t,int);
extern int mpfr_div_si(mpfr_t,const mpfr_t,long,int),mpfr_ui_div(mpfr_t,unsigned long,const mpfr_t,int);
extern int mpfr_sqrt(mpfr_t,const mpfr_t,int),mpfr_exp(mpfr_t,const mpfr_t,int);
extern int mpfr_cos(mpfr_t,const mpfr_t,int),mpfr_sin(mpfr_t,const mpfr_t,int);
extern int mpfr_neg(mpfr_t,const mpfr_t,int),mpfr_const_pi(mpfr_t,int);
extern int __gmpfr_fprintf(FILE*,const char*,...);
extern int mpfr_cmp(const mpfr_t,const mpfr_t),mpfr_cmp_si(const mpfr_t,long),mpfr_printf(const char*,...);
#define DN 3
#define UP 2
#define P 1024
#define MAXN 192
#define FAILURE(...) do{fprintf(stderr,__VA_ARGS__);fprintf(stderr,"\n");exit(2);}while(0)
typedef struct {mpfr_t l,h;} Iv;
static void init(Iv *a){mpfr_init2(a->l,P);mpfr_init2(a->h,P);}
static void clear(Iv *a){mpfr_clear(a->l);mpfr_clear(a->h);}
static void val(Iv *a,long n){mpfr_set_si(a->l,n,DN);mpfr_set_si(a->h,n,UP);}
static void rat(Iv *a,long p,long q){mpfr_set_si(a->l,p,DN);mpfr_set_si(a->h,p,UP);mpfr_div_si(a->l,a->l,q,DN);mpfr_div_si(a->h,a->h,q,UP);}
static void str(Iv *a,const char*s){if(mpfr_set_str(a->l,s,10,DN)||mpfr_set_str(a->h,s,10,UP))FAILURE("invalid number %s",s);}
static void ends(Iv *a,const char*l,const char*h){if(mpfr_set_str(a->l,l,10,DN)||mpfr_set_str(a->h,h,10,UP))FAILURE("invalid bounds");}
static void cp(Iv *a,const Iv*b){mpfr_set(a->l,b->l,DN);mpfr_set(a->h,b->h,UP);}
static void add(Iv*c,const Iv*a,const Iv*b){mpfr_add(c->l,a->l,b->l,DN);mpfr_add(c->h,a->h,b->h,UP);}
static void sub(Iv*c,const Iv*a,const Iv*b){mpfr_sub(c->l,a->l,b->h,DN);mpfr_sub(c->h,a->h,b->l,UP);}
static void neg(Iv*c,const Iv*a){mpfr_neg(c->l,a->h,DN);mpfr_neg(c->h,a->l,UP);}
static void prod(Iv*c,const Iv*a,const Iv*b){
 mpfr_t low[4],high[4];const mpfr_t *aa[2]={&a->l,&a->h},*bb[2]={&b->l,&b->h};
 for(int i=0;i<4;i++){mpfr_init2(low[i],P);mpfr_init2(high[i],P);int j=i/2,k=i%2;mpfr_mul(low[i],*aa[j],*bb[k],DN);mpfr_mul(high[i],*aa[j],*bb[k],UP);}
 int il=0,ih=0;for(int i=1;i<4;i++){if(mpfr_cmp(low[i],low[il])<0)il=i;if(mpfr_cmp(high[i],high[ih])>0)ih=i;}
 mpfr_set(c->l,low[il],DN);mpfr_set(c->h,high[ih],UP);
 for(int i=0;i<4;i++){mpfr_clear(low[i]);mpfr_clear(high[i]);}
}
static void divi(Iv*c,const Iv*a,const Iv*b){
 if(mpfr_cmp_si(b->l,0)<=0&&mpfr_cmp_si(b->h,0)>=0)FAILURE("interval division by zero");
 Iv re;init(&re);mpfr_ui_div(re.l,1,b->h,DN);mpfr_ui_div(re.h,1,b->l,UP);prod(c,a,&re);clear(&re);
}
static void mulsi(Iv*c,const Iv*a,long i){Iv t;init(&t);val(&t,i);prod(c,a,&t);clear(&t);}
static void divsi(Iv*c,const Iv*a,long i){if(i<=0)FAILURE("divsi requires positive denominator");mpfr_div_si(c->l,a->l,i,DN);mpfr_div_si(c->h,a->h,i,UP);}
static void sqrti(Iv*c,const Iv*a){if(mpfr_cmp_si(a->l,0)<0)FAILURE("sqrt negative");mpfr_sqrt(c->l,a->l,DN);mpfr_sqrt(c->h,a->h,UP);}
static void expi(Iv*c,const Iv*a){mpfr_exp(c->l,a->l,DN);mpfr_exp(c->h,a->h,UP);}
/* Globally Lipschitz 1 on R: endpoint cos/sin ball enlarged by width. */
static void trigi(Iv*c,const Iv*a,int sine){
 mpfr_t l0,l1,h0,h1,w;mpfr_init2(l0,P);mpfr_init2(l1,P);mpfr_init2(h0,P);mpfr_init2(h1,P);mpfr_init2(w,P);
 if(sine){mpfr_sin(l0,a->l,DN);mpfr_sin(l1,a->h,DN);mpfr_sin(h0,a->l,UP);mpfr_sin(h1,a->h,UP);}
 else {mpfr_cos(l0,a->l,DN);mpfr_cos(l1,a->h,DN);mpfr_cos(h0,a->l,UP);mpfr_cos(h1,a->h,UP);}
 mpfr_sub(w,a->h,a->l,UP);
 mpfr_sub(c->l,(mpfr_cmp(l0,l1)<0?l0:l1),w,DN);
 mpfr_add(c->h,(mpfr_cmp(h0,h1)>0?h0:h1),w,UP);
 mpfr_clear(l0);mpfr_clear(l1);mpfr_clear(h0);mpfr_clear(h1);mpfr_clear(w);
}
static void cosi(Iv*c,const Iv*a){trigi(c,a,0);}static void sini(Iv*c,const Iv*a){trigi(c,a,1);}
static void printiv(const char*label,const Iv*a){printf("%s [",label);mpfr_printf("%.36RDe",a->l);printf(", ");mpfr_printf("%.36RUe",a->h);printf("]\n");}
static int signi(const Iv*a){return mpfr_cmp_si(a->l,0)>0?1:mpfr_cmp_si(a->h,0)<0?-1:0;}
/* Legendre recurrence and exact derivative enclosure */
static void leg(int n,const Iv*x,Iv*out,Iv*der){
 Iv p0,p1,p2,t0,t1,t2,xx,den;Iv *v[]={&p0,&p1,&p2,&t0,&t1,&t2,&xx,&den};for(int i=0;i<8;i++)init(v[i]);
 val(&p0,1);cp(&p1,x);
 if(n==0){cp(out,&p0);val(der,0);}
 else{
  for(int k=1;k<n;k++){
   prod(&t0,x,&p1);mulsi(&t0,&t0,2*k+1);mulsi(&t1,&p0,k);sub(&t2,&t0,&t1);divsi(&p2,&t2,k+1);
   cp(&p0,&p1);cp(&p1,&p2);
  }
  cp(out,&p1);
  prod(&t0,x,&p1);sub(&t1,&t0,&p0);mulsi(&t0,&t1,n);prod(&xx,x,x);val(&den,1);sub(&t2,&xx,&den);divi(der,&t0,&t2);
 }
 for(int i=0;i<8;i++)clear(v[i]);
}
typedef struct{Iv x,w;} Node;
static void load_nodes(const char*file,int n,long center_n,long center_d,long radius_n,long radius_d,Node*out){
 FILE*f=fopen(file,"r");if(!f)FAILURE("cannot open %s",file);int declared=0;if(fscanf(f,"%d",&declared)!=1||declared!=n)FAILURE("node file mismatch %s",file);
 Iv center,rad,one,t0,t1,t2,p0,p1,der,rootbox,weight,r2;Iv*tmp[]={&center,&rad,&one,&t0,&t1,&t2,&p0,&p1,&der,&rootbox,&weight,&r2};for(int k=0;k<12;k++)init(tmp[k]);rat(&center,center_n,center_d);rat(&rad,radius_n,radius_d);val(&one,1);
 Iv sum;init(&sum);val(&sum,0);
 mpfr_t last_right; mpfr_init2(last_right,P); mpfr_set_si(last_right,-2,DN);
 char lo[192],hi[192];
 for(int i=0;i<n;i++){
  if(fscanf(f,"%191s %191s",lo,hi)!=2)FAILURE("bad root line %d",i);
  ends(&rootbox,lo,hi);
  if(mpfr_cmp(rootbox.l,rootbox.h)>=0)FAILURE("reversed interval %d",i);
  if(i>0 && mpfr_cmp(last_right,rootbox.l)>=0)FAILURE("unordered or overlapping roots %d",i);
  mpfr_set(last_right,rootbox.h,UP);
  /* Sign reversal certified with directed-rounding interval polynomial evaluation. */
  Iv xl,xh;init(&xl);init(&xh);mpfr_set(xl.l,rootbox.l,DN);mpfr_set(xl.h,rootbox.l,UP);mpfr_set(xh.l,rootbox.h,DN);mpfr_set(xh.h,rootbox.h,UP);
  leg(n,&xl,&p0,&der);leg(n,&xh,&p1,&der);
  if(signi(&p0)==0||signi(&p1)==0||signi(&p0)==signi(&p1)){printiv("p(left)",&p0);printiv("p(right)",&p1);FAILURE("root %d not sign-bracketed for n=%d",i,n);}
  clear(&xl);clear(&xh);
  leg(n,&rootbox,&p0,&der);
  prod(&t0,&rootbox,&rootbox);sub(&t1,&one,&t0);prod(&t2,&der,&der);prod(&r2,&t1,&t2);
  val(&t0,2);divi(&weight,&t0,&r2);
  if(mpfr_cmp_si(weight.l,0)<=0)FAILURE("nonpositive weight");
  init(&out[i].x);init(&out[i].w);
  prod(&t0,&rad,&rootbox);add(&out[i].x,&center,&t0);prod(&out[i].w,&rad,&weight);
  add(&t0,&sum,&out[i].w);cp(&sum,&t0);
 }
 char extra[3];if(fscanf(f,"%2s",extra)==1)FAILURE("extra data in node file");fclose(f);
 mulsi(&t0,&rad,2);
 if(mpfr_cmp(sum.l,t0.h)>0||mpfr_cmp(sum.h,t0.l)<0)FAILURE("weights do not enclose 2*radius");
 printf("verified %d GL root brackets and positive weights: %s\n",n,file);
 for(int k=0;k<12;k++)clear(tmp[k]);clear(&sum);mpfr_clear(last_right);
}
static void free_nodes(Node*out,int n){for(int i=0;i<n;i++){clear(&out[i].x);clear(&out[i].w);}}
static void cosh_pos(Iv *out,const Iv *u){Iv ep,en,m,t;init(&ep);init(&en);init(&m);init(&t);expi(&ep,u);neg(&m,u);expi(&en,&m);add(&t,&ep,&en);divsi(out,&t,2);clear(&ep);clear(&en);clear(&m);clear(&t);}

/* Rigorous outward enclosures of ideal Chebyshev-Lobatto node values.
   Uniform Gauss 192 Bessel error <=3000*(10/13)^384+1e-100<1e-38,
   for z>=18/5, real r<7. See rigorous derivation in companion file.
   Node positions are themselves enclosed by interval cosine evaluation.
*/
int main(int argc,char **argv){
 if(argc!=5)FAILURE("usage: %s FIRST_NORM_INDEX LAST_NORM_INDEX NODES.csv COEFFS.csv",argv[0]);
 int fi=atoi(argv[1]),la=atoi(argv[2]);
 const int norms[]={1,2,4,5,8,9,10,13,16,17,18,20,25,26,29,32,34,36,37,40,41,45,49,50};
 const int N=sizeof(norms)/sizeof(norms[0]);
 if(fi<0||la> N||la<=fi)FAILURE("Bad norm slice; supply 0 <= FIRST < LAST <=24");
 FILE *out=fopen(argv[3],"w");if(!out)FAILURE("cannot open output");
 FILE *coeff=fopen(argv[4],"w");if(!coeff)FAILURE("cannot open coefficient output");
 fprintf(out,"# Picard real-r radial Chebyshev-Lobatto node intervals, degree=120, interval MPFR 512bits\n");
 fprintf(coeff,"# n,k,c_f_lo,c_f_hi,c_fprime_lo,c_fprime_hi (standard Chebyshev coefficient convention)\n");
 fprintf(out,"# n,j,ylo,yhi,flo,fhi,fylo,fyhi ; errors enlarged +-1e-38 for K and Kprime\n");
 Node *b=calloc(192,sizeof(Node));if(!b)FAILURE("allocation failed");
 load_nodes("gl_n192.txt",192,3,1,3,1,b);
 Iv r,pi,twopi,one,z,iy,a,y,t,ang,eps,K,Kp,expv,tmp,tmp2,term,negz,f,fy,rootn,nn,center,rad;
 Iv *v[]={&r,&pi,&twopi,&one,&z,&iy,&a,&y,&t,&ang,&eps,&K,&Kp,&expv,&tmp,&tmp2,&term,&negz,&f,&fy,&rootn,&nn,&center,&rad};
 for(size_t i=0;i<sizeof(v)/sizeof(*v);i++)init(v[i]);
 rat(&r,6622119348586011L,1000000000000000L); /* exact decimal 6.622119348586011 */
 mpfr_const_pi(pi.l,DN);mpfr_const_pi(pi.h,UP);mulsi(&twopi,&pi,2);
 rat(&center,17,15);rat(&rad,8,15);str(&eps,"1e-38");
 Iv cu[192],wc[192];
 for(int j=0;j<192;j++){
   init(&cu[j]);init(&wc[j]);cosh_pos(&cu[j],&b[j].x);
   prod(&tmp,&r,&b[j].x);cosi(&t,&tmp);prod(&wc[j],&t,&b[j].w);
 }
 Iv nodeF[121],nodeFy[121];
 for(int j=0;j<=120;j++){init(&nodeF[j]);init(&nodeFy[j]);}
 for(int index=fi;index<la;index++){
   val(&nn,norms[index]);sqrti(&rootn,&nn);prod(&a,&twopi,&rootn);
   for(int j=0;j<=120;j++){
     /* cos(pi*j/120) gives exact mathematical Lobatto node; interval-enclosed. */
     mulsi(&tmp,&pi,j);divsi(&ang,&tmp,120);cosi(&t,&ang);
     prod(&tmp,&rad,&t);add(&y,&center,&tmp);prod(&z,&a,&y);
     if(mpfr_cmp_si(y.l,0)<=0)FAILURE("Invalid radial node");
     val(&K,0);val(&Kp,0);
     for(int k=0;k<192;k++){
       prod(&tmp,&z,&cu[k]);neg(&negz,&tmp);expi(&expv,&negz);
       prod(&term,&expv,&wc[k]);add(&tmp,&K,&term);cp(&K,&tmp);
       prod(&tmp2,&term,&cu[k]);add(&tmp,&Kp,&tmp2);cp(&Kp,&tmp);
     }
     mpfr_sub(K.l,K.l,eps.h,DN);mpfr_add(K.h,K.h,eps.h,UP);
     neg(&tmp,&Kp);cp(&Kp,&tmp);
     mpfr_sub(Kp.l,Kp.l,eps.h,DN);mpfr_add(Kp.h,Kp.h,eps.h,UP);
     prod(&f,&y,&K);prod(&tmp,&z,&Kp);add(&fy,&K,&tmp);
     cp(&nodeF[j],&f);cp(&nodeFy[j],&fy);
     fprintf(out,"%d,%d,",norms[index],j);
     __gmpfr_fprintf(out,"%.48RDe,%.48RUe,%.48RDe,%.48RUe,%.48RDe,%.48RUe\n",y.l,y.h,f.l,f.h,fy.l,fy.h);
   }
   Iv sumF,sumDy,carg,cosv,contrib,prodval,dst;
   Iv *work[]={&sumF,&sumDy,&carg,&cosv,&contrib,&prodval,&dst};
   for(int l=0;l<7;l++)init(work[l]);
   for(int k=0;k<=120;k++){
     val(&sumF,0);val(&sumDy,0);
     for(int j=0;j<=120;j++){
       mulsi(&carg,&pi,k*j);divsi(&carg,&carg,120);cosi(&cosv,&carg);
       prod(&contrib,&nodeF[j],&cosv);
       if(j==0||j==120)divsi(&contrib,&contrib,2);
       add(&tmp,&sumF,&contrib);cp(&sumF,&tmp);
       prod(&contrib,&nodeFy[j],&cosv);
       if(j==0||j==120)divsi(&contrib,&contrib,2);
       add(&tmp,&sumDy,&contrib);cp(&sumDy,&tmp);
     }
     mulsi(&tmp,&sumF,2);divsi(&sumF,&tmp,120);
     mulsi(&tmp,&sumDy,2);divsi(&sumDy,&tmp,120);
     if(k==0||k==120){divsi(&tmp,&sumF,2);cp(&sumF,&tmp);divsi(&tmp,&sumDy,2);cp(&sumDy,&tmp);}
     fprintf(coeff,"%d,%d,",norms[index],k);
     __gmpfr_fprintf(coeff,"%.48RDe,%.48RUe,%.48RDe,%.48RUe\n",sumF.l,sumF.h,sumDy.l,sumDy.h);
   }
   for(int l=0;l<7;l++)clear(work[l]);
   printf("Certified radial n=%d, 121 Chebyshev nodal intervals + both coefficient tables\n",norms[index]);fflush(stdout);
 }
 for(int j=0;j<=120;j++){clear(&nodeF[j]);clear(&nodeFy[j]);}
 for(int j=0;j<192;j++){clear(&cu[j]);clear(&wc[j]);}free_nodes(b,192);free(b);
 for(size_t i=0;i<sizeof(v)/sizeof(*v);i++)clear(v[i]);
 fclose(out);fclose(coeff);return 0;
}
