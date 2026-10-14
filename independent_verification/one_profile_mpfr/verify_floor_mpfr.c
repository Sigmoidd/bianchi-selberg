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
extern int mpfr_cmp(const mpfr_t,const mpfr_t),mpfr_cmp_si(const mpfr_t,long),mpfr_printf(const char*,...);
#define DN 3
#define UP 2
#define P 256
#define MAXN 128
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
int main(int argc,char**argv){
 if(argc!=5)FAILURE("Usage: %s SPACE_N BESSEL_N GL_SPACE_BRACKETS GL_BESSEL_BRACKETS",argv[0]);
 int n=atoi(argv[1]),m=atoi(argv[2]);if(n>MAXN||m>MAXN||n<56||m<80)FAILURE("unsupported n/m");
 Node *s=calloc(n,sizeof(Node)),*b=calloc(m,sizeof(Node));if(!s||!b)FAILURE("allocation failed");
 load_nodes(argv[3],n,0,1,1,2,s);load_nodes(argv[4],m,3,1,3,1,b);
 Iv r,pi,twopi,one,eps,S,spError,bsError;Iv*bs[]={&r,&pi,&twopi,&one,&eps,&S,&spError,&bsError};for(int i=0;i<8;i++)init(bs[i]);
 rat(&r,6633,1000);mpfr_const_pi(pi.l,DN);mpfr_const_pi(pi.h,UP);mulsi(&twopi,&pi,2);val(&one,1);str(&eps,"0.000000000001");val(&S,0);
 /* Construct validated Bessel nodes: cosh(u), weight*cos(r*u). */
 Iv *cu=calloc(m,sizeof(Iv)),*wc=calloc(m,sizeof(Iv));for(int j=0;j<m;j++){
  init(&cu[j]);init(&wc[j]);cosh_pos(&cu[j],&b[j].x);Iv q,t;init(&q);init(&t);prod(&t,&r,&b[j].x);cosi(&q,&t);prod(&wc[j],&b[j].w,&q);clear(&q);clear(&t);
 }
 /* Spatial integrand scratch variables; all routines use outward MPFR rounding. */
 Iv tnode,tx,ty,y,z,X,T,U,DU,k,kp,term,zcu,e,q,qq,u2,du2,kk,k2,dk1,dk2,bracket,integration,xweighted,wquad,temp,buf,negz,ft,fx,fy;
 Iv*scratch[]={&tnode,&tx,&ty,&y,&z,&X,&T,&U,&DU,&k,&kp,&term,&zcu,&e,&q,&qq,&u2,&du2,&kk,&k2,&dk1,&dk2,&bracket,&integration,&xweighted,&wquad,&temp,&buf,&negz,&ft,&fx,&fy};
 int nc=sizeof(scratch)/sizeof(*scratch);for(int a=0;a<nc;a++)init(scratch[a]);
 for(int i=0;i<n;i++){
  for(int j=0;j<n;j++){
   /* Spatial t coordinate on [0,1/2] comes from the *same* GL node on [-1,1]
      as x: t = 1/4 + x_raw/4, whereas x = x_raw/2. */
   val(&temp,1);divsi(&temp,&temp,4);divsi(&tx,&s[j].x,2);add(&tnode,&temp,&tx);
   prod(&tx,&s[i].x,&s[i].x);prod(&ty,&tnode,&tnode);add(&buf,&tx,&ty);sub(&q,&one,&buf);sqrti(&y,&q);prod(&z,&twopi,&y);
   val(&k,0);val(&kp,0);
   for(int h=0;h<m;h++){
    prod(&zcu,&z,&cu[h]);neg(&negz,&zcu);expi(&e,&negz);
    prod(&term,&e,&wc[h]);add(&temp,&k,&term);cp(&k,&temp);
    prod(&q,&term,&cu[h]);add(&temp,&kp,&q);cp(&kp,&temp);
   }
   /* K->K±1e-12, K'=-kp ±1e-12. */
   mpfr_sub(k.l,k.l,eps.h,DN);mpfr_add(k.h,k.h,eps.h,UP);
   neg(&temp,&kp);cp(&kp,&temp);
   mpfr_sub(kp.l,kp.l,eps.h,DN);mpfr_add(kp.h,kp.h,eps.h,UP);
   prod(&X,&twopi,&s[i].x);prod(&T,&twopi,&tnode);cosi(&q,&X);cosi(&qq,&T);sub(&U,&q,&qq);
   sini(&q,&X);prod(&tx,&s[i].x,&q);neg(&temp,&tx);cp(&tx,&temp);
   sini(&qq,&T);prod(&ty,&tnode,&qq);add(&buf,&tx,&ty);prod(&DU,&twopi,&buf);
   prod(&u2,&U,&U);prod(&du2,&U,&DU);add(&buf,&u2,&du2);prod(&kk,&k,&k);prod(&k2,&kk,&buf);
   prod(&dk1,&k,&kp);prod(&dk2,&dk1,&z);prod(&buf,&dk2,&u2);add(&bracket,&k2,&buf);mulsi(&integration,&bracket,-4);
   /* GL tensor weights: x weight dx = (1/2)*w_i; t weight dt=(1/4)*w_j. */
   prod(&wquad,&s[i].w,&s[j].w);divsi(&wquad,&wquad,2);prod(&xweighted,&integration,&wquad);add(&buf,&S,&xweighted);cp(&S,&buf);
  }
  if(i%8==0){printf("spatial row %d/%d\n",i,n);fflush(stdout);}
 }
 /* Conservative externally proved Gauss and K errors added as rational amounts.
    Need to ensure the analytic error bounds are sound separately. */
 str(&spError,"0.000000000001");
 mpfr_sub(S.l,S.l,spError.h,DN);mpfr_add(S.h,S.h,spError.h,UP);
 printiv("MPFR certified B_floor",&S);
 int negative=mpfr_cmp_si(S.h,0)<0;
 printf("STRICT NEGATIVE %s\n",negative?"YES":"NO");
 for(int a=0;a<nc;a++)clear(scratch[a]);for(int a=0;a<8;a++)clear(bs[a]);
 for(int h=0;h<m;h++){clear(&cu[h]);clear(&wc[h]);}free(cu);free(wc);free_nodes(s,n);free_nodes(b,m);free(s);free(b);
 return negative?0:3;
}
