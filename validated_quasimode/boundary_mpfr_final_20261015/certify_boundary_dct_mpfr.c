/*
 * Full 128x128 Picard floor-J DCT with outward-rounded MPFR interval arithmetic.
 * Standalone except for the exact MPFR Bessel polynomial coefficient CSV and
 * exact rational 21-term Fourier seed list. Never uses a numerical DCT midpoint
 * as a proof. Exit code zero ONLY for a strict coefficient l1 bound.
 * Build: cc -O2 -std=c11 -o certify_boundary_dct_mpfr certify_boundary_dct_mpfr.c -l:libmpfr.so.6 -lgmp -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__has_include) && __has_include(<mpfr.h>)
#include <mpfr.h>
#else
 typedef long mpfr_prec_t;typedef long mpfr_exp_t;typedef int mpfr_sign_t;
 typedef struct{mpfr_prec_t _mpfr_prec;mpfr_sign_t _mpfr_sign;mpfr_exp_t _mpfr_exp;unsigned long *_mpfr_d;}__mpfr_struct;
 typedef __mpfr_struct mpfr_t[1];
#endif
extern void mpfr_init2(mpfr_t,long);extern void mpfr_clear(mpfr_t);
extern int mpfr_set(mpfr_t,const mpfr_t,int),mpfr_set_str(mpfr_t,const char*,int,int),mpfr_set_si(mpfr_t,long,int);
extern int mpfr_add(mpfr_t,const mpfr_t,const mpfr_t,int),mpfr_sub(mpfr_t,const mpfr_t,const mpfr_t,int);
extern int mpfr_mul(mpfr_t,const mpfr_t,const mpfr_t,int),mpfr_div_ui(mpfr_t,const mpfr_t,unsigned long,int);
extern int mpfr_mul_si(mpfr_t,const mpfr_t,long,int),mpfr_neg(mpfr_t,const mpfr_t,int);
extern int mpfr_sqrt(mpfr_t,const mpfr_t,int),mpfr_cos(mpfr_t,const mpfr_t,int),mpfr_sin(mpfr_t,const mpfr_t,int),mpfr_const_pi(mpfr_t,int);
extern int mpfr_cmp(const mpfr_t,const mpfr_t),mpfr_cmp_si(const mpfr_t,long);
extern int mpfr_printf(const char*,...);
#define D 3
#define U 2
#define PREC 256
#define KMAX 121
#define MAXN 128
#define MAXM 22
#define DIE(...) do{fprintf(stderr,__VA_ARGS__);fprintf(stderr,"\n");exit(2);}while(0)
typedef struct {mpfr_t l,h;} I;
static void init(I*x){mpfr_init2(x->l,PREC);mpfr_init2(x->h,PREC);}
static void clr(I*x){mpfr_clear(x->l);mpfr_clear(x->h);}
static void set(I*x,long a){mpfr_set_si(x->l,a,D);mpfr_set_si(x->h,a,U);}
static void cpy(I*x,const I*y){mpfr_set(x->l,y->l,D);mpfr_set(x->h,y->h,U);}
static void readi(I*x,const char*l,const char*h){if(mpfr_set_str(x->l,l,10,D)||mpfr_set_str(x->h,h,10,U))DIE("bad interval decimal");}
static void add(I*z,const I*x,const I*y){mpfr_add(z->l,x->l,y->l,D);mpfr_add(z->h,x->h,y->h,U);}
static void sub(I*z,const I*x,const I*y){mpfr_sub(z->l,x->l,y->h,D);mpfr_sub(z->h,x->h,y->l,U);}
static void neg(I*z,const I*x){mpfr_neg(z->l,x->h,D);mpfr_neg(z->h,x->l,U);}
static void mul(I*z,const I*x,const I*y){
 mpfr_t l[4],h[4];const mpfr_t *a[2]={&x->l,&x->h},*b[2]={&y->l,&y->h};
 for(int k=0;k<4;k++){mpfr_init2(l[k],PREC);mpfr_init2(h[k],PREC);int i=k/2,j=k%2;mpfr_mul(l[k],*a[i],*b[j],D);mpfr_mul(h[k],*a[i],*b[j],U);}
 int lo=0,hi=0;for(int k=1;k<4;k++){if(mpfr_cmp(l[k],l[lo])<0)lo=k;if(mpfr_cmp(h[k],h[hi])>0)hi=k;}
 mpfr_set(z->l,l[lo],D);mpfr_set(z->h,h[hi],U);
 for(int k=0;k<4;k++){mpfr_clear(l[k]);mpfr_clear(h[k]);}
}
static void mulsi(I*z,const I*x,long k){if(k>=0){mpfr_mul_si(z->l,x->l,k,D);mpfr_mul_si(z->h,x->h,k,U);}else{mpfr_mul_si(z->l,x->h,k,D);mpfr_mul_si(z->h,x->l,k,U);}}
static void divui(I*z,const I*x,unsigned long k){mpfr_div_ui(z->l,x->l,k,D);mpfr_div_ui(z->h,x->h,k,U);}
static void sqrti(I*z,const I*x){if(mpfr_cmp_si(x->l,0)<=0)DIE("sqrt interval invalid");mpfr_sqrt(z->l,x->l,D);mpfr_sqrt(z->h,x->h,U);}
static void trigi(I*z,const I*x,int sine){
 // sin/cos globally 1-Lipschitz: endpoint values expanded by input width
 mpfr_t p,q,t,u,w;mpfr_init2(p,PREC);mpfr_init2(q,PREC);mpfr_init2(t,PREC);mpfr_init2(u,PREC);mpfr_init2(w,PREC);
 if(sine){mpfr_sin(p,x->l,D);mpfr_sin(q,x->h,D);mpfr_sin(t,x->l,U);mpfr_sin(u,x->h,U);}
 else{mpfr_cos(p,x->l,D);mpfr_cos(q,x->h,D);mpfr_cos(t,x->l,U);mpfr_cos(u,x->h,U);}
 mpfr_sub(w,x->h,x->l,U);mpfr_sub(z->l,(mpfr_cmp(p,q)<0?p:q),w,D);mpfr_add(z->h,(mpfr_cmp(t,u)>0?t:u),w,U);
 mpfr_clear(p);mpfr_clear(q);mpfr_clear(t);mpfr_clear(u);mpfr_clear(w);
}
static void cosi(I*z,const I*x){trigi(z,x,0);}static void sini(I*z,const I*x){trigi(z,x,1);}
static I *newiv(size_t n){I*a=malloc(n*sizeof(I));if(!a)DIE("alloc");for(size_t i=0;i<n;i++)init(a+i);return a;}
static void freeiv(I*a,size_t n){for(size_t i=0;i<n;i++)clr(a+i);free(a);}
static void piiv(I*p){mpfr_const_pi(p->l,D);mpfr_const_pi(p->h,U);}
static void pr(const char*label,const I*x){printf("%s: [",label);mpfr_printf("%.30RDe",x->l);printf(", ");mpfr_printf("%.30RUe",x->h);printf("]\n");}

typedef struct{int a,b,n;long num,den;} Mode;
static Mode mode[MAXM];static int nm;
static I coeff[51][KMAX], deriv[51][KMAX];static int exists[51];
static void load_terms(const char*name){FILE*f=fopen(name,"r");if(!f)DIE("cannot read terms %s",name);
 while(nm<MAXM&&fscanf(f,"%d%d%d%ld%ld",&mode[nm].a,&mode[nm].b,&mode[nm].n,&mode[nm].num,&mode[nm].den)==5)nm++;
 fclose(f);if(nm!=21)DIE("expected exactly 21 modes, got %d",nm);
}
static void load_radial(const char*name){FILE*f=fopen(name,"r");if(!f)DIE("cannot open coefficients");
 char line[2048],s1[256],s2[256],s3[256],s4[256];int n,k,count=0;
 while(fgets(line,sizeof(line),f)){
  if(line[0]=='#')continue;
  if(sscanf(line,"%d,%d,%255[^,],%255[^,],%255[^,],%255s",&n,&k,s1,s2,s3,s4)!=6)DIE("invalid csv line");
  if(n<1||n>50||k<0||k>=KMAX)DIE("invalid radial index");
  if(!exists[n]){for(int v=0;v<KMAX;v++){init(&coeff[n][v]);init(&deriv[n][v]);}exists[n]=1;}
  readi(&coeff[n][k],s1,s2);readi(&deriv[n][k],s3,s4);count++;
 }
 fclose(f);if(count!=24*KMAX)DIE("incomplete radial table: %d",count);
 for(int i=0;i<nm;i++)if(!exists[mode[i].n])DIE("missing n=%d",mode[i].n);
}
static void rat(I*z,long a,long b){set(z,a);divui(z,z,b);}
static void cheb_table(const I*x,I*T){I a,b,c;init(&a);init(&b);init(&c);
 set(&T[0],1);cpy(&T[1],x);
 for(int k=2;k<KMAX;k++){
  mul(&a,x,&T[k-1]);mulsi(&b,&a,2);sub(&c,&b,&T[k-2]);cpy(&T[k],&c);
 }
 clr(&a);clr(&b);clr(&c);
}
static void eval_radial(const I*T,int n,I*F,I*G){I a,b,c,d;init(&a);init(&b);init(&c);init(&d);set(F,0);set(G,0);
 for(int k=0;k<KMAX;k++){
  mul(&a,&coeff[n][k],&T[k]);add(&b,F,&a);cpy(F,&b);
  mul(&c,&deriv[n][k],&T[k]);add(&d,G,&c);cpy(G,&d);
 }
 clr(&a);clr(&b);clr(&c);clr(&d);
}
static I nodes_x[MAXN],nodes_t[MAXN];
static I CX[8][MAXN],SX[8][MAXN],CT[8][MAXN],ST[8][MAXN];
static I cm[MAXN][MAXN];
static I samples[MAXN][MAXN];
static I partial[MAXN][MAXN];
static void initcoords(int N){I pi,tmp,ang,node,two_pi;init(&pi);init(&tmp);init(&ang);init(&node);init(&two_pi);piiv(&pi);mulsi(&two_pi,&pi,2);
 for(int i=0;i<N;i++){
  init(&nodes_x[i]);init(&nodes_t[i]);
  mulsi(&tmp,&pi,2*i+1);divui(&ang,&tmp,2*N);cosi(&node,&ang);
  divui(&nodes_x[i],&node,2);
  set(&tmp,1);add(&ang,&tmp,&node);divui(&nodes_t[i],&ang,4);
  for(int a=0;a<8;a++){
   init(&CX[a][i]);init(&SX[a][i]);init(&CT[a][i]);init(&ST[a][i]);
   mulsi(&tmp,&two_pi,a);mul(&ang,&tmp,&nodes_x[i]);cosi(&CX[a][i],&ang);sini(&SX[a][i],&ang);
   mul(&ang,&tmp,&nodes_t[i]);cosi(&CT[a][i],&ang);sini(&ST[a][i],&ang);
  }
 }
 for(int i=0;i<N;i++)for(int j=0;j<N;j++){
  init(&cm[i][j]);init(&samples[i][j]);init(&partial[i][j]);
  mulsi(&tmp,&pi,(2*i+1)*j);divui(&ang,&tmp,2*N);cosi(&cm[i][j],&ang);
 }
 clr(&pi);clr(&tmp);clr(&ang);clr(&node);clr(&two_pi);
}
static void Jnode(int i,int j,I*out){
 I x2,t2,s,z,y,tau,one,zero,F[51],G[51], T[KMAX],aa,bb,cc,dd,ee,ff,gg,hh,tmp,tmp2,rx,rt,angle,twopi,term,c;
 I *v[]={&x2,&t2,&s,&z,&y,&tau,&one,&zero,&aa,&bb,&cc,&dd,&ee,&ff,&gg,&hh,&tmp,&tmp2,&rx,&rt,&angle,&twopi,&term,&c};
 for(size_t k=0;k<sizeof(v)/sizeof(*v);k++)init(v[k]);for(int k=0;k<KMAX;k++)init(&T[k]);for(int n=0;n<=50;n++)if(exists[n]){init(&F[n]);init(&G[n]);}
 set(&one,1);set(&zero,0);piiv(&angle);mulsi(&twopi,&angle,2);
 mul(&x2,&nodes_x[i],&nodes_x[i]);mul(&t2,&nodes_t[j],&nodes_t[j]);add(&s,&x2,&t2);sub(&z,&one,&s);sqrti(&y,&z);
 rat(&aa,17,15);sub(&bb,&y,&aa);mulsi(&cc,&bb,15);divui(&tau,&cc,8);
 cheb_table(&tau,T);
 for(int n=1;n<=50;n++)if(exists[n])eval_radial(T,n,&F[n],&G[n]);
 set(out,0);
 for(int k=0;k<nm;k++){
  const Mode*m=mode+k;int a=m->a,b=m->b,n=m->n;
  // T = cos(ax)cos(bt)-cos(bx)cos(at)
  mul(&aa,&CX[a][i],&CT[b][j]);mul(&bb,&CX[b][i],&CT[a][j]);sub(&cc,&aa,&bb);
  // Tx = 2pi *(-a sin(ax) cos(bt) + b sin(bx) cos(at))
  mul(&aa,&SX[a][i],&CT[b][j]);mulsi(&bb,&aa,-a);
  mul(&aa,&SX[b][i],&CT[a][j]);mulsi(&dd,&aa,b);add(&ee,&bb,&dd);mul(&rx,&twopi,&ee);
  // Tt = 2pi *(-b cos(ax) sin(bt) + a cos(bx) sin(at))
  mul(&aa,&CX[a][i],&ST[b][j]);mulsi(&bb,&aa,-b);
  mul(&aa,&CX[b][i],&ST[a][j]);mulsi(&dd,&aa,a);add(&ee,&bb,&dd);mul(&rt,&twopi,&ee);
  mul(&aa,&nodes_x[i],&rx);mul(&bb,&nodes_t[j],&rt);add(&dd,&aa,&bb);
  mul(&aa,&F[n],&dd);mul(&bb,&y,&G[n]);mul(&ee,&bb,&cc);add(&gg,&aa,&ee);
  rat(&c,m->num,m->den);mul(&hh,&c,&gg);mul(&tmp,&y,&hh);sub(&tmp2,out,&tmp);cpy(out,&tmp2);
 }
 for(int n=0;n<=50;n++)if(exists[n]){clr(&F[n]);clr(&G[n]);}
 for(int k=0;k<KMAX;k++)clr(&T[k]);for(size_t k=0;k<sizeof(v)/sizeof(*v);k++)clr(v[k]);
}
int main(int argc,char**argv){
 if(argc<4)DIE("usage: %s radial_coeffs.csv terms.txt N [limit]",argv[0]);int N=atoi(argv[3]);if(N<2||N>128)DIE("N invalid");
 load_terms(argv[2]);load_radial(argv[1]);initcoords(N);
 I s,t,p,tmp,l1,maxwidth,negative,zero;I *vv[]={&s,&t,&p,&tmp,&l1,&maxwidth,&negative,&zero};for(int i=0;i<8;i++)init(vv[i]);set(&l1,0);set(&maxwidth,0);set(&zero,0);
 for(int i=0;i<N;i++){
  for(int j=0;j<N;j++){
   Jnode(i,j,&samples[i][j]);
   sub(&s,&samples[i][j],&zero);
   mpfr_sub(t.l,s.h,s.l,U);mpfr_set(t.h,t.l,U);
   if(mpfr_cmp(t.h,maxwidth.h)>0)cpy(&maxwidth,&t);
  }
  if(i%16==0){fprintf(stderr,"MPFR certified J rows: %d/%d\n",i,N);fflush(stderr);}
 }
 pr("max interval width among J node values",&maxwidth);
 // DCT of exact nodal intervals: 2/N for nonzero frequency, 1/N at DC.
 for(int k=0;k<N;k++)for(int j=0;j<N;j++){
  set(&s,0);
  for(int i=0;i<N;i++){
   mul(&t,&samples[i][j],&cm[i][k]);add(&p,&s,&t);cpy(&s,&p);
  }
  mulsi(&p,&s,k?2:1);divui(&partial[k][j],&p,N);
 }
 for(int k=0;k<N;k++)for(int l=0;l<N;l++){
  set(&s,0);
  for(int j=0;j<N;j++){
   mul(&t,&partial[k][j],&cm[j][l]);add(&p,&s,&t);cpy(&s,&p);
  }
  mulsi(&p,&s,l?2:1);divui(&tmp,&p,N);
  // abs interval upper enclosure
  if(mpfr_cmp_si(tmp.l,0)>=0){mpfr_add(l1.h,l1.h,tmp.h,U);mpfr_add(l1.l,l1.l,tmp.l,D);}
  else if(mpfr_cmp_si(tmp.h,0)<=0){mpfr_sub(l1.h,l1.h,tmp.l,U);mpfr_sub(l1.l,l1.l,tmp.h,D);}
  else{
   // only need rigorous upper (and trivial lower zero)
   I ab;init(&ab);mpfr_neg(ab.l,tmp.l,U);mpfr_set(ab.h,ab.l,U);
   mpfr_t *best=(mpfr_cmp(ab.h,tmp.h)>0?&ab.h:&tmp.h);
   mpfr_add(l1.h,l1.h,*best,U);clr(&ab);
  }
 }
 pr("outward DCT coefficient l1 (incl radial uncertainty)",&l1);
 I gate;init(&gate);readi(&gate,"3.5e-12","3.5e-12");int pass=mpfr_cmp(l1.h,gate.l)<0;
 printf("CERTIFICATE_DCT_L1_%s\n",pass?"PASS":"FAIL");
 clr(&gate);for(int i=0;i<8;i++)clr(vv[i]);return pass?0:3;
}
