/* MPFR-checkable explicit constants for the analytic quadrature theorem.
   Reuses the standalone directed interval primitives, not mpmath.iv.
   Build: gcc -O2 -std=c11 -o audit_constants_mpfr audit_constants_mpfr.c -l:libmpfr.so.6 -lgmp -lm
*/
#define main discarded_numerical_main
#include "verify_floor_mpfr.c"
#undef main
static void lt(const char*msg,const Iv*a,const Iv*b){
 if(mpfr_cmp(a->h,b->l)>=0)FAILURE("FAILED %s",msg);
 printf("PASS: %s\n",msg);
}
static void gt(const char*msg,const Iv*a,const Iv*b){
 if(mpfr_cmp(a->l,b->h)<=0)FAILURE("FAILED %s",msg);
 printf("PASS: %s\n",msg);
}
static void pw(Iv*out,Iv*base,int n){Iv acc,tmp;init(&acc);init(&tmp);val(&acc,1);for(int i=0;i<n;i++){prod(&tmp,&acc,base);cp(&acc,&tmp);}cp(out,&acc);clear(&acc);clear(&tmp);}
int main(void){
 Iv rho,inv,v,a,r,pi,t0,t1,t2,t3,bound,coshval,check,expo,t_small,kval,kprime,twopi;
 Iv*vs[]={&rho,&inv,&v,&a,&r,&pi,&t0,&t1,&t2,&t3,&bound,&coshval,&check,&expo,&t_small,&kval,&kprime,&twopi};
 for(int i=0;i<18;i++)init(vs[i]);
 rat(&rho,13,10);val(&t0,1);divi(&inv,&t0,&rho);rat(&r,6633,1000);mpfr_const_pi(pi.l,DN);mpfr_const_pi(pi.h,UP);mulsi(&twopi,&pi,2);
 /* Im u <= 3/2(13/10 -10/13); Re u <= 3+3/2(13/10+10/13) */
 sub(&t1,&rho,&inv);mulsi(&t2,&t1,3);divsi(&v,&t2,2);
 add(&t1,&rho,&inv);mulsi(&t2,&t1,3);divsi(&t3,&t2,2);val(&t0,3);add(&a,&t0,&t3);
 prod(&t0,&r,&v);cosh_pos(&coshval,&t0);val(&bound,102);lt("|cos(ru)| <= cosh(r Im u) < 102",&coshval,&bound);
 cosh_pos(&coshval,&a);val(&bound,226);lt("|cosh(u)| < 226",&coshval,&bound);
 cosi(&t1,&v);mulsi(&t2,&t1,4);rat(&bound,8,3);gt("4 cos(Im u) > 8/3",&t2,&bound);
 rat(&t0,-8,3);expi(&expo,&t0);rat(&bound,1,14);lt("exp(-8/3) < 1/14",&expo,&bound);
 rat(&t0,102*226,14);val(&bound,2000);lt("102 *226 /14 <2000",&t0,&bound);
 /* Floor ellipses: width Im x <=6/35, Im t <=3/35 */
 rat(&t0,6,35);prod(&t1,&twopi,&t0);cosh_pos(&coshval,&t1);rat(&bound,33,20);lt("cosh(12 pi /35)<33/20",&coshval,&bound);
 rat(&t0,3,35);prod(&t1,&twopi,&t0);cosh_pos(&coshval,&t1);rat(&bound,6,5);lt("cosh(6 pi /35)<6/5",&coshval,&bound);
 rat(&t0,2235,4900);rat(&bound,1,4);gt("Re(1-x²-t²)>1/4",&t0,&bound);
 mulsi(&t0,&pi,4);val(&bound,13);lt("4pi<13",&t0,&bound);
 /* K bound exp(-pi)/sqrt(2) < 1/20 */
 neg(&t0,&pi);expi(&t1,&t0);val(&t0,2);sqrti(&t2,&t0);divi(&kval,&t1,&t2);rat(&bound,1,20);lt("exp(-pi)/sqrt(2)<1/20 bounds |K|",&kval,&bound);
 /* K' bound exp(-pi) sqrt(pi/(2(pi-1)))<1/20 */
 val(&t0,1);sub(&t2,&pi,&t0);mulsi(&t3,&t2,2);divi(&t2,&pi,&t3);sqrti(&t3,&t2);prod(&kprime,&t1,&t3);rat(&bound,1,20);lt("exp(-pi)*sqrt(pi/(2(pi-1)))<1/20 bounds |K'|",&kprime,&bound);
 /* 2pi (3/5) (33/20+6/5) < 15 bounds |DU| */
 rat(&t0,33,20);rat(&t1,6,5);add(&t2,&t0,&t1);rat(&t0,3,5);prod(&t3,&t2,&t0);prod(&t2,&t3,&twopi);val(&bound,15);lt("|DU|<15",&t2,&bound);
 /* 4*[ (1/400)*(16+60) + 13*(1/400)*16] < 100 */
 rat(&t0,4*(76+208),400);val(&bound,100);lt("floor integrand |F|<100",&t0,&bound);
 rat(&rho,10,13);pw(&t0,&rho,159);mulsi(&t1,&t0,160000);str(&t0,"1e-100");add(&t2,&t1,&t0);str(&bound,"1e-12");lt("Bessel Gauss error n>=80<1e-12",&t2,&bound);
 rat(&rho,5,7);pw(&t0,&rho,111);mulsi(&t1,&t0,10000);str(&bound,"1e-12");lt("tensor spatial Gauss error n>=56<1e-12",&t1,&bound);
 puts("ALL MPFR ANALYTIC CONSTANT INEQUALITIES PASS");
 for(int i=0;i<18;i++)clear(vs[i]);return 0;
}