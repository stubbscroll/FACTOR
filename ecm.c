#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

gmp_randstate_t gmprand;

/* ECM from "factorization and primality testing" (david m bressoud) */
/* sample factorizations in main take 1634 seconds */

void sub_2i(mpz_t r,mpz_t s,mpz_t n,mpz_t a,mpz_t b,mpz_t X2,mpz_t Z2) {
	mpz_t t,u,v;
	mpz_init(t); mpz_init(u); mpz_init(v);
	/* calculate X2 */
	mpz_mul(t,r,r); mpz_mul(u,s,s); mpz_mul(u,u,a);
	mpz_sub(t,t,u); mpz_mod(t,t,n);
	mpz_mul(t,t,t);
	mpz_powm_ui(u,s,3,n); mpz_mul(u,u,r); mpz_mul(u,u,b); mpz_mul_ui(u,u,8);
	mpz_sub(t,t,u); mpz_mod(v,t,n);
	/* calculate Z2 */
	mpz_powm_ui(t,r,3,n); mpz_mul(u,a,r); mpz_mul(u,u,s); mpz_mul(u,u,s);
	mpz_add(t,t,u);
	mpz_powm_ui(u,s,3,n); mpz_mul(u,u,b); mpz_add(t,t,u); mpz_mod(t,t,n);
	mpz_mul(t,s,t); mpz_mul_ui(t,t,4); mpz_mod(Z2,t,n); mpz_set(X2,v);
	mpz_clear(t); mpz_clear(u); mpz_clear(v);
}

void sub_2i_plus(mpz_t r,mpz_t s,mpz_t u,mpz_t v,mpz_t X,mpz_t Z,mpz_t n,mpz_t a,mpz_t b,mpz_t U1,mpz_t U2) {
	mpz_t t,t1,t2,u1;
	mpz_init(t); mpz_init(t1); mpz_init(t2); mpz_init(u1);
	/* obtain U1 */
	mpz_mul(t1,r,u); mpz_mul(t,a,s); mpz_mul(t,t,v); mpz_sub(t1,t1,t); mpz_mod(t1,t1,n);
	mpz_mul(t2,r,v); mpz_mul(t,s,u); mpz_add(t2,t2,t); mpz_mul(t2,t2,v);
	mpz_mul(t2,t2,s); mpz_mul(t2,t2,b); mpz_mod(t2,t2,n);
	mpz_mul(t1,t1,t1); mpz_mul_ui(t2,t2,4); mpz_sub(t1,t1,t2); mpz_mul(t1,t1,Z);
	mpz_mod(u1,t1,n);
	/* obtain U2 */
	mpz_mul(t,u,s); mpz_mul(t1,r,v); mpz_sub(t,t,t1); mpz_mod(t,t,n);
	mpz_mul(t,t,t); mpz_mul(t,t,X); mpz_mod(U2,t,n);
	mpz_set(U1,u1);
	mpz_clear(t); mpz_clear(t1); mpz_clear(t2); mpz_clear(u1);
}

/* multiply (X,Z) with k in y^2=x^3+ax+b mod n */
void nextvalues(mpz_t X,mpz_t Z,int p,mpz_t n,mpz_t a,mpz_t b) {
	mpz_t X1,Z1,X2,Z2,U1,U2,t;
	int len=0,c[32],i;
	mpz_init(U1); mpz_init(U2); mpz_init(t);
	mpz_init_set(X1,X); mpz_init_set(Z1,Z); mpz_init(X2); mpz_init(Z2);
	while(p) c[len++]=p&1,p>>=1;
	sub_2i(X,Z,n,a,b,X2,Z2);
	/* len-1 or len-2? */
	for(i=len-2;i>=0;i--) {
		sub_2i_plus(X1,Z1,X2,Z2,X,Z,n,a,b,U1,U2);
		if(!c[i]) {
			sub_2i(X1,Z1,n,a,b,t,Z1);
			mpz_set(X1,t); mpz_set(X2,U1); mpz_set(Z2,U2);
		} else {
			sub_2i(X2,Z2,n,a,b,t,Z2);
			mpz_set(X2,t); mpz_set(X1,U1); mpz_set(Z1,U2);
		}
	}
	mpz_set(X,X1); mpz_set(Z,Z1);
	mpz_clear(U1); mpz_clear(U2); mpz_clear(t);
	mpz_clear(X1); mpz_clear(Z1); mpz_clear(X2); mpz_clear(Z2);
}

/* elliptic curve method! returns 1 and factor in out if factor is found,
   otherwise return 0, maxb is max prime, maxc is number of curves to test */
int ECM(mpz_t n,mpz_t out,int maxb,int maxc) {
	int r=0,p;
	mpz_t X,Y,Z,a,b,g,t;
	mpz_init(X); mpz_init(Y); mpz_init(Z); mpz_init(a); mpz_init(b); mpz_init(g); mpz_init(t);
	gmp_printf("start ecm on %Zd with b1 %d curves %d\n",n,maxb,maxc);
	while(maxc--) {
	newcurve:
		/* pick random curve */
		mpz_urandomm(X,gmprand,n); mpz_urandomm(Y,gmprand,n); mpz_urandomm(a,gmprand,n);
		/* determine b=(y^2-x^3-ax)%n */
		mpz_mul(b,Y,Y); mpz_powm_ui(t,X,3,n); mpz_sub(b,b,t);
		mpz_mul(t,a,X); mpz_sub(b,b,t); mpz_mod(b,b,n);
//		gmp_printf("curve a %Zd\n      b %Zd\n      x %Zd\n      y %Zd\n",a,b,X,Y);
		/* check gcd(4a^3+27b^2,n) */
		mpz_powm_ui(g,a,3,n); mpz_mul_si(g,g,4); mpz_mul(t,b,b);
		mpz_mul_si(t,t,27); mpz_add(g,g,t); mpz_gcd(g,g,n);
		if(!mpz_cmp(g,n)) goto newcurve;
		if(mpz_cmp_ui(g,1)>0) {
			/* factor found! */
			r=1;
			mpz_set(a,g);
			goto end;
		}
		mpz_set_ui(Z,1);
		for(p=2;p<=maxb;p++) {
			nextvalues(X,Z,p,n,a,b);
			if(!(p&15)) {
				mpz_gcd(g,Z,n);
				if(!mpz_cmp(g,n)) goto fail;
				if(mpz_cmp_ui(g,1)>0) {
					r=1;
					mpz_set(out,g);
					goto end;
				}
			}
		}
		/* loop done, check again */
		mpz_gcd(g,Z,n);
		if(mpz_cmp_ui(g,1)>0) {
			r=1;
			mpz_set(out,g);
			goto end;
		}
	fail:;
	}
end:
	mpz_clear(X); mpz_clear(Z); mpz_clear(a); mpz_clear(b); mpz_clear(g); mpz_clear(t);
	return r;
}

void try(mpz_t n) {
	mpz_t a;
	mpz_init(a);
	if(ECM(n,a,2000,25)) goto done;
	if(ECM(n,a,11000,90)) goto done;
	if(ECM(n,a,50000,300)) goto done;
	if(ECM(n,a,250000,700)) goto done;
	if(ECM(n,a,1000000,1800)) goto done;
	puts("no factor found");
	goto fail;
done:
	gmp_printf("found factor %Zd\n",a);
fail:
	mpz_clear(a);
}



int main() {
	mpz_t n;
	gmp_randinit_mt(gmprand);
	mpz_init(n);
	mpz_set_str(n,"189029013605764030727921585951",10);
	try(n);
	mpz_set_str(n,"74520163184103070906530082210517",10);
	try(n);
	mpz_set_str(n,"523221436353855391814506581063557",10);
	try(n);
	mpz_set_str(n,"1564875138070655023123959837084599",10);
	try(n);
	mpz_set_str(n,"78325683705012095897299536068804821",10);
	try(n);
	mpz_set_str(n,"228264844518616987380835399399539853",10);
	try(n);
	mpz_set_str(n,"7511663247147032357037656316584448877",10);
	try(n);
	mpz_set_str(n,"25348924873403921164412907702279733193",10);
	try(n);
	mpz_set_str(n,"208105107011856763735887399456439331987",10);
	try(n);
	mpz_set_str(n,"3565260354721980199129400248402571306803",10);
	try(n);
/*	mpz_set_str(n,"32160137412888834732051225949878741400809992284289",10);
	try(n);
	mpz_set_str(n,"160967735740568108627966290684899321608893044314961348169843",10);
	try(n);
	mpz_set_str(n,"216564934649977779183332104119550719684705171673087005634079598195092857334543",10);
	try(n);
*/
	mpz_clear(n);
	return 0;
}
