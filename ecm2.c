#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

/* ECM from "prime numbers - a computational perspective" (crandall, pomerance) */
/* sample factorizations in main take 469 seconds */

gmp_randstate_t gmprand;

#define MAXP 100000000
typedef long long ll;

char sieve[MAXP];
int prime[5761456];
int primes;

void createsieve() {
	int i,j,q;
	memset(sieve,1,sizeof(sieve));
	q=sqrt(MAXP);
	for(sieve[0]=sieve[1]=0,i=2;i<=q;i++) if(sieve[i]) for(j=i*i;j<MAXP;j+=i) sieve[j]=0;
}

void genprimes() {
	int i;
	for(primes=i=0;i<MAXP;i++) if(sieve[i]) prime[primes++]=i;
}

/* global variables for ECM */
struct {
	mpz_t n,C;
	/* temp variables */
	mpz_t t,t1,t2,t3,t4;
	mpz_t U,V,T,W;
} ecm;

void ECMaddh(mpz_t X1,mpz_t Z1,mpz_t X2,mpz_t Z2,mpz_t X,mpz_t Z,mpz_t Xo,mpz_t Zo) {
	/* Xo = Z*(X1*X2 - Z1*Z2)^2 */
	mpz_mul(ecm.t1,X1,X2); mpz_submul(ecm.t1,Z1,Z2); mpz_mul(ecm.t1,ecm.t1,ecm.t1);
	mpz_mul(ecm.t1,ecm.t1,Z);
	/* Zo = X*(X1*Z2 - X2*Z1)^2 */
	mpz_mul(ecm.t2,X1,Z2); mpz_submul(ecm.t2,X2,Z1); mpz_mul(ecm.t2,ecm.t2,ecm.t2);
	mpz_mul(Zo,ecm.t2,X); mpz_mod(Zo,Zo,ecm.n);
	mpz_mod(Xo,ecm.t1,ecm.n);
}

/* double X,Y and put result in X2,Z2 */
void ECMdoubleh(mpz_t X,mpz_t Z,mpz_t X2,mpz_t Z2) {
	/* t1=X*X, t2=Z*Z */
	mpz_mul(ecm.t1,X,X); mpz_mul(ecm.t2,Z,Z);
	/* Z2 = 4*Z*(X*X*X + C*X*X*Z + X*Z*Z), build inner stuff in t3 */
	mpz_mul(ecm.t3,ecm.t1,X); mpz_mul(ecm.t4,ecm.t1,ecm.C); mpz_mul(ecm.t4,ecm.t4,Z);
	mpz_add(ecm.t3,ecm.t3,ecm.t4); mpz_mul(ecm.t4,ecm.t2,X);
	mpz_add(ecm.t3,ecm.t3,ecm.t4);
	mpz_mul(ecm.t3,Z,ecm.t3); mpz_mul_ui(Z2,ecm.t3,4); mpz_mod(Z2,Z2,ecm.n);
	/* X2=(X*X-Z*Z)^2 */
	mpz_sub(X2,ecm.t1,ecm.t2); mpz_mul(X2,X2,X2); mpz_mod(X2,X2,ecm.n);
}

void ECMmultiply(mpz_t X,mpz_t Z,int p,mpz_t X2,mpz_t Z2) {
	int b;
	if(p<2) puts("error not implemented"),exit(0);
	if(p==2) return ECMdoubleh(X,Z,X2,Z2);
	mpz_set(ecm.U,X); mpz_set(ecm.V,Z);
	ECMdoubleh(X,Z,ecm.T,ecm.W);
	for(b=30;;b--) if(p&(1<<b)) break;
	for(b--;b>=0;b--) {
		if(p&(1<<b)) {
			ECMaddh(ecm.T,ecm.W,ecm.U,ecm.V,X,Z,ecm.U,ecm.V);
			ECMdoubleh(ecm.T,ecm.W,ecm.T,ecm.W);
		} else {
			ECMaddh(ecm.U,ecm.V,ecm.T,ecm.W,X,Z,ecm.T,ecm.W);
			ECMdoubleh(ecm.U,ecm.V,ecm.U,ecm.V);
		}
	}
	if(p&1) return ECMaddh(ecm.U,ecm.V,ecm.T,ecm.W,X,Z,X2,Z2);
	ECMdoubleh(ecm.U,ecm.V,X2,Z2);
}

#define D 100

/* faster ECM from "prime numbers - a computational perspective" (crandall,
   pomerance), algorithm 7.4.4 */
/* return 1 and factor in out if factor is found, otherwise return 0.
   B1 is max prime (must be even), maxc is number of curves to test */
int ECM(mpz_t n,mpz_t out,int B1,int maxc) {
	long long q;
	int r=0,B2=100*B1,i,B,j,delta;
	mpz_t sigma,u,v,Qx,Qz,g,t,t1;
	mpz_t Sx[D],Sz[D],beta[D],Tx,Tz,Rx,Rz,alpha;
	mpz_init(t); mpz_init(t1);
	mpz_init(sigma); mpz_init(u); mpz_init(v); mpz_init(Qx); mpz_init(Qz); mpz_init(g);
	/* global variables */
	mpz_init_set(ecm.n,n); mpz_init(ecm.C);
	mpz_init(ecm.U); mpz_init(ecm.V); mpz_init(ecm.T); mpz_init(ecm.W);
	mpz_init(ecm.t); mpz_init(ecm.t1); mpz_init(ecm.t2); mpz_init(ecm.t3);
	gmp_printf("ecm %Zd B1 %d maxc %d d %d\n",n,B1,maxc,D);
	while(maxc--) {
		/* choose random curve */
		do {
			mpz_urandomm(sigma,gmprand,n); /* sigma between 6 and n-1 */
		} while(mpz_cmp_si(sigma,6)<0);
		/* u=(sigma^2-5) mod n */
		mpz_mul(u,sigma,sigma); mpz_sub_ui(u,u,5); mpz_mod(u,u,n);
		/* v=4*sigma mod n */
		mpz_mul_si(v,sigma,4); mpz_mod(v,v,n);
		/* C=((u-v)^3)(3u+v)/(4u^3v)-2) mod n */
		mpz_sub(t,v,u); mpz_powm_ui(t,t,3,n);
		mpz_mul_ui(t1,u,3); mpz_add(t1,t1,v); mpz_mul(t,t,t1); mpz_mod(t,t,n);
		mpz_powm_ui(t1,u,3,n); mpz_mul(t1,t1,v); mpz_mul_si(t1,t1,4); mpz_mod(t1,t1,n);
		mpz_invert(t1,t1,n); mpz_mul(t,t,t1); mpz_sub_ui(t,t,2); mpz_mod(ecm.C,t,n);
		/* Qx=u^3 mod n, Qy=v^3 mod n */
		mpz_powm_ui(Qx,u,3,n); mpz_powm_ui(Qz,v,3,n);
		/* perform stage 1 */
		for(j=0;j<primes && prime[j]<B1;j++) {
			for(q=1;q<=B1;q*=prime[j]) ECMmultiply(Qx,Qz,prime[j],Qx,Qz);
		}
		mpz_gcd(g,Qz,n);
		if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) {
			r=1;
			mpz_set(out,g);
			goto end;
		}
		/* perform stage 2 */
		for(i=0;i<D;i++) mpz_init(Sx[i]),mpz_init(Sz[i]),mpz_init(beta[i]);
		mpz_init(Tx); mpz_init(Tz); mpz_init(Rx); mpz_init(Rz); mpz_init(alpha);
		ECMdoubleh(Qx,Qz,Sx[0],Sz[0]);
		ECMdoubleh(Sx[0],Sz[0],Sx[1],Sz[1]);
		for(i=1;i<=D;i++) {
			if(i>2) ECMaddh(Sx[i-2],Sz[i-2],Sx[0],Sz[0],Sx[i-3],Sz[i-3],Sx[i-1],Sz[i-1]);
			mpz_mul(beta[i-1],Sx[i-1],Sz[i-1]); mpz_mod(beta[i-1],beta[i-1],n);
		}
		mpz_set_ui(g,1);
		B=B1-1;
		ECMmultiply(Qx,Qz,B-2*D,Tx,Tz);
		ECMmultiply(Qx,Qz,B,Rx,Rz);
		for(i=B;i<B2;i+=2*D) {
			mpz_mul(alpha,Rx,Rz); mpz_mod(alpha,alpha,n);
			for(;j<primes && prime[j]<=i+2*D;j++) {
				delta=(prime[j]-i)/2-1;
				if(delta>=D || delta<0) printf("error p %d i %d delta %d error\n",prime[j],i,delta),exit(0);
				/* g = g*( (Rx-Sx[delta])*(Rz+Sz[delta])-alpha+beta[delta]) mod n */
				mpz_sub(t,Rx,Sx[delta]); mpz_add(t1,Rz,Sz[delta]);
				mpz_mul(t,t1,t1); mpz_sub(t,t,alpha); mpz_add(t,t,beta[delta]);
				mpz_mul(g,g,t); mpz_mod(g,g,n);
			}
			ECMaddh(Rx,Rz,Sx[D-1],Sz[D-1],Tx,Tz,t,t1);
			mpz_set(Tx,Rx); mpz_set(Tz,Rz); mpz_set(Rx,t); mpz_set(Rz,t1);
		}
		mpz_gcd(g,g,n);
		if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) r=1,mpz_set(out,g);
		for(i=0;i<D;i++) mpz_clear(Sx[i]),mpz_clear(Sz[i]),mpz_clear(beta[i]);
		mpz_clear(Tx); mpz_clear(Tz); mpz_clear(Rx); mpz_clear(Rz); mpz_clear(alpha);
		if(r) goto end;
	}
end:
	mpz_clear(ecm.n); mpz_clear(ecm.C);
	mpz_clear(ecm.t); mpz_clear(ecm.t1); mpz_clear(ecm.t2); mpz_clear(ecm.t3);
	mpz_clear(ecm.U); mpz_clear(ecm.V); mpz_clear(ecm.T); mpz_clear(ecm.W);
	mpz_clear(sigma); mpz_clear(u); mpz_clear(v); mpz_clear(Qx); mpz_clear(Qz); mpz_clear(g);
	mpz_clear(t); mpz_clear(t1);
	return r;
}
#undef D

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
	createsieve();
	genprimes();
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
	try(n);*/
	mpz_clear(n);
	return 0;
}
