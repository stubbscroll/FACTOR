#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

/* quadratic sieve!
   this is a basic version without anything */
/* for improvements, read:
   http://mersennewiki.org/index.php/Self-Initializing_Quadratic_Sieve
   "prime numbers - a computational perspective" (crandall, pomerance)

   code examples: (implementation of siqs)
   http://gforge.inria.fr/projects/ecm/
   http://codes-sources.commentcamarche.net/source/25981-self-initializing-quadratic-sieve
   http://www.alpertron.com.ar/ECM.HTM
*/

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

unsigned int powmod(unsigned int n,unsigned int k,unsigned int mod) {
	unsigned int ans=1;
	while(k) {
		if(k&1) ans=(unsigned long long)ans*n%mod;
		k>>=1;
		n=(unsigned long long)n*n%mod;
	}
	return ans;
}

int rand31() {
	return (rand()&32767)+((rand()&32767)<<15)+((rand()&1)<<30);
}

/* find square root of a modulo p (p prime) */
int sqrtmod(int a,int p) {
	int p8,alpha,i;
	int x,c,s,n,b,J,r2a,r,ret;
	mpz_t za,zp;
	if(p==2) return a&1;
	mpz_init_set_si(za,a);
	mpz_init_set_si(zp,p);
	if(mpz_jacobi(za,zp)!=1) { /* no square root */
		ret=0;
		goto end;
	}
	p8=p&7;
	if(p8==3 || p8==5 || p8==7) {
		if((p8&3)==3) {
			ret=powmod(a,(p+1)/4,p);
			goto end;
		}
		x=powmod(a,(p+3)/8,p);
		c=(ll)x*x%p;
		ret=c==a?x:(ll)x*powmod(2,(p-1)/4,p)%p;
		goto end;
	}
	alpha=0;
	s=p-1;
	while(!(s&1)) s>>=1,alpha++;
	r=powmod(a,(s+1)/2,p);
	r2a=(ll)r*powmod(a,(s+1)/2-1,p)%p;
	do {
		n=rand31()%(p-2)+2;
		mpz_set_si(za,n);
	} while(mpz_jacobi(za,zp)!=-1);
	b=powmod(n,s,p);
	J=0;
	for(i=0;i<alpha-1;i++) {
		c=powmod(b,2*J,p);
		c=(ll)r2a*c%p;
		c=powmod(c,1<<(alpha-i-2),p);
		if(c==p-1) J+=1<<i;
	}
	ret=(ll)r*powmod(b,J,p)%p;
end:
	mpz_clear(zp);
	mpz_clear(za);
	return ret;
}

/* number of extra relations in linear algebra */
#define EXTRAREL 10
/* cache block, experiment with this */
#define BLOCKSIZE 30000
/* slack in logarithm for trial division */
/* 13 is good for 42 digits, 1.5% of trial divisions are non-smooth, but misses
   some smooth numbers */
#define LGSLACK 13

/* global variables for quadratic sieve */
struct {
	/* basic info */
	mpz_t n;                /* number to factorize */
	/* factor base */
	int B;                  /* factor base bound */
	int *p;                 /* primes in factor base */
	int *a;                 /* root of n mod p */
	int *lg;                /* integer logarithm base 2 */
	int fn;                 /* number of primes in factor base */
	/* sieve */
	short sieve[BLOCKSIZE]; /* the sieve */
	mpz_t *rel;             /* smooth numbers */
	int rn;                 /* number of smooth numbers */
	/* matrix */
	unsigned long long **m;
	/* final assembly */
	int *ev;                /* cumulative exponent vector */
} qs;

#define SETBIT(a,b) qs.m[(a)][(b)>>6]|=(1ULL<<(b&63))
#define XORBIT(a,b) qs.m[(a)][(b)>>6]^=(1ULL<<(b&63))
#define ISSET(a,b) (qs.m[(a)][(b)>>6]&(1ULL<<(b&63)))

/* gaussian elimination mod 2 on bitmasks, A is n*m, b is n*o */
/* a is a malloced array of pointers, each a[i] is of size
   sizeof(unsigned long long)*(m+o+63)/64 */
/* return 0: no solutions, 1: one solution, 2: free variables */
int bitgauss64(int n,int m,int o) {
	int i,j,k,z=m+o,c=0,fri=0,bz=(z+63)>>6;
	unsigned long long t;
	/* process each column */
	for(i=0;i<m;i++) {
		/* TODO check words instead of bits */
		for(j=c;j<n;j++) if(ISSET(j,i)) break;
		if(j==n) { fri=1; continue; }
		/* swap? */
		if(j>c)  for(k=0;k<bz;k++) {
			t=qs.m[j][k],qs.m[j][k]=qs.m[c][k],qs.m[c][k]=t;
		}
		/* subtract multiples of this row */
		for(j=0;j<n;j++) if(j!=c && ISSET(j,i)) {
			for(k=0;k<bz;k++) qs.m[j][k]^=qs.m[c][k];
		}
		c++;
	}
	/* detect no solution: rows with 0=b and b!=0 */
	for(i=0;i<n;i++) {
		/* TODO make bit-efficient solution later */
		for(j=0;j<m;j++) if(ISSET(i,j)) goto ok;
		for(;j<z;j++) if(ISSET(i,j)) return 0;
	ok:;
	}
	return 1+fri;
}

void QSgenfactorbase() {
	int i;
	mpz_t t;
	mpz_init(t);
	qs.fn=0;
	/* generate factor base! p must be a quadratic residue of n */
	for(qs.fn=1,i=0;i<primes && prime[i]<qs.B;i++) {
		mpz_set_ui(t,prime[i]);
		if(mpz_jacobi(qs.n,t)>0) qs.fn++;
	}
	/* allocate memory while we're at it */
	if(!(qs.p=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(qs.a=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(qs.lg=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(qs.rel=malloc(sizeof(mpz_t)*(qs.fn+EXTRAREL)))) puts("out of memory"),exit(1);
	for(i=0;i<qs.fn;i++) mpz_init(qs.rel[i]);
	if(!(qs.m=malloc(sizeof(unsigned long long *)*(qs.fn+EXTRAREL)))) puts("out of memory"),exit(1);
	for(i=0;i<qs.fn+EXTRAREL;i++) if(!(qs.m[i]=calloc(((qs.fn+EXTRAREL+63)/64),sizeof(unsigned long long)))) puts("out of memory"),exit(1);
	qs.fn=0;
	qs.p[qs.fn++]=-1;
	for(i=0;i<primes && prime[i]<qs.B;i++) {
		mpz_set_ui(t,prime[i]);
		if(mpz_jacobi(qs.n,t)>0) {
			qs.p[qs.fn]=prime[i];
			/* find square root a = n (mod p) */
			qs.a[qs.fn]=sqrtmod(mpz_fdiv_ui(qs.n,prime[i]),prime[i]);
			qs.lg[qs.fn]=log(prime[i])/log(2)+.5;
			qs.fn++;
		}
	}
	printf("factor base bound %d primes %d\n",qs.B,qs.fn);
	mpz_clear(t);
}

int correct,false,sieves;

void QStrialdiv(mpz_t start) {
	mpz_t t,t1,t2;
	int lim,i,j,lo,hi,mid,k;
	mpz_init_set(t,start); mpz_init(t1); mpz_init(t2);
	/* find center value in interval, t=start+BLOCKSIZE/2 */
	mpz_add_ui(t,t,BLOCKSIZE/2);
	/* take lg(t^2-n) */
	mpz_mul(t,t,t);
	mpz_sub(t,t,qs.n);
	mpz_abs(t,t);
	lim=0.5+log(mpz_get_d(t))/log(2);
	for(i=0;i<BLOCKSIZE;i++) if(qs.sieve[i]>=lim-LGSLACK) {
		mpz_set(t,start);
		mpz_add_ui(t,t,i);
		mpz_mul(t,t,t);
		mpz_sub(t,t,qs.n);
		mpz_set(t2,t);
		/* put -1 in exponent vector if negative */
		if(mpz_cmp_ui(t,0)<0) SETBIT(0,qs.rn);
		mpz_abs(t,t);
		for(j=1;j<qs.fn;j++) {
			if(qs.p[j]<46340) {
				if(mpz_cmp_ui(t,qs.p[j]*qs.p[j])<0) break;
			} else {
				mpz_set_ui(t1,qs.p[j]);
				mpz_mul(t1,t1,t1);
				if(mpz_cmp(t,t1)<0) break;
			}
			while(!mpz_fdiv_ui(t,qs.p[j])) {
				mpz_fdiv_q_ui(t,t,qs.p[j]);
				/* put factor in exponent vector */
				XORBIT(j,qs.rn);
			}
		}
		if(mpz_cmp_ui(t,1)>0 && mpz_cmp_ui(t,qs.p[qs.fn-1])<0) {
			/* find index of remaining factor in list */
			lo=j+1; hi=qs.fn; k=mpz_get_ui(t);
			while(lo<hi) {
				mid=lo+(hi-lo)/2;
				if(qs.p[mid]<k) lo=mid+1;
				else hi=mid;
			}
			/* put factor in exponent vector */
			XORBIT(lo,qs.rn);
			mpz_set_ui(t,1);
		}
		if(!mpz_cmp_ui(t,1)) {
			mpz_set(qs.rel[qs.rn],start);
			mpz_add_ui(qs.rel[qs.rn],qs.rel[qs.rn],i);
			qs.rn++;
			correct++;
		//gmp_printf("add smooth %Zd %d/%d sieve-lg %d lim %d\n",t2,qs.rn,qs.fn+EXTRAREL,qs.sieve[i],lim);
		}
		else false++;
		if(qs.rn==qs.fn+EXTRAREL) break;
	}
	mpz_clear(t); mpz_clear(t1); mpz_clear(t2);
}

/* TODO possible improvements: calculate pfrontm from pfrontp (same wich
   pback) and eliminate the 2 m-arrays */
void QSsieve() {
	int *pfrontp,*pfrontm,*pbackp,*pbackm;
	int i;
	mpz_t xfront,xback,t;
	if(!(pfrontp=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(pfrontm=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(pbackp=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(pbackm=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	mpz_init(xfront); mpz_init(xback); mpz_init(t);
	/* start sieving from ceil(sqrt(n)), in both directions */
	mpz_sqrt(xfront,qs.n);
	mpz_add_ui(xfront,xfront,1);
	mpz_set(xback,xfront);
	mpz_sub_ui(xback,xback,BLOCKSIZE);
	qs.rn=0;
	/* find offsets for front and back sieves, both roots */
	for(i=1;i<qs.fn;i++) {
		pfrontp[i]=qs.a[i]-mpz_fdiv_ui(xfront,qs.p[i]);
		while(pfrontp[i]<0) pfrontp[i]+=qs.p[i];
		pfrontm[i]=-qs.a[i]-mpz_fdiv_ui(xfront,qs.p[i]);
		while(pfrontm[i]<0) pfrontm[i]+=qs.p[i];
		pbackp[i]=pfrontp[i]-qs.p[i];
		pbackm[i]=pfrontm[i]-qs.p[i];
		if(pfrontp[i]<0) printf("error");
		if(pfrontm[i]<0) printf("error");
	}
	/* sieve! */
	correct=false=sieves=0;
	do {
		/* forward */
		for(i=0;i<BLOCKSIZE;i++) qs.sieve[i]=0;
		i=1;
		if(qs.p[i]==2) {
			while(pfrontp[i]<BLOCKSIZE) qs.sieve[pfrontp[i]]+=qs.lg[i],pfrontp[i]+=qs.p[i];
			pfrontp[i]-=BLOCKSIZE;
			i++;
		}
		for(;i<qs.fn;i++) {
			while(pfrontp[i]<BLOCKSIZE) qs.sieve[pfrontp[i]]+=qs.lg[i],pfrontp[i]+=qs.p[i];
			pfrontp[i]-=BLOCKSIZE;
			while(pfrontm[i]<BLOCKSIZE) qs.sieve[pfrontm[i]]+=qs.lg[i],pfrontm[i]+=qs.p[i];
			pfrontm[i]-=BLOCKSIZE;
		}
		/* find smooth numbers (positive) */
		if(qs.rn<qs.fn+EXTRAREL) QStrialdiv(xfront);
		mpz_add_ui(xfront,xfront,BLOCKSIZE);
		/* backward */
		for(i=0;i<BLOCKSIZE;i++) qs.sieve[i]=0;
		i=1;
		if(qs.p[i]==2) {
			pbackp[i]+=BLOCKSIZE;
			while(pbackp[i]>=0) qs.sieve[pbackp[i]]+=qs.lg[i],pbackp[i]-=qs.p[i];
			i++;
		}
		for(;i<qs.fn;i++) {
			pbackp[i]+=BLOCKSIZE;
			while(pbackp[i]>=0) qs.sieve[pbackp[i]]+=qs.lg[i],pbackp[i]-=qs.p[i];
			pbackm[i]+=BLOCKSIZE;
			while(pbackm[i]>=0) qs.sieve[pbackm[i]]+=qs.lg[i],pbackm[i]-=qs.p[i];
		}
		/* find smooth numbers (negative) */
		if(qs.rn<qs.fn+EXTRAREL) QStrialdiv(xback);
		mpz_sub_ui(xback,xback,BLOCKSIZE);
		sieves++;
		if(sieves%1000==0) printf("[%d] ",qs.rn);
	} while(qs.rn<qs.fn+EXTRAREL);
	printf("%d rel %d fail trial division %d sieve blocks\n",correct,false,sieves);
	free(pfrontp); free(pfrontm); free(pbackp); free(pbackm);
	mpz_clear(xfront); mpz_clear(xback); mpz_clear(t);
}

/* build final exponent vector by trial division of x^2-n */
void QSbuildev(mpz_t x) {
	mpz_t t,t1;
	int i,lo,hi,mid,k;
	mpz_init_set(t,x); mpz_init(t1);
	mpz_mul(t,t,t);
	mpz_sub(t,t,qs.n);
	if(mpz_cmp_ui(t,0)<0) {
		qs.ev[0]++;
		mpz_abs(t,t);
	}
	for(i=1;i<qs.fn;i++) {
		if(qs.p[i]<46340) {
			if(mpz_cmp_ui(t,qs.p[i]*qs.p[i])<0) break;
		} else {
			mpz_set_ui(t1,qs.p[i]);
			mpz_mul(t1,t1,t1);
			if(mpz_cmp(t,t1)<0) break;
		}
		while(!mpz_fdiv_ui(t,qs.p[i])) {
			mpz_fdiv_q_ui(t,t,qs.p[i]);
			qs.ev[i]++;
		}
	}
	if(mpz_cmp_ui(t,1)>0 && mpz_cmp_ui(t,qs.p[qs.fn-1])<0) {
		/* find index of remaining factor in list */
		lo=i+1; hi=qs.fn; k=mpz_get_ui(t);
		while(lo<hi) {
			mid=lo+(hi-lo)/2;
			if(qs.p[mid]<k) lo=mid+1;
			else hi=mid;
		}
		mpz_set_ui(t,1);
		qs.ev[lo]++;
	}
	if(mpz_cmp_ui(t,1)) {
		gmp_printf("sanity error, remaining factor %Zd\n",t);
		exit(1);
	}
	mpz_clear(t); mpz_clear(t1);
}

/* examine each solution vector and try to get factor */
int QSroot(mpz_t a) {
	int i,j,k,r=0,f;
	char *freevar,*v;
	mpz_t x,y;
	mpz_init(x); mpz_init(y);
	if(!(qs.ev=malloc(sizeof(int)*qs.rn))) puts("out of memory"),exit(1);
	/* find all free variables. variable i is free if there is no row having
	   its first 1-element in column i */
	if(!(freevar=malloc(qs.rn))) puts("out of memory"),exit(1);
	if(!(v=malloc(qs.rn))) puts("out of memory"),exit(1);
	for(i=0;i<qs.rn;i++) freevar[i]=1;
	for(i=0;i<qs.fn;i++) {
		for(j=0;j<qs.rn;j++) if(ISSET(i,j)) {
			freevar[j]=0;
			break;
		}
	}
	for(f=0;f<qs.fn;f++) if(freevar[f]) {
		/* set free variable i to 1 and the others to 0 */
		for(i=0;i<qs.rn;i++) v[i]=i==f;
		/* solution vector by back-substitution! set the first 1-element to
		   the xor of the others */
		for(i=qs.fn-1;i>=0;i--) {
			for(j=0;j<qs.rn;j++) if(ISSET(i,j)) goto ok;
			continue;
		ok:
			for(k=j++;j<qs.rn;j++) if(ISSET(i,j) && v[j]) v[k]^=1;
		}
		/* v[i]=1 means that i-th relation is part of the solution */
		/* take square root of left side, the product of x^2 */
		mpz_set_ui(x,1);
		for(i=0;i<qs.rn;i++) if(v[i]) mpz_mul(x,x,qs.rel[i]),mpz_mod(x,x,qs.n);
		/* take square root of right side, the product of (x^2-n) */
		for(i=0;i<qs.rn;i++) qs.ev[i]=0;
		/* we didn't want to spend lots of memory storing the factorization of
		   each x^2-n, so trial divide again */
		for(i=0;i<qs.rn;i++) if(v[i]) QSbuildev(qs.rel[i]);
		mpz_set_ui(y,1);
		/* multiply half the exponents */
		for(i=0;i<qs.rn;i++) for(j=0;j<qs.ev[i];j+=2) {
			mpz_mul_si(y,y,qs.p[i]);
			mpz_mod(y,y,qs.n);
		}
		gmp_printf("free var %d:\n  x = %Zd\n  y = %Zd\n",f,x,y);
		/* get factor */
		mpz_sub(x,x,y);
		mpz_gcd(x,x,qs.n);
		gmp_printf("  gcd %Zd\n",x);
		if(mpz_cmp_ui(x,1)>0 && mpz_cmp(x,qs.n)<0) {
			mpz_set(a,x);
			r=1;
			break;
		}
	}
	free(v);
	free(freevar);
	free(qs.ev);
	return r;
}

/* quadratic sieve! */
/* warning, don't invoke on even numbers or powers */
int QS(mpz_t n,mpz_t a) {
	double L=mpz_get_d(n);
	int r,i;
	/* the multiplier is tweakable! */
	qs.B=(int)exp(0.5*sqrt(log(L)*log(log(L))))*1.9;
	mpz_init_set(qs.n,n);
	QSgenfactorbase();
	QSsieve();
	bitgauss64(qs.fn,qs.rn,0);
	r=QSroot(a);
	for(i=0;i<qs.rn;i++) mpz_clear(qs.rel[i]);
	free(qs.p); free(qs.a); free(qs.lg); free(qs.rel);
	for(i=0;i<qs.rn;i++) free(qs.m[i]);
	free(qs.m);
	mpz_clear(qs.n);
	return r;
}

void try(mpz_t n) {
	mpz_t a;
	mpz_init(a);
	if(QS(n,a)) {
		gmp_printf("factor %Zd\n",a);
		goto end;
	}
	puts("no factor found");
end:
	mpz_clear(a);
}

int main() {
	mpz_t n;
	createsieve();
	genprimes();
	mpz_init(n);
	mpz_set_str(n,"174224571863520493293247799005065324265471",10);
	try(n);
	mpz_set_str(n,"216564934649977779183332104119550719684705171673087005634079598195092857334543",10);
	try(n);
	mpz_clear(n);
	return 0;
}
