#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <sys/time.h>

/* quadratic sieve!
   basic version (from qs.c) with large prime variation.
   smooth numbers with one large prime less than B2 are collected, and
   sieving stops as soon as the number of full relations plus the number
   of relations obtained from merged partial relations equals the factor
   base size plus a small number. we maintain a count of merged partial
   relations using a hash table that keeps count of the number of occurrences
   of each large prime.

   problems:
   - how to detect smooth numbers with 1 large prime in sieve? how to avoid
     an excessive amount of trial division?
*/

double gettime() {
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+t.tv_usec/1000000.;
}

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

void mpz_set_ull(mpz_t m,unsigned long long a) {
	mpz_import(m,1,-1,sizeof(unsigned long long),0,0,&a);
}

/* WARNING, don't call on m > maxull */
/* TODO check for m > maxull */
ll mpz_get_ull(mpz_t m) {
	long long r;
	mpz_export(&r,0,-1,sizeof(long long),0,0,m);
	return r;
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
#define EXTRAREL 20
/* cache block, experiment with this */
#define BLOCKSIZE 30000

/* structures for partial relations */
struct partial2_s {
	mpz_t rel;               /* smooth number with one large prime */
	unsigned long long *v;   /* exponent vector mod 2 (only factor base primes) */
};

struct partial_s {
	long long p;             /* large prime */
	int p2ix;                /* index into partial array */
};

/* global variables for quadratic sieve */
struct {
	/* basic info */
	mpz_t n;                 /* number to factorize */
	/* factor base */
	int B1;                  /* upper bound for primes in factor base */
	long long B2;            /* upper bound for large prime */
	mpz_t Big2;              /* large prime bound as bigint */
	int LGSLACK;             /* base-2 log slack to accept for trial division */
	int *p;                  /* primes in factor base */
	int *a;                  /* root of n mod p */
	int *lg;                 /* integer logarithm base 2 */
	int fn;                  /* number of primes in factor base */
	/* full relations */
	short sieve[BLOCKSIZE];  /* the sieve */
	mpz_t *rel;              /* smooth numbers, this includes the "first half" of
	                            merged partial relations (actually not the smooth
	                            number itself, but the x resulting in x^2-n) */
	int rn;                  /* number of relations */
	int frn;                 /* number of full relations */
	/* partial relations */
	struct partial_s *part;  /* partial relations, can be sorted, points into
	                            part2 where most of the actual data is */
	struct partial2_s *part2;/* partial relations: smooth number (with large
	                            prime) and exponent vector */
	int prn;                 /* number of partial relations */
	int maxprn;              /* number of partial relations allocated */
	/* merged partial relations consist of 2 almost smooth numbers. the first is
	   in *rel[i], the other is *rel2[i-f] where f is the number of full
	   relations. */
	mpz_t *rel2;             /* smooth numbers, second half of merged partial
	                            relations (actually, the x giving x^2-n) */
	mpz_t *lp;               /* common large prime of merged relation */
	/* hash table for counting occurrences of large primes */
	long long *hkey;         /* keys (large prime) in hash table */
	/* TODO consider making hcount uchar if memory use becomes significant */
	int *hcount;             /* value (count) in hash table */
	long long maxhash;       /* number of elements in hash table */
	int mcnt;                /* number of relations formed from partials */
	/* matrix */
	unsigned long long **m;
	/* final assembly */
	int *ev;                 /* cumulative exponent vector */
} qs;

/* find index in hash table for element key. if it doesn't exist, create it */
long long QSgethash(long long key) {
	long long ix=key%qs.maxhash;
	while(qs.hkey[ix] && qs.hkey[ix]!=key) {
		ix++;
		if(ix==qs.maxhash) ix=0;
	}
	if(!qs.hkey[ix]) qs.hkey[ix]=key;
	return ix;
}

int QScompp(const void *A,const void *B) {
	const struct partial_s *a=A,*b=B;
	if(a->p<b->p) return -1;
	if(a->p>b->p) return 1;
	return 0;
}

#define SETBITVAL(a,b,v) qs.m[(a)][(b)>>6]=(qs.m[(a)][(b)>>6]&~(1ULL<<((b)&63)))|((v>0)<<((b)&63))
#define SETBIT(a,b) qs.m[(a)][(b)>>6]|=(1ULL<<((b)&63))
#define CLRBIT(a,b) qs.m[(a)][(b)>>6]&=~(1ULL<<((b)&63))
#define XORBIT(a,b) qs.m[(a)][(b)>>6]^=(1ULL<<((b)&63))
#define ISSET(a,b) (qs.m[(a)][(b)>>6]&(1ULL<<((b)&63)))

/* gaussian elimination mod 2 on bitmasks, A is n*m, b is n*o */
/* a is a malloced array of pointers, each a[i] is of size
   sizeof(unsigned long long)*(m+o+63)/64 */
/* TODO optimize this later */
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

int QSgenfactorbase() {
	int i,j;
	mpz_t t;
	mpz_init(t);
	qs.fn=0;
	/* generate factor base! p must be a quadratic residue of n */
	for(qs.fn=2,i=1;i<primes && prime[i]<qs.B1;i++) {
		mpz_set_ui(t,prime[i]);
		if(mpz_jacobi(qs.n,t)>0) qs.fn++;
	}
	/* allocate memory! */
	if(!(qs.p=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(qs.a=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(qs.lg=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(qs.rel=malloc(sizeof(mpz_t)*(qs.fn+EXTRAREL)))) puts("out of memory"),exit(1);
	for(i=0;i<qs.fn+EXTRAREL;i++) mpz_init(qs.rel[i]);
	if(!(qs.m=malloc(sizeof(unsigned long long *)*(qs.fn)))) puts("out of memory"),exit(1);
	for(i=0;i<qs.fn;i++) if(!(qs.m[i]=calloc(((qs.fn+EXTRAREL+63)/64),sizeof(unsigned long long)))) puts("out of memory"),exit(1);
	qs.fn=0;
	qs.p[qs.fn++]=-1;
	qs.p[qs.fn]=2; qs.lg[qs.fn]=1; qs.a[qs.fn++]=1;
	for(i=1;i<primes && prime[i]<qs.B1;i++) {
		mpz_set_ui(t,prime[i]);
		j=mpz_jacobi(qs.n,t);
		/* in the exceedingly rare event that the prime divides n: return factor */
		if(!j) return prime[i];
		if(j>0) {
			qs.p[qs.fn]=prime[i];
			/* find square root a = n (mod p) */
			qs.a[qs.fn]=sqrtmod(mpz_fdiv_ui(qs.n,prime[i]),prime[i]);
			qs.lg[qs.fn]=log(prime[i])/log(2)+.5;
			qs.fn++;
		}
	}
	printf("  factor base bound %d primes %d\n",qs.B1,qs.fn);
	mpz_clear(t);
	return 0;
}

int correct,false,sieves;

/* dir: 1 if positive x^2-n, -1 if negative */
void QStrialdiv(mpz_t start,int dir) {
	mpz_t t,t1,t2;
	long long big,ix,iy;
	int lim,i,j,lo,hi,mid,k;
	mpz_init_set(t,start); mpz_init(t1); mpz_init(t2);
	/* find lowest value in interval */
	if(dir<0) mpz_add_ui(t,t,BLOCKSIZE-1);
	/* take lg(t^2-n) */
	mpz_mul(t,t,t);
	mpz_sub(t,t,qs.n);
	mpz_abs(t,t);
	lim=0.5+log(mpz_get_d(t))/log(2);
	/* TODO find good way to modify LGSLACK for large prime */
	for(i=0;i<BLOCKSIZE;i++) if(qs.sieve[i]>=lim-qs.LGSLACK) {
		mpz_set(t,start);
		mpz_add_ui(t,t,i);
		mpz_mul(t,t,t);
		mpz_sub(t,t,qs.n);
		mpz_set(t2,t);
		/* put -1 in exponent vector if negative */
		if(mpz_cmp_ui(t,0)<0) SETBIT(0,qs.rn);
		mpz_abs(t,t);
		for(j=1;j<qs.fn;j++) {
			/* break if remainder < next prime^2 */
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
			lo=j; hi=qs.fn; k=mpz_get_ui(t);
			while(lo<hi) {
				mid=lo+(hi-lo)/2;
				if(qs.p[mid]<k) lo=mid+1;
				else hi=mid;
			}
			if(k!=qs.p[lo]) printf("sanity error, remainder %d found %d\n",k,qs.p[lo]);
			/* put factor in exponent vector */
			if(lo<0 || lo>=qs.fn) printf("lo %d out of bounds\n",lo);
			XORBIT(lo,qs.rn);
			mpz_set_ui(t,1);
		} else if(mpz_cmp_ui(t,1)>0 && mpz_cmp(t,qs.Big2)<1) {
			/* the remainder is a large prime within our B2 bounds */
			if(qs.prn==qs.maxprn) puts("ERROR, exhausted large prime storage"),exit(1);
			big=mpz_get_ull(t);
			qs.part[qs.prn].p2ix=qs.prn;
			qs.part[qs.prn].p=big;
			mpz_init_set(qs.part2[qs.prn].rel,start);
			mpz_add_ui(qs.part2[qs.prn].rel,start,i);
			if(!(qs.part2[qs.prn].v=calloc((qs.fn+63)/64,sizeof(unsigned long long)))) puts("out of memory"),exit(1);
			for(j=0;j<qs.fn;j++) if(ISSET(j,qs.rn)) qs.part2[qs.prn].v[j>>6]|=(1ULL<<(j&63));
			qs.prn++;
			/* update count */
			if(++qs.hcount[QSgethash(big)]>1) qs.mcnt++;
			/* don't set t to 1 to invoke fail routine on purpose */
//			gmp_printf("add large prime %Zd of %Zd (sieve-ix %Zd) count %d ix %d\n",t,t2,qs.part2[qs.prn-1].rel,qs.mcnt,qs.prn-1);
//			for(j=0;j<qs.fn;j++) if(qs.part2[qs.prn-1].v[j>>6]&(1ULL<<(j&63))) printf("1"); else printf("0"); printf("\n");
		}
		if(!mpz_cmp_ui(t,1)) {
			mpz_set(qs.rel[qs.rn],start);
			mpz_add_ui(qs.rel[qs.rn],qs.rel[qs.rn],i);
			qs.rn++;
			correct++;
//			gmp_printf("add smooth %Zd %d/%d sieve-lg %d lim %d\n",t2,qs.rn,qs.fn+EXTRAREL,qs.sieve[i],lim);
//			for(j=0;j<qs.fn;j++) if(ISSET(j,qs.rn-1)) printf("1"); else printf("0"); printf("\n");
		} else {
			/* factorization failed, clear vector */
			for(j=0;j<qs.fn;j++) CLRBIT(j,qs.rn);
			false++;
		}
		/* break if we have enough relations */
		if(qs.rn+qs.mcnt==qs.fn+EXTRAREL) {
			/* if we have enough partial relations, merge them and end sieving */
			if(qs.mcnt) {
				/* sort on largest prime */
				qsort(qs.part,qs.prn,sizeof(qs.part[0]),QScompp);
				/* allocate storage for merged smooth numbers */
				if(!(qs.rel2=malloc(qs.mcnt*sizeof(mpz_t)))) puts("out of memory"),exit(1);
				if(!(qs.lp=malloc(qs.mcnt*sizeof(mpz_t)))) puts("out of memory"),exit(1);
				for(i=0;i<qs.mcnt;i++) mpz_init(qs.rel2[i]),mpz_init(qs.lp[i]);
				i=qs.rn;
				/* for each prime with count>=2: merge 0-1, 0-2, ..., 0-(cnt-1) */
				for(ix=0;ix<qs.prn;) {
					for(iy=ix+1;iy<qs.prn && qs.part[ix].p==qs.part[iy].p;iy++) {
						/* merge ix and iy */
						mpz_set(qs.rel[i],qs.part2[qs.part[ix].p2ix].rel);
						mpz_set(qs.rel2[i-qs.rn],qs.part2[qs.part[iy].p2ix].rel);
						mpz_set_ull(qs.lp[i-qs.rn],qs.part[ix].p);
						for(j=0;j<qs.fn;j++) if(qs.part2[qs.part[ix].p2ix].v[j>>6]&(1ULL<<(j&63)))
							XORBIT(j,i);
						for(j=0;j<qs.fn;j++) if(qs.part2[qs.part[iy].p2ix].v[j>>6]&(1ULL<<(j&63)))
							XORBIT(j,i);
/*						gmp_printf("merge %d (%d %d) (%Zd %Zd)\n",(int)qs.part[ix].p,qs.part[ix].p2ix,qs.part[iy].p2ix,qs.part2[qs.part[ix].p2ix].rel,qs.part2[qs.part[iy].p2ix].rel);
						for(j=0;j<qs.fn;j++) if(qs.part2[qs.part[ix].p2ix].v[j>>6]&(1ULL<<(j&63))) printf("1"); else printf("0");printf("\n");
						for(j=0;j<qs.fn;j++) if(qs.part2[qs.part[iy].p2ix].v[j>>6]&(1ULL<<(j&63))) printf("1"); else printf("0");printf("\n");
						for(j=0;j<qs.fn;j++) if(ISSET(j,i)) printf("1"); else printf("0"); printf("\n");*/
						i++;
//						printf("done %d/%d\n",i,qs.fn+EXTRAREL);
					}
					ix=iy;
				}
				/* adjust the meaning of some variables... */
			}
			qs.frn=qs.rn;
			qs.rn+=qs.mcnt;
			break;
		}
	}
	mpz_clear(t); mpz_clear(t1); mpz_clear(t2);
}

/* TODO possible improvements: calculate pfrontm from pfrontp (same with
   pback) and eliminate the 2 m-arrays and trade memory access for
   calculation */
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
	printf("  ");
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
		if(qs.rn+qs.mcnt<qs.fn+EXTRAREL) QStrialdiv(xfront,1);
		mpz_add_ui(xfront,xfront,BLOCKSIZE);
		/* backward */
		if(mpz_cmp_ui(xback,0)<0) continue;
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
		if(qs.rn<qs.fn+EXTRAREL) QStrialdiv(xback,-1);
		mpz_sub_ui(xback,xback,BLOCKSIZE);
		sieves++;
		if(sieves%1000000==0) printf("[%d] ",qs.rn);
	} while(qs.rn+qs.mcnt<qs.fn+EXTRAREL);
	printf("%d fr %d pr %d mr %d fail trial division %d sieve blocks\n",correct,qs.prn,qs.rn-qs.frn,false,sieves);
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
		lo=i; hi=qs.fn; k=mpz_get_ui(t);
		while(lo<hi) {
			mid=lo+(hi-lo)/2;
			if(qs.p[mid]<k) lo=mid+1;
			else hi=mid;
		}
		if(k!=qs.p[lo]) printf("sanity error, remainder %d found %d\n",k,qs.p[lo]);
		mpz_set_ui(t,1);
		qs.ev[lo]++;
	}
	if(mpz_cmp(t,qs.Big2)>=0) {
		gmp_printf("sanity error, remaining factor %Zd\n",t);
		exit(1);
	}
	mpz_clear(t); mpz_clear(t1);
}

/* examine each solution vector and try to get factor */
int QSroot(mpz_t a) {
	int i,j,k,r=0,f,tried=0;
	char *freevar,*v;
	mpz_t x,y,t;
	mpz_init(x); mpz_init(y); mpz_init(t);
	if(!(qs.ev=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
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
	for(f=0;f<qs.rn;f++) if(freevar[f]) {
		tried++;
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
		for(i=0;i<qs.rn;i++) if(v[i]) {
			mpz_mul(x,x,qs.rel[i]),mpz_mod(x,x,qs.n);
			if(i>=qs.frn) mpz_mul(x,x,qs.rel2[i-qs.frn]);
		}
		/* take square root of right side, the product of (x^2-n) */
		for(i=0;i<qs.fn;i++) qs.ev[i]=0;
		/* we didn't want to spend lots of memory storing the factorization of
		   each x^2-n, so trial divide again */
		for(i=0;i<qs.rn;i++) if(v[i]) QSbuildev(qs.rel[i]);
		for(i=0;i<qs.rn-qs.frn;i++) if(v[i+qs.frn]) QSbuildev(qs.rel2[i]);
		mpz_set_ui(y,1);
		/* multiply half the exponents */
		for(i=0;i<qs.fn;i++) for(j=0;j<qs.ev[i];j+=2) {
			mpz_mul_si(y,y,qs.p[i]);
			mpz_mod(y,y,qs.n);
		}
		/* multiply in the large primes */
		for(i=0;i<qs.rn-qs.frn;i++) if(v[i+qs.frn]) mpz_mul(y,y,qs.lp[i]),mpz_mod(y,y,qs.n);
		/* get factor */
//		gmp_printf("try gcd %Zd %Zd =>",x,y);
		mpz_sub(x,x,y);
		mpz_gcd(x,x,qs.n);
//		gmp_printf(" %Zd\n",x);
		if(mpz_cmp_ui(x,1)>0 && mpz_cmp(x,qs.n)<0) {
			mpz_set(a,x);
			printf("  found factor after %d nullvectors\n",tried);
			r=1;
			break;
		}
	}
	free(v);
	free(freevar);
	free(qs.ev);
	mpz_clear(x); mpz_clear(y); mpz_clear(t);
	return r;
}

/* quadratic sieve! */
/* warning, don't invoke on even integers or powers */
int QS(mpz_t n,mpz_t a) {
	double L=mpz_get_d(n);
	int r,i,d=mpz_sizeinbase(n,10);
	/* TODO tweak B1 and B2 */
	qs.B1=(int)exp(0.5*sqrt(log(L)*log(log(L))));
	if(d<=50) qs.B1*=1.4;
	else if(d>=68) qs.B1*=0.9;
	else qs.B1*=1.4-((d-50)/18*0.5);
	qs.B1=qs.B1*0.5+10;
	qs.B2=(qs.B1)*100;
	/* ensure that B2<B1*B1 */
	if(qs.B2>(long long)qs.B1*qs.B1) qs.B2=(long long)qs.B1*qs.B1-1;
	mpz_init(qs.Big2);
	mpz_set_ull(qs.Big2,qs.B2);
	/* the following formula of LGSLACK found by experimentation */
	qs.LGSLACK=13;
	if(d>=33) qs.LGSLACK+=(d-33)/7;
	/* additional slack for large prime */
//	qs.LGSLACK+=log(qs.B1)/log(2);
	/* allocate hash table. TODO measure fill degree later */
	qs.maxhash=qs.B2/log(qs.B2)-qs.B1/log(qs.B1);
	if(!(qs.hkey=calloc(qs.maxhash,sizeof(long long)))) puts("out of memory"),exit(1);
	if(!(qs.hcount=calloc(qs.maxhash,sizeof(int)))) puts("out of memory"),exit(1);
	qs.mcnt=0;
	mpz_init_set(qs.n,n);
	if((i=QSgenfactorbase())) {
		/* factor base prime divides n */
		mpz_set_si(a,i);
		r=1;
		goto done;
	}
	qs.maxprn=qs.fn*10+10000;
	qs.prn=0;
	if(!(qs.part=malloc(qs.maxprn*sizeof(struct partial_s)))) puts("out of memory"),exit(1);
	if(!(qs.part2=malloc(qs.maxprn*sizeof(struct partial2_s)))) puts("out of memory"),exit(1);
	QSsieve();
	bitgauss64(qs.fn,qs.rn,0);
	r=QSroot(a);
	/* dealloc stuff */
	if(qs.mcnt) {
		for(i=0;i<qs.mcnt;i++) mpz_clear(qs.rel2[i]),mpz_clear(qs.lp[i]);
		free(qs.rel2); free(qs.lp);
	}
	for(i=0;i<qs.prn;i++) free(qs.part2[i].v);
	for(i=0;i<qs.prn;i++) mpz_clear(qs.part2[i].rel);
	free(qs.part2);
	free(qs.part);
done:
	free(qs.hcount);
	free(qs.hkey);
	for(i=0;i<qs.rn;i++) mpz_clear(qs.rel[i]);
	free(qs.p); free(qs.a); free(qs.lg); free(qs.rel);
	for(i=0;i<qs.fn;i++) free(qs.m[i]);
	free(qs.m);
	mpz_clear(qs.Big2);
	mpz_clear(qs.n);
	return r;
}

void try(char *s) {
	double start,end;
	mpz_t n,a;
	mpz_init_set_str(n,s,10); mpz_init(a);
	gmp_printf("%Zd (%d):\n",n,strlen(s));
	start=gettime();
	if(QS(n,a)) goto done;
	end=gettime()-start;
	printf("  %.3f s no factor found\n",end);
	goto fail;
done:
	end=gettime()-start;
	if(end<0) end=0;
	gmp_printf("  %.3f s, found %Zd (%c) * ",end,a,mpz_probab_prime_p(a,200)?'P':'C');
	mpz_fdiv_q(a,n,a);
	gmp_printf("(%c)\n",mpz_probab_prime_p(a,200)?'P':'C');
fail:
	mpz_clear(n); mpz_clear(a);
}

int main() {
	mpz_t n;
	char s[1000];
	mpz_init(n);
	createsieve();
	genprimes();
	while(scanf("%999s",s)==1) try(s);
	mpz_clear(n);
	return 0;
}
