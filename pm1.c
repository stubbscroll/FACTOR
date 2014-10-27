#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <sys/time.h>

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

/* pollard p-1 with stage 2 with large prime. hope that n has a factor p such
   that p-1 is B1-smooth with at most one factor larger than B1 and not larger
   than B2. return 0 if it fails, otherwise return 1 and return factor in a. */
int pollardp1(mpz_t n,mpz_t a,int B1,int B2,int maxc) {
	int b,c,r=0,q;
	mpz_t m,g;
	mpz_init(m); mpz_init(g);
	for(c=2;c<2+maxc;c++) {
		mpz_set_ui(m,c);
		for(b=2;b<=B1;b++) {
			mpz_powm_ui(m,m,b,n);
			/* check for factor periodically */
			if(!(b&1023)) {
				mpz_sub_ui(g,m,1);
				mpz_gcd(g,g,n);
				if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) {
					r=1;
					mpz_set(a,g);
					goto end;
				}
			}
		}
		/* we're done, check again */
		mpz_sub_ui(g,m,1);
		mpz_gcd(g,g,n);
		if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) {
			r=1;
			mpz_set(a,g);
			goto end;
		}
		/* stage 2! for each prime p between B1 and B2, check p*m */
		if(!(b&1)) b++;
		for(q=0;b<=B2;b+=2) if(sieve[b]) {
			mpz_powm_ui(m,m,b-q,n); q=b;
			mpz_sub_ui(g,m,1);
			mpz_gcd(g,g,n);
			if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) {
				r=1;
				mpz_set(a,g);
				goto end;
			}
		}
	}
end:
	mpz_clear(m); mpz_clear(g);
	return r;
}

int B1[]={500000,2000000,5000000,10000000,99999999};
int B2[]={1000000,5000000,20000000,99999999,99999999};

void try(char *s) {
	double start,end;
	int i;
	mpz_t n,a;
	mpz_init_set_str(n,s,10); mpz_init(a);
	gmp_printf("%Zd (%d):\n",n,strlen(s));
	start=gettime();
	for(i=0;i<5;i++) {
		printf("  try B1=%d B2=%d\n",B1[i],B2[i]);
		if(pollardp1(n,a,B1[i],B2[i],2)) goto ok;
	}
	puts("  no factor found");
	goto fail;
ok:
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
	while(scanf("%999s",s)==1) try(s);
	mpz_clear(n);
	return 0;
}
