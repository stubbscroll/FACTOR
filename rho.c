#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <sys/time.h>

double gettime() {
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+t.tv_usec/1000000.;
}

/* pollard-rho! when a factor is found, return it in a. the routine will run
   until a factor is found. it is the caller's responsibility to call this
   function on composite values that aren't too huge. */
void pollardrho(mpz_t n,mpz_t a) {
	int c,j;
	mpz_t max,range,k,x1,x2,product,temp,g;
	mpz_init(max); mpz_init(range); mpz_init(k); mpz_init(x1);
	mpz_init(x2); mpz_init(product); mpz_init(temp); mpz_init(g);
	/* try x^2+c for increasing c */
	for(c=1;;c++) {
		/* birthday! probability 0.5 to find a factor p after 1.179*sqrt(p)
		   iterations. set the max number of iterations slightly higher before
		   trying another c */
		mpz_mul_ui(max,n,2);
		mpz_root(max,max,4);
		mpz_set_ui(range,1);
		mpz_set_ui(x1,2);
		mpz_set_ui(x2,4); mpz_add_ui(x2,x2,c); mpz_mod(x2,x2,n);
		mpz_set_ui(product,1);
		for(j=0;mpz_cmp_ui(max,0);mpz_sub_ui(max,max,1)) {
			for(mpz_set_ui(k,0);mpz_cmp(k,range)<0;mpz_add_ui(k,k,1)) {
				mpz_powm_ui(x2,x2,2,n); mpz_add_ui(x2,x2,c);
				if(mpz_cmp(x2,n)>-1) mpz_sub(x2,x2,n);
				mpz_sub(temp,x1,x2);
				mpz_abs(temp,temp);
				mpz_mul(product,product,temp); mpz_mod(product,product,n);
				j++;
				if(!(j&7)) {
					mpz_gcd(g,product,n);
					if(!mpz_cmp(g,n)) goto failed;
					if(mpz_cmp_si(g,1)>0) {
						mpz_set(a,g);
						goto done;
					}
					mpz_set_ui(product,1);
				}
			}
			mpz_set(x1,x2);
			mpz_mul_2exp(range,range,1);
			for(mpz_set_ui(k,0);mpz_cmp(k,range)<0;mpz_add_ui(k,k,1)) {
				mpz_powm_ui(x2,x2,2,n); mpz_add_ui(x2,x2,c);
				if(mpz_cmp(x2,n)>-1) mpz_sub(x2,x2,n);
			}
		}
	failed:;
	}
done:
	mpz_clear(max); mpz_clear(range); mpz_clear(k); mpz_clear(x1);
	mpz_clear(x2); mpz_clear(product); mpz_clear(temp); mpz_clear(g);
}

void try(char *s) {
	double start,end;
	mpz_t n,a;
	mpz_init_set_str(n,s,10); mpz_init(a);
	gmp_printf("%Zd (%d):\n",n,strlen(s));
	start=gettime();
	pollardrho(n,a);
	end=gettime()-start;
	if(end<0) end=0;
	gmp_printf("  %.3f s, found %Zd (%c) * ",end,a,mpz_probab_prime_p(a,200)?'P':'C');
	mpz_fdiv_q(a,n,a);
	gmp_printf("(%c)\n",mpz_probab_prime_p(a,200)?'P':'C');
	mpz_clear(n); mpz_clear(a);
}

int main() {
	mpz_t n;
	char s[1000];
	mpz_init(n);
	while(scanf("%999s",s)==1) try(s);
	mpz_clear(n);
	return 0;
}
