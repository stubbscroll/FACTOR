some factorization algorithms

rho.c            pollard rho with brent's improvements.
pm1.c            pollard p-1 with one large prime.
ecm.c            basic ECM with no stage 2 implemented from "factorization
                 and primality testing" by bressoud.
ecm2.c           inversionless ECM with montgomery representation and stage 2
                 with one large prime implemented from "prime numbers - a
                 computational perspective" by crandall, pomerance (algorithm
                 7.4.4). around 3x faster than ecm.c, but around 1.5x slower
                 than alpertron (measured on a very small set of numbers).
qs.c             basic quadratic sieve, the only improvement is negative numbers.
                 using polynomial x^2-n with sieve interval centered at 
                 x=ceil(sqrt(n)), and gaussian elimination with dense matrix.
qs2.c            qs.c with large prime variation.

in.txt           file with large composites, some have more than 2 prime
                 factors. all numbers are factors of mersenne numbers. none
                 of the numbers have prime factors less than 1 million.

improvements to quadratic sieve will come later. i want to add large primes,
two large primes, multiple polynomials, self-initialization and better linear
algebra (block lanczos or block wiedemann or structured gaussian elimination).
number field sieve will probably not be added here, as i probably won't be able
to factor numbers above the crossover point where nfs beats qs. check
https://github.com/stubbscroll/nfs instead.

usage: run program and enter number(s) to factorize to stdin
