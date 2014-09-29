some factorization algorithms

rho.c            pollard rho with brent's improvements.
p1.c             pollard p-1 with one large prime.
ecm.c            basic ECM implemented from "factorization and primality
                 testing" by bressoud.
ecm2.c           improved ECM implemented from "prime numbers - a computational
                 perspective" by crandall, pomerance (algorithm 7.4.4). around 3x
                 faster than ecm.3, but around 1.5x slower than alpertron.
qs.c             basic quadratic sieve, the only improvement is negative numbers.
                 using polynomial x^2-n and gaussian elimination with dense
                 matrix.

in.txt           file with large composites, some have more than 2 prime
                 factors. all numbers are factors of mersenne numbers.

improvements of quadratic sieve will come later. i want to add large primes,
multiple polynomials, self-initialization and better linear algebra (block
lanczos or block wiedemann). number field sieve will probably not be added
here, check https://github.com/stubbscroll/nfs instead.
