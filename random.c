/*	Psoudo-random number generator
 *
 *	Uses xorshift1024 from “An experimental exploration of Marsaglia’s xorshift
 *	generators, scrambled”, Sebastiano Vigna
 */

#include "random.h"

uint64_t rand_u64 (randctx * const st) {
	uint64_t s0 = st->s[ st->p ];
	uint64_t s1 = st->s[ st->p = ( st->p + 1 ) & (XORSHIFT_S-1) ];
	s1 ^= s1 << 31; // a
	s1 ^= s1 >> 11; // b
	s0 ^= s0 >> 30; // c
	return ( st->s[ st->p ] = s0 ^ s1 ) * 1181783497276652981LL;
}

/*	Generate random double in [0,1]
 */
double rand_d01 (randctx * const st) {
	return (double) rand_u64 (st) / (double) UINT64_MAX;
}

/*	Generate random double in [-1,1]
 */
double rand_d11 (randctx * const st) {
	return rand_d01 (st) * 2.0 - 1.0;
}

/*	Random boolean value in {0,1}
 */
int rand_bool (randctx * const st) {
	return rand_u64 (st) & 1;
}

/*	Generate random uint64_t with Intel’s rdrand instruction
 */
static uint64_t rand64 () {
	unsigned long long rand;
	while (!__builtin_ia32_rdrand64_step (&rand));
	return rand;
}

/*	Seed rng with rdrand
 */
void rand_seed (randctx * const st) {
	/* seed with high-quality randomness */
	for (unsigned char i = 0; i < XORSHIFT_S; i++) {
		st->s[i] = rand64 ();
	}
}

#if 0
uint64_t xorshift_step (randctx * const st) {
	uint64_t x = st->s[0];
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c
	st->s[0] = x;
	return x * 2685821657736338717LL;
}
#endif

