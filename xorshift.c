/*	“An experimental exploration of Marsaglia’s xorshift generators,
 *	scrambled”, Sebastiano Vigna */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <unistd.h>

#include "xorshift.h"

uint64_t xorshift_step (randctx * const st) {
	uint64_t s0 = st->s[ st->p ];
	uint64_t s1 = st->s[ st->p = ( st->p + 1 ) & (XORSHIFT_S-1) ];
	s1 ^= s1 << 31; // a
	s1 ^= s1 >> 11; // b
	s0 ^= s0 >> 30; // c
	return ( st->s[ st->p ] = s0 ^ s1 ) * 1181783497276652981LL;
}

static uint64_t rand64 () {
	uint64_t rand;
	while (!__builtin_ia32_rdrand64_step (&rand));
	return rand;
}

void xorshift_seed (randctx * const st) {
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

