/*	Psoudo-random number generator
 *
 *	Uses xorshift1024 from “An experimental exploration of Marsaglia’s xorshift
 *	generators, scrambled”, Sebastiano Vigna
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <unistd.h>

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

/*	Seed rng with rdrand
 */
void rand_seed (randctx * const st) {
	int fd = open ("/dev/urandom", O_RDONLY);
	assert (fd != -1);
	int ret = read (fd, &st->s, sizeof (st->s));
	assert (ret != -1);
	close (fd);
	st->p = 0;
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

