#pragma once

#include <stdint.h>

#define XORSHIFT_S 16

typedef struct {
	uint64_t s[XORSHIFT_S];
	int p;
} randctx;

void rand_seed (randctx * const st);
uint64_t rand_u64 (randctx * const st);
double rand_d01 (randctx * const st);
double rand_d11 (randctx * const st);
int rand_bool (randctx * const st);

#define vlen(x) (sizeof(x)/sizeof(*x))
#define rand_distrib(st,v) ((v)[rand_u64(st)%vlen(v)])
#define rand_mod(st,max) (rand_u64(st)%(max))
                                   
