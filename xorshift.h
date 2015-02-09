#pragma once

#include <stdint.h>

#define XORSHIFT_S 16

typedef struct {
	uint64_t s[XORSHIFT_S];
	int p;
} randctx;

uint64_t xorshift_step (randctx * const st);
void xorshift_seed (randctx * const st);

