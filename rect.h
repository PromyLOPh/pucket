#pragma once
#include <stdbool.h>

#include "vector.h"

typedef struct {
	/* bucket width/height */
	uint2 dim;
	double4 *data;
	unsigned long int badvals, samples;
} bucket;

void bucket_init (bucket * const b, const uint2 dim);
bool bucket_deserialize (bucket * const b, const char *file);
void bucket_serialize (bucket * const b, const char *file);

bool render_bucket (flam3_genome * const genome, bucket * const bucket,
		const double timelimit);
void render_image (const flam3_genome * const genome, const bucket * const b,
		void * const out, const unsigned int bytes_per_channel);

