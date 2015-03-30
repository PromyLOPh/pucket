/*
    Copyright (C) 1992-2009 Spotworks LLC

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <stdbool.h>
#include "vector.h"
#include "flam3.h"

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

