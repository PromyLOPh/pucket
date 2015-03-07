/*
    FLAM3 - cosmic recursive fractal flames
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

#include <assert.h>
#include <stdlib.h>
#include <omp.h>

#include "private.h"
#include "variations.h"
#include "palettes.h"
#include "math.h"
#include "rect.h"

typedef struct {
	double timelimit;
	unsigned int sub_batch_size, fuse;
	unsigned short *xform_distrib;

	/* camera stuff */
	double ws0, wb0s0, hs1, hb1s1; /* shortcuts for indexing */
	double bounds[4]; /* Corner coords of viewable area */
	double2 rot[3]; /* Rotation transformation */
	double ppux, ppuy;
} render_constants;

/*	Lookup color [0,1]
 */
static double4 color_palette_lookup (const double color,
		const color_palette_mode mode, const flam3_palette_entry * const map,
		const unsigned int map_count) {
	assert (color >= 0.0 && color <= 1.0);

	switch (mode) {
		case PALETTE_MODE_LINEAR: {
			const double ix = color * map_count;
			const double bottomix = floor (ix);
			const double frac = ix - bottomix;
			const unsigned int intix = bottomix;

			if (intix == map_count-1) {
				return vector_d4 (map[intix].color);
			} else {
				const double4 c1 = vector_d4 (map[intix].color);
				const double4 c2 = vector_d4 (map[intix+1].color);
				return c1 * (1.0-frac) + c2 * frac;
			}
			break;
		}

		case PALETTE_MODE_STEP: {
			const unsigned int intix = nearbyint (color * (map_count-1));
			return vector_d4 (map[intix].color);
			break;
		}

		default:
			assert (0);
			break;
	}
}

static void iter_thread (flam3_genome * const input_genome,
		bucket * const bucket, const render_constants * const c,
		volatile bool * const stopped) {
	randctx rc;
	rand_seed (&rc);

	flam3_genome genome;
	flam3_copy (&genome, input_genome);

	double4 *iter_storage;
	int ret = posix_memalign ((void **) &iter_storage, sizeof (*iter_storage),
			c->sub_batch_size * sizeof (*iter_storage));
	assert (ret == 0);
	assert (iter_storage != NULL);  

	const double starttime = omp_get_wtime ();

	do {
		/* Seed iterations */
		const double4 start = (double4) {
				rand_d11(&rc),
				rand_d11(&rc),
				rand_d01(&rc),
				rand_d01(&rc),
				};

		/* Execute iterations */
		const unsigned long badcount = flam3_iterate(&genome,
				c->sub_batch_size, c->fuse, start, iter_storage,
				c->xform_distrib, &rc);

#pragma omp critical
		{
			/* Add the badcount to the counter */
			bucket->badvals += badcount;
			bucket->samples += c->sub_batch_size;

			/* Put them in the bucket accumulator */
			for (unsigned int j = 0; j < c->sub_batch_size; j++) {
				double4 p = iter_storage[j];

				if (genome.rotate != 0.0) {
					const double2 p01 = (double2) { p[0], p[1] };
					const double2 rotatedp = apply_affine (p01, c->rot);
					p[0] = rotatedp[0];
					p[1] = rotatedp[1];
				}

				/* Skip if out of bounding box or invisible */
				if (p[0] >= c->bounds[0] && p[1] >= c->bounds[1] &&
						p[0] <= c->bounds[2] && p[1] <= c->bounds[3] &&
						p[3] > 0) {
					const size_t ix = (int)(c->ws0 * p[0] - c->wb0s0) + bucket->dim[0] * (int)(c->hs1 * p[1] - c->hb1s1);
#if HAVE_BUILTIN_PREFETCH
					/* prefetch for reading (0) with no locality (0). This (partially)
					 * hides the load latency for the += operation at the end of this
					 * block */
					__builtin_prefetch (&bucket->data[ix], 0, 0);
#endif

					double4 interpcolor = color_palette_lookup (p[2],
							genome.palette_mode, input_genome->palette,
							256);

					const double logvis = p[3];
					if (logvis != 1.0) {
						interpcolor *= logvis;
					}

					bucket->data[ix] += interpcolor;
				}
			}
		}
#pragma omp master
		{
			if (omp_get_wtime () - starttime > c->timelimit) {
				*stopped = true;
			}
		}
	} while (!(*stopped));

	free (iter_storage);
}

/*	Perform clipping
 */
static double4 clip (const double4 in, const double g, const double linrange,
		const double highpow, const double vibrancy) {
	double alpha, ls;

	if (in[3] <= 0.0) {
		alpha = 0.0;
		ls = 0.0;
	} else {
		alpha = flam3_calc_alpha (in[3], g, linrange);
		ls = vibrancy * alpha / in[3];
		alpha = clamp (alpha, 0.0, 1.0);
	}

	double4 newrgb = flam3_calc_newrgb (in, ls, highpow);
	newrgb += (1.0-vibrancy) * pow_d4 (in, g);
	if (alpha > 0.0) {
		newrgb /= alpha;
	} else {
		newrgb = (double4) {0, 0, 0, 0};
	}
	newrgb[3] = alpha;
	newrgb = clamp_d4 (newrgb, 0.0, 1.0);

	return newrgb;
}

void bucket_init (bucket * const b, const uint2 dim) {
	memset (b, 0, sizeof (*b));
	b->dim = dim;

	size_t size = dim[0] * dim[1] * sizeof (*b->data);
	int ret = posix_memalign ((void **) &b->data, sizeof (*b->data), size);
	assert (ret == 0);
	assert (b->data != NULL);

	memset (b->data, 0, size);
}

/* just a random 32 bit value */
#define BUCKET_CACHE_IDENT 0x252007d2

/*	Read bucket from file
 */
bool bucket_deserialize (bucket * const b, const char *file) {
	FILE *fd = fopen (file, "r");
	if (fd == NULL) {
		return false;
	}

	uint32_t ident;
	size_t ret = fread (&ident, sizeof (ident), 1, fd);
	assert (ret == 1);
	assert (ident == BUCKET_CACHE_IDENT);

	uint32_t w, h;
	ret = fread (&w, sizeof (w), 1, fd);
	assert (ret == 1);
	ret = fread (&h, sizeof (h), 1, fd);
	assert (ret == 1);
	assert (b->dim[0] == w && b->dim[1] == h);

	uint64_t samples, badvals;
	ret = fread (&samples, sizeof (samples), 1, fd);
	assert (ret == 1);
	ret = fread (&badvals, sizeof (badvals), 1, fd);
	assert (ret == 1);
	b->samples = samples;
	b->badvals = badvals;

	ret = fread (b->data, sizeof (*b->data), w*h, fd);
	assert (ret == w*h);

	fclose (fd);

	return true;
}

/*	Write bucket into a file
 */
void bucket_serialize (bucket * const b, const char *file) {
	FILE *fd = fopen (file, "w");
	assert (fd != NULL);

	uint32_t ident = BUCKET_CACHE_IDENT;
	fwrite (&ident, sizeof (ident), 1, fd);

	assert (sizeof (b->dim[0]) >= sizeof (uint32_t));
	fwrite (&b->dim[0], sizeof (uint32_t), 1, fd);
	fwrite (&b->dim[1], sizeof (uint32_t), 1, fd);

	assert (sizeof (b->samples) >= sizeof (uint64_t));
	assert (sizeof (b->badvals) >= sizeof (uint64_t));
	fwrite (&b->samples, sizeof (uint64_t), 1, fd);
	fwrite (&b->badvals, sizeof (uint64_t), 1, fd);

	fwrite (b->data, sizeof (*b->data), b->dim[0]*b->dim[1], fd);

	fclose (fd);
}

static void compute_camera (const flam3_genome * const genome,
		const bucket * const bucket, render_constants * const c) {
	assert (genome != NULL);
	assert (bucket != NULL);
	assert (c != NULL);

	double corner0, corner1;

	const double scale = pow(2.0, genome->zoom);

	c->ppux = genome->pixels_per_unit * scale;
	c->ppuy = c->ppux;
	//ppux /=  spec->pixel_aspect_ratio;
	corner0 = genome->center[0] - bucket->dim[0] / c->ppux / 2.0;
	corner1 = genome->center[1] - bucket->dim[1] / c->ppuy / 2.0;
	c->bounds[0] = corner0;
	c->bounds[1] = corner1;
	c->bounds[2] = corner0 + bucket->dim[0] / c->ppux;
	c->bounds[3] = corner1 + bucket->dim[1] / c->ppuy;
	const double size[2] = {1.0 / (c->bounds[2] - c->bounds[0]),
							1.0 / (c->bounds[3] - c->bounds[1])};
	rotate_center ((double2) { genome->rot_center[0], genome->rot_center[1] },
			genome->rotate, c->rot);
	c->ws0 = bucket->dim[0] * size[0];
	c->wb0s0 = c->ws0 * c->bounds[0];
	c->hs1 = bucket->dim[1] * size[1];
	c->hb1s1 = c->hs1 * c->bounds[1];
}

bool render_bucket (flam3_genome * const genome, bucket * const bucket,
		const double timelimit) {
	assert (bucket != NULL);
	assert (genome != NULL);

	int ret = prepare_precalc_flags(genome);
	assert (ret == 0);

	render_constants c = {
			.fuse = 100,
			.sub_batch_size = 10000,
			.xform_distrib = flam3_create_xform_distrib(genome),
			.timelimit = timelimit,
			};
	assert (c.xform_distrib != NULL);

	/* compute camera */
	compute_camera (genome, bucket, &c);

	bool stopped = false;
#pragma omp parallel shared(stopped)
	iter_thread (genome, bucket, &c, &stopped);

	free (c.xform_distrib);

	return true;
}

void render_image (const flam3_genome * const genome, const bucket * const bucket,
		void * const out, const unsigned int bytes_per_channel) {
	assert (genome != NULL);
	assert (bucket != NULL);
	assert (bucket->data != NULL);

	const unsigned int pixels = bucket->dim[0] * bucket->dim[1];
	const unsigned int channels = 4;

	/* XXX: copied from above */
	const double scale = pow(2.0, genome->zoom);
	const double ppux = genome->pixels_per_unit * scale;
	const double ppuy = ppux;

	const double sample_density = (double) bucket->samples / (double) pixels;
	const double g = 1.0 / genome->gamma;
	const double linrange = genome->gam_lin_thresh;
	const double vibrancy = genome->vibrancy;
	/* XXX: the original formula has a factor 268/256 in here, not sure why */
	const double k1 = genome->contrast * genome->brightness;
	const double area = (double) pixels / (ppux * ppuy);
	const double k2 = 1.0 / (genome->contrast * area * sample_density);
	const double highpow = genome->highlight_power;

#pragma omp parallel for
	for (unsigned int i = 0; i < pixels; i++) {
		double4 t = bucket->data[i];

		const double ls = (k1 * log(1.0 + t[3] * k2))/t[3];

		t = t * ls;
		t = clip (t, g, linrange, highpow, vibrancy);

		const double maxval = (1 << (bytes_per_channel*8)) - 1;
		t = nearbyint_d4 (t * maxval);

		switch (bytes_per_channel) {
			case 2: {
				uint16_t * const p = &((uint16_t *) out)[channels * i];
				p[0] = t[0];
				p[1] = t[1];
				p[2] = t[2];
				p[3] = t[3];
				break;
			}

			case 1: {
				uint8_t * const p = &((uint8_t *) out)[channels * i];
				p[0] = t[0];
				p[1] = t[1];
				p[2] = t[2];
				p[3] = t[3];
				break;
			}

			default:
				assert (0);
				break;
		}
	}
}
