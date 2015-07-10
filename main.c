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

#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <argp.h>

#include "random.h"
#include "img.h"
#include "rect.h"
#include "math.h"
#include "genome.h"
#include "palettes_builtin.h"

#define streq(a,b) (strcmp (a, b) == 0)

const char *argp_program_version = PACKAGE "-" VERSION;

typedef struct {
	unsigned int bpc;
	float scale, time;
	char *cache;
} render_arguments;

static error_t parse_render_opt (int key, char *arg,
		struct argp_state * const state) {
	render_arguments * const arguments = state->input;
	switch (key) {
		case 'b': {
			int i = atoi (arg);
			if (i == 8 || i == 16) {
				arguments->bpc = i;
			} else {
				argp_error (state, "Bits per channel must be 8 or 16");
			}
			break;
		}

		case 't': {
			float i = atof (arg);
			if (i <= 0) {
				argp_error (state, "Time must be > 0");
			} else {
				arguments->time = i;
			}
			break;
		}

		case 's':
			arguments->scale = atof (arg);
			if (arguments->scale <= 0.0) {
				argp_error (state, "Scale must be > 0");
			}
			break;

		case 'c':
			arguments->cache = strdup (arg);
			break;

		case ARGP_KEY_ARG:
			if (state->arg_num > 0) {
				return ARGP_ERR_UNKNOWN;
			}
			break;

		case ARGP_KEY_END:
			break;

		default:
			return ARGP_ERR_UNKNOWN;
			break;
	}
	return 0;
}

static void do_render (const render_arguments * const arguments) {
	randctx rc;

	rand_seed(&rc);

	int ncps;
	flam3_genome * const cps = flam3_parse_xml2 (STDIN_FILENO,
			flam3_defaults_on, &ncps, &rc);
	if (cps == NULL) {
		fprintf(stderr,"error reading genomes from file\n");
		exit(1);
	}
	assert (ncps == 1);

	flam3_genome * const genome = &cps[0];

	genome->height *= arguments->scale;
	genome->width *= arguments->scale;
	genome->pixels_per_unit *= arguments->scale;

	const unsigned int bytes_per_channel = arguments->bpc/8;
	const unsigned int channels = 4;
	const size_t this_size = channels * genome->width * genome->height *
			bytes_per_channel;
	void *image = (void *) calloc(this_size, sizeof(char));

	bucket bucket;
	bucket_init (&bucket, (uint2) { genome->width, genome->height });
	if (arguments->cache != NULL) {
		bucket_deserialize (&bucket, arguments->cache);
	}

	render_bucket (genome, &bucket, arguments->time);

	if (arguments->cache != NULL) {
		bucket_serialize (&bucket, arguments->cache);
	}

	fprintf (stderr, "%lu samples, %lu bad\n",
			bucket.samples, bucket.badvals);
	render_image (genome, &bucket, image, bytes_per_channel);

	write_png (stdout, image, genome->width, genome->height,
			bytes_per_channel);
}

static void print_genome (flam3_genome * const genome) {
	printf("<pick version=\"" PACKAGE "-" VERSION "\">\n");
	flam3_print (stdout, genome, NULL);
	printf("</pick>\n");
}

typedef struct {
	const char *palette;
	unsigned int width, height, max_xforms, max_var;
	signed int max_symmetry, min_symmetry;
	double post_likelihood, final_likelihood, symmetry_likelihood;
} random_arguments;

static error_t parse_random_opt (int key, char *arg,
		struct argp_state * const state) {
	random_arguments * const arguments = state->input;
	switch (key) {
		case 'h': {
			int i = atoi (arg);
			if (i <= 0) {
				argp_error (state, "Height must be > 0");
			} else {
				arguments->height = i;
			}
			break;
		}

		case 'w': {
			int i = atoi (arg);
			if (i <= 0) {
				argp_error (state, "Width must be > 0");
			} else {
				arguments->width = i;
			}
			break;
		}

		case 'x': {
			int i = atoi (arg);
			if (i <= 0) {
				argp_error (state, "Max xforms must be > 0");
			} else {
				arguments->max_xforms = i;
			}
			break;
		}

		case 'v': {
			int i = atoi (arg);
			if (i <= 0) {
				argp_error (state, "Max variations must be > 0");
			} else {
				arguments->max_var = i;
			}
			break;
		}

		case ARGP_KEY_ARG:
			if (state->arg_num > 0) {
				return ARGP_ERR_UNKNOWN;
			}
			break;

		case ARGP_KEY_END:
			break;

		default:
			return ARGP_ERR_UNKNOWN;
			break;
	}

	return 0;
}

#define GOLDEN_RATIO (1.618033988749894848204586834)
#define GOLDEN_RATIO_M1 (GOLDEN_RATIO-1.0)
#define GOLDEN_RATIO_DIV (GOLDEN_RATIO_M1/GOLDEN_RATIO)
static double golden_bit (randctx * const rc) {
	return rand_bool (rc) ? GOLDEN_RATIO_DIV : GOLDEN_RATIO_M1;
}

static void adjust_bounding_box (flam3_genome * const genome, randctx * const rc) {
	double bmin[2], bmax[2];
	flam3_estimate_bounding_box(genome, 0.01, 100000, bmin, bmax, rc);
	if (rand_d01(rc) < 0.3) {
		genome->center[0] = (bmin[0] + bmax[0]) / 2.0;
		genome->center[1] = (bmin[1] + bmax[1]) / 2.0;
	} else {
		double mix0, mix1;
		if (rand_bool(rc)) {
			mix0 = golden_bit(rc) + rand_d11(rc)/5;
			mix1 = golden_bit(rc);
		} else if (rand_bool(rc)) {
			mix0 = golden_bit(rc);
			mix1 = golden_bit(rc) + rand_d11(rc)/5;
		} else {
			mix0 = golden_bit(rc) + rand_d11(rc)/5;
			mix1 = golden_bit(rc) + rand_d11(rc)/5;
		}
		genome->center[0] = mix0 * bmin[0] + (1-mix0)*bmax[0];
		genome->center[1] = mix1 * bmin[1] + (1-mix1)*bmax[1];
	}
	const double zoomin = rand_d01 (rc) * 5;
	genome->pixels_per_unit = genome->width / (bmax[0] - bmin[0]) * zoomin;
}

static void do_random (const random_arguments * const arguments) {
	randctx rc;
	rand_seed(&rc);

	flam3_genome genome;
	clear_cp (&genome,flam3_defaults_on);

	genome.hue_rotation = rand_mod(&rc, 8) ? 0.0 : rand_d01(&rc);
	const palette * const p = palette_random (&builtin_palettes, &rc);
	assert (p != NULL);
	palette_copy (p, &genome.palette);
	palette_rotate_hue (&genome.palette, genome.hue_rotation);
	genome.interpolation = flam3_interpolation_linear;
	genome.rotate = rand_d01 (&rc) * 360.0;

	unsigned int nxforms = rand_mod (&rc, arguments->max_xforms) + 1;
	flam3_add_xforms(&genome,nxforms,0,0);
	/* Add a final xform 15% of the time */
	const bool add_final = rand_d01(&rc) < arguments->final_likelihood;
	if (add_final) {
		flam3_add_xforms(&genome,1,0,1);
		++nxforms;
	}

	/* Loop over xforms */
	assert (nxforms > 0);
	for (unsigned int i = 0; i < nxforms; i++) {
		flam3_xform * const xform = &genome.xform[i];
		const bool add_post = rand_d01 (&rc) < arguments->post_likelihood;
		xform_rand (xform, add_post, arguments->max_var, &rc);
		xform->density = 1.0 / nxforms;
		xform->animate = 1.0;
	}

	/* Randomly add symmetry (but not if we've already added a final xform) */
	if (rand_d01(&rc) < arguments->symmetry_likelihood && !add_final) {
		assert (arguments->max_symmetry >= arguments->min_symmetry);
		unsigned int symrange = arguments->max_symmetry - arguments->min_symmetry + 1;
		flam3_add_symmetry(&genome, rand_mod (&rc, symrange) + arguments->min_symmetry);
	}

	/* random resets genome, adjust before finding appropriate bbox */
	genome.width = arguments->width;
	genome.height = arguments->height;

	adjust_bounding_box (&genome, &rc);

	print_genome (&genome);
}

typedef struct {
	int method;
	unsigned int symmetry;
} mutate_arguments;

static error_t parse_mutate_opt (int key, char *arg,
		struct argp_state * const state) {
	mutate_arguments * const arguments = state->input;
	switch (key) {
		case 'm':
			if (arg == NULL) {
				arguments->method = MUTATE_NOT_SPECIFIED;
			} else if (streq (arg, "all-vars")) {
				arguments->method = MUTATE_ALL_VARIATIONS;
			} else if (streq(arg,"one-xform")) {
				arguments->method = MUTATE_ONE_XFORM_COEFS;
			} else if (streq(arg,"add-symmetry")) {
				arguments->method = MUTATE_ADD_SYMMETRY;
			} else if (streq(arg,"post-xforms")) {
				arguments->method = MUTATE_POST_XFORMS;
			} else if (streq(arg,"color-palette")) {
				arguments->method = MUTATE_COLOR_PALETTE;
			} else if (streq(arg,"delete-xform")) {
				arguments->method = MUTATE_DELETE_XFORM;
			} else if (streq(arg,"all-coefs")) {
				arguments->method = MUTATE_ALL_COEFS;
			} else {
				argp_error (state, "Unknown method %s", arg);
			}
			break;

		case ARGP_KEY_ARG:
			if (state->arg_num > 0) {
				return ARGP_ERR_UNKNOWN;
			}
			break;

		case ARGP_KEY_END:
			break;

		default:
			return ARGP_ERR_UNKNOWN;
			break;
	}

	return 0;
}

static void do_mutate (const mutate_arguments * const arguments) {
	randctx rc;

	rand_seed(&rc);

	int ncps;
	flam3_genome * const cps = flam3_parse_xml2 (STDIN_FILENO,
			flam3_defaults_on, &ncps, &rc);
	if (cps == NULL) {
		fprintf(stderr,"error reading genomes from file\n");
		exit(1);
	}
	assert (ncps == 1);

	flam3_genome * const genome = &cps[0];

	int ivars = 0;
	const double speed = 1.0;
	flam3_mutate (genome, arguments->method, &ivars, 1, arguments->symmetry,
			speed, &builtin_palettes, &rc);

	print_genome (genome);
}

typedef struct {
	int method;
} cross_arguments;

static error_t parse_cross_opt (int key, char *arg,
		struct argp_state * const state) {
	mutate_arguments * const arguments = state->input;
	switch (key) {
		case 'm':
			if (arg == NULL) {
				arguments->method = CROSS_NOT_SPECIFIED;
			} else if (streq(arg,"union")) {
				arguments->method = CROSS_UNION;
			} else if (streq(arg,"interpolate")) {
				arguments->method = CROSS_INTERPOLATE;
			} else if (streq(arg,"alternate")) {
				arguments->method = CROSS_ALTERNATE;
			} else {
				argp_error (state, "Unknown method %s", arg);
			}
			break;

		case ARGP_KEY_ARG:
			if (state->arg_num > 0) {
				return ARGP_ERR_UNKNOWN;
			}
			break;

		case ARGP_KEY_END:
			break;

		default:
			return ARGP_ERR_UNKNOWN;
			break;
	}

	return 0;
}

static void do_cross (const cross_arguments * const arguments) {
	randctx rc;

	rand_seed(&rc);

	int ncps;
	flam3_genome * const cps = flam3_parse_xml2 (STDIN_FILENO,
			flam3_defaults_on, &ncps, &rc);
	if (cps == NULL) {
		fprintf(stderr,"error reading genomes from file\n");
		exit(1);
	}
	assert (ncps == 2);

	flam3_genome * const genome_a = &cps[0], * const genome_b = &cps[1];
	flam3_genome genome_out;

	flam3_cross (genome_a, genome_b, &genome_out, arguments->method, &rc);

	print_genome (&genome_out);
}

typedef struct {
	unsigned int tries, resolution;
	bool change_palette;
	double time;
} improvecolors_arguments;

static error_t parse_improvecolors_opt (int key, char *arg,
		struct argp_state * const state) {
	improvecolors_arguments * const arguments = state->input;

	return 0;
}

static void do_improvecolors (const improvecolors_arguments * const arguments) {
	randctx rc;

	rand_seed(&rc);

	int ncps;
	flam3_genome * const cps = flam3_parse_xml2 (STDIN_FILENO,
			flam3_defaults_on, &ncps, &rc);
	if (cps == NULL) {
		fprintf(stderr,"error reading genomes from file\n");
		exit(1);
	}
	assert (ncps == 1);

	flam3_improve_colors (&cps[0], arguments->tries, arguments->change_palette,
			arguments->resolution, arguments->time, &builtin_palettes, &rc);

	print_genome (&cps[0]);
}

typedef struct {
	float time;
} interpolate_arguments;

static error_t parse_interpolate_opt (int key, char *arg,
		struct argp_state * const state) {
	interpolate_arguments * const arguments = state->input;
	switch (key) {
		case 't': {
			assert (arg != NULL);
			float i = atof (arg);
			if (i >= 0.0 && i <= 1.0) {
				arguments->time = i;
			} else {
				argp_error (state, "Time must be between 0 and 1");
			}	
			break;
		}

		case ARGP_KEY_ARG:
			if (state->arg_num > 0) {
				return ARGP_ERR_UNKNOWN;
			}
			break;

		case ARGP_KEY_END:
			break;

		default:
			return ARGP_ERR_UNKNOWN;
			break;
	}

	return 0;
}

static void do_interpolate (const interpolate_arguments * const arguments) {
	randctx rc;

	rand_seed(&rc);

	int ncps;
	flam3_genome * const cps = flam3_parse_xml2 (STDIN_FILENO,
			flam3_defaults_on, &ncps, &rc);
	if (cps == NULL) {
		fprintf(stderr,"error reading genomes from file\n");
		exit(1);
	}
	assert (ncps == 2);
	cps[0].time = 0.0;
	cps[1].time = 1.0;

	flam3_genome genome_out;
	flam3_interpolate (cps, 2, arguments->time, 0, &genome_out);

	print_genome (&genome_out);
}

static void show_help (const char * const argv0) {
	const char *progname = strrchr (argv0, (int) '/');
	if (progname == NULL) {
		progname = argv0;
	} else {
		++progname;
	}
	fprintf (stderr,
			"Usage: %s cross [OPTION...]\n"
			"   Or: %s improvecolors [OPTION...]\n"
			"   Or: %s interpolate [OPTION...]\n"
			"   Or: %s mutate [OPTION...]\n"
			"   Or: %s random [OPTION...]\n"
			"   Or: %s render [OPTION...]\n",
			progname, progname, progname, progname, progname);
}

int main (int argc, char **argv) {
	if (argc < 2) {
		show_help (argv[0]);
		return EXIT_FAILURE;
	}

	const char * const command = argv[1];

	if (streq (command, "cross")) {
		const struct argp_option options[] = {
				{"method", 'm', "XXX", OPTION_ARG_OPTIONAL, "Cross method" },
				{ 0 },
				};
		const char doc[] = PACKAGE "-cross -- a fractal flame renderer";
		const struct argp argp = {
				.options = options, .parser = parse_cross_opt,
				.args_doc = NULL, .doc = doc, .children = NULL
				};

		cross_arguments arguments = {
				.method = CROSS_NOT_SPECIFIED,
				};

		argp_parse (&argp, argc, argv, 0, NULL, &arguments);
		do_cross (&arguments);
	} else if (streq (command, "improvecolors")) {
		const struct argp_option options[] = {
				{ 0 },
				};
		const char doc[] = PACKAGE "-improvecolors -- improve flame colors";
		const struct argp argp = {
				.options = options, .parser = parse_improvecolors_opt,
				.args_doc = NULL, .doc = doc, .children = NULL
				};

		improvecolors_arguments arguments = {
				.tries = 100,
				.resolution = 10,
				.change_palette = true,
				.time = 1.0,
				};

		argp_parse (&argp, argc, argv, 0, NULL, &arguments);
		do_improvecolors (&arguments);
	} else if (streq (command, "interpolate")) {
		const struct argp_option options[] = {
				{"time", 't', "float", OPTION_ARG_OPTIONAL, "Time step (0.5)" },
				{ 0 },
				};
		const char doc[] = PACKAGE "-interpolate -- fractal flame interpolation";
		const struct argp argp = {
				.options = options, .parser = parse_interpolate_opt,
				.args_doc = NULL, .doc = doc, .children = NULL
				};

		interpolate_arguments arguments = {
				.time = 0.5,
				};

		argp_parse (&argp, argc, argv, 0, NULL, &arguments);
		do_interpolate (&arguments);
#if 0
	} else if (streq (command, "mutate")) {
		const struct argp_option options[] = {
				{"method", 'm', "XXX", OPTION_ARG_OPTIONAL, "Mutation method" },
				{ 0 },
				};
		const char doc[] = PACKAGE "-mutate -- a fractal flame renderer";
		const struct argp argp = {
				.options = options, .parser = parse_mutate_opt,
				.args_doc = NULL, .doc = doc, .children = NULL
				};

		mutate_arguments arguments = {
				.method = MUTATE_NOT_SPECIFIED,
				.symmetry = 0,
				};

		argp_parse (&argp, argc, argv, 0, NULL, &arguments);
		do_mutate (&arguments);
#endif
	} else if (streq (command, "random")) {
		/* generate random genome */
		const struct argp_option options[] = {
				{"height", 'h', "pixels", 0, "Output flame height (1000)" },
				{"width", 'w', "pixels", 0, "Output flame width (1000)" },
				{"max-xforms", 'x', "number", 0, "Max number of xforms (6)" },
				{"max-var", 'v', "number", 0, "Max number of variations per xform (unlimited)" },
				{ 0 },
				};
		const char doc[] = PACKAGE "-random -- a fractal flame generator";
		const struct argp argp = {
				.options = options, .parser = parse_random_opt,
				.args_doc = NULL, .doc = doc, .children = NULL
				};

		random_arguments arguments = {
				.palette = "flam3-palettes.xml",
				.width = 1000,
				.height = 1000,
				.max_xforms = 6,
				.max_var = flam3_nvariations,
				.post_likelihood = 0.4,
				.final_likelihood = 0.15,
				.symmetry_likelihood = 0.25,
				.min_symmetry = -6,
				.max_symmetry = 6,
				};

		argp_parse (&argp, argc, argv, 0, NULL, &arguments);
		do_random (&arguments);
	} else if (streq (command, "render")) {
		/* render flame to image file */
		const struct argp_option options[] = {
				{"scale", 's', "factor", 0, "Scale image dimensions by factor (1.0)" },
				{"bpc", 'b', "8|16", 0, "Bits per channel of output image (8)" },
				{"time", 't', "seconds", 0, "Rendering time" },
				{"cache", 'c', "path", 0, "Cache file" },
				{ 0 },
				};
		const char doc[] = PACKAGE "-render -- a fractal flame renderer";
		const struct argp argp = {
				.options = options, .parser = parse_render_opt,
				.args_doc = NULL, .doc = doc, .children = NULL,
				};

		render_arguments arguments = {
				.bpc = 8,
				.scale = 1.0,
				.time = 1.0,
				.cache = NULL,
				};

		argp_parse (&argp, argc, argv, 0, NULL, &arguments);
		do_render (&arguments);
	} else {
		show_help (argv[0]);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

