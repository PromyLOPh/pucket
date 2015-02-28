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
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <argp.h>

#include "random.h"
#include "private.h"
#include "img.h"

#define streq(a,b) (strcmp (a, b) == 0)

const char *argp_program_version =
  "vlam3-pre";

typedef struct {
	bool verbose;
	unsigned int threads, bpc, quality;
	float scale;
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

		case 'q': {
			int i = atoi (arg);
			if (i < 1) {
				argp_error (state, "Quality must be >= 1");
			} else {
				arguments->quality = i;
			}
			break;
		}

		case 's':
			arguments->scale = atof (arg);
			if (arguments->scale <= 0.0) {
				argp_error (state, "Scale must be > 0");
			}
			break;

		case 't': {
			int i = atoi (arg);
			if (i <= 0) {
				argp_error (state, "Threads must be > 0");
			} else {
				arguments->threads = i;
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

static void do_render (const render_arguments * const arguments) {
	randctx rc;

	rand_seed(&rc);

	int ncps;
	flam3_genome * const cps = flam3_parse_from_file (stdin, NULL,
			flam3_defaults_on, &ncps, &rc);
	if (cps == NULL) {
		fprintf(stderr,"error reading genomes from file\n");
		exit(1);
	}
	assert (ncps == 1);

	flam3_genome * const genome = &cps[0];

	/* Force ntemporal_samples to 1 for -render */
	genome->ntemporal_samples = 1;
	genome->sample_density = arguments->quality;
	genome->height *= arguments->scale;
	genome->width *= arguments->scale;
	genome->pixels_per_unit *= arguments->scale;

	flam3_frame f;
	f.genomes = genome;
	f.ngenomes = 1;
	f.verbose = arguments->verbose;
	f.time = 0.0;
	f.pixel_aspect_ratio = 1.0;
	f.progress = 0;
	f.nthreads = arguments->threads;
	f.earlyclip = 0;
	f.sub_batch_size = 10000;
	f.bytes_per_channel = arguments->bpc / 8;

	const unsigned int channels = 4;
	const size_t this_size = channels * genome->width * genome->height *
			f.bytes_per_channel;
	void *image = (void *) calloc(this_size, sizeof(char));

	stat_struct stats;
	if (flam3_render (&f, image, flam3_field_both, &stats)) {
		fprintf(stderr,"error rendering image: aborting.\n");
		exit(1);
	}

	flam3_img_comments fpc;
	write_png (stdout, image, genome->width, genome->height, &fpc,
			f.bytes_per_channel);
}

static void print_genome (flam3_genome * const genome) {
	printf("<pick version=\"FLAM3-%s\">\n", flam3_version());
	flam3_print (stdout, genome, NULL, flam3_dont_print_edits);
	printf("</pick>\n");
}

typedef struct {
	int symmetry;
	const char *palette;
	unsigned int width, height;
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
	genome->rot_center[0] = genome->center[0];
	genome->rot_center[1] = genome->center[1];
	genome->pixels_per_unit = genome->width / (bmax[0] - bmin[0]);
}

static void do_random (const random_arguments * const arguments) {
	randctx rc;
	rand_seed(&rc);

	flam3_genome genome = { .edits = NULL };
	int ivars = flam3_variation_random;
	flam3_random (&genome, &ivars, 1, arguments->symmetry, 0, &rc);

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
	flam3_genome * const cps = flam3_parse_from_file (stdin, NULL,
			flam3_defaults_on, &ncps, &rc);
	if (cps == NULL) {
		fprintf(stderr,"error reading genomes from file\n");
		exit(1);
	}
	assert (ncps == 1);

	flam3_genome * const genome = &cps[0];

	int ivars = flam3_variation_random;
	const double speed = 1.0;
	flam3_mutate (genome, arguments->method, &ivars, 1, arguments->symmetry,
			speed, &rc);

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
	flam3_genome * const cps = flam3_parse_from_file (stdin, NULL,
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

#if 0
static void do_improvecolors () {
	flam3_improve_colors(&cp_orig, 100, 0, 10, &rc);
}

static void do_interpolate () {
	for (ftime = first_frame; ftime <= last_frame; ftime += 1) {
		iscp=0;
		for (i=0;i<ncp;i++) {
			if ( ftime==cp[i].time ) {
				flam3_copy(&interpolated, &(cp[i]) );
				iscp=1;
			}
		}
		if (iscp==0) {
			flam3_interpolate(cp, ncp, (double)ftime, stagger, &interpolated);
			for (i=0;i<ncp;i++) {
				if ( ftime==cp[i].time-1 ) {
					iscp=1;
				}
			}
			if (iscp==0)
				interpolated.interpolation_type = flam3_inttype_linear;
		}

		if (templ) flam3_apply_template(&interpolated, templ);
		gprint(&interpolated, 1);
	}
}
#endif

static void show_help (const char * const argv0) {
	const char *progname = strrchr (argv0, (int) '/');
	if (progname == NULL) {
		progname = argv0;
	} else {
		++progname;
	}
	fprintf (stderr,
			"Usage: %s cross [OPTION...]\n"
			"   Or: %s mutate [OPTION...]\n"
			"   Or: %s random [OPTION...]\n"
			"   Or: %s render [OPTION...]\n",
			progname, progname, progname, progname);
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
		const char doc[] = "vlame3-cross -- a fractal flame renderer";
		const struct argp argp = {
				.options = options, .parser = parse_cross_opt,
				.args_doc = NULL, .doc = doc, .children = NULL
				};

		cross_arguments arguments = {
				.method = CROSS_NOT_SPECIFIED,
				};

		argp_parse (&argp, argc, argv, 0, NULL, &arguments);
		do_cross (&arguments);
	} else if (streq (command, "mutate")) {
		const struct argp_option options[] = {
				{"method", 'm', "XXX", OPTION_ARG_OPTIONAL, "Mutation method" },
				{ 0 },
				};
		const char doc[] = "vlame3-mutate -- a fractal flame renderer";
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
	} else if (streq (command, "random")) {
		/* generate random genome */
		const struct argp_option options[] = {
				{"height", 'h', "pixels", 0, "Output flame height" },
				{"width", 'w', "pixels", 0, "Output flame width" },
				{ 0 },
				};
		const char doc[] = "vlame3-random -- a fractal flame generator";
		const struct argp argp = {
				.options = options, .parser = parse_random_opt,
				.args_doc = NULL, .doc = doc, .children = NULL
				};

		random_arguments arguments = {
				.symmetry = 0,
				.palette = "flam3-palettes.xml",
				.width = 1000,
				.height = 1000,
				};

		argp_parse (&argp, argc, argv, 0, NULL, &arguments);
		do_random (&arguments);
	} else if (streq (command, "render")) {
		/* render flame to image file */
		const struct argp_option options[] = {
				{"threads", 't', "num", 0, "Number of threads (auto)" },
				{"scale", 's', "factor", 0, "Scale image dimensions by factor (1.0)" },
				{"bpc", 'b', "8|16", 0, "Bits per channel of output image (8)" },
				{"quality", 'q', "num", 0, "Average samples per pixel (100)" },
				{"width", 'w', "pixels", 0, "Output image width" },
				{"height", 'h', "pixels", 0, "Output image height" },
				{ 0 },
				};
		const char doc[] = "vlame3-render -- a fractal flame renderer";
		const struct argp argp = {
				.options = options, .parser = parse_render_opt,
				.args_doc = NULL, .doc = doc, .children = NULL,
				};

		render_arguments arguments = {
				.threads = flam3_count_nthreads(),
				.bpc = 8,
				.scale = 1.0,
				.quality = 100,
				.verbose = true,
				};

		argp_parse (&argp, argc, argv, 0, NULL, &arguments);
		do_render (&arguments);
	} else {
		show_help (argv[0]);
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

