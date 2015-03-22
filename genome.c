/*
    Copyright (C) 1992-2009 Spotworks LLC
	              2015 pucket contributors

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

#include <math.h>
#include <string.h>

#include "genome.h"
#include "math.h"
#include "flam3.h"
#include "random.h"

/* Generate random params for parametric variations */
static void random_xform_param (flam3_xform * const xform, const unsigned int i,
		randctx * const rc) {
	switch (i) {
		case VAR_BLOB:
			xform->blob_low = 0.2 + 0.5 * rand_d01(rc);
			xform->blob_high = 0.8 + 0.4 * rand_d01(rc);
			xform->blob_waves = (int)(2 + 5 * rand_d01(rc));
			break;

		case VAR_PDJ:
			xform->pdj_a = 3.0 * rand_d11(rc);
			xform->pdj_b = 3.0 * rand_d11(rc);
			xform->pdj_c = 3.0 * rand_d11(rc);
			xform->pdj_d = 3.0 * rand_d11(rc);
			break;

		case VAR_FAN2:
			xform->fan2_x = rand_d11(rc);
			xform->fan2_y = rand_d11(rc);
			break;

		case VAR_RINGS2:
			xform->rings2_val = 2*rand_d01(rc);
			break;

		case VAR_PERSPECTIVE:
			xform->perspective_angle = rand_d01(rc);
			xform->perspective_dist = 2*rand_d01(rc) + 1.0;
			break;

		case VAR_JULIAN:
			xform->julian_power = (int)(5*rand_d01(rc) + 2);
			xform->julian_dist = 1.0;
			break;

		case VAR_JULIASCOPE:
			xform->juliascope_power = (int)(5*rand_d01(rc) + 2);
			xform->juliascope_dist = 1.0;
			break;

		case VAR_RADIAL_BLUR:
			xform->radial_blur_angle = (2 * rand_d01(rc) - 1);
			break;

		case VAR_PIE:
			xform->pie_slices = (int) 10.0*rand_d01(rc);
			xform->pie_thickness = rand_d01(rc);
			xform->pie_rotation = 2.0 * M_PI * rand_d11(rc);
			break;

		case VAR_NGON:
			xform->ngon_sides = (int) rand_d01(rc)* 10 + 3;
			xform->ngon_power = 3*rand_d01(rc) + 1;
			xform->ngon_circle = 3*rand_d01(rc);
			xform->ngon_corners = 2*rand_d01(rc)*xform->ngon_circle;
			break;

		case VAR_CURL:
			xform->curl_c1 = rand_d01(rc);
			xform->curl_c2 = rand_d01(rc);
			break;

		case VAR_RECTANGLES:
			xform->rectangles_x = rand_d01(rc);
			xform->rectangles_y = rand_d01(rc);
			break;

		case VAR_DISC2:
			xform->disc2_rot = 0.5 * rand_d01(rc);
			xform->disc2_twist = 0.5 * rand_d01(rc);
			break;

		case VAR_SUPER_SHAPE:
			xform->super_shape_rnd = rand_d01(rc);
			xform->super_shape_m = (int) rand_d01(rc)*6;
			xform->super_shape_n1 = rand_d01(rc)*40;
			xform->super_shape_n2 = rand_d01(rc)*20;
			xform->super_shape_n3 = xform->super_shape_n2;
			xform->super_shape_holes = 0.0;
			break;

		case VAR_FLOWER:
			xform->flower_petals = 4 * rand_d01(rc);
			xform->flower_holes = rand_d01(rc);
			break;

		case VAR_CONIC:
			xform->conic_eccentricity = rand_d01(rc);
			xform->conic_holes = rand_d01(rc);
			break;

		case VAR_PARABOLA:
			xform->parabola_height = 0.5 + rand_d01(rc);
			xform->parabola_width = 0.5 + rand_d01(rc);
			break;

		case VAR_BENT2:
			xform->bent2_x = 3*(-0.5 + rand_d01(rc));
			xform->bent2_y = 3*(-0.5 + rand_d01(rc));
			break;

		case VAR_BIPOLAR:
			xform->bipolar_shift = 2.0 * rand_d01(rc) - 1;
			break;

		case VAR_CELL:
			xform->cell_size = 2.0 * rand_d01(rc) + 0.5;
			break;

		case VAR_CPOW:
			xform->cpow_r = 3.0 * rand_d01(rc);
			xform->cpow_i = rand_d01(rc) - 0.5;
			xform->cpow_power = (int)(5.0 * rand_d01(rc));
			break;

		case VAR_CURVE:
			xform->curve_xamp = 5 * (rand_d01(rc)-.5);
			xform->curve_yamp = 4 * (rand_d01(rc)-.5);
			xform->curve_xlength = 2 * (rand_d01(rc)+.5);
			xform->curve_ylength = 2 * (rand_d01(rc)+.5);
			break;

		case VAR_ESCHER:
			xform->escher_beta = M_PI * rand_d11(rc);
			break;

		case VAR_LAZYSUSAN:
			xform->lazysusan_x = 2.0*rand_d11(rc);
			xform->lazysusan_y = 2.0*rand_d11(rc);
			xform->lazysusan_spin = M_PI*rand_d11(rc);
			xform->lazysusan_space = 2.0*rand_d11(rc);
			xform->lazysusan_twist = 2.0*rand_d11(rc);
			break;

		case VAR_MODULUS:
			xform->modulus_x = rand_d11(rc);
			xform->modulus_y = rand_d11(rc);
			break;

		case VAR_OSCILLOSCOPE:
			xform->oscope_separation = 1.0 + rand_d11(rc);
			xform->oscope_frequency = M_PI * rand_d11(rc);
			xform->oscope_amplitude = 1.0 + 2 * rand_d01(rc);
			xform->oscope_damping = rand_d01(rc);
			break;

		case VAR_POPCORN2:
			xform->popcorn2_x = 0.2 * rand_d01(rc);
			xform->popcorn2_y = 0.2 * rand_d01(rc);
			xform->popcorn2_c = 5 * rand_d01(rc);
			break;

		case VAR_SEPARATION:
			xform->separation_x = 1 + rand_d11(rc);
			xform->separation_y = 1 + rand_d11(rc);
			xform->separation_xinside = rand_d11(rc);
			xform->separation_yinside = rand_d11(rc);
			break;

		case VAR_SPLIT:
			xform->split_xsize = rand_d11(rc);
			xform->split_ysize = rand_d11(rc);
			break;

		case VAR_SPLITS:
			xform->splits_x = rand_d11(rc);
			xform->splits_y = rand_d11(rc);
			break;

		case VAR_STRIPES:
			xform->stripes_space = rand_d01(rc);
			xform->stripes_warp = 5*rand_d01(rc);
			break;

		case VAR_WEDGE:
			xform->wedge_angle = M_PI*rand_d01(rc);
			xform->wedge_hole = 0.5*rand_d11(rc);
			xform->wedge_count = floor(5*rand_d01(rc))+1;
			xform->wedge_swirl = rand_d01(rc);
			break;

		case VAR_WEDGE_JULIA:
			xform->wedge_julia_power = (int)(5*rand_d01(rc) + 2);
			xform->wedge_julia_dist = 1.0;
			xform->wedge_julia_count = (int)(3*rand_d01(rc) + 1);
			xform->wedge_julia_angle = M_PI * rand_d01(rc);
			break;

		case VAR_WEDGE_SPH:
			xform->wedge_sph_angle = M_PI*rand_d01(rc);
			xform->wedge_sph_hole = 0.5*rand_d11(rc);
			xform->wedge_sph_count = floor(5*rand_d01(rc))+1;
			xform->wedge_sph_swirl = rand_d01(rc);
			break;

		case VAR_WHORL:
			xform->whorl_inside = rand_d01(rc);
			xform->whorl_outside = rand_d01(rc);
			break;

		case VAR_WAVES2:
			xform->waves2_scalex = 0.5 + rand_d01(rc);
			xform->waves2_scaley = 0.5 + rand_d01(rc);
			xform->waves2_freqx = 4 * rand_d01(rc);
			xform->waves2_freqy = 4 * rand_d01(rc);
			break;         

		case VAR_AUGER:
			xform->auger_sym = 0;
			xform->auger_weight = 0.5 + rand_d01(rc)/2.0;
			xform->auger_freq = floor(5*rand_d01(rc))+1;
			xform->auger_scale = rand_d01(rc);
			break;         

		case VAR_FLUX:
			xform->flux_spread = 0.5 + rand_d01(rc)/2.0;
			break;

		case VAR_MOBIUS:
			xform->mobius_re_a = rand_d11(rc);
			xform->mobius_im_a = rand_d11(rc);
			xform->mobius_re_b = rand_d11(rc);
			xform->mobius_im_b = rand_d11(rc);
			xform->mobius_re_c = rand_d11(rc);
			xform->mobius_im_c = rand_d11(rc);
			xform->mobius_re_d = rand_d11(rc);
			xform->mobius_im_d = rand_d11(rc);
			break;
	}
}

/*	Fill xform with random values
 */
void xform_rand (flam3_xform * const xform, const bool add_post,
		const unsigned int max_var, randctx * const rc) {
	assert (xform != NULL);
	assert (rc != NULL);

	/* XXX: the original code alternates between 0/1 for every xform */
	xform->color = rand_bool (rc) ? 1.0 : 0.0;

	for (unsigned int j = 0; j < 3; j++) {
		for (unsigned int k = 0; k < 2; k++) {
			xform->c[j][k] = rand_d11(rc);
			if (add_post) {
				xform->post[j][k] = rand_d11(rc);
			} else {
				xform->post[j][k] = (double)(k==j);
			}
		}
	}

	memset (xform->var, 0, sizeof (*xform->var));
	for (unsigned int i = 0; i < max_var; i++) {
		const unsigned int v = rand_mod (rc, flam3_nvariations);
		double w;
		do {
			w = rand_d01 (rc);
		} while (w == 0.0);
		xform->var[v] += w;
		random_xform_param (xform, v, rc);

		/* small number of variations is more likely */
		if (rand_bool (rc)) {
			break;
		}
	}

	/* Normalize weights to 1.0 total. */
	normalize (xform->var, flam3_nvariations);
}

