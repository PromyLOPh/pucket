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

#include <stdio.h>
#include <libxml/parser.h>

#define flam3_palette_random       (-1)
#define flam3_palette_interpolated (-2)

#define flam3_defaults_on          (1)
#define flam3_defaults_off         (0)

#define flam3_name_len    64

extern char *flam3_variation_names[];

#define flam3_nvariations 99

#define flam3_interpolation_linear 0
#define flam3_interpolation_smooth 1

#define flam3_inttype_linear 0
#define flam3_inttype_log    1
#define flam3_inttype_compat 2 /* Linear and old behaviour */
#define flam3_inttype_older  3 /* rotate padded xforms     */

typedef enum {
	PALETTE_MODE_STEP = 0,
	PALETTE_MODE_LINEAR = 1,
} color_palette_mode;

#define VAR_LINEAR   0
#define VAR_SINUSOIDAL  1
#define VAR_SPHERICAL  2
#define VAR_SWIRL 3
#define VAR_HORSESHOE  4
#define VAR_POLAR 5
#define VAR_HANDKERCHIEF 6
#define VAR_HEART 7
#define VAR_DISC 8
#define VAR_SPIRAL 9
#define VAR_HYPERBOLIC 10
#define VAR_DIAMOND 11
#define VAR_EX 12
#define VAR_JULIA 13
#define VAR_BENT 14
#define VAR_WAVES 15
#define VAR_FISHEYE 16
#define VAR_POPCORN 17
#define VAR_EXPONENTIAL 18
#define VAR_POWER 19
#define VAR_COSINE 20
#define VAR_RINGS 21
#define VAR_FAN 22
#define VAR_BLOB 23
#define VAR_PDJ 24
#define VAR_FAN2 25
#define VAR_RINGS2 26
#define VAR_EYEFISH 27
#define VAR_BUBBLE 28
#define VAR_CYLINDER 29
#define VAR_PERSPECTIVE 30
#define VAR_NOISE 31
#define VAR_JULIAN 32
#define VAR_JULIASCOPE 33
#define VAR_BLUR 34
#define VAR_GAUSSIAN_BLUR 35
#define VAR_RADIAL_BLUR 36
#define VAR_PIE 37
#define VAR_NGON 38
#define VAR_CURL 39
#define VAR_RECTANGLES 40
#define VAR_ARCH 41
#define VAR_TANGENT 42
#define VAR_SQUARE 43
#define VAR_RAYS 44
#define VAR_BLADE 45
#define VAR_SECANT2 46
#define VAR_TWINTRIAN 47
#define VAR_CROSS 48
#define VAR_DISC2 49
#define VAR_SUPER_SHAPE 50
#define VAR_FLOWER 51
#define VAR_CONIC 52
#define VAR_PARABOLA 53
#define VAR_BENT2 54
#define VAR_BIPOLAR 55
#define VAR_BOARDERS 56
#define VAR_BUTTERFLY 57
#define VAR_CELL 58
#define VAR_CPOW 59
#define VAR_CURVE 60
#define VAR_EDISC 61
#define VAR_ELLIPTIC 62
#define VAR_ESCHER 63
#define VAR_FOCI 64
#define VAR_LAZYSUSAN 65
#define VAR_LOONIE 66
#define VAR_PRE_BLUR 67
#define VAR_MODULUS 68
#define VAR_OSCILLOSCOPE 69
#define VAR_POLAR2 70
#define VAR_POPCORN2 71
#define VAR_SCRY 72
#define VAR_SEPARATION 73
#define VAR_SPLIT 74
#define VAR_SPLITS 75
#define VAR_STRIPES 76
#define VAR_WEDGE 77
#define VAR_WEDGE_JULIA 78
#define VAR_WEDGE_SPH 79
#define VAR_WHORL 80
#define VAR_WAVES2 81
#define VAR_EXP 82
#define VAR_LOG 83
#define VAR_SIN 84
#define VAR_COS 85
#define VAR_TAN 86
#define VAR_SEC 87
#define VAR_CSC 88
#define VAR_COT 89
#define VAR_SINH 90
#define VAR_COSH 91
#define VAR_TANH 92
#define VAR_SECH 93
#define VAR_CSCH 94
#define VAR_COTH 95
#define VAR_AUGER 96
#define VAR_FLUX 97
#define VAR_MOBIUS 98

#include "vector.h"
#include "random.h"

typedef struct xform {
   double var[flam3_nvariations];   /* interp coefs between variations */
   double2 c[3];      /* the coefs to the affine part of the function */
   double2 post[3];   /* the post transform */
   double density;      /* probability that this function is chosen. 0 - 1 */
   double color;     /* color coords for this function. 0 - 1 */
   double color_speed;  /* scaling factor on color added to current iteration */
   double animate;      /* whether or not this xform rotates (in sheep) >0 means stationary */
   double opacity;   /* 0=invisible, 1=totally visible */
   double vis_adjusted; /* adjusted visibility for better transitions */
   
   int padding;/* Set to 1 for padding xforms */
   double wind[2]; /* winding numbers */

   int precalc_angles_flag;
   int precalc_atan_xy_flag;
   int precalc_atan_yx_flag;
   double has_preblur;
   int has_post;

   /* Params for new parameterized variations */
   /* Blob */
   double blob_low;
   double blob_high;
   double blob_waves;

   /* PDJ */
   double pdj_a;
   double pdj_b;
   double pdj_c;
   double pdj_d;

   /* Fan2 */
   double fan2_x;
   double fan2_y;

   /* Rings2 */
   double rings2_val;

   /* Perspective */
   double perspective_angle;
   double perspective_dist;

   /* Julia_N */
   double julian_power;
   double julian_dist;

   /* Julia_Scope */
   double juliascope_power;
   double juliascope_dist;

   /* Radial_Blur */
   double radial_blur_angle;

   /* Pie */
   double pie_slices;
   double pie_rotation;
   double pie_thickness;

   /* Ngon */
   double ngon_sides;
   double ngon_power;
   double ngon_circle;
   double ngon_corners;

   /* Curl */
   double curl_c1;
   double curl_c2;

   /* Rectangles */
   double rectangles_x;
   double rectangles_y;

   /* Disc 2 */
   double disc2_rot;
   double disc2_twist;

   /* Supershape */
   double super_shape_rnd;
   double super_shape_m;
   double super_shape_n1;
   double super_shape_n2;
   double super_shape_n3;
   double super_shape_holes;
   
   /* Flower */
   double flower_petals;
   double flower_holes;
   
   /* Conic */
   double conic_eccentricity;
   double conic_holes;
   
   /* Parabola */
   double parabola_height;
   double parabola_width;
   
   /* Bent2 */
   double bent2_x;
   double bent2_y;
   
   /* Bipolar */
   double bipolar_shift;
   
   /* Cell */
   double cell_size;
   
   /* Cpow */
   double cpow_r;
   double cpow_i;
   double cpow_power; /* int in apo */
   
   /* Curve */
   double curve_xamp,curve_yamp;
   double curve_xlength,curve_ylength;
   
   /* Escher */
   double escher_beta;
   
   /* Lazysusan */
   double lazysusan_spin;
   double lazysusan_space;
   double lazysusan_twist;
   double lazysusan_x, lazysusan_y;
   
   /* Modulus */
   double modulus_x, modulus_y;
   
   /* Oscope */
   double oscope_separation;
   double oscope_frequency;
   double oscope_amplitude;
   double oscope_damping;
   
   /* Popcorn2 */
   double popcorn2_x, popcorn2_y, popcorn2_c;
   
   /* Separation */
   double separation_x, separation_xinside;
   double separation_y, separation_yinside;
   
   /* Split */
   double split_xsize;
   double split_ysize;
   
   /* Splits */
   double splits_x,splits_y;
   
   /* Stripes */
   double stripes_space;
   double stripes_warp;
   
   /* Wedge */
   double wedge_angle, wedge_hole;
   double wedge_count, wedge_swirl;
   
   /* Wedge_Julia */
   double wedge_julia_angle;
   double wedge_julia_count;
   double wedge_julia_power;
   double wedge_julia_dist;
   
   /* Wedge_Sph */
   double wedge_sph_angle, wedge_sph_count;
   double wedge_sph_hole, wedge_sph_swirl;
   
   /* Whorl */
   double whorl_inside, whorl_outside;
   
   /* Waves2 */
   double waves2_freqx, waves2_scalex;
   double waves2_freqy, waves2_scaley;
   
   /* Auger */
   double auger_sym, auger_weight;
   double auger_freq, auger_scale;

   /* Flux */
   double flux_spread;

   /* Mobius */
   double mobius_re_a, mobius_im_a;
   double mobius_re_b, mobius_im_b;
   double mobius_re_c, mobius_im_c;
   double mobius_re_d, mobius_im_d;
      
   /* If perspective is used, precalculate these values */
   /* from the _angle and _dist                         */
   double persp_vsin;
   double persp_vfcos;

   /* If Julia_N is used, precalculate these values */
   double julian_rN;
   double julian_cn;

   /* If Julia_Scope is used, precalculate these values */
   double juliascope_rN;
   double juliascope_cn;
   
   /* if Wedge_Julia, precalculate */
   double wedgeJulia_rN;
   double wedgeJulia_cn;
   double wedgeJulia_cf;

   /* If Radial_Blur is used, precalculate these values */
   double radialBlur_spinvar;
   double radialBlur_zoomvar;

   /* Precalculate these values for waves */
   double2 waves_d2;

   /* If disc2 is used, precalculate these values */
   double disc2_sinadd;
   double disc2_cosadd;
   double disc2_timespi;

   /* If supershape is used, precalculate these values */
   double super_shape_pm_4;
   double super_shape_pneg1_n1;

   int num_active_vars;
   double active_var_weights[flam3_nvariations];
   int varFunc[flam3_nvariations];
} flam3_xform;

#include "palettes.h"

typedef struct {
   char flame_name[flam3_name_len+1]; /* 64 chars plus a null */
   double time;
   int interpolation;
   int interpolation_type;
   int num_xforms;
   int final_xform_index;
   int final_xform_enable;
   flam3_xform *xform;
   
   /* Xaos implementation */
   double **chaos;
   int chaos_enable;
   
   int genome_index;                   /* index into source file */
   int symmetry;                /* 0 means none */
   palette palette;
   char *input_image;           /* preview/temporary! */
   double brightness;           /* 1.0 = normal */
   double contrast;             /* 1.0 = normal */
   double gamma;
   double highlight_power;
   int  width, height;          /* of the final image */
   double center[2];             /* of camera */
   double rotate;                /* camera */
   double vibrancy;              /* blend between color algs (0=old,1=new) */
   double hue_rotation;          /* applies to cmap, 0-1 */
   double zoom;                  /* effects ppu, sample density, scale */
   double pixels_per_unit;       /* vertically */

   /* Small-gamma linearization threshold */
   double gam_lin_thresh;

   /* for cmap_interpolated hack */
   int palette_index0;
   double hue_rotation0;
   int palette_index1;
   double hue_rotation1;
   double palette_blend;

   int palette_mode;
} flam3_genome;

typedef struct {
	int from;
	int to;
	double scalar;
} flam3_chaos_entry;

/* xform manipulation */

void flam3_add_xforms(flam3_genome *cp, int num_to_add, int interp_padding, int final_flag);
void flam3_delete_xform(flam3_genome *thiscp, int idx_to_delete);
void flam3_copy(flam3_genome *dest, const flam3_genome * const src);
void flam3_copyx(flam3_genome * restrict dest, flam3_genome * restrict src,
		int dest_std_xforms, int dest_final_xform);
void flam3_copy_params(flam3_xform * restrict dest, flam3_xform * restrict src, int varn);

unsigned short* flam3_create_xform_distrib(flam3_genome *cp);
void clear_cp(flam3_genome *cp, int def_flag);

/* samples is array nsamples*4 long of x,y,color triples.
   using (samples[0], samples[1]) as starting XY point and
   (samples[2], samples[3]) as starting color coordinate,
   perform fuse iterations and throw them away, then perform
   nsamples iterations and save them in the samples array */
int flam3_iterate(flam3_genome *cp, int n, int fuse, const double4 in, double4 *samples, const unsigned short *xform_distrib, randctx *rc);

/* genomes is array ngenomes long, with times set and in ascending order.
   interpolate to the requested time and return in result */
void flam3_interpolate(flam3_genome *genomes, int ngenomes, double time, double stagger, flam3_genome *result);

/* print genome to given file with extra_attributes if not NULL */
void flam3_print(FILE *f, flam3_genome *g, char *extra_attributes);

void flam3_random(flam3_genome *cp, const unsigned int max_xform,
		const palette_collection * const pc, randctx * const rc);

void flam3_mutate(flam3_genome *cp, int mutate_mode, int *ivars, int ivars_n, int sym, double speed, const palette_collection * const pc, randctx *rc);
void flam3_cross(flam3_genome *cp0, flam3_genome *cp1, flam3_genome *out, int cross_mode, randctx *rc);

/* return NULL in case of error */
flam3_genome *flam3_parse_xml2(const int, int default_flag, int *ncps, randctx * const);

void flam3_add_symmetry(flam3_genome *cp, int sym, randctx * const rc);

void flam3_improve_colors(flam3_genome *g, int ntries, int change_palette, int color_resolution, const palette_collection * const pc, randctx * const rc);
int flam3_estimate_bounding_box(flam3_genome *g, double eps, int nsamples,
             double *bmin, double *bmax, randctx *rc);

typedef struct {
   double         pixel_aspect_ratio;    /* width over height of each pixel */
   flam3_genome  *genomes;
   int            ngenomes;
   int            bytes_per_channel;
   int            earlyclip;
   double         time;
   int            (*progress)(void *, double, int, double);
   void          *progress_parameter;
   int           nthreads;
   int           sub_batch_size;
} flam3_frame;

/* Mutation modes */
#define MUTATE_NOT_SPECIFIED   -1
#define MUTATE_ALL_VARIATIONS  0
#define MUTATE_ONE_XFORM_COEFS 1
#define MUTATE_ADD_SYMMETRY    2
#define MUTATE_POST_XFORMS     3
#define MUTATE_COLOR_PALETTE   4
#define MUTATE_DELETE_XFORM    5
#define MUTATE_ALL_COEFS       6

/* Cross modes */
#define CROSS_NOT_SPECIFIED   -1
#define CROSS_UNION           0
#define CROSS_INTERPOLATE     1  
#define CROSS_ALTERNATE       2

typedef struct {
	unsigned int consec, lastxf;
	double4 p;
	const unsigned short *xform_distrib;
	const flam3_genome *genome;
} iterator;

void iterator_init (iterator * const iter, const flam3_genome * const genome,
		const unsigned short * const xform_distrib, randctx * const rc);
bool iterator_step (iterator * const iter, double4 * const ret, randctx * const rc);
