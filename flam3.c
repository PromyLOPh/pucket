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

#include "rect.h"
#include "img.h"
#include "build/config.h"
#include "variations.h"
#include "interpolation.h"
#include "parser.h"
#include "palettes.h"
#include "random.h"
#include "math.h"
#include <string.h>
#include <locale.h>
#include <assert.h>

#define CHOOSE_XFORM_GRAIN 16384
#define CHOOSE_XFORM_GRAIN_M1 16383

static void flam3_print_xform(FILE *f, flam3_xform *x, int final_flag, int numstd, double *chaos_row);

unsigned short *flam3_create_xform_distrib(flam3_genome *cp) {

   /* Xform distrib is created in this function             */   
   int numrows;
   int dist_row,i;
   unsigned short *xform_distrib;
   
   numrows = cp->num_xforms - (cp->final_xform_index>=0) + 1;
   xform_distrib = calloc(numrows*CHOOSE_XFORM_GRAIN,sizeof(unsigned short));
   
   /* First, set up the first row of the xform_distrib (raw weights) */
   flam3_create_chaos_distrib(cp, -1, xform_distrib);
   
   /* Check for non-unity chaos */
   cp->chaos_enable = 1 - flam3_check_unity_chaos(cp);
   
   if (cp->chaos_enable) {
   
      /* Now set up a row for each of the xforms */
      dist_row = 0;
      for (i=0;i<cp->num_xforms;i++) {
      
         if (cp->final_xform_index == i)
            continue;
         else
            dist_row++;
         
         if (flam3_create_chaos_distrib(cp, i, &(xform_distrib[CHOOSE_XFORM_GRAIN*(dist_row)]))) {
            free(xform_distrib);
            return(NULL);
         }
      }
   }
   
   return(xform_distrib);
}

int flam3_check_unity_chaos(flam3_genome *cp) {

   int i,j;
   int num_std;
   int unity=1;
   num_std = cp->num_xforms - (cp->final_xform_index >= 0);
   
   for (i=0;i<num_std;i++) {
      for (j=0;j<num_std;j++) {
         if ( fabs(cp->chaos[i][j]-1.0) > EPS)
            unity=0;
      }
   }
   
   return(unity);
}

int flam3_create_chaos_distrib(flam3_genome *cp, int xi, unsigned short *xform_distrib) {

   /* Xform distrib is a preallocated array of CHOOSE_XFORM_GRAIN chars */
   /* address of array is passed in, contents are modified              */
   double t,r,dr;
   int i,j;
   int num_std;
   
   //fprintf(stdout,"storing at %ld\n",xform_distrib);
   
   num_std = cp->num_xforms - (cp->final_xform_index >= 0);
   
   dr = 0.0;
   for (i = 0; i < num_std; i++) {
      double d = cp->xform[i].density;
                     
      if (xi>=0)
         d *= cp->chaos[xi][i];
         
      //fprintf(stdout,"%f ",d);
      if (d < 0.0) {
         fprintf(stderr, "xform weight must be non-negative, not %g.\n", d);
         return(1);
      }         
      
      dr += d;
   }
   
   //fprintf(stdout,"dr=%f\n",dr);
   
   if (dr == 0.0) {
      fprintf(stderr, "cannot iterate empty flame.\n");
      return(1);
   }
   
   dr = dr / CHOOSE_XFORM_GRAIN;

   j = 0;
   t = cp->xform[0].density;
   if (xi>=0)
     t *= cp->chaos[xi][0];
   r = 0.0;
   for (i = 0; i < CHOOSE_XFORM_GRAIN; i++) {
      while (r >= t) {
         j++;

         if (xi>=0)
            t += cp->xform[j].density*cp->chaos[xi][j];
         else
            t += cp->xform[j].density;
         
      }
      //fprintf(stdout,"%d ",j);
      xform_distrib[i] = j;
      r += dr;
   }
   //fprintf(stdout,"\n---\n");
   
   return(0);
}

void iterator_init (iterator * const iter, const flam3_genome * const genome,
		const unsigned short * const xform_distrib, randctx * const rc) {
	iter->consec = 0;
	iter->lastxf = 0;
	iter->genome = genome;
	iter->xform_distrib = xform_distrib;
	iter->p = (double4) {
			rand_d11(rc),
			rand_d11(rc),
			rand_d01(rc),
			rand_d01(rc),
			};
}

/*	xform_precalc must be called once for each xform before running this
 *	function
 */
bool iterator_step (iterator * const iter, double4 * const ret, randctx * const rc) {
	const flam3_genome * const genome = iter->genome;
	unsigned int fn;
	double4 q;

	if (genome->chaos_enable)
		fn = iter->xform_distrib[ iter->lastxf*CHOOSE_XFORM_GRAIN + (rand_u64(rc) & CHOOSE_XFORM_GRAIN_M1)];
	else
		fn = iter->xform_distrib[ rand_u64(rc) & CHOOSE_XFORM_GRAIN_M1 ];

	if (apply_xform(genome, fn, iter->p, &q, rc)>0) {
		++iter->consec;
		if (iter->consec < 5) {
			iter->p = q;
			return false;
		} else
			iter->consec = 0;
	} else
		iter->consec = 0;

	/* Store the last used transform */
	iter->lastxf = fn+1;

	iter->p = q;

	if (genome->final_xform_enable == 1) {
		if (genome->xform[genome->final_xform_index].opacity==1 || 
				rand_d01(rc)<genome->xform[genome->final_xform_index].opacity) {
			apply_xform(genome, genome->final_xform_index, iter->p, &q, rc);
			/* Keep the opacity from the original xform */
			q = (double4) { q[0], q[1], q[2], iter->p[3] };
		}
	}

	*ret = q;

	return true;
}

/*
 * create a control point that interpolates between the control points
 * passed in CPS.  CPS must be sorted by time.
 */
void flam3_interpolate(flam3_genome cps[], int ncps,
             double time, double stagger, flam3_genome *result) {
   int i1, i2;
   double c[2];
   flam3_genome cpi[4];
   int smoothflag = 0;

   if (1 == ncps) {
      flam3_copy(result, &(cps[0]));
      return;
   }
      
   if (cps[0].time >= time) {
      i1 = 0;
      i2 = 1;
   } else if (cps[ncps - 1].time <= time) {
      i1 = ncps - 2;
      i2 = ncps - 1;
   } else {
      i1 = 0;
      while (cps[i1].time < time)
         i1++;

      i1--;
      i2 = i1 + 1;

   }

   c[0] = (cps[i2].time - time) / (cps[i2].time - cps[i1].time);
   c[1] = 1.0 - c[0];

   memset(cpi, 0, 4*sizeof(flam3_genome));
   
   /* To interpolate the xforms, we will make copies of the source cps  */
   /* and ensure that they both have the same number before progressing */
   if (flam3_interpolation_linear == cps[i1].interpolation) {
       flam3_align(&cpi[0], &cps[i1], 2);
       smoothflag = 0;
     
   } else {
       if (0 == i1) {
          fprintf(stderr, "error: cannot use smooth interpolation on first segment.\n");
          fprintf(stderr, "reverting to linear interpolation.\n");
          flam3_align(&cpi[0], &cps[i1], 2);
          smoothflag = 0;
       }

       if (ncps-1 == i2) {
          fprintf(stderr, "error: cannot use smooth interpolation on last segment.\n");
          fprintf(stderr, "reverting to linear interpolation.\n");
          flam3_align(&cpi[0], &cps[i1], 2);
          smoothflag = 0;
       }

       flam3_align(&cpi[0], &cps[i1-1], 4);
       smoothflag = 1;
   }
   
   /* Clear the destination cp */
   clear_cp(result, 1);
      
   if (cpi[0].final_xform_index >= 0) {
      flam3_add_xforms(result, cpi[0].num_xforms-1, 0, 0);
      flam3_add_xforms(result, 1, 0, 1);
   } else
      flam3_add_xforms(result, cpi[0].num_xforms, 0, 0);
   

   result->time = time;
   result->interpolation = flam3_interpolation_linear;
   result->interpolation_type = cpi[0].interpolation_type;

   if (!smoothflag) {
       flam3_interpolate_n(result, 2, cpi, c, stagger);
   } else {
       interpolate_catmull_rom(cpi, c[1], result);
       clear_cp(&(cpi[2]),0);
       clear_cp(&(cpi[3]),0);
   }
   
   clear_cp(&(cpi[0]),0);
   clear_cp(&(cpi[1]),0);

}

void flam3_copy_params(flam3_xform * restrict dest, flam3_xform * restrict src, int varn) {

   /* We only want to copy param var coefs for this one */
   if (varn==VAR_BLOB) {
      /* Blob */
      dest->blob_low = src->blob_low;
      dest->blob_high = src->blob_high;
      dest->blob_waves = src->blob_waves;
   } else if (varn==VAR_PDJ) {
      /* PDJ */
      dest->pdj_a = src->pdj_a;
      dest->pdj_b = src->pdj_b;
      dest->pdj_c = src->pdj_c;
      dest->pdj_d = src->pdj_d;
   } else if (varn==VAR_FAN2) {
      /* Fan2 */
      dest->fan2_x = src->fan2_x;
      dest->fan2_y = src->fan2_y;
   } else if (varn==VAR_RINGS2) {
      /* Rings2 */
      dest->rings2_val = src->rings2_val;
   } else if (varn==VAR_PERSPECTIVE) {
      /* Perspective */
      dest->perspective_angle = src->perspective_angle;
      dest->perspective_dist = src->perspective_dist;
      dest->persp_vsin = src->persp_vsin;
      dest->persp_vfcos = src->persp_vfcos;
   } else if (varn==VAR_JULIAN) {
      /* Julia_N */
      dest->julian_power = src->julian_power;
      dest->julian_dist = src->julian_dist;
      dest->julian_rN = src->julian_rN;
      dest->julian_cn = src->julian_cn;
   } else if (varn==VAR_JULIASCOPE) {
      /* Julia_Scope */
      dest->juliascope_power = src->juliascope_power;
      dest->juliascope_dist = src->juliascope_dist;
      dest->juliascope_rN = src->juliascope_rN;
      dest->juliascope_cn = src->juliascope_cn;
   } else if (varn==VAR_RADIAL_BLUR) {
      /* Radial Blur */
      dest->radial_blur_angle = src->radial_blur_angle;
   } else if (varn==VAR_PIE) {
      /* Pie */
      dest->pie_slices = src->pie_slices;
      dest->pie_rotation = src->pie_rotation;
      dest->pie_thickness = src->pie_thickness;
   } else if (varn==VAR_NGON) {
      /* Ngon */
      dest->ngon_sides = src->ngon_sides;
      dest->ngon_power = src->ngon_power;
      dest->ngon_corners = src->ngon_corners;
      dest->ngon_circle = src->ngon_circle;
   } else if (varn==VAR_CURL) {
      /* Curl */
      dest->curl_c1 = src->curl_c1;
      dest->curl_c2 = src->curl_c2;
   } else if (varn==VAR_RECTANGLES) {
      /* Rect */
      dest->rectangles_x = src->rectangles_x;
      dest->rectangles_y = src->rectangles_y;
   } else if (varn==VAR_DISC2) {
      /* Disc2 */
      dest->disc2_rot = src->disc2_rot;
      dest->disc2_twist = src->disc2_twist;
   } else if (varn==VAR_SUPER_SHAPE) {
      /* Supershape */
      dest->super_shape_rnd = src->super_shape_rnd;
      dest->super_shape_m = src->super_shape_m;
      dest->super_shape_n1 = src->super_shape_n1;
      dest->super_shape_n2 = src->super_shape_n2;
      dest->super_shape_n3 = src->super_shape_n3;
      dest->super_shape_holes = src->super_shape_holes;
   } else if (varn==VAR_FLOWER) {
      /* Flower */
      dest->flower_petals = src->flower_petals;
      dest->flower_petals = src->flower_petals;
   } else if (varn==VAR_CONIC) {
      /* Conic */
      dest->conic_eccentricity = src->conic_eccentricity;
      dest->conic_holes = src->conic_holes;
   } else if (varn==VAR_PARABOLA) {
      /* Parabola */
      dest->parabola_height = src->parabola_height;
      dest->parabola_width = src->parabola_width;
   } else if (varn==VAR_BENT2) {
      /* Bent2 */
      dest->bent2_x = src->bent2_x;
      dest->bent2_y = src->bent2_y;
   } else if (varn==VAR_BIPOLAR) {
      /* Bipolar */
      dest->bipolar_shift = src->bipolar_shift;
   } else if (varn==VAR_CELL) {
      /* Cell */
      dest->cell_size = src->cell_size;
   } else if (varn==VAR_CPOW) {
      /* Cpow */
      dest->cpow_i = src->cpow_i;
      dest->cpow_r = src->cpow_r;
      dest->cpow_power = src->cpow_power;
   } else if (varn==VAR_CURVE) {
      /* Curve */
      dest->curve_xamp = src->curve_xamp;
      dest->curve_yamp = src->curve_yamp;
      dest->curve_xlength = src->curve_xlength;
      dest->curve_ylength = src->curve_ylength;
   } else if (varn==VAR_ESCHER) {
      /* Escher */
      dest->escher_beta = src->escher_beta;
   } else if (varn==VAR_LAZYSUSAN) {
      /* Lazysusan */
      dest->lazysusan_x = src->lazysusan_x;
      dest->lazysusan_y = src->lazysusan_y;
      dest->lazysusan_spin = src->lazysusan_spin;
      dest->lazysusan_space = src->lazysusan_space;
      dest->lazysusan_twist = src->lazysusan_twist;
   } else if (varn==VAR_MODULUS) {
      /* Modulus */
      dest->modulus_x = src->modulus_x;
      dest->modulus_y = src->modulus_y;
   } else if (varn==VAR_OSCILLOSCOPE) {
      /* Oscope */
      dest->oscope_separation = src->oscope_separation;
      dest->oscope_frequency = src->oscope_frequency;
      dest->oscope_amplitude = src->oscope_amplitude;
      dest->oscope_damping = src->oscope_damping;
   } else if (varn==VAR_POPCORN2) {
      /* Popcorn2 */
      dest->popcorn2_x = src->popcorn2_x;
      dest->popcorn2_y = src->popcorn2_y;
      dest->popcorn2_c = src->popcorn2_c;
   } else if (varn==VAR_SEPARATION) {
      /* Separation */
      dest->separation_x = src->separation_x;
      dest->separation_y = src->separation_y;
      dest->separation_xinside = src->separation_xinside;
      dest->separation_yinside = src->separation_yinside;
   } else if (varn==VAR_SPLIT) {
      /* Split */
      dest->split_xsize = src->split_xsize;
      dest->split_ysize = src->split_ysize;
   } else if (varn==VAR_SPLITS) {
      /* Splits */
      dest->splits_x = src->splits_x;
      dest->splits_y = src->splits_y;
   } else if (varn==VAR_STRIPES) {
      /* Stripes */
      dest->stripes_space = src->stripes_space;
      dest->stripes_warp = src->stripes_warp;
   } else if (varn==VAR_WEDGE) {
      /* Wedge */
      dest->wedge_angle = src->wedge_angle;
      dest->wedge_hole = src->wedge_hole;
      dest->wedge_count = src->wedge_count;
      dest->wedge_swirl = src->wedge_swirl;
   } else if (varn==VAR_WEDGE_JULIA) {
      /* Wedge_Julia */
      dest->wedge_julia_angle = src->wedge_julia_angle;
      dest->wedge_julia_count = src->wedge_julia_count;
      dest->wedge_julia_power = src->wedge_julia_power;
      dest->wedge_julia_dist = src->wedge_julia_dist;
      dest->wedgeJulia_cf = src->wedgeJulia_cf;
      dest->wedgeJulia_cn = src->wedgeJulia_cn;
      dest->wedgeJulia_rN = src->wedgeJulia_rN;
   } else if (varn==VAR_WEDGE_SPH) {
      /* Wedge_sph */
      dest->wedge_sph_angle = src->wedge_sph_angle;
      dest->wedge_sph_hole = src->wedge_sph_hole;
      dest->wedge_sph_count = src->wedge_sph_count;
      dest->wedge_sph_swirl = src->wedge_sph_swirl;      
   } else if (varn==VAR_WHORL) {
      /* whorl */
      dest->whorl_inside = src->whorl_inside;
      dest->whorl_outside = src->whorl_outside;
   } else if (varn==VAR_WAVES2) {
      /* waves2 */
      dest->waves2_scalex = src->waves2_scalex;
      dest->waves2_scaley = src->waves2_scaley;
      dest->waves2_freqx = src->waves2_freqx;
      dest->waves2_freqy = src->waves2_freqy;
   } else if (varn==VAR_AUGER) {
      /* auger */
      dest->auger_sym = src->auger_sym;
      dest->auger_weight = src->auger_weight;
      dest->auger_freq = src->auger_freq;
      dest->auger_scale = src->auger_scale;
   } else if (varn==VAR_FLUX) {
      /* flux */
      dest->flux_spread = src->flux_spread;
   } else if (varn==VAR_MOBIUS) {
      /* mobius */
      dest->mobius_re_a = src->mobius_re_a;
      dest->mobius_re_b = src->mobius_re_b;
      dest->mobius_re_c = src->mobius_re_c;
      dest->mobius_re_d = src->mobius_re_d;
      dest->mobius_im_a = src->mobius_im_a;
      dest->mobius_im_b = src->mobius_im_b;
      dest->mobius_im_c = src->mobius_im_c;
      dest->mobius_im_d = src->mobius_im_d;
   }
}

/* Xform support functions */
void flam3_add_xforms(flam3_genome *thiscp, int num_to_add, int interp_padding, int final_flag) {

   int i,j;
   int old_num = thiscp->num_xforms;
   int oldstd,numstd;
   flam3_xform tmp;
   
   oldstd = thiscp->num_xforms - (thiscp->final_xform_index >= 0);
   
   /* !!! must make sure that if final_flag is specified, we don't already have a final xform! !!! */

//   if (thiscp->num_xforms > 0)
      thiscp->xform = (flam3_xform *)realloc(thiscp->xform, (thiscp->num_xforms + num_to_add) * sizeof(flam3_xform));
//   else
//      thiscp->xform = (flam3_xform *)malloc(num_to_add * sizeof(flam3_xform));

   thiscp->num_xforms += num_to_add;

   /* Initialize all the new xforms */
   initialize_xforms(thiscp, old_num);

   /* Set the padding flag for the new xforms */
   if (interp_padding) {
      for (i = old_num ; i < thiscp->num_xforms ; i++)
         thiscp->xform[i].padding=1;
   }
   
   /* If the final xform is not the last xform in the list, make it so */
   if (thiscp->final_xform_index >= 0 && thiscp->final_xform_index != thiscp->num_xforms-1) {
      tmp = thiscp->xform[thiscp->final_xform_index];
      for (i=thiscp->final_xform_index; i < thiscp->num_xforms-1; i++)
         thiscp->xform[i] = thiscp->xform[i+1];
      
      thiscp->final_xform_index = thiscp->num_xforms-1;
      thiscp->xform[thiscp->final_xform_index] = tmp;
   }
   
   if (final_flag) {
      /* Set the final xform index */
      thiscp->final_xform_enable = 1;
      thiscp->final_xform_index = thiscp->num_xforms-1;
   } else {
      /* Handle the chaos array */
      numstd = thiscp->num_xforms - (thiscp->final_xform_index>=0);
      
      /* Pad existing rows */
      for (i=0;i<oldstd;i++) {
         thiscp->chaos[i] = realloc(thiscp->chaos[i], numstd * sizeof(double));
         for (j=oldstd; j<numstd; j++)
            thiscp->chaos[i][j] = 1.0;
      }
      
      /* Add new rows */
      thiscp->chaos = realloc(thiscp->chaos,numstd * sizeof(double *));
      for (i=oldstd; i<numstd; i++) {
         thiscp->chaos[i] = malloc(numstd * sizeof(double));
         for (j=0;j<numstd;j++)
            thiscp->chaos[i][j] = 1.0;
      }
   }
}

void flam3_delete_xform(flam3_genome *thiscp, int idx_to_delete) {

   int i,j;
   int num_std = thiscp->num_xforms - (thiscp->final_xform_index >= 0);

   if (thiscp->final_xform_index != idx_to_delete) {
      /* We're going to delete the nth std xform. */

      /* Delete the nth_std row of the chaos array */
      free(thiscp->chaos[idx_to_delete]);

      /* Shift the pointers down one */
      for (i=idx_to_delete+1;i<num_std;i++)
         thiscp->chaos[i-1] = thiscp->chaos[i];
      
      /* Realloc the pointer array */
      thiscp->chaos = realloc(thiscp->chaos,(num_std-1)*sizeof(double *));
      num_std--;
      
      /* Loop over all of the rows and remove the nth_std element from them */
      for (i=0;i<num_std;i++) {
         for (j=idx_to_delete+1;j<num_std+1;j++) {
            thiscp->chaos[i][j-1] = thiscp->chaos[i][j];
         }
         /* Realloc the vector to have one less element */
         thiscp->chaos[i] = realloc(thiscp->chaos[i],num_std*sizeof(double));
         
      }
   }      
   
   /* Handle the final xform index */
   if (thiscp->final_xform_index == idx_to_delete) {
      thiscp->final_xform_index = -1;
      thiscp->final_xform_enable = 0;
   } else if (thiscp->final_xform_index > idx_to_delete) {
      thiscp->final_xform_index--;
   }

   /* Move all of the xforms down one - this does not require manual motion xform adjustment */
   for (i=idx_to_delete; i<thiscp->num_xforms-1; i++)
      thiscp->xform[i] = thiscp->xform[i+1];

   thiscp->num_xforms--;

   /* Reduce the memory storage by one xform */
   thiscp->xform = (flam3_xform *)realloc(thiscp->xform, sizeof(flam3_xform) * thiscp->num_xforms);
   
}

void flam3_copy_xform(flam3_xform *dest, flam3_xform *src) {
   memcpy (dest, src, sizeof (*dest));
}

/* Copy one control point to another */
void flam3_copy(flam3_genome *dest, const flam3_genome * const src) {
   
   int i,ii;
   int numstd;

   /* If there are any xforms in dest before the copy, clean them up */
   clear_cp(dest, 1);

   /* Copy main contents of genome */
   memcpy(dest, src, sizeof(flam3_genome));

   /* Only the pointer to the xform was copied, not the actual xforms. */
   /* We need to create new xform memory storage for this new cp       */
   /* This goes for chaos, too.                                        */
   dest->num_xforms = 0;
   dest->final_xform_index = -1;
   dest->xform = NULL;
   dest->chaos = NULL;

   /* Add the standard xforms first */
   numstd = src->num_xforms-(src->final_xform_index>=0);
   flam3_add_xforms(dest, numstd, 0, 0);
   for (i=0;i<numstd;i++)
      flam3_copy_xform(&dest->xform[i], &src->xform[i]);
      
   /* Add the final x if it's present */
   if (src->final_xform_index>=0) {
      i = src->final_xform_index;
      flam3_add_xforms(dest, 1, 0, 1);
      ii = dest->final_xform_index;
      flam3_copy_xform(&dest->xform[ii],&src->xform[i]);      
   }
   
   /* Also, only the pointer to the chaos array was copied.
    * We have to take care of that as well.                 */   
   for (i=0;i<numstd;i++)
      memcpy(dest->chaos[i],src->chaos[i], numstd * sizeof(double));
      
   palette_copy (&src->palette, &dest->palette);
}

void flam3_copyx(flam3_genome * restrict dest, flam3_genome * restrict src,
		int dest_std_xforms, int dest_final_xform) {

   int i,numsrcstd;
   
   /* If there are any xforms in dest before the copy, clean them up */
   clear_cp(dest, 1);

   /* Copy main contents of genome */
   memcpy(dest, src, sizeof(flam3_genome));

   /* Only the pointer to the xform was copied, not the actual xforms. */
   /* We need to create new xform memory storage for this new cp       */
   /* This goes for chaos, too.                                        */
   dest->num_xforms = 0;
   dest->xform = NULL;
   dest->chaos = NULL;
   dest->final_xform_index = -1;

   /* Add the padded standard xform list */
   /* Set the pad to 1 for these */
   flam3_add_xforms(dest, dest_std_xforms, 1, 0);

   numsrcstd = src->num_xforms - (src->final_xform_index >= 0);

   for(i=0;i<numsrcstd;i++) {

      /* When we copy the old xform, the pad is set to 0 */
      flam3_copy_xform(&dest->xform[i],&src->xform[i]);

      /* Copy the initial chaos from the src - the rest are already 1 */
      memcpy(dest->chaos[i], src->chaos[i], numsrcstd*sizeof(double));
      
   }   
   
   /* Add the final xform if necessary */
   if (dest_final_xform > 0) {
      flam3_add_xforms(dest, dest_final_xform, 1, 1);

	  flam3_xform * const xf = &dest->xform[dest->num_xforms-1];

      if (src->final_xform_enable > 0) {
      
         i = src->final_xform_index;
         
         flam3_copy_xform(xf, &src->xform[i]);
         
      } else {
         /* Interpolated-against final xforms need animate & color_speed set to 0.0 */
         xf->animate=0.0;
         xf->color_speed=0.0;
      }

   } else {
      dest->final_xform_index = -1;
      dest->final_xform_enable = 0;
   }

}

void clear_cp(flam3_genome *cp, int default_flag) {
	memset (cp, 0, sizeof (*cp));

    cp->gamma = 4.0;
    cp->vibrancy = 1.0;
    cp->contrast = 1.0;
    cp->brightness = 4.0;
    cp->pixels_per_unit = 50;
    cp->interpolation = flam3_interpolation_linear;

    if (default_flag==flam3_defaults_on) {
       /* If defaults are on, set to reasonable values */
       cp->highlight_power = -1.0;
       cp->width = 100;
       cp->height = 100;
       cp->gam_lin_thresh = 0.01;
       cp->interpolation_type = flam3_inttype_log;
       cp->palette_mode = PALETTE_MODE_STEP;

    } else {
       /* Defaults are off, so set to UN-reasonable values. */
       cp->highlight_power = -1.0;
       cp->zoom = 999999999;
       cp->width = -1;
       cp->height = -1;
       cp->gam_lin_thresh = -1;
//       cp->motion_exp = -999;
       cp->interpolation_type = -1;
       cp->palette_mode = -1;
    }

    cp->final_xform_index = -1;
}

flam3_genome *flam3_parse_xml2(const int fd, int default_flag, int *ncps,
		randctx * const rc) {

   xmlDocPtr doc; /* Parsed XML document tree */
   xmlNode *rootnode;
   int i;
   int loc_all_ncps=0;
   flam3_genome *loc_all_cp=NULL;
   char* locale = NULL;
   char* lorig  = setlocale(LC_NUMERIC, NULL);

   doc = xmlReadFd (fd, NULL, NULL, XML_PARSE_NONET);
   /* Check for errors */
   if (doc==NULL) {
      fprintf(stderr, "Failed to parse \n");
      return NULL;
   }

   /* What is the root node of the document? */
   rootnode = xmlDocGetRootElement(doc);
   
   // force use of "C" locale when writing reals.
   // first save away the current settings.
   if (lorig == NULL)
      fprintf(stderr, "error: couldn't get current locale\n");
   else {
      int slen = strlen(lorig) + 1;
      locale = (char*)malloc(slen);
      if (locale != NULL)
         memcpy(locale, lorig, slen);
   }
   if (setlocale(LC_NUMERIC, "C") == NULL)
      fprintf(stderr, "error: couldn't set C locale\n");

   /* Have to use &loc_all_cp since the memory gets allocated in scan_for_flame_nodes */
   scan_for_flame_nodes(rootnode, default_flag,&loc_all_cp,&loc_all_ncps, rc);

   // restore locale
   if (locale != NULL) {
      if (setlocale(LC_NUMERIC, locale) == NULL)
         fprintf(stderr, "error: couldn't replace locale settings\n");
      free(locale);
   }
   
   xmlFreeDoc(doc);

   *ncps = loc_all_ncps;
   
   /* Check to see if the first control point or the second-to-last */
   /* control point has interpolation="smooth".  This is invalid    */
   /* and should be reset to linear (with a warning).               */
   if (loc_all_ncps>=1) {
      if (loc_all_cp[0].interpolation == flam3_interpolation_smooth) {
         fprintf(stderr,"Warning: smooth interpolation cannot be used for first segment.\n"
                        "         switching to linear.\n");
         loc_all_cp[0].interpolation = flam3_interpolation_linear;
      }
   }
   
   if (loc_all_ncps>=2) {
      if (loc_all_cp[(loc_all_ncps)-2].interpolation == flam3_interpolation_smooth) {
         fprintf(stderr,"Warning: smooth interpolation cannot be used for last segment.\n"
                        "         switching to linear.\n");
         loc_all_cp[loc_all_ncps-2].interpolation = flam3_interpolation_linear;
      }
   }
   
   /* Finally, ensure that consecutive 'rotate' parameters never exceed */
   /* a difference of more than 180 degrees (+/-) for interpolation.    */
   /* An adjustment of +/- 360 degrees is made until this is true.      */
   if (*ncps>1) {
   
      for (i=1;i<*ncps;i++) {

         /* Only do this adjustment if we're not in compat mode */
         if (flam3_inttype_compat != loc_all_cp[i-1].interpolation_type
        && flam3_inttype_older != loc_all_cp[i-1].interpolation_type) {
      
            while (loc_all_cp[i].rotate < loc_all_cp[i-1].rotate-180)
               loc_all_cp[i].rotate += 360;
            
            while (loc_all_cp[i].rotate > loc_all_cp[i-1].rotate+180)
               loc_all_cp[i].rotate -= 360;
         }
      }
   }
   
   //Note that concurrent calls to flam3, if in parallel, potentially segfault
   //if this function is called.  technically it's required but it doesn't
   //leak memory continuously.
   //xmlCleanupParser();
   
   return loc_all_cp;
}

void flam3_print(FILE *f, flam3_genome *cp, char *extra_attributes) {
   int i,numstd;

   // force use of "C" locale when writing reals.
   // first save away the current settings.
   char* locale = NULL;
   char* lorig  = setlocale(LC_NUMERIC, NULL);
   if (lorig == NULL)
      fprintf(stderr, "error: couldn't get current locale\n");
   else {
      int slen = strlen(lorig) + 1;
      locale = (char*)malloc(slen);
      if (locale != NULL)
      memcpy(locale, lorig, slen);
   }
   if (setlocale(LC_NUMERIC, "C") == NULL)
      fprintf(stderr, "error: couldn't set C locale\n");

   
   fprintf(f, "<flame version=\"" PACKAGE "-" VERSION "\" time=\"%g\"", cp->time);
   
   if (cp->flame_name[0]!=0)
      fprintf(f, " name=\"%s\"",cp->flame_name);
   
   fprintf(f, " size=\"%d %d\"", cp->width, cp->height);
   fprintf(f, " center=\"%g %g\"", cp->center[0], cp->center[1]);
   fprintf(f, " scale=\"%g\"", cp->pixels_per_unit);

   if (cp->zoom != 0.0)
      fprintf(f, " zoom=\"%g\"", cp->zoom);

   fprintf(f, " rotate=\"%g\"", cp->rotate);

   fprintf(f, " brightness=\"%g\"", cp->brightness);
   fprintf(f, " gamma=\"%g\"", cp->gamma);
   
   fprintf(f, " highlight_power=\"%g\"", cp->highlight_power);
      
   fprintf(f, " vibrancy=\"%g\"", cp->vibrancy);
   fprintf(f, " gamma_threshold=\"%g\"", cp->gam_lin_thresh);
   
   if (PALETTE_MODE_STEP == cp->palette_mode)
      fprintf(f, " palette_mode=\"step\"");
   else if (PALETTE_MODE_LINEAR == cp->palette_mode)
      fprintf(f, " palette_mode=\"linear\"");

   if (flam3_interpolation_linear != cp->interpolation)
       fprintf(f, " interpolation=\"smooth\"");
       
   if (flam3_inttype_linear == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"linear\"");
   else if (flam3_inttype_log == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"log\"");
   else if (flam3_inttype_compat == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"old\"");
   else if (flam3_inttype_older == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"older\"");

   if (extra_attributes)
      fprintf(f, " %s", extra_attributes);

   fprintf(f, ">\n");

   if (cp->symmetry)
      fprintf(f, "   <symmetry kind=\"%d\"/>\n", cp->symmetry);
   
   numstd = cp->num_xforms - (cp->final_xform_index>=0);
   
   for (i = 0; i < cp->num_xforms; i++) {
   
      if (i==cp->final_xform_index)   
         flam3_print_xform(f, &cp->xform[i], 1, numstd, NULL);
      else 
         flam3_print_xform(f, &cp->xform[i], 0, numstd, cp->chaos[i]);

   }

   for (i = 0; i < cp->palette.count; i++) {
	  double4 rgba = cp->palette.color[i] * 255.0;

      fprintf(f, "   ");
      
      if (rgba[3] ==255.0) {
      
		 fprintf(f, "<color index=\"%d\" rgb=\"%.6g %.6g %.6g\"/>", i, rgba[0],
				 rgba[1], rgba[2]);
      } else {
		 fprintf(f, "   <color index=\"%d\" rgba=\"%.6g %.6g %.6g %.6g\"/>", i,
				 rgba[0], rgba[1], rgba[2], rgba[3]);
      }
//      if (i%4 == 3)
         fprintf(f, "\n");
         
   }

   fprintf(f, "</flame>\n");

   if (locale != NULL) {
      if (setlocale(LC_NUMERIC, locale) == NULL)
         fprintf(stderr, "error: couldn't restore locale settings\n");
      free(locale);
   }
}

#define PRINTNON(p) do { if (x->p != 0.0) fprintf(f, #p "=\"%f\" ",x->p); } while(0)

static void flam3_print_xform(FILE *f, flam3_xform *x, int final_flag,
		int numstd, double *chaos_row) {
	if (final_flag)
		fprintf(f, "   <finalxform ");
	else
		fprintf(f, "   <xform weight=\"%g\" ", x->density);

	fprintf(f, "color=\"%g\" ", x->color);

	fprintf(f, "color_speed=\"%g\" ", x->color_speed);

	if (!final_flag)
		fprintf(f, "animate=\"%g\" ", x->animate);

	for (unsigned int j = 0; j < flam3_nvariations; j++) {
		const double v = x->var[j];
		if (0.0 != v) {
			fprintf(f, "%s=\"%g\" ", flam3_variation_names[j], v);
			switch (j) {
				case VAR_BLOB:
					fprintf(f, "blob_low=\"%g\" ", x->blob_low);
					fprintf(f, "blob_high=\"%g\" ", x->blob_high);
					fprintf(f, "blob_waves=\"%g\" ", x->blob_waves);
					break;

				case VAR_PDJ:
					fprintf(f, "pdj_a=\"%g\" ", x->pdj_a);
					fprintf(f, "pdj_b=\"%g\" ", x->pdj_b);
					fprintf(f, "pdj_c=\"%g\" ", x->pdj_c);
					fprintf(f, "pdj_d=\"%g\" ", x->pdj_d);
					break;

				case VAR_FAN2:
					fprintf(f, "fan2_x=\"%g\" ", x->fan2_x);
					fprintf(f, "fan2_y=\"%g\" ", x->fan2_y);
					break;

				case VAR_RINGS2:
					fprintf(f, "rings2_val=\"%g\" ", x->rings2_val);
					break;

				case VAR_PERSPECTIVE:
					fprintf(f, "perspective_angle=\"%g\" ", x->perspective_angle);
					fprintf(f, "perspective_dist=\"%g\" ", x->perspective_dist);
					break;

				case VAR_JULIAN:
					fprintf(f, "julian_power=\"%g\" ", x->julian_power);
					fprintf(f, "julian_dist=\"%g\" ", x->julian_dist);
					break;

				case VAR_JULIASCOPE:
					fprintf(f, "juliascope_power=\"%g\" ", x->juliascope_power);
					fprintf(f, "juliascope_dist=\"%g\" ", x->juliascope_dist);
					break;

				case VAR_RADIAL_BLUR:
					fprintf(f, "radial_blur_angle=\"%g\" ", x->radial_blur_angle);
					break;

				case VAR_PIE:
					fprintf(f, "pie_slices=\"%g\" ", x->pie_slices);
					fprintf(f, "pie_rotation=\"%g\" ", x->pie_rotation);
					fprintf(f, "pie_thickness=\"%g\" ", x->pie_thickness);
					break;

				case VAR_NGON:
					fprintf(f, "ngon_sides=\"%g\" ", x->ngon_sides);
					fprintf(f, "ngon_power=\"%g\" ", x->ngon_power);
					fprintf(f, "ngon_corners=\"%g\" ", x->ngon_corners);
					fprintf(f, "ngon_circle=\"%g\" ", x->ngon_circle);
					break;

				case VAR_CURL:
					fprintf(f, "curl_c1=\"%g\" ", x->curl_c1);
					fprintf(f, "curl_c2=\"%g\" ", x->curl_c2);
					break;

				case VAR_RECTANGLES:
					fprintf(f, "rectangles_x=\"%g\" ", x->rectangles_x);
					fprintf(f, "rectangles_y=\"%g\" ", x->rectangles_y);
					break;

				case VAR_DISC2:
					fprintf(f, "disc2_rot=\"%g\" ", x->disc2_rot);
					fprintf(f, "disc2_twist=\"%g\" ", x->disc2_twist);
					break;

				case VAR_SUPER_SHAPE:
					fprintf(f, "super_shape_rnd=\"%g\" ", x->super_shape_rnd);
					fprintf(f, "super_shape_m=\"%g\" ", x->super_shape_m);
					fprintf(f, "super_shape_n1=\"%g\" ", x->super_shape_n1);
					fprintf(f, "super_shape_n2=\"%g\" ", x->super_shape_n2);
					fprintf(f, "super_shape_n3=\"%g\" ", x->super_shape_n3);
					fprintf(f, "super_shape_holes=\"%g\" ", x->super_shape_holes);
					break;

				case VAR_FLOWER:
					fprintf(f, "flower_petals=\"%g\" ", x->flower_petals);
					fprintf(f, "flower_holes=\"%g\" ", x->flower_holes);
					break;

				case VAR_CONIC:
					fprintf(f, "conic_eccentricity=\"%g\" ", x->conic_eccentricity);
					fprintf(f, "conic_holes=\"%g\" ", x->conic_holes);
					break;

				case VAR_PARABOLA:
					fprintf(f, "parabola_height=\"%g\" ", x->parabola_height);
					fprintf(f, "parabola_width=\"%g\" ", x->parabola_width);
					break;

				case VAR_BENT2:
					fprintf(f, "bent2_x=\"%g\" ", x->bent2_x);
					fprintf(f, "bent2_y=\"%g\" ", x->bent2_y);
					break;

				case VAR_BIPOLAR:
					fprintf(f, "bipolar_shift=\"%g\" ", x->bipolar_shift);
					break;

				case VAR_CELL:
					fprintf(f, "cell_size=\"%g\" ", x->cell_size);
					break;

				case VAR_CPOW:
					fprintf(f, "cpow_i=\"%g\" ", x->cpow_i);
					fprintf(f, "cpow_r=\"%g\" ", x->cpow_r);
					fprintf(f, "cpow_power=\"%g\" ", x->cpow_power);
					break;

				case VAR_CURVE:
					fprintf(f, "curve_xamp=\"%g\" ", x->curve_xamp);
					fprintf(f, "curve_yamp=\"%g\" ", x->curve_yamp);
					fprintf(f, "curve_xlength=\"%g\" ", x->curve_xlength);
					fprintf(f, "curve_ylength=\"%g\" ", x->curve_ylength);
					break;

				case VAR_ESCHER:
					fprintf(f, "escher_beta=\"%g\" ", x->escher_beta);
					break;

				case VAR_LAZYSUSAN:
					fprintf(f, "lazysusan_x=\"%g\" ", x->lazysusan_x);
					fprintf(f, "lazysusan_y=\"%g\" ", x->lazysusan_y);
					fprintf(f, "lazysusan_spin=\"%g\" ", x->lazysusan_spin);
					fprintf(f, "lazysusan_space=\"%g\" ", x->lazysusan_space);
					fprintf(f, "lazysusan_twist=\"%g\" ", x->lazysusan_twist);
					break;

				case VAR_MODULUS:
					fprintf(f, "modulus_x=\"%g\" ", x->modulus_x);
					fprintf(f, "modulus_y=\"%g\" ", x->modulus_y);
					break;

				case VAR_OSCILLOSCOPE:
					fprintf(f, "oscilloscope_separation=\"%g\" ", x->oscope_separation);
					fprintf(f, "oscilloscope_frequency=\"%g\" ", x->oscope_frequency);
					fprintf(f, "oscilloscope_amplitude=\"%g\" ", x->oscope_amplitude);
					fprintf(f, "oscilloscope_damping=\"%g\" ", x->oscope_damping);
					break;

				case VAR_POPCORN2:
					fprintf(f, "popcorn2_x=\"%g\" ", x->popcorn2_x);
					fprintf(f, "popcorn2_y=\"%g\" ", x->popcorn2_y);
					fprintf(f, "popcorn2_c=\"%g\" ", x->popcorn2_c);
					break;

				case VAR_SPLIT:
					fprintf(f, "split_xsize=\"%g\" ", x->split_xsize);
					fprintf(f, "split_ysize=\"%g\" ", x->split_ysize);
					break;

				case VAR_SPLITS:
					fprintf(f, "splits_x=\"%g\" ", x->splits_x);
					fprintf(f, "splits_y=\"%g\" ", x->splits_y);
					break;

				case VAR_STRIPES:
					fprintf(f, "stripes_space=\"%g\" ", x->stripes_space);
					fprintf(f, "stripes_warp=\"%g\" ", x->stripes_warp);
					break;

				case VAR_WEDGE:
					fprintf(f, "wedge_angle=\"%g\" ", x->wedge_angle);
					fprintf(f, "wedge_hole=\"%g\" ", x->wedge_hole);
					fprintf(f, "wedge_count=\"%g\" ", x->wedge_count);
					fprintf(f, "wedge_swirl=\"%g\" ", x->wedge_swirl);
					break;

				case VAR_WEDGE_JULIA:
					fprintf(f, "wedge_julia_angle=\"%g\" ", x->wedge_julia_angle);
					fprintf(f, "wedge_julia_count=\"%g\" ", x->wedge_julia_count);
					fprintf(f, "wedge_julia_power=\"%g\" ", x->wedge_julia_power);
					fprintf(f, "wedge_julia_dist=\"%g\" ", x->wedge_julia_dist);
					break;

				case VAR_WEDGE_SPH:
					fprintf(f, "wedge_sph_angle=\"%g\" ", x->wedge_sph_angle);
					fprintf(f, "wedge_sph_hole=\"%g\" ", x->wedge_sph_hole);
					fprintf(f, "wedge_sph_count=\"%g\" ", x->wedge_sph_count);
					fprintf(f, "wedge_sph_swirl=\"%g\" ", x->wedge_sph_swirl);
					break;

				case VAR_WHORL:
					fprintf(f, "whorl_inside=\"%g\" ", x->whorl_inside);
					fprintf(f, "whorl_outside=\"%g\" ", x->whorl_outside);
					break;

				case VAR_WAVES2:
					fprintf(f, "waves2_scalex=\"%g\" ", x->waves2_scalex);
					fprintf(f, "waves2_scaley=\"%g\" ", x->waves2_scaley);
					fprintf(f, "waves2_freqx=\"%g\" ", x->waves2_freqx);
					fprintf(f, "waves2_freqy=\"%g\" ", x->waves2_freqy);
					break;

				case VAR_AUGER:
					fprintf(f, "auger_sym=\"%g\" ", x->auger_sym);
					fprintf(f, "auger_weight=\"%g\" ", x->auger_weight);
					fprintf(f, "auger_freq=\"%g\" ", x->auger_freq);
					fprintf(f, "auger_scale=\"%g\" ", x->auger_scale);
					break;

				case VAR_FLUX:
					fprintf(f, "flux_spread=\"%g\" ", x->flux_spread);
					break;

				case VAR_MOBIUS:
					fprintf(f, "mobius_re_a=\"%g\" ", x->mobius_re_a);
					fprintf(f, "mobius_im_a=\"%g\" ", x->mobius_im_a);
					fprintf(f, "mobius_re_b=\"%g\" ", x->mobius_re_b);
					fprintf(f, "mobius_im_b=\"%g\" ", x->mobius_im_b);
					fprintf(f, "mobius_re_c=\"%g\" ", x->mobius_re_c);
					fprintf(f, "mobius_im_c=\"%g\" ", x->mobius_im_c);
					fprintf(f, "mobius_re_d=\"%g\" ", x->mobius_re_d);
					fprintf(f, "mobius_im_d=\"%g\" ", x->mobius_im_d);
					break;

				case VAR_SEPARATION:
					fprintf(f, "separation_x=\"%g\" ", x->separation_x);
					fprintf(f, "separation_y=\"%g\" ", x->separation_y);
					fprintf(f, "separation_xinside=\"%g\" ", x->separation_xinside);
					fprintf(f, "separation_yinside=\"%g\" ", x->separation_yinside);
					break;

				default:
					/* pass */
					break;
			}
		}
	}

	fprintf(f, "coefs=\"");
	for (unsigned int j = 0; j < 3; j++) {
		if (j) fprintf(f, " ");
		fprintf(f, "%g %g", x->c[j][0], x->c[j][1]);
	}
	fprintf(f, "\"");

	if (!id_matrix(x->post)) {
		fprintf(f, " post=\"");
		for (unsigned int j = 0; j < 3; j++) {
			if (j) fprintf(f, " ");
			fprintf(f, "%g %g", x->post[j][0], x->post[j][1]);
		}
		fprintf(f, "\"");
	}

	if (!final_flag) {

		/* Print out the chaos row for this xform */
		int numcols = numstd;

		while (numcols > 0 && chaos_row[numcols-1]==1.0)
			numcols--;

		if (numcols>0) {
			fprintf(f, " chaos=\"");
			for (unsigned int j=0;j<numcols;j++)
				fprintf(f, "%g ",chaos_row[j]);
			fprintf(f, "\"");
		}
	}

	fprintf(f, " opacity=\"%g\" />\n",x->opacity);
}

static double round6(double x) {
  x *= 1e6;
  if (x < 0) x -= 1.0;
  return 1e-6*(int)(x+0.5);
}

static int compare_xforms(const void *av, const void *bv) {
   flam3_xform *a = (flam3_xform *) av;
   flam3_xform *b = (flam3_xform *) bv;
   double aa[2][2];
   double bb[2][2];
   double ad, bd;

   aa[0][0] = a->c[0][0];
   aa[0][1] = a->c[0][1];
   aa[1][0] = a->c[1][0];
   aa[1][1] = a->c[1][1];
   bb[0][0] = b->c[0][0];
   bb[0][1] = b->c[0][1];
   bb[1][0] = b->c[1][0];
   bb[1][1] = b->c[1][1];
   ad = det_matrix(aa);
   bd = det_matrix(bb);

   if (a->color_speed > b->color_speed) return 1;
   if (a->color_speed < b->color_speed) return -1;
   if (a->color_speed) {
      if (ad < 0) return -1;
      if (bd < 0) return 1;
      ad = atan2(a->c[0][0], a->c[0][1]);
      bd = atan2(b->c[0][0], b->c[0][1]);
   }

   if (ad < bd) return -1;
   if (ad > bd) return 1;
   return 0;
}

/* sym=2 or more means rotational
   sym=1 means identity, ie no symmetry
   sym=0 means pick a random symmetry (maybe none)
   sym=-1 means bilateral (reflection)
   sym=-2 or less means rotational and reflective
*/
void flam3_add_symmetry(flam3_genome *cp, int sym, randctx * const rc) {
   int i, j, k;
   double a;
   int result = 0;

   if (0 == sym) {
      static int sym_distrib[] = {
         -4, -3,
         -2, -2, -2,
         -1, -1, -1,
         2, 2, 2,
         3, 3,
         4, 4,
      };
      if (rand_bool(rc)) {
         sym = rand_distrib(rc, sym_distrib);
      } else if (rand_mod(rc, 32)) {
         sym = rand_mod(rc, 13)-6;
      } else {
         sym = rand_mod(rc, 51)-25;
      }
   }

   if (1 == sym || 0 == sym) return;

   cp->symmetry = sym;

   if (sym < 0) {

      i = cp->num_xforms;
      if (cp->final_xform_enable)
         i -= 1;

      flam3_add_xforms(cp,1,0,0);

	  flam3_xform * const xf = &cp->xform[i];

      xf->density = 1.0;
      xf->color_speed = 0.0;
      xf->animate = 0.0;
      xf->var[0] = 1.0;
      for (j = 1; j < flam3_nvariations; j++)
         xf->var[j] = 0;
      xf->color = 1.0;
      xf->c[0][0] = -1.0;
      xf->c[0][1] = 0.0;
      xf->c[1][0] = 0.0;
      xf->c[1][1] = 1.0;
      xf->c[2][0] = 0.0;
      xf->c[2][1] = 0.0;

      result++;
      sym = -sym;
   }

   a = 2*M_PI/sym;

   for (k = 1; k < sym; k++) {

      i = cp->num_xforms;
      if (cp->final_xform_enable)
         i -= 1;

      flam3_add_xforms(cp, 1, 0,0);

	  flam3_xform * const xf = &cp->xform[i];

      xf->density = 1.0;
      xf->color_speed = 0.0;
      xf->animate = 0.0;
      xf->var[0] = 1.0;
      for (j = 1; j < flam3_nvariations; j++)
         xf->var[j] = 0;
      xf->color = (sym<3) ? 0.0 : ((k-1.0)/(sym-2.0));
      xf->c[0][0] = round6(cos(k*a));
      xf->c[0][1] = round6(sin(k*a));
      xf->c[1][0] = round6(-cp->xform[i].c[0][1]);
      xf->c[1][1] = cp->xform[i].c[0][0];
      xf->c[2][0] = 0.0;
      xf->c[2][1] = 0.0;

      result++;
   }   
   
   qsort((char *) &cp->xform[cp->num_xforms-result], result,
      sizeof(flam3_xform), compare_xforms);
      
}

void flam3_cross(flam3_genome *cp0, flam3_genome *cp1, flam3_genome *out, int cross_mode, randctx *rc) {

   int i,j, rb;

   if (cross_mode == CROSS_NOT_SPECIFIED) {
   
      double s = rand_d01(rc);
      
      if (s < 0.1)
         cross_mode = CROSS_UNION;
      else if (s < 0.2)
         cross_mode = CROSS_INTERPOLATE;
      else
         cross_mode = CROSS_ALTERNATE;

   }
   
   if (cross_mode == CROSS_UNION) {
   
      flam3_xform mycopy;
   
      /* Make a copy of cp0 */
      flam3_copy(out, cp0);
      
      for (j=0;j<cp1->num_xforms;j++) {
         /* Skip over the final xform, if it's present.    */
         /* Default behavior keeps the final from parent0. */
         if (cp1->final_xform_index == j)		     
            continue;

         i = out->num_xforms;
         if (out->final_xform_enable)
            i -= 1;

         flam3_add_xforms(out, 1, 0, 0);
         flam3_copy_xform(&out->xform[i],&cp1->xform[j]);
      }
      
      /* Put the final xform last (if there is one) */
      /* We do not need to do complicated xform copies here since we're just moving them around */
      if (out->final_xform_index >= 0) {
         mycopy = out->xform[out->final_xform_index];
         out->xform[out->final_xform_index] = out->xform[out->num_xforms-1];
         out->xform[out->num_xforms-1] = mycopy;
         out->final_xform_index = out->num_xforms-1;
      }
      
      
   } else if (cross_mode == CROSS_INTERPOLATE) {
   
      /* linearly interpolate somewhere between the two */
      flam3_genome parents[2];
      double t = rand_d01(rc);

      memset(parents, 0, 2*sizeof(flam3_genome));

      flam3_copy(&(parents[0]), cp0);
      flam3_copy(&(parents[1]), cp1);

      parents[0].time = 0.0;
      parents[1].time = 1.0;
      flam3_interpolate(parents, 2, t, 0, out);
      
      clear_cp(&parents[0],flam3_defaults_on);
      clear_cp(&parents[1],flam3_defaults_on);
      
   } else {
   
      /* alternate mode */
      int got0, got1, used_parent;
      char *trystr;

      trystr = calloc(4 * (cp0->num_xforms + cp1->num_xforms), sizeof(char));

      /* each xform comes from a random parent, possible for an entire parent to be excluded */
      do {

         trystr[0] = 0;
         got0 = got1 = 0;
         rb = rand_bool(rc);

         /* Copy the parent, sorting the final xform to the end if it's present. */
         if (rb)
            flam3_copyx(out, cp1, cp1->num_xforms - (cp1->final_xform_index > 0), cp1->final_xform_enable);
         else
            flam3_copyx(out, cp0, cp0->num_xforms - (cp0->final_xform_index > 0), cp0->final_xform_enable);

         used_parent = rb;

         /* Only replace non-final xforms */

         for (i = 0; i < out->num_xforms - out->final_xform_enable; i++) {
            rb = rand_bool(rc);

            /* Replace xform if bit is 1 */
            if (rb==1) {
               if (used_parent==0) {
                  if (i < cp1->num_xforms && cp1->xform[i].density > 0) {
                     flam3_copy_xform(&out->xform[i],&cp1->xform[i]);
                     got1 = 1;
                  } else {
                     got0 = 1;
                  }
               } else {
                  if (i < cp0->num_xforms && cp0->xform[i].density > 0) {
                     flam3_copy_xform(&out->xform[i],&cp0->xform[i]);
                     got0 = 1;
                  } else {
                     got1 = 1;
                  }
               }
            } else {
               if (used_parent)
                  got1 = 1;
               else
                  got0 = 1;
            }

         }
             
         if (used_parent==0 && cp0->final_xform_enable)
            got0 = 1;
         else if (used_parent==1 && cp1->final_xform_enable)
            got1 = 1;
               
      } while ((i > 1) && !(got0 && got1));

      free(trystr);
   }
   
   /* reset color coords */
   for (i = 0; i < out->num_xforms; i++) {
      out->xform[i].color = i&1;
   }
               
   /* Potentially genetically cross the two colormaps together */
   if (rand_d01(rc) < 0.4) {
                              
      /* Select the starting parent */
      int startParent=rand_bool(rc);
      int ci;
                        
      /* Loop over the entries, switching to the other parent 1% of the time */
      for (ci=0;ci<out->palette.count;ci++) {
         if (rand_d01(rc)<.01) {
            startParent = 1-startParent;
         }
                     
         out->palette.color[ci] = startParent ? cp1->palette.color[ci] :
		 		cp0->palette.color[ci];
      }
   }

}

#if 0
void flam3_mutate(flam3_genome *cp, int mutate_mode, int *ivars, int ivars_n, int sym, double speed, const palette_collection * const pc, randctx *rc) {

   double randselect;
   flam3_genome mutation;
   int i,j,done;
   
   /* If mutate_mode = -1, choose a random mutation mode */
   if (mutate_mode == MUTATE_NOT_SPECIFIED) {
   
      randselect = rand_d01(rc);
      
      if (randselect < 0.1)
         mutate_mode = MUTATE_ALL_VARIATIONS;
      else if (randselect < 0.3)
         mutate_mode = MUTATE_ONE_XFORM_COEFS;
      else if (randselect < 0.5)
         mutate_mode = MUTATE_ADD_SYMMETRY;
      else if (randselect < 0.6)
         mutate_mode = MUTATE_POST_XFORMS;
      else if (randselect < 0.7)
         mutate_mode = MUTATE_COLOR_PALETTE;
      else if (randselect < 0.8)
         mutate_mode = MUTATE_DELETE_XFORM;
      else
         mutate_mode = MUTATE_ALL_COEFS;

   }
   
   memset(&mutation, 0, sizeof(flam3_genome));
      
   if (mutate_mode == MUTATE_ALL_VARIATIONS) {
   
      do {
         /* Create a random flame, and use the variations */
         /* to replace those in the original              */
         flam3_random(&mutation, cp->num_xforms, pc, rc);
         for (i = 0; i < cp->num_xforms; i++) {
            for (j = 0; j < flam3_nvariations; j++) {
               if (cp->xform[i].var[j] != mutation.xform[i].var[j]) {
               
                  /* Copy the new var weights */
                  cp->xform[i].var[j] = mutation.xform[i].var[j];

                  /* Copy parameters for this variation only */
                  flam3_copy_params(&(cp->xform[i]),&(mutation.xform[i]),j);

                  done = 1;
               }
            }
         }
      } while (!done);
      
   } else if (mutate_mode == MUTATE_ONE_XFORM_COEFS) {
   
      int modxf;
   
      /* Generate a 2-xform random */
      flam3_random(&mutation, ivars, ivars_n, sym, 2, pc, rc);
      
      /* Which xform do we mutate? */
      modxf = rand_mod(rc, cp->num_xforms);
      
      /* if less than 3 xforms, then change only the translation part */
      if (2 >= cp->num_xforms) {
         for (j = 0; j < 2; j++)
            cp->xform[modxf].c[2][j] = mutation.xform[0].c[2][j];
      } else {
         for (i = 0; i < 3; i++)
            for (j = 0; j < 2; j++)
               cp->xform[modxf].c[i][j] = mutation.xform[0].c[i][j];
      }
      
   } else if (mutate_mode == MUTATE_ADD_SYMMETRY) {
   
      flam3_add_symmetry(cp, 0, rc);
      
   } else if (mutate_mode == MUTATE_POST_XFORMS) {
   
      int b = 1 + rand_mod(rc, 6);
      int same = rand_mod(rc, 4); /* 25% chance of using the same post for all of them */
      
      for (i = 0; i < cp->num_xforms; i++) {
         int copy = (i > 0) && same;

         if (copy) { /* Copy the post from the first xform to the rest of them */
            for (j = 0; j < 3; j++) {
               cp->xform[i].post[j][0] = cp->xform[0].post[j][0];
               cp->xform[i].post[j][1] = cp->xform[0].post[j][1];
            }

         } else {

            if (b&1) { /* 50% chance */
            
               double f = M_PI * rand_d11(rc);
               double t[2][2];

               t[0][0] = (cp->xform[i].c[0][0] * cos(f) + cp->xform[i].c[0][1] * -sin(f));
               t[0][1] = (cp->xform[i].c[0][0] * sin(f) + cp->xform[i].c[0][1] * cos(f));
               t[1][0] = (cp->xform[i].c[1][0] * cos(f) + cp->xform[i].c[1][1] * -sin(f));
               t[1][1] = (cp->xform[i].c[1][0] * sin(f) + cp->xform[i].c[1][1] * cos(f));

               cp->xform[i].c[0][0] = t[0][0];
               cp->xform[i].c[0][1] = t[0][1];
               cp->xform[i].c[1][0] = t[1][0];
               cp->xform[i].c[1][1] = t[1][1];

               f *= -1.0;

               t[0][0] = (cp->xform[i].post[0][0] * cos(f) + cp->xform[i].post[0][1] * -sin(f));
               t[0][1] = (cp->xform[i].post[0][0] * sin(f) + cp->xform[i].post[0][1] * cos(f));
               t[1][0] = (cp->xform[i].post[1][0] * cos(f) + cp->xform[i].post[1][1] * -sin(f));
               t[1][1] = (cp->xform[i].post[1][0] * sin(f) + cp->xform[i].post[1][1] * cos(f));

               cp->xform[i].post[0][0] = t[0][0];
               cp->xform[i].post[0][1] = t[0][1];
               cp->xform[i].post[1][0] = t[1][0];
               cp->xform[i].post[1][1] = t[1][1];

            }

            if (b&2) { /* 33% chance */
            
               double f = 0.2 + rand_d01(rc);
               double g = 0.2 + rand_d01(rc);

               if (rand_bool(rc))
                  f = 1.0 / f;
               
               if (rand_bool(rc))
                  g = f;
               else {               
                  if (rand_bool(rc))
                     g = 1.0 / g;
               }

               cp->xform[i].c[0][0] /= f;
               cp->xform[i].c[0][1] /= f;
               cp->xform[i].c[1][1] /= g;
               cp->xform[i].c[1][0] /= g;
               cp->xform[i].post[0][0] *= f;
               cp->xform[i].post[1][0] *= f;
               cp->xform[i].post[0][1] *= g;
               cp->xform[i].post[1][1] *= g;
            }

            if (b&4) { /* 16% chance */

               double f = rand_d11(rc);
               double g = rand_d11(rc);

               cp->xform[i].c[2][0] -= f;
               cp->xform[i].c[2][1] -= g;
               cp->xform[i].post[2][0] += f;
               cp->xform[i].post[2][1] += g;
            }
         }
      }
   } else if (mutate_mode == MUTATE_COLOR_PALETTE) {
   
      double s = rand_d01(rc);

      if (s < 0.4) { /* randomize xform color coords */
      
         flam3_improve_colors(cp, 100, 0, 10, pc, rc);
         
      } else if (s < 0.8) { /* randomize xform color coords and palette */
      
         flam3_improve_colors(cp, 25, 1, 10, pc, rc);
         
      } else { /* randomize palette only */
		 const palette * const p = palette_random (pc, rc);
		 assert (p != NULL);
		 palette_copy (p, &cp->palette);
		 palette_rotate_hue (&cp->palette, cp->hue_rotation);
      }
   } else if (mutate_mode == MUTATE_DELETE_XFORM) {
   
      int nx = rand_mod(rc, cp->num_xforms);

      if (cp->num_xforms > 1)
         flam3_delete_xform(cp,nx);

   } else { /* MUTATE_ALL_COEFS */ 
   
      int x;
      flam3_random(&mutation, ivars, ivars_n, sym, cp->num_xforms, pc, rc);

      /* change all the coefs by a fraction of the random */
      for (x = 0; x < cp->num_xforms; x++) {
         for (i = 0; i < 3; i++) {
            for (j = 0; j < 2; j++) {
               cp->xform[x].c[i][j] += speed * mutation.xform[x].c[i][j];

            }
         }
         /* Eventually, we can mutate the parametric variation coefs here. */
      }
   }
   
   clear_cp(&mutation,flam3_defaults_on);

}
#endif

static int sort_by_x(const void *av, const void *bv) {
    double4 a = *((double4 *) av);
    double4 b = *((double4 *) bv);
    if (a[0] < b[0]) return -1;
    if (a[0] > b[0]) return 1;
    return 0;
}

static int sort_by_y(const void *av, const void *bv) {
    double4 a = *((double4 *) av);
    double4 b = *((double4 *) bv);
    if (a[1] < b[1]) return -1;
    if (a[1] > b[1]) return 1;
    return 0;
}

/*
 * find a 2d bounding box that does not enclose eps of the fractal density
 * in each compass direction.
 */
int flam3_estimate_bounding_box(flam3_genome *cp, double eps, int maxsamples,
             double *bmin, double *bmax, randctx *rc) {
   int i;
   int low_target, high_target;
   double min[2], max[2];
   double4 *points;
   int bv;
   unsigned short *xform_distrib;

   if (maxsamples <= 0) maxsamples = 10000;

   int ret = posix_memalign ((void **) &points, sizeof (*points), sizeof(*points) * maxsamples);
   assert (ret == 0 && points != NULL);

   if (prepare_precalc_flags(cp))
      return(-1);
   xform_distrib = flam3_create_xform_distrib(cp);
   if (xform_distrib==NULL)
      return(-1);
   for (unsigned int i = 0; i < cp->num_xforms; i++) {
	   xform_precalc (&cp->xform[i]);
   }

   iterator iter;
   iterator_init (&iter, cp, xform_distrib, rc);

   /* throw away fuse steps */
   for (unsigned int i = 0; i < 20; i++) {
	   double4 p;
	   iterator_step (&iter, &p, rc);
   }

   /* actual iterations */
   unsigned int samples = 0;
   for (unsigned int i = 0; i < maxsamples; i++) {
	   if (iterator_step (&iter, &points[samples], rc)) {
		   ++samples;
	   }
   }

   free(xform_distrib);
      
   if ( bv/(double)samples > eps )
      eps = 3*bv/(double)samples;
   
   if ( eps > 0.3 )
      eps = 0.3;
      
   low_target = (int)(samples * eps);
   high_target = samples - low_target;

   
   min[0] = min[1] =  1e10;
   max[0] = max[1] = -1e10;

   for (i = 0; i < samples; i++) {
      const double4 p = points[i];
      if (p[0] < min[0]) min[0] = p[0];
      if (p[1] < min[1]) min[1] = p[1];
      if (p[0] > max[0]) max[0] = p[0];
      if (p[1] > max[1]) max[1] = p[1];
   }

   if (low_target == 0) {
      bmin[0] = min[0];
      bmin[1] = min[1];
      bmax[0] = max[0];
      bmax[1] = max[1];
      free(points);
      return(bv);
   }

   qsort(points, samples, sizeof(double4), sort_by_x);
   bmin[0] = points[low_target][0];
   bmax[0] = points[high_target][0];

   qsort(points, samples, sizeof(double4), sort_by_y);
   bmin[1] = points[low_target][1];
   bmax[1] = points[high_target][1];
   free(points);
   
   return(bv);
}

