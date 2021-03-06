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

#include "interpolation.h"
#include "palettes.h"

#include <assert.h>

#define INTERP(x)  do { result->x = 0.0; \
   for (k = 0; k < ncp; k++) result->x += c[k] * cpi[k].x; } while(0)

#define INTERI(x)  do { double tt = 0.0; \
   for (k = 0; k < ncp; k++) tt += c[k] * cpi[k].x; \
   result->x = (int)rint(tt); } while(0)

int id_matrix(double2 s[3]) {
   return
      (s[0][0] == 1.0) &&
      (s[0][1] == 0.0) &&
      (s[1][0] == 0.0) &&
      (s[1][1] == 1.0) &&
      (s[2][0] == 0.0) &&
      (s[2][1] == 0.0);
}

static void clear_matrix(double2 m[3]) {
   const double2 zero = (double2) { 0.0, 0.0 };
   m[0] = zero;
   m[1] = zero;
   m[2] = zero;
}

static void sum_matrix(double s, const double2 m1[3], double2 m2[3]) {

   m2[0] += s * m1[0];
   m2[1] += s * m1[1];
   m2[2] += s * m1[2];
}


#if 0
void interpolate_cmap(flam3_palette cmap, double blend,
                      int index0, double hue0, int index1, double hue1,
					  randctx * const rc) {
                 
   flam3_palette p0,p1;
   int i, j, rcode;

   rcode = flam3_get_palette(index0, p0, hue0, rc);
   if (rcode<0)
      fprintf(stderr,"unable to retrieve palette %d, setting to white\n", index0);
   rcode = flam3_get_palette(index1, p1, hue1, rc);
   if (rcode<0)
      fprintf(stderr,"unable to retrieve palette %d, setting to white\n", index1);

   for (i = 0; i < 256; i++) {
      double4 t, s;
	  double t4, s4;
    
      s = rgb2hsv(vector_d4 (p0[i].color));
      t = rgb2hsv(vector_d4 (p1[i].color));
      
      s[3] = p0[i].color[3];
      t[3] = p1[i].color[3];
      
      s4 = p0[i].index;
      t4 = p1[i].index;
    
      for (j = 0; j < 4; j++)
         t[j] = ((1.0-blend) * s[j]) + (blend * t[j]);
      t4 = ((1.0-blend) * s4) + (blend * t4);
         
	  const double4 c = hsv2rgb(t);
      cmap[i].color[0] = c[0];
      cmap[i].color[1] = c[1];
      cmap[i].color[2] = c[2];
      cmap[i].color[3] = t[3];
      cmap[i].index = t4;
   }
}
#endif

static void interp_and_convert_back(double *c, int ncps, int xfi, double cxang[4][2], 
                             double cxmag[4][2], double cxtrn[4][2],double store_array[3][2]) {

   int i,col;
   
   double accang[2],accmag[2];
   double expmag;
   int accmode[2];
   
   accang[0] = 0.0;
   accang[1] = 0.0;
   accmag[0] = 0.0;
   accmag[1] = 0.0;

   accmode[0]=accmode[1]=0;
   
   /* accumulation mode defaults to logarithmic, but in special */
   /* cases we want to switch to linear accumulation            */
   for (col=0; col<2; col++) {
      for (i=0; i<ncps; i++) {
         if (log(cxmag[i][col])<-10)
            accmode[col]=1; // Mode set to linear interp
      }
   }
   
   for (i=0; i<ncps; i++) {
      for (col=0; col<2; col++) {
      
         accang[col] += c[i] * cxang[i][col];
         
         if (accmode[col]==0)
            accmag[col] += c[i] * log(cxmag[i][col]);
         else 
            accmag[col] += c[i] * (cxmag[i][col]);
            
         /* translation is ready to go */
         store_array[2][col] += c[i] * cxtrn[i][col];
      }
   }
   
   /* Convert the angle back to rectangular */
   for (col=0;col<2;col++) {
      if (accmode[col]==0)
         expmag = exp(accmag[col]);
      else
         expmag = accmag[col];
      
      store_array[col][0] = expmag * cos(accang[col]);
      store_array[col][1] = expmag * sin(accang[col]);
   }
   
}

static void convert_linear_to_polar(flam3_genome *cp, int ncps, int xfi, int cflag, 
                             double cxang[4][2], double cxmag[4][2], double cxtrn[4][2]) {

   double c1[2],d,t,refang;
   int col,k;
   int zlm[2];

   for (k=0; k<ncps;k++) {

      /* Establish the angles and magnitudes for each component */
      /* Keep translation linear */
      zlm[0]=zlm[1]=0;
      for (col=0;col<2;col++) {
      
         if (cflag==0) {
            c1[0] = cp[k].xform[xfi].c[col][0];
            c1[1] = cp[k].xform[xfi].c[col][1];
            t = cp[k].xform[xfi].c[2][col];            
         } else {
            c1[0] = cp[k].xform[xfi].post[col][0];
            c1[1] = cp[k].xform[xfi].post[col][1];
            t = cp[k].xform[xfi].post[2][col];
         }
         
         cxang[k][col] = atan2(c1[1],c1[0]);
         cxmag[k][col] = sqrt(c1[0]*c1[0] + c1[1]*c1[1]);
         
         if (cxmag[k][col]== 0.0)
            zlm[col]=1;
         
         cxtrn[k][col] = t;
      }
      
      if (zlm[0]==1 && zlm[1]==0)
         cxang[k][0] = cxang[k][1];
      else if (zlm[0]==0 && zlm[1]==1)
         cxang[k][1] = cxang[k][0];
      
   }
   
   /* Make sure the rotation is the shorter direction around the circle */
   /* by adjusting each angle in succession, and rotate clockwise if 180 degrees */
   for (col=0; col<2; col++) {
      for (k=1;k<ncps;k++) {

         /* Adjust angles differently if we have an asymmetric case */   
         if (cp[k].xform[xfi].wind[col]>0 && cflag==0) {

            /* Adjust the angles to make sure that it's within wind:wind+2pi */
            refang = cp[k].xform[xfi].wind[col] - 2*M_PI;

            /* Make sure both angles are within [refang refang+2*pi] */
            while(cxang[k-1][col] < refang)
                  cxang[k-1][col] += 2*M_PI;
            
            while(cxang[k-1][col] > refang + 2*M_PI)
                  cxang[k-1][col] -= 2*M_PI;
                  
            while(cxang[k][col] < refang)
                  cxang[k][col] += 2*M_PI;
            
            while(cxang[k][col] > refang + 2*M_PI)
                  cxang[k][col] -= 2*M_PI;

         } else {

            /* Normal way of adjusting angles */
            d = cxang[k][col]-cxang[k-1][col];
      
            /* Adjust to avoid the -pi/pi discontinuity */
            if (d > M_PI+EPS)
               cxang[k][col] -= 2*M_PI;
            else if (d < -(M_PI-EPS) ) /* Forces clockwise rotation at 180 */
               cxang[k][col] += 2*M_PI;
         }
      }
   }
}

void interpolate_catmull_rom(flam3_genome cps[], double t, flam3_genome *result) {
   double t2 = t * t;
   double t3 = t2 * t;
   double cmc[4];

   cmc[0] = (2*t2 - t - t3) / 2;
   cmc[1] = (3*t3 - 5*t2 + 2) / 2;
   cmc[2] = (4*t2 - 3*t3 + t) / 2;
   cmc[3] = (t3 - t2) / 2;

   flam3_interpolate_n(result, 4, cps, cmc, 0);
}

static double smoother(double t) {
  return 3*t*t - 2*t*t*t;
}

static double get_stagger_coef(double t, double stagger_prc, int num_xforms, int this_xform) {

   /* max_stag is the spacing between xform start times if stagger_prc = 1.0 */
   double max_stag = (double)(num_xforms-1)/num_xforms;
   
   /* scale the spacing by stagger_prc */
   double stag_scaled = stagger_prc * max_stag;

   /* t ranges from 1 to 0 (the contribution of cp[0] to the blend) */
   /* the first line below makes the first xform interpolate first */
   /* the second line makes the last xform interpolate first */
   double st = stag_scaled * (num_xforms - 1 - this_xform) / (num_xforms-1);
//   double st = stag_scaled * (this_xform) / (num_xforms-1);
   double et = st + (1-stag_scaled);
   
//   printf("t=%f xf:%d st=%f et=%f : : %f\n",t,this_xform,st,et,smoother((t-st)/(1-stag_scaled)));
   
   if (t <= st)
      return (0);
   else if (t >= et)
      return (1);
   else
      return ( smoother((t-st)/(1-stag_scaled)) );

}
   

/* all cpi and result must be aligned (have the same number of xforms,
   and have final xform in the same slot) */
void flam3_interpolate_n(flam3_genome *result, int ncp,
          flam3_genome *cpi, double *c, double stagger) {
   int i, j, k, numstd;
   
	  /* HSV palette interpolation */
   
      for (i = 0; i < cpi[0].palette.count; i++) {
		 double4 s, t;
         int alpha1 = 1;

         s[0] = s[1] = s[2] = s[3] = 0.0;
         
         for (k = 0; k < ncp; k++) {
		    assert (cpi[k].palette.count == cpi[0].palette.count);
            t = rgb2hsv(cpi[k].palette.color[i]);
            for (j = 0; j < 3; j++)
               s[j] += c[k] * t[j];
            
            s[3] += c[k] * cpi[k].palette.color[i][3];
            if (cpi[k].palette.color[i][3] != 1.0)
               alpha1 = 0;
         }

         if (alpha1 == 1)
            s[3] = 1.0;
       
		 const double4 ret_color = hsv2rgb(s);
		 palette_add (&result->palette, (double4) { ret_color[0], ret_color[1],
				                    ret_color[2], s[3] });
       
      }

   result->symmetry = 0;
   result->palette_mode = cpi[0].palette_mode;

   result->interpolation_type = cpi[0].interpolation_type;
   INTERP(brightness);
   INTERP(contrast);
   INTERP(highlight_power);
   INTERP(gamma);
   INTERP(vibrancy);
   INTERP(hue_rotation);
   INTERI(width);
   INTERI(height);
   INTERP(center[0]);
   INTERP(center[1]);
   INTERP(pixels_per_unit);
   INTERP(zoom);
   INTERP(rotate);
   INTERP(gam_lin_thresh);
   
   /* Interpolate the chaos array */
   numstd = cpi[0].num_xforms - (cpi[0].final_xform_index >= 0);
   for (i=0;i<numstd;i++) {
      for (j=0;j<numstd;j++) {
         INTERP(chaos[i][j]);
         if (result->chaos[i][j]<0) result->chaos[i][j]=0;
         //chaos can be > 1
         //if (result->chaos[i][j]>1) result->chaos[i][j]=1.0;
      }
   }

   /* Interpolate each xform */
   for (i = 0; i < cpi[0].num_xforms; i++) {
   
      double csave[2];     
      double td;
      int all_id;
      int nx = cpi[0].num_xforms-(cpi[0].final_xform_index>=0);
      
      if (ncp==2 && stagger>0 && i!=cpi[0].final_xform_index) {
         csave[0] = c[0];
         csave[1] = c[1];
         c[0] = get_stagger_coef(csave[0],stagger,nx,i);
         c[1] = 1.0-c[0];
      }
      
      
      INTERP(xform[i].density);
      td = result->xform[i].density;
      result->xform[i].density = (td < 0.0) ? 0.0 : td;
      INTERP(xform[i].color);
      if (result->xform[i].color<0) result->xform[i].color=0;
      if (result->xform[i].color>1) result->xform[i].color=1;

      INTERP(xform[i].color_speed);
      if (result->xform[i].color_speed<0) result->xform[i].color_speed=0;
      if (result->xform[i].color_speed>1) result->xform[i].color_speed=1;
      
      INTERP(xform[i].opacity);      
      INTERP(xform[i].animate);
      INTERP(xform[i].blob_low);
      INTERP(xform[i].blob_high);
      INTERP(xform[i].blob_waves);
      INTERP(xform[i].pdj_ac[0]);
      INTERP(xform[i].pdj_bd[0]);
      INTERP(xform[i].pdj_ac[1]);
      INTERP(xform[i].pdj_bd[1]);
      INTERP(xform[i].fan2_x);
      INTERP(xform[i].fan2_y);
      INTERP(xform[i].rings2_val);
      INTERP(xform[i].perspective_angle);
      INTERP(xform[i].perspective_dist);
      INTERP(xform[i].julian_power);
      INTERP(xform[i].julian_dist);
      INTERP(xform[i].juliascope_power);
      INTERP(xform[i].juliascope_dist);
      INTERP(xform[i].radial_blur_angle);
      INTERP(xform[i].pie_slices);
      INTERP(xform[i].pie_rotation);
      INTERP(xform[i].pie_thickness);
      INTERP(xform[i].ngon_sides);
      INTERP(xform[i].ngon_power);
      INTERP(xform[i].ngon_circle);
      INTERP(xform[i].ngon_corners);
      INTERP(xform[i].curl_c1);
      INTERP(xform[i].curl_c2);
      INTERP(xform[i].rectangles_x);
      INTERP(xform[i].rectangles_y);
      INTERP(xform[i].disc2_rot);
      INTERP(xform[i].disc2_twist);
      INTERP(xform[i].super_shape_rnd);
      INTERP(xform[i].super_shape_m);
      INTERP(xform[i].super_shape_n1);
      INTERP(xform[i].super_shape_n2);
      INTERP(xform[i].super_shape_n3);
      INTERP(xform[i].super_shape_holes);
      INTERP(xform[i].flower_petals);
      INTERP(xform[i].flower_holes);
      INTERP(xform[i].conic_eccentricity);
      INTERP(xform[i].conic_holes);
      INTERP(xform[i].parabola_height);
      INTERP(xform[i].parabola_width);
      INTERP(xform[i].bent2_x);
      INTERP(xform[i].bent2_y);
      INTERP(xform[i].bipolar_shift);
      INTERP(xform[i].cell_size);
      INTERP(xform[i].cpow_r);
      INTERP(xform[i].cpow_i);
      INTERP(xform[i].cpow_power);
      INTERP(xform[i].curve_xamp);
      INTERP(xform[i].curve_yamp);
      INTERP(xform[i].curve_xlength);
      INTERP(xform[i].curve_ylength);
      INTERP(xform[i].escher_beta);
      INTERP(xform[i].lazysusan_x);
      INTERP(xform[i].lazysusan_y);
      INTERP(xform[i].lazysusan_twist);
      INTERP(xform[i].lazysusan_space);
      INTERP(xform[i].lazysusan_spin);
      INTERP(xform[i].modulus_x);
      INTERP(xform[i].modulus_y);
      INTERP(xform[i].oscope_separation);
      INTERP(xform[i].oscope_frequency);
      INTERP(xform[i].oscope_amplitude);
      INTERP(xform[i].oscope_damping);
      INTERP(xform[i].popcorn2_x);
      INTERP(xform[i].popcorn2_y);
      INTERP(xform[i].popcorn2_c);
      INTERP(xform[i].separation_x);
      INTERP(xform[i].separation_xinside);
      INTERP(xform[i].separation_y);
      INTERP(xform[i].separation_yinside);
      INTERP(xform[i].split_xsize);
      INTERP(xform[i].split_ysize);
      INTERP(xform[i].splits_x);
      INTERP(xform[i].splits_y);
      INTERP(xform[i].stripes_space);
      INTERP(xform[i].stripes_warp);
      INTERP(xform[i].wedge_angle);
      INTERP(xform[i].wedge_hole);
      INTERP(xform[i].wedge_count);
      INTERP(xform[i].wedge_swirl);
      INTERP(xform[i].wedge_julia_angle);
      INTERP(xform[i].wedge_julia_count);
      INTERP(xform[i].wedge_julia_power);
      INTERP(xform[i].wedge_julia_dist);
      INTERP(xform[i].wedge_sph_angle);
      INTERP(xform[i].wedge_sph_hole);
      INTERP(xform[i].wedge_sph_count);
      INTERP(xform[i].wedge_sph_swirl);
      INTERP(xform[i].whorl_inside);
      INTERP(xform[i].whorl_outside);
      INTERP(xform[i].waves2_scalex);
      INTERP(xform[i].waves2_scaley);
      INTERP(xform[i].waves2_freqx);
      INTERP(xform[i].waves2_freqy);
      INTERP(xform[i].auger_sym);
      INTERP(xform[i].auger_weight);
      INTERP(xform[i].auger_freq);
      INTERP(xform[i].auger_scale);
      INTERP(xform[i].flux_spread);
      INTERP(xform[i].mobius_re_a);
      INTERP(xform[i].mobius_im_a);
      INTERP(xform[i].mobius_re_b);
      INTERP(xform[i].mobius_im_b);
      INTERP(xform[i].mobius_re_c);
      INTERP(xform[i].mobius_im_c);
      INTERP(xform[i].mobius_re_d);
      INTERP(xform[i].mobius_im_d);
      INTERP(xform[i].asteria_alpha);
      INTERP(xform[i].bcollide_num);
      INTERP(xform[i].bcollide_a);

      for (j = 0; j < flam3_nvariations; j++)
         INTERP(xform[i].var[j]);

      if (flam3_inttype_log == cpi[0].interpolation_type) {
         double cxmag[4][2];  // XXX why only 4? should be ncp
         double cxang[4][2];
         double cxtrn[4][2];

         /* affine part */
         clear_matrix(result->xform[i].c);
         convert_linear_to_polar(cpi,ncp,i,0,cxang,cxmag,cxtrn);
         interp_and_convert_back(c, ncp, i, cxang, cxmag, cxtrn,result->xform[i].c);

         /* post part */
         all_id = 1;
         for (k=0; k<ncp; k++)
            all_id &= id_matrix(cpi[k].xform[i].post);
         
         clear_matrix(result->xform[i].post);
         if (all_id) {
            result->xform[i].post[0][0] = 1.0;
            result->xform[i].post[1][1] = 1.0;
         } else {
            convert_linear_to_polar(cpi,ncp,i,1,cxang,cxmag,cxtrn);
            interp_and_convert_back(c, ncp, i, cxang, cxmag, cxtrn,result->xform[i].post);
         }         
         
      } else {

         /* Interpolate c matrix & post */
         clear_matrix(result->xform[i].c);
         clear_matrix(result->xform[i].post);
         all_id = 1;
         for (k = 0; k < ncp; k++) {
            sum_matrix(c[k], cpi[k].xform[i].c, result->xform[i].c);
            sum_matrix(c[k], cpi[k].xform[i].post, result->xform[i].post);

            all_id &= id_matrix(cpi[k].xform[i].post);

         }
         if (all_id) {
            clear_matrix(result->xform[i].post);
            result->xform[i].post[0][0] = 1.0;
            result->xform[i].post[1][1] = 1.0;
         }
      }
      
      if (ncp==2 && stagger>0 && i!=cpi[0].final_xform_index) {
         c[0] = csave[0];
         c[1] = csave[1];
      }
      
   }
   
}

static void establish_asymmetric_refangles(flam3_genome *cp, int ncps) {

   int k, xfi, col;
   
   double cxang[4][2],d,c1[2];

   for (xfi=0; xfi<cp[0].num_xforms; xfi++) {
   
      /* Final xforms don't rotate regardless of their symmetry */
      if (cp[0].final_xform_enable==1 && xfi==cp[0].final_xform_index)
         continue;

      for (k=0; k<ncps;k++) {

         /* Establish the angle for each component */
         /* Should potentially functionalize */
         for (col=0;col<2;col++) {
         
            c1[0] = cp[k].xform[xfi].c[col][0];
            c1[1] = cp[k].xform[xfi].c[col][1];
            
            cxang[k][col] = atan2(c1[1],c1[0]);
         }
      }
      
      for (k=1; k<ncps; k++) {
     
         for (col=0;col<2;col++) {

            int sym0,sym1;
            int padsymflag;

            d = cxang[k][col]-cxang[k-1][col];

            /* Adjust to avoid the -pi/pi discontinuity */
            if (d > M_PI+EPS)
            cxang[k][col] -= 2*M_PI;
            else if (d < -(M_PI-EPS) )
            cxang[k][col] += 2*M_PI;

            /* If this is an asymmetric case, store the NON-symmetric angle    */
            /* Check them pairwise and store the reference angle in the second */
            /* to avoid overwriting if asymmetric on both sides                */
            padsymflag = 0;
         
            sym0 = (cp[k-1].xform[xfi].animate==0 || (cp[k-1].xform[xfi].padding==1 && padsymflag));
            sym1 = (cp[k].xform[xfi].animate==0 || (cp[k].xform[xfi].padding==1 && padsymflag));

            if ( sym1 && !sym0 )
               cp[k].xform[xfi].wind[col] = cxang[k-1][col] + 2*M_PI;
            else if ( sym0 && !sym1 )
               cp[k].xform[xfi].wind[col] = cxang[k][col] + 2*M_PI;

         }
      }
   }
}

void flam3_align(flam3_genome *dst, flam3_genome *src, int nsrc) {
   int i, tfx, tnx, max_nx = 0, max_fx = 0;
   int already_aligned=1;
   int xf,j;
   int ii,fnd;
   double normed;
   int usethisone;
   
   usethisone = (nsrc/2) - 1;
   
   max_nx = src[0].num_xforms - (src[0].final_xform_index >= 0);
   max_fx = src[0].final_xform_enable;
   
   for (i = 1; i < nsrc; i++) {
      tnx = src[i].num_xforms - (src[i].final_xform_index >= 0);
      if (max_nx != tnx) {
         already_aligned = 0;
         if (tnx > max_nx) max_nx = tnx;
      }
      
      tfx = src[i].final_xform_enable;
      if (max_fx != tfx) {
         already_aligned = 0;
         max_fx |= tfx;
      }
   }

   /* Pad the cps to equal xforms */
   for (i = 0; i < nsrc; i++) {
      flam3_copyx(&dst[i], &src[i], max_nx, max_fx);
   }

   /* Skip if this genome is compatibility mode */
   if (dst[usethisone].interpolation_type == flam3_inttype_compat ||
         dst[usethisone].interpolation_type == flam3_inttype_older)
      return;

      
   /* Check to see if there's a parametric variation present in one xform   */
   /* but not in an aligned xform.  If this is the case, use the parameters */
   /* from the xform with the variation as the defaults for the blank one.  */
   
   /* All genomes will have the same number of xforms at this point */
   /* num = max_nx + max_fx */
   for (i = 0; i<nsrc; i++) {


      for (xf = 0; xf<max_nx+max_fx; xf++) {
                  
         /* Loop over the variations to see which of them are set to 0 */
         /* Note that there are no parametric variations < 23 */
         for (j = 23; j < flam3_nvariations; j++) {
         
              if (dst[i].xform[xf].var[j]==0) {
            
                 if (i>0) {
                              
                    /* Check to see if the prior genome's xform is populated */
                    if (dst[i-1].xform[xf].var[j] != 0) {
                  
                       /* Copy the prior genome's parameters and continue */
                       flam3_copy_params(&(dst[i].xform[xf]), &(dst[i-1].xform[xf]), j);
                       continue;
                    }

                 } else if (i<nsrc-1) {

                    /* Check to see if the next genome's xform is populated */
                    if (dst[i+1].xform[xf].var[j] != 0) {
                  
                       /* Copy the next genome's parameters and continue */
                       flam3_copy_params(&(dst[i].xform[xf]), &(dst[i+1].xform[xf]), j);
                       continue;
                    }
                 }
              }
          } /* variations */

          if (dst[i].xform[xf].padding == 1 && !already_aligned) {
         
             /* This is a new xform.  Let's see if we can choose a better 'identity' xform. */
             /* Check the neighbors to see if any of these variations are used: */
             /* rings2, fan2, blob, perspective, julian, juliascope, ngon, curl, super_shape, split */
             /* If so, we can use a better starting point for these */
            
             /* Remove linear from the list */
             dst[i].xform[xf].var[0] = 0.0;
            
             /* Look through all of the 'companion' xforms to see if we get a match on any of these */
             fnd=0;

             /* Only do the next substitution for log interpolation */
             if ( (i==0 && dst[i].interpolation_type == flam3_inttype_log)
                  || (i>0 && dst[i-1].interpolation_type==flam3_inttype_log) ) {

             for (ii=-1; ii<=1; ii+=2) {

                /* Skip if out of bounds */
                if (i+ii<0 || i+ii>=nsrc)
                   continue;
                  
                /* Skip if this is also padding */
                if (dst[i+ii].xform[xf].padding==1)
                   continue;

                /* Spherical / Ngon (trumps all others due to holes)       */
                /* Interpolate these against a 180 degree rotated identity */
                /* with weight -1.                                         */
                /* Added JULIAN/JULIASCOPE to get rid of black wedges      */
                if (dst[i+ii].xform[xf].var[VAR_SPHERICAL]>0 ||
                      dst[i+ii].xform[xf].var[VAR_NGON]>0 || 
                      dst[i+ii].xform[xf].var[VAR_JULIAN]>0 || 
                      dst[i+ii].xform[xf].var[VAR_JULIASCOPE]>0 ||
                      dst[i+ii].xform[xf].var[VAR_POLAR]>0 ||
                      dst[i+ii].xform[xf].var[VAR_WEDGE_SPH]>0 ||
                      dst[i+ii].xform[xf].var[VAR_WEDGE_JULIA]>0) {
                 
                   dst[i].xform[xf].var[VAR_LINEAR] = -1.0;
                   /* Set the coefs appropriately */
                   dst[i].xform[xf].c[0][0] = -1.0;
                   dst[i].xform[xf].c[0][1] = 0.0;
                   dst[i].xform[xf].c[1][0] = 0.0;
                   dst[i].xform[xf].c[1][1] = -1.0;
                   dst[i].xform[xf].c[2][0] = 0.0;
                   dst[i].xform[xf].c[2][1] = 0.0;               
                   fnd=-1;
                }
             }

             }

             if (fnd==0) {

                for (ii=-1; ii<=1; ii+=2) {

                   /* Skip if out of bounds */
                   if (i+ii<0 || i+ii>=nsrc)
                      continue;
                     
                   /* Skip if also padding */
                   if (dst[i+ii].xform[xf].padding==1)
                      continue;

                   /* Rectangles */
                   if (dst[i+ii].xform[xf].var[VAR_RECTANGLES]>0) {
                      dst[i].xform[xf].var[VAR_RECTANGLES] = 1.0;
                      dst[i].xform[xf].rectangles_x = 0.0;
                      dst[i].xform[xf].rectangles_y = 0.0;
                      fnd++;
                   }

                   /* Rings 2 */
                   if (dst[i+ii].xform[xf].var[VAR_RINGS2]>0) {
                      dst[i].xform[xf].var[VAR_RINGS2] = 1.0;
                      dst[i].xform[xf].rings2_val = 0.0;
                      fnd++;
                   }
                  
                   /* Fan 2 */
                   if (dst[i+ii].xform[xf].var[VAR_FAN2]>0) {
                      dst[i].xform[xf].var[VAR_FAN2] = 1.0;
                      dst[i].xform[xf].fan2_x = 0.0;
                      dst[i].xform[xf].fan2_y = 0.0;
                      fnd++;
                   }
               
                   /* Blob */
                   if (dst[i+ii].xform[xf].var[VAR_BLOB]>0) {
                      dst[i].xform[xf].var[VAR_BLOB] = 1.0;
                      dst[i].xform[xf].blob_low = 1.0;
                      dst[i].xform[xf].blob_high = 1.0;
                      dst[i].xform[xf].blob_waves = 1.0;
                      fnd++;
                   }
               
                   /* Perspective */
                   if (dst[i+ii].xform[xf].var[VAR_PERSPECTIVE]>0) {
                      dst[i].xform[xf].var[VAR_PERSPECTIVE] = 1.0;
                      dst[i].xform[xf].perspective_angle = 0.0;
                      /* Keep the perspective distance as-is */
                      fnd++;
                   }
               
                   /* Curl */
                   if (dst[i+ii].xform[xf].var[VAR_CURL]>0) {
                      dst[i].xform[xf].var[VAR_CURL] = 1.0;
                      dst[i].xform[xf].curl_c1 = 0.0;
                      dst[i].xform[xf].curl_c2 = 0.0;
                      fnd++;
                   }

                   /* Super-Shape */
                   if (dst[i+ii].xform[xf].var[VAR_SUPER_SHAPE]>0) {
                      dst[i].xform[xf].var[VAR_SUPER_SHAPE] = 1.0;
                      /* Keep supershape_m the same */
                      dst[i].xform[xf].super_shape_n1 = 2.0;
                      dst[i].xform[xf].super_shape_n2 = 2.0;
                      dst[i].xform[xf].super_shape_n3 = 2.0;
                      dst[i].xform[xf].super_shape_rnd = 0.0;
                      dst[i].xform[xf].super_shape_holes = 0.0;
                      fnd++;
                   }
                }
             }

             /* If we didn't have any matches with those, */
             /* try the affine ones, fan and rings        */
             if (fnd==0) {
            
                for (ii=-1; ii<=1; ii+=2) {

                   /* Skip if out of bounds */
                   if (i+ii<0 || i+ii>=nsrc)
                      continue;                  

                   /* Skip if also a padding xform */
                   if (dst[i+ii].xform[xf].padding==1)
                      continue;
                     
                   /* Fan */
                   if (dst[i+ii].xform[xf].var[VAR_FAN]>0) {
                      dst[i].xform[xf].var[VAR_FAN] = 1.0;
                      fnd++;
                   }

                   /* Rings */
                   if (dst[i+ii].xform[xf].var[VAR_RINGS]>0) {
                      dst[i].xform[xf].var[VAR_RINGS] = 1.0;
                      fnd++;
                   }

                }
               
                if (fnd>0) {
                   /* Set the coefs appropriately */
                   dst[i].xform[xf].c[0][0] = 0.0;
                   dst[i].xform[xf].c[0][1] = 1.0;
                   dst[i].xform[xf].c[1][0] = 1.0;
                   dst[i].xform[xf].c[1][1] = 0.0;
                   dst[i].xform[xf].c[2][0] = 0.0;
                   dst[i].xform[xf].c[2][1] = 0.0;               
                }
             }
                                          
             /* If we still have no matches, switch back to linear */
             if (fnd==0)

                dst[i].xform[xf].var[VAR_LINEAR] = 1.0;

             else if (fnd>0) {

                /* Otherwise, go through and normalize the weights. */
                normed = 0.0;
                for (j = 0; j < flam3_nvariations; j++)
                   normed += dst[i].xform[xf].var[j];
                  
                for (j = 0; j < flam3_nvariations; j++)
                   dst[i].xform[xf].var[j] /= normed;

             }         
          }
       } /* xforms */
   } /* genomes */
                              
}

