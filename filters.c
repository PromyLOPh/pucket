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

#include <math.h>
#include <assert.h>

#include "filters.h"


/*
 * filter function definitions
 * from Graphics Gems III code
 * and ImageMagick resize.c
 */


double flam3_spatial_support[flam3_num_spatialfilters] = {

   1.5, /* gaussian */
   1.0, /* hermite */
   0.5, /* box */
   1.0, /* triangle */
   1.5, /* bell */
   2.0, /* b spline */
   2.0, /* mitchell */
   1.0, /* blackman */
   2.0, /* catrom */
   1.0, /* hanning */
   1.0, /* hamming */
   3.0, /* lanczos3 */
   2.0, /* lanczos2 */
   1.5  /* quadratic */
};

double flam3_hermite_filter(double t) {
   /* f(t) = 2|t|^3 - 3|t|^2 + 1, -1 <= t <= 1 */
   if(t < 0.0) t = -t;
   if(t < 1.0) return((2.0 * t - 3.0) * t * t + 1.0);
   return(0.0);
}

double flam3_box_filter(double t) {
   if((t > -0.5) && (t <= 0.5)) return(1.0);
   return(0.0);
}

double flam3_triangle_filter(double t) {
   if(t < 0.0) t = -t;
   if(t < 1.0) return(1.0 - t);
   return(0.0);
}

double flam3_bell_filter(double t) {
   /* box (*) box (*) box */
   if(t < 0) t = -t;
   if(t < .5) return(.75 - (t * t));
   if(t < 1.5) {
      t = (t - 1.5);
      return(.5 * (t * t));
   }
   return(0.0);
}

double flam3_b_spline_filter(double t) {

   /* box (*) box (*) box (*) box */
   double tt;

   if(t < 0) t = -t;
   if(t < 1) {
      tt = t * t;
      return((.5 * tt * t) - tt + (2.0 / 3.0));
   } else if(t < 2) {
      t = 2 - t;
      return((1.0 / 6.0) * (t * t * t));
   }
   return(0.0);
}

double flam3_sinc(double x) {
   x *= M_PI;
   if(x != 0) return(sin(x) / x);
   return(1.0);
}

double flam3_blackman_filter(double x) {
  return(0.42+0.5*cos(M_PI*x)+0.08*cos(2*M_PI*x));
}

double flam3_catrom_filter(double x) {
  if (x < -2.0)
    return(0.0);
  if (x < -1.0)
    return(0.5*(4.0+x*(8.0+x*(5.0+x))));
  if (x < 0.0)
    return(0.5*(2.0+x*x*(-5.0-3.0*x)));
  if (x < 1.0)
    return(0.5*(2.0+x*x*(-5.0+3.0*x)));
  if (x < 2.0)
    return(0.5*(4.0+x*(-8.0+x*(5.0-x))));
  return(0.0);
}

double flam3_mitchell_filter(double t) {
   double tt;

   tt = t * t;
   if(t < 0) t = -t;
   if(t < 1.0) {
      t = (((12.0 - 9.0 * flam3_mitchell_b - 6.0 * flam3_mitchell_c) * (t * tt))
         + ((-18.0 + 12.0 * flam3_mitchell_b + 6.0 * flam3_mitchell_c) * tt)
         + (6.0 - 2 * flam3_mitchell_b));
      return(t / 6.0);
   } else if(t < 2.0) {
      t = (((-1.0 * flam3_mitchell_b - 6.0 * flam3_mitchell_c) * (t * tt))
         + ((6.0 * flam3_mitchell_b + 30.0 * flam3_mitchell_c) * tt)
         + ((-12.0 * flam3_mitchell_b - 48.0 * flam3_mitchell_c) * t)
         + (8.0 * flam3_mitchell_b + 24 * flam3_mitchell_c));
      return(t / 6.0);
   }
   return(0.0);
}

double flam3_hanning_filter(double x) {
  return(0.5+0.5*cos(M_PI*x));
}

double flam3_hamming_filter(double x) {
  return(0.54+0.46*cos(M_PI*x));
}

double flam3_lanczos3_filter(double t) {
   if(t < 0) t = -t;
   if(t < 3.0) return(flam3_sinc(t) * flam3_sinc(t/3.0));
   return(0.0);
}

double flam3_lanczos2_filter(double t) {
   if(t < 0) t = -t;
   if(t < 2.0) return(flam3_sinc(t) * flam3_sinc(t/2.0));
   return(0.0);
}

double flam3_gaussian_filter(double x) {
  return(exp((-2.0*x*x))*sqrt(2.0/M_PI));
}

double flam3_quadratic_filter(double x) {
  if (x < -1.5)
    return(0.0);
  if (x < -0.5)
    return(0.5*(x+1.5)*(x+1.5));
  if (x < 0.5)
    return(0.75-x*x);
  if (x < 1.5)
    return(0.5*(x-1.5)*(x-1.5));
  return(0.0);
}

double flam3_spatial_filter(int knum, double x) {

   if (knum==0)
      return flam3_gaussian_filter(x);
   else if (knum==1)
      return flam3_hermite_filter(x);
   else if (knum==2)
      return flam3_box_filter(x);
   else if (knum==3)
      return flam3_triangle_filter(x);
   else if (knum==4)
      return flam3_bell_filter(x);
   else if (knum==5)
      return flam3_b_spline_filter(x);
   else if (knum==6)
      return flam3_mitchell_filter(x);
   else if (knum==7)
      return flam3_sinc(x)*flam3_blackman_filter(x);
   else if (knum==8)
      return flam3_catrom_filter(x);
   else if (knum==9)
      return flam3_sinc(x)*flam3_hanning_filter(x);
   else if (knum==10)
      return flam3_sinc(x)*flam3_hamming_filter(x);
   else if (knum==11)
      return flam3_lanczos3_filter(x)*flam3_sinc(x/3.0);   
   else if (knum==12)
      return flam3_lanczos2_filter(x)*flam3_sinc(x/2.0);
   else if (knum==13)
      return flam3_quadratic_filter(x);
   assert (0);
}

int normalize_vector(double *v, int n) {
   double t = 0.0;
   int i;
   for (i = 0; i < n; i++)
      t += v[i];
   if (0.0 == t) return 1;
   t = 1.0 / t;
   for (i = 0; i < n; i++)
      v[i] *= t;
   return 0;
}

double flam3_create_temporal_filter(int numsteps, int filter_type, double filter_exp, double filter_width,
                                    double **temporal_filter, double **temporal_deltas) {

   double maxfilt = 0.0;
   double sumfilt = 0.0;
   double slpx,halfsteps;
   double *deltas, *filter;
   
   int i;

   /* Allocate memory - this must be freed in the calling routine! */   
   deltas = (double *)malloc(numsteps*sizeof(double));
   filter = (double *)malloc(numsteps*sizeof(double));
   
   /* Deal with only one step */
   if (numsteps==1) {
      deltas[0] = 0;
      filter[0] = 1.0;
      *temporal_deltas = deltas;
      *temporal_filter = filter;
      return(1.0);
   }
      
   /* Define the temporal deltas */   
   for (i = 0; i < numsteps; i++)
      deltas[i] = ((double)i /(double)(numsteps - 1) - 0.5)*filter_width;
      
   /* Define the filter coefs */
   if (flam3_temporal_exp == filter_type) {

      for (i=0; i < numsteps; i++) {

         if (filter_exp>=0)
            slpx = ((double)i+1.0)/numsteps;
         else
            slpx = (double)(numsteps - i)/numsteps;

         /* Scale the color based on these values */
         filter[i] = pow(slpx,fabs(filter_exp));
         
         /* Keep the max */
         if (filter[i]>maxfilt)
            maxfilt = filter[i];
      }

   } else if (flam3_temporal_gaussian == filter_type) {

      halfsteps = numsteps/2.0;
      for (i=0; i < numsteps; i++) {
      
         /* Gaussian */
         filter[i] = flam3_spatial_filter(flam3_gaussian_kernel,
                           flam3_spatial_support[flam3_gaussian_kernel]*fabs(i - halfsteps)/halfsteps);
         /* Keep the max */
         if (filter[i]>maxfilt)
            maxfilt = filter[i];
      }
      
   } else { // (flam3_temporal_box)

      for (i=0; i < numsteps; i++)
         filter[i] = 1.0;
         
	   maxfilt = 1.0;
	   
   }

   /* Adjust the filter so that the max is 1.0, and */
   /* calculate the K2 scaling factor  */
   for (i=0;i<numsteps;i++) {
      filter[i] /= maxfilt;
      sumfilt += filter[i];
   }
         
   sumfilt /= numsteps;
   
   *temporal_deltas = deltas;
   *temporal_filter = filter;
   
   return(sumfilt);
}                                     
 
