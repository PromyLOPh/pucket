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

#include "variations.h"
#include "interpolation.h" 

#ifdef HAVE_AMDLIBM
#define REPLACE_WITH_AMDLIBM
#include <amdlibm.h>
#endif

#define badvalue(x) (((x)!=(x))||((x)>1e10)||((x)<-1e10))

typedef struct {
   double precalc_atan, precalc_sina;  /* Precalculated, if needed */
   double precalc_cosa, precalc_sqrt;
   double precalc_sumsq,precalc_atanyx;

   flam3_xform *xform; /* For the important values */

   /* Pointer to the RNG state */
   randctx *rc;

} flam3_iter_helper;

char *flam3_variation_names[1+flam3_nvariations] = {
  "linear",
  "sinusoidal",
  "spherical",
  "swirl",
  "horseshoe",
  "polar",
  "handkerchief",
  "heart",
  "disc",
  "spiral",
  "hyperbolic",
  "diamond",
  "ex",
  "julia",
  "bent",
  "waves",
  "fisheye",
  "popcorn",
  "exponential",
  "power",
  "cosine",
  "rings",
  "fan",
  "blob",
  "pdj",
  "fan2",
  "rings2",
  "eyefish",
  "bubble",
  "cylinder",
  "perspective",
  "noise",
  "julian",
  "juliascope",
  "blur",
  "gaussian_blur",
  "radial_blur",
  "pie",
  "ngon",
  "curl",
  "rectangles",
  "arch",
  "tangent",
  "square",
  "rays",
  "blade",
  "secant2",
  "twintrian",
  "cross",
  "disc2",
  "super_shape",
  "flower",
  "conic",
  "parabola",
  "bent2",
  "bipolar",
  "boarders",
  "butterfly",
  "cell",
  "cpow",
  "curve",
  "edisc",
  "elliptic",
  "escher",
  "foci",
  "lazysusan",
  "loonie",
  "pre_blur",
  "modulus",
  "oscilloscope",
  "polar2",
  "popcorn2",
  "scry",
  "separation",
  "split",
  "splits",
  "stripes",
  "wedge",
  "wedge_julia",
  "wedge_sph",
  "whorl",
  "waves2",
  "exp",
  "log",
  "sin",
  "cos",
  "tan",
  "sec",
  "csc",
  "cot",
  "sinh",
  "cosh",
  "tanh",
  "sech",
  "csch",
  "coth",
  "auger",
  "flux",
  "mobius",
  0
};

/*
 * VARIATION FUNCTIONS
 * must be of the form void (void *, double)
 */
static double2 var0_linear (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* linear */
   /* nx = tx;
      ny = ty;
      p[0] += v * nx;
      p[1] += v * ny; */

   return weight * in;
}

static double2 var1_sinusoidal (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* sinusoidal */
   /* nx = sin(tx);
      ny = sin(ty);
      p[0] += v * nx;
      p[1] += v * ny; */

   return weight * (double2) {sin(in[0]), sin(in[1])};
}

static double2 var2_spherical (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* spherical */
   /* double r2 = tx * tx + ty * ty + 1e-6;
      nx = tx / r2;
      ny = ty / r2;
      p[0] += v * nx;
      p[1] += v * ny; */

   const double r2 = weight / ( f->precalc_sumsq + EPS);

   return r2 * in;
}

static double2 var3_swirl (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* swirl */
   /* double r2 = tx * tx + ty * ty;    /k here is fun
      double c1 = sin(r2);
      double c2 = cos(r2);
      nx = c1 * tx - c2 * ty;
      ny = c2 * tx + c1 * ty;
      p[0] += v * nx;
      p[1] += v * ny; */

   double r2 = f->precalc_sumsq;
   double c1,c2;
   
   sincos(r2,&c1,&c2);
//   double c1 = sin(r2);
//   double c2 = cos(r2);
   const double2 n = (double2) {
                     c1 * in[0] - c2 * in[1],
                     c2 * in[0] + c1 * in[1],
					 };

   return weight * n;
}

static double2 var4_horseshoe (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* horseshoe */
   /* a = atan2(tx, ty);
      c1 = sin(a);
      c2 = cos(a);
      nx = c1 * tx - c2 * ty;
      ny = c2 * tx + c1 * ty;
      p[0] += v * nx;
      p[1] += v * ny;  */

   const double r = weight / (f->precalc_sqrt + EPS);

   return r * (double2) {
          (in[0] - in[1]) * (in[0] + in[1]),
          2.0 * in[0] * in[1],
		  };
}

static double2 var5_polar (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* polar */
   /* nx = atan2(tx, ty) / M_PI;
      ny = sqrt(tx * tx + ty * ty) - 1.0;
      p[0] += v * nx;
      p[1] += v * ny; */

   const double2 n = (double2) {
                     f->precalc_atan * M_1_PI,
                     f->precalc_sqrt - 1.0,
					 };

   return weight * n;
}

static double2 var6_handkerchief (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* folded handkerchief */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      p[0] += v * sin(a+r) * r;
      p[1] += v * cos(a-r) * r; */

   const double a = f->precalc_atan;
   const double r = f->precalc_sqrt;

   return weight * r * (double2) { sin(a+r), cos(a-r) };
}

static double2 var7_heart (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* heart */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      a *= r;
      p[0] += v * sin(a) * r;
      p[1] += v * cos(a) * -r; */

   const double a = f->precalc_sqrt * f->precalc_atan;
   double ca,sa;
   const double r = weight * f->precalc_sqrt;
   
   sincos(a,&sa,&ca);

   return (double2) { r * sa, (-r) * ca };
}

static double2 var8_disc (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* disc */
   /* nx = tx * M_PI;
      ny = ty * M_PI;
      a = atan2(nx, ny);
      r = sqrt(nx*nx + ny*ny);
      p[0] += v * sin(r) * a / M_PI;
      p[1] += v * cos(r) * a / M_PI; */

   const double a = f->precalc_atan * M_1_PI;
   const double r = M_PI * f->precalc_sqrt;
   double sr,cr;
   sincos(r,&sr,&cr);

   return weight * (double2) { sr, cr } * a;
}

static double2 var9_spiral (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* spiral */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty) + 1e-6;
      p[0] += v * (cos(a) + sin(r)) / r;
      p[1] += v * (sin(a) - cos(r)) / r; */

   const double r = f->precalc_sqrt + EPS;
   const double r1 = weight/r;
   double sr,cr;
   sincos(r,&sr,&cr);

   return r1 * (double2) { f->precalc_cosa + sr, f->precalc_sina - cr };
}

static double2 var10_hyperbolic (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* hyperbolic */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty) + 1e-6;
      p[0] += v * sin(a) / r;
      p[1] += v * cos(a) * r; */

   const double r = f->precalc_sqrt + EPS;

   return weight * (double2) { f->precalc_sina / r, f->precalc_cosa * r};
}

static double2 var11_diamond (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* diamond */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      p[0] += v * sin(a) * cos(r);
      p[1] += v * cos(a) * sin(r); */

   const double r = f->precalc_sqrt;
   double sr,cr;
   sincos(r,&sr,&cr);

   return weight * (double2) { f->precalc_sina * cr, f->precalc_cosa * sr };
}

static double2 var12_ex (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* ex */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      n0 = sin(a+r);
      n1 = cos(a-r);
      m0 = n0 * n0 * n0 * r;
      m1 = n1 * n1 * n1 * r;
      p[0] += v * (m0 + m1);
      p[1] += v * (m0 - m1); */

   const double a = f->precalc_atan;
   const double r = f->precalc_sqrt;

   const double2 n = (double2) { sin(a+r), cos(a-r) };

   const double2 m = n * n * n * r;

   return weight * (double2) { m[0] + m[1], m[0] - m[1] };
}

static double2 var13_julia (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* julia */
   /* a = atan2(tx, ty)/2.0;
      if (flam3_random_bit()) a += M_PI;
      r = pow(tx*tx + ty*ty, 0.25);
      nx = r * cos(a);
      ny = r * sin(a);
      p[0] += v * nx;
      p[1] += v * ny; */

   double r;
   double a = 0.5 * f->precalc_atan;
   double sa,ca;

   if (rand_bool(f->rc))
      a += M_PI;

   r = weight * sqrt(f->precalc_sqrt);
   
   sincos(a,&sa,&ca);

   return r * (double2) { ca, sa };
}

static double2 var14_bent (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* bent */
   /* nx = tx;
      ny = ty;
      if (nx < 0.0) nx = nx * 2.0;
      if (ny < 0.0) ny = ny / 2.0;
      p[0] += v * nx;
      p[1] += v * ny; */

   return weight * in * (double2) {
                          in[0] < 0.0 ? 2.0 : 1.0,
						  in[1] < 0.0 ? 0.5 : 1.0,
						  };
}

static double2 var15_waves (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* waves */
   /* dx = coef[2][0];
      dy = coef[2][1];
      nx = tx + coef[1][0]*sin(ty/((dx*dx)+EPS));
      ny = ty + coef[1][1]*sin(tx/((dy*dy)+EPS));
      p[0] += v * nx;
      p[1] += v * ny; */

   const double2 c1 = f->xform->c[1];
   const double2 inswap = (double2) { in[1], in[0] };
   const double2 a = inswap * f->xform->waves_d2;

   const double2 n = in + c1 * (double2) { sin(a[0]), sin(a[1]), };

   return weight * n;
}

static double2 var16_fisheye (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* fisheye */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      r = 2 * r / (r + 1);
      nx = r * cos(a);
      ny = r * sin(a);
      p[0] += v * nx;
      p[1] += v * ny; */

   double r = f->precalc_sqrt;

   r = 2.0 * weight / (r+1.0);

   /* XXX this seems to be wrong */
   const double2 tswap = (double2) { in[1], in[0] };

   return r * tswap;
}

static double2 var17_popcorn (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* popcorn */
   /* dx = tan(3*ty);
      dy = tan(3*tx);
      nx = tx + coef[2][0] * sin(dx);
      ny = ty + coef[2][1] * sin(dy);
      p[0] += v * nx;
      p[1] += v * ny; */

   const double dx = tan(3.0*in[1]);
   const double dy = tan(3.0*in[0]);

   const double2 n = in + f->xform->c[2] * (double2) { sin(dx), sin(dy) };

   return weight * n;
}

static double2 var18_exponential (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* exponential */
   /* dx = exp(tx-1.0);
      dy = M_PI * ty;
      nx = cos(dy) * dx;
      ny = sin(dy) * dx;
      p[0] += v * nx;
      p[1] += v * ny; */

   double dx = weight * exp(in[0] - 1.0);
   double dy = M_PI * in[1];
   double sdy,cdy;
   
   sincos(dy,&sdy,&cdy);
   
   return dx * (double2) { cdy, sdy };
}

static double2 var19_power (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* power */
   /* a = atan2(tx, ty);
      sa = sin(a);
      r = sqrt(tx*tx + ty*ty);
      r = pow(r, sa);
      nx = r * precalc_cosa;
      ny = r * sa;
      p[0] += v * nx;
      p[1] += v * ny; */

   const double r = weight * pow(f->precalc_sqrt, f->precalc_sina);

   return r * (double2) { f->precalc_cosa, f->precalc_sina };
}

static double2 var20_cosine (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* cosine */
   /* nx = cos(tx * M_PI) * cosh(ty);
      ny = -sin(tx * M_PI) * sinh(ty);
      p[0] += v * nx;
      p[1] += v * ny; */

   const double a = in[0] * M_PI;
   double sa,ca;
   
   sincos(a,&sa,&ca);

   const double2 n = (double2) { ca * cosh(in[1]), -sa * sinh(in[1]) };

   return weight * n;
}

static double2 var21_rings (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* rings */
   /* dx = coef[2][0];
      dx = dx * dx + EPS;
      r = sqrt(tx*tx + ty*ty);
      r = fmod(r + dx, 2*dx) - dx + r*(1-dx);
      a = atan2(tx, ty);
      nx = cos(a) * r;
      ny = sin(a) * r;
      p[0] += v * nx;
      p[1] += v * ny; */

   const double dx = f->xform->c[2][0] * f->xform->c[2][0] + EPS;
   double r = f->precalc_sqrt;
   r = weight * (fmod(r+dx, 2*dx) - dx + r * (1 - dx));

   return r * (double2) { f->precalc_cosa,  f->precalc_sina };
}

static double2 var22_fan (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* fan */
   /* dx = coef[2][0];
      dy = coef[2][1];
      dx = M_PI * (dx * dx + EPS);
      dx2 = dx/2;
      a = atan(tx,ty);
      r = sqrt(tx*tx + ty*ty);
      a += (fmod(a+dy, dx) > dx2) ? -dx2 : dx2;
      nx = cos(a) * r;
      ny = sin(a) * r;
      p[0] += v * nx;
      p[1] += v * ny; */

   const double dx = M_PI * (f->xform->c[2][0] * f->xform->c[2][0] + EPS);
   const double dy = f->xform->c[2][1];
   const double dx2 = 0.5 * dx;

   double a = f->precalc_atan;
   const double r = weight * f->precalc_sqrt;
   double sa,ca;

   a += (fmod(a+dy,dx) > dx2) ? -dx2 : dx2;
   sincos(a,&sa,&ca);

   return r * (double2) { ca, sa };
}

static double2 var23_blob (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* blob */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      r = r * (bloblow + (blobhigh-bloblow) * (0.5 + 0.5 * sin(blobwaves * a)));
      nx = sin(a) * r;
      ny = cos(a) * r;

      p[0] += v * nx;
      p[1] += v * ny; */

   double r = f->precalc_sqrt;
   const double a = f->precalc_atan;
   const double bdiff = f->xform->blob_high - f->xform->blob_low;

   r = r * (f->xform->blob_low +
            bdiff * (0.5 + 0.5 * sin(f->xform->blob_waves * a)));

   return weight * (double2) { f->precalc_sina, f->precalc_cosa } * r;
}

static double2 var24_pdj (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* pdj */
   /* nx1 = cos(pdjb * tx);
      nx2 = sin(pdjc * tx);
      ny1 = sin(pdja * ty);
      ny2 = cos(pdjd * ty);

      p[0] += v * (ny1 - nx1);
      p[1] += v * (nx2 - ny2); */

   const double2 a = (double2) {
                     sin(f->xform->pdj_a * in[1]),
					 sin(f->xform->pdj_c * in[0]) };
   const double2 b = (double2) {
                     cos(f->xform->pdj_b * in[0]),
					 cos(f->xform->pdj_d * in[1]) };

   return weight * (a - b);
}

static double2 var25_fan2 (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* fan2 */
   /* a = precalc_atan;
      r = precalc_sqrt;

      dy = fan2y;
      dx = M_PI * (fan2x * fan2x + EPS);
      dx2 = dx / 2.0;

      t = a + dy - dx * (int)((a + dy)/dx);

      if (t > dx2)
         a = a - dx2;
      else
         a = a + dx2;

      nx = sin(a) * r;
      ny = cos(a) * r;

      p[0] += v * nx;
      p[1] += v * ny; */

   double dy = f->xform->fan2_y;
   double dx = M_PI * (f->xform->fan2_x * f->xform->fan2_x + EPS);
   double dx2 = 0.5 * dx;
   double a = f->precalc_atan;
   double sa,ca;
   double r = weight * f->precalc_sqrt;

   double t = a + dy - dx * (int)((a + dy)/dx);

   if (t>dx2)
      a = a-dx2;
   else
      a = a+dx2;
      
   sincos(a,&sa,&ca);

   return r * (double2) { sa, ca };
}

static double2 var26_rings2 (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* rings2 */
   /* r = precalc_sqrt;
      dx = rings2val * rings2val + EPS;
      r += dx - 2.0*dx*(int)((r + dx)/(2.0 * dx)) - dx + r * (1.0-dx);
      nx = precalc_sina * r;
      ny = precalc_cosa * r;
      p[0] += v * nx;
      p[1] += v * ny; */

   double r = f->precalc_sqrt;
   double dx = f->xform->rings2_val * f->xform->rings2_val + EPS;

   r += -2.0*dx*(int)((r+dx)/(2.0*dx)) + r * (1.0-dx);

   return weight * (double2) { f->precalc_sina, f->precalc_cosa } * r;
}

static double2 var27_eyefish (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* eyefish */
   /* r = 2.0 * v / (precalc_sqrt + 1.0);
      p[0] += r*tx;
      p[1] += r*ty; */

   const double r = (weight * 2.0) / (f->precalc_sqrt + 1.0);

   return r * in;
}

static double2 var28_bubble (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* bubble */

   const double r = weight / (0.25 * (f->precalc_sumsq) + 1);

   return r * in;
}

static double2 var29_cylinder (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* cylinder (01/06) */

   return weight * (double2) { sin(in[0]), in[1] };
}

static double2 var30_perspective (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* perspective (01/06) */

   const double t = 1.0 / (f->xform->perspective_dist - in[1] * f->xform->persp_vsin);

   return weight * (double2) { f->xform->perspective_dist, f->xform->persp_vfcos } * in * t;
}

static double2 var31_noise (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* noise (03/06) */

   double tmpr, sinr, cosr, r;

   tmpr = rand_d01(f->rc) * 2 * M_PI;
   sincos(tmpr,&sinr,&cosr);

   r = weight * rand_d01(f->rc);

   return in * r * (double2) { cosr, sinr };
}

static double2 var32_juliaN_generic (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* juliaN (03/06) */

   int t_rnd = trunc((f->xform->julian_rN)*rand_d01(f->rc));
   
   double tmpr = (f->precalc_atanyx + 2 * M_PI * t_rnd) / f->xform->julian_power;

   double r = weight * pow(f->precalc_sumsq, f->xform->julian_cn);
   double sina, cosa;
   sincos(tmpr,&sina,&cosa);

   return r * (double2) { cosa, sina };
}

static double2 var33_juliaScope_generic (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* juliaScope (03/06) */

   int t_rnd = trunc((f->xform->juliascope_rN) * rand_d01(f->rc));

   double tmpr, r;
   double sina, cosa;

   if ((t_rnd & 1) == 0)
      tmpr = (2 * M_PI * t_rnd + f->precalc_atanyx) / f->xform->juliascope_power;
   else
      tmpr = (2 * M_PI * t_rnd - f->precalc_atanyx) / f->xform->juliascope_power;

   sincos(tmpr,&sina,&cosa);

   r = weight * pow(f->precalc_sumsq, f->xform->juliascope_cn);

   return r * (double2) { cosa, sina };
}

static double2 var34_blur (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* blur (03/06) */

   double tmpr, sinr, cosr, r;

   tmpr = rand_d01(f->rc) * 2 * M_PI;
   sincos(tmpr,&sinr,&cosr);

   r = weight * rand_d01(f->rc);

   return r * (double2) { cosr, sinr };
}

static double2 var35_gaussian (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* gaussian (09/06) */

   double ang, r, sina, cosa;

   ang = rand_d01(f->rc) * 2 * M_PI;
   sincos(ang,&sina,&cosa);

   r = weight * ( rand_d01(f->rc) + rand_d01(f->rc)
                   + rand_d01(f->rc) + rand_d01(f->rc) - 2.0 );

   return r * (double2) { cosa, sina };
}

static double2 var36_radial_blur (const double2 in, const flam3_iter_helper * const f, double weight) {
   /* radial blur (09/06) */
   /* removed random storage 6/07 */

   double rndG, ra, rz, tmpa, sa, ca;

   /* Get pseudo-gaussian */
   rndG = weight * (rand_d01(f->rc) + rand_d01(f->rc)
                   + rand_d01(f->rc) + rand_d01(f->rc) - 2.0);

   /* Calculate angle & zoom */
   ra = f->precalc_sqrt;
   tmpa = f->precalc_atanyx + f->xform->radialBlur_spinvar*rndG;
   sincos(tmpa,&sa,&ca);
   rz = f->xform->radialBlur_zoomvar * rndG - 1;

   return ra * (double2) { ca, sa } + rz * in;
}

static double2 var37_pie(const double2 in, const flam3_iter_helper * const f, double weight) {
   /* pie by Joel Faber (June 2006) */

   double a, r, sa, ca;
   int sl;

   sl = (int) (rand_d01(f->rc) * f->xform->pie_slices + 0.5);
   a = f->xform->pie_rotation +
       2.0 * M_PI * (sl + rand_d01(f->rc) * f->xform->pie_thickness) / f->xform->pie_slices;
   r = weight * rand_d01(f->rc);
   sincos(a,&sa,&ca);

   return r * (double2) { ca, sa };
}

static double2 var38_ngon(const double2 in, const flam3_iter_helper * const f, double weight) {
   /* ngon by Joel Faber (09/06) */

   double r_factor,theta,phi,b, amp;

   r_factor = pow(f->precalc_sumsq, f->xform->ngon_power/2.0);

   theta = f->precalc_atanyx;
   b = 2*M_PI/f->xform->ngon_sides;

   phi = theta - (b*floor(theta/b));
   if (phi > b/2)
      phi -= b;

   amp = f->xform->ngon_corners * (1.0 / (cos(phi) + EPS) - 1.0) + f->xform->ngon_circle;
   amp /= (r_factor + EPS);

   return weight * in * amp;
}

static double2 var39_curl(const double2 in, const flam3_iter_helper * const f, double weight)
{
    double re = 1.0 + f->xform->curl_c1 * in[0] + f->xform->curl_c2 * (in[0] * in[0] - in[1] * in[1]);
    double im = f->xform->curl_c1 * in[1] + 2.0 * f->xform->curl_c2 * in[0] * in[1];

    double r = weight / (re*re + im*im);

    double2 tswap = (double2) { in[1], -in[0] };

    return (in * re + tswap * im) * r;
}

static double2 var40_rectangles(const double2 in, const flam3_iter_helper * const f, double weight)
{
    return weight * (double2) {
	                f->xform->rectangles_x==0 ? in[0] : ((2 * floor(in[0] / f->xform->rectangles_x) + 1) * f->xform->rectangles_x - in[0]),
                    f->xform->rectangles_y==0 ? in[1] : ((2 * floor(in[1] / f->xform->rectangles_y) + 1) * f->xform->rectangles_y - in[1]),
					};
}

static double2 var41_arch(const double2 in, const flam3_iter_helper * const f, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Arch;
   var
     sinr, cosr: double;
   begin
     SinCos(random * vars[29]*pi, sinr, cosr);
     FPx := FPx + sinr*vars[29];
     FPy := FPy + sqr(sinr)/cosr*vars[29];
   end;
   */
   
   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   double ang = rand_d01(f->rc) * weight * M_PI;
   double sinr,cosr;
   sincos(ang,&sinr,&cosr);

   return weight * (double2) { sinr, (sinr*sinr)/cosr };

}

static double2 var42_tangent(const double2 in, const flam3_iter_helper * const f, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Tangent;
   begin
     FPx := FPx + vars[30] * (sin(FTx)/cos(FTy));
     FPy := FPy + vars[30] * (sin(FTy)/cos(FTy));
   end;
   */

   return weight * (double2) { sin(in[0])/cos(in[1]), tan(in[1]) };

}

static double2 var43_square(const double2 in, const flam3_iter_helper * const f, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.SquareBlur;
   begin
     FPx := FPx + vars[31] * (random - 0.5);
     FPy := FPy + vars[31] * (random - 0.5);
   end;
   */

   return weight * ((double2) {
                   rand_d01(f->rc),
                   rand_d01(f->rc),
				   } - 0.5);

}

static double2 var44_rays(const double2 in, const flam3_iter_helper * const f, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Rays;
   var
     r, sinr, cosr, tgr: double;
   begin
     SinCos(random * vars[32]*pi, sinr, cosr);
     r := vars[32] / (sqr(FTx) + sqr(FTy) + EPS);
     tgr := sinr/cosr;
     FPx := FPx + tgr * (cos(FTx)*vars[32]) * r;
     FPy := FPy + tgr * (sin(FTy)*vars[32]) * r;
   end;
   */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   double ang = weight * rand_d01(f->rc) * M_PI;
   double r = weight / (f->precalc_sumsq + EPS);
   double tanr = weight * tan(ang) * r;


   return tanr * (double2) { cos(in[0]), sin(in[1]) };

}

static double2 var45_blade(const double2 in, const flam3_iter_helper * const f, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Blade;
   var
     r, sinr, cosr: double;
   begin
     r := sqrt(sqr(FTx) + sqr(FTy))*vars[33];
     SinCos(r*random, sinr, cosr);
     FPx := FPx + vars[33] * FTx * (cosr + sinr);
     FPy := FPy + vars[33] * FTx * (cosr - sinr);
   end;
   */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   double r = rand_d01(f->rc) * weight * f->precalc_sqrt;
   double sinr,cosr;
   
   sincos(r,&sinr,&cosr);

   return weight * in[0] * (cosr + (double2) { sinr, -sinr });

}

static double2 var46_secant2(const double2 in, const flam3_iter_helper * const f, double weight)
{
   /* Intended as a 'fixed' version of secant */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   double r = weight * f->precalc_sqrt;
   double cr = cos(r);
   double icr = 1.0/cr;

   return weight * (double2) { in[0], cr<0 ? (icr + 1) : (icr - 1) };
}

static double2 var47_twintrian(const double2 in, const flam3_iter_helper * const f, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.TwinTrian;
   var
     r, diff, sinr, cosr: double;
   begin
     r := sqrt(sqr(FTx) + sqr(FTy))*vars[35];
     SinCos(r*random, sinr, cosr);
     diff := Math.Log10(sinr*sinr)+cosr;
     FPx := FPx + vars[35] * FTx * diff;
     FPy := FPy + vars[35] * FTx * (diff - (sinr*pi));
   end;
   */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   double r = rand_d01(f->rc) * weight * f->precalc_sqrt;
   double sinr,cosr,diff;
   
   sincos(r,&sinr,&cosr);
   diff = log10(sinr*sinr)+cosr;
   
   if (badvalue(diff))
      diff = -30.0;      

   
   return weight * in[0] * (diff - (double2) { 0.0, sinr*M_PI });

}

static double2 var48_cross(const double2 in, const flam3_iter_helper * const f, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Cross;
   var
     r: double;
   begin
     r := vars[36]*sqrt(1/(sqr(sqr(FTx)-sqr(FTy))+EPS));
     FPx := FPx + FTx * r;
     FPy := FPy + FTy * r;
   end;
   */

   double s = in[0]*in[0] - in[1]*in[1];
   double r = weight * sqrt(1.0 / (s*s+EPS));

   return r * in;
}

static double2 var49_disc2(const double2 in, const flam3_iter_helper * const f, double weight)
{
   /* Z+ variation Jan 07
   c := vvar/PI;
   k := rot*PI;
     sinadd := Sin(add);
     cosadd := Cos(add);
   cosadd := cosadd - 1;
   if (add > 2*PI) then begin
     cosadd := cosadd * (1 + add - 2*PI);
     sinadd := sinadd * (1 + add - 2*PI)
   end
   else if (add < -2*PI) then begin
     cosadd := cosadd * (1 + add + 2*PI);
     sinadd := sinadd * (1 + add + 2*PI)
   end
   end;
   procedure TVariationDisc2.CalcFunction;
   var
     r, sinr, cosr: extended;
   begin
     SinCos(k * (FTx^+FTy^), sinr, cosr);   //rot*PI
     r := c * arctan2(FTx^, FTy^); //vvar/PI
     FPx^ := FPx^ + (sinr + cosadd) * r;
     FPy^ := FPy^ + (cosr + sinadd) * r;
   */

   double r,t,sinr, cosr;

   t = f->xform->disc2_timespi * (in[0] + in[1]);
   sincos(t,&sinr,&cosr);
   r = weight * f->precalc_atan / M_PI;

   return r * (double2) {
              sinr + f->xform->disc2_cosadd,
			  cosr + f->xform->disc2_sinadd,
			  };

}

static double2 var50_supershape(const double2 in, const flam3_iter_helper * const f, double weight) {

   double theta;
   double t1,t2,r;
   double st,ct;
   double myrnd;

   theta = f->xform->super_shape_pm_4 * f->precalc_atanyx + M_PI_4;
   
   sincos(theta,&st,&ct);

   t1 = fabs(ct);
   t1 = pow(t1,f->xform->super_shape_n2);

   t2 = fabs(st);
   t2 = pow(t2,f->xform->super_shape_n3);
   
   myrnd = f->xform->super_shape_rnd;

   r = weight * ( (myrnd*rand_d01(f->rc) + (1.0-myrnd)*f->precalc_sqrt) - f->xform->super_shape_holes) 
      * pow(t1+t2,f->xform->super_shape_pneg1_n1) / f->precalc_sqrt;

   return r * in;
}

static double2 var51_flower(const double2 in, const flam3_iter_helper * const f, double weight) {
    /* cyberxaos, 4/2007 */
    /*   theta := arctan2(FTy^, FTx^);
         r := (random-holes)*cos(petals*theta);
         FPx^ := FPx^ + vvar*r*cos(theta);
         FPy^ := FPy^ + vvar*r*sin(theta);*/
 
    double theta = f->precalc_atanyx;
    double r = weight * (rand_d01(f->rc) - f->xform->flower_holes) * 
                    cos(f->xform->flower_petals*theta) / f->precalc_sqrt;

	return r * in;
}
    
static double2 var52_conic(const double2 in, const flam3_iter_helper * const f, double weight) {
    /* cyberxaos, 4/2007 */
    /*   theta := arctan2(FTy^, FTx^);
         r :=  (random - holes)*((eccentricity)/(1+eccentricity*cos(theta)));
         FPx^ := FPx^ + vvar*r*cos(theta);
         FPy^ := FPy^ + vvar*r*sin(theta); */
 
    double ct = in[0] / f->precalc_sqrt;
    double r = weight * (rand_d01(f->rc) - f->xform->conic_holes) * 
                    f->xform->conic_eccentricity / (1 + f->xform->conic_eccentricity*ct) / f->precalc_sqrt;

	return r * in;
}

static double2 var53_parabola(const double2 in, const flam3_iter_helper * const f, double weight) {
    /* cyberxaos, 4/2007 */
    /*   r := sqrt(sqr(FTx^) + sqr(FTy^));
         FPx^ := FPx^ + parabola_height*vvar*sin(r)*sin(r)*random;  
         FPy^ := FPy^ + parabola_width*vvar*cos(r)*random; */
 
    double r = f->precalc_sqrt;
    double sr,cr;
    
    sincos(r,&sr,&cr);
    
	return weight * (double2) {
                    f->xform->parabola_height * sr*sr * rand_d01(f->rc),
                    f->xform->parabola_width * cr * rand_d01(f->rc),
					};
    
}      

static double2 var54_bent2 (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Bent2 in the Apophysis Plugin Pack */   
   
   return weight * in * (double2) {
                          in[0] < 0.0 ? f->xform->bent2_x : 1.0,
						  in[1] < 0.0 ? f->xform->bent2_y : 1.0,
						  };
}

static double2 var55_bipolar (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Bipolar in the Apophysis Plugin Pack */   
   
   double x2y2 = f->precalc_sumsq;
   double t = x2y2+1;
   double x2 = 2*in[0];
   double ps = -M_PI_2 * f->xform->bipolar_shift;
   double y = 0.5 * atan2(2.0 * in[1], x2y2 - 1.0) + ps;
   
   if (y > M_PI_2)
       y = -M_PI_2 + fmod(y + M_PI_2, M_PI);
   else if (y < -M_PI_2)
       y = M_PI_2 - fmod(M_PI_2 - y, M_PI);

   return weight * M_2_PI * (double2) {
                   0.25 * log ( (t+x2) / (t-x2) ),
                   y,
				   };
}

static double2 var56_boarders (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Boarders in the Apophysis Plugin Pack */   
   
   double2 round = (double2) { rint(in[0]), rint(in[1]) };
   double2 offset = in - round;
    
   if (rand_d01(f->rc) >= 0.75) {
	  return weight*(offset*0.5 + round);
   } else {
      
      if (fabs(offset[0]) >= fabs(offset[1])) {
         
         if (offset[0] >= 0.0) {
			return weight*(offset*0.5 + round + 0.25 * (double2) { 1.0, offset[1] / offset[0] });
         } else {
		    return weight*(offset*0.5 + round - 0.25 * (double2) { 1.0, offset[1] / offset[0] });
         }
         
      } else {
         if (offset[1] >= 0.0) {
            return weight*(offset*0.5 + round + (double2) { offset[0]/offset[1]*0.25, 0.25 });
         } else {
		    return weight*(offset*0.5 + round - (double2) { offset[0]/offset[1]*0.25, 0.25 });
         }
      }
   }
}

static double2 var57_butterfly (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Butterfly in the Apophysis Plugin Pack */   
   
   /* wx is weight*4/sqrt(3*pi) */
   double wx = weight*1.3029400317411197908970256609023;
   
   double y2 = in[1]*2.0;
   double r = wx*sqrt(fabs(in[1] * in[0])/(EPS + in[0]*in[0] + y2*y2));
   
   return r * (double2) { in[0], y2 };
   
}

static double2 var58_cell (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Cell in the Apophysis Plugin Pack */   

   double inv_cell_size = 1.0/f->xform->cell_size;
    
   /* calculate input cell */
   double2 b = (double2) {
               floor(in[0]*inv_cell_size),
               floor(in[1]*inv_cell_size),
			   };

   /* Offset from cell origin */
   double2 a = in - b*f->xform->cell_size;
   
   /* interleave cells */
   if (b[1] >= 0.0) {
      if (b[0] >= 0.0) {
		 b *= 2.0;
      } else {
		 b = (double2) { -1.0, 1.0 } * (b*2.0 + (double2) { 1.0, 0.0 });
      }
   } else {
      if (b[0] >= 0.0) {
	     b = (double2) { 1.0, -1.0 } * (b*2.0 + (double2) { 0.0, 1.0 });
      } else {
	     b = -(2.0*b+1.0);
      }
   }
   
   return weight * (a + b*f->xform->cell_size) * (double2) { 1.0, -1.0 };
}

static double2 var59_cpow (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Cpow in the Apophysis Plugin Pack */   

   double a = f->precalc_atanyx;
   double lnr = 0.5 * log(f->precalc_sumsq);
   double va = 2.0 * M_PI / f->xform->cpow_power;
   double vc = f->xform->cpow_r / f->xform->cpow_power;
   double vd = f->xform->cpow_i / f->xform->cpow_power;
   double ang = vc*a + vd*lnr + va*floor(f->xform->cpow_power*rand_d01(f->rc));
   double sa,ca;
   
   double m = weight * exp(vc * lnr - vd * a);
   
   sincos(ang,&sa,&ca);
   
   return m * (double2) { ca, sa };
   
}

static double2 var60_curve (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Curve in the Apophysis Plugin Pack */   
   
   double pc_xlen = f->xform->curve_xlength*f->xform->curve_xlength;
   double pc_ylen = f->xform->curve_ylength*f->xform->curve_ylength;
   
   if (pc_xlen<1E-20) pc_xlen = 1E-20;
   
   if (pc_ylen<1E-20) pc_ylen = 1E-20;

   return weight * (in + (double2) {
                   f->xform->curve_xamp * exp(-in[1]*in[1]/pc_xlen),
                   f->xform->curve_yamp * exp(-in[0]*in[0]/pc_ylen),
				   });
      
}

static double2 var61_edisc (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Edisc in the Apophysis Plugin Pack */   
   
   double tmp = f->precalc_sumsq + 1.0;
   double tmp2 = 2.0 * in[0];
   double r1 = sqrt(tmp+tmp2);
   double r2 = sqrt(tmp-tmp2);
   double xmax = (r1+r2) * 0.5;
   double a1 = log(xmax + sqrt(xmax - 1.0));
   double a2 = -acos(in[0]/xmax);
   double w = weight / 11.57034632;
   double snv,csv,snhu,cshu;
   
   sincos(a1,&snv,&csv);
   
   snhu = sinh(a2);
   cshu = cosh(a2);
   
   if (in[1] > 0.0) snv = -snv;
   
   return w * (double2) { cshu * csv, snhu * snv };
   
}

static double2 var62_elliptic (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Elliptic in the Apophysis Plugin Pack */

   double tmp = f->precalc_sumsq + 1.0;
   double x2 = 2.0 * in[0];
   double xmax = 0.5 * (sqrt(tmp+x2) + sqrt(tmp-x2));
   double a = in[0] / xmax;
   double b = 1.0 - a*a;
   double ssx = xmax - 1.0;
   double w = weight / M_PI_2;
   
   if (b<0)
      b = 0;
   else
      b = sqrt(b);
      
   if (ssx<0)
      ssx = 0;
   else
      ssx = sqrt(ssx);
      
   return w * (double2) {
              atan2(a,b),
			  (in[1] > 0.0 ? 1.0 : -1.0) * log(xmax + ssx),
			  };
   
}

static double2 var63_escher (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Escher in the Apophysis Plugin Pack */
   
   double seb,ceb;
   double vc,vd;
   double m,n;
   double sn,cn;

   double a = f->precalc_atanyx;
   double lnr = 0.5 * log(f->precalc_sumsq);

   sincos(f->xform->escher_beta,&seb,&ceb);
   
   vc = 0.5 * (1.0 + ceb);
   vd = 0.5 * seb;

   m = weight * exp(vc*lnr - vd*a);
   n = vc*a + vd*lnr;
   
   sincos(n,&sn,&cn);
   
   return m * (double2) { cn, sn };

}

static double2 var64_foci (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Foci in the Apophysis Plugin Pack */

   double expx = exp(in[0]) * 0.5;
   double expnx = 0.25 / expx;
   double sn,cn,tmp;
   
   sincos(in[1],&sn,&cn);
   tmp = weight/(expx + expnx - cn);
   
   return tmp * (double2) { expx - expnx, sn };
      
}

static double2 var65_lazysusan (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Lazysusan in the Apophysis Plugin Pack */
   
   double2 b = in + (double2) { -f->xform->lazysusan_x, f->xform->lazysusan_y };
   double r = sqrt(b[0]*b[0] + b[1]*b[1]);
   double sina, cosa;
   
   if (r<weight) {
      double a = atan2(b[1],b[0]) + f->xform->lazysusan_spin +
                 f->xform->lazysusan_twist*(weight-r);
      sincos(a,&sina,&cosa);
      r = weight * r;
      
	  return r * (double2) { cosa, sina } + (double2) {
	                                        f->xform->lazysusan_x,
											-f->xform->lazysusan_y,
											};
   } else {
      
      r = weight * (1.0 + f->xform->lazysusan_space / r);
      
	  return r * b + (double2) {
                     f->xform->lazysusan_x,
					 -f->xform->lazysusan_y,
					 };
   
   }
      
}

static double2 var66_loonie (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Loonie in the Apophysis Plugin Pack */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */
   
   double r2 = f->precalc_sumsq;
   double w2 = weight*weight;
   
   if (r2 < w2) {
      double r = weight * sqrt(w2/r2 - 1.0);
	  return r * in;
   } else {
	  return weight * in;
   }
         
}

static double2 var67_pre_blur (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* pre-xform: PreBlur (Apo 2.08) */
   
   /* Get pseudo-gaussian */
   double rndG = weight * (rand_d01(f->rc) + rand_d01(f->rc)
                   + rand_d01(f->rc) + rand_d01(f->rc) - 2.0);
   double rndA = rand_d01(f->rc) * 2.0 * M_PI;
   double sinA,cosA;
   
   sincos(rndA,&sinA,&cosA);
   
   /* Note: original coordinate changed */
   return rndG * (double2) { cosA, sinA };

}

static double2 var68_modulus (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* Modulus in the Apophysis Plugin Pack */
   
   double xr = 2*f->xform->modulus_x;
   double yr = 2*f->xform->modulus_y;
   double a, b;
   
   if (in[0] > f->xform->modulus_x)
      a = (-f->xform->modulus_x + fmod(in[0] + f->xform->modulus_x, xr));
   else if (in[0] < -f->xform->modulus_x)
      a = ( f->xform->modulus_x - fmod(f->xform->modulus_x - in[0], xr));
   else
      a = in[0];
      
   if (in[1] > f->xform->modulus_y)
      b = (-f->xform->modulus_y + fmod(in[1] + f->xform->modulus_y, yr));
   else if (in[1] < -f->xform->modulus_y)
      b = ( f->xform->modulus_y - fmod(f->xform->modulus_y - in[1], yr));
   else
      b = in[1];
         
   return weight * (double2) { a, b };
}

static double2 var69_oscope (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* oscilloscope from the apophysis plugin pack */
   
   double tpf = 2 * M_PI * f->xform->oscope_frequency;
   double t;
   
   if (f->xform->oscope_damping == 0.0)
      t = f->xform->oscope_amplitude * cos(tpf*in[0]) + f->xform->oscope_separation;
   else {
      t = f->xform->oscope_amplitude * exp(-fabs(in[0])*f->xform->oscope_damping)
          * cos(tpf*in[0]) + f->xform->oscope_separation;
   }
   
   return (double2) { 1.0, fabs(in[1]) <= t ? -1.0 : 1.0 } * weight * in;
}

static double2 var70_polar2 (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* polar2 from the apophysis plugin pack */
   
   double p2v = weight / M_PI;
   
   return (double2) { p2v * f->precalc_atan, p2v/2.0 * log(f->precalc_sumsq) };
}

static double2 var71_popcorn2 (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* popcorn2 from the apophysis plugin pack */
   
   return weight * (in + (double2) {
                   f->xform->popcorn2_x * sin(tan(in[1]*f->xform->popcorn2_c)),
				   f->xform->popcorn2_y * sin(tan(in[0]*f->xform->popcorn2_c)),
				   });

}

static double2 var72_scry (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* scry from the apophysis plugin pack */
   /* note that scry does not multiply by weight, but as the */
   /* values still approach 0 as the weight approaches 0, it */
   /* should be ok                                           */ 

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */
   
   double t = f->precalc_sumsq;
   double r = 1.0 / (f->precalc_sqrt * (t + 1.0/(weight+EPS)));
   
   return r * in;

}

static double2 var73_separation (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* separation from the apophysis plugin pack */

   const double sx2 = f->xform->separation_x * f->xform->separation_x;
   const double sy2 = f->xform->separation_y * f->xform->separation_y;
   
   const double2 a = (double2) {
                     sqrt(in[0]*in[0] + sx2),
                     sqrt(in[1]*in[1] + sy2),
					 };
   const double2 b = (double2) {
                     in[0]*f->xform->separation_xinside,
                     in[1]*f->xform->separation_yinside,
					 };
   const double2 sign = (double2) {
                         in[0] > 0.0 ? 1.0 : -1.0,
						 in[1] > 0.0 ? 1.0 : -1.0,
						 };
   const double2 bsign = (double2) {
                         in[0] > 0.0 ? -1.0 : 1.0,
						 in[1] > 0.0 ? -1.0 : 1.0,
                         };

   return sign * (weight * (a + bsign * b));
}

static double2 var74_split (const double2 in, const flam3_iter_helper * const f, double weight) {
   
   /* Split from apo plugins pack */

   return weight * (double2) {
                   cos(in[1]*f->xform->split_ysize*M_PI) >= 0 ? 1.0 : -1.0,
                   cos(in[0]*f->xform->split_xsize*M_PI) >= 0 ? 1.0 : -1.0,
				   } * in;
}

static double2 var75_splits (const double2 in, const flam3_iter_helper * const f, double weight) {
   
   /* Splits from apo plugins pack */

   return weight * (in + (double2) {
                           in[0] >= 0 ? f->xform->splits_x : -f->xform->splits_x,
                           in[1] >= 0 ? f->xform->splits_y : -f->xform->splits_y,
						   });
}

static double2 var76_stripes (const double2 in, const flam3_iter_helper * const f, double weight) {
   
   /* Stripes from apo plugins pack */

   double roundx,offsetx;
   
   roundx = floor(in[0] + 0.5);
   offsetx = in[0] - roundx;

   return weight * (double2) {
                   offsetx*(1.0-f->xform->stripes_space)+roundx,
                   in[1] + offsetx*offsetx*f->xform->stripes_warp,
				   };

}

static double2 var77_wedge (const double2 in, const flam3_iter_helper * const f, double weight) {
   
   /* Wedge from apo plugins pack */

   double r = f->precalc_sqrt;
   double a = f->precalc_atanyx + f->xform->wedge_swirl * r;
   double c = floor( (f->xform->wedge_count * a + M_PI)*M_1_PI*0.5);
   
   double comp_fac = 1 - f->xform->wedge_angle*f->xform->wedge_count*M_1_PI*0.5;
   double sa, ca;
   
   a = a * comp_fac + c * f->xform->wedge_angle;
   
   sincos(a,&sa,&ca);

   r = weight * (r + f->xform->wedge_hole);
   
   return r * (double2) { ca, sa };

}

static double2 var78_wedge_julia (const double2 in, const flam3_iter_helper * const f, double weight) {

   /* wedge_julia from apo plugin pack */

   double r = weight * pow(f->precalc_sumsq, f->xform->wedgeJulia_cn);
   int t_rnd = (int)((f->xform->wedgeJulia_rN)*rand_d01(f->rc));
   double a = (f->precalc_atanyx + 2 * M_PI * t_rnd) / f->xform->wedge_julia_power;
   double c = floor( (f->xform->wedge_julia_count * a + M_PI)*M_1_PI*0.5 );
   double sa,ca;
   
   a = a * f->xform->wedgeJulia_cf + c * f->xform->wedge_julia_angle;
   
   sincos(a,&sa,&ca);

   return r * (double2) { ca, sa };
}

static double2 var79_wedge_sph (const double2 in, const flam3_iter_helper * const f, double weight) {
   
   /* Wedge_sph from apo plugins pack */

   double r = 1.0/(f->precalc_sqrt+EPS);
   double a = f->precalc_atanyx + f->xform->wedge_sph_swirl * r;
   double c = floor( (f->xform->wedge_sph_count * a + M_PI)*M_1_PI*0.5);
   
   double comp_fac = 1 - f->xform->wedge_sph_angle*f->xform->wedge_sph_count*M_1_PI*0.5;
   double sa, ca;
   
   a = a * comp_fac + c * f->xform->wedge_sph_angle;

   sincos(a,&sa,&ca);   
   r = weight * (r + f->xform->wedge_sph_hole);
   
   return r * (double2) { ca, sa };

}

static double2 var80_whorl (const double2 in, const flam3_iter_helper * const f, double weight) {
   
   /* whorl from apo plugins pack */
   
   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   double r = f->precalc_sqrt;
   double a,sa,ca;

   if (r<weight)
      a = f->precalc_atanyx + f->xform->whorl_inside/(weight-r);
   else
      a = f->precalc_atanyx + f->xform->whorl_outside/(weight-r);
   
   sincos(a,&sa,&ca);
   
   return weight * r * (double2) { ca, sa };

}

static double2 var81_waves2 (const double2 in, const flam3_iter_helper * const f, double weight) {
   
   /* waves2 from Joel F */
   
   return weight * (in + (double2) {
                   f->xform->waves2_scalex*sin(in[1] * f->xform->waves2_freqx),
                   f->xform->waves2_scaley*sin(in[0] * f->xform->waves2_freqy),
				   });

}

/* complex vars by cothe */
/* exp log sin cos tan sec csc cot sinh cosh tanh sech csch coth */

static double2 var82_exp (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Exponential EXP
   double expe = exp(in[0]);
   double expcos,expsin;
   sincos(in[1],&expsin,&expcos);
   return weight * expe * (double2) { expcos, expsin };
}
        
static double2 var83_log (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Natural Logarithm LOG
   // needs precalc_atanyx and precalc_sumsq
   return weight * (double2) { 0.5 * log(f->precalc_sumsq), f->precalc_atanyx };
}

static double2 var84_sin (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Sine SIN
   double sinsin,sinacos,sinsinh,sincosh;
   sincos(in[0],&sinsin,&sinacos);
   sinsinh = sinh(in[1]);
   sincosh = cosh(in[1]);
   return weight * (double2) { sinsin * sincosh, sinacos * sinsinh };
}

static double2 var85_cos (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Cosine COS
   double cossin,coscos,cossinh,coscosh;
   sincos(in[0],&cossin,&coscos);
   cossinh = sinh(in[1]);
   coscosh = cosh(in[1]);
   return weight * (double2) { coscos * coscosh, -1.0 * cossin * cossinh };
}

static double2 var86_tan (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Tangent TAN
   double tansin,tancos,tansinh,tancosh;
   double tanden;
   sincos(2*in[0],&tansin,&tancos);
   tansinh = sinh(2.0*in[1]);
   tancosh = cosh(2.0*in[1]);
   tanden = 1.0/(tancos + tancosh);
   return weight * tanden * (double2) { tansin, tansinh };
}

static double2 var87_sec (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Secant SEC
   double secsin,seccos,secsinh,seccosh;
   double secden;
   sincos(in[0],&secsin,&seccos);
   secsinh = sinh(in[1]);
   seccosh = cosh(in[1]);
   secden = 2.0/(cos(2*in[0]) + cosh(2*in[1]));
   return weight * secden * (double2) { seccos * seccosh, secsin * secsinh };
}

static double2 var88_csc (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Cosecant CSC
   double cscsin,csccos,cscsinh,csccosh;
   double cscden;
   sincos(in[0],&cscsin,&csccos);
   cscsinh = sinh(in[1]);
   csccosh = cosh(in[1]);
   cscden = 2.0/(cosh(2.0*in[1]) - cos(2.0*in[0]));
   return weight * cscden * (double2) { cscsin * csccosh, -1.0 * csccos * cscsinh };
}

static double2 var89_cot (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Cotangent COT
   double cotsin,cotcos,cotsinh,cotcosh;
   double cotden;
   sincos(2.0*in[0],&cotsin,&cotcos);
   cotsinh = sinh(2.0*in[1]);
   cotcosh = cosh(2.0*in[1]);
   cotden = 1.0/(cotcosh - cotcos);
   return weight * cotden * (double2) { cotsin, -1.0 * cotsinh };
}

static double2 var90_sinh (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Hyperbolic Sine SINH
   double sinhsin,sinhcos,sinhsinh,sinhcosh;
   sincos(in[1],&sinhsin,&sinhcos);
   sinhsinh = sinh(in[0]);
   sinhcosh = cosh(in[0]);
   return weight * (double2) { sinhsinh * sinhcos, sinhcosh * sinhsin };
}

static double2 var91_cosh (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Hyperbolic Cosine COSH
   double coshsin,coshcos,coshsinh,coshcosh;
   sincos(in[1],&coshsin,&coshcos);
   coshsinh = sinh(in[0]);
   coshcosh = cosh(in[0]);
   return weight * (double2) { coshcosh * coshcos, coshsinh * coshsin };
}

static double2 var92_tanh (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Hyperbolic Tangent TANH
   double tanhsin,tanhcos,tanhsinh,tanhcosh;
   double tanhden;
   sincos(2.0*in[1],&tanhsin,&tanhcos);
   tanhsinh = sinh(2.0*in[0]);
   tanhcosh = cosh(2.0*in[0]);
   tanhden = 1.0/(tanhcos + tanhcosh);
   return weight * tanhden * (double2) { tanhsinh, tanhsin };
}

static double2 var93_sech (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Hyperbolic Secant SECH
   double sechsin,sechcos,sechsinh,sechcosh;
   double sechden;
   sincos(in[1],&sechsin,&sechcos);
   sechsinh = sinh(in[0]);
   sechcosh = cosh(in[0]);
   sechden = 2.0/(cos(2.0*in[1]) + cosh(2.0*in[0]));
   return weight * sechden * (double2) { sechcos * sechcosh, -1.0 * sechsin * sechsinh };
}

static double2 var94_csch (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Hyperbolic Cosecant CSCH
   double cschsin,cschcos,cschsinh,cschcosh;
   double cschden;
   sincos(in[1],&cschsin,&cschcos);
   cschsinh = sinh(in[0]);
   cschcosh = cosh(in[0]);
   cschden = 2.0/(cosh(2.0*in[0]) - cos(2.0*in[1]));
   return weight * cschden * (double2) { cschsinh * cschcos, -1.0 * cschcosh * cschsin };
}

static double2 var95_coth (const double2 in, const flam3_iter_helper * const f, double weight) {
   //Hyperbolic Cotangent COTH
   double cothsin,cothcos,cothsinh,cothcosh;
   double cothden;
   sincos(2.0*in[1],&cothsin,&cothcos);
   cothsinh = sinh(2.0*in[0]);
   cothcosh = cosh(2.0*in[0]);
   cothden = 1.0/(cothcosh - cothcos);
   return weight * cothden * (double2) { cothsinh, cothsin };
}

static double2 var96_auger (const double2 in, const flam3_iter_helper * const f, double weight) {

    // Auger, by Xyrus01
    double s = sin(f->xform->auger_freq * in[0]);
    double t = sin(f->xform->auger_freq * in[1]);
    double dy = in[1] + f->xform->auger_weight*(f->xform->auger_scale*s/2.0 + fabs(in[1])*s);
    double dx = in[0] + f->xform->auger_weight*(f->xform->auger_scale*t/2.0 + fabs(in[0])*t);

	return weight * (double2) { in[0] + f->xform->auger_sym*(dx-in[0]), dy };
}

static double2 var97_flux (const double2 in, const flam3_iter_helper * const f, double weight) {

    // Flux, by meckie
    double xpw = in[0] + weight;
    double xmw = in[0] - weight;
    double avgr = weight * (2 + f->xform->flux_spread) * sqrt( sqrt(in[1]*in[1] + xpw*xpw) / sqrt(in[1]*in[1] + xmw*xmw));
    double avga = ( atan2(in[1], xmw) - atan2(in[1],xpw) ) * 0.5;

    double s = sin(avga);
    double c = cos(avga);

    return avgr * (double2) { cos(avga), sin(avga) };
}

static double2 var98_mobius (const double2 in, const flam3_iter_helper * const f, double weight) {

    // Mobius, by eralex
    double re_u, im_u, re_v, im_v, rad_v;

    re_u = f->xform->mobius_re_a * in[0] - f->xform->mobius_im_a * in[1] + f->xform->mobius_re_b;
    im_u = f->xform->mobius_re_a * in[1] + f->xform->mobius_im_a * in[0] + f->xform->mobius_im_b;
    re_v = f->xform->mobius_re_c * in[0] - f->xform->mobius_im_c * in[1] + f->xform->mobius_re_d;
    im_v = f->xform->mobius_re_c * in[1] + f->xform->mobius_im_c * in[0] + f->xform->mobius_im_d;

    rad_v = weight / (re_v*re_v + im_v*im_v);

	return rad_v * (double2) { (re_u*re_v + im_u*im_v), (im_u*re_v - re_u*im_v) };
}
    

/* Precalc functions */

static void perspective_precalc(flam3_xform *xf) {
   double ang = xf->perspective_angle * M_PI / 2.0;
   xf->persp_vsin = sin(ang);
   xf->persp_vfcos = xf->perspective_dist * cos(ang);
}

static void juliaN_precalc(flam3_xform *xf) {
   xf->julian_rN = fabs(xf->julian_power);
   xf->julian_cn = xf->julian_dist / (double)xf->julian_power / 2.0;
}

static void wedgeJulia_precalc(flam3_xform *xf) {
   xf->wedgeJulia_cf = 1.0 - xf->wedge_julia_angle * xf->wedge_julia_count * M_1_PI * 0.5;
   xf->wedgeJulia_rN = fabs(xf->wedge_julia_power);
   xf->wedgeJulia_cn = xf->wedge_julia_dist / xf->wedge_julia_power / 2.0;
}

static void juliaScope_precalc(flam3_xform *xf) {
   xf->juliascope_rN = fabs(xf->juliascope_power);
   xf->juliascope_cn = xf->juliascope_dist / (double)xf->juliascope_power / 2.0;
}

static void radial_blur_precalc(flam3_xform *xf) {
   sincos(xf->radial_blur_angle * M_PI / 2.0,
             &xf->radialBlur_spinvar, &xf->radialBlur_zoomvar);
}

static void waves_precalc(flam3_xform *xf) {
   const double2 d = xf->c[2];

   xf->waves_d2 = 1.0/(d * d + EPS);
}

static void disc2_precalc(flam3_xform *xf) {
   double add = xf->disc2_twist;
   double k;

   xf->disc2_timespi = xf->disc2_rot * M_PI;

   sincos(add,&xf->disc2_sinadd,&xf->disc2_cosadd);
   xf->disc2_cosadd -= 1;

   if (add > 2 * M_PI) {
      k = (1 + add - 2*M_PI);
      xf->disc2_cosadd *= k;
      xf->disc2_sinadd *= k;
   }

   if (add < -2 * M_PI) {
      k = (1 + add + 2*M_PI);
      xf->disc2_cosadd *= k;
      xf->disc2_sinadd *= k;
   }
}

static void supershape_precalc(flam3_xform *xf) {
   xf->super_shape_pm_4 = xf->super_shape_m / 4.0;
   xf->super_shape_pneg1_n1 = -1.0 / xf->super_shape_n1;
}

void xform_precalc(flam3_genome *cp, int xi) {

   perspective_precalc(&(cp->xform[xi]));
   juliaN_precalc(&(cp->xform[xi]));
   juliaScope_precalc(&(cp->xform[xi]));
   radial_blur_precalc(&(cp->xform[xi]));
   waves_precalc(&(cp->xform[xi]));
   disc2_precalc(&(cp->xform[xi]));
   supershape_precalc(&(cp->xform[xi]));
   wedgeJulia_precalc(&(cp->xform[xi]));   
}   

int prepare_precalc_flags(flam3_genome *cp) {

   double d;
   int i,j,totnum;

   /* Loop over valid xforms */
   for (i = 0; i < cp->num_xforms; i++) {
      d = cp->xform[i].density;
      if (d < 0.0) {
         fprintf(stderr, "xform %d weight must be non-negative, not %g.\n",i,d);
         return(1);
      }

      if (i != cp->final_xform_index && d == 0.0)
         continue;

      totnum = 0;

      cp->xform[i].vis_adjusted = adjust_percentage(cp->xform[i].opacity);

      cp->xform[i].precalc_angles_flag=0;
      cp->xform[i].precalc_atan_xy_flag=0;
      cp->xform[i].precalc_atan_yx_flag=0;
      cp->xform[i].has_preblur=0;
      cp->xform[i].has_post = !(id_matrix(cp->xform[i].post));


      for (j = 0; j < flam3_nvariations; j++) {

         if (cp->xform[i].var[j]!=0) {

            cp->xform[i].varFunc[totnum] = j;
            cp->xform[i].active_var_weights[totnum] = cp->xform[i].var[j];

            if (j==VAR_POLAR) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_HANDKERCHIEF) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_HEART) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_DISC) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_SPIRAL) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_HYPERBOLIC) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_DIAMOND) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_EX) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_JULIA) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_POWER) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_RINGS) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_FAN) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_BLOB) {
               cp->xform[i].precalc_atan_xy_flag=1;
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_FAN2) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_RINGS2) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_JULIAN) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_JULIASCOPE) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_RADIAL_BLUR) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_NGON) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_DISC2) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_SUPER_SHAPE) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_FLOWER) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_CONIC) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_CPOW) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_ESCHER) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_PRE_BLUR) {
               cp->xform[i].has_preblur=cp->xform[i].var[j];
            } else if (j==VAR_POLAR2) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_WEDGE) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_WEDGE_JULIA) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_WEDGE_SPH) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_WHORL) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_LOG) {
               cp->xform[i].precalc_atan_yx_flag=1;
            }
            
            totnum++;
         }
      }

      cp->xform[i].num_active_vars = totnum;

   }
   
   return(0);
}

/*	Apply affine coordinate transformation
 */
static double2 apply_affine (const double2 in, const double2 matrix[3]) {
	return matrix[0] * in[0] + matrix[1] * in[1] + matrix[2];
}

static double sum(const double2 in) {
	return in[0] + in[1];
}

int apply_xform(flam3_genome *cp, int fn, const double4 p, double4 *q_ret, randctx *rc)
{
   flam3_iter_helper f;
   int var_n;
   double next_color,s,s1;
   double weight;

   f.rc = rc;

   s1 = cp->xform[fn].color_speed;

   const double2 q23 = (double2) {
         s1 * cp->xform[fn].color + (1.0-s1) * p[2],
         cp->xform[fn].vis_adjusted,
		 };

   //fprintf(stderr,"%d : %f %f %f\n",fn,cp->xform[fn].c[0][0],cp->xform[fn].c[1][0],cp->xform[fn].c[2][0]);

   const double2 t = apply_affine ((double2) { p[0], p[1] }, cp->xform[fn].c);

   /* Pre-xforms go here, and modify the f.tx and f.ty values */
   if (cp->xform[fn].has_preblur!=0.0)
      var67_pre_blur(t, &f, cp->xform[fn].has_preblur);

   /* Always calculate sumsq and sqrt */
   f.precalc_sumsq = sum(t*t);
   f.precalc_sqrt = sqrt(f.precalc_sumsq);

   /* Check to see if we can precalculate any parts */
   /* Precalculate atanxy, sin, cos */
   if (cp->xform[fn].precalc_atan_xy_flag > 0) {
      f.precalc_atan = atan2(t[0],t[1]);
   }
   
   if (cp->xform[fn].precalc_angles_flag > 0) {
      f.precalc_sina = t[0] / f.precalc_sqrt;
      f.precalc_cosa = t[1] / f.precalc_sqrt;
   }

   /* Precalc atanyx */
   if (cp->xform[fn].precalc_atan_yx_flag > 0) {
      f.precalc_atanyx = atan2(t[1],t[0]);
   }

   f.xform = &(cp->xform[fn]);

   double2 accum = (double2) {0.0, 0.0};
   for (var_n=0; var_n < cp->xform[fn].num_active_vars; var_n++) {
      
      weight = cp->xform[fn].active_var_weights[var_n];

      switch (cp->xform[fn].varFunc[var_n]) {
 
         case (VAR_LINEAR):
                accum += var0_linear(t, &f, weight); break;               
         case (VAR_SINUSOIDAL):
                accum += var1_sinusoidal(t, &f, weight); break;
         case (VAR_SPHERICAL):
                accum += var2_spherical(t, &f, weight); break;
         case (VAR_SWIRL):
                accum += var3_swirl(t, &f, weight); break;
         case (VAR_HORSESHOE):
                accum += var4_horseshoe(t, &f, weight); break;               
         case (VAR_POLAR): 
                accum += var5_polar(t, &f, weight); break;
         case (VAR_HANDKERCHIEF):
                accum += var6_handkerchief(t, &f, weight); break;               
         case (VAR_HEART):
                accum += var7_heart(t, &f, weight); break;               
         case (VAR_DISC):
                accum += var8_disc(t, &f, weight); break;               
         case (VAR_SPIRAL):
                accum += var9_spiral(t, &f, weight); break;               
         case (VAR_HYPERBOLIC):
                accum += var10_hyperbolic(t, &f, weight); break;               
         case (VAR_DIAMOND):
                accum += var11_diamond(t, &f, weight); break;               
         case (VAR_EX):
                accum += var12_ex(t, &f, weight); break;               
         case (VAR_JULIA): 
                accum += var13_julia(t, &f, weight); break;               
         case (VAR_BENT):
                accum += var14_bent(t, &f, weight); break;
         case (VAR_WAVES):
                accum += var15_waves(t, &f, weight); break;
         case (VAR_FISHEYE): 
                accum += var16_fisheye(t, &f, weight); break;
         case (VAR_POPCORN):
                accum += var17_popcorn(t, &f, weight); break;
         case (VAR_EXPONENTIAL):
                accum += var18_exponential(t, &f, weight); break;
         case (VAR_POWER): 
                accum += var19_power(t, &f, weight); break;               
         case (VAR_COSINE):
                accum += var20_cosine(t, &f, weight); break;
         case (VAR_RINGS):
                accum += var21_rings(t, &f, weight); break;
         case (VAR_FAN):
                accum += var22_fan(t, &f, weight); break;
         case (VAR_BLOB):
                accum += var23_blob(t, &f, weight); break;               
         case (VAR_PDJ):
                accum += var24_pdj(t, &f, weight); break;
         case (VAR_FAN2):
                accum += var25_fan2(t, &f, weight); break;
         case (VAR_RINGS2): 
                accum += var26_rings2(t, &f, weight); break;              
         case (VAR_EYEFISH): 
                accum += var27_eyefish(t, &f, weight); break;             
         case (VAR_BUBBLE):
                accum += var28_bubble(t, &f, weight); break;
         case (VAR_CYLINDER):
                accum += var29_cylinder(t, &f, weight); break;
         case (VAR_PERSPECTIVE):
                accum += var30_perspective(t, &f, weight); break;
         case (VAR_NOISE):
                accum += var31_noise(t, &f, weight); break;
         case (VAR_JULIAN): 
                accum += var32_juliaN_generic(t, &f, weight); break;            
         case (VAR_JULIASCOPE):
                accum += var33_juliaScope_generic(t, &f, weight);break;
         case (VAR_BLUR):
                accum += var34_blur(t, &f, weight); break;
         case (VAR_GAUSSIAN_BLUR):
                accum += var35_gaussian(t, &f, weight); break;
         case (VAR_RADIAL_BLUR):
                accum += var36_radial_blur(t, &f, weight); break;
         case (VAR_PIE):
                accum += var37_pie(t, &f, weight); break;
         case (VAR_NGON):
                accum += var38_ngon(t, &f, weight); break;          
         case (VAR_CURL):
                accum += var39_curl(t, &f, weight); break;
         case (VAR_RECTANGLES):
                accum += var40_rectangles(t, &f, weight); break;
         case (VAR_ARCH):
                accum += var41_arch(t, &f, weight); break;
         case (VAR_TANGENT):
                accum += var42_tangent(t, &f, weight); break;
         case (VAR_SQUARE):
                accum += var43_square(t, &f, weight); break;
         case (VAR_RAYS):
                accum += var44_rays(t, &f, weight); break;
         case (VAR_BLADE): 
                accum += var45_blade(t, &f, weight); break;              
         case (VAR_SECANT2): 
                accum += var46_secant2(t, &f, weight); break;               
         case (VAR_TWINTRIAN): 
                accum += var47_twintrian(t, &f, weight); break;               
         case (VAR_CROSS):
                accum += var48_cross(t, &f, weight); break;
         case (VAR_DISC2):
                accum += var49_disc2(t, &f, weight); break;            
         case (VAR_SUPER_SHAPE):
                accum += var50_supershape(t, &f, weight); break;
         case (VAR_FLOWER):
                accum += var51_flower(t, &f, weight); break;            
         case (VAR_CONIC):
                accum += var52_conic(t, &f, weight); break;
         case (VAR_PARABOLA): 
                accum += var53_parabola(t, &f, weight); break;              
         case (VAR_BENT2): 
                accum += var54_bent2(t, &f, weight); break;              
         case (VAR_BIPOLAR): 
                accum += var55_bipolar(t, &f, weight); break;              
         case (VAR_BOARDERS): 
                accum += var56_boarders(t, &f, weight); break;              
         case (VAR_BUTTERFLY): 
                accum += var57_butterfly(t, &f, weight); break;              
         case (VAR_CELL): 
                accum += var58_cell(t, &f, weight); break;              
         case (VAR_CPOW): 
                accum += var59_cpow(t, &f, weight); break;              
         case (VAR_CURVE): 
                accum += var60_curve(t, &f, weight); break;              
         case (VAR_EDISC): 
                accum += var61_edisc(t, &f, weight); break;              
         case (VAR_ELLIPTIC): 
                accum += var62_elliptic(t, &f, weight); break;              
         case (VAR_ESCHER): 
                accum += var63_escher(t, &f, weight); break;              
         case (VAR_FOCI): 
                accum += var64_foci(t, &f, weight); break;              
         case (VAR_LAZYSUSAN): 
                accum += var65_lazysusan(t, &f, weight); break;              
         case (VAR_LOONIE): 
                accum += var66_loonie(t, &f, weight); break;              
         case (VAR_MODULUS): 
                accum += var68_modulus(t, &f, weight); break;              
         case (VAR_OSCILLOSCOPE): 
                accum += var69_oscope(t, &f, weight); break;              
         case (VAR_POLAR2): 
                accum += var70_polar2(t, &f, weight); break;              
         case (VAR_POPCORN2): 
                accum += var71_popcorn2(t, &f, weight); break;              
         case (VAR_SCRY): 
                accum += var72_scry(t, &f, weight); break;              
         case (VAR_SEPARATION): 
                accum += var73_separation(t, &f, weight); break;              
         case (VAR_SPLIT):
                accum += var74_split(t, &f, weight); break;
         case (VAR_SPLITS):
                accum += var75_splits(t, &f, weight); break;
         case (VAR_STRIPES):
                accum += var76_stripes(t, &f, weight); break;
         case (VAR_WEDGE):
                accum += var77_wedge(t, &f, weight); break;
         case (VAR_WEDGE_JULIA):
                accum += var78_wedge_julia(t, &f, weight); break;
         case (VAR_WEDGE_SPH):
                accum += var79_wedge_sph(t, &f, weight); break;
         case (VAR_WHORL):
                accum += var80_whorl(t, &f, weight); break;
         case (VAR_WAVES2):
                accum += var81_waves2(t, &f, weight); break;
         case (VAR_EXP):
                accum += var82_exp(t, &f, weight); break;
         case (VAR_LOG):
                accum += var83_log(t, &f, weight); break;
         case (VAR_SIN):
                accum += var84_sin(t, &f, weight); break;
         case (VAR_COS):
                accum += var85_cos(t, &f, weight); break;
         case (VAR_TAN):
                accum += var86_tan(t, &f, weight); break;
         case (VAR_SEC):
                accum += var87_sec(t, &f, weight); break;
         case (VAR_CSC):
                accum += var88_csc(t, &f, weight); break;
         case (VAR_COT):
                accum += var89_cot(t, &f, weight); break;
         case (VAR_SINH):
                accum += var90_sinh(t, &f, weight); break;
         case (VAR_COSH):
                accum += var91_cosh(t, &f, weight); break;
         case (VAR_TANH):
                accum += var92_tanh(t, &f, weight); break;
         case (VAR_SECH):
                accum += var93_sech(t, &f, weight); break;
         case (VAR_CSCH):
                accum += var94_csch(t, &f, weight); break;
         case (VAR_COTH):
                accum += var95_coth(t, &f, weight); break;
         case (VAR_AUGER):
                accum += var96_auger(t, &f, weight); break;
         case (VAR_FLUX):
                accum += var97_flux(t, &f, weight); break;
         case (VAR_MOBIUS):
                accum += var98_mobius(t, &f, weight); break;
      }

   }
   double2 q01;
   /* apply the post transform */
   if (cp->xform[fn].has_post) {
      q01 = apply_affine (accum, cp->xform[fn].post);
   } else {
      q01 = accum;
   }

   /* Check for badvalues and return randoms if bad */
   if (badvalue(q01[0]) || badvalue(q01[1])) {
      *q_ret = (double4) { rand_d11(rc), rand_d11(rc), q23[0], q23[1] };
      return(1);
   } else {
	  *q_ret = (double4) { q01[0], q01[1], q23[0], q23[1] };
      return(0);
   }

}

void initialize_xforms(flam3_genome *thiscp, int start_here) {

   int i,j;
   for (i = start_here ; i < thiscp->num_xforms ; i++) {
      thiscp->xform[i].padding = 0;
      thiscp->xform[i].density = 0.0;
      thiscp->xform[i].color_speed = 0.5;
      thiscp->xform[i].animate = 1.0;
      thiscp->xform[i].color = i&1;
      thiscp->xform[i].opacity = 1.0;
      thiscp->xform[i].var[0] = 1.0;
      thiscp->xform[i].motion_freq = 0;
      thiscp->xform[i].motion_func = 0;
      thiscp->xform[i].num_motion = 0;
      thiscp->xform[i].motion = NULL;
      for (j = 1; j < flam3_nvariations; j++)
         thiscp->xform[i].var[j] = 0.0;
      thiscp->xform[i].c[0][0] = 1.0;
      thiscp->xform[i].c[0][1] = 0.0;
      thiscp->xform[i].c[1][0] = 0.0;
      thiscp->xform[i].c[1][1] = 1.0;
      thiscp->xform[i].c[2][0] = 0.0;
      thiscp->xform[i].c[2][1] = 0.0;
      thiscp->xform[i].post[0][0] = 1.0;
      thiscp->xform[i].post[0][1] = 0.0;
      thiscp->xform[i].post[1][0] = 0.0;
      thiscp->xform[i].post[1][1] = 1.0;
      thiscp->xform[i].post[2][0] = 0.0;
      thiscp->xform[i].post[2][1] = 0.0;
      thiscp->xform[i].wind[0] = 0.0;
      thiscp->xform[i].wind[1] = 0.0;
      thiscp->xform[i].blob_low = 0.0;
      thiscp->xform[i].blob_high = 1.0;
      thiscp->xform[i].blob_waves = 1.0;
      thiscp->xform[i].pdj_a = 0.0;
      thiscp->xform[i].pdj_b = 0.0;
      thiscp->xform[i].pdj_c = 0.0;
      thiscp->xform[i].pdj_d = 0.0;
      thiscp->xform[i].fan2_x = 0.0;
      thiscp->xform[i].fan2_y = 0.0;
      thiscp->xform[i].rings2_val = 0.0;
      thiscp->xform[i].perspective_angle = 0.0;
      thiscp->xform[i].perspective_dist = 0.0;
      thiscp->xform[i].persp_vsin = 0.0;
      thiscp->xform[i].persp_vfcos = 0.0;
      thiscp->xform[i].radial_blur_angle = 0.0;
      thiscp->xform[i].disc2_rot = 0.0;
      thiscp->xform[i].disc2_twist = 0.0;
      thiscp->xform[i].disc2_sinadd = 0.0;
      thiscp->xform[i].disc2_cosadd = 0.0;
      thiscp->xform[i].disc2_timespi = 0.0;
      thiscp->xform[i].flower_petals = 0.0;
      thiscp->xform[i].flower_holes = 0.0;
      thiscp->xform[i].parabola_height = 0.0;
      thiscp->xform[i].parabola_width = 0.0;
      thiscp->xform[i].bent2_x = 1.0;
      thiscp->xform[i].bent2_y = 1.0;
      thiscp->xform[i].bipolar_shift = 0.0;
      thiscp->xform[i].cell_size = 1.0;
      thiscp->xform[i].cpow_r = 1.0;
      thiscp->xform[i].cpow_i = 0.0;
      thiscp->xform[i].cpow_power = 1.0;
      thiscp->xform[i].curve_xamp = 0.0;
      thiscp->xform[i].curve_yamp = 0.0;
      thiscp->xform[i].curve_xlength = 1.0;
      thiscp->xform[i].curve_ylength = 1.0;
      thiscp->xform[i].escher_beta = 0.0;
      thiscp->xform[i].lazysusan_space = 0.0;
      thiscp->xform[i].lazysusan_twist = 0.0;
      thiscp->xform[i].lazysusan_spin = 0.0;
      thiscp->xform[i].lazysusan_x = 0.0;
      thiscp->xform[i].lazysusan_y = 0.0;
      thiscp->xform[i].modulus_x = 0.0;
      thiscp->xform[i].modulus_y = 0.0;
      thiscp->xform[i].oscope_separation = 1.0;
      thiscp->xform[i].oscope_frequency = M_PI;
      thiscp->xform[i].oscope_amplitude = 1.0;
      thiscp->xform[i].oscope_damping = 0.0;
      thiscp->xform[i].popcorn2_c = 0.0;
      thiscp->xform[i].popcorn2_x = 0.0;
      thiscp->xform[i].popcorn2_y = 0.0;
      thiscp->xform[i].separation_x = 0.0;
      thiscp->xform[i].separation_xinside = 0.0;
      thiscp->xform[i].separation_y = 0.0;
      thiscp->xform[i].separation_yinside = 0.0;
      thiscp->xform[i].split_xsize = 0.0;
      thiscp->xform[i].split_ysize = 0.0;
      thiscp->xform[i].splits_x = 0.0;
      thiscp->xform[i].splits_y = 0.0;
      thiscp->xform[i].stripes_space = 0.0;
      thiscp->xform[i].stripes_warp = 0.0;
      thiscp->xform[i].wedge_angle = 0.0;
      thiscp->xform[i].wedge_hole = 0.0;
      thiscp->xform[i].wedge_count = 1.0;
      thiscp->xform[i].wedge_swirl = 0.0;
      thiscp->xform[i].wedge_sph_angle = 0.0;
      thiscp->xform[i].wedge_sph_hole = 0.0;
      thiscp->xform[i].wedge_sph_count = 1.0;
      thiscp->xform[i].wedge_sph_swirl = 0.0;

      thiscp->xform[i].wedge_julia_power = 1.0;
      thiscp->xform[i].wedge_julia_dist = 0.0;
      thiscp->xform[i].wedge_julia_count = 1.0;
      thiscp->xform[i].wedge_julia_angle = 0.0;
      thiscp->xform[i].wedgeJulia_cf = 0.0;
      thiscp->xform[i].wedgeJulia_cn = 0.5;
      thiscp->xform[i].wedgeJulia_rN = 1.0;
      thiscp->xform[i].whorl_inside = 0.0;
      thiscp->xform[i].whorl_outside = 0.0;
      
      thiscp->xform[i].waves2_scalex = 0.0;       
      thiscp->xform[i].waves2_scaley = 0.0;       
      thiscp->xform[i].waves2_freqx = 0.0;       
      thiscp->xform[i].waves2_freqy = 0.0;  
      
      thiscp->xform[i].auger_freq = 1.0;
      thiscp->xform[i].auger_weight = 0.5;
      thiscp->xform[i].auger_sym = 0.0;
      thiscp->xform[i].auger_scale = 1.0;     

      thiscp->xform[i].flux_spread = 0.0;
       
      thiscp->xform[i].julian_power = 1.0;
      thiscp->xform[i].julian_dist = 1.0;
      thiscp->xform[i].julian_rN = 1.0;
      thiscp->xform[i].julian_cn = 0.5;
      thiscp->xform[i].juliascope_power = 1.0;
      thiscp->xform[i].juliascope_dist = 1.0;
      thiscp->xform[i].juliascope_rN = 1.0;
      thiscp->xform[i].juliascope_cn = 0.5;
      thiscp->xform[i].radialBlur_spinvar = 0.0;
      thiscp->xform[i].radialBlur_zoomvar = 1.0;
      thiscp->xform[i].pie_slices = 6.0;
      thiscp->xform[i].pie_rotation = 0.0;
      thiscp->xform[i].pie_thickness = 0.5;
      thiscp->xform[i].ngon_sides = 5;
      thiscp->xform[i].ngon_power = 3;
      thiscp->xform[i].ngon_circle = 1;
      thiscp->xform[i].ngon_corners = 2;
      thiscp->xform[i].curl_c1 = 1.0;
      thiscp->xform[i].curl_c2 = 0.0;
      thiscp->xform[i].rectangles_x = 1.0;
      thiscp->xform[i].rectangles_y = 1.0;
      thiscp->xform[i].amw_amp = 1.0;
      thiscp->xform[i].super_shape_rnd = 0.0;
      thiscp->xform[i].super_shape_m = 0.0;
      thiscp->xform[i].super_shape_n1 = 1.0;
      thiscp->xform[i].super_shape_n2 = 1.0;
      thiscp->xform[i].super_shape_n3 = 1.0;
      thiscp->xform[i].super_shape_holes = 0.0;
      thiscp->xform[i].conic_eccentricity = 1.0;
      thiscp->xform[i].conic_holes = 0.0;

      thiscp->xform[i].mobius_re_a = 0.0;
      thiscp->xform[i].mobius_re_b = 0.0;
      thiscp->xform[i].mobius_re_c = 0.0;
      thiscp->xform[i].mobius_re_d = 0.0;
      thiscp->xform[i].mobius_im_a = 0.0;
      thiscp->xform[i].mobius_im_b = 0.0;
      thiscp->xform[i].mobius_im_c = 0.0;
      thiscp->xform[i].mobius_im_d = 0.0;
   }
}
