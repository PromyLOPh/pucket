/*
    Copyright (C) 1992-2009 Spotworks LLC
    Copyright (C) 2015 pucket contributors

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

#include <math.h>
#include <assert.h>

#include "build/config.h"

#ifdef HAVE_AMDLIBM
#define REPLACE_WITH_AMDLIBM
#include <amdlibm.h>
#undef nearbyint
#undef floor
#endif

#define clamp(a,min,max) (a > max ? max : (a < min ? min : a))

/*	Apply affine coordinate transformation
 */
inline double2 apply_affine (const double2 in, const double2 matrix[3]) {
	return matrix[0] * in[0] + matrix[1] * in[1] + matrix[2];
}

/*	Create affine rotation matrix, angle in degree
 */
inline void rotate (const double angle, double2 matrix[3]) {
	double s, c;
	sincos (angle * 2.0 * M_PI / 360.0, &s, &c);
	matrix[0] = (double2) { c, s };
	matrix[1] = (double2) { -s, c };
	matrix[2] = (double2) { 0.0, 0.0 };
}

/*	Create affine translation matrix
 */
inline void translate (const double2 xy, double2 matrix[3]) {
	matrix[0] = (double2) { 1.0, 0.0 };
	matrix[1] = (double2) { 0.0, 1.0 };
	matrix[2] = xy;
}

/*	Create affine scaling matrix
 */
inline void scale (const double2 xy, double2 matrix[3]) {
	matrix[0] = (double2) { xy[0], 0.0 };
	matrix[1] = (double2) { 0.0, xy[1] };
	matrix[2] = (double2) { 0.0, 0.0 };
}

/*	Multiply two affine matrices a, b and store the result in c.
 *
 *	The last row of each matrix is assumed to be 0, 0, 1.
 */
inline void matrixmul (const double2 a[3], const double2 b[3], double2 c[3]) {
	c[0] = a[0] * b[0][0] + a[1] * b[0][1];
	c[1] = a[0] * b[1][0] + a[1] * b[1][1];
	c[2] = a[0] * b[2][0] + a[1] * b[2][1] + a[2];
}

/*	Affine matrix that transforms rect from (x1, y1, x2, y2) into rect to
 */
inline void translate_rect (const double4 from, const double4 to,
		double2 matrix[3]) {
	const double2 from_edge = (double2) { from[0], from[1] },
			to_edge = (double2) { to[0], to[1] };
	/* first align one of A and B’s edges */
	double2 translate_edge[3];
	translate (to_edge - from_edge, translate_edge);
	/* then scale it up or down */
	double2 scale_rect[3];
	scale ((double2) { (to[2] - to[0])/(from[2] - from[0]),
			(to[3] - to[1])/(from[3] - from[1])}, scale_rect);
	/* the result is scale*translate (i.e. translate first) */
	matrixmul (scale_rect, translate_edge, matrix);
}

/* Create rotation around center. Note that matrix multiplication is
 * right-associative, thus A*B*C == A*(B*C) */
inline void rotate_center (const double2 center, const double angle, double2 out[3]) {
	double2 rot[3], trans_a[3], trans_b[3], tmp[3];
	translate (-1.0 * center, trans_a);
	rotate (angle, rot);
	translate (center, trans_b);
	matrixmul (rot, trans_a, tmp);
	matrixmul (trans_b, tmp, out);
}

inline double sum(const double2 in) {
	return in[0] + in[1];
}

inline void normalize (double * const a, const size_t n) {
	double sum = 0.0;
	for (unsigned int j = 0; j < n; j++) {
		sum += a[j];
	}
	assert (sum > 0.0);

	for (unsigned int j = 0; j < n; j++) {
		a[j] /= sum;
	}
}

#define max(a,b) ((a) > (b) ? (a) : (b))

/*	Vector wrapping function, could be replaced by true vector functions later
 */

inline double4 clamp_d4 (const double4 in, const double min, const double max) {
	return (double4) {
			clamp (in[0], min, max),
			clamp (in[1], min, max),
			clamp (in[2], min, max),
			clamp (in[3], min, max),
			};
}

inline double4 pow_d4 (const double4 in, double exp) {
	return (double4) {
			pow (in[0], exp),
			pow (in[1], exp),
			pow (in[2], exp),
			pow (in[3], exp),
			};
}

inline double4 nearbyint_d4 (const double4 in) {
	return (double4) {
			nearbyint (in[0]),
			nearbyint (in[1]),
			nearbyint (in[2]),
			nearbyint (in[3]),
			};
}
