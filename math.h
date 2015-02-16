/*
    FLAM3 - cosmic recursive fractal flames
    Copyright (C) 1992-2009 Spotworks LLC
    Copyright (C) 2015 vlam3 contributors

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

#include "build/config.h"

#ifdef HAVE_AMDLIBM
#define REPLACE_WITH_AMDLIBM
#include <amdlibm.h>
#endif

#define clamp(a,min,max) (a > max ? max : (a < min ? min : a))

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
