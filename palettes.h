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

#include <stdlib.h>
#include <stdbool.h>

#include "vector.h"

/* One palette */
typedef struct {
	/* number of items allocated */
	size_t allocated;
	/* number of items in palette */
	size_t count;
	double4 *color;
} palette;

/* A collection/array of palettes */
typedef struct {
	size_t count;
	palette *p;
} palette_collection;

double4 rgb2hsv(double4);
double4 hsv2rgb(double4);

#include "flam3.h"

void palette_add (palette * const p, const double4 c);
const palette *palette_random (const palette_collection * const pc, randctx * const rc);
void palette_copy (const palette * restrict const src, palette * restrict const dest);
void palette_rotate_hue (palette * const p, double rot);
bool palette_read_collection (const char * const filename,
		palette_collection * const pc);

