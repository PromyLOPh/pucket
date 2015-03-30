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

#include "math.h"
#include "flam3.h"

int id_matrix(double2 s[3]);
double det_matrix(double s[2][2]);

#if 0
void interpolate_cmap(flam3_palette cmap, double blend,
              int index0, double hue0, int index1, double hue1, randctx * const);
#endif

void interpolate_catmull_rom(flam3_genome cps[], double t, flam3_genome *result);
void flam3_interpolate_n(flam3_genome *result, int ncp, flam3_genome *cpi, double *c, double stagger);
void flam3_align(flam3_genome *dst, flam3_genome *src, int nsrc);
