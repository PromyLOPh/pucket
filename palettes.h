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

#ifndef palettes_included
#define palettes_included


typedef struct {
    int number;
    char name[flam3_name_len];
    unsigned char colors[256][3];
} lib_palette;

#include "vector.h"

double4 rgb2hsv(double4);
double4 hsv2rgb(double4);

double flam3_calc_alpha(double density, double gamma, double linrange);
double4 flam3_calc_newrgb(double4 cbuf, double ls, double highpow);

#endif

