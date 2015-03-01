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

#ifndef private_included
#define private_included

#include "flam3.h"
#include "build/config.h"
#include <stdlib.h>

#include <ctype.h>
#include <time.h>
#include <string.h>
#include <libxml/parser.h>

#include <unistd.h>
#include <libgen.h>

#include <math.h>

#define EPS (1e-10)
#define CMAP_SIZE 256
#define CMAP_SIZE_M1 255
#define rbit() (flam3_random_bit())
#define flam3_variation_none   (-1)
#define max_specified_vars     (100)
#define vlen(x) (sizeof(x)/sizeof(*x))

double flam3_sinc(double x);

#define  flam3_mitchell_b   (1.0 / 3.0)
#define  flam3_mitchell_c   (1.0 / 3.0)


#endif
