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

#ifndef variations_included
#define variations_included

#include "private.h"

void xform_precalc(flam3_genome *cp, int xi);
int prepare_precalc_flags(flam3_genome *);

int apply_xform(flam3_genome *cp, int fn, const double4 p, double4 *, randctx *rc);
void initialize_xforms(flam3_genome *thiscp, int start_here);
#endif
