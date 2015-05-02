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

#include "flam3.h"

void xform_precalc (flam3_xform * const xform);
int prepare_precalc_flags(flam3_genome *);

int apply_xform(const flam3_xform * const xf, const unsigned int fn,
		const double4 p, double4 *q_ret, randctx * const rc);
void initialize_xforms(flam3_genome *thiscp, int start_here);

