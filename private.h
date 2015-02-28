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

#include <pthread.h>

#define EPS (1e-10)
#define CMAP_SIZE 256
#define CMAP_SIZE_M1 255
#define rbit() (flam3_random_bit())
#define flam3_variation_none   (-1)
#define max_specified_vars     (100)
#define vlen(x) (sizeof(x)/sizeof(*x))

/* Structures for passing parameters to iteration threads */
typedef struct {
   unsigned short *xform_distrib;    /* Distribution of xforms based on weights */
   flam3_frame *spec; /* Frame contains timing information */
   double bounds[4]; /* Corner coords of viewable area */
   double2 rot[3]; /* Rotation transformation */
   double size[2];
   int width, height; /* buffer width/height */
   double ws0, wb0s0, hs1, hb1s1; /* shortcuts for indexing */
   flam3_palette_entry *dmap; /* palette */
   double color_scalar; /* <1.0 if non-uniform motion blur is set */
   double4 *buckets; /* Points to the first accumulator */
   double badvals; /* accumulates all badvalue resets */
   double batch_size;
   int aborted, cmap_size;
   time_t *progress_timer;
   time_t *progress_timer_history;
   double *progress_history;
   int *progress_history_mark;
   /* mutex for bucket accumulator */
   pthread_mutex_t bucket_mutex;
   
} flam3_iter_constants;



typedef struct {
   flam3_genome cp; /* Full copy of genome for use by the thread */
   int first_thread;
   int timer_initialize;
   flam3_iter_constants *fic; /* Constants for render */
} flam3_thread_helper;

double flam3_sinc(double x);

#define  flam3_mitchell_b   (1.0 / 3.0)
#define  flam3_mitchell_c   (1.0 / 3.0)


#endif
