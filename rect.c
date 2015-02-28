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

#include <assert.h>
#include <stdlib.h>

#include "private.h"
#include "variations.h"
#include "palettes.h"
#include "math.h"
#include "rect.h"

/* allow this many iterations for settling into attractor */
#define FUSE_27 15
#define FUSE_28 100

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
   /* mutex for bucket accumulator */
   pthread_mutex_t bucket_mutex;
} flam3_iter_constants;

typedef struct {
   flam3_genome cp; /* Full copy of genome for use by the thread */
   flam3_iter_constants *fic; /* Constants for render */
   /* thread number */
   size_t i;
} flam3_thread_helper;

/*	Lookup color [0,1]
 */
static double4 color_palette_lookup (const double color,
		const color_palette_mode mode, const flam3_palette map,
		const unsigned int map_count) {
	assert (color >= 0.0 && color <= 1.0);

	switch (mode) {
		case PALETTE_MODE_LINEAR: {
			const double ix = color * map_count;
			const double bottomix = floor (ix);
			const double frac = ix - bottomix;
			const unsigned int intix = bottomix;

			if (intix == map_count-1) {
				return map[intix].color;
			} else {
				return map[intix].color * (1.0-frac) +
					map[intix+1].color * frac;
			}
			break;
		}

		case PALETTE_MODE_STEP: {
			const unsigned int intix = nearbyint (color * map_count);
			return map[intix].color;
			break;
		}

		default:
			assert (0);
			break;
	}
}

static void *iter_thread(void *fth) {
   double sub_batch;
   int j;
   flam3_thread_helper *fthp = (flam3_thread_helper *)fth;
   flam3_iter_constants *ficp = fthp->fic;
   int SBS = ficp->spec->sub_batch_size;
   int fuse;
   int cmap_size = ficp->cmap_size;
   double4 *iter_storage;
   randctx rc;

   rand_seed (&rc);

   int ret = posix_memalign ((void **) &iter_storage, sizeof (*iter_storage),
		   SBS * sizeof (*iter_storage));
   assert (ret == 0);
   assert (iter_storage != NULL);  

   fuse = (ficp->spec->earlyclip) ? FUSE_28 : FUSE_27;

   for (sub_batch = 0; sub_batch < ficp->batch_size; sub_batch+=SBS) {
      int sub_batch_size, badcount;
      /* sub_batch is double so this is sketchy */
      sub_batch_size = (sub_batch + SBS > ficp->batch_size) ?
                           (ficp->batch_size - sub_batch) : SBS;

      /* Seed iterations */
      const double4 start = (double4) {
	                        rand_d11(&rc),
                            rand_d11(&rc),
                            rand_d01(&rc),
                            rand_d01(&rc),
							};

      /* Execute iterations */
      badcount = flam3_iterate(&(fthp->cp), sub_batch_size, fuse, start, iter_storage, ficp->xform_distrib, &rc);

      /* Lock mutex for access to accumulator */
      pthread_mutex_lock(&ficp->bucket_mutex);

      /* Add the badcount to the counter */
      ficp->badvals += badcount;

      /* Put them in the bucket accumulator */
      for (j = 0; j < sub_batch_size; j++) {
         double4 p = iter_storage[j];

         if (fthp->cp.rotate != 0.0) {
		 	const double2 p01 = (double2) { p[0], p[1] };
		 	const double2 rotatedp = apply_affine (p01, ficp->rot);
		 	p[0] = rotatedp[0];
		 	p[1] = rotatedp[1];
         }

		 /* Skip if out of bounding box or invisible */
         if (p[0] >= ficp->bounds[0] && p[1] >= ficp->bounds[1] &&
		     p[0] <= ficp->bounds[2] && p[1] <= ficp->bounds[3] &&
			 p[3] > 0) {
			const size_t ix = (int)(ficp->ws0 * p[0] - ficp->wb0s0) + ficp->width * (int)(ficp->hs1 * p[1] - ficp->hb1s1);
#if HAVE_BUILTIN_PREFETCH
			/* prefetch for reading (0) with no locality (0). This (partially)
			 * hides the load latency for the += operation at the end of this
			 * block */
			__builtin_prefetch (&ficp->buckets[ix], 0, 0);
#endif

			double4 interpcolor = color_palette_lookup (p[2],
					fthp->cp.palette_mode, ficp->dmap, cmap_size);

            const double logvis = p[3];
            if (logvis != 1.0) {
			   interpcolor *= logvis;
            }

            ficp->buckets[ix] += interpcolor;

         }
      }
      
      /* Release mutex */
      pthread_mutex_unlock(&ficp->bucket_mutex);

   }

   free (iter_storage);
   return NULL;
}

/*	Perform clipping
 */
static double4 clip (const double4 in, const double g, const double linrange,
		const double highpow, const double vibrancy) {
	double alpha, ls;

	if (in[3] <= 0.0) {
		alpha = 0.0;
		ls = 0.0;
	} else {
		alpha = flam3_calc_alpha (in[3], g, linrange);
		ls = vibrancy * alpha / in[3];
		alpha = clamp (alpha, 0.0, 1.0);
	}

	double4 newrgb = flam3_calc_newrgb (in, ls, highpow);
	newrgb += (1.0-vibrancy) * pow_d4 (in, g);
	if (alpha > 0.0) {
		newrgb /= alpha;
	} else {
		newrgb = (double4) {0, 0, 0, 0};
	}
	newrgb[3] = alpha;
	newrgb = clamp_d4 (newrgb, 0.0, 1.0);

	return newrgb;
}

int render_parallel (flam3_frame *spec, void *out, stat_struct *stats) {
   long nbuckets;
   int i, j, k;
   double ppux=0, ppuy=0;
   int image_width, image_height;    /* size of the image to produce */
   int out_width;
   int bytes_per_channel = spec->bytes_per_channel;
   double highpow;
   flam3_palette dmap;
   double vibrancy = 0.0;
   double gamma = 0.0;
   int vib_gam_n = 0;
   flam3_genome cp;
   unsigned short *xform_distrib;
   flam3_iter_constants fic;
   flam3_thread_helper *fth;
   pthread_attr_t pt_attr;
   pthread_t *myThreads=NULL;
   int thi;
   int cmap_size;
   
   fic.badvals = 0;
   fic.aborted = 0;

   stats->num_iters = 0;

   /* correct for apophysis's use of 255 colors in the palette rather than all 256 */
   cmap_size = 256;

   memset(&cp,0, sizeof(flam3_genome));

   /* interpolate and get a control point                      */
   flam3_interpolate(spec->genomes, spec->ngenomes, spec->time, 0, &cp);
   highpow = cp.highlight_power;

   /* Initialize the thread helper structures */
   fth = (flam3_thread_helper *)calloc(spec->nthreads,sizeof(flam3_thread_helper));
   for (i=0;i<spec->nthreads;i++)
      fth[i].cp.final_xform_index=-1;
      
   /* Set up the output image dimensions, adjusted for scanline */   
   const unsigned int channels = 4;
   image_width = cp.width;
   out_width = image_width;
   image_height = cp.height;

   /* Allocate the space required to render the image */
   fic.height = image_height;
   fic.width  = image_width;

   nbuckets = (long)fic.width * (long)fic.height;

   double4 *buckets;
   int ret = posix_memalign ((void **) &buckets, sizeof (*buckets),
                             nbuckets * sizeof (*buckets));
   assert (ret == 0);
   assert (buckets != NULL);
   double4 *accumulate;
   ret = posix_memalign ((void **) &accumulate, sizeof (*accumulate),
                         nbuckets * sizeof (*accumulate));
   assert (ret == 0);
   assert (accumulate != NULL);

   memset(accumulate, 0, sizeof(*accumulate) * nbuckets);


   /* Batch loop - outermost */
   {
      double sample_density=0.0;
      double k1, area, k2;

      memset(buckets, 0, sizeof(*buckets) * nbuckets);

      {

         /* Get the xforms ready to render */
         if (prepare_precalc_flags(&cp)) {
            fprintf(stderr,"prepare xform pointers returned error: aborting.\n");
            return(1);
         }
         xform_distrib = flam3_create_xform_distrib(&cp);
         if (xform_distrib==NULL) {
            fprintf(stderr,"create xform distrib returned error: aborting.\n");
            return(1);
         }

         /* compute the colormap entries.                             */
         /* the input colormap is 256 long with entries from 0 to 1.0 */
         for (j = 0; j < CMAP_SIZE; j++) {
            dmap[j].index = cp.palette[(j * 256) / CMAP_SIZE].index / 256.0;
            for (k = 0; k < 4; k++)
               dmap[j].color[k] = cp.palette[(j * 256) / CMAP_SIZE].color[k];
         }

         /* compute camera */
         {
            double corner0, corner1;
            double scale;

            if (cp.sample_density <= 0.0) {
              fprintf(stderr,
                 "sample density (quality) must be greater than zero,"
                 " not %g.\n", cp.sample_density);
              return(1);
            }

            scale = pow(2.0, cp.zoom);
            sample_density = cp.sample_density * scale * scale;

            ppux = cp.pixels_per_unit * scale;
            ppuy = ppux;
            ppux /=  spec->pixel_aspect_ratio;
            corner0 = cp.center[0] - image_width / ppux / 2.0;
            corner1 = cp.center[1] - image_height / ppuy / 2.0;
            fic.bounds[0] = corner0;
            fic.bounds[1] = corner1;
            fic.bounds[2] = corner0 + image_width  / ppux;
            fic.bounds[3] = corner1 + image_height / ppuy;
            fic.size[0] = 1.0 / (fic.bounds[2] - fic.bounds[0]);
            fic.size[1] = 1.0 / (fic.bounds[3] - fic.bounds[1]);
			rotate_center ((double2) { cp.rot_center[0], cp.rot_center[1] },
						   cp.rotate, fic.rot);
            fic.ws0 = fic.width * fic.size[0];
            fic.wb0s0 = fic.ws0 * fic.bounds[0];
            fic.hs1 = fic.height * fic.size[1];
            fic.hb1s1 = fic.hs1 * fic.bounds[1];

         }

         /* number of samples is based only on the output image size */
         double nsamples = sample_density * image_width * image_height;
         
         /* how many of these samples are rendered in this loop? */
         double batch_size = nsamples;

         stats->num_iters += batch_size;
                  
         /* Fill in the iter constants */
         fic.xform_distrib = xform_distrib;
         fic.spec = spec;
         fic.batch_size = batch_size / (double)spec->nthreads;
         fic.cmap_size = cmap_size;

         fic.dmap = (flam3_palette_entry *)dmap;
         fic.buckets = (void *)buckets;
         
         /* Initialize the thread helper structures */
         for (thi = 0; thi < spec->nthreads; thi++) {
            fth[thi].fic = &fic;
            fth[thi].i = thi;
            flam3_copy(&(fth[thi].cp),&cp);
         }

         /* Let's make some threads */
         myThreads = (pthread_t *)malloc(spec->nthreads * sizeof(pthread_t));

         pthread_mutex_init(&fic.bucket_mutex, NULL);

         pthread_attr_init(&pt_attr);
         pthread_attr_setdetachstate(&pt_attr,PTHREAD_CREATE_JOINABLE);

         for (thi=0; thi <spec->nthreads; thi ++)
            pthread_create(&myThreads[thi], &pt_attr, (void *)iter_thread, (void *)(&(fth[thi])));

         pthread_attr_destroy(&pt_attr);

         /* Wait for them to return */
         for (thi=0; thi < spec->nthreads; thi++)
            pthread_join(myThreads[thi], NULL);

         pthread_mutex_destroy(&fic.bucket_mutex);
         
         free(myThreads);
         
         /* Free the xform_distrib array */
         free(xform_distrib);
             
         if (fic.aborted) {
            goto done;
         }

         vibrancy += cp.vibrancy;
         gamma += cp.gamma;
         vib_gam_n++;

      }

	  /* XXX: the original formula has a factor 268/256 in here, not sure why */
      k1 = cp.contrast * cp.brightness;
      area = image_width * image_height / (ppux * ppuy);
      k2 = 1.0 / (cp.contrast * area * sample_density);
#if 0
      printf("iw=%d,ih=%d,ppux=%f,ppuy=%f\n",image_width,image_height,ppux,ppuy);
      printf("contrast=%f, brightness=%f, PREFILTER=%d\n",
        cp.contrast, cp.brightness, PREFILTER_WHITE);
      printf("area = %f, WHITE_LEVEL=%d, sample_density=%f\n",
        area, WHITE_LEVEL, sample_density);
      printf("k1=%f,k2=%15.12f\n",k1,k2);
#endif

      for (j = 0; j < fic.height; j++) {
         for (i = 0; i < fic.width; i++) {
			const double4 c = buckets[i + j * fic.width];

            if (0.0 == c[3])
               continue;

            const double ls = (k1 * log(1.0 + c[3] * k2))/c[3];

            accumulate[i + j * fic.width] += c * ls;
         }
      }

   }

   /* filter the accumulation buffer down into the image */
   if (1) {
      int x, y;
      const double g = 1.0 / (gamma / vib_gam_n);

      double linrange = cp.gam_lin_thresh;

      vibrancy /= vib_gam_n;
      
      /* If we're in the early clip mode, perform this first step to  */
      /* apply the gamma correction and clipping before the spat filt */
      
      if (spec->earlyclip) {
         for (j = 0; j < fic.height; j++) {
            for (i = 0; i < fic.width; i++) {
               const double4 in = accumulate[i + j*fic.width];
			   accumulate[i + j*fic.width] = clip (in, g, linrange, highpow,
					   vibrancy);
            }
         }
      }

      /* Apply the spatial filter */
      y = 0;
      for (j = 0; j < image_height; j++) {
         x = 0;
         for (i = 0; i < image_width; i++) {
			double4 t = accumulate[x + y * fic.width];

            /* The old way, spatial filter first and then clip after gamma */
            if (!spec->earlyclip) {
			   t = clip (t, g, linrange, highpow, vibrancy);
            }

			const double maxval = (1 << (bytes_per_channel*8)) - 1;
			t = nearbyint_d4 (t * maxval);

			if (bytes_per_channel == 2) {
				uint16_t * const p = &((uint16_t *) out)[channels * (i + j * out_width)];
				p[0] = t[0];
				p[1] = t[1];
				p[2] = t[2];
				p[3] = t[3];
			} else if (bytes_per_channel == 1) {
				uint8_t * const p = &((uint8_t *) out)[channels * (i + j * out_width)];
				p[0] = t[0];
				p[1] = t[1];
				p[2] = t[2];
				p[3] = t[3];
			} else {
				assert (0);
			}

            x += 1;
         }
         y += 1;
      }
   }

 done:

   stats->badvals = fic.badvals;

   free(buckets);
   free(accumulate);
   /* We have to clear the cps in fth first */
   for (thi = 0; thi < spec->nthreads; thi++) {
      clear_cp(&(fth[thi].cp),0);
   }   
   free(fth);
   clear_cp(&cp,0);

   return(0);

}
