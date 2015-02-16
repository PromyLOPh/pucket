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
#include "filters.h"
#include "variations.h"
#include "palettes.h"
#include "math.h"

/*
 * for batch
 *   generate de filters
 *   for temporal_sample_batch
 *     interpolate
 *     compute colormap
 *     for subbatch
 *       compute samples
 *       buckets += cmap[samples]
 *   accum += time_filter[temporal_sample_batch] * log[buckets] * de_filter
 * image = filter(accum)
 */

/* allow this many iterations for settling into attractor */
#define FUSE_27 15
#define FUSE_28 100

static void iter_thread(void *fth) {
   double sub_batch;
   int j;
   flam3_thread_helper *fthp = (flam3_thread_helper *)fth;
   flam3_iter_constants *ficp = fthp->fic;
   struct timespec pauset;
   int SBS = ficp->spec->sub_batch_size;
   int fuse;
   int cmap_size = ficp->cmap_size;
   int cmap_size_m1 = ficp->cmap_size-1;

   double eta = 0.0;
   
   fuse = (ficp->spec->earlyclip) ? FUSE_28 : FUSE_27;

   pauset.tv_sec = 0;
   pauset.tv_nsec = 100000000;


   if (fthp->timer_initialize) {
   	*(ficp->progress_timer) = 0;
   	memset(ficp->progress_timer_history,0,64*sizeof(time_t));
   	memset(ficp->progress_history,0,64*sizeof(double));
   	*(ficp->progress_history_mark) = 0;
   }
   
   for (sub_batch = 0; sub_batch < ficp->batch_size; sub_batch+=SBS) {
      int sub_batch_size, badcount;
      time_t newt = time(NULL);
      /* sub_batch is double so this is sketchy */
      sub_batch_size = (sub_batch + SBS > ficp->batch_size) ?
                           (ficp->batch_size - sub_batch) : SBS;
                           
      if (fthp->first_thread && newt != *(ficp->progress_timer)) {
         double percent = 100.0 *
             ((((sub_batch / (double) ficp->batch_size) + ficp->temporal_sample_num)
             / ficp->ntemporal_samples) + ficp->batch_num)/ficp->nbatches;
         int old_mark = 0;
         int ticker;

         if (ficp->spec->verbose)
            fprintf(stderr, "\rchaos: %5.1f%%", percent);
            
         *(ficp->progress_timer) = newt;
         if (ficp->progress_timer_history[*(ficp->progress_history_mark)] &&
                ficp->progress_history[*(ficp->progress_history_mark)] < percent)
            old_mark = *(ficp->progress_history_mark);

         if (percent > 0) {
            eta = (100 - percent) * (*(ficp->progress_timer) - ficp->progress_timer_history[old_mark])
                  / (percent - ficp->progress_history[old_mark]);

            if (ficp->spec->verbose) {
               ticker = (*(ficp->progress_timer)&1)?':':'.';
               if (eta < 1000)
                  ticker = ':';
               if (eta > 100)
                  fprintf(stderr, "  ETA%c %.1f minutes", ticker, eta / 60);
               else
                  fprintf(stderr, "  ETA%c %ld seconds ", ticker, (long) ceil(eta));
               fprintf(stderr, "              \r");
               fflush(stderr);
            }
         }

         ficp->progress_timer_history[*(ficp->progress_history_mark)] = *(ficp->progress_timer);
         ficp->progress_history[*(ficp->progress_history_mark)] = percent;
         *(ficp->progress_history_mark) = (*(ficp->progress_history_mark) + 1) % 64;
      }

      /* Custom progress function */
      if (ficp->spec->progress) {
         if (fthp->first_thread) {
         
            int rv;
         
            /* Recalculate % done, as the other calculation only updates once per second */
            double percent = 100.0 *
                ((((sub_batch / (double) ficp->batch_size) + ficp->temporal_sample_num)
                / ficp->ntemporal_samples) + ficp->batch_num)/ficp->nbatches;
                
            rv = (*ficp->spec->progress)(ficp->spec->progress_parameter, percent, 0, eta);
            
            if (rv==2) { /* PAUSE */
               
               time_t tnow = time(NULL);
               time_t tend;
               int lastpt;
               
               ficp->aborted = -1;
               
               do {
				   nanosleep(&pauset,NULL);
                  rv = (*ficp->spec->progress)(ficp->spec->progress_parameter, percent, 0, eta);
               } while (rv==2);
               
               /* modify the timer history to compensate for the pause */
               tend = time(NULL)-tnow;
               
               ficp->aborted = 0;

               for (lastpt=0;lastpt<64;lastpt++) {
                  if (ficp->progress_timer_history[lastpt]) {
                      ficp->progress_timer_history[lastpt] += tend;
                  }
               }
               
            }
                  
            if (rv==1) { /* ABORT */
				   ficp->aborted = 1;
               pthread_exit((void *)0);
            }
         } else {
            if (ficp->aborted<0) {

            do {
               nanosleep(&pauset,NULL);
            } while (ficp->aborted==-1);
            }
            if (ficp->aborted>0) pthread_exit((void *)0);
         }
      }

      /* Seed iterations */
      const double4 start = (double4) {
	                        rand_d11(&(fthp->rc)),
                            rand_d11(&(fthp->rc)),
                            rand_d01(&(fthp->rc)),
                            rand_d01(&(fthp->rc)),
							};

      /* Execute iterations */
      badcount = flam3_iterate(&(fthp->cp), sub_batch_size, fuse, start, fthp->iter_storage, ficp->xform_distrib, &(fthp->rc));

      /* Lock mutex for access to accumulator */
      pthread_mutex_lock(&ficp->bucket_mutex);

      /* Add the badcount to the counter */
      ficp->badvals += badcount;

      /* Put them in the bucket accumulator */
      for (j = 0; j < sub_batch_size; j++) {
         double p0, p1, p00, p11;
         double dbl_index0,dbl_frac;
         double4 interpcolor;
         int color_index0;
         const double4 p = fthp->iter_storage[j];

         if (fthp->cp.rotate != 0.0) {
            p00 = p[0] - fthp->cp.rot_center[0];
            p11 = p[1] - fthp->cp.rot_center[1];
            p0 = p00 * ficp->rot[0][0] + p11 * ficp->rot[0][1] + fthp->cp.rot_center[0];
            p1 = p00 * ficp->rot[1][0] + p11 * ficp->rot[1][1] + fthp->cp.rot_center[1];
         } else {
            p0 = p[0];
            p1 = p[1];
         }

         if (p0 >= ficp->bounds[0] && p1 >= ficp->bounds[1] && p0 <= ficp->bounds[2] && p1 <= ficp->bounds[3]) {

            double logvis=1.0;

            /* Skip if invisible */
            if (p[3]==0)
               continue;
            else
               logvis = p[3];
            
            dbl_index0 = p[2] * cmap_size;
            color_index0 = (int) (dbl_index0);
            
            if (flam3_palette_mode_linear == fthp->cp.palette_mode) {
               if (color_index0 < 0) {
                  color_index0 = 0;
                  dbl_frac = 0;
               } else if (color_index0 >= cmap_size_m1) {
                  color_index0 = cmap_size_m1-1;
                  dbl_frac = 1.0;
               } else {
                  /* interpolate between color_index0 and color_index0+1 */
                  dbl_frac = dbl_index0 - (double)color_index0;
               }
                        
               interpcolor = ficp->dmap[color_index0].color * (1.0-dbl_frac) + 
                                    ficp->dmap[color_index0+1].color * dbl_frac;
            } else { /* Palette mode step */
            
               if (color_index0 < 0) {
                  color_index0 = 0;
               } else if (color_index0 >= cmap_size_m1) {
                  color_index0 = cmap_size_m1;
               }
                        
               interpcolor = ficp->dmap[color_index0].color;
            }

            if (p[3]!=1.0) {
			   interpcolor *= logvis;
            }

            ficp->buckets[(int)(ficp->ws0 * p0 - ficp->wb0s0) + ficp->width * (int)(ficp->hs1 * p1 - ficp->hb1s1)] += interpcolor;

         }
      }
      
      /* Release mutex */
      pthread_mutex_unlock(&ficp->bucket_mutex);

   }
     pthread_exit((void *)0);
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

int render_rectangle(flam3_frame *spec, void *out,
			     int field, stat_struct *stats) {
   long nbuckets;
   int i, j, k, batch_num, temporal_sample_num;
   double nsamples, batch_size;
   double *filter, *temporal_filter, *temporal_deltas, *batch_filter;
   double ppux=0, ppuy=0;
   int image_width, image_height;    /* size of the image to produce */
   int out_width;
   int filter_width=0;
   int bytes_per_channel = spec->bytes_per_channel;
   int oversample;
   double highpow;
   int nbatches;
   int ntemporal_samples;
   flam3_palette dmap;
   int gutter_width;
   double vibrancy = 0.0;
   double gamma = 0.0;
   int vib_gam_n = 0;
   time_t progress_began=0;
   int verbose = spec->verbose;
   flam3_genome cp;
   unsigned short *xform_distrib;
   flam3_iter_constants fic;
   flam3_thread_helper *fth;
   pthread_attr_t pt_attr;
   pthread_t *myThreads=NULL;
   int thi;
   time_t tstart,tend;   
   double sumfilt;
   int cmap_size;
   
   /* Per-render progress timers */
   time_t progress_timer=0;
   time_t progress_timer_history[64];
   double progress_history[64];
   int progress_history_mark = 0;

   tstart = time(NULL);

   fic.badvals = 0;
   fic.aborted = 0;

   stats->num_iters = 0;

   /* correct for apophysis's use of 255 colors in the palette rather than all 256 */
   cmap_size = 256;

   memset(&cp,0, sizeof(flam3_genome));

   /* interpolate and get a control point                      */
   flam3_interpolate(spec->genomes, spec->ngenomes, spec->time, 0, &cp);
   oversample = cp.spatial_oversample;
   highpow = cp.highlight_power;
   nbatches = cp.nbatches;
   ntemporal_samples = cp.ntemporal_samples;

   if (nbatches < 1) {
       fprintf(stderr, "nbatches must be positive, not %d.\n", nbatches);
       return(1);
   }

   if (oversample < 1) {
       fprintf(stderr, "oversample must be positive, not %d.\n", oversample);
       return(1);
   }

   /* Initialize the thread helper structures */
   fth = (flam3_thread_helper *)calloc(spec->nthreads,sizeof(flam3_thread_helper));
   for (i=0;i<spec->nthreads;i++)
      fth[i].cp.final_xform_index=-1;
      
   /* Set up the output image dimensions, adjusted for scanline */   
   const unsigned int channels = 4;
   image_width = cp.width;
   out_width = image_width;
   if (field) {
      image_height = cp.height / 2;
      
      if (field == flam3_field_odd)
         out = (unsigned char *)out + channels * bytes_per_channel * out_width;
         
      out_width *= 2;
   } else
      image_height = cp.height;


   /* Spatial Filter kernel creation */
   filter_width = flam3_create_spatial_filter(spec, field, &filter);
   
   /* handle error */
   if (filter_width<0) {
      fprintf(stderr,"flam3_create_spatial_filter returned error: aborting\n");
      return(1);
   }
      
   /* note we must free 'filter' at the end */

   /* batch filter */
   /* may want to revisit this at some point */
   batch_filter = (double *) malloc(sizeof(double) * nbatches);
   for (i=0; i<nbatches; i++)
      batch_filter[i]=1.0/(double)nbatches;

   /* temporal filter - we must free temporal_filter and temporal_deltas at the end */
   sumfilt = flam3_create_temporal_filter(nbatches*ntemporal_samples, 
                                          cp.temporal_filter_type,
                                          cp.temporal_filter_exp,
                                          cp.temporal_filter_width,
                                          &temporal_filter, &temporal_deltas);
                                                                                    

   /*
      the number of additional rows of buckets we put at the edge so
      that the filter doesn't go off the edge
   */
   gutter_width = (filter_width - oversample) / 2;

   /* Allocate the space required to render the image */
   fic.height = oversample * image_height + 2 * gutter_width;
   fic.width  = oversample * image_width  + 2 * gutter_width;

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
   double4 ** const iter_storage = malloc (spec->nthreads * sizeof (*iter_storage));
   assert (iter_storage != NULL);
   for (size_t i = 0; i < spec->nthreads; i++) {
      ret = posix_memalign ((void **) &iter_storage[i],
	                        sizeof (*iter_storage[i]),
							spec->sub_batch_size * sizeof (*iter_storage[i]));
	  assert (ret == 0);
	  assert (iter_storage[i] != NULL);
   }

   if (verbose) {
      fprintf(stderr, "chaos: ");
      progress_began = time(NULL);
   }

   memset(accumulate, 0, sizeof(*accumulate) * nbuckets);


   /* Batch loop - outermost */
   for (batch_num = 0; batch_num < nbatches; batch_num++) {
      double sample_density=0.0;
      double k1, area, k2;

      memset(buckets, 0, sizeof(*buckets) * nbuckets);

      /* Temporal sample loop */
      for (temporal_sample_num = 0; temporal_sample_num < ntemporal_samples; temporal_sample_num++) {

         double temporal_sample_time;
         double color_scalar = temporal_filter[batch_num*ntemporal_samples + temporal_sample_num];

         temporal_sample_time = spec->time +
            temporal_deltas[batch_num*ntemporal_samples + temporal_sample_num];

         /* Interpolate and get a control point */
         flam3_interpolate(spec->genomes, spec->ngenomes, temporal_sample_time, 0, &cp);

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
               dmap[j].color[k] = cp.palette[(j * 256) / CMAP_SIZE].color[k] * color_scalar;
         }

         /* compute camera */
         if (1) {
            double t0, t1, shift=0.0, corner0, corner1;
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
            ppuy = field ? (ppux / 2.0) : ppux;
            ppux /=  spec->pixel_aspect_ratio;
            switch (field) {
               case flam3_field_both: shift =  0.0; break;
               case flam3_field_even: shift = -0.5; break;
               case flam3_field_odd:  shift =  0.5; break;
            }
            shift = shift / ppux;
            t0 = (double) gutter_width / (oversample * ppux);
            t1 = (double) gutter_width / (oversample * ppuy);
            corner0 = cp.center[0] - image_width / ppux / 2.0;
            corner1 = cp.center[1] - image_height / ppuy / 2.0;
            fic.bounds[0] = corner0 - t0;
            fic.bounds[1] = corner1 - t1 + shift;
            fic.bounds[2] = corner0 + image_width  / ppux + t0;
            fic.bounds[3] = corner1 + image_height / ppuy + t1 + shift;
            fic.size[0] = 1.0 / (fic.bounds[2] - fic.bounds[0]);
            fic.size[1] = 1.0 / (fic.bounds[3] - fic.bounds[1]);
            fic.rot[0][0] = cos(cp.rotate * 2 * M_PI / 360.0);
            fic.rot[0][1] = -sin(cp.rotate * 2 * M_PI / 360.0);
            fic.rot[1][0] = -fic.rot[0][1];
            fic.rot[1][1] = fic.rot[0][0];
            fic.ws0 = fic.width * fic.size[0];
            fic.wb0s0 = fic.ws0 * fic.bounds[0];
            fic.hs1 = fic.height * fic.size[1];
            fic.hb1s1 = fic.hs1 * fic.bounds[1];

         }

         /* number of samples is based only on the output image size */
         nsamples = sample_density * image_width * image_height;
         
         /* how many of these samples are rendered in this loop? */
         batch_size = nsamples / (nbatches * ntemporal_samples);

         stats->num_iters += batch_size;
                  
         /* Fill in the iter constants */
         fic.xform_distrib = xform_distrib;
         fic.spec = spec;
         fic.batch_size = batch_size / (double)spec->nthreads;
         fic.temporal_sample_num = temporal_sample_num;
         fic.ntemporal_samples = ntemporal_samples;
         fic.batch_num = batch_num;
         fic.nbatches = nbatches;
         fic.cmap_size = cmap_size;

         fic.dmap = (flam3_palette_entry *)dmap;
         fic.color_scalar = color_scalar;
         fic.buckets = (void *)buckets;
         
         /* Need a timer per job */
         fic.progress_timer = &progress_timer;
         fic.progress_timer_history = &(progress_timer_history[0]);
         fic.progress_history = &(progress_history[0]);
         fic.progress_history_mark = &progress_history_mark;

         /* Initialize the thread helper structures */
         for (thi = 0; thi < spec->nthreads; thi++) {

            /* Create a new state for this thread */
			rand_seed (&fth[thi].rc);

            if (0==thi) {

               fth[thi].first_thread=1;
               if (0==batch_num && 0==temporal_sample_num)
               	fth[thi].timer_initialize = 1;
               else
               	fth[thi].timer_initialize = 0;
               	
            } else {
               fth[thi].first_thread=0;
	         	fth[thi].timer_initialize = 0;
            }

            fth[thi].iter_storage = iter_storage[thi];
            fth[thi].fic = &fic;
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
            if (verbose) fprintf(stderr, "\naborted!\n");
            goto done;
         }

         vibrancy += cp.vibrancy;
         gamma += cp.gamma;
         vib_gam_n++;

      }

	  /* XXX: the original formula has a factor 268/256 in here, not sure why */
      k1 = cp.contrast * cp.brightness * batch_filter[batch_num];
      area = image_width * image_height / (ppux * ppuy);
      k2 = (oversample * oversample * nbatches) /
             (cp.contrast * area * sample_density * sumfilt);
#if 0
      printf("iw=%d,ih=%d,ppux=%f,ppuy=%f\n",image_width,image_height,ppux,ppuy);
      printf("contrast=%f, brightness=%f, PREFILTER=%d, temporal_filter=%f\n",
        cp.contrast, cp.brightness, PREFILTER_WHITE, temporal_filter[batch_num]);
      printf("oversample=%d, nbatches=%d, area = %f, WHITE_LEVEL=%d, sample_density=%f\n",
        oversample, nbatches, area, WHITE_LEVEL, sample_density);
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

   if (verbose) {
     fprintf(stderr, "\rchaos: 100.0%%  took: %ld seconds   \n", time(NULL) - progress_began);
     fprintf(stderr, "filtering...");
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
            int ii, jj;
			double4 t = (double4) { 0.0, 0.0, 0.0, 0.0 };

            for (ii = 0; ii < filter_width; ii++) {
               for (jj = 0; jj < filter_width; jj++) {
                  const double k = filter[ii + jj * filter_width];
                  const double4 ac = accumulate[x+ii + (y+jj)*fic.width];
                  
				  t += k * ac;
               }
            }

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

            x += oversample;
         }
         y += oversample;
      }
   }

 done:

   stats->badvals = fic.badvals;

   free(temporal_filter);
   free(temporal_deltas);
   free(batch_filter);
   free(filter);
   free(buckets);
   free(accumulate);
   /* We have to clear the cps in fth first */
   for (thi = 0; thi < spec->nthreads; thi++) {
      clear_cp(&(fth[thi].cp),0);
	  free (iter_storage[thi]);
   }   
   free (iter_storage);
   free(fth);
   clear_cp(&cp,0);

   tend = time(NULL);
   stats->render_seconds = (int)(tend-tstart);
   
   return(0);

}
