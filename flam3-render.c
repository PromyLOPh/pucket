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

#include <limits.h>


#include "private.h"
#include "img.h"

int print_progress(void *foo, double fraction, int stage, double eta) {
    fprintf(stderr, "stage=%s progress=%g eta=%g\n", stage?"filtering":"chaos", fraction, eta);
  return 0;
}

int main(int argc, char **argv) {
   flam3_frame f;
   char *ai;
   flam3_genome *cps;
   int ncps;
   int i;
   void *image=NULL;
   FILE *fp;
   char fname[256];
   size_t this_size;
   char *prefix = args("prefix", "");
   char *out = args("out", NULL);
   const char *format = "png";
   int verbose = argi("verbose", 1);
   int bits = argi("bits", 33);
   int bpc = argi("bpc",8);
   int transparency = argi("transparency", 0);
   char *inf = getenv("in");
   double qs = argf("qs", 1.0);
   double ss = argf("ss", 1.0);
   double pixel_aspect = argf("pixel_aspect", 1.0);
   int sub_batch_size = argi("sub_batch_size",10000);
   int name_enable = argi("name_enable",0);
   int num_threads = argi("nthreads",0);
   int earlyclip = argi("earlyclip",0);
   FILE *in;
   unsigned int channels;
   long start_time = (long)time(0);
   flam3_img_comments fpc;
   stat_struct stats;
   char numiter_string[64];
   char badval_string[64];
   char rtime_string[64];
   randctx rc;

   if (1 != argc) {
     docstring();
     exit(0);
   }

   /* Init random number generators */
   rand_seed(&rc);

   /* Set the number of threads */
   if (num_threads==0) {
      num_threads = flam3_count_nthreads();
      if (verbose > 1)
         fprintf(stderr,"Automatically detected %d core(s)...\n",num_threads);
   } else{
      if (verbose)
         fprintf(stderr,"Manually specified %d thread(s)...\n",num_threads);
   }


   channels = strcmp(format, "png") ? 3 : 4;

   /* Check for 16-bit-per-channel processing */
   if ( (16 == bpc) && (strcmp(format,"png") != 0)) {
	fprintf(stderr,"Support for 16 bpc images is only present for the png format.\n");
	exit(1);
   } else if (bpc != 8 && bpc != 16) {
	fprintf(stderr,"Unexpected bpc specified (%d)\n",bpc);
	exit(1);
   }
   
   if (pixel_aspect <= 0.0) {
     fprintf(stderr, "pixel aspect ratio must be positive, not %g.\n",
        pixel_aspect);
     exit(1);
   }

   if (inf)
     in = fopen(inf, "rb");
   else
     in = stdin;
   if (NULL == in) {
     perror(inf);
     exit(1);
   }

   cps = flam3_parse_from_file(in, inf, flam3_defaults_on, &ncps, &rc);
   if (NULL == cps) {
     fprintf(stderr,"error reading genomes from file\n");
     exit(1);
   }
   
   if (inf)
      fclose(in);

   for (i = 0; i < ncps; i++) {
      /* Force ntemporal_samples to 1 for -render */
      cps[i].ntemporal_samples = 1;
      cps[i].sample_density *= qs;
      cps[i].height = (int)(cps[i].height * ss);
      cps[i].width = (int)(cps[i].width * ss);
      cps[i].pixels_per_unit *= ss;
      if (cps[i].height<=0 || cps[i].width<=0) {
         fprintf(stderr,"output image has dimension <=0, aborting.\n");
         exit(1);
      }
   }

   if (out && (ncps > 1)) {
      fprintf(stderr, "hqi-flame: warning: writing multiple images "
      "to one file.  all but last will be lost.\n");
   }


   for (i = 0; i < ncps; i++) {
      if (verbose && ncps > 1) {
         fprintf(stderr, "flame = %d/%d ", i+1, ncps);
      }

//      f.temporal_filter_radius = 0.0;
      f.genomes = &cps[i];
      f.ngenomes = 1;
      f.verbose = verbose;
      f.bits = bits;
      f.time = 0.0;
      f.pixel_aspect_ratio = pixel_aspect;
      f.progress = 0;//print_progress;
      f.nthreads = num_threads;
      f.earlyclip = earlyclip;
      f.sub_batch_size = sub_batch_size;
      
      if (16==bpc)
         f.bytes_per_channel = 2;
      else
         f.bytes_per_channel = 1;
         

      this_size = (size_t)channels * (size_t)cps[i].width 
                  * (size_t)cps[i].height * f.bytes_per_channel;
	  image = (void *) calloc(this_size, sizeof(char));

      if (verbose && ncps > 1) {
         fprintf(stderr, "\n");
      }
      cps[i].ntemporal_samples = 1;
      if (flam3_render(&f, image, flam3_field_both, channels, transparency, &stats)) {
         fprintf(stderr,"error rendering image: aborting.\n");
         exit(1);
      }

      if (NULL != out) {
         strcpy(fname,out);
      } else if (name_enable && cps[i].flame_name[0]>0) {
         sprintf(fname, "%s.%s",cps[i].flame_name,format);
      } else {
         sprintf(fname, "%s%05d.%s", prefix, i, format);
      }
      if (verbose) {
         fprintf(stderr, "writing %s...", fname);
      }
      fp = fopen(fname, "wb");
      if (NULL == fp) {
         perror(fname);
         exit(1);
      }

      /* Generate temp file with genome */
      fpc.genome = flam3_print_to_string(f.genomes);
      
      sprintf(badval_string,"%g",stats.badvals/(double)stats.num_iters);
      fpc.badvals = badval_string;
      sprintf(numiter_string,"%g",(double)stats.num_iters);
      fpc.numiters = numiter_string;
      sprintf(rtime_string,"%d",stats.render_seconds);
      fpc.rtime = rtime_string;

      write_png(fp, image, cps[i].width, cps[i].height, &fpc, f.bytes_per_channel);
      /* Free string */
      free(fpc.genome);

      fclose(fp);

      if (verbose) {
         fprintf(stderr, "done.\n");
      }
   }
   if (verbose && (ncps > 1)) {
      long total_time = (long)time(0) - start_time;

      if (total_time > 100)
         fprintf(stderr, "total time = %.1f minutes\n", total_time / 60.0);
      else
         fprintf(stderr, "total time = %ld seconds\n", total_time);
   }
   
   for (i=0;i<ncps;i++) {
   
      xmlFreeDoc(cps[i].edits);
      clear_cp(&cps[i],0);
   
   }
   free(cps);
   
   free(image);
   return 0;
}
