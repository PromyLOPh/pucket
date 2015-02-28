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

#include "private.h"
#include "palettes.h"

lib_palette *the_palettes = NULL;
int npalettes;

static void parse_palettes(xmlNode *node) {
   xmlAttrPtr attr;
   char *val;
   lib_palette *pal;
   int hex_error;
    
   while (node) {
      if (node->type == XML_ELEMENT_NODE && !xmlStrcmp(node->name, (const xmlChar *)"palette")) {
         attr = node->properties;
	      pal = &the_palettes[npalettes];
	      memset(pal, 0, sizeof(lib_palette));

         while (attr) {
            val = (char *) xmlGetProp(node, attr->name);
            if (!xmlStrcmp(attr->name, (const xmlChar *)"data")) {
               int count = 256;
               int c_idx = 0;
               int r,g,b;
               int col_count = 0;
               int sscanf_ret;
               
               c_idx=0;
               col_count = 0;
               hex_error = 0;
               
               do {
                  sscanf_ret = sscanf((char *)&(val[c_idx]),"00%2x%2x%2x",&r,&g,&b);
                  if (sscanf_ret != 3) {
                     fprintf(stderr,"error:  problem reading hexadecimal color data '%8s'\n",&val[c_idx]);
                     hex_error = 1;
                     break;
                  }
			         c_idx += 8;
                  while (isspace( (int)val[c_idx]))
                     c_idx++;

                  pal->colors[col_count][0] = (double)r/255.0;
                  pal->colors[col_count][1] = (double)g/255.0;
                  pal->colors[col_count][2] = (double)b/255.0;
                  
                  col_count++;
               } while (col_count<count);
            } else if (!xmlStrcmp(attr->name, (const xmlChar *)"number")) {
               pal->number = atoi(val);
            } else if (!xmlStrcmp(attr->name, (const xmlChar *)"name")) {
               strncpy(pal->name, val, flam3_name_len);
               pal->name[flam3_name_len-1] = 0;
            }
            
            xmlFree(val);
            attr = attr->next;
         }
         
         if (hex_error == 0) {
            npalettes++;
            the_palettes = realloc(the_palettes, (1 + npalettes) * sizeof(lib_palette));
         }
      } else
         parse_palettes(node->children);
         
	   node = node->next;
   }
}

static int init_palettes(char *filename) {
    FILE *fp;
    xmlDocPtr doc;
    xmlNode *rootnode;
    int i, c, slen = 5000;
    char *s;

    fp = fopen(filename, "rb");
    if (NULL == fp) {
        fprintf(stderr, "flam3: could not open palette file ");
        perror(filename);
        return(-1);
    }
   
   /* Incrementally read XML file into a string */
   s = malloc(slen);
   i = 0;
   do {
      c = getc(fp);
      if (EOF == c) {
         if (ferror(fp)) {
            perror(filename);
            return(-1);
         }
         break;
      }
      s[i++] = c;
      if (i == slen-1) {
         slen *= 2;
         s = realloc(s, slen);
      }
   } while (1);

   fclose(fp);
   s[i] = 0;
   
   doc = xmlReadMemory(s, (int)strlen(s), filename, NULL, XML_PARSE_NONET);
   if (NULL == doc) {
       fprintf(stderr, "error parsing %s (%s).\n", filename, s);
       return(-1);
   }
   rootnode = xmlDocGetRootElement(doc);
   the_palettes = malloc(sizeof(lib_palette));
   npalettes = 0;
   parse_palettes(rootnode);
   xmlFreeDoc(doc);
   
   free(s);
   xmlCleanupParser();
   return(1);
}

int flam3_get_palette(int n, flam3_palette c, double hue_rotation, randctx * const rc) {
   int cmap_len = 256;
   int idx, i, j, rcode;

   // set palette to all white in case there are problems
   for (i = 0; i < cmap_len; i++) {
   	c[i].index = i;
      for (j = 0; j < 4; j++)
         c[i].color[j] = 1.0;
   }

   if (NULL == the_palettes) {   
      rcode = init_palettes("flam3-palettes.xml");
      if (rcode<0) {
         fprintf(stderr,"error reading xml palette file, setting to all white\n");
         return(-1);
      }
   }
   
   if (flam3_palette_random == n)
      n = the_palettes[rand_mod(rc, npalettes)].number;

   for (idx = 0; idx < npalettes; idx++) {
      
      if (n == the_palettes[idx].number) {
      	/* Loop over elements of cmap */
	      for (i = 0; i < cmap_len; i++) {
            int ii = (i * 256) / cmap_len;
            double4 rgb, hsv;
            
             rgb = (double4) {
			       the_palettes[idx].colors[ii][0],
				   the_palettes[idx].colors[ii][1],
				   the_palettes[idx].colors[ii][2],
				   1.0,
				   };
	       
	         hsv = rgb2hsv(rgb);
	         hsv[0] += hue_rotation * 6.0;
	         rgb = hsv2rgb(hsv);
	         
	         c[i].index = i;
	       
		     c[i].color = rgb;
		         
		     c[i].color[3] = 1.0;
         }
	   
	      return n;
      }
   }
   
   fprintf(stderr, "warning: palette number %d not found, using white.\n", n);

   return(-1);
}

/* rgb 0 - 1,
   h 0 - 6, s 0 - 1, v 0 - 1 */
double4 rgb2hsv(double4 rgb)
 {
  double rd, gd, bd, h, s, v, max, min, del, rc, gc, bc;

  rd = rgb[0];
  gd = rgb[1];
  bd = rgb[2];

  /* compute maximum of rd,gd,bd */
  if (rd>=gd) { if (rd>=bd) max = rd;  else max = bd; }
         else { if (gd>=bd) max = gd;  else max = bd; }

  /* compute minimum of rd,gd,bd */
  if (rd<=gd) { if (rd<=bd) min = rd;  else min = bd; }
         else { if (gd<=bd) min = gd;  else min = bd; }

  del = max - min;
  v = max;
  if (max != 0.0) s = (del) / max;
             else s = 0.0;

  h = 0;
  if (s != 0.0) {
    rc = (max - rd) / del;
    gc = (max - gd) / del;
    bc = (max - bd) / del;

    if      (rd==max) h = bc - gc;
    else if (gd==max) h = 2 + rc - bc;
    else if (bd==max) h = 4 + gc - rc;

    if (h<0) h += 6;
  }

  return (double4) { h, s, v, rgb[3] };
}


/* h 0 - 6, s 0 - 1, v 0 - 1
   rgb 0 - 1 */
double4 hsv2rgb(double4 hsv)
{
   double h = hsv[0], s = hsv[1], v = hsv[2];
  int    j;
  double rd, gd, bd;
  double f, p, q, t;

   while (h >= 6.0) h = h - 6.0;
   while (h <  0.0) h = h + 6.0;
   j = (int) floor(h);
   f = h - j;
   p = v * (1-s);
   q = v * (1 - (s*f));
   t = v * (1 - (s*(1 - f)));
   
   switch (j) {
    case 0:  rd = v;  gd = t;  bd = p;  break;
    case 1:  rd = q;  gd = v;  bd = p;  break;
    case 2:  rd = p;  gd = v;  bd = t;  break;
    case 3:  rd = p;  gd = q;  bd = v;  break;
    case 4:  rd = t;  gd = p;  bd = v;  break;
    case 5:  rd = v;  gd = p;  bd = q;  break;
    default: rd = v;  gd = t;  bd = p;  break;
   }

   return (double4) { rd, gd, bd, hsv[3] };
}

double flam3_calc_alpha(double density, double gamma, double linrange) {

   double dnorm = density;
   double funcval = pow(linrange, gamma);
   double frac,alpha;
   
   if (dnorm>0) {
      if (dnorm < linrange) {
         frac = dnorm/linrange;
         alpha = (1.0-frac) * dnorm * (funcval / linrange) + frac * pow(dnorm,gamma);
      } else
         alpha = pow(dnorm,gamma);
   } else
      alpha = 0;
      
   return(alpha);
}
          
double4 flam3_calc_newrgb(double4 cbuf, double ls, double highpow) {

   int rgbi;
   double newls,lsratio;
   double a, maxa=-1.0, maxc=0;
   double adjhlp;
   
   if (ls==0.0 || (cbuf[0]==0.0 && cbuf[1]==0.0 && cbuf[2]==0.0)) {
      return (double4) { 0, 0, 0, 0 };
   }
   
   /* Identify the most saturated channel */
   for (rgbi=0;rgbi<3;rgbi++) {
      a = ls * (cbuf[rgbi]);
      if (a>maxa) {
         maxa = a;
         maxc = cbuf[rgbi];
      }
   }
      
   /* If a channel is saturated and we have a non-negative highlight power */
   /* modify the color to prevent hue shift                                */
   if (maxa > 1.0 && highpow>=0.0) {
      newls = 1.0/maxc;
      lsratio = pow(newls/ls,highpow);

      /* Calculate the max-value color (ranged 0 - 1) */
	  double4 newrgb = newls*(cbuf);

      /* Reduce saturation by the lsratio */
      double4 newhsv = rgb2hsv(newrgb);
      newhsv[1] *= lsratio;

      return hsv2rgb(newhsv);
   } else {
      newls = 1.0/maxc;
      adjhlp = -highpow;
      if (adjhlp>1)
         adjhlp=1;
      if (maxa<=1.0)
         adjhlp=1.0;

	  /* Calculate the max-value color (ranged 0 - 1) interpolated with the old
	   * behaviour */
		
      return ((1.0-adjhlp)*newls + adjhlp*ls)*(cbuf);
   }
}

static int random_xform(flam3_genome *g, int excluded, randctx * const rc) {
   int ntries = 0;
   while (ntries++ < 100) {
      int i = rand_mod(rc, g->num_xforms);
      if (g->xform[i].density > 0.0 && i != excluded)
         return i;
   }
   return -1;
}


static double try_colors(flam3_genome *g, int color_resolution) {
    int *hist;
    int i, hits, res = color_resolution;
    int res3 = res * res * res;
    flam3_frame f;
    unsigned char *image, *p;
    flam3_genome saved;
    double scalar;
    int pixtotal;
    stat_struct stats;

    memset(&saved, 0, sizeof(flam3_genome));

    flam3_copy(&saved, g);

    g->sample_density = 1;
    
    /* Scale the image so that the total number of pixels is ~10000 */
    pixtotal = g->width * g->height;    
    scalar = sqrt( 10000.0 / (double)pixtotal);
    g->width *= scalar;
    g->height *= scalar;
    g->pixels_per_unit *= scalar;      
    
//    g->width = 100; // XXX keep aspect ratio
//    g->height = 100;
//    g->pixels_per_unit = 50;

    f.bytes_per_channel=1;
    f.verbose = 0;
    f.genomes = g;
    f.ngenomes = 1;
    f.earlyclip = 1;
    f.pixel_aspect_ratio = 1.0;
    f.progress = 0;
    f.nthreads = 1;
    f.sub_batch_size = 10000;
        
    image = (unsigned char *) calloc(g->width * g->height, 3);
    if (flam3_render(&f, image, flam3_field_both, &stats)) {
       fprintf(stderr,"Error rendering test image for trycolors.  Aborting.\n");
       return(-1);
    }

    hist = calloc(sizeof(int), res3);
    p = image;
    for (i = 0; i < g->height * g->width; i++) {
       hist[(p[0] * res / 256) +
            (p[1] * res / 256) * res +
            (p[2] * res / 256) * res * res]++;
       p += 3;
    }

    if (0) {
       int j, k;
       for (i = 0; i < res; i++) {
          fprintf(stderr, "\ni=%d: \n", i);
          for (j = 0; j < res; j++) {
             for (k = 0; k < res; k++) {
                fprintf(stderr, " %5d", hist[i * res * res + j * res + k]);
             }
             fprintf(stderr, "\n");
          }
       }
    }

    hits = 0;
    for (i = 0; i < res3; i++) {
       if (hist[i]) hits++;
    }

    free(hist);
    free(image);

    g->sample_density = saved.sample_density;
    g->width = saved.width;
    g->height = saved.height;
    g->pixels_per_unit = saved.pixels_per_unit;

    /* Free xform storage */
    clear_cp(&saved,flam3_defaults_on);

    return (double) (hits / res3);
}

static void change_colors(flam3_genome *g, int change_palette, randctx * const rc) {
   int i;
   int x0, x1;
   if (change_palette) {
      g->hue_rotation = 0.0;
      g->palette_index = flam3_get_palette(flam3_palette_random, g->palette, 0.0, rc);
      if (g->palette_index < 0)
         fprintf(stderr,"error retrieving random palette, setting to all white\n");
   }
   for (i = 0; i < g->num_xforms; i++) {
      g->xform[i].color = rand_d01(rc);
   }
   x0 = random_xform(g, -1, rc);
   x1 = random_xform(g, x0, rc);
   if (x0 >= 0 && rand_bool(rc)) g->xform[x0].color = 0.0;
   if (x1 >= 0 && rand_bool(rc)) g->xform[x1].color = 1.0;
}

void flam3_improve_colors(flam3_genome *g, int ntries, int change_palette, int color_resolution, randctx * const rc) {
   int i;
   double best, b;
   flam3_genome best_genome;

   memset(&best_genome, 0, sizeof(flam3_genome));

   best = try_colors(g, color_resolution);
   if (best<0) {
      fprintf(stderr,"error in try_colors, skipping flam3_improve_colors\n");
      return;
   }

   flam3_copy(&best_genome,g);
   for (i = 0; i < ntries; i++) {
      change_colors(g, change_palette, rc);
      b = try_colors(g, color_resolution);
      if (b < 0) {
         fprintf(stderr,"error in try_colors, aborting tries\n");
         break;
      }      
      if (b > best) {
         best = b;
         flam3_copy(&best_genome,g);
      }
   }

   flam3_copy(g,&best_genome);
   clear_cp(&best_genome,flam3_defaults_on);
}

