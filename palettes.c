/*
    Copyright (C) 1992-2009 Spotworks LLC
	Copyright (C) 2015 pucket contributors

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
#include <string.h>

#include "math.h"
#include "palettes.h"
#include "rect.h"

/*	Add color to the end of palette
 */
void palette_add (palette * const p, const double4 c) {
	++p->count;
	if (p->count > p->allocated || p->color == NULL) {
		double4 *newcolors;
		if (p->allocated == 0) {
			/* flam3’s default palette size */
			p->allocated = 256;
		} else {
			p->allocated *= 2;
		}
		int ret = posix_memalign ((void **) &newcolors, sizeof (*newcolors),
				sizeof (*newcolors) * p->allocated);
		assert (ret == 0 && newcolors != NULL);
		memcpy (newcolors, p->color,
				(p->count-1)*sizeof (*p->color));
		free (p->color);
		p->color = newcolors;
	}
	p->color[p->count-1] = c;
}

static void parse_palettes(xmlNode *node, palette_collection * const pc) {
	/* node can be NULL */
	assert (pc != NULL);

	while (node) {
		if (node->type == XML_ELEMENT_NODE &&
				!xmlStrcmp(node->name, (const xmlChar *)"palette")) {
			xmlAttrPtr attr = node->properties;

			++pc->count;
			pc->p = realloc (pc->p, pc->count * sizeof (*pc->p));
			assert (pc->p != NULL);

			palette * const p = &pc->p[pc->count-1];
			memset (p, 0, sizeof (*p));

			while (attr) {
				char *val = (char *) xmlGetProp(node, attr->name);
				if (!xmlStrcmp(attr->name, (const xmlChar *)"data")) {
					int c_idx = 0;

					while (val[c_idx] != '\0') {
						if (isspace( (int)val[c_idx])) {
							c_idx++;
						} else {
							unsigned int r,g,b;
							const int sscanf_ret = sscanf((char *)&(val[c_idx]),
									"00%2x%2x%2x",&r,&g,&b);
							assert (sscanf_ret == 3);

							palette_add (p, (double4) { (double) r/255.0,
									(double) g/255.0, (double) b/255.0, 1.0 });

							c_idx += 8;
						}
					}
				}

				xmlFree(val);
				attr = attr->next;
			}
		} else
			parse_palettes(node->children, pc);

		node = node->next;
	}
}

/*	Read palette collection from file
 */
bool palette_read_collection (const char * const filename,
		palette_collection * const pc) {
	xmlDocPtr doc;
	xmlNode *rootnode;

	if ((doc = xmlReadFile (filename, NULL, XML_PARSE_NONET)) == NULL) {
		return false;
	}

	memset (pc, 0, sizeof (*pc));

	rootnode = xmlDocGetRootElement(doc);
	parse_palettes (rootnode, pc);
	xmlFreeDoc(doc);

	return true;
}

const palette *palette_random (const palette_collection * const pc,
							   randctx * const rc) {
	size_t i = rand_mod (rc, pc->count);
	return &pc->p[i];
}

void palette_copy (const palette * restrict const src, palette * restrict const dest) {
	assert (src != NULL);
	assert (dest != NULL);

	*dest = *src;
	if (src->name != NULL) {
		dest->name = strdup (src->name);
	}
	if (src->color != NULL) {
		int ret = posix_memalign ((void **) &dest->color, sizeof (*dest->color),
				sizeof (*dest->color) * dest->count);
		assert (ret == 0 && dest->color != NULL);
		memcpy (dest->color, src->color, dest->count * sizeof (*dest->color));
	}
}

void palette_rotate_hue (palette * const p, double rot) {
	for (unsigned int i = 0; i < p->count; i++) {
		double4 rgb = p->color[i];
		double4 hsv = rgb2hsv(rgb);
		hsv[0] += rot * 6.0;
		rgb = hsv2rgb(hsv);
		p->color[i] = rgb;
	}
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

static int random_xform(flam3_genome *g, int excluded, randctx * const rc) {
   int ntries = 0;
   while (ntries++ < 100) {
      int i = rand_mod(rc, g->num_xforms);
      if (g->xform[i].density > 0.0 && i != excluded)
         return i;
   }
   return -1;
}

static double try_colors(flam3_genome *g, unsigned int color_resolution,
		double timelimit) {
	assert (g != NULL);
	assert (color_resolution > 0);

	int *hist;
	unsigned int res = color_resolution, res3 = res * res * res;
	unsigned char *image, *p;
	flam3_genome saved;
	int pixtotal;

	memset(&saved, 0, sizeof(flam3_genome));

	flam3_copy(&saved, g);

	/* Scale the image so that the total number of pixels is ~10000 */
	pixtotal = g->width * g->height;
	const double scalar = sqrt( 10000.0 / (double)pixtotal);
	g->width *= scalar;
	g->height *= scalar;
	g->pixels_per_unit *= scalar;

	const unsigned int bytes_per_channel=1;
	const unsigned int channels = 4;

	image = (unsigned char *) calloc(g->width * g->height, bytes_per_channel * channels);

	bucket bucket;
	bucket_init (&bucket, (uint2) { g->width, g->height });
	render_bucket (g, &bucket, timelimit);
	render_image (g, &bucket, image, bytes_per_channel);

	hist = calloc(sizeof(int), res3);
	p = image;
	for (unsigned int i = 0; i < g->height * g->width; i++) {
		hist[(p[0] * res / 256) +
			(p[1] * res / 256) * res +
			(p[2] * res / 256) * res * res]++;
		p += channels;
	}

	unsigned int hits = 0;
	for (unsigned int i = 0; i < res3; i++) {
		if (hist[i]) hits++;
	}

	free(hist);
	free(image);

	g->width = saved.width;
	g->height = saved.height;
	g->pixels_per_unit = saved.pixels_per_unit;

	/* Free xform storage */
	clear_cp(&saved,flam3_defaults_on);

	return (double) hits / (double) res3;
}

static void change_colors(flam3_genome *g, bool change_palette,
		const palette_collection * const pc, randctx * const rc) {
   int i;
   int x0, x1;
   if (change_palette) {
      g->hue_rotation = 0.0;
	  const palette * const p = palette_random (pc, rc);
	  assert (p != NULL);
	  palette_copy (p, &g->palette);
   }
   for (i = 0; i < g->num_xforms; i++) {
      g->xform[i].color = rand_d01(rc);
   }
   x0 = random_xform(g, -1, rc);
   x1 = random_xform(g, x0, rc);
   if (x0 >= 0 && rand_bool(rc)) g->xform[x0].color = 0.0;
   if (x1 >= 0 && rand_bool(rc)) g->xform[x1].color = 1.0;
}

void flam3_improve_colors(flam3_genome *g, unsigned int ntries,
		bool change_palette,
		unsigned int color_resolution, double timelimit,
		const palette_collection * const pc,
		randctx * const rc) {
	const double trytime = timelimit/(double) ntries;
	flam3_genome best_genome;

	memset(&best_genome, 0, sizeof(flam3_genome));

	double best = try_colors(g, color_resolution, trytime);
	assert (best >= 0.0);

	flam3_copy(&best_genome,g);
	for (unsigned int i = 0; i < ntries; i++) {
		change_colors(g, change_palette, pc, rc);
		const double b = try_colors(g, color_resolution, trytime);
		assert (b >= 0.0);
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

