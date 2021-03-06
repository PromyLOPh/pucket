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

#include "parser.h"
#include "interpolation.h"
#include <errno.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

static int parse_flame_element(xmlNode *flame_node, flam3_genome *loc_current_cp,
                        randctx * const rc);
static int parse_xform_xml(xmlNode *chld_node,flam3_xform *this_xform, int *num_xaos, 
                    flam3_chaos_entry **xaos, int numstd, int motionxf);

static int flam3_conversion_failed;

static int flam3_atoi(char *nstr) {

   /* Note that this is NOT thread-safe, but simplifies things significantly. */
   int res;
   char *endp;

   /* Reset errno */
   errno=0;

   /* Convert the string using strtol */
   res = strtol(nstr, &endp, 10);

   /* Check errno & return string */
   if (endp!=nstr+strlen(nstr)) {
      flam3_conversion_failed = 1;
      fprintf(stderr,"flam3_atoi : Error converting :%s: extra chars\n",nstr);
   }
   if (errno) {
      flam3_conversion_failed = 1;
      fprintf(stderr,"flam3_atoi : Error converting :%s:\n",nstr);
   }
   return(res);
}

static double flam3_atof(char *nstr) {

   /* Note that this is NOT thread-safe, but simplifies things significantly. */
   double res;
   char *endp;

   /* Reset errno */
   errno=0;

   /* Convert the string using strtod */
   res = strtod(nstr, &endp);
    
   /* Check errno & return string */
   if (endp!=nstr+strlen(nstr)) {
      flam3_conversion_failed = 1;
      fprintf(stderr,"flam3_atof: Error converting :%s: extra chars\n",nstr);
   }
   if (errno) {
      flam3_conversion_failed = 1;
      fprintf(stderr,"flam3_atof: Error converting :%s:\n",nstr);
   }
   return(res);
}

#define flam3_variation_none   (-1)

static int var2n(const char *s) {
   int i;
   
   for (i = 0; i < flam3_nvariations; i++)
      if (!strcmp(s, flam3_variation_names[i])) return i;
      
   return flam3_variation_none;
}

#if 0
int flam3_parse_hexformat_colors(char *colstr, flam3_genome *cp, int numcolors, int chan) {

   int c_idx=0;
   int col_count=0;
   int r,g,b,a;
   int sscanf_ret;
   char tmps[2];
   int skip = (int)fabs(chan);
   
   /* Strip whitespace prior to first color */
   while (isspace( (int)colstr[c_idx]))
      c_idx++;

   do {

      /* Parse an RGB triplet at a time... */
      a = 255;
      if (chan==3)
         sscanf_ret = sscanf(&(colstr[c_idx]),"%2x%2x%2x",&r,&g,&b);
      else if (chan==-4)
         sscanf_ret = sscanf(&(colstr[c_idx]),"00%2x%2x%2x",&r,&g,&b);
      else // chan==4
         sscanf_ret = sscanf(&(colstr[c_idx]),"%2x%2x%2x%2x",&r,&g,&b,&a);

      if ((chan!=4 && sscanf_ret != 3) || (chan==4 && sscanf_ret != 4)) {
         fprintf(stderr, "Error:  Problem reading hexadecimal color data.\n");
         return(1);
      }

      c_idx += 2*skip;

      while (isspace( (int)colstr[c_idx]))
         c_idx++;

      cp->palette[col_count].color[0] = r / 255.0;
      cp->palette[col_count].color[1] = g / 255.0;
      cp->palette[col_count].color[2] = b / 255.0;
      cp->palette[col_count].color[3] = a / 255.0;
      cp->palette[col_count].index = col_count;

      col_count++;

   } while (col_count<numcolors);
   
   if (sscanf(&(colstr[c_idx]),"%1s",tmps)>0) {
      fprintf(stderr,"error: extra data at end of hex color data '%s'\n",&(colstr[c_idx]));
      return(1);
   }
   
   return(0);
}
#endif

#if 0
int flam3_interp_missing_colors(flam3_genome *cp) {

    /* Check for a non-full palette */
    int minix,maxix;
    int colorli,colorri;
    int wrapmin,wrapmax;
    int intl, intr;
    int str,enr;
    int i,j,k;
    double prcr;
      
    for (i=0; i<256; i++) {
        if (cp->palette[i].index >= 0) {
            minix = i;
            break;
        }
    }
    
    if (i==256) {
        /* No colors.  Set all indices properly. */
        for (i=0;i<256;i++)
            cp->palette[i].index = i;            
        return(1);
    }
    
    wrapmin = minix + 256;
      
    for (i=255;i>=0;i--) {
        if (cp->palette[i].index >= 0) {
            maxix = i;
            break;
        }
    }
      
    wrapmax = maxix - 256;
      
    /* Loop over the indices looking for negs */
    i = 0;
    while(i<256) {

        if (cp->palette[i].index < 0) {
            /* Start of a range of negs */
            str = i;
            intl = i-1;
            colorli = intl;
            while (cp->palette[i].index<0 && i<256) {
                enr = i;
                intr = i+1;
                colorri = intr;
                i++;
            }

            if (intl==-1) {
                intl = wrapmax;
                colorli = maxix;
            }

            if (intr==256) {
                intr = wrapmin;
                colorri = minix;
            }
                
            for (j=str;j<=enr;j++) {

                prcr = (j-intl)/(double)(intr-intl);
                
                for (k=0;k<=3;k++)
                    cp->palette[j].color[k] = cp->palette[colorli].color[k] * (1.0-prcr) + cp->palette[colorri].color[k] * prcr;
                    
                cp->palette[j].index = j;
            }

            i = colorri+1;
        } else
            i ++;
    }

    return(0);
}
#endif


void scan_for_flame_nodes(xmlNode *cur_node, int default_flag, flam3_genome **all_cps, int *all_ncps, randctx * const rc) {

   xmlNode *this_node = NULL;
   flam3_genome loc_current_cp;
   flam3_genome *genome_storage = *all_cps; /* To simplify semantics */
   size_t f3_storage;
   int pfe_success;
   
   memset(&loc_current_cp,0,sizeof(flam3_genome));

   /* Loop over this level of elements */
   for (this_node=cur_node; this_node; this_node = this_node->next) {

      /* Check to see if this element is a <flame> element */
      if (this_node->type == XML_ELEMENT_NODE && !xmlStrcmp(this_node->name, (const xmlChar *)"flame")) {

         /* This is a flame element.  Parse it. */
         clear_cp(&loc_current_cp, default_flag);

         pfe_success = parse_flame_element(this_node,&loc_current_cp, rc);
         
         if (pfe_success>0) {
            fprintf(stderr,"error parsing flame element\n");
            all_cps = NULL; /* leaks memory but terminates */
            /* !!! free all_cp properly !!! */
            *all_ncps = 0;
            return;
         }

         /* Copy this cp into the array */
         f3_storage = (1+*all_ncps)*sizeof(flam3_genome);
         genome_storage = realloc(genome_storage, f3_storage);
         
         /* Must set value of pointer to new storage location */
         *all_cps = genome_storage;

         /* Clear out the realloc'd memory */
         memset(&(genome_storage[*all_ncps]),0,sizeof(flam3_genome));

         loc_current_cp.genome_index = *all_ncps;

         flam3_copy(&(genome_storage[*all_ncps]), &loc_current_cp);
         (*all_ncps) ++;         

      } else {
         /* Check all of the children of this element */
         scan_for_flame_nodes(this_node->children, default_flag, all_cps, all_ncps, rc);
      }
   }
   
   /* Clear the cp (frees allocated memory) */
   clear_cp(&loc_current_cp, default_flag);
   
}


static int parse_flame_element(xmlNode *flame_node, flam3_genome *loc_current_cp,
                        randctx * const rc) {
   flam3_genome *cp = loc_current_cp;
   xmlNode *chld_node;
   xmlAttrPtr att_ptr, cur_att;
   int solo_xform=-1;
   char *att_str;
   int num_std_xforms=-1;
   char tmps[2];
   int i;
   flam3_xform tmpcpy;
   flam3_chaos_entry *xaos=NULL;
   int num_xaos=0;

   /* Reset the conversion error flag */
   /* NOT threadsafe                  */
   flam3_conversion_failed=0;

   /* Store this flame element in the current cp */
   
   /* The top level element is a flame element. */
   /* Read the attributes of it and store them. */
   att_ptr = flame_node->properties;

   if (att_ptr==NULL) {
      fprintf(stderr, "Error : <flame> element has no attributes.\n");
      return(1);
   }

   memset(cp->flame_name,0,flam3_name_len+1);

   for (cur_att = att_ptr; cur_att; cur_att = cur_att->next) {

       att_str = (char *) xmlGetProp(flame_node,cur_att->name);
       
      /* Compare attribute names */
      if (!xmlStrcmp(cur_att->name, (const xmlChar *)"time")) {
         cp->time = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"interpolation")) {
     if (!strcmp("linear", att_str)) {
         cp->interpolation = flam3_interpolation_linear;
     } else if  (!strcmp("smooth", att_str)) {
         cp->interpolation = flam3_interpolation_smooth;
     } else {
         fprintf(stderr, "warning: unrecognized interpolation type %s.\n", att_str);
     }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"interpolation_space") ||
                 !xmlStrcmp(cur_att->name, (const xmlChar *)"interpolation_type")) {
      
         if (!strcmp("linear", att_str))
            cp->interpolation_type = flam3_inttype_linear;
         else if (!strcmp("log", att_str))
            cp->interpolation_type = flam3_inttype_log;
         else if (!strcmp("old", att_str))
            cp->interpolation_type = flam3_inttype_compat;
         else if (!strcmp("older", att_str))
            cp->interpolation_type = flam3_inttype_older;
         else
            fprintf(stderr,"warning: unrecognized interpolation_type %s.\n",att_str);
     
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"name")) {
         strncpy(cp->flame_name, att_str, flam3_name_len);
         i = (int)strlen(cp->flame_name)-1;
         while(i-->0) {
            if (isspace(cp->flame_name[i]))
               cp->flame_name[i] = '_';
         }

      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"size")) {
         if (sscanf(att_str, "%d %d%1s", &cp->width, &cp->height, tmps) != 2) {
            fprintf(stderr,"error: invalid size attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"center")) {
         if (sscanf(att_str, "%lf %lf%1s", &cp->center[0], &cp->center[1], tmps) != 2) {
            fprintf(stderr,"error: invalid center attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"scale")) {
         cp->pixels_per_unit = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rotate")) {
         cp->rotate = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"zoom")) {
         cp->zoom = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"palette_mode")) {
         if (!strcmp("step", att_str))
            cp->palette_mode = PALETTE_MODE_STEP;
         else if (!strcmp("linear", att_str))
            cp->palette_mode = PALETTE_MODE_LINEAR;
         else
            fprintf(stderr,"warning: unrecognized palette mode %s.  Using step.\n",att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"brightness")) {
         cp->brightness = flam3_atof(att_str);
/*      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"contrast")) {
         cp->contrast = flam3_atof(att_str);*/
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"gamma")) {
         cp->gamma = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"highlight_power")) {
         cp->highlight_power = atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"vibrancy")) {
         cp->vibrancy = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"hue")) {
         cp->hue_rotation = fmod(flam3_atof(att_str), 1.0);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"gamma_threshold")) {
         cp->gam_lin_thresh = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"soloxform")) {
         solo_xform = flam3_atof(att_str);
      }

      xmlFree(att_str);

   }

   /* Finished with flame attributes.  Now look at children of flame element. */
   for (chld_node=flame_node->children; chld_node; chld_node = chld_node->next) {

      /* Is this a color node? */
      if (!xmlStrcmp(chld_node->name, (const xmlChar *)"color")) {
         int index = -1;
         double r=0.0,g=0.0,b=0.0,a=0.0;

         /* Loop through the attributes of the color element */
         att_ptr = chld_node->properties;

         if (att_ptr==NULL) {
            fprintf(stderr,"Error:  No attributes for color element.\n");
            return(1);
         }

         for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

            att_str = (char *) xmlGetProp(chld_node,cur_att->name);

            a = 255.0;

            if (!xmlStrcmp(cur_att->name, (const xmlChar *)"index")) {
               index = flam3_atoi(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rgb")) {
               if (sscanf(att_str, "%lf %lf %lf%1s", &r, &g, &b, tmps) != 3) {
                  fprintf(stderr,"error: invalid rgb attribute '%s'\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rgba")) {
               if (sscanf(att_str, "%lf %lf %lf %lf%1s", &r, &g, &b, &a, tmps) != 4) {
                  fprintf(stderr,"error: invalid rgba attribute '%s'\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"a")) {
               if (sscanf(att_str, "%lf%1s", &a, tmps) != 1) {
                  fprintf(stderr,"error: invalid a attribute '%s'\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
            } else {
               fprintf(stderr,"Error:  Unknown color attribute '%s'\n",cur_att->name);
               xmlFree(att_str);
               return(1);
            }

            xmlFree(att_str);
         }
 
		assert (cp->palette.count == index);
		const double4 c = (double4) { r/255.0, g/255.0, b/255.0, a/255.0 };
		palette_add (&cp->palette, c);

#if 0
      } else if (!xmlStrcmp(chld_node->name, (const xmlChar *)"colors")) {

         int count = 0;

         /* Loop through the attributes of the colors element */
         att_ptr = chld_node->properties;

         if (att_ptr==NULL) {
            fprintf(stderr,"Error: No attributes for colors element.\n");
            return(1);
         }

         for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

            att_str = (char *) xmlGetProp(chld_node,cur_att->name);

            if (!xmlStrcmp(cur_att->name, (const xmlChar *)"count")) {
               count = flam3_atoi(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"data")) {
               if (flam3_parse_hexformat_colors(att_str, cp, count, -4) > 0) {
                  fprintf(stderr,"error parsing hexformatted colors\n");
                  xmlFree(att_str);
                  return(1);
               }
            } else {
               fprintf(stderr,"Error:  Unknown color attribute '%s'\n",cur_att->name);
               xmlFree(att_str);
               return(1);
            }

            xmlFree(att_str);
         }

      } else if (!xmlStrcmp(chld_node->name, (const xmlChar *)"palette")) {

         /* This could be either the old form of palette or the new form */
         /* Make sure BOTH are not specified, otherwise either are ok    */
         int numcolors=0;
         int numbytes=0;
         int old_format=0;
         int new_format=0;
         int index0, index1;
         double hue0, hue1;
         double blend = 0.5;
         index0 = index1 = flam3_palette_random;
         hue0 = hue1 = 0.0;

         /* Loop through the attributes of the palette element */
         att_ptr = chld_node->properties;

         if (att_ptr==NULL) {
            fprintf(stderr,"Error:  No attributes for palette element.\n");
            return(1);
         }

         for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

            att_str = (char *) xmlGetProp(chld_node,cur_att->name);

            if (!xmlStrcmp(cur_att->name, (const xmlChar *)"index0")) {
               old_format++;
               index0 = flam3_atoi(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"index1")) {
               old_format++;
               index1 = flam3_atoi(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"hue0")) {
               old_format++;
               hue0 = flam3_atof(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"hue1")) {
               old_format++;
               hue1 = flam3_atof(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"blend")) {
               old_format++;
               blend = flam3_atof(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"count")) {
               new_format++;
               numcolors = flam3_atoi(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"format")) {
               new_format++;
               if (!strcmp(att_str,"RGB"))
                  numbytes=3;
               else if (!strcmp(att_str,"RGBA"))
                  numbytes=4;
               else {
                  fprintf(stderr,"Error: Unrecognized palette format string (%s)\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
            } else {
               fprintf(stderr,"Error:  Unknown palette attribute '%s'\n",cur_att->name);
               xmlFree(att_str);
               return(1);
            }

            xmlFree(att_str);
         }

         /* Old or new format? */
         if (new_format>0 && old_format>0) {
            fprintf(stderr,"Error: mixing of old and new palette tag syntax not allowed.\n");
            return(1);
         }

         if (old_format>0)
            interpolate_cmap(cp->palette, blend, index0, hue0, index1, hue1, rc);
         else {

            char *pal_str;

            /* Read formatted string from contents of tag */

            pal_str = (char *) xmlNodeGetContent(chld_node);

            if (flam3_parse_hexformat_colors(pal_str, cp, numcolors, numbytes) > 0) {
               fprintf(stderr,"error reading hexformatted colors\n");
               xmlFree(pal_str);
               return(1);
            }

            xmlFree(pal_str);
         }
#endif
      } else if (!xmlStrcmp(chld_node->name, (const xmlChar *)"symmetry")) {

         int kind=0;
         int bef,aft;

         /* Loop through the attributes of the symmetry element */
         att_ptr = chld_node->properties;

         if (att_ptr==NULL) {
            fprintf(stderr,"Error:  No attributes for symmetry element.\n");
            return(1);
         }

         for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

            att_str = (char *) xmlGetProp(chld_node,cur_att->name);

            if (!xmlStrcmp(cur_att->name, (const xmlChar *)"kind")) {
               kind = flam3_atoi(att_str);
            } else {
               fprintf(stderr,"Error:  Unknown symmetry attribute '%s'\n",cur_att->name);
               xmlFree(att_str);
               return(1);
            }

            xmlFree(att_str);
         }

         bef = cp->num_xforms;
         flam3_add_symmetry(cp, kind);
         aft = cp->num_xforms;
         num_std_xforms += (aft-bef);
         

      } else if (!xmlStrcmp(chld_node->name, (const xmlChar *)"xform") ||
                  !xmlStrcmp(chld_node->name, (const xmlChar *)"finalxform")) {

         int xf = cp->num_xforms;
         
         if (!xmlStrcmp(chld_node->name, (const xmlChar *)"finalxform")) {

            if (cp->final_xform_index >=0) {
               fprintf(stderr,"Error:  Cannot specify more than one final xform.\n");
               return(1);
            }

            flam3_add_xforms(cp, 1, 0, 1);
            cp->xform[xf].var[0]=0.0;
            cp->final_xform_index = xf;
            /* Now, if present, the xform enable defaults to on */
            cp->final_xform_enable = 1;
            
         } else {

            /* Add one to the counter */
            flam3_add_xforms(cp, 1, 0, 0);
            
            /* If there was already a final xform, we have to change xf to point to the second to last xform */
            if (cp->final_xform_index>=0)
               xf--;
               
            cp->xform[xf].var[0]=0.0;
            num_std_xforms++;
                     
         }
         
         if (parse_xform_xml(chld_node, &(cp->xform[xf]), &num_xaos, &xaos, num_std_xforms, 0) != 0)
            return(1);
            
         if (cp->final_xform_index == xf && cp->xform[xf].density != 0.0) {
            fprintf(stderr,"Error: Final xforms should not have weight specified.\n");
            return(1);
         }
      }
   } /* Done parsing flame element. */
   
   num_std_xforms++;

   for (i=0;i<num_std_xforms;i++) {
      
      /* Adjust opacity with solo xform setting */
      if (solo_xform>=0 && i!=solo_xform)
         cp->xform[i].opacity = 0.0;

   }
      
   /* Set the chaos array entries with the values in the xaos list */
   for (i=0;i<num_xaos;i++)
      cp->chaos[xaos[i].from][xaos[i].to] = xaos[i].scalar;
      
   free(xaos);
   
   /* If there is a final xform in this cp, move it to the end of the list */   
   if (cp->final_xform_index >=0 && cp->final_xform_index != (cp->num_xforms-1)) {
      /* Make a copy of the final xform */
      tmpcpy = cp->xform[cp->final_xform_index];
      
      /* Move each other xform up one */
      for (i=cp->final_xform_index+1;i<cp->num_xforms;i++)
         cp->xform[i-1] = cp->xform[i];
         
      /* Put the final at the end */
      cp->xform[cp->num_xforms-1] = tmpcpy;
      
      cp->final_xform_index = cp->num_xforms - 1;
   }
   
   /* Check for bad parse */
   if (flam3_conversion_failed) {
      fprintf(stderr,"error: parsing a double or int attribute's value.\n");
      return(1);
   }
   
   return(0);

}

static int parse_xform_xml(xmlNode *chld_node,flam3_xform *this_xform, int *num_xaos, 
                    flam3_chaos_entry **xaos, int numstd, int motionxf) {

   xmlAttrPtr att_ptr, cur_att;
   char *att_str, *cpy;
   char tmps[2];
   int j,k;

   /* Loop through the attributes of the xform element */
   att_ptr = chld_node->properties;

   if (att_ptr==NULL) {
      fprintf(stderr,"Error: No attributes for element.\n");
      return(1);
   }

   for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

      att_str = (char *) xmlGetProp(chld_node,cur_att->name);

      cpy = att_str;
      if (!xmlStrcmp(cur_att->name, (const xmlChar *)"weight")) {
         this_xform->density = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"symmetry")) {
         /* Deprecated.  Set both color_speed and animate to this value. */
         this_xform->color_speed = (1.0-flam3_atof(att_str))/2.0;
         this_xform->animate = flam3_atof(att_str)>0 ? 0 : 1;
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"color_speed")) {
         this_xform->color_speed = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"animate")) {
         this_xform->animate = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"color")) {
         double tmpc1;
         this_xform->color = 0.0;
         /* Try two coords first */
         if (sscanf(att_str, "%lf %lf%1s", &this_xform->color, &tmpc1, tmps) != 2) {
            /* Try one color */
            if (sscanf(att_str, "%lf%1s", &this_xform->color,tmps) != 1) {
               fprintf(stderr,"Error: malformed xform color attribute '%s'\n",att_str);
               xmlFree(att_str);
               return(1);
            }
         }            
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"var1")) {
         for (j=0; j < flam3_nvariations; j++) {
            this_xform->var[j] = 0.0;
         }
         j = flam3_atoi(att_str);

         if (j < 0 || j >= flam3_nvariations) {
            fprintf(stderr,"Error:  Bad variation (%d)\n",j);
            xmlFree(att_str);
            return(1);
         }

         this_xform->var[j] = 1.0;
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"var")) {
         for (j=0; j < flam3_nvariations; j++) {
            char *cpy2;
            errno=0;
            this_xform->var[j] = strtod(cpy, &cpy2);
            if (errno != 0 || cpy==cpy2) {
               fprintf(stderr,"error: bad value in var attribute '%s'\n",att_str);
               xmlFree(att_str);
               return(1);
            }
            cpy=cpy2;
         }

         if (cpy != att_str+strlen(att_str)) {
            fprintf(stderr,"error: extra chars at the end of var attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"chaos")) {      
         /* Chaos scalars */
         
         char *tok;
         double scal;
         int toi=0;
         
         if (motionxf==1) {
            fprintf(stderr,"error: motion element cannot have a chaos attribute.\n");
            xmlFree(att_str);
            return(1);
         }
         
         /* The att string contains at least one value, delimited by a space */
         tok = strtok(cpy," ");
         while (tok!=NULL) {
            scal = flam3_atof(tok);
            
            /* Skip 1.0 entries */
            if (scal==1.0) {
               toi++;
               tok = strtok(NULL," ");
               continue;
            }
            
            /* Realloc the xaos list */
            *xaos = realloc((*xaos),(*num_xaos+1) * sizeof(flam3_chaos_entry));

            /* Populate the xaos list */
            (*xaos)[*num_xaos].from = numstd;
            (*xaos)[*num_xaos].to = toi;
            (*xaos)[*num_xaos].scalar = scal;
            toi++;
            (*num_xaos)++;
                                          
            /* Get the next token */
            tok = strtok(NULL," ");
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"plotmode")) {

         if (motionxf==1) {
            fprintf(stderr,"error: motion element cannot have a plotmode attribute.\n");
            xmlFree(att_str);
            return(1);
         }
         
         if (!strcmp("off", att_str))
            this_xform->opacity = 0.0;
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"opacity")) {
         this_xform->opacity = flam3_atof(att_str);
         
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"coefs")) {
         for (k=0; k<3; k++) {
            for (j=0; j<2; j++) {
               char *cpy2;
               errno = 0;
               this_xform->c[k][j] = strtod(cpy, &cpy2);
               if (errno != 0 || cpy==cpy2) {
                  fprintf(stderr,"error: bad value in coefs attribute '%s'\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
               cpy=cpy2;
            }
         }
         if (cpy != att_str+strlen(att_str)) {
            fprintf(stderr,"error: extra chars at the end of coefs attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"post")) {
         for (k = 0; k < 3; k++) {
            for (j = 0; j < 2; j++) {
               char *cpy2;
               errno = 0;
               this_xform->post[k][j] = strtod(cpy, &cpy2);
               if (errno != 0 || cpy==cpy2) {
                  fprintf(stderr,"error: bad value in post attribute '%s'\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
               cpy=cpy2;
            }
         }
         if (cpy != att_str+strlen(att_str)) {
            fprintf(stderr,"error: extra chars at end of post attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"blob_low")) {
         this_xform->blob_low = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"blob_high")) {
         this_xform->blob_high = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"blob_waves")) {
         this_xform->blob_waves = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pdj_a")) {
         this_xform->pdj_ac[0] = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pdj_b")) {
         this_xform->pdj_bd[0] = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pdj_c")) {
         this_xform->pdj_ac[1] = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pdj_d")) {
         this_xform->pdj_bd[1] = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"fan2_x")) {
         this_xform->fan2_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"fan2_y")) {
         this_xform->fan2_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rings2_val")) {
         this_xform->rings2_val = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"perspective_angle")) {
         this_xform->perspective_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"perspective_dist")) {
         this_xform->perspective_dist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"julian_power")) {
         this_xform->julian_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"julian_dist")) {
         this_xform->julian_dist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"juliascope_power")) {
         this_xform->juliascope_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"juliascope_dist")) {
         this_xform->juliascope_dist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"radial_blur_angle")) {
         this_xform->radial_blur_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pie_slices")) {
         this_xform->pie_slices = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pie_rotation")) {
         this_xform->pie_rotation = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pie_thickness")) {
         this_xform->pie_thickness = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"ngon_sides")) {
         this_xform->ngon_sides = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"ngon_power")) {
         this_xform->ngon_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"ngon_circle")) {
         this_xform->ngon_circle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"ngon_corners")) {
         this_xform->ngon_corners = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curl_c1")) {
         this_xform->curl_c1 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curl_c2")) {
         this_xform->curl_c2 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rectangles_x")) {
         this_xform->rectangles_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rectangles_y")) {
         this_xform->rectangles_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"disc2_rot")) {
         this_xform->disc2_rot = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"disc2_twist")) {
         this_xform->disc2_twist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_rnd")) {
         this_xform->super_shape_rnd = flam3_atof(att_str);
         /* Limit to [0,1] */
         if (this_xform->super_shape_rnd<0)
            this_xform->super_shape_rnd=0;
         else if (this_xform->super_shape_rnd>1)
            this_xform->super_shape_rnd=1;
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_m")) {
         this_xform->super_shape_m = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_n1")) {
         this_xform->super_shape_n1 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_n2")) {
         this_xform->super_shape_n2 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_n3")) {
         this_xform->super_shape_n3 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_holes")) {
         this_xform->super_shape_holes = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"flower_petals")) {
         this_xform->flower_petals = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"flower_holes")) {
         this_xform->flower_holes = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"conic_eccentricity")) {
         this_xform->conic_eccentricity = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"conic_holes")) {
         this_xform->conic_holes = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"parabola_height")) {
         this_xform->parabola_height = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"parabola_width")) {
         this_xform->parabola_width = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bent2_x")) {
         this_xform->bent2_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bent2_y")) {
         this_xform->bent2_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bipolar_shift")) {
         this_xform->bipolar_shift = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"cell_size")) {
         this_xform->cell_size = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"cpow_i")) {
         this_xform->cpow_i = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"cpow_r")) {
         this_xform->cpow_r = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"cpow_power")) {
         this_xform->cpow_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curve_xamp")) {
         this_xform->curve_xamp = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curve_yamp")) {
         this_xform->curve_yamp = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curve_xlength")) {
         this_xform->curve_xlength = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curve_ylength")) {
         this_xform->curve_ylength = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"escher_beta")) {
         this_xform->escher_beta = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_x")) {
         this_xform->lazysusan_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_y")) {
         this_xform->lazysusan_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_spin")) {
         this_xform->lazysusan_spin = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_space")) {
         this_xform->lazysusan_space = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_twist")) {
         this_xform->lazysusan_twist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"modulus_x")) {
         this_xform->modulus_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"modulus_y")) {
         this_xform->modulus_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscilloscope_separation")) {
         this_xform->oscope_separation = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscilloscope_frequency")) {
         this_xform->oscope_frequency = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscilloscope_amplitude")) {
         this_xform->oscope_amplitude = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscilloscope_damping")) {
         this_xform->oscope_damping = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscope_separation")) {
         this_xform->oscope_separation = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscope_frequency")) {
         this_xform->oscope_frequency = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscope_amplitude")) {
         this_xform->oscope_amplitude = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscope_damping")) {
         this_xform->oscope_damping = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"popcorn2_x")) {
         this_xform->popcorn2_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"popcorn2_y")) {
         this_xform->popcorn2_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"popcorn2_c")) {
         this_xform->popcorn2_c = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"separation_x")) {
         this_xform->separation_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"separation_xinside")) {
         this_xform->separation_xinside = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"separation_y")) {
         this_xform->separation_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"separation_yinside")) {
         this_xform->separation_yinside = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"split_xsize")) {
         this_xform->split_xsize = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"split_ysize")) {
         this_xform->split_ysize = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"splits_x")) {
         this_xform->splits_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"splits_y")) {
         this_xform->splits_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"stripes_space")) {
         this_xform->stripes_space = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"stripes_warp")) {
         this_xform->stripes_warp = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_angle")) {
         this_xform->wedge_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_hole")) {
         this_xform->wedge_hole = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_count")) {
         this_xform->wedge_count = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_swirl")) {
         this_xform->wedge_swirl = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_julia_angle")) {
         this_xform->wedge_julia_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_julia_count")) {
         this_xform->wedge_julia_count = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_julia_power")) {
         this_xform->wedge_julia_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_julia_dist")) {
         this_xform->wedge_julia_dist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_sph_angle")) {
         this_xform->wedge_sph_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_sph_hole")) {
         this_xform->wedge_sph_hole = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_sph_count")) {
         this_xform->wedge_sph_count = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_sph_swirl")) {
         this_xform->wedge_sph_swirl = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"whorl_inside")) {
         this_xform->whorl_inside = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"whorl_outside")) {
         this_xform->whorl_outside = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"waves2_scalex")) {
         this_xform->waves2_scalex = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"waves2_scaley")) {
         this_xform->waves2_scaley = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"waves2_freqx")) {
         this_xform->waves2_freqx = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"waves2_freqy")) {
         this_xform->waves2_freqy = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"auger_freq")) {
         this_xform->auger_freq = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"auger_weight")) {
         this_xform->auger_weight = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"auger_sym")) {
         this_xform->auger_sym = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"auger_scale")) {
         this_xform->auger_scale = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"flux_spread")) {
         this_xform->flux_spread = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"Re_A") || !xmlStrcmp(cur_att->name, (const xmlChar *)"mobius_re_a")) {
         this_xform->mobius_re_a = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"Re_B") || !xmlStrcmp(cur_att->name, (const xmlChar *)"mobius_re_b")) {
         this_xform->mobius_re_b = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"Re_C") || !xmlStrcmp(cur_att->name, (const xmlChar *)"mobius_re_c")) {
         this_xform->mobius_re_c = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"Re_D") || !xmlStrcmp(cur_att->name, (const xmlChar *)"mobius_re_d")) {
         this_xform->mobius_re_d = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"Im_A") || !xmlStrcmp(cur_att->name, (const xmlChar *)"mobius_im_a")) {
         this_xform->mobius_im_a = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"Im_B") || !xmlStrcmp(cur_att->name, (const xmlChar *)"mobius_im_b")) {
         this_xform->mobius_im_b = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"Im_C") || !xmlStrcmp(cur_att->name, (const xmlChar *)"mobius_im_c")) {
         this_xform->mobius_im_c = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"Im_D") || !xmlStrcmp(cur_att->name, (const xmlChar *)"mobius_im_d")) {
         this_xform->mobius_im_d = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"asteria_alpha")) {
         this_xform->asteria_alpha = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bcollide_num")) {
         this_xform->bcollide_num = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bcollide_a")) {
         this_xform->bcollide_a = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bmod_radius")) {
         this_xform->bmod_radius = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bmod_distance")) {
         this_xform->bmod_distance = flam3_atof(att_str);
      } else {
         int v = var2n((char *) cur_att->name);
         if (v != flam3_variation_none)
            this_xform->var[v] = flam3_atof(att_str);
         else
            fprintf(stderr,"Warning: unrecognized variation %s.  Ignoring.\n",(char *)cur_att->name);
      }


      xmlFree(att_str);
   }
   return(0);
}

