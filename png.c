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


#include <stdio.h>
#include <stdlib.h>
#include <png.h>
#include <setjmp.h>

#include "build/config.h"
#include "img.h"
#include "flam3.h"
#include "private.h"

#include <arpa/inet.h>

void write_png(FILE *file, void *image, int width, int height, flam3_img_comments *fpc, int bpc) {
  png_structp  png_ptr;
  png_infop    info_ptr;
  size_t i;
  unsigned short testbe = 1;
  void **rows = malloc(sizeof(void *) * height);

  for (i = 0; i < height; i++)
    rows[i] = (unsigned char *)image + i * width * 4 * bpc;
      
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
				    NULL, NULL, NULL);
  info_ptr = png_create_info_struct(png_ptr);

  if (setjmp(png_jmpbuf(png_ptr))) {
     fclose(file);
     png_destroy_write_struct(&png_ptr, &info_ptr);
     perror("writing file");
     return;
  }
  png_init_io(png_ptr, file);

  png_set_IHDR(png_ptr, info_ptr, width, height, 8*bpc,
	       PNG_COLOR_TYPE_RGBA,
	       PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_BASE,
	       PNG_FILTER_TYPE_BASE);
	       
  png_write_info(png_ptr, info_ptr);

  /* Must set this after the write_info */
  if (2==bpc && testbe != htons(testbe)) {
     png_set_swap(png_ptr);
  }

  png_write_image(png_ptr, (png_bytepp) rows);
  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  free(rows);

}

