/*=========================================================/
/  makePNG module for PX425 2021 assignment 4.             /
/  Writes a PNG of the simulation grid                     /
/                                                          /
/  Original code created by N. Hine  - November 2021       /
/  (based on previous code by D. Quigley)                  /                              /
/=========================================================*/

#include <png.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "makePNG.h"

void writepng(char *filename, int **grid, int width, int height){
  /*==================================================================/
  / Function to plot discrete data for upto 7 values into a PNG image /
  / Values outside of this range will generate errors. Note that -1   /
  / is mapped to black.                                               /
  /-------------------------------------------------------------------/
  / Original B&W version by S. Brown - University of Warwick          /
  / Modified to use colour scales by G. Enstone - also Warwick        /
  /==================================================================*/
  int x,y;

  int pixels_per_block = 2; 
  
  /* Open output file */
  FILE *fp = fopen(filename, "wb");
  png_byte **col_pointers;
  if (!fp) bork("Couldn't open %s for writing.\n",filename);

  /* Set colour scale */
  const int palette_size    = 3;    
  const int num_colours     = 3;
  const int pixel_bit_depth = 8;
  png_color palette[palette_size];
  png_color colour;

  int colors[][3] =   { {71,46,230},    /* blue */
                        {0,0,0},
                       {230,230,46} } ; /* yellow */


  int icol;
  for (icol=0;icol<palette_size;icol++) {
    colour = (png_color){colors[icol][0],colors[icol][1],colors[icol][2]};
    palette[icol] = colour;
  }

  
  /* Set up the PNG file */
  png_structp png_ptr = png_create_write_struct
    (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) bork("Couldn't allocate PNG.\n");
  
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) bork("Couldn't allocate PNG info.\n");
  
  if (setjmp (png_jmpbuf (png_ptr))) bork("PNG error: long_jump\n");
  
  /* Set image attributes. */
  png_set_IHDR (png_ptr,
		info_ptr,
		width*pixels_per_block,
		height*pixels_per_block,
                pixel_bit_depth, 
		PNG_COLOR_TYPE_PALETTE,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
  
  /*png_set_packing(png_ptr); */
  png_set_PLTE(png_ptr, info_ptr, palette, palette_size);

  
  /* Initialise PNG columns */
  col_pointers = png_malloc (png_ptr, sizeof(png_byte *)*height*pixels_per_block);	

  int xp=0; /* x pixel */
  int yp=0; /* y pixel */
  
  int xs=0; /* index within a block */
  int ys=0;

  png_byte block_colour;

  for (xp=0;xp<width*pixels_per_block;xp++) {
    col_pointers[xp]=png_malloc(png_ptr, sizeof(png_byte)*width*pixels_per_block);
  }

  int corner_x;
  int corner_y;

  for (y=0;y<height;y++) {
  
    corner_y = y*pixels_per_block;       
  
    for (x=0;x<width;x++) {

      corner_x = x*pixels_per_block;
      
      /* Interior */
      block_colour = (png_byte)((grid[height-y-1][x]+1)%num_colours);      
      for (xs=0;xs<pixels_per_block;xs++){
	xp = corner_x + xs;
	for (ys=0;ys<pixels_per_block;ys++){
	  yp = corner_y + ys;
	  col_pointers[yp][xp] = block_colour;
	}
      }
      
    }
  }

  
  /* Output code */
  png_init_io(png_ptr, fp);
  png_set_rows(png_ptr, info_ptr, col_pointers);
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_PACKING, NULL);
  png_write_end(png_ptr,info_ptr);


  
  /* Tidy up */
  fclose(fp);
  for (y=0;y<height;y++){
    png_free(png_ptr, col_pointers[y]);
  }
  png_free(png_ptr, col_pointers);
  png_destroy_write_struct(&png_ptr,&info_ptr);
}


void bork(char *msg,...){
  va_list args;
  va_start(args,msg);
  vfprintf(stderr,msg,args);
  exit(EXIT_FAILURE);
}


