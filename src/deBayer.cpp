#include "pgm_io.h"

/* separate img_bayer interleaved Bayer matrix into red, blue, and green pixel planes:
   rgrgr <- row 0;  unread pixels in lower case
   gbgbg
   rgRGR
   gbGBG
   rgRGR
   ... where each input pixel has only one color component.
   Note that first and last columns and rows are ignored...???
   .. with the first red value  from row 2, column 2,
      the first blue pixel from row 3, column 3
      and first green pixels from r2, column 3 and row 3, column 2.
   .. then red is written to imgR row 1 column 1, leaving a black pixel border.
   Green pixel plane has a 2-pixel wide black border.
 */
void deBayer_char(image_char &img_bayer, image_char &imgR, image_char &imgG, image_char &imgB)
{
	printf("de-Bayer into separate red, green, blue planes... ");
	int wiRB = imgR->xsize, heRB = imgR->ysize;
	unsigned char *red, *blue, *green;

	for (int i = 1; i < wiRB-1; i++)
	{
		for (int j = 1; j < heRB-1; j++)
		{
			red = img_bayer->data + i*2+j*2*img_bayer->xsize;
			imgR->data[i + j * wiRB] = *red;

			blue = img_bayer->data + i * 2 + 1 + (j * 2 + 1) * img_bayer->xsize;
			imgB->data[i + j * wiRB] = *blue;

			green = img_bayer->data + i*2+1+j*2*img_bayer->xsize;
			imgG->data[i*2+1+j*2*imgG->xsize] = *green;
			green = img_bayer->data + i * 2 + (j * 2 + 1) * img_bayer->xsize;
			imgG->data[i * 2 + (j * 2 + 1) * imgG->xsize] = *green;

			imgG->data[i * 2 + j * 2 * imgG->xsize] = 
			(
			 2 +  img_bayer->data[i * 2 + 1 + j * 2 * img_bayer->xsize]
			 + img_bayer->data[i * 2 - 1 + j * 2 * img_bayer->xsize]
			 + img_bayer->data[i * 2 + (j * 2 + 1) * img_bayer->xsize]
			 + img_bayer->data[i * 2 + (j * 2 - 1) * img_bayer->xsize]
			) >> 2;
			imgG->data[i * 2 + 1 + (j * 2 + 1) * imgG->xsize] = 
			(
			 2 + img_bayer->data[i * 2 + 1 + j * 2 * img_bayer->xsize]
		 	 + img_bayer->data[i * 2 + (j * 2 + 1) * img_bayer->xsize]
			 + img_bayer->data[i * 2 + 2 + (j * 2 + 1) * img_bayer->xsize]
			 + img_bayer->data[i * 2 + 1 + (j * 2 + 2) * img_bayer->xsize]
			) >> 2;
		}
	}
	printf("done\n");
}

void deBayer(image_char &img_bayer, image_double &imgR, image_double &imgG, image_double &imgB)
{
	printf("de-Bayer into separate red, green, blue planes... ");
	int wiRB = imgR->xsize, heRB = imgR->ysize;
	unsigned char *red, *blue, *green;

	for (int i = 1; i < wiRB-1; i++)
	{
		for (int j = 1; j < heRB-1; j++)
		{
			red = img_bayer->data + i*2+j*2*img_bayer->xsize;
			imgR->data[i + j * wiRB] = *red;

			blue = img_bayer->data + i * 2 + 1 + (j * 2 + 1) * img_bayer->xsize;
			imgB->data[i + j * wiRB] = *blue;

			green = img_bayer->data + i*2+1+j*2*img_bayer->xsize;
			imgG->data[i*2+1+j*2*imgG->xsize] = *green;
			green = img_bayer->data + i * 2 + (j * 2 + 1) * img_bayer->xsize;
			imgG->data[i * 2 + (j * 2 + 1) * imgG->xsize] = *green;

			imgG->data[i * 2 + j * 2 * imgG->xsize] = 0.25 *
			(
			   img_bayer->data[i * 2 + 1 + j * 2 * img_bayer->xsize]
			 + img_bayer->data[i * 2 - 1 + j * 2 * img_bayer->xsize]
			 + img_bayer->data[i * 2 + (j * 2 + 1) * img_bayer->xsize]
			 + img_bayer->data[i * 2 + (j * 2 - 1) * img_bayer->xsize]
			);
			imgG->data[i * 2 + 1 + (j * 2 + 1) * imgG->xsize] = 0.25 *
			(
			   img_bayer->data[i * 2 + 1 + j * 2 * img_bayer->xsize]
		 	 + img_bayer->data[i * 2 + (j * 2 + 1) * img_bayer->xsize]
			 + img_bayer->data[i * 2 + 2 + (j * 2 + 1) * img_bayer->xsize]
			 + img_bayer->data[i * 2 + 1 + (j * 2 + 2) * img_bayer->xsize]
			);
		}
	}
	printf("done\n");
}
