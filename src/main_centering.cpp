/* Chromatic aberration correction prototype.
Copyright (C) 2014 Victoria Rudakova <vicrucann@gmail.com>
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#define _CRT_SECURE_NO_DEPRECATE	// fopen(), fscanf() warnings; must be macro;  unavoidable VCR101
#include "centers.h"
#include "abberation.h"
#include "distortion.h"
#include "pgm_io.h"
#include "spline.h"
#include "correction.h"
#include "misc.h"

#include <ios>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

typedef unsigned int uint;

#include "gnuplot2file.h"

static void average_image(image_double &img_avg, image_double img) {
	int w = img->xsize;
	int h = img->ysize;
	new_image_double(img_avg, h, w);
	double *v0 = img->data, *v1 = v0;
	double *d = img_avg->data;
	for (double *dx = d + w; d < dx; d++)
		*d = *v1++;
	double *v2 = v1 + w;
	for (int v = 1; v < h-1; v++) {
		*d++ = *v1; 
		for (double *dx = d + w - 2; d < dx; d++) {
			double pix = *v0 + v0[1] + v0[2] + *v1 + v1[1] + v1[2] + *v2 + v2[1] + v2[2];
			*d = pix/9;
			v0++; v1++; v2++;
		}
		*d++ = v1[1]; v1+=2; v0+=2; v2+=2;
	}
	for (double *dx = d + w; d < dx; d++)
		*d = *v1++;
}

template <typename T>
int initial_tache(image_double I, vector<T>& h, T& rayon, bool color, T x, T y) {
	int COL_IMA = I->xsize;
	int LIG_IMA = I->ysize;
	int j = (int)x;
	int i = (int)y;
	int d = (int)(2*rayon);
	if (2*d+1 > LIG_IMA)
		d=(LIG_IMA-1)/2;
	if(2*d+1 > COL_IMA)
		d=(COL_IMA-1)/2;
	if(i<d)
		i=d+1;
	if(i>LIG_IMA-1-d)
		i=LIG_IMA-2-d;
	if(j<d)
		j=d+1;
	if(j>COL_IMA-1-d)
		j=COL_IMA-2-d;
	int val_haut=0;
	int val_bas=255;
	for (int k = -d; k <= d; k++) {
		for (int l = -d; l <= d; l++) {
			T lum = I->data[j+l+(i+k)*COL_IMA];
			if (lum > val_haut)
				val_haut = (int)lum;
			else if (lum<val_bas)
				val_bas=(int)lum;
		} }
	T seuil = 0;
	if (!color)
		seuil = val_bas + (val_haut - val_bas)/3 * 2;
	else if (color)
		seuil = val_bas + (val_haut - val_bas)/3;

	matrix<T> tab;
	tab.init(0, 2*d+1, 2*d+1);
	int label = 1;

	for (int k = -d+1; k <= d; k++){
		for(int l = -d+1; l <= d-1; l++){
			T lum= I->data[j+l+(i+k)*COL_IMA];
			if (lum < seuil) {
				int imin=l;
				int imax=l+1;
				while( (I->data[j+imax+(i+k)*COL_IMA] <= seuil) && (imax <= d-1) ){
					imax++; }
				int vallab=0;
				for(int m = imin; m <= imax; m++){
					if (tab(k+d-1,m+d) != 0){
						vallab = (int)tab(k+d-1,m+d);
					}
				}
				if (vallab == 0){
					vallab=label;
					label++;
				}
				for(int m = imin; m <= imax; m++){
					tab(k+d,m+d)=vallab;
				}
				l=imax;
			}
		}
	}
	matrix<T> bary;
	bary.init(0, label, 4);
	for(int k = -d; k <= d; k++){
		for(int l = -d; l <= d; l++){
			if(tab(k+d,l+d)!=0){
				T lum= I->data[j+l+(i+k)*COL_IMA];
				bary((int)tab(k+d,l+d),0)+=(255-lum)*(j+l);
				bary((int)tab(k+d,l+d),1)+=(255-lum)*(i+k);
				bary((int)tab(k+d,l+d),2)+=(255-lum);
				bary((int)tab(k+d,l+d),3)++;
			}
		}
	}
	int distmin=100;
	int labelmin=0;
	for(int k = 1; k < label; k++){
		T dist = std::sqrt( (bary(k,0)/bary(k,2)-x) * (bary(k,0)/bary(k,2)-x)+
			(bary(k,1)/bary(k,2)-y) * (bary(k,1)/bary(k,2)-y));
		if(dist < distmin && bary(k,3) > 25 ){ /* 25 = surface min*/
			distmin=(int)dist;
			labelmin=k;
		}
	}
	if(labelmin == 0) {
		printf("spot is too small (<=25 pixels)\n");
		return 1;
	}
	x=bary(labelmin,0)/bary(labelmin,2);
	y=bary(labelmin,1)/bary(labelmin,2);
	T sx2 = 0, sy2 = 0, sxy = 0, ss = 0;
	for(int k = -d; k <= d; k++){
		for(int l = -d; l <= d; l++){
			if(tab(k+d,l+d) == labelmin){
				sx2+=(j+l-x)*(j+l-x);
				sy2+=(i+k-y)*(i+k-y);
				sxy+=(i+k-y)*(j+l-x);
				ss++;
			}
		}
	}
	T lambda1 = ((sx2+sy2)/ss + std::sqrt(( (sx2+sy2)*(sx2+sy2)+4*(sxy*sxy-sx2*sy2)))/ss)/2.0;
	T lambda2 = ((sx2*sy2-sxy*sxy)/(ss*ss))/lambda1;
	rayon=std::sqrt(lambda1)*2;
	h[0] = 1.0/rayon; 								/* lambda1 */
	h[1] = (std::sqrt(lambda1/lambda2))/rayon;		/* lambda2 */
	h[2] = std::atan2(sx2/ss-lambda1, -sxy/ss);		/* alpha */
	h[3] = x;			/* tu */
	h[4] = y;			/* tv */
	h[5] = 0.25;		/* circle radius 1 */
	h[6] = -2.0;	 	/* pente */
	h[7] = 0.25;		/* rayon cercle 2 */
	h[8] = val_haut; 	/* val_haut */
	h[9] = val_bas; 	/* val_bas */
	h[10] = 1.0; 		/* position step */
	return 0;
}

template <typename T>
void trgtDataCalc(vector<T> &trgData, image_double img_avg, T cx, T cy, T delta) {
	int xbegin = (int)(cx + 0.5 - delta);
	int xend = (int)(cx + 0.5 + delta);
	int ybegin = (int)(cy + 0.5 - delta);
	int yend = (int)(cy + 0.5 + delta);
	int nerr = (yend-ybegin+1)*(xend-xbegin+1);
	trgData.init(0, nerr);
	int wi = img_avg->xsize;
	int he = img_avg->ysize;
	int idx = 0;
	for (int v = ybegin; v <= yend; v++) {
		for (int u = xbegin; u <= xend; u++) {
			if (u >= 1 && u <= wi-2 && v >= 1 && v <= he-2)
				trgData[idx] = img_avg->data[u+v*img_avg->xsize];
			else
				trgData[idx] = 255; // for 'tache noire'
			idx++;
		}
	}
}

template <typename T>
T centerLMA(image_double sub_img, bool clr, T& centerX, T& centerY)
{
	image_double img_avg; average_image(img_avg, sub_img);
	int w = sub_img->xsize;
	int h = sub_img->ysize;
	T cx = w/2, cy = h/2, radi = 0.4*w;
	vector<T> P(11);
	initial_tache(sub_img, P, radi, clr, cx, cy);
	vector<T> trgData; trgtDataCalc<T>(trgData, img_avg, P[3], P[4], radi*2);
	LMTacheC<T> ellipseLMA(img_avg, P[3], P[4], radi*2, clr, w, h);
	T rmse = ellipseLMA.minimize(P, trgData, 0.001);
	free_image_double(img_avg);
	//T lambda1 = P[0]; T lambda2 = P[1]; T theta = P[2];
	centerX = P[3];
	centerY = P[4];
	return rmse;
}

template <typename T>
void takeSubImg(image_double &img, image_char IMG, T cx, T cy, T radi, int& x0, int& y0)
{
	int size = (int)(2.5 * radi);
	x0 = (int)(cx - 0.5*size);
	y0 = (int)(cy - 0.5*size);
	int i, j, jx = IMG->ysize - y0, ix = IMG->xsize - x0;
	if (size < jx)
		jx = size;
	if (size < ix)
		ix = size;
	new_image_double(img, size, size);
	double *d = img->data;
	for (j = 0; j < -y0; j++) 
		for (i = 0; i < size; i++)
			*d++ = 255;
	for (; j <  jx; j++) {
		for (i = 0; i < -x0; i++)
			*d++ = 255;
		unsigned char *r = IMG->data + x0+i+(y0+j)*IMG->xsize;
		for (; i < ix; i++)
			*d++ = *r++;
		for (; i < size; i++)
			*d++ = 255;
	}
	for (; j < size; j++)
		for (i = 0; i < size; i++)
			*d++ = 255;
}

template <typename T>
void takeSubImg(image_double &img, image_double IMG, T cx, T cy, T radi, int& x0, int& y0)
{
	int size = (int)(2.5 * radi);
	x0 = (int)(cx - 0.5*size);
	y0 = (int)(cy - 0.5*size);
	int i, j, jx = IMG->ysize - y0, ix = IMG->xsize - x0;
	if (size < jx)
		jx = size;
	if (size < ix)
		ix = size;
	new_image_double(img, size, size);
	double *d = img->data;
	for (j = 0; j < -y0; j++) 
		for (i = 0; i < size; i++)
			*d++ = 255;
	for (; j <  jx; j++) {
		for (i = 0; i < -x0; i++)
			*d++ = 255;
		double *r = IMG->data + x0+i+(y0+j)*IMG->xsize;
		for (; i < ix; i++)
			*d++ = *r++;
		for (; i < size; i++)
			*d++ = 255;
	}
	for (; j < size; j++)
		for (i = 0; i < size; i++)
			*d++ = 255;
}

template <typename T>
bool loadKeypts(const char* fname, std::vector<CCStats>& cc_green, std::vector<CCStats>& cc_red)
{
	std::ifstream f(fname);
	while( f.good() ) {
		std::string str;
		std::getline(f, str);
		if( f.good() ) {
			std::istringstream s(str);
			CCStats red{}, green{};
			s >> green.centerX >> green.centerY >> red.centerX >> red.centerY;
			if(!s.fail() )
			{
				cc_green.push_back(green);
				cc_red.push_back(red);
			}
		}
	}

	return true;
}

template <typename T>
void image_rotate_left(image_double &res, image_double img)
{
	new_image_double_ini(res, img->ysize, img->xsize, 0);
	int idx = 0;
	for (int x = img->xsize-1; x >= 0; x--)
	{
		for (int y = 0; y < img->ysize; y++)
		{
			res->data[idx] = img->data[x+y*img->xsize];
			idx++;
		}
	}
	return res;
}

template <typename T>
void image_rotate_right(image_double &res, image_double img)
{
	new_image_double_ini(res, img->ysize, img->xsize, 0);
	int idx = 0;
	for (int x = 0; x < img->xsize; x++)
	{
		for (int y = img->ysize-1; y >= 0; y--)
		{
			res->data[idx] = img->data[x+y*img->xsize];
			idx++;
		}
	}
	return res;
}

template <typename T>
void binarization(image_double &imgbiR, image_double &imgbiG, image_double &imgbiB,
	image_char &imgR, image_char &imgG, image_char &imgB,
	T threR, T threG, T threB)
{
	int wiRB = imgR->xsize;
	int heRB = imgR->ysize;
	int Gcols = imgG->xsize;
	int Grows = imgG->ysize;
	T red, blue, green;
	for (int k, i = 0; i < wiRB; i++) {
		for (int j = 0; j < heRB; j++) {
			red = imgR->data[k = i+j*wiRB];
			if (red <= threR) imgbiR->data[k] = 0;

			blue = imgB->data[k];
			if (blue <= threB) imgbiB->data[k] = 0;

			if (Gcols == wiRB && Grows == heRB) {
				green = imgG->data[k = i+j*wiRB];
				if (green <= threG) imgbiG->data[k] = 0;	}
			else {
				green = imgG->data[k = i*2+j*2*imgbiG->xsize];
				if (green <= threG) imgbiG->data[k] = 0;
				k++;
				green = imgG->data[k];
				if (green <= threG) imgbiG->data[k] = 0;
				green = imgG->data[k = i*2+(j*2+1)*imgbiG->xsize];
				if (green <= threG) imgbiG->data[k] = 0;
				k++;
				green = imgG->data[k];
				if (green <= threG) imgbiG->data[k] = 0; }
		}
	}
}

template <typename T>
void img_extremas(image_char &img, T &min, T &max) {
	size_t xsz = img->xsize, ysz = img->ysize;
	// start at [4][4]
	unsigned char* yp, *p = img->data + 4 + 4 * img->xsize;

	max = min = *p;
	// scan 8 fewer rows and columns
	ysz -= 8; xsz -= 8;
	for (yp = p + img->xsize * ysz; p < yp; p += 7)
		for (unsigned char *xp = p + xsz; p < xp; p++) {
			if (*p > max)
				max = *p;
			else if (*p < min)
				min = *p;
		}
	return;
}


template <typename T>
void circle_redefine(image_double &imgR, image_char &imgG, image_double &imgB,
	vector<T> &xR, vector<T> &yR, vector<T> &rR,
	vector<T> &xGr, vector<T> &yGr, vector<T> &rG,
	vector<T> &xB, vector<T> &yB, vector<T> &rB,
	vector<T> &xGb, vector<T> &yGb,
    T scale, bool clr, int ntaches, bool green = true)
{
	printf("\nLevenberg-Marquardt damped least squares center redefinition for R,B channels... \n");
	for (int i = 0; i < ntaches; i++) {
		int x0R = 0, y0R = 0, x0G = 0, y0G = 0, x0B = 0, y0B = 0;
		image_double sub_imgR{}; takeSubImg(sub_imgR, imgR, xR[i], yR[i], rR[i], x0R, y0R);
		image_double sub_imgB{}; takeSubImg(sub_imgB, imgB, xB[i], yB[i], rB[i], x0B, y0B);

		T cxR=0, cyR=0, cxB=0, cyB=0;
		centerLMA<T>(sub_imgR, clr, cxR, cyR);	// sets center cxR, cyR
		centerLMA<T>(sub_imgB, clr, cxB, cyB);

		free_image_double(sub_imgR);
		free_image_double(sub_imgB);

		xR[i] = scale*(x0R + cxR);
		yR[i] = scale*(y0R + cyR);
		xB[i] = scale*(x0B + cxB);
		yB[i] = scale*(y0B + cyB);

		if (green) {
			image_double sub_imgG{}; takeSubImg(sub_imgG, imgG, xGr[i], yGr[i], rG[i], x0G, y0G);
			T cxG=0, cyG=0;
			centerLMA<T>(sub_imgG, clr, cxG, cyG);
			free_image_double(sub_imgG);
			xGr[i] = x0G + cxG;
			yGr[i] = y0G + cyG;
			xGb[i] = x0G + cxG;
			yGb[i] = y0G + cyG; }

		double percent = ((double)i / (double)ntaches)*100;
		if (!(i % (int)(0.2*ntaches)))
		std::cout << int(percent)+1 << "%" << std::flush;
//			printf("%i%c", (int)percent+1, '%');
		else if (!(i % (int)(0.04*ntaches)))
		 std::cout << "." << std::flush;
//			printf(".");
	}
}

template <typename T>
void keypnts_circle(image_char &imgR, image_char &imgG, image_char &imgB,
	vector<T> &xR, vector<T> &yR, vector<T> &rR,
	vector<T> &xGr, vector<T> &yGr, vector<T> &rG,
	vector<T> &xB, vector<T> &yB, vector<T> &rB,
	vector<T> &xGb, vector<T> &yGb,
	T scale, bool clr)
{
	int wiRB = imgR->xsize, heRB = imgR->ysize;
	int Gcols = (int)(wiRB*scale), Grows = (int)(heRB*scale);

	T maxR = 0, maxG = 0, maxB = 0;
	T minR = 255, minG = 255, minB = 255;
	img_extremas(imgR, minR, maxR);
	img_extremas(imgG, minG, maxG);
	img_extremas(imgB, minB, maxB);
	T threR = 0.5 * (maxR-minR);
	T threG = 0.54* (maxG-minG);
	T threB = 0.4 * (maxB-minB);

	image_double imgbiR; new_image_double_ini(imgbiR, wiRB, heRB, 255);
	image_double imgbiG; new_image_double_ini(imgbiG, Gcols, Grows, 255);
	image_double imgbiB; new_image_double_ini(imgbiB, wiRB, heRB, 255);

	binarization(imgbiR, imgbiG, imgbiB, imgR, imgG, imgB, threR, threG, threB);
	//write_pgm_image_double(imgbiB, FOLDER "b.pgm");

	printf("\nfinding connected components:");
	std::vector<CCStats> ccstatsR, ccstatsG, ccstatsB;
	CC(ccstatsR, imgbiR, 'R'); printf(" number = %u ", (uint)ccstatsR.size());
	CC(ccstatsG, imgbiG, 'G'); printf(" number = %u ", (uint)ccstatsG.size());
	CC(ccstatsB, imgbiB, 'B'); printf(" number = %u ", (uint)ccstatsB.size());

	assert(ccstatsR.size() == ccstatsG.size() && ccstatsG.size() == ccstatsB.size());
	printf("\nRGB spot centers initialization is done;  begin matching... ");

	int ntaches = (int)ccstatsG.size();
	xR.init(1, ntaches); yR.init(1, ntaches); rR.init(1, ntaches);
	xB.init(1, ntaches); yB.init(1, ntaches); rB.init(1, ntaches);
	xGr.init(1, ntaches); yGr.init(1, ntaches); rG.init(1, ntaches);
	xGb.init(1, ntaches); yGb.init(1, ntaches);
	xR = -1; xB = -1;
	yR = -1; yB = -1;
	for (int i = 0; i < ntaches; i++) {
		T xg = ccstatsG[i].centerX;
		T yg = ccstatsG[i].centerY;
		int idxB = -1, idxR = findMatch(xg, yg, ccstatsR, scale);
		if (0 <= idxR)
		{
			idxB = findMatch(xg, yg, ccstatsB, scale);
			if (0 > idxB)
				printf("blue findMatch(%lf, %lf) fail\n", xg, yg);
		}
		else printf("red findMatch(%lf, %lf) fail\n", xg, yg);
		if (0 > idxB || 0 > idxR)
		{
			i--; ntaches--;	// forfeit
			continue;
		}

		xGr[i] = xg; yGr[i] = yg;
		xGb[i] = xg; yGb[i] = yg;
		rG[i] = 0.5*(ccstatsG[i].radius1+ccstatsG[i].radius2);
		xR[i] = ccstatsR[idxR].centerX;
		yR[i] = ccstatsR[idxR].centerY;
		rR[i] = 0.5*(ccstatsR[idxR].radius1+ccstatsR[idxR].radius2);
		xB[i] = ccstatsB[idxB].centerX;
		yB[i] = ccstatsB[idxB].centerY;
		rB[i] = 0.5*(ccstatsB[idxB].radius1+ccstatsB[idxB].radius2);
	}
	printf("done.\n");

//	circle_redefine(imgR, imgG, imgB, xR, yR, rR, xGr, yGr, rG, xB, yB, rB, xGb, yGb, scale, clr, ntaches);

	free_image_double(imgbiR);
	free_image_double(imgbiG);
	free_image_double(imgbiB);
}

template <typename T>
void print_RMSE(vector<T> &xR,  vector<T> &yR,	// red circle centers
				vector<T> &xGr, vector<T> &yGr,	// green circle centers
				vector<T> &xB,  vector<T> &yB,	// blue circle centers
				vector<T> &xGb, vector<T> &yGb)
{
	printf("\nData stats calculation... ");
	int ntachesr = xR.size();
	int ntachesb = xB.size();
	T RMSE_red_dist = 0, RMSE_blue_dist = 0;
	T mean_red_dist = 0, mean_blue_dist = 0, stddev_red_dist = 0, stddev_blue_dist = 0, dispR_dist = 0, dispB_dist = 0;
	for (int i = 0; i < ntachesr; i++)
	{
		T dispR = std::sqrt((xR[i]-xGr[i])*(xR[i]-xGr[i]) + (yR[i]-yGr[i])*(yR[i]-yGr[i]));
		if (dispR_dist < dispR) dispR_dist = dispR;
		RMSE_red_dist += dispR*dispR;
		mean_red_dist += dispR;
	}
	for (int i = 0; i < ntachesb; i++)
	{
		T dispB = std::sqrt((xB[i]-xGb[i])*(xB[i]-xGb[i]) + (yB[i]-yGb[i])*(yB[i]-yGb[i]));
		if (dispB_dist < dispB) dispB_dist = dispB;
		RMSE_blue_dist += dispB*dispB;
		mean_blue_dist += dispB;
	}
	mean_red_dist /= ntachesr;
	mean_blue_dist /= ntachesb;
	for (int i = 0; i < ntachesr; i++)
		stddev_red_dist += (std::sqrt((xR[i]-xGr[i])*(xR[i]-xGr[i]) + (yR[i]-yGr[i])*(yR[i]-yGr[i])) - mean_red_dist) *
			(std::sqrt((xR[i]-xGr[i])*(xR[i]-xGr[i]) + (yR[i]-yGr[i])*(yR[i]-yGr[i])) - mean_red_dist);
	for (int i = 0; i < ntachesb; i++)
		stddev_blue_dist += (std::sqrt((xB[i]-xGb[i])*(xB[i]-xGb[i]) + (yB[i]-yGb[i])*(yB[i]-yGb[i])) - mean_blue_dist) *
			(std::sqrt((xB[i]-xGb[i])*(xB[i]-xGb[i]) + (yB[i]-yGb[i])*(yB[i]-yGb[i])) - mean_blue_dist);
	printf(" done.\n");

	printf("RMSE R-G & B-G are:	%f	%f \n", std::sqrt(RMSE_red_dist/ntachesr), std::sqrt(RMSE_blue_dist/ntachesb));
	printf("mean R-G & B-G are:	%f	%f \n", mean_red_dist, mean_blue_dist);
	printf("stdd R-G & B-G are:	%f	%f \n", std::sqrt(stddev_red_dist/ntachesr), std::sqrt(stddev_blue_dist/ntachesb));
	printf("maxd R-G & B-G are:	%f	%f \n", dispR_dist, dispB_dist);
}

template <typename T>
void keypnts2file(const char* fnameXYdist,
	vector<T> &xR, vector<T> &yR, vector<T> &xGr, vector<T> &yGr,	// red, green centers
	vector<T> &xB, vector<T> &yB, vector<T> &xGb, vector<T> &yGb)	// blue, (redundant, identical) green centers
{
	printf("Saving corrected keypoints to file... ");
	FILE *pfile_dist;
	pfile_dist = fopen(fnameXYdist, "wt");
	if(pfile_dist == NULL) {
		printf("cannot open file %s.\n", fnameXYdist);
		exit(1);
	}

	int lenR = xR.size();
	int lenB = xB.size();
	int lenmax = std::max(lenR, lenB);
	int lenmin = std::min(lenR, lenB);
	char *ffmt = "%f %f %f %f %f %f %f %f\n";

	for (int i = 0; i < lenmax; i++)
	{
		if (i < lenmin)
		{
			fprintf(pfile_dist, ffmt, xR[i], yR[i], xGr[i], yGr[i], xGb[i], yGb[i], xB[i], yB[i]);
		}
		else
		{
			if (lenR < lenB)
				fprintf(pfile_dist, ffmt, 0.f, 0.f, 0.f, 0.f, xGb[i], yGb[i], xB[i], yB[i]);
			else
				fprintf(pfile_dist, ffmt, xR[i], yR[i], xGr[i], yGr[i], 0.f, 0.f, 0.f, 0.f);
		}
	}
	fclose(pfile_dist);
	printf("done.\n");
}

template <typename T>
void get_polynom(vector<T>& xF, vector<T>& yF, vector<T>& xGf, vector<T>& yGf,
	vector<T>& paramsXF, vector<T>& paramsYF, int degX, int degY, T xp, T yp)
{
	printf("Obtaining correction polynomials... ");
	//vector<T> polyB = getParamsCorrection(xB, yB, xGr, yGr, degX, degY, xp, yp);
	vector<T> polyF = getParamsCorrection(xF, yF, xGf, yGf, degX, degY, xp, yp);
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	paramsXF = polyF.copyRef(0, sizex-1);
	paramsYF = polyF.copyRef(sizex, sizex+sizey-1);
	//paramsXB = polyB.copyRef(0, sizex-1);
	//paramsYB = polyB.copyRef(sizex, sizex+sizey-1);
	printf("done.\n");
}

template <typename T>
void correct_channel(image_char &imgF, image_double &imgFz, vector<T> &paramsXF, vector<T> &paramsYF,
	int spline_order, int degX, int degY, T xp, T yp, int Gcols, int Grows, T scale)
{
	printf("calculating channel correction... ");
	image_double imgD {};
	prepare_spline(imgF, imgD, spline_order);		// spline.h
	for (int i = 0; i < Gcols; i++) {
		for (int j = 0; j < Grows; j++) {
			T p1=0, p2=0;
			undistortPixel(p1, p2, paramsXF, paramsYF, i, j, xp, yp, degX, degY);
			// +0.5 to compensate -0.5 in interpolation function
			T clr = interpolate_image_double(imgD, spline_order, p1/scale+0.5, p2/scale+0.5);
			if (clr < 0) clr = 0; 
			else if (clr > 255) clr = 255;
			imgFz->data[i+j*imgFz->xsize] = clr;
		}
		double percent = ((double)i / (double)Gcols)*100;
		if (!(i % (int)(0.2*Gcols))) printf("%i%c", (int)percent+1, '%');
		else if (!(i % (int)(0.04*Gcols))) printf(".");
	}
	free_image_double(imgD);
	printf("done.\n");
}

template <typename T>
void circuit(int argc, char ** argv, bool clr, bool test = false)
{
	char* fnameRGB = argv[1];
	char* fnameR = argv[2];
	char* fnameG = argv[3];
	char* fnameB = argv[4];
	char* fnameXYdist = argv[5];
	char* fnameXYcorr = argv[6];
	int nImgs = (argc-7)/4;

	T scale = 2;
	image_char img_bayer{};  read_pgm_image_char(img_bayer, fnameRGB);
	// // Uncomment below to replace the line above, make sure you know which direction to rotate (left or right)
	//image_char img_bayer2;  read_pgm_image_char(img_bayer2, fnameRGB);
	//image_char img_bayer; image_rotate_right<T>(img_bayer, img_bayer2); free_image_char(img_bayer2);
	int wi = img_bayer->xsize, he = img_bayer->ysize;
	int wiRB = wi/2, heRB = he/2;
	int Gcols = (int)(wiRB*scale), Grows = (int)(heRB*scale);
	image_char imgR; new_image_char_ini(imgR, wiRB, heRB, 255);
	image_char imgG; new_image_char_ini(imgG, Gcols, Grows, 255);
	image_char imgB; new_image_char_ini(imgB, wiRB, heRB, 255);
	deBayer_char(img_bayer, imgR, imgG, imgB);
	vector<T> xR, yR, xGr, yGr, xB, yB, xGb, yGb, rR, rG, rB;
	keypnts_circle<T>(imgR, imgG, imgB, xR, yR, rR, xGr, yGr, rG, xB, yB, rB, xGb, yGb, scale, clr);
	//keypnts_sift<T>(imgR, imgG, imgB, xR, yR, xGr, yGr, xB, yB, xGb, yGb, scale, clr);
	keypnts2file(fnameXYdist, xR, yR, xGr, yGr, xB, yB, xGb, yGb);
	print_RMSE(xR, yR, xGr, yGr, xB, yB, xGb, yGb);

	vector<T> paramsXR, paramsYR, paramsXB, paramsYB;
	int degX = 11, degY = 11;
	T xp = (T)imgG->xsize/2+0.2, yp = (T)imgG->ysize/2+0.2;
	get_polynom<T>(xR, yR, xGr, yGr, paramsXR, paramsYR, degX, degY, xp, yp);
	get_polynom<T>(xB, yB, xGb, yGb, paramsXB, paramsYB, degX, degY, xp, yp);

	int spline_order = 3;
	image_double imgRz; new_image_double_ini(imgRz, Gcols, Grows, 255);
	image_double imgBz; new_image_double_ini(imgBz, Gcols, Grows, 255);
	printf("Red ");
	correct_channel<T>(imgR, imgRz, paramsXR, paramsYR, spline_order, degX, degY, xp, yp, Gcols, Grows, scale);
	printf("Blue ");
	correct_channel<T>(imgB, imgBz, paramsXB, paramsYB, spline_order, degX, degY, xp, yp, Gcols, Grows, scale);

	printf("\nSaving images to file... \n");
	write_pgm_image_double(imgRz, fnameR);
	write_pgm_image_char(imgG, fnameG);
	write_pgm_image_double(imgBz, fnameB);

	bool green_proc = false;
	vector<T> rB_scale = rB*scale;
	vector<T> rR_scale = rR*scale;
	circle_redefine<T>(imgRz, imgG, imgBz, xR, yR, rR_scale,
                       xGr, yGr, rG, xB, yB, rB_scale, xGb, yGb, 1, clr, xR.size(), green_proc);
	//keypnts_sift<T>(imgRz, imgG, imgBz, xR, yR, xGr, yGr, xB, yB, xGb, yGb, 1, clr);
	keypnts2file(fnameXYcorr, xR, yR, xGr, yGr, xB, yB, xGb, yGb);
	print_RMSE(xR, yR, xGr, yGr, xB, yB, xGb, yGb);

	printf("\nCorrecting blue and red channels for other input images... \n");
	for (int i = 0; i < nImgs; i++)
	{
		image_char imgn_bayer{}; read_pgm_image_char(imgn_bayer, argv[7 + i * 4 + 0]);
		image_char Rin; new_image_char_ini(Rin, wiRB, heRB, 255);
		image_char Gin; new_image_char_ini(Gin, Gcols, Grows, 255);
		image_char Bin; new_image_char_ini(Bin, wiRB, heRB, 255);
		//separate the channels
		deBayer_char(imgn_bayer, Rin, Gin, Bin);
		// measure test image RMSE if necessary
		vector<T> xnR, ynR, xnGr, ynGr, xnB, ynB, xnGb, ynGb, rnR, rnG, rnB;
		if (test) {
			//keypnts_sift<T>(Rin, Gin, Bin, xnR, ynR, xnGr, ynGr, xnB, ynB, xnGb, ynGb, scale, clr);
			keypnts_circle<T>(Rin, Gin, Bin, xnR, ynR, rnR, xnGr, ynGr, rnG, xnB, ynB, rnB, xnGb, ynGb, scale, clr);
			print_RMSE(xnR, ynR, xnGr, ynGr, xnB, ynB, xnGb, ynGb);
		}
		// perform the correction
		image_double Rout; new_image_double_ini(Rout, Gcols, Grows, 255);
		image_double Bout; new_image_double_ini(Bout, Gcols, Grows, 255);
		printf("Red ");
		correct_channel<T>(Rin, Rout, paramsXR, paramsYR, spline_order, degX, degY, xp, yp, Gcols, Grows, scale);
		printf("Blue ");
		correct_channel<T>(Bin, Bout, paramsXB, paramsYB, spline_order, degX, degY, xp, yp, Gcols, Grows, scale);
		// save corrected images to files
		printf("\nSaving images to file... \n");
		write_pgm_image_double(Rout, argv[7+i*4+1]);
		write_pgm_image_char(Gin, argv[7+i*4+2]);
		write_pgm_image_double(Bout, argv[7+i*4+3]);
		// if its test image, measure RMSE
		if (test) {
			//keypnts_sift<T>(Rout, Gin, Bout, xnR, ynR, xnGr, ynGr, xnB, ynB, xnGb, ynGb, 1, clr);
			bool green_proc = false;
			vector<T> rnR_scale = rnR*scale;
			vector<T> rnB_scale = rnB*scale;
			circle_redefine<T>(Rout, Gin, Bout, xnR, ynR,
                               rnR_scale, xnGr, ynGr, rnG, xnB, ynB, rnB_scale, xnGb, ynGb, 1, clr, xnR.size(), green_proc);
			print_RMSE(xnR, ynR, xnGr, ynGr, xnB, ynB, xnGb, ynGb);
		}
		// free memory
		free_image_char(Rin); free_image_char(Gin); free_image_char(Bin);
		free_image_double(Rout); free_image_double(Bout);
	}
	// free memory
	free_image_char(img_bayer);
	free_image_char(imgR); free_image_char(imgG); free_image_char(imgB);
	free_image_double(imgRz); free_image_double(imgBz);
}

template <typename T>
void printMono(FILE* pfile, T mono, int degX, int degY) {
	if (degX != 0 && degY != 0) {
		if (mono > 0)	fprintf(pfile, "+ %.16g * x^%i * y^%i\n", mono, degX, degY);
		else			fprintf(pfile, "- %.16g * x^%i * y^%i\n", -1*mono, degX, degY); }
	else if (degX == 0 && degY == 0) {
		if (mono > 0)	fprintf(pfile, "+ %.16g\n", mono);
		else			fprintf(pfile, "- %.16g\n", -1*mono); }
	else if (degX == 0) {
		if (mono > 0)	fprintf(pfile, "+ %.16g * y^%i\n", mono, degY);
		else			fprintf(pfile, "- %.16g * y^%i\n", -1*mono, degY); }
	else {
		if (mono > 0)	fprintf(pfile, "+ %.16g * x^%i\n", mono, degX);
		else			fprintf(pfile, "- %.16g * x^%i\n", -1*mono, degX); }
}

template <typename T>
void save_poly(char* fname, vector<T>& paramsX, vector<T>& paramsY, const int degX, const int degY)
{
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	//vector<T> paramsX = poly_params.copyRef(0, sizex-1);
	//vector<T> paramsY = poly_params.copyRef(sizex, sizex+sizey-1);
	FILE *pfile;
	pfile = fopen(fname, "wt");
	if(pfile == NULL)
	{
		printf("cannot open file %s.\n", fname);
		return;
	}
	int idx = 0;
	fprintf(pfile, "# polyX(x,y): \n");
	for (int i = degX; i >= 0; i--) {
		for (int j = 0; j <= i; j++) {
			printMono(pfile, paramsX[idx], i-j, j);
			idx++; 	} }
	idx = 0;
	fprintf(pfile, "# polyY(x,y): \n");
	for (int i = degY; i >= 0; i--) {
		for (int j = 0; j <= i; j++) {
			printMono(pfile, paramsY[idx], i-j, j);
			idx++; } }
	fclose(pfile);
}

template <typename T>
void polyEstimation(int argc, char ** argv, bool clr) {
	printf("Polynomial estimation... \n");
	char* fnameRGB = argv[1];
	char* fnamePolyR = argv[2]; 
	char* fnamePolyB = argv[3]; 
	image_char imgR{}, imgG{}, imgB{};

	read_pnm_char(imgR, imgG, imgB, fnameRGB);

	if (5 == argc)
	{
		printf(" done;  writing %s", argv[4]);
		write_ppm_image_char(imgR, imgG, imgB, argv[4]);
	}

	vector<T> xR, yR, xGr, yGr, xB, yB, xGb, yGb, rR, rG, rB;
	keypnts_circle<T>(imgR, imgG, imgB, xR, yR, rR, xGr, yGr, rG, xB, yB, rB, xGb, yGb, 2, clr);
	// solved in gnuplot2file():
	// dxR = cc[0][0] + cc[0][1]*x + cc[0][2]*y + cc[0][3]*x*x + cc[0][4]*y*y + cc[0][5]*x*x*x
	//	 + cc[0][6]*y*y*y + cc[0][7]*x*y + cc[0][8]*x*x*y + cc[0][9]*x*y*y;
	double cc[4][10] = { 0 }, *coef[4] = { cc[0], cc[1], cc[2], cc[3]};
	gnuplot2file("Before_redefine", xR, yR, xGr, yGr, xB, yB, imgR, imgG, coef);
//	gnuplot2file("After_redefine", xR, yR, xGr, yGr, xB, yB, imgR, imgG, coef);
//	keypnts2file(FOLDER "keypnts.p", xR, yR, xGr, yGr, xB, yB, xGb, yGb);
//	print_RMSE(xR, yR, xGr, yGr, xB, yB, xGb, yGb);
	exit(0);

	vector<T> paramsXR, paramsYR, paramsXB, paramsYB;
//	int degX = 5, degY = 5;
	int degX = 11, degY = 11;
	T xp = (T)imgG->xsize/2+0.2, yp = (T)imgG->ysize/2+0.2;
	get_polynom<T>(xR, yR, xGr, yGr, paramsXR, paramsYR, degX, degY, xp, yp);
	get_polynom<T>(xB, yB, xGb, yGb, paramsXB, paramsYB, degX, degY, xp, yp);

	save_poly(fnamePolyR, paramsXR, paramsYR, degX, degY);
	save_poly(fnamePolyB, paramsXB, paramsYB, degX, degY);

	free_image_char(imgR); free_image_char(imgG); free_image_char(imgB);
}

/* Pop out one character from char array */
static char popchar( char *c, int idx, int size) {
	char res = c[idx];
	for (int i = idx; i < size; i++)
		if (i > idx) c[i-1] = c[i];
	return res;
}

static int ok(int ret)
{
	return (0 == ret || EOF == ret) ? 1 : 0;
}

/* Difines the degee values for X and Y polinomials from file */
static void read_degree(FILE *pfile, int& degX, int& degY) {
	degX = 0; degY = 0;
	char buffer[250], sign[2] = { 0 };
	double coef = 0;
	const char* grid = "#", *star = "*" , *plus = "+", *minu = "-";
	bool flagX = false, flagY = false;
	while (!feof(pfile)) {
		/* expects each monimial to have a form of: "+/- coef * x^deg1 * y^deg2 " */
		/* deg1 and/or deg2 can be zeros, leaving just a coefficent value */
		if (!ok(fscanf(pfile, "%s", sign)))
		{
			printf("read_degree():  bad sign\n");
			break;
		}
		/* if sign "#" is met - it's a comment, we may skip it. */
		if (strcmp(sign, grid) != 0) {
			char mono1[5], mono2[5];
			if (!ok(fscanf(pfile, "%lf", &coef)))
				break;
			long currPos = ftell (pfile);
			if (!ok(fscanf(pfile, "%s", sign)))
				break;
			if (strcmp(sign, star) != 0) {
				if (!feof(pfile)) {
					fseek(pfile, currPos, SEEK_SET);
					assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} }
			else {
				if (!ok(fscanf(pfile, "%s", mono1)))
					break;
				currPos = ftell (pfile);
				if (!ok(fscanf(pfile, "%s", sign)))
					break;
				if (strcmp(sign, star) != 0) {
					if (!feof(pfile)) {
						fseek(pfile, currPos, SEEK_SET);
						assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0);
					}
				}
				else if (!ok(fscanf(pfile, "%s", mono2)))
					break;
			}
			/* pop-out "x^" and "y^" so that to calculate degrees. */
			popchar(mono1, 0, 5); popchar(mono1, 0, 5);
			popchar(mono2, 0, 5); popchar(mono2, 0, 5);
			int tmpdeg = std::max(atoi(mono1), atoi(mono2));
			/* there are two polynomials in file: degree of X is defined by flagX, second - flagY. */
			/* assumed first poly is for X direction, second is for Y. */
			if (flagX)
			{
				if (tmpdeg > degX)
					degX = tmpdeg; }
				else if (tmpdeg > degY)
					degY = tmpdeg;		
		} else {
			if (!ok(fscanf(pfile, "%s", buffer)))
				break;
			if (degX == 0 && degY == 0)
				flagX = true;
			else { flagX = false; flagY = true; }
		}
	}
	fseek(pfile, 0, SEEK_SET);
}

/* Finds an index of coefTerm[] vector based on given x and y degrees. */
static int coefIdx(int degree, int x, int y) {
	int a1 = degree + 1;
	int n = std::abs(x+y-degree)+1;
	int an = x+y+1;
	int Sn = n * (a1 + an);
	Sn /= 2;
	return Sn-an+y; }

/* Reads the poly coefficients and insert them into vector of coefficients - coefTerm[]. */
/* Also returns the degrees for each polynomial. */
template <typename T>
vector<T> read_poly(char* fname, int& degX, int& degY) {
	FILE *pfile;
	pfile = fopen(fname, "r");
	if(pfile == NULL) printf("unable to open file %s.\n", fname);
	/* get the degrees for each polynomial - need it to know the size of vectors paramsX and paramsY. */
	read_degree(pfile, degX, degY);
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	vector<T> paramsX(sizex), paramsY(sizey);
	char buffer[500], sign[2] = { '\0' };
	const char* grid = "#", *star = "*", *plus = "+", *minu = "-";
	bool flagX = false, flagY = false;
	while (!feof(pfile)) {
		/* reading is done by the same manner as in "read_degree(pfile, degX, degY);" */
		if (!ok(fscanf(pfile, "%s", sign)))
			break;
		if (strncmp(sign, grid, 2) != 0 && (strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 ) ) {
			char mono1[5], mono2[5];
			double coef = 0;
			if (!ok(fscanf(pfile, "%lf", &coef)))
				break;
			long currPos = ftell (pfile);
			/* save the coef sign. */
			if (strcmp(sign, minu) == 0) coef *= -1;
			if (!ok(fscanf(pfile, "%s", sign)))
				break;
			if (strcmp(sign, star) != 0) {
				if (!feof(pfile)) {
					fseek(pfile, currPos, SEEK_SET);
					assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} 	}
			else {
				if (!ok(fscanf(pfile, "%s", mono1)))
					break;
				currPos = ftell (pfile);
				if (!ok(fscanf(pfile, "%s", sign)))
					break;
				if (strcmp(sign, star) != 0) {
					if (!feof(pfile)) {
						fseek(pfile, currPos, SEEK_SET);
						assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0);
					}
				}
				else if (!ok(fscanf(pfile, "%s", mono2)))
					break;
			}
			char x = popchar(mono1, 0, 5); popchar(mono1, 0, 5);
			char y = popchar(mono2, 0, 5); popchar(mono2, 0, 5);
			/* see which degree belongs to which variable. */
			int tmpdegX = 0, tmpdegY = 0;
			if ( x == 'x')		tmpdegX = atoi(mono1);
			else if ( x == 'y')	tmpdegY = atoi(mono1);
			if ( y == 'y') 	tmpdegY = atoi(mono2);
			/* save the coef to accoring paramsX/Y; fist poly in file belongs to paramsX, second - paramsY. */
			if (flagX && !flagY) {
				int idx = coefIdx(degX, tmpdegX, tmpdegY);
				paramsX[idx] = coef; }
			else {
				int idx = coefIdx(degY, tmpdegX, tmpdegY);
				paramsY[idx] = coef; }
		}
		else {
			if (!ok(fscanf(pfile, "%s", buffer)))
				break;
			if (!flagX) flagX = true;
			else flagY = true; }
	}
	fclose(pfile);
	vector<T> poly_params(sizex+sizey);
	/* copy the paramsX and paramsY to one vector. */
	/* later we will be able to separate them since the degrees for each poly are known. */
	for (int k = 0; k < sizex; k++) poly_params[k] = paramsX[k];
	for (int k = 0; k < sizey; k++) poly_params[k+sizex] = paramsY[k];
	return poly_params;
}

template <typename T>
void aberCorrection(int argc, char ** argv, bool clr)
{
	printf("\nAberration correction... \n");
	char* fnameRGB = argv[1];
	char* fnamePolyR = argv[2];
	char* fnamePolyB = argv[3];
	char *fnameR = argv[4], *fnameG, *fnameB;
	if (7 == argc)
	{
		fnameG = argv[5];
		fnameB = argv[6];
	}
	int degX = 11, degY = 11;
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	image_char Rin{}, Gin{}, Bin{};

	read_pnm_char(Rin, Gin, Bin, fnameRGB);
	uint Gcols = Gin ? Gin->xsize : 0, Grows = Gin ? Gin->ysize : 0;

	vector<T> paramsR = read_poly<T>(fnamePolyR, degX, degY);
	vector<T> paramsB = read_poly<T>(fnamePolyB, degX, degY);
	vector<T> paramsXR = paramsR.copyRef(0, sizex-1);
	vector<T> paramsYR = paramsR.copyRef(sizex, sizex+sizey-1);
	vector<T> paramsXB = paramsB.copyRef(0, sizex-1);
	vector<T> paramsYB = paramsB.copyRef(sizex, sizex+sizey-1);

//	printf("Gcols = %d;  Grows = %d, Gin->xsize = %d, Gin->ysize = %d for %s\n",
//			Gcols, Grows, Gin->xsize, Gin->ysize, fnameRGB);
//	Gcols = 5634;  Grows = 3752, Gin->xsize = 5634, Gin->ysize = 3752 for ../../../../data/_MG_7626.pgm

	int spline_order = 3;
	image_double Rout; new_image_double_ini(Rout, Gcols, Grows, 255);
	image_double Bout; new_image_double_ini(Bout, Gcols, Grows, 255);
	T xp = 0.2, yp = 0.2;
	if (Gin) {
		xp += Gin->xsize/2;
		yp += Gin->ysize/2;
 	}	else error("aberCorrection:  Null Gin");

	printf("Red  ");
	correct_channel<T>(Rin, Rout, paramsXR, paramsYR, spline_order, degX, degY, xp, yp, Gcols, Grows, 2);
	printf("Blue ");
	correct_channel<T>(Bin, Bout, paramsXB, paramsYB, spline_order, degX, degY, xp, yp, Gcols, Grows, 2);

	printf("\nSaving images to file... \n");
	if (7 == argc)
	{
		write_pgm_image_double(Rout, argv[4]);
		write_pgm_image_char(Gin, argv[5]);
		write_pgm_image_double(Bout, argv[6]);
	}
	else write_ppm_image_double(Rout, Gin, Bout, argv[4]);

	free_image_char(Rin); free_image_char(Gin); free_image_char(Bin);
	free_image_double(Rout); free_image_double(Bout);
}

int main(int argc, char ** argv)
{
	bool clr = false; // deals with black circles on white background
	bool test = false; // true if the image to correct is a test image to measure the correction RMSE

	if (1 == argc)
	{
		const char * foo[] = { argv[0], "../../../../data/_MG_7626.pgm",
								"../../../../data/_MG_7626_polyR.txt", "../../../../data/_MG_7626_polyB.txt",
								FOLDER "_MG_7626R.pgm", FOLDER "_MG_7626G.pgm", FOLDER "_MG_7626B.pgm" };
		return test_read_matrix(foo[1]);
//		foo[1] = FOLDER "uncorrected.ppm";
		foo[2] = FOLDER "BayerIMG_7626_polyR.txt";
		foo[3] = FOLDER "BayerIMG_7626_polyB.txt";
//		foo[4] = FOLDER "BayerFromPPM_7626.ppm";
//		foo[4] = FOLDER "BayerFromPGM_7626.ppm";
		foo[4] = FOLDER "BayerFromPGM_7626R.pgm";
		foo[5] = FOLDER "BayerFromPGM_7626G.pgm";
		foo[6] = FOLDER "BayerFromPGM_7626B.pgm";
//		foo[5] = foo[6] = "";

		printf("Polynomial estimation:\n");
//		printf("%s %s %s %s\n", foo[0], foo[1], foo[2], foo[3]);
//		printf("%s %s %s %s %s\n", foo[0], foo[1], foo[2], foo[3], foo[4]);
		printf("%s %s %s %s %s %s %s\n", foo[0], foo[1], foo[2], foo[3], foo[4], foo[5], foo[6]);
		polyEstimation<double>(7, (char**)foo, clr);
		return 0;

		printf("CA Polynomial correction:\n");
		printf("%s %s %s %s %s %s %s\n", foo[0], foo[1], foo[2], foo[3], foo[4], foo[5], foo[6]);
		aberCorrection<double>(5, (char**)foo, clr);
		return 0;
	}

	else if (argc > 7)	// runs all circuit, change settings inside
		circuit<double>(argc, argv, clr, test);

	else if (argc == 4 || 5 == argc)	// estimates and saves polynomial; optionally writes uncorrected PPM
		polyEstimation<double>(argc, argv, clr);

	else if (argc == 7)	// reads image and poly, corrects input and saves corrected channels separately
		aberCorrection<double>(argc, argv, clr);

	else {
		printf("Program usage: \n");
		printf("Polynomial estimation:\n");
		printf("CAcorr fname_raw_calib.pgm fname_poly_red.txt fname_poly_blue.txt \n\n");
		printf("CA correction using estimated polynomial:\n");
		printf("CAcorr fname_raw.pgm fname_poly_red.txt fname_poly_blue.txt "
				"fname_raw_red_corr.pgm fname_raw_green_corr.pgm fname_raw_blue_corr.pgm \n\n");
		printf("Running all circuit (polynomial estimation - image correction):\n");
		printf("CAcorr fname_raw_calib.pgm fname_raw_calib_red_corr.pgm fname_raw_calib_green_corr.pgm"
		" fname_raw_calib_blue_corr.pgm fname_raw_calib_keyp_dist.txt fname_raw_calib_keyp_corr.txt "
		"[fname_img_n.pgm fname_img_n_red_corr.pgm fname_img_n_green_corr.pgm fname_img_n_blue_corr.pgm, ...]\n\n");
	}

	return 0;
}
