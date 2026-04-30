## [Cimg open-source C++ library for image processing](https://www.cimg.eu/)
Cimg can [delegate file I/O](https://stackoverflow.com/a/47375364) to [ImageMagick](https://imagemagick.org/download/) `convert` or `magick`
- tick the box entitled `Install legacy commands` or similar, during ImageMagick installation.
### Bicubic resampling - not bundled with Cimg
- [OscarCasadoLorenzo / BicubicInterpolation](https://github.com/OscarCasadoLorenzo/BicubicInterpolation)
	- lacks explicit Windows support

### instead of Cimg, convert non-PPM images to PPM using ImageMagick's `magick.exe`
- [execute another .exe, in Windows from C++](https://stackoverflow.com/questions/17703721/how-to-execute-another-exe-from-a-c-program-in-windows)
	- [`CreateProcessA()` function](https://learn.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-createprocessa)
- or just run a .bat file, using magick to convert to/from PPM
	- [ImageMagick Usage under Windows](https://usage.imagemagick.org/windows/)

### [Two dimensional bicubic interpolation implementation in C++](https://codereview.stackexchange.com/questions/263289/two-dimensional-bicubic-interpolation-implementation-in-c)
- [Bi-Cubic resize implementation in C++](https://stackoverflow.com/questions/17640173/implementation-of-bi-cubic-resize)

## CorrCA `interpolate_spline()` and `prepare_spline()`
core spline interpolation functions used for CorrCA chromatic aberration correction.  

## prepare_spline()

`prepare_spline()` pre-processes an image to enable high-quality spline interpolation 
by computing spline coefficients through recursive filtering [2](#0-1) .

**Key aspects:**
- Only processes images with spline order >= 3 (returns true for lower orders)  
- Initializes poles of the associated z-filter based on spline order   
- Applies 1D inverse spline filtering horizontally (lines) and vertically (columns)   

## interpolate_spline()

`interpolate_spline()` samples an image at sub-pixel coordinates  
using various interpolation methods.  

**Supported interpolation orders:**
- Order 0: Nearest neighbor (pixel replication)  
- Order 1: Bilinear interpolation  
- Order -3: Bicubic Keys' function (Catmull-Rom)  
- Order 3: Cubic spline  
- Orders 5,7,9,11: Higher-order splines  

**Usage in correction pipeline:**
The functions work together in `correct_channel()`
- `prepare_spline()` is called once per channel before warping, then
- `interpolate_spline()` is called for each pixel during the correction loop.  

### Notes
- `prepare_spline()` modifies the input image in-place to store spline coefficients
- `interpolate_spline()` handles boundary conditions by returning 0 for out-of-bounds coordinates
- The default Keys parameter (-0.5) in `interpolate_spline()` corresponds to Catmull-Rom spline
- These functions are essential for sub-pixel accurate image warping chromatic aberration correction

<details><summary><b>Citations</b></summary>

**File:** src/spline.h (L51-55)
```text
bool prepare_spline(image_double& im, int order);
bool interpolate_spline(image_double& im, int order,
                        double x, double y,
                        double& out,
                        double paramKeys=-.5f);
```

**File:** src/spline.cpp (L180-196)
```cpp
bool prepare_spline(image_double& im, int order)
{
    if(order < 3)
        return true;

    // Init poles of associated z-filter
    double z[5];
    if(! fill_poles(z, order))
        return false;
    int npoles = order/2;

	for(int y = 0; y < im->ysize; y++) // Filter on lines
		invspline1D(im->data+y*im->xsize, 1, im->xsize, z, npoles);
	for(int x = 0; x < im->xsize; x++) // Filter on columns
		invspline1D(im->data+x, 1*im->xsize, im->ysize, z, npoles);
    return true;
}
```

**File:** src/spline.cpp (L312-386)
```cpp
bool interpolate_spline( image_double& im, int order, double x, double y, double& out, double paramKeys)
{
	double  cx[12],cy[12];

	/* CHECK ORDER */
	if(order != 0 && order != 1 && order != -3 &&
		order != 3 && order != 5 && order != 7 && order != 9 && order != 11)
		return false;

	double ak[13];
	if(order > 3)
		init_splinen(ak, order);

	bool bInside = false;
	/* INTERPOLATION */
	if(order == 0) { /* zero order interpolation (pixel replication) */
		int xi = (int)floor((double)x);
		int yi = (int)floor((double)y);
		bInside = valid_image_double(im, xi, yi);//im.valid(xi, yi);
		if(bInside) {
			double p = im->data[xi+yi*im->xsize]; //im.pixel(xi, yi);
			out = p;
		}
	} else { /* higher order interpolations */
		bInside = (x>=0.0f && x<=(double)im->xsize && y>=0.0f && y<=(double)im->ysize);
		if(bInside) {
			x -= 0.5f; y -= 0.5f;
			int xi = (x<0)? -1: (int)x;
			int yi = (y<0)? -1: (int)y;
			double ux = x - (double)xi;
			double uy = y - (double)yi;
			switch(order)  {
			case 1: /* first order interpolation (bilinear) */
				cx[0] = ux; cx[1] = 1.0f-ux;
				cy[0] = uy; cy[1] = 1.0f-uy;
				break;
			case -3: /* third order interpolation (bicubic Keys' function) */
				keys(cx, ux, paramKeys);
				keys(cy, uy, paramKeys);
				break;
			case 3: /* spline of order 3 */
				spline3(cx, ux);
				spline3(cy, uy);
				break;
			default: /* spline of order >3 */
				splinen(cx, ux, ak, order);
				splinen(cy, uy, ak, order);
				break;
			}
			int n2 = (order==-3)? 2: (order+1)/2;
			int n1 = 1-n2;
			/* this test saves computation time */
			if(valid_image_double(im, xi+n1, yi+n1) && valid_image_double(im, xi+n2,yi+n2)) {
				out = 0.0f;
				for(int dy = n1; dy <= n2; dy++) {
					int adrs = (xi+n1) + (yi+dy) * im->xsize;
					for(int dx = n1; dx <= n2; dx++) {
						double f = im->data[adrs];
						out += cy[n2-dy]*cx[n2-dx] * f;
						adrs++;
					}
				}

			} else
				out = 0.0f;
			for(int dy = n1; dy <= n2; dy++)
				for(int dx = n1; dx <= n2; dx++) {
					double v = 0.0f; // the image is not infinite, there is no data outside
					out += cy[n2-dy]*cx[n2-dx]*v;
				}

		}
	}
	return bInside;
}
```

**File:** src/main_centering.cpp (L575-580)
```cpp
	prepare_spline(imgF, spline_order);
	for (int i = 0; i < wiG; i++) {
		for (int j = 0; j < heG; j++) {
			T p1=0, p2=0;
			undistortPixel(p1, p2, paramsXF, paramsYF, i, j, xp, yp, degX, degY);
			T clr = interpolate_image_double(imgF, spline_order, p1/scale+0.5, p2/scale+0.5); // +0.5 to compensate -0.5 in interpolation function
```

</details>
