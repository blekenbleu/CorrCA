## Chromatic aberration (CA) correction software (C/C++)
from [corrCA-prototype](https://github.com/vicrucann/corrCA-prototype)

### See branch [PPM](https://github.com/blekenbleu/CorrCA/tree/PPM)
Branch `main` expects `fname_raw_calib.pgm` to have Bayer-matrixed RGB *pixels*,  
unlike PPM (AKA [portable pixmap](https://en.wikipedia.org/wiki/Netpbm)), which have RGB components (subpixels) for each pixel.  
The [PPM branch](https://github.com/blekenbleu/CorrCA/tree/PPM) adds `read_ppm_image_double()` to handle `.ppm` files.
- `chromaberrat` renders a green color plane to the same size as Bayer-matrixed input,  
	interpolating to replace red and blue pixels
- `chromaberrat` renders red and blue color planes to *half* Bayer-matrixed input width and height.
- `write_ppm_image_double()` outputs PPM files sized to red and blue color planes, subsampling green pixels.
- for compatibility, `read_ppm_image_double()` interpolates PPM green plane  
	back to twice width and height of red and blue planes.

### Input parameters
- [Debugging with command-line parameters in Visual Studio](https://stackoverflow.com/questions/298708/debugging-with-command-line-parameters-in-visual-studio)

Depending on argument count, `chromaberrat` performs different tasks:

##### Three input arguments 
* `fname_raw_calib.pgm fname_poly_red.txt fname_poly_blue.txt`  
* Based on raw calibration image, the program estimates correction polynomials  
	for red and bleu channels and saves them to chosen txt files  
* **EXAMPLE**: `data/_MG_7626.pgm data/_MG_7626_polyR.txt data/_MG_7626_polyB.txt`  

##### Six input arguments  
* `fname_raw.pgm fname_poly_red.txt fname_poly_blue.txt fname_raw_red_corr.pgm fname_raw_green_corr.pgm fname_raw_blue_corr.pgm`  
* Read a raw image that is needed to be corrected, reads correction polynomials  
	and performs the correction of the image; three corrected channels are saved separately  
* **EXAMPLE**: `data/_MG_7628.pgm data/_MG_7626_polyR.txt data/_MG_7626_polyB.txt data/_MG_7628_R_corr.pgm data/_MG_7628_G_corr.pgm data/_MG_7628_B_corr.pgm`  

##### More than six input arguments  
* `fname_raw_calib.pgm fname_raw_calib_red_corr.pgm fname_raw_calib_green_corr.pgm fname_raw_calib_blue_corr.pgm
 fname_raw_calib_keyp_dist.txt fname_raw_calib_keyp_corr.txt
 [fname_img_n.pgm fname_img_n_red_corr.pgm fname_img_n_green_corr.pgm fname_img_n_blue_corr.pgm, ...]`  
- Runs all the circuit: first it builds the correction polynomials based on calibration raw image,  
	it corrects the calibration image and all other image arguments with those polynomials.  
	In txt files, it saves the geometrical coordinates of the calibration image keypoints in all the channels  
	- before and after the correction.  
	Those txt files can be used later for visualization tests in Matlab (i.e., misalignment field and histograms).
	Correction polynomials are not saved; for those, use  **Three input arguments**.  
* **EXAMPLE**: `data/_MG_7626.pgm data/_MG_7626_R_corr.pgm data/_MG_7626_G_corr.pgm data/_MG_7626_B_corr.pgm
	 data/_MG_7626_dist.txt data/_MG_7626_corr.txt data/_MG_7628.pgm data/_MG_7628_R_corr.pgm data/_MG_7628_G_corr.pgm data/_MG_7628_B_corr.pgm`  

#### *29 Apr 2026*
- build with Visual Studio 2022
  - by Open -> Folder
- process `data/_MG_7626.pgm` by default (simplify Visual Studio debugging)
  - 'LMA center redefinition' runs slowly...
  - calculates x,y polynomial coefficients to `x^11`, `y^11`;
	- should check coefficients and stop when exceedingly small...

#### *30 Apr 2026*
- employ [DeepWiki](https://docs.devin.ai/work-with-devin/deepwiki)

#### *2 May 2026* PPM branch
- inplement and visually debug `read_pgm_image_double()`, `write_pgm_image_double()`

#### *3 May 2026* PPM branch
- binary `write_ppm_image_double()`, `write_pgm_image_double()`  

## Supported Image Format (from [DeepWiki](https://deepwiki.com/blekenbleu/CorrCA))
Portable Pixmap (PPM) image files had not been supported in this codebase,  
which specifically worked with Portable Gray Map (PGM) format files. 

- **Input files**: Raw Bayer-matrix calibration and correction images use `.pgm` extension (e.g., `_MG_7626.pgm`)  
- **Output files**: Corrected R,G,B channel images saved as separate PGM files  
	(e.g., `_R_corr.pgm`, `_G_corr.pgm`, `_B_corr.pgm`)  

### Implementation Details

Image I/O uses:
- `read_pgm_image_double()` for reading PGM files  
- `write_pgm_image_double()` for writing PGM files  
- `read_ppm_image_double()` for reading PPM files  
- `write_ppm_image_double()` for writing PPM files  

`read_pgm_image_double()` reads Bayer-patterned images  
	&emsp; *where color components are interleaved in a single grayscale channel*.

### Notes

The codebase reads raw Bayer pattern images where R, G, and B pixels were interleaved in PGM files.  
These are separated into R,G,B planes for processing, with PGM format for all image file I/O operations.  

DeepWiki pages to explore:
- [Project Layout and Data Files](https://deepwiki.com/blekenbleu/CorrCA/1.2-project-layout-and-data-files)
- [Glossary](https://deepwiki.com/blekenbleu/CorrCA/5-glossary)

<details><summary><b>Citations</b></summary>

**File:** [README.md](#three-input-arguments)
```markdown
* `fname_raw_calib.pgm fname_poly_red.txt fname_poly_blue.txt`  
* Based on raw calibration image, the program estimates correction polynomials for red and bleu channels and saves them to chosen txt files  
* **EXAMPLE**: `data/_MG_7626.pgm data/_MG_7626_polyR.txt data/_MG_7626_polyB.txt`  

###### Six input arguments  
* `fname_raw.pgm fname_poly_red.txt fname_poly_blue.txt fname_raw_red_corr.pgm fname_raw_green_corr.pgm fname_raw_blue_corr.pgm`
* Read a raw image that is needed to be corrected, reads correction polynomials and performs the correction of the image; three corrected channels are saved separately
* **EXAMPLE**: `data/_MG_7628.pgm data/_MG_7626_polyR.txt data/_MG_7626_polyB.txt data/_MG_7628_R_corr.pgm data/_MG_7628_G_corr.pgm data/_MG_7628_B_corr.pgm`
```

**File:** [src/main_centering.cpp](https://github.com/blekenbleu/CorrCA/blob/3f73695786a4ddc9241167a208be13184f5618cb/src/main_centering.cpp#L604)
```cpp
	image_double img_bayer = read_pgm_image_double(fnameRGB);
```

**File:** [src/main_centering.cpp](https://github.com/blekenbleu/CorrCA/blob/3f73695786a4ddc9241167a208be13184f5618cb/src/main_centering.cpp#L634)
```cpp
	write_pgm_image_double(imgRz, fnameR);
	write_pgm_image_double(imgG, fnameG);
	write_pgm_image_double(imgBz, fnameB);
```

</details>

## To process white circles on black background
Change the `clr` parameter from `false` to `true`:  

#### 1. Change line 952 in `src/main_centering.cpp`:
```cpp
bool clr = true; // deals with white circles on black background
```

#### 2. Thresholding Logic
The `clr` parameter affects thresholding in two key functions:

**In `initial_tache()` (lines 73-77):**
- When `clr = false` (black on white): `seuil = val_bas + (val_haut - val_bas)/3 * 2;` 
- When `clr = true` (white on black): `seuil = val_bas + (val_haut - val_bas)/3;`

**`binarization()` (lines 281-313)**
checks if pixel values are `<= threshold` and sets them to 0
With `clr = true`, this logic inverts to properly detect white circles.

#### 3. Function Calls
The `clr` parameter is passed through the processing chain:
- `circuit()` → `keypnts_circle()` → `initial_tache()`
- `circle_redefine()` → `centerLMA()` → `initial_tache()`

### Notes
The `clr` parameter controls the entire detection pipeline for "tache" (calibration spots).  
Setting it to `true` inverts the thresholding logic throughout the system  
to handle white circles on black backgrounds instead of the default black circles on white backgrounds.

<details><summary>Citations</summary>

**File:** src/main_centering.cpp (L45-46)
```cpp
int initial_tache(image_double I, vector<T>& h, T& rayon, bool color, T x, T y) {
	int COL_IMA = I->xsize; 
```

**File:** src/main_centering.cpp (L74-75)
```cpp
	if (!color)
		seuil = val_bas + (val_haut - val_bas)/3 * 2;
```

**File:** src/main_centering.cpp (L76-77)
```cpp
	else if (color)
		seuil = val_bas + (val_haut - val_bas)/3;
```

**File:** src/main_centering.cpp (L293-293)
```cpp
			if (red <= threR) imgbiR->data[i+j*wiRB] = 0;
```

**File:** src/main_centering.cpp (L344-344)
```cpp
		centerLMA<T>(sub_imgR, clr, cxR, cyR);
```

**File:** src/main_centering.cpp (L616-616)
```cpp
	keypnts_circle<T>(imgR, imgG, imgB, xR, yR, rR, xGr, yGr, rG, xB, yB, rB, xGb, yGb, scale, clr);
```

**File:** src/main_centering.cpp (L952-952)
```cpp
	bool clr = false; // deals with black circles on white background
```

</details>

## Dynamic Circle Detection and Count Determination in CorrCA
CorrCA dynamically determines circle count  
through connected component analysis rather than using hard-coded values.  
Circles are detected by binarizing images [1a],  
finding connected components [1b],  
filtering by size (>180 pixels) and shape (compactness 0.7-1.3, radius >8) [1c,1d],  
then storing valid circles [1e].  
The resulting count drives the entire processing pipeline [2a,2b,2c,2d]  
and populates coordinate arrays for matching across color channels [3a,3b,3c,3d].

<details>

#### 1. Circle Detection Pipeline - Connected Component Analysis
How the system dynamically discovers and counts calibration circles from image content
#### 1a. Binarize image channels (`main_centering.cpp:399`)
Convert color channels to binary for connected component detection
```text
binarization(imgbiR, imgbiG, imgbiB, imgR, imgG, imgB, threR, threG, threB);
```
#### 1b. Detect connected components in red channel (`main_centering.cpp:404`)
Find all connected components and store their statistics
```text
CC(ccstatsR, imgbiR, 'R'); printf(" number = %i ", ccstatsR.size());
```
#### 1c. Filter by minimum pixel count (`abberation.h:156`)
Only consider components with more than 180 pixels
```text
if (npix > 180) {
```
#### 1d. Filter by shape and size (`abberation.h:159`)
Accept only circular shapes with radius > 8 and compactness 0.7-1.3
```text
if (std::min(stats.radius1, stats.radius2) > 8 && compactness < 1.3 && compactness > 0.7)
```
#### 1e. Add valid circle to results (`abberation.h:161`)
Store statistics for circles that pass all filters
```text
ccstats.push_back(stats);
```
#### 2. Circle Count Allocation and Processing
How the dynamically detected circle count drives the processing pipeline
#### 2a. Report detected circle counts (`main_centering.cpp:408`)
Display the number of circles found in each color channel
```text
printf("\nnumber of connected components per channel R=%i G=%i B=%i\n", ccstatsR.size(), ccstatsG.size(), ccstatsB.size());
```
#### 2b. Validate equal circle counts (`main_centering.cpp:409`)
Ensure all channels detected the same number of circles
```text
assert(ccstatsR.size() == ccstatsG.size() && ccstatsG.size() == ccstatsB.size());
```
#### 2c. Extract circle count (`main_centering.cpp:413`)
Store the dynamic circle count for processing
```text
int ntaches = ccstatsG.size();
```
#### 2d. Allocate arrays based on circle count (`main_centering.cpp:414`)
Create coordinate arrays sized to the detected circle count
```text
xR = xR.ones(ntaches); yR = yR.ones(ntaches); rR = rR.ones(ntaches);
```
#### 3. Circle Matching and Data Structure Population
How the system matches circles across color channels and populates processing arrays
#### 3a. Iterate through all detected circles (`main_centering.cpp:420`)
Process each circle found in the green channel
```text
for (int i = 0; i < ntaches; i++) {
```
#### 3b. Find matching red circle (`main_centering.cpp:423`)
Locate corresponding circle in red channel
```text
int idxR = findMatch(xg, yg, ccstatsR, scale);
```
#### 3c. Find matching blue circle (`main_centering.cpp:424`)
Locate corresponding circle in blue channel
```text
int idxB = findMatch(xg, yg, ccstatsB, scale);
```
#### 3d. Store red circle coordinates (`main_centering.cpp:429`)
Save matched circle center coordinates
```text
xR[i] = ccstatsR[idxR].centerX; yR[i] = ccstatsR[idxR].centerY;
```

</details>

## Reducing polynomial fit precision, e.g. coefficient count
To reduce polynomial fit precision and coefficient count,  
lower `deg[XY]` from the default value of 11.  
`deg[XY]` is hardcoded in three main functions in `src/main_centering.cpp`.

### Key Changes Required

#### 1. Modify `deg[XY]` in `circuit()` function
Change line 622:
```cpp
int degX = 5, degY = 5;  // Reduced from 11
```

#### 2. Modify `deg[XY]` in `polyEstimation()` function  
Change line 756:
```cpp
int degX = 5, degY = 5;  // Reduced from 11
```

#### 3. Modify `deg[XY]` in `aberCorrection()` function
Change line 913:
```cpp
int degX = 5, degY = 5;  // Reduced from 11
```

### Coefficient Count Calculation

The number of coefficients is `(deg + 1) * (deg + 2) / 2`  
for each polynomial dimension:

- **Degree 11**: 78 coefficients per dimension (156 total)
- **Degree 5**: 21 coefficients per dimension (42 total)
- **Degree 3**: 10 coefficients per dimension (20 total)

### Impact on Processing

Lowering `deg[XY]` affects:
1. **Polynomial estimation** in `get_polynom()`
2. **Channel correction** in `correct_channel()`
3. **File I/O** when reading/writing polynomial files

### Notes
Reducing `deg[XY]` decreases correction accuracy and computations.  
Experiment with different values to balance precision and performance for specific use cases.

<details><summary><b>Citations</b></summary>

**File:** src/main_centering.cpp (L555-568)
```cpp
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
```

**File:** src/main_centering.cpp (L571-590)
```cpp
void correct_channel(image_double& imgF, image_double& imgFz, vector<T>& paramsXF, vector<T>& paramsYF, 
	int spline_order, int degX, int degY, T xp, T yp, int wiG, int heG, int scale)
{
	printf("\nCorrected channel is being calculated... \n");
	prepare_spline(imgF, spline_order);
	for (int i = 0; i < wiG; i++) {
		for (int j = 0; j < heG; j++) {
			T p1=0, p2=0;
			undistortPixel(p1, p2, paramsXF, paramsYF, i, j, xp, yp, degX, degY);
			T clr = interpolate_image_double(imgF, spline_order, p1/scale+0.5, p2/scale+0.5); // +0.5 to compensate -0.5 in interpolation function
			if (clr < 0) clr = 0; 
			if (clr > 255) clr = 255;
			imgFz->data[i+j*imgFz->xsize] = clr;
		}
		double percent = ((double)i / (double)wiG)*100;
		if (!(i % (int)(0.2*wiG))) printf("%i%c", (int)percent+1, '%');
		else if (!(i % (int)(0.04*wiG))) printf(".");
	}
	printf("done.\n");
}
```

**File:** src/main_centering.cpp (L837-900)
```cpp
vector<T> read_poly(char* fname, int& degX, int& degY) {
	FILE *pfile;
	pfile = fopen(fname, "r");
	if(pfile == NULL) printf("unable to open file %s.\n", fname);
	/* get the degrees for each polynomial - need it to know the size of vectors paramsX and paramsY. */
	read_degree(pfile, degX, degY);
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	vector<T> paramsX(sizex), paramsY(sizey);
	char buffer[500], sign[2];
	const char* grid = "#", *star = "*", *plus = "+", *minu = "-";
	bool flagX = false, flagY = false;
	while (!feof(pfile)) {
		/* reading is done by the same manner as in "read_degree(pfile, degX, degY);" */
		fscanf(pfile, "%s", sign);
		if (strcmp(sign, grid) != 0 && (strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 ) ) {
			char mono1[5], mono2[5];
			double coef = 0;
			fscanf(pfile, "%lf", &coef);
			long currPos = ftell (pfile);
			/* save the coef sign. */
			if (strcmp(sign, minu) == 0) coef *= -1;
			fscanf(pfile, "%s", sign);
			if (strcmp(sign, star) != 0) {
				if (!feof(pfile)) {
					fseek(pfile, currPos, SEEK_SET);
					assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} 	}
			else {
				fscanf(pfile, "%s", mono1);
				currPos = ftell (pfile);
				fscanf(pfile, "%s", sign);
				if (strcmp(sign, star) != 0) {
					if (!feof(pfile)) {
						fseek(pfile, currPos, SEEK_SET);
						assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} 	}
				else fscanf(pfile, "%s", mono2); }
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
			fscanf(pfile, "%s", buffer);
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
```
</details>
