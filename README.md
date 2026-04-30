## Chromatic aberration (CA) correction software (C/C++)
from [corrCA-prototype](https://github.com/vicrucann/corrCA-prototype)

#### Input parameters

Depending on the number of paramters, the program performs different tasks:

###### Three input arguments 
* `fname_raw_calib.pgm fname_poly_red.txt fname_poly_blue.txt`  
* Based on raw calibration image, the program estimates correction polynomials for red and bleu channels and saves them to chosen txt files  
* **EXAMPLE**: `data/_MG_7626.pgm data/_MG_7626_polyR.txt data/_MG_7626_polyB.txt`  

###### Six input arguments  
* `fname_raw.pgm fname_poly_red.txt fname_poly_blue.txt fname_raw_red_corr.pgm fname_raw_green_corr.pgm fname_raw_blue_corr.pgm`  
* Read a raw image that is needed to be corrected, reads correction polynomials and performs the correction of the image; three corrected channels are saved separately  
* **EXAMPLE**: `data/_MG_7628.pgm data/_MG_7626_polyR.txt data/_MG_7626_polyB.txt data/_MG_7628_R_corr.pgm data/_MG_7628_G_corr.pgm data/_MG_7628_B_corr.pgm`  

###### More than six input arguments  
* `fname_raw_calib.pgm fname_raw_calib_red_corr.pgm fname_raw_calib_green_corr.pgm fname_raw_calib_blue_corr.pgm fname_raw_calib_keyp_dist.txt fname_raw_calib_keyp_corr.txt [fname_img_n.pgm fname_img_n_red_corr.pgm fname_img_n_green_corr.pgm fname_img_n_blue_corr.pgm, ...]`  
* Runs all the circuit: first it builds the correction polynomials based on calibration raw image, it corrects the calibration image and then it corrects all the images that were taken with the same camera settings. In txt files, it saves the geometrical coordinates of the calibration image keypoints in all the channels - before and after the correction. Those txt files can be used later for visualization tests in Matlab (i.e., misalignment field and histograms). The correction polynomials are not saved, for this, use point *three-input-arguments run*.  
* **EXAMPLE**: `data/_MG_7626.pgm data/_MG_7626_R_corr.pgm data/_MG_7626_G_corr.pgm data/_MG_7626_B_corr.pgm data/_MG_7626_dist.txt data/_MG_7626_corr.txt data/_MG_7628.pgm data/_MG_7628_R_corr.pgm data/_MG_7628_G_corr.pgm data/_MG_7628_B_corr.pgm`  

#### *29 Apr 2026*
- build with Visual Studio 2022
  - by Open -> Folder
- process `data/_MG_7626.pgm` by default (simplify Visual Studio debugging)
  - 'LMA center redefinition' runs slowly...
  - calculates x,y polynomial coefficients to `x^11`, `y^11`;
	- should check coefficients and stop when exceedingly small...

#### *30 Apr 2026*
- employ [DeepWiki](https://docs.devin.ai/work-with-devin/deepwiki)


## Supported Image Format (from [DeepWiki](https://deepwiki.com/blekenbleu/CorrCA))
Portable Pixmap (PPM) image files are not supported in this codebase,  
which specifically works with Portable Gray Map (PGM) format files. 

CorrCA supports *only* PGM (Portable Gray Map) format for both input and output operations:

- **Input files**: Raw calibration and correction images use `.pgm` extension (e.g., `_MG_7626.pgm`)  
- **Output files**: Corrected channel images are saved as separate PGM files (e.g., `_R_corr.pgm`, `_G_corr.pgm`, `_B_corr.pgm`)  

### Implementation Details

Image I/O is implemented through:
- `read_pgm_image_double()` for reading PGM files  
- `write_pgm_image_double()` for writing PGM files  

It processes Bayer-patterned images *where color components are interleaved in a single grayscale channel*.

### Notes

The codebase processes raw Bayer pattern images where R, G, and B components are interleaved in a single channel.  
This is separated into individual color channels during processing, maintaining PGM format for all file operations.  
There are no references to PPM format support anywhere in the source code or documentation.

DeepWiki pages to explore:
- [Project Layout and Data Files](https://deepwiki.com/blekenbleu/CorrCA/1.2-project-layout-and-data-files)
- [Glossary](https://deepwiki.com/blekenbleu/CorrCA/5-glossary)

#### Citations

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

#### Citations

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


