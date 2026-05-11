# Modified Radial Symmetry
CorCA generally fails for polynomial degrees less than 11, because  
- polynomial models are generally problematic
  - bicubic splines are typically employed
- Optical lens aberration are often nearly radially symmetric
  - better to convert XY distances to radii from some center
- typical lens correction software does not include centering
  - it also does not consider ellipsoidal distortion

Even if images are captured from perfectly centered optics,  
removing Bayer matrix from raw images may leave blue and red planes  
*slightly shifted from green*.  
Subjects may be other than perfectly perpendicular to lens axes,  
confounding transverse CA with longitudinal.  

### workaround
- first, perform linear regression for each red and blue delta column vs green x, y
 - determine most significant constant coefficients from those models for red and blue  
   **This is a likely center for radial symmetry**  
 - for most significant x and y coefficients, *calculate their ratio*  
   **This is a likely ellipsoidal distortion**
- rescaling x or y should reduce ellipsoidal distortion
- calculate radii from regression center to x and y independent variables
- perform linear regression for red and blue delta columns based on calculated radii
  - consider [pseudoinverse](https://jabrekh.github.io/assets/materials/LA_Spring2023_L15-QRInverse.pdf) to speed regressions
	- in that paper, The set of real numbers is denoted by &#x211D;.
    - QR decomposition by the Gram-Schmidt process is inherently numerically unstable,  
      but relatively easy to implement.
    - [GitHub pseudoinverse `source.cpp` comparing QR Factorization and SVD](https://github.com/ggory15/Pseudo_inverse_using_QR_Factorizaion)
		- uses [Eigen](https://libeigen.gitlab.io/);&nbsp; latest stable release is Eigen 5.0.0.
			- [Building Eigen 3.3.8 with Visual Studio](https://github.com/Qannaf/Building-PCL-with-Visual-Studio-from-source-on-Windows/blob/main/Eigen3.3.8.md)
			- [Configure the EIGEN 3.4.0 library in Visual Studio](https://www.programmersought.com/article/568011439760/)
			- ["`this->solve(rhs)` is more efficient and numerically stable than `this->pseudoInverse()*rhs`"](https://libeigen.gitlab.io/eigen/docs-5.0/classEigen_1_1CompleteOrthogonalDecomposition.html)
  - [iterative Conjugate Gradient Linear System Solver](https://github.com/MB-VGU-projects/cpp-linalg-regression-1)
