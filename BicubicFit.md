[*back*](README.md#9-may-2026-ppm-branch)
## Chromatic Aberration [Bicubic Spline](https://maxcandocia.com/article/2024/Oct/14/sampled-bicubic-spline-fitting/),  hybrid Linear Regression Models
Instead of high order polynomials, consider fitting
## a 2x2 Bicubic Spline  
for which Max Candocia made [proof-of-concept python and R code available on GitHub](https://github.com/mcandocia/bicubic-spline-fitting)

Instead of trying to convert that to C,  
instead try a [multiple linear regression model](https://github.com/oscar8880/Multiple-Linear-Regression),
then replace its linear kernel with a bicubic:  
 &emsp; `CAxy = a + b*x + c*x*x + d*x*x*x + e*y + f*y*y + g*y*y*y + h*y*x`  

While that is more complex than [typical radially symmetric lens Transverse Chromatic Aberration](https://lensfun.github.io/calibration-tutorial/lens-tca.html) models:  
 &emsp; `CAr = a * r^4 + b * r^3 + c * r^2 + v * r`

.. it anticipates photomicography issues
- possible stage tilt, introducing some LoCA
- imperfect, misaligned and uncentered optics

## a hybrid model with Linear regression, then polynomial radius fit
using e.g. [JASP](https://github.com/blekenbleu/Multiple-Linear-Regression/blob/CA/JASP.md) for linear regression
- should better match photomicography optical physics and geometry
- linear regression to sort radius center and possible x-y scaling distortion
- polynomial fit to radius
