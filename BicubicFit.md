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
```
dxR = 0.826 - 1.614 xG - 0.017 yG
dyR = 0.589 - 0.018 xG - 1.148 yG
dxB = -0.992 + 0.113 xG - 0.023 yG
dyB = -1.049 - 0.031 xG + 0.132 yG
```
[xG,yG center estimation discussion](data/Gxyfit.txt)
- polynomial fit to radius
