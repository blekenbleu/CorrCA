[*back*](README.md#9-may-2026-ppm-branch)
## [Bicubic Spline Fitting](https://maxcandocia.com/article/2024/Oct/14/sampled-bicubic-spline-fitting/)
Instead of high order polynomials, consider fitting a 2x2 Bicubic Spline
Max Candocia made [proof-of-concept python and R code available on GitHub](https://github.com/mcandocia/bicubic-spline-fitting)

Instead of trying to convert that to C,  
instead try a [multiple linear regression model](https://github.com/oscar8880/Multiple-Linear-Regression),
then replace its linear kernel with a bicubic:  
 &emsp; `CAxy = a + b*x + c*x*x + d*x*x*x + e*y + f*y*y + g*y*y*y + h*y*x`  

While that is more complex than [typical radially symmetric lens Transverse Chromatic Aberration](https://lensfun.github.io/calibration-tutorial/lens-tca.html) models:  
 &emsp; `CAr = a * r^4 + b * r^3 + c * r^2 + v * r`

.. it anticipates photomicography issues
- possible stage tilt, introducing some LoCA
- imperfect, misaligned and uncentered optics

