[*back*](README.md#9-may-2026-ppm-branch)
## Chromatic Aberration [Bicubic](https://maxcandocia.com/article/2024/Oct/14/sampled-bicubic-spline-fitting/) Polynomial Linear Regression Models
Max Candocia made [BiCubic **Spline** proof-of-concept python and R code available on GitHub](https://github.com/mcandocia/bicubic-spline-fitting)

Rather than trying to convert that to C,  
instead try a [multiple linear regression model](https://github.com/oscar8880/Multiple-Linear-Regression),
with a bicubic polynomial kernel:  
 &emsp; `CAxy = a + b*x + c*x*x + d*x*x*x + e*y + f*y*y + g*y*y*y`  

While more complex than [typical radially symmetric lens Transverse Chromatic Aberration](https://lensfun.github.io/calibration-tutorial/lens-tca.html) model:  
 &emsp; `CAr = a * r^4 + b * r^3 + c * r^2 + v * r`

.. it anticipates photomicography issues
- possible stage tilt, introducing some LoCA
- imperfect, misaligned and uncentered optics

## Linear regression, polynomial radial vs Cartesian fit
using e.g. [JASP](https://github.com/blekenbleu/Multiple-Linear-Regression/blob/CA/JASP.md) for linear regression
- Cartesian may better match photomicography optical physics and geometry
- linear regression to first estimate radius center and possible x-y scaling distortion
```
dxR = 0.826 - 1.614 xG - 0.017 yG
dyR = 0.589 - 0.018 xG - 1.148 yG
dxB = -0.992 + 0.113 xG - 0.023 yG
dyB = -1.049 - 0.031 xG + 0.132 yG
```
[xG,yG center estimation discussion](data/Gxyfit.txt)

### radially symmetric red-green and blue-green center difference data fit
JASP linear regression has radius-squared more significant than radius or radius-cubed.  
That simplifies solving with unknown x,y center A,B:  
`radius_squared r2 = (x-A)**2 + (y-B)**2 = x**2 - 2Ax + A**2 + y**2 - 2Bx + B**2`  

With unknown A and B, then instead fitting polynomial model coefficients:  
`polynomial xy2 = C0 + C1*x + C2*y + C3*x*x + C4*y*y + C5*x*x*x + C6*y*y*y`  
.. should simultaneously solve for radial centers and other x,y distortions,  
as is supported by statistics, with `r2` models based on A,B center guesstimates for dyR:

| Model  | RMSE  | F     | Tintercept| Std Error |
|--------|-------|-------|-----------|-----------|
| [dxBr2](data/dxBr2.png)| 0.103 | 88.87 | -27.68    | 0.094 |
| [dxBxy2](data/dxBxy2.png) | 0.093 | 128.1 | -26.85    | 0.024 |
| [dxRr2](data/dxRr2.png)  | 0.078 | 6223  | 20.48     | 0.071 |
| [dxRxy2](data/dxRxy2.png) | 0.074 | 8465  | 46.65     | 0.015 |
| [dyBr2](data/dyBr2.png)  | 0.101 | 36.34 | -16.76    | 0.092 |
| [dyBxy2](data/dyBxy2.png) | 0.100 | 54.36 | -59.25    | 0.017 |
| [dyRr2](data/dyRr2.png)  | 0.075 | 3524  | 13.48     | 0.068 |
| [dyRxy2](data/dyRxy2.png) | 0.075 | 4410  | 44.86     | 0.013 |

.. with only `dxBxy2 Tintercept` worse than corresponding r2
