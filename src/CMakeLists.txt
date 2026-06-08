cmake_minimum_required(VERSION 3.10)

project(Aberration_Correction)

add_executable(CAcorr distortion.cpp image.cpp main_centering.cpp misc.cpp
	numerics.cpp pgm_io.cpp spline.cpp vector.cpp deBayer.cpp regress.cpp
	report.cpp reduceCA.cpp matrix_correction.cpp matrix_pnm.cpp)
