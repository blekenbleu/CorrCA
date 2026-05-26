cmake_minimum_required(VERSION 3.10)

project(Aberration_Correction)

add_executable(CAcorr abberation.h centers.h correction.h distortion.cpp
	 distortion.h image.cpp image.h LM.h main_centering.cpp matrix.cpp array.h
	 misc.h misc.cpp numerics.cpp numerics.h pgm_io.cpp pgm_io.h spline.cpp
	 spline.h vector.cpp t_test.h correct.cpp deBayer.cpp)
