cmake_minimum_required(VERSION 3.10)

project(Aberration_Correction)

add_executable(CAcorr src/distortion.cpp src/image.cpp src/main_centering.cpp src/misc.cpp
	src/numerics.cpp src/pgm_io.cpp src/spline.cpp src/vector.cpp src/deBayer.cpp src/regress.cpp
	src/report.cpp src/reduceCA.cpp src/matrix_correction.cpp src/matrix_pnm.cpp
	src/test_pnm_matrix.cpp src/pseudoinverse.cpp)
