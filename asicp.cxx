#include <iostream>

#include "asicp.h"
#include "asopa.h"

void asicp(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	   double threshold, size_t max_iterations,
	   Eigen::Matrix3d &R, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	   double &FRE)
{
	if(X.cols() != 3 || Y.cols() != 3) {
		std::cerr << "X and Y must have be 3 row matrices" << std::endl;
		return;	
	}

	size_t Xn = X.rows();
	size_t Yn = Y.rows();
	
	R = Eigen::Matrix3d::Identity();
	A = Eigen::Matrix3d::Identity();

	Eigen::Matrix3d B;
	Eigen::Vector3d B_mag;

	FRE = 0.0f;
	double FRE_orig = 0.0f;
	Eigen::MatrixXd FRE_mag;

	Eigen::Matrix3d AA;
	Eigen::Matrix3d RR;
	
	
	for(size_t i=0; i < max_iterations; i++) {
	
		
	}
}
