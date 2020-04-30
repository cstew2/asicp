#pragma once
#include <Eigen/Dense>

/*
  Use Principal Component Analysis to solve for the 
  approximite scaling relationship between 
  the two points sets.
 */

int pca_scales(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	       Eigen::Matrix3d &R,
	       Eigen::Matrix3d &A);
