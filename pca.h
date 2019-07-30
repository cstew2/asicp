#pragma once
#include <Eigen/Dense>

/*
  Use Principal Component Analysis to solve for the 
  approximite scaling and rotation relationship between 
  the two points sets.
 */

int pca_registration(Eigen::MatrixXd X, Eigen::MatrixXd Y,
		     Eigen::Matrix3d &Q,
		     Eigen::Matrix3d &A,
		     Eigen::Vector3d &t,
		     double &RMSE);
