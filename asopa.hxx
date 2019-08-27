#pragma once
#include <Eigen/Dense>

/*
  The solution is based on:
  Mohammed Bennani Dosse, Jos Ten Berge
  "Anisotropic Orthogonal Procrustes Analysis"
  Journal of Classification
  June 2010, Volume 27, pp 111-128

  Ansiotropically Scaled Orthogonal Procrustes Algorithm
  Finds 3D Transformation from X to Y with known correspondence

  Inputs: 3xn Matrix           - X
          3xn Matrix           - Y
          Accuracy of estimate - threshold

  Outputs: 3x3 Rotation Matrix         - Q,
           3x3 Scaling Matrix          - A
           1x3 Translation vector      - t
	   Root Mean Square Error      - RMSE
	       
*/

int asopa(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	  double threshold,
	  Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	  double &RMSE);
