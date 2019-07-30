#pragma once
#include <Eigen/Dense>

/*
  The solution is based on:
  Mohammed Bennani Dosse and Jos Ten Berge (2010),
  "Anisotropic Orthogonal Procrustes Analysis"
  Journal of Classification 27:111-128
    
  Ansiotropically Scaled Orthogonal Procrustes Algorithm:
  Finds 3D Transformation from X to Y. Note that this is
  algorithm finds the prescaling solution, i.e. the 
  scaling happens before the rotation. This algorithm
  uses the majorisation principle to solve for the rotation 
  matrix not block relaxation.
      
  Inputs: 3xn Matrix                 - X
          3xn Matrix                 - Y
	  Accuracy of transformation - threshold

  Outputs: 3x3 Rotation Matrix         - R
           3x3 Scaling Matrix          - A
           1x3 Translation vector      - t
           Root Mean Square Error      - RMSE
*/

int asopa(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	  double threshold, bool estimate,
	  Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	  double &RMSE);
