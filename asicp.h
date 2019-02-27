#pragma once
#include <Eigen/Dense>

 /*
      The solution is based on:
      Elvis C.S. Chen, A. Jonathan McLeod, John S.H. Baxter, Terry M. Peters,
      "Registration of 3D shapes under anisotropic scaling"
      International Journal of Computer Assisted Radiology and Surgery
      June 2015, Volume 10, Issue 6, pp 867–878

      Ansiotropically Scaled Orthogonal Procrustes Algorithm
      Finds 3D Transformation from X to Y
      Inputs: 3xn Matrix - X
              3xn Matrix - Y
	      Error threshold to halt iteration - threshold
	      Upper bound on number of iterations - max_iterations
      Outputs: 3x3 Rotation Matrix - R,
               3x3 Scaling Matrix  - A
	       1x3 Translation vector - t
	       
    */

void asicp(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	   double threshold, size_t max_iterations,
	   Eigen::Matrix3d &R, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	   double &FRE);
