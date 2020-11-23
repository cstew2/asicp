#pragma once
#include <Eigen/Dense>

/*
  The solution is based on:
  Elvis C.S. Chen, A. Jonathan McLeod, John S.H. Baxter, Terry M. Peters,
  "Registration of 3D shapes under anisotropic scaling"
  International Journal of Computer Assisted Radiology and Surgery
  June 2015, Volume 10, Issue 6, pp 867–878

  Ansiotropically Scaled Iterative Closest Point Algorithm
  Finds the 3D transformation from X to Y

  Inputs: 3xn Matrix           - X
          3xn Matrix           - Y
          Accuracy of estimate - threshold
          Maximum interations  - max_iterations
	  Esimates given       - estimate

  Outputs: 3x3 Rotation Matrix                               - Q,
           3x3 Scaling Matrix                                - A
           1x3 Translation vector                            - t
	   Indices of Y that are cloeset to points in X  - indices
	   Root Mean Square Error                            - RMSE
	       
*/


int asicp(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	  double threshold, size_t max_iterations, double asopa_threshold,
	  bool estimate, int rotations,
	  Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	  std::vector<size_t> &indices, double &RMSE);

int asicp_rot(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	      double threshold, size_t max_iterations, double asopa_threshold,
	      Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	      std::vector<size_t> &indices, double &RMSE);
