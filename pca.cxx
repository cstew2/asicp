#include <iostream>

#include "pca.hxx"

int pca_scales(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	       Eigen::Matrix3d &R,
	       Eigen::Matrix3d &A)
{
	size_t Xn = X.cols();
	size_t Yn = Y.cols();
	
	A = Eigen::Matrix3d::Identity();
       		
	//translate input by centroids (data has a mean of 0)
	Eigen::MatrixXd X_trans = X.colwise() - X.rowwise().mean();
	Eigen::MatrixXd Y_trans = Y.colwise() - Y.rowwise().mean();

	//rotate X by R
	Eigen::MatrixXd X_rot = R * X_trans;
	
	//computer cross covariance matrix
	//C= X'Y/n-1 (use un-scaled for esimating scaling factors?)
	Eigen::Matrix3d Sx = 1/(double)Xn * (X_rot * X_rot.transpose());
	Eigen::Matrix3d Sy = 1/(double)Yn * (Y_trans * Y_trans.transpose());
	
	//solve for eigenvectors and eigenvalues
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
	es.compute(Sx);
	Eigen::Matrix3d X_evec = es.eigenvectors();
	Eigen::Vector3d X_eval = es.eigenvalues();
	es.compute(Sy);
	Eigen::Matrix3d Y_evec = es.eigenvectors();
	Eigen::Vector3d Y_eval = es.eigenvalues();
	
	A.diagonal() = Eigen::Vector3d(X_eval(0)/Y_eval(0),
				       X_eval(1)/Y_eval(1),
				       X_eval(2)/Y_eval(2));

	//std::cout << X_eval(0) << " " << Y_eval(0) << " = " << X_eval(0)/Y_eval(0) << std::endl;
	//std::cout << X_eval(1) << " " << Y_eval(1) << " = " << X_eval(1)/Y_eval(1) << std::endl;
	//std::cout << X_eval(2) << " " << Y_eval(2) << " = " << X_eval(2)/Y_eval(2) << std::endl;

	//std::cout << X_eval(0)/Y_eval(0)+X_eval(1)/Y_eval(1)+X_eval(2)/Y_eval(2) << std::endl;
	
	return 0;
}
