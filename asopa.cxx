#include <iostream>
#include "asopa.h"

void asopa(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	   double threshold,
	   Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	    double &FRE, Eigen::MatrixXd &FRE_mag)
{
	//check input dimensions
	if(X.cols() != 3 || Y.cols() != 3) {
		std::cerr << "X and Y must be 3xn matrices" << std::endl;
		return;	
	}
	if(Y.rows() != X.rows()) {
		std::cerr << "X and Y are required to have the same number of points" << std::endl;
		return;
	}

	//number of points
	size_t n = X.rows();
	
	//find centroids
	Eigen::Vector3d X_centroid = X.colwise().mean();
	Eigen::Vector3d Y_centroid = Y.colwise().mean();

	//translate input by centroids
	Eigen::MatrixXd X_trans = X.rowwise() - X_centroid.transpose();
	Eigen::MatrixXd Y_trans = Y.rowwise() - Y_centroid.transpose();
	
	//normalise translated source input
	Eigen::MatrixXd X_norm = X_trans.colwise().normalized();
	
	//Find the cross covariance matrix
	Eigen::Matrix3d B = Y_trans.transpose() * X_norm;

	//use SVD to decompose B such that B=USV^T
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(B, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d U = svd.matrixU();
	Eigen::Matrix3d V = svd.matrixV().transpose();

	//create D where the last element is the determinant of UV
     	Eigen::Matrix3d D(3, 3);
	D = Eigen::Matrix3d::Identity();
	D(2, 2) = (U*V).determinant();

	//Find the rotation
	Q = V * D * U;
	std::cout << "Q:" << std::endl << Q << std::endl << std::endl;

	//Calculate the FRE for the given rotation
	Eigen::MatrixXd FRE_vect(n, 3);
	FRE_vect = X_norm * Q.transpose() - Y_trans;
	FRE = sqrt(FRE_vect.squaredNorm()/n);

	//starting FRE value to compare to
	double FRE_orig = 2.0 * (FRE + threshold);
	
	//majorisation to solve for rotation
       	while(fabs(FRE_orig - FRE) > threshold) {
		//recompute SVD values
		Eigen::Matrix3d QB = Q.transpose() * B;
		Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
		I(0,0) = QB(0,0);
		I(1,1) = QB(1,1);
		I(2,2) = QB(2,2);
		
		svd.compute(B * I, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		V = svd.matrixV().transpose();
		//recompute rotation
		D(2, 2) = (U*V).determinant();
		Q = U * D * V;
		
		//recompute FRE value
		FRE_vect = X_norm * Q.transpose() - Y_trans;
		FRE_orig = FRE;
		FRE = 0.0f;
		FRE = sqrt(FRE_vect.squaredNorm()/n);
	}

	//calculate final scaling
	B = Y_trans.transpose() * X_trans;
	U = B.transpose() * Q;
	V = X_trans.transpose() * X_trans;
	A = Eigen::Matrix3d::Zero();
	A(0, 0) = U(0, 0)/V(0, 0);
	A(1, 1) = U(1, 1)/V(1, 1);
	A(2, 2) = U(2, 2)/V(2, 2);

	//calculate final FRE values
	FRE_vect = Y_trans - (X_trans * A * Q.transpose());
	FRE_mag = FRE_vect.rowwise().squaredNorm();
	FRE = sqrt(FRE_vect.squaredNorm()/n);

	//calculate final translation
	t =  Y_centroid - (Q * A * X_centroid);
}
