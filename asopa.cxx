#include <iostream>
#include "asopa.hxx"

int asopa(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	  double threshold,
	  Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	  double &RMSE)
{
	//check input dimensions
	if(X.rows() != 3 || Y.rows() != 3) {
		std::cerr << "X and Y must be 3xn matrices" << std::endl;
		return -1;
	}
	if(X.cols() == 0 || Y.cols() == 0) {
		std::cerr << "No points given" << std::endl;
		return -1;
	}
	if(Y.cols() != X.cols()) {
		std::cerr << "X and Y are required to have the same number of points" << std::endl;
		return - 1;
	}
	if(X.cols() == 1 || Y.cols() == 1) {
		Q = Eigen::Matrix3d::Identity();
		A = Eigen::Matrix3d::Identity();
		t = (Y.col(0) - X.col(0)).col(0);
		return 0;
	}
	
	Q = Eigen::Matrix3d::Identity();
	A = Eigen::Matrix3d::Identity();
	t = Eigen::Vector3d::Zero();
	
	//number of points
	size_t n = X.cols();
	
	//find centroids
	Eigen::Vector3d X_centroid = X.rowwise().mean();
	Eigen::Vector3d Y_centroid = Y.rowwise().mean();

	//translate input by centroids
	Eigen::MatrixXd X_trans = X.colwise() - X_centroid;
	Eigen::MatrixXd Y_trans = Y.colwise() - Y_centroid;
	
	//normalise translated source input
	Eigen::MatrixXd X_norm = X_trans.rowwise().normalized();
	
	//Find the cross covariance matrix
	Eigen::Matrix3d B = Y_trans * X_norm.transpose();

	//use SVD to decompose B such that B=USV^T
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(B, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d U = svd.matrixU();
	Eigen::Matrix3d V = svd.matrixV().transpose();

	//create D where the last element is the determinant of UV
     	Eigen::Matrix3d D(3, 3);
	D = Eigen::Matrix3d::Identity();
	D(2, 2) = (U*V).determinant();

	//Find the rotation
	Q = U * D * V;

	//Calculate the RMSE for the given rotation
	Eigen::MatrixXd RMSE_vect(3, n);
	RMSE_vect = Y_trans - Q.transpose() * X_norm;
	RMSE = sqrt(RMSE_vect.squaredNorm()/n);

	//starting RMSE value to compare to
	double RMSE_orig = 2.0 * (RMSE + threshold);

	//majorisation to solve for rotation
       	while(fabs(RMSE_orig - RMSE) > threshold) {
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
		
		//recompute RMSE value
		RMSE_vect = Y_trans - Q.transpose() * X_norm;
		RMSE_orig = RMSE;
		RMSE = 0.0f;
		RMSE = sqrt(RMSE_vect.squaredNorm()/n);
	}

	//calculate final scaling
	B = Y_trans * X_trans.transpose();
	U = B.transpose() * Q;
	V = X_trans * X_trans.transpose();
	A = Eigen::Matrix3d::Zero();
	A(0, 0) = U(0, 0)/V(0, 0);
	A(1, 1) = U(1, 1)/V(1, 1);
	A(2, 2) = U(2, 2)/V(2, 2);

	//calculate final translation
	t =  Y_centroid - (Q * (A * X_centroid));
	
	//calculate final RMSE values
	RMSE_vect = Y - ((Q * (A * X)).colwise() + t);
	RMSE = sqrt(RMSE_vect.squaredNorm()/n);
	
	return 0;
}
