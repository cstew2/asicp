#include <iostream>

#include "asicp.h"
#include "asopa.h"

int asicp(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	  double threshold, size_t max_iterations,
	  Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	  double &FRE)
{
	if(X.rows() != 3 || Y.rows() != 3) {
		std::cerr << "X and Y must be column matrices with 3 rows" << std::endl;
		return -1;	
	}

	//size of inputs
	size_t Xn = X.cols();
	size_t Yn = Y.cols();

	//find centroids of source and target
	Eigen::Vector3d X_centroid = X.rowwise().mean();
	Eigen::Vector3d Y_centroid = Y.rowwise().mean();

	//translate input by centroids
	Eigen::MatrixXd X_trans = X.colwise() - X_centroid;
	Eigen::MatrixXd Y_trans = Y.colwise() - Y_centroid;
	
	//current and previous iteration of X coords
	Eigen::MatrixXd X_curr(3, Xn);
	Eigen::MatrixXd X_prev(3, Xn);
	Eigen::Vector3d x_scale(3);
	
	//ordered closest points in Y_trans to points in X_curr
	Eigen::MatrixXd Y_close(3, Yn);
	Eigen::MatrixXd Y_scale(3, Yn);

	Eigen::MatrixXd::Index Y_index;
	
	//Set translation to identity
	Q = Eigen::Matrix3d::Identity();
	A = Eigen::Matrix3d::Identity();
	t = Eigen::Vector3d::Zero();

	//transformation between iterations
	Eigen::Matrix3d Q_mod;
	Eigen::Matrix3d A_mod;
	Eigen::Vector3d t_mod;

	Q_mod = Eigen::Matrix3d::Identity();
	A_mod = Eigen::Matrix3d::Identity();
	t_mod = Eigen::Vector3d::Zero();
	
	double FRE_prev = sqrt((Y_trans - X_trans).squaredNorm()/Yn);
	
	double FRE_asopa = 0.0;
	
	for(size_t j=0; j < max_iterations; j++) {
		//current transformation
		X_curr = (Q * (A * X_trans));

		//find correspondence
		for(size_t i=0; i < Yn; i++) {
			//scale target
			Y_scale(0,i) = Y_trans(0,i)/sqrt(A(0,0));
			Y_scale(1,i) = Y_trans(1,i)/sqrt(A(1,1));
			Y_scale(2,i) = Y_trans(2,i)/sqrt(A(2,2));
		}

		for(size_t i=0; i < Xn; i++) {
			//scale source
			x_scale(0) = Y_trans(0,i)/sqrt(A(0,0));
			x_scale(1) = Y_trans(1,i)/sqrt(A(1,1));
			x_scale(2) = Y_trans(2,i)/sqrt(A(2,2));
			
			//calculate shortest mahalanobis distance to find
			//pair of points between X and Y
			(Y_scale.colwise() - x_scale).colwise().squaredNorm().minCoeff(&Y_index);
			Y_close.col(i) = Y_trans.col(i);
		}
		
		//find transformation
		asopa(X_curr, Y_close, threshold, Q_mod, A_mod, t_mod, FRE_asopa);
		Q = Q * Q_mod;
		A = A * A_mod;
		
		//calculate error

		
		//set prev to current
		X_prev = X_curr;
	}

	//calculate transform using centroids with rotation and scaling
	t =  Y_centroid - (Q * (A * X_centroid));
	
	return 0;
}

void initial_scales(Eigen::MatrixXd X, Eigen::MatrixXd Y,
		    Eigen::Matrix3d Q,
		    double &scale)
{
	size_t Xn = Y.cols();
	size_t Yn = Y.cols();

	//find centroids of each point set
	Eigen::Vector3d X_centroid = X.rowwise().mean();
	Eigen::Vector3d Y_centroid = Y.rowwise().mean();

	//translate input by centroids
	Eigen::MatrixXd X_trans = X.colwise() - X_centroid;
	Eigen::MatrixXd Y_trans = Y.colwise() - Y_centroid;

	Eigen::MatrixXd X_t = Q * X_trans;
		
	Eigen::MatrixXd A = 1.0/(double)Xn * (X_t * X_t.transpose());
	Eigen::MatrixXd B = 1.0/(double)Yn * (Y_trans * Y_trans.transpose());

	
}
