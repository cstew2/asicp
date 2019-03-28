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

	size_t Xn = X.cols();
	size_t Yn = Y.cols();
	
	Eigen::Vector3d X_centroid = X.rowwise().mean();
	Eigen::Vector3d Y_centroid = Y.rowwise().mean();

	//translate input by centroids
	Eigen::MatrixXd X_trans = X.colwise() - X_centroid;
	Eigen::MatrixXd Y_trans = Y.colwise() - Y_centroid;
	
	//find initial rotation and scaling factors
	Q = Eigen::Matrix3d::Identity();
	A = Eigen::Matrix3d::Identity();
	
	//transformed source and target for corespondence
	Eigen::MatrixXd X_t(3, Xn);
	//scaled row of X for finding corespondence
	Eigen::Vector3d x_s(3);
	//scaled Y points for finding corespondence
	Eigen::MatrixXd Y_s(3, Yn);
	//index of closest point in Y to point in x
	Eigen::MatrixXd::Index index;
	//closest points of Y to points in X
	Eigen::MatrixXd Y_close(3, Xn);

	//asopa variables
	Eigen::Vector3d tt(3);
	double asopa_FRE = 0.0f;
	Eigen::MatrixXd asopa_FRE_mag = Eigen::MatrixXd();

	//error variables
	FRE = 0.0f;
	Eigen::MatrixXd FRE_vect(Xn, 3);
	
	for(size_t i=0; i < max_iterations; i++) {
		//current transformed X points
		X_t = (Q * (A * X)).colwise() + t;
		//find corrispodence of points of X to points of Y
		//scale all Y points
		for(size_t i=0; i < Yn; i++) {
			Y_s(0, i) = Y_trans(0, i) * Q(0, 0);
			Y_s(1, i) = Y_trans(1, i) * Q(1, 1);
			Y_s(2, i) = Y_trans(2, i) * Q(2, 2);
		}

		//find mahalanobis distance between X and Y points
		for(size_t i=0; i < Xn; i++) {
			//scale X points
			x_s(0) = X_t(0, i) * Q(0, 0);
			x_s(1) = X_t(1, i) * Q(1, 1);
			x_s(2) = X_t(2, i) * Q(2, 2);
			(Y_s.colwise() - x_s).colwise().squaredNorm().minCoeff(&index);
			Y_close.col(i) = Y.col(index);
		}
		//find registration
		asopa(X_t, Y_close, threshold, Q, A, tt, asopa_FRE);		
		
		//calculate FRE
		FRE_vect = Y_close - Q * (A * X_t);
		FRE = sqrt(FRE_vect.squaredNorm()/Xn);

		//check if FRE is low enough to end
		if(FRE < threshold) {
			break;
		}
		
	}

	//find translation
	t =  Y_centroid - (Q * (A * X_centroid));
	
	//find final error
	FRE_vect = Y_close - ((Q * (A * X)).colwise() + t);
	FRE = sqrt(FRE_vect.squaredNorm()/Xn);
	
	return 0;
}

