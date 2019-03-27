#include <iostream>

#include "asicp.h"
#include "asopa.h"

int asicp(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	   double threshold, size_t max_iterations,
	   Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::RowVector3d &t,
	   double &FRE)
{
	if(X.cols() != 3 || Y.cols() != 3) {
		std::cerr << "X and Y must have be 3 row matrices" << std::endl;
		return -1;	
	}

	size_t Xn = X.rows();
	size_t Yn = Y.rows();
	
	Eigen::Vector3d X_centroid = X.colwise().mean();
	Eigen::Vector3d Y_centroid = Y.colwise().mean();

	//translate input by centroids
	Eigen::MatrixXd X_trans = X.rowwise() - X_centroid.transpose();
	Eigen::MatrixXd Y_trans = Y.rowwise() - Y_centroid.transpose();
	
	//find initial rotation and scaling factors
	Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d S = Eigen::Matrix3d::Identity();
	
	//transformed source and target for corespondence
	Eigen::MatrixXd X_t = X_trans;
	//scaled row of X for finding corespondence
	Eigen::RowVector3d x_s;
	//scaled Y points for finding corespondence
	Eigen::MatrixXd Y_s;
	//index of closest point in Y to point in x
	Eigen::MatrixXd::Index index;
	//closest points of Y to points in X
	Eigen::MatrixXd Y_close;

	//asopa variables
	Eigen::RowVector3d tt(3);
	double asopa_FRE = 0.0f;
	Eigen::MatrixXd asopa_FRE_mag = Eigen::MatrixXd();

	//error variables
	FRE = 0.0f;
	Eigen::MatrixXd FRE_vect(Xn, 3);
	
	for(size_t i=0; i < max_iterations; i++) {
		//current transformed X points
		X_t = (Q * A * X_trans.transpose()).transpose().rowwise() + t;
		//find corrispodence of points of X to points of Y
		//scale all Y points
		for(size_t i=0; i < Yn; i++) {
			Y_s.row(i) = Y_trans.row(i) * S.diagonal().transpose();
		}

		//find mahalanobis distance between X and Y points
		for(size_t i=0; i < Xn; i++) {
			//scale X points
		        //x_s = X_t.row(i) * S.diagonal().transpose();
			//(Y_s.rowwise() - x_s).rowwise().squaredNorm().minCoeff(&index);
			//Y_close.row(i) = Y.row(index);
		}
		
		//find registration
		//asopa(X_t, Y_trans, threshold, R, S, tt, asopa_FRE);		
		
		//calculate FRE
		//FRE_vect = Y_close - X_t * A * Q;
		//FRE = sqrt(FRE_vect.squaredNorm()/Xn);

		//check if FRE is low enough to end
		if(FRE < threshold) {
			break;
		}
		
	}
	
	//find final error
	//FRE_vect = Y_close - ((X_trans * A) * Q);
	//FRE = sqrt(FRE_vect.squaredNorm()/Xn);
	
	t =  Y_centroid - (Q * (A * X_centroid));
	
	return 0;
}

