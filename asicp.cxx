#include <iostream>

#include <nanoflann.hpp>

#include "asicp.hxx"
#include "asopa.hxx"
#include "pca.hxx"

int asicp(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	  double threshold, size_t max_iterations, double asopa_threshold, bool estimate,
	  Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	  double &RMSE)
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
	
	//current iteration of X points
	Eigen::MatrixXd X_curr(3, Xn);
	//scaled X points from X_curr
	Eigen::MatrixXd X_scale(3, Xn);
	
	//ordered closest points in Y_trans to points in X_curr
	Eigen::MatrixXd Y_close(3, Yn);
	//scaled Y points
	Eigen::MatrixXd Y_scale(3, Yn);

	//set up kd-tree
	typedef nanoflann::KDTreeEigenMatrixAdaptor<
		Eigen::MatrixXd, -1, nanoflann::metric_L2>
		kd_tree_t;
	
	const int max_leaves = 10;

	if(!estimate) {
		Q = Eigen::Matrix3d::Identity();
		A = Eigen::Matrix3d::Identity();
		t = Eigen::Vector3d::Zero();
	}
	else {
		pca_registration(X, Y, Q, A, t, RMSE);
	}
	
	//transformation between iterations
	Eigen::Matrix3d Q_mod;
	Eigen::Matrix3d A_mod;
	Eigen::Vector3d t_mod;

	//initialise kd-tree of Y
	
	double RMSE_prev = sqrt((Y_trans - X_trans).squaredNorm()/Yn);
	double RMSE_asopa = 0.0;
	
	for(size_t j=0; j < max_iterations; j++) {
		std::cout << "iteration " << j << std::endl;
		std::cout << "Q" << std::endl << Q << std::endl;
		std::cout << "A" << std::endl << A << std::endl;
		std::cout << "t" << std::endl << t << std::endl;
		
		//current transformation
		X_curr = (Q * (A * X_trans));

		//find correspondence
		for(size_t i=0; i < Yn; i++) {
			//scale target
			Y_scale(0,i) = Y_trans(0,i)/sqrt(A(0,0));
			Y_scale(1,i) = Y_trans(1,i)/sqrt(A(1,1));
			Y_scale(2,i) = Y_trans(2,i)/sqrt(A(2,2));
		}
		kd_tree_t kd_tree(Yn, std::cref(Y_scale), max_leaves);
	        kd_tree.index->buildIndex();
	        
		for(size_t i=0; i < Xn; i++) {
			//scale source
			X_scale(0, i) = Y_trans(0,i)/sqrt(A(0,0));
			X_scale(1, i) = Y_trans(1,i)/sqrt(A(1,1));
			X_scale(2, i) = Y_trans(2,i)/sqrt(A(2,2));
		}
		const size_t num_results = 3;
		std::vector<size_t> ret_indexes(num_results);
		std::vector<double> out_dists_sqr(num_results);

		for(size_t i=0; i < Xn; i++) {
			nanoflann::KNNResultSet<double> resultSet(num_results);
			resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
			kd_tree.index->findNeighbors(resultSet, &X_scale(i),
						     nanoflann::SearchParams(1));
		}
		
		//find transformation
		asopa(X_curr, Y_close, asopa_threshold, Q_mod, A_mod, t_mod, RMSE_asopa);
		Q = Q * Q_mod;
		A = A * A_mod;
		
		//calculate error
		RMSE = sqrt((Y_close - (Q * (A * X_curr))).squaredNorm()/Xn);
		if(abs(RMSE_prev - RMSE) < threshold) {
			break;
		}
		
		//set prev to current
		RMSE_prev = RMSE;
	}
	
	//calculate transform using centroids with rotation and scaling
	t =  Y_centroid - (Q * (A * X_centroid));

	//calculate final error
	RMSE = sqrt((Y - ((Q * (A * X)).colwise() + t)).squaredNorm()/Xn);

	return 0;
}

