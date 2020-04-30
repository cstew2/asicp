#include <iostream>

#include <nanoflann.hpp>

#include "asicp.hxx"
#include "asopa.hxx"
#include "pca.hxx"
#include "rotation.hxx"

int asicp(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	  double threshold, size_t max_iterations, double asopa_threshold,
	  bool estimate, int rotations,
	  Eigen::Matrix3d &Q, Eigen::Matrix3d &A, Eigen::Vector3d &t,
	  double &RMSE)
{
	Eigen::Matrix3d best_Q;

	if(estimate) {
		asicp_rot(X, Y, threshold, max_iterations,
			  asopa_threshold, Q, A, t, RMSE);
		return 0;
	}

	Eigen::Matrix3d Q_best(3,3);
	Eigen::Matrix3d A_best(3,3);
	Eigen::Vector3d t_best(3);
	double RMSE_best = std::numeric_limits<double>::max();

	std::vector<Eigen::Quaterniond> rots = get_rots(rotations);

	//go through discrete subgroup of SO(3)
	for(int i=0; i < rots.size(); i++) {
		std::cout << "Iteration:" << i << std::endl;
		std::cout << rots[i].coeffs() << std::endl << std::endl;
		
		//get rotation
		Q = rots[i].toRotationMatrix();
				
		//std::cout << "R:" << std::endl << Q << std::endl;
		//std::cout << "A:" << std::endl << A << std::endl;
		
		//go through each scale
		for(int j=0; j < 3; j++) {
			
			//calculate scales with current rotation
			pca_scales(X, Y, Q, A);

			asicp_rot(X, Y, threshold, max_iterations, asopa_threshold, Q, A, t, RMSE);
			
			if(RMSE < RMSE_best) {
				Q_best = Q;
				A_best = A;
				t_best = t;
			}
			
			//change scaling order
			double temp = A(0,0);
			A(0,0) = A(1,1);
			A(1,1) = A(2,2);
			A(2,2) = temp;
		}
	}

	Q = Q_best;
	A = A_best;
	t = t_best;
	
	return 0;
}

int asicp_rot(Eigen::MatrixXd X, Eigen::MatrixXd Y,
	      double threshold, size_t max_iterations, double asopa_threshold,
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

	Eigen::Matrix3d A_orig = A;
	
	//initialise kd-tree
	typedef nanoflann::KDTreeEigenMatrixAdaptor<
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
		3,
		nanoflann::metric_L2,
		false>
		kd_tree_t;
	
	double RMSE_prev = sqrt((Y_trans - X_trans).squaredNorm()/Yn);
	double RMSE_asopa = 0.0;
	
	for(size_t j=0; j < max_iterations; j++) {
		//std::cout << "iteration " << j << " error delta: " << std::abs(RMSE_prev - RMSE) << std::endl;
		//std::cout << "Q" << std::endl << Q << std::endl;
		//std::cout << "A" << std::endl << A << std::endl;
		//std::cout << "t" << std::endl << t << std::endl;

		//bound scale values
		if(A(0, 0) < 0.9 * A_orig(0, 0)) {
			A(0, 0) = 1.1 * A_orig(0, 0);
		}
		if(A(0, 0) > 1.1 * A_orig(0, 0)) {
			A(0, 0) = 0.9 * A_orig(0, 0);
		}
		
		if(A(1, 1) < 0.9 * A_orig(1, 1)) {
			A(1, 1) = 1.1 * A_orig(1, 1);
		}
		if(A(1, 1) > 1.1 * A_orig(1, 1)) {
			A(1, 1) = 0.9 * A_orig(1, 1);
		}
		
		if(A(2, 2) < 0.9 * A_orig(2, 2)) {
			A(2, 2) = 1.1 * A_orig(2, 2);
		}
		if(A(2, 2) > 1.1 * A_orig(2, 2)) {
			A(2, 2) = 0.9 * A_orig(2, 2);
		}
		
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
			X_scale(0, i) = X_curr(0,i)/sqrt(A(0,0));
			X_scale(1, i) = X_curr(1,i)/sqrt(A(1,1));
			X_scale(2, i) = X_curr(2,i)/sqrt(A(2,2));
		}

		//initialize kd tree with scaled Y points
		kd_tree_t kd_tree(3, Y_scale, 10);
		kd_tree.index->buildIndex();
		for(size_t i=0; i < Xn; i++) {
			size_t ret_index;
			double out_dist_sqr;
			nanoflann::KNNResultSet<double> result_set(1);
			result_set.init(&ret_index, &out_dist_sqr);
			double query[3];
			query[0] = X_scale(0,i);
			query[1] = X_scale(1,i);
			query[2] = X_scale(2,i);
			//find closest points in the scaled Y points to the current transformed points of X
			kd_tree.index->findNeighbors(result_set,
						     &query[0],
						     nanoflann::SearchParams(10));
			Y_close.col(i) = X_curr.col(ret_index);
		}
       		
		//find transformation
		asopa(X_curr, Y_close, asopa_threshold, Q, A, t, RMSE_asopa);      		
		
		//set previous error to current
		RMSE_prev = RMSE;

		//calculate error
		RMSE = sqrt((Y_close - (Q * (A * X_curr))).squaredNorm()/Xn);
		if(std::abs(RMSE_prev - RMSE) < threshold) {
			break;
		}
	}
	
	//calculate transform using centroids with rotation and scaling
	t =  Y_centroid - (Q * (A * X_centroid));

	//calculate final error
	RMSE = sqrt((Y - ((Q * (A * X)).colwise() + t)).squaredNorm()/Xn);

	return 0;
}

