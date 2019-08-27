#include <iostream>
#include <Eigen/Dense>

#include "asopa.hxx"
#include "asicp.hxx"
#include "pca.hxx"

void test_asopa(void)
{
	size_t n = 29;

	srand(time(NULL));
	
	Eigen::Matrix3d R = (Eigen::AngleAxisd(M_PI/3.19, Eigen::Vector3d::UnitX()).matrix() *
			     Eigen::AngleAxisd(3.99*M_PI/2.8, Eigen::Vector3d::UnitY()).matrix() *
			     Eigen::AngleAxisd(5.12*M_PI/0.4, Eigen::Vector3d::UnitZ()).matrix());
	Eigen::Matrix3d S = Eigen::Vector3d::Random(3).asDiagonal() * 22.2;
	Eigen::Vector3d l = Eigen::RowVector3d::Random(3) * 10.7;
	
	Eigen::MatrixXd X = Eigen::MatrixXd::Random(3, n) * 102.4;
	Eigen::MatrixXd Y = (R * (S * X)).colwise() +  l;

	Eigen::Matrix3d Q(3,3);
	Eigen::Matrix3d A(3,3);
	Eigen::Vector3d t(3);

	double threshold = 1e-9f;
	double RMSE = 0.0f;

	if(!asopa(X, Y, threshold, Q, A, t, RMSE)) {
		if(n < 30) {
			std::cout << "Initial configurations:" << std::endl;
			std::cout << "X" << std::endl << X << std::endl << std::endl
				  << "Y" << std::endl << Y << std::endl << std::endl;
			
			std::cout << "Q*A*X+t:" << std::endl
				  << (Q * (A * X)).colwise() + t
				  << std::endl << std::endl;
		}
	
		std::cout << "Actual Transformation from X -> Y" << std::endl;
		std::cout << "R:" << std::endl << R << std::endl << std::endl;
		std::cout << "S:" << std::endl << S << std::endl << std::endl;
		std::cout << "l:" << std::endl << l << std::endl << std::endl;	

		std::cout << "Calculated Transformation from X -> Y" << std::endl;
		std::cout << "Q:" << std::endl << Q << std::endl << std::endl;
		std::cout << "A:" << std::endl << A << std::endl << std::endl;
		std::cout << "t:" << std::endl << t << std::endl << std::endl;
		std::cout << "RMSE: " << RMSE << std::endl << std::endl;

		std::cout << "Error statistics" << std::endl;
		std::cout << "QA-RS Error: " << sqrt((Q*A-R*S).squaredNorm()/3*3) << std::endl;
		std::cout << "t-l Error: " << sqrt((t-l).squaredNorm()/3) << std::endl;

		std::cout << "Actual Error of Y - (QAX+t)" << std::endl;
		std::cout << "RMSE:" <<
			sqrt((Y - ((Q * (A * X)).colwise() + t)).squaredNorm()/n)
			  << std::endl;
	}
	
}

void test_asicp(void)
{
	size_t n = 29;
	
	srand(time(NULL));
	
	Eigen::Matrix3d R = (Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitX()).matrix() *
			     Eigen::AngleAxisd(M_PI/6, Eigen::Vector3d::UnitY()).matrix() *
			     Eigen::AngleAxisd(M_PI/3, Eigen::Vector3d::UnitZ()).matrix());
	Eigen::Matrix3d S = Eigen::Vector3d::Random(3).asDiagonal() * 1.2;
	Eigen::Vector3d l = Eigen::Vector3d::Random(3) * 2.2;
	
	Eigen::MatrixXd X = Eigen::MatrixXd::Random(3, n) * 10.4;
	Eigen::MatrixXd Y = (R * (S * X)).colwise() +  l;
	Eigen::MatrixXd Y_mod = Y;
	
	Eigen::Matrix3d Q;
	Eigen::Matrix3d A;
	Eigen::Vector3d t;

	Q = Eigen::Matrix3d::Identity();
	A = Eigen::Matrix3d::Identity();
	t = Eigen::Vector3d::Zero();
	
	double threshold = 1e-9;
	double asopa_threshold = 1e-32;
	double RMSE = 0.0f;
	size_t maximum_iterations = 1000;

	if(!asicp(X, Y, threshold, maximum_iterations, asopa_threshold, true, Q, A, t, RMSE)) {
		if(n < 30) {
			std::cout << "Initial configurations:" << std::endl;
			std::cout << "X" << std::endl << X << std::endl << std::endl
				  << "Y" << std::endl << Y << std::endl << std::endl;
			
			std::cout << "Q*A*X+t:" << std::endl
				  << (Q * (A * X)).colwise() + t
				  << std::endl << std::endl;
		}
	
		std::cout << "Actual Transformation from X -> Y" << std::endl;
		std::cout << "R:" << std::endl << R << std::endl << std::endl;
		std::cout << "S:" << std::endl << S << std::endl << std::endl;
		std::cout << "l:" << std::endl << l << std::endl << std::endl;	

		std::cout << "Calculated Transformation from X -> Y" << std::endl;
		std::cout << "Q:" << std::endl << Q << std::endl << std::endl;
		std::cout << "A:" << std::endl << A << std::endl << std::endl;
		std::cout << "t:" << std::endl << t << std::endl << std::endl;
		std::cout << "RMSE: " << RMSE << std::endl << std::endl;

		std::cout << "Error statistics" << std::endl;
		std::cout << "Q-R Error: " << sqrt((Q*A-R*S).squaredNorm()/3*3) << std::endl;
		std::cout << "t-l Error: " << sqrt((t-l).squaredNorm()/3) << std::endl;
		
		std::cout << "Actual Error of Y - (QAX+t)" << std::endl;
		std::cout << "RMSE:" <<
			sqrt((Y - ((Q * (A * X)).colwise() + t)).squaredNorm()/n)
			  << std::endl;
	}
}

void test_pca(void)
{
	size_t n = 29;

	srand(time(NULL));
	
	Eigen::Matrix3d R = (Eigen::AngleAxisd(M_PI/3.19, Eigen::Vector3d::UnitX()).matrix() *
			     Eigen::AngleAxisd(3.99*M_PI/2.8, Eigen::Vector3d::UnitY()).matrix() *
			     Eigen::AngleAxisd(5.12*M_PI/0.4, Eigen::Vector3d::UnitZ()).matrix());
	/*
	R <<    1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0;
	*/	

	Eigen::Vector3d s = Eigen::RowVector3d::Random(3) * 10.7;
	s << 1.0, 1.0, 1.0;
	Eigen::Matrix3d S = s.asDiagonal();
	
	Eigen::Vector3d l = Eigen::RowVector3d::Random(3) * 10.7;
	
	Eigen::MatrixXd X = Eigen::MatrixXd::Random(3, n) * 102.4;
	Eigen::MatrixXd Y = (R * (S * X)).colwise() +  l;

	Eigen::Matrix3d Q(3,3);
	Eigen::Matrix3d A(3,3);
	Eigen::Vector3d t(3);

	double RMSE = 0.0f;

	if(!pca_registration(X, Y, Q, A, t, RMSE)) {
		if(n < 30) {
			std::cout << "Initial configurations:" << std::endl;
			std::cout << "X" << std::endl << X << std::endl << std::endl
				  << "Y" << std::endl << Y << std::endl << std::endl;
			
			std::cout << "Q*A*X+t:" << std::endl
				  << (Q * (A * X)).colwise() + t
				  << std::endl << std::endl;
		}
	
		std::cout << "Actual Transformation from X -> Y" << std::endl;
		std::cout << "R:" << std::endl << R << std::endl << std::endl;
		std::cout << "S:" << std::endl << S << std::endl << std::endl;
		std::cout << "l:" << std::endl << l << std::endl << std::endl;	

		std::cout << "Calculated Transformation from X -> Y" << std::endl;
		std::cout << "Q:" << std::endl << Q << std::endl << std::endl;
		std::cout << "A:" << std::endl << A << std::endl << std::endl;
		std::cout << "t:" << std::endl << t << std::endl << std::endl;
		std::cout << "RMSE: " << RMSE << std::endl << std::endl;

		std::cout << "Error statistics" << std::endl;
		std::cout << "QA-RS Error: " << sqrt((Q*A-R*S).squaredNorm()/3*3) << std::endl;
		std::cout << "t-l Error: " << sqrt((t-l).squaredNorm()/3) << std::endl;

		std::cout << "Actual Error of Y - (QAX+t)" << std::endl;
		std::cout << "RMSE:" <<
			sqrt((Y - ((Q * (A * X)).colwise() + t)).squaredNorm()/n)
			  << std::endl;
	}
}
