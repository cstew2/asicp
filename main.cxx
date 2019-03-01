#include <iostream>
#include <Eigen/Dense>

#include "asopa.h"
#include "asicp.h"

void test_asopa(void)
{
	size_t n = 4;

	srand(time(NULL));
	
	Eigen::Matrix3d R = (Eigen::AngleAxisd(M_PI/3.19, Eigen::Vector3d::UnitX()).matrix() *
			     Eigen::AngleAxisd(3.99*M_PI/2.8, Eigen::Vector3d::UnitY()).matrix() *
			     Eigen::AngleAxisd(5.12*M_PI/0.4, Eigen::Vector3d::UnitZ()).matrix());
	Eigen::Matrix3d S = Eigen::Vector3d::Random(3).asDiagonal() * 20.2;
	Eigen::RowVector3d l = Eigen::RowVector3d::Random(3) * 10.7;
	
	Eigen::MatrixXd X = Eigen::MatrixXd::Random(n, 3) * 102.4;
	Eigen::MatrixXd Y = (R * S * X.transpose()).transpose().rowwise() +  l;

	Eigen::Matrix3d Q(3,3);
	Eigen::Matrix3d A(3,3);
	Eigen::RowVector3d t(3);

	double threshold = 1e-32f;
	double FRE = 0.0f;
	Eigen::MatrixXd FRE_mag;

	if(!asopa(X, Y, threshold, Q, A, t, FRE, FRE_mag)) {
		if(n < 50) {
			std::cout << "Initial configurations:" << std::endl;
			std::cout << "X" << std::endl << X << std::endl << std::endl
				  << "Y" << std::endl << Y << std::endl << std::endl;
			
			std::cout << "Q*A*X+t:" << std::endl
				  << (Q * A * X.transpose()).transpose().rowwise() + t
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
		std::cout << "FRE: " << FRE << std::endl << std::endl;

		std::cout << "Actual Error of Y - (QAX+t)" << std::endl;
		std::cout << "FRE:" <<
			sqrt((Y - ((Q * A * X.transpose()).transpose().rowwise() + t)).squaredNorm()/n)
			  << std::endl;
	}

	
	
}

void test_asicp(void)
{
	Eigen::MatrixXd X(6, 3);
	X << 0,0,1,
		2,0,1,
		1,0,0,
		0,1,0,
		0,0,0,
		1,0,4;
		
	       
		
	Eigen::MatrixXd Y(4, 3);
	Y << 3,0,0,
		1,0,0,
		2,0,0,
		0,0,0;

	Eigen::Matrix3d Q(3,3);
	Eigen::Matrix3d A(3,3);
	Eigen::Vector3d t(3);

	double threshold = 1e-9;
	double FRE = 0.0f;
	size_t maximum_iterations = 2000;
	
	std::cout << "X" << std::endl << X << std::endl << std::endl
		  << "Y" << std::endl << Y << std::endl;
	asicp(X, Y, threshold, maximum_iterations, Q, A, t, FRE);

	std::cout << "Q:" << std::endl << Q << std::endl << std::endl;
	std::cout << "A:" << std::endl << A << std::endl << std::endl;
	std::cout << "t:" << std::endl << t << std::endl << std::endl;
	std::cout << "FRE: " << FRE << std::endl << std::endl;
}

int main(void)
{
	std::cout << "Testing ASOPA" << std::endl;
	test_asopa();

	/*
	std::cout << "Testing ASICP" << std::endl;
	test_asicp();
	*/
}

