#include <iostream>
#include <Eigen/Dense>

#include "asopa.h"
#include "asicp.h"

void test_asopa(void)
{
	size_t n = 20;

	srand(time(NULL));
	
	Eigen::Matrix3d R = (Eigen::AngleAxisd(M_PI/3.19, Eigen::Vector3d::UnitX()).matrix() *
			     Eigen::AngleAxisd(3.99*M_PI/2.8, Eigen::Vector3d::UnitY()).matrix() *
			     Eigen::AngleAxisd(5.12*M_PI/0.4, Eigen::Vector3d::UnitZ()).matrix());
	Eigen::Matrix3d S = Eigen::Vector3d::Random(3).asDiagonal() * 20.2;
	Eigen::Vector3d l = Eigen::RowVector3d::Random(3) * 10.7;
	
	Eigen::MatrixXd X = Eigen::MatrixXd::Random(3, n) * 102.4;
	Eigen::MatrixXd Y = (R * (S * X)).colwise() +  l;

	Eigen::Matrix3d Q(3,3);
	Eigen::Matrix3d A(3,3);
	Eigen::Vector3d t(3);

	double threshold = 1e-9f;
	double FRE = 0.0f;

	if(!asopa(X, Y, threshold, Q, A, t, FRE)) {
		if(n < 10) {
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
		std::cout << "FRE: " << FRE << std::endl << std::endl;

		std::cout << "Error statistics" << std::endl;
		std::cout << "Q-R Error: " << sqrt((Q*A-R*S).squaredNorm()/3*3) << std::endl;
		std::cout << "t-l Error: " << sqrt((t-l).squaredNorm()/3) << std::endl;
		
		std::cout << "Q-R Error: " << fabs(sqrt((Q*A).squaredNorm()/3*3) - sqrt((R*S).squaredNorm()/3*3)) << std::endl;
		std::cout << "t-l Error: " << fabs(sqrt((t).squaredNorm()/3) - sqrt((l).squaredNorm()/3)) << std::endl;
		
		std::cout << "Actual Error of Y - (QAX+t)" << std::endl;
		std::cout << "FRE:" <<
			sqrt((Y - ((Q * (A * X)).colwise() + t)).squaredNorm()/n)
			  << std::endl;
	}
	
}

void test_asicp(void)
{
	size_t n = 200;
	
	srand(time(NULL));
	
	Eigen::Matrix3d R = (Eigen::AngleAxisd(M_PI/3.19, Eigen::Vector3d::UnitX()).matrix() *
			     Eigen::AngleAxisd(3.99*M_PI/2.8, Eigen::Vector3d::UnitY()).matrix() *
			     Eigen::AngleAxisd(5.12*M_PI/0.4, Eigen::Vector3d::UnitZ()).matrix());
	
	Eigen::Matrix3d S = Eigen::Vector3d::Random(3).asDiagonal() * 10.2;
	Eigen::Vector3d l = Eigen::Vector3d::Random(3) * 10.7;
	
	Eigen::MatrixXd X = Eigen::MatrixXd::Random(3, n) * 102.4;
	Eigen::MatrixXd Y = (R * (S * X)).colwise() +  l;
	Eigen::MatrixXd Y_mod(3, Y.cols()/2);

	for(size_t i=0; i < n/2; i++) {
		Y_mod.col(i) = Y.col(i*2);
	}

	Eigen::Matrix3d Q(3,3);
	Eigen::Matrix3d A(3,3);
	Eigen::Vector3d t(3);

	double threshold = 1e-9;
	double FRE = 0.0f;
	size_t maximum_iterations = 20;

	if(!asicp(X, Y_mod, threshold, maximum_iterations, Q, A, t, FRE)) {
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
		std::cout << "FRE: " << FRE << std::endl << std::endl;

		std::cout << "Error statistics" << std::endl;
		std::cout << "Q-R Error: " << sqrt((Q*A-R*S).squaredNorm()/3*3) << std::endl;
		std::cout << "t-l Error: " << sqrt((t-l).squaredNorm()/3) << std::endl;
		
		std::cout << "Q-R Error: " << fabs(sqrt((Q*A).squaredNorm()/3*3) - sqrt((R*S).squaredNorm()/3*3)) << std::endl;
		std::cout << "t-l Error: " << fabs(sqrt((t).squaredNorm()/3) - sqrt((l).squaredNorm()/3)) << std::endl;
		
		std::cout << "Actual Error of Y - (QAX+t)" << std::endl;
		std::cout << "FRE:" <<
			sqrt((Y - ((Q * (A * X)).colwise() + t)).squaredNorm()/n)
			  << std::endl;
	}
}

int main(void)
{
	//std::cout << "Testing ASOPA" << std::endl;
	//test_asopa();
	
	std::cout << std::endl << std::endl;

	std::cout << "Testing ASICP" << std::endl;
	test_asicp();

	return 0;
}

