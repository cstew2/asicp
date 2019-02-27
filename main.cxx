#include <iostream>
#include <Eigen/Dense>

#include "asopa.h"
#include "asicp.h"

void test_asopa(void)
{
	
	Eigen::MatrixXd X(4, 3);
	X << 0,0,1,
		1,0,0,
		0,1,0,
		0,0,0;
	       
		
	Eigen::MatrixXd Y(4, 3);
	Y << 3,0,0,
		1,0,0,
		2,0,0,
		0,0,0;

	Eigen::Matrix3d Q(3,3);
	Eigen::Matrix3d A(3,3);
	Eigen::Vector3d t(3);

	double threshold = 1e-9f;
	double FRE = 0.0f;
	Eigen::MatrixXd FRE_mag;
	
	std::cout << "X" << std::endl << X << std::endl << std::endl
		  << "Y" << std::endl << Y << std::endl;
	asopa(X, Y, threshold, Q, A, t, FRE, FRE_mag);

	std::cout << "Q:" << std::endl << Q << std::endl << std::endl;
	std::cout << "A:" << std::endl << A << std::endl << std::endl;
	std::cout << "t:" << std::endl << t << std::endl << std::endl;
	std::cout << "FRE: " << FRE << std::endl << std::endl;
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

