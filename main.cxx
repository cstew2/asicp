#include <iostream>

#include "test.h"

int main(void)
{
	std::cout << "Testing ASOPA" << std::endl;
	test_asopa();
	
	std::cout << std::endl << std::endl;

	std::cout << "Testing ASICP" << std::endl;
	test_asicp();

	std::cout << std::endl << std::endl;

	std::cout << "Testing PCA" << std::endl;
	test_asicp();

	
	return 0;
}

