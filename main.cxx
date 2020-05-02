#include <iostream>
#include <string.h>

#include "main.hxx"
#include "test.hxx"

int run_tests(void)
{
	/*
	std::cout << "Testing ASOPA" << std::endl;
	test_asopa();
	
	std::cout << std::endl << std::endl;

	*/
	std::cout << "Testing ASICP" << std::endl;
	test_asicp();

	std::cout << std::endl << std::endl;
	
	return  0;
}

int print_help(void)
{
	std::cout << prog_name << " usage: asicp [OPTIONS]..."
		  << std::endl << std::endl
		  << "Options" << std::endl
		  << "-f \t args: \"source file\" \"destination file\""    << std::endl
		  << "-h \t prints help" << std::endl
		  << "-t \t runs test" << std::endl;
	return 0;
}

int main(int argc, char **argv)
{
	if(argc == 1) {
		print_help();	
	}
	
	std::string filename_src;
	std::string filename_dst;

	for(int i=0; i < argc; i++) {
		if(!strcmp(argv[i], "-h")) {
			print_help();
		}
		else if(!strcmp(argv[i], "-t")) {
			run_tests();	
		}
		else if(!strcmp(argv[i], "-f")) {
			filename_src = argv[++i];
			filename_dst = argv[++i];
		}
	}

	
	
	return 0;
}

