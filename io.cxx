#include "io.h"

#include <iostream>
#include <fstream>

void read_ply_points(char *filename)
{
	std::ofstream fp;
	fp.open(filename, std::ios::in);

	if(fp.is_open()) {
		
	}
}

void write_points(char *filename, char *data)
{
	
}
