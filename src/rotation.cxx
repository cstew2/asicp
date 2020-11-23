#include <iostream>
#include <cmath>

#include <Eigen/Dense>

#include "rotation.hxx"

std::vector<Eigen::Quaterniond> get_rots(int n)
{
	std::vector<Eigen::Quaterniond> q;

	std::vector<Eigen::Quaterniond> q_i;
	std::vector<Eigen::Quaterniond> q_j;
	std::vector<Eigen::Quaterniond> q_k;
	
	q.push_back(Eigen::Quaternion<double>(0, 0, 0, 0));
	q.push_back(Eigen::Quaternion<double>(0, 1, 0, 0));
	q.push_back(Eigen::Quaternion<double>(0, 0, 1, 0));
	q.push_back(Eigen::Quaternion<double>(0, 0, 0, 1));
	
	//generate for each axis
	for(double i=1; i < n; i++) {	
		double a = cos(0.5*(M_PI*(i/n)));
		double r = sqrt(1-a*a);
		q_i.push_back(Eigen::Quaternion<double>( r,  a,  0,  0));
		q_i.push_back(Eigen::Quaternion<double>( r, -a,  0,  0));
		q_i.push_back(Eigen::Quaternion<double>(-r, -a,  0,  0));
		q_i.push_back(Eigen::Quaternion<double>(-r,  a,  0,  0));
		
		q_j.push_back(Eigen::Quaternion<double>( r,  0,  a,  0));
		q_j.push_back(Eigen::Quaternion<double>( r,  0, -a,  0));
		q_j.push_back(Eigen::Quaternion<double>(-r,  0, -a,  0));
		q_j.push_back(Eigen::Quaternion<double>(-r,  0,  a,  0));
		
		q_k.push_back(Eigen::Quaternion<double>( r,  0,  0,  a));
		q_k.push_back(Eigen::Quaternion<double>( r,  0,  0, -a));
		q_k.push_back(Eigen::Quaternion<double>(-r,  0,  0, -a));
		q_k.push_back(Eigen::Quaternion<double>(-r,  0,  0,  a));
	}

	for(int i=1; i < n; i++) {
		for(int j=1; j < n; j++) {
			q.push_back(q_i[i]);
			q.push_back(q_j[i]);
			q.push_back(q_k[i]);
				
			q.push_back(q_i[i] * q_j[j]);
			q.push_back(q_i[i] * q_k[j]);
			
			q.push_back(q_j[i] * q_i[j]);
			q.push_back(q_j[i] * q_k[j]);

			q.push_back(q_k[i] * q_i[j]);
			q.push_back(q_k[i] * q_j[j]);
		}
	}

	//should think of a better way to do this
	for(int i=0; i < q.size(); i++) {
		for(int j=0; j < q.size(); j++) {
			if(q[i].coeffs() == q[j].coeffs() && i!=j) {
			//std::cout << q[i].coeffs() << std::endl << q[j].coeffs() << std::endl << std::endl;
				q.erase(q.begin() + j);
			}
		}
	}
	
	return q;
}
