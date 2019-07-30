#include "pca.h"

int pca_registration(Eigen::MatrixXd X, Eigen::MatrixXd Y,
		     Eigen::Matrix3d &Q,
		     Eigen::Matrix3d &A,
		     Eigen::Vector3d &t,
		     double &RMSE)
{
	size_t Xn = X.cols();
	size_t Yn = Y.cols();

	//remove noise/outliers?
	
	//translate input by centroids (data has a mean of 0)
	Eigen::MatrixXd X_trans = X.colwise() - X.rowwise().mean();
	Eigen::MatrixXd Y_trans = Y.colwise() - Y.rowwise().mean();
	
	//computer cross covariance matrix
	//C= X'Y/n-1 (use un-scaled for esimating scaling factors?)
	Eigen::Matrix3d Sx = 1/(Xn-1) * (X * X.transpose());
	Eigen::Matrix3d Sy = 1/(Yn-1) * (Y * Y.transpose());
	
	//solve for eigenvectors and eigenvalues
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
	es.compute(Sx);
	Eigen::Matrix3d X_evec = es.eigenvectors();
	Eigen::Vector3d X_eval = es.eigenvalues();
	
	es.compute(Sy);
	Eigen::Matrix3d Y_evec = es.eigenvectors();
	Eigen::Vector3d Y_eval = es.eigenvalues();

	//could try the all possible arrangements of PCs not just in order of largest
	
	Eigen::Matrix3d PC_0 = Eigen::Quaterniond().setFromTwoVectors
		(Y_evec.col(0), X_evec.col(0)).toRotationMatrix();

	Eigen::Matrix3d PC_1 = Eigen::Quaterniond().setFromTwoVectors
		(Y_evec.col(1), X_evec.col(1)).toRotationMatrix();
	
	Eigen::Matrix3d PC_2 = Eigen::Quaterniond().setFromTwoVectors
		(Y_evec.col(2), X_evec.col(2)).toRotationMatrix();

	Q = PC_0 * PC_1 * PC_2;

	Sx = (X * X.transpose());
	Sy = (Y * Y.transpose());

	es.compute(Sx);
	X_evec = es.eigenvectors();
	X_eval = es.eigenvalues();
	
	es.compute(Sy);
	Y_evec = es.eigenvectors();
	Y_eval = es.eigenvalues();
	
	A(0) = (Y_evec.col(0).squaredNorm())/(X_evec.col(0).squaredNorm());
	A(1) = (Y_evec.col(1).squaredNorm())/(X_evec.col(1).squaredNorm());
	A(2) = (Y_evec.col(2).squaredNorm())/(X_evec.col(2).squaredNorm());

	t =  Y.rowwise().mean() - (Q * (A * X.rowwise().mean()));

	RMSE = sqrt((Y - ((Q * (A * X)).colwise() + t)).squaredNorm()/Xn);
}
