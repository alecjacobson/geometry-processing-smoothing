#include "cotmatrix.h"
#include "math.h"
#include <Eigen/Sparse>
#include <unordered_map>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
	int m = F.rows();
	double cot = 0.0;
	double sum = 0;
	double h = 0;
	double sin = 0;
	double cos = 0;
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(12*m);
	// Could just use cos(x)^2 + sin(x)^2 = 1 to compute cot(x) from cos(x). 
	// This way has less sqrt()
	for (int t = 0; t < m; ++t) {
		int i = F(t,0);
		int j = F(t,1);
		int k = F(t,2);
		double ei  = l(t,0);
		double ej  = l(t,1);
		double eij = l(t,2);

		cos = (eij*eij - ej*ej - ei*ei)/ (-2*ej*ei);
		sin = sqrt(1-cos*cos);
		cot = cos/sin;
		tripletList.push_back(Eigen::Triplet<double>(i, j, 0.5*cot));
		tripletList.push_back(Eigen::Triplet<double>(j, i, 0.5*cot));
		tripletList.push_back(Eigen::Triplet<double>(i, i, -0.5*cot));
		tripletList.push_back(Eigen::Triplet<double>(j, j, -0.5*cot));

		cos = (ei*ei - eij*eij - ej*ej) / (-2*eij*ej);
		sin = sqrt(1 - cos*cos);
		cot = cos / sin;
		tripletList.push_back(Eigen::Triplet<double>(j, k, 0.5*cot));
		tripletList.push_back(Eigen::Triplet<double>(k, j, 0.5*cot));
		tripletList.push_back(Eigen::Triplet<double>(j, j, -0.5*cot));
		tripletList.push_back(Eigen::Triplet<double>(k, k, -0.5*cot));

		cos = (ej*ej - ei*ei - eij*eij) / (-2*ei*eij);
		sin = sqrt(1 - cos*cos);
		cot = cos / sin;
		tripletList.push_back(Eigen::Triplet<double>(i, k, 0.5*cot));
		tripletList.push_back(Eigen::Triplet<double>(k, i, 0.5*cot));
		tripletList.push_back(Eigen::Triplet<double>(i, i, -0.5*cot));
		tripletList.push_back(Eigen::Triplet<double>(k, k, -0.5*cot));
		
	}

	L.setFromTriplets(tripletList.begin(),tripletList.end()); // sums up triplets with same (i,j) keys
}

