#include "cotmatrix.h"
#include <math.h>

double cot(double a, double b, double c) {
	double s = (a + b + c) / 2.0;
	double A = sqrt(s * (s - a) * (s - b) * (s - c));
	return (b*b + c*c - a*a) / (4.0 * A);
}


void cotmatrix(
	const Eigen::MatrixXd & l,
	const Eigen::MatrixXi & F,
	Eigen::SparseMatrix<double> & L)
{
	int n = F.maxCoeff()+1;
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(9 * F.rows());
	
	for (int i = 0; i < F.rows(); i++) {
		double a = l(i, 0);
		double b = l(i, 1);
		double c = l(i, 2);
		// case i == j
		tripletList.push_back(T(F(i, 0), F(i, 0), -cot(b, c, a)/2.0));
		tripletList.push_back(T(F(i, 0), F(i, 0), -cot(c, a, b)/2.0));
		tripletList.push_back(T(F(i, 1), F(i, 1), -cot(a, b, c)/2.0));
		tripletList.push_back(T(F(i, 1), F(i, 1), -cot(c, a, b)/2.0));
		tripletList.push_back(T(F(i, 2), F(i, 2), -cot(a, b, c)/2.0));
		tripletList.push_back(T(F(i, 2), F(i, 2), -cot(b, c, a)/2.0));
		// case i != j
		tripletList.push_back(T(F(i, 0), F(i, 2), cot(b, c, a)/2.0));
		tripletList.push_back(T(F(i, 0), F(i, 1), cot(c, a, b)/2.0));
		tripletList.push_back(T(F(i, 1), F(i, 2), cot(a, b, c)/2.0));
		tripletList.push_back(T(F(i, 1), F(i, 0), cot(c, a, b)/2.0));
		tripletList.push_back(T(F(i, 2), F(i, 1), cot(a, b, c)/2.0));
		tripletList.push_back(T(F(i, 2), F(i, 0), cot(b, c, a)/2.0));
	}

	L.resize(n, n);
	L.setFromTriplets(tripletList.begin(), tripletList.end());
}

