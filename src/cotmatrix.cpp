#include "cotmatrix.h"
#include <iostream>
#include <math.h>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
	L.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);
	L.setZero();
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(6 * F.rows());

	for (int f = 0; f < F.rows(); f++) {
		int i = F(f, 0); int j = F(f, 1); int k = F(f, 2);
		double a = l(f, 0); double b = l(f, 1); double c = l(f, 2);
		double s = (a + b + c) / 2.0;
		double A = sqrt(s * (s-a) * (s-b) * (s-c));
		double cot;
		// Consider edge ij
		cot = (pow(a, 2) + pow(b, 2) - pow(c, 2)) / (4 * A);
		triplets.push_back({ i, j, 0.5 * cot });
		triplets.push_back({ j, i, 0.5 * cot });

		// Consider edge jk
		cot = (pow(b, 2) + pow(c, 2) - pow(a, 2)) / (4 * A);
		triplets.push_back({ j, k, 0.5 * cot });
		triplets.push_back({ k, j, 0.5 * cot });

		// Consider edge ik
		cot = (pow(a, 2) + pow(c, 2) - pow(b, 2)) / (4 * A);
		triplets.push_back({ i, k, 0.5 * cot });
		triplets.push_back({ k, i, 0.5 * cot });
	}
	L.setFromTriplets(triplets.begin(), triplets.end());

	// Fill in the diagonal entries
	Eigen::SparseMatrix<double> Ld(F.maxCoeff() + 1, F.maxCoeff() + 1);
	Ld.setZero();
	std::vector<Eigen::Triplet<double>> tripletsd;
	tripletsd.reserve(L.rows());
	for (int f = 0; f < L.rows(); f++) {
		tripletsd.push_back({ f, f,  L.row(f).sum() });
	}

	Ld.setFromTriplets(tripletsd.begin(), tripletsd.end());

	L = L - Ld;
}

