#include "cotmatrix.h"
#include <iostream>
#include <vector>
#include <cmath>

typedef Eigen::Triplet<double> T;

// Given edge lengths of a triangle ABC
// Return the cotan of angle A
double my_cot(double a, double b, double c) {
	// Heron's Formula
	double s = (a + b + c) / 2.0;
	double area = sqrt(s * (s - a) * (s - b) * (s - c));

	// Law of Sine and Cosine
	double sinA = 2.0 * area / (b * c);
	double cosA = (a*a - b*b - c*c) / (-2.0 * b * c);

	return cosA / sinA;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
	int nV = F.maxCoeff() + 1;
	std::vector<T> tripletList;

	for(int fi = 0; fi < F.rows(); fi++) {
		for(int e = 0; e < 3; e++) {
			int i = F(fi, e); 
			int j = F(fi, (e + 1) % 3);

			double a = l(fi, (e + 2) % 3);
			double b = l(fi, (e + 0) % 3);
			double c = l(fi, (e + 1) % 3);

			double cotA = my_cot(a, b, c);

			double Lij = 0.5 * cotA;

			tripletList.push_back(T(i, j, Lij)); // Duplicates will be summed
			tripletList.push_back(T(j, i, Lij));

			// Diagonal
			tripletList.push_back(T(i, i, -Lij));
			tripletList.push_back(T(j, j, -Lij));
		}
	}

	L.resize(nV, nV);
	L.setFromTriplets(tripletList.begin(), tripletList.end());
}

