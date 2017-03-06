#include "cotmatrix.h"

void cotmatrix(
	const Eigen::MatrixXd & l,
	const Eigen::MatrixXi & F,
	Eigen::SparseMatrix<double> & L)
{
	// Add your code here
	auto num_verts = F.maxCoeff() + 1;


	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(F.rows() * 3 * 4);

	for (int i = 0; i < l.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			//compute cotangent for edge F(i, j), F(i, j+1 %3)
			//Have triangle edge lengths a, b, c, where a is the edge we're looking at.

			double b = l(i, j);
			double c = l(i, (j + 1) % 3);
			double a = l(i, (j + 2) % 3);

			double cos_alpha = (b*b + c*c - a*a) / (2 * b*c);
			double s = (a + b + c) / 2;
			double A = std::sqrt(s*(s - a)*(s - b)*(s - c));
			double sin_alpha = a * 2 * A / (a*b*c);
			double entry = cos_alpha / sin_alpha * 0.5;

			int fi = F(i, j);
			int fj = F(i, (j + 1) % 3);

			triplets.emplace_back(fi, fj, entry);
			triplets.emplace_back(fj, fi, entry);
			triplets.emplace_back(fi, fi, -entry);
			triplets.emplace_back(fj, fj, -entry);
		}
	}

	L.resize(num_verts, num_verts);
	L.setFromTriplets(triplets.begin(), triplets.end());
}

