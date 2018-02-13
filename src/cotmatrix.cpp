#include "cotmatrix.h"

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
	L.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);

	float e_i, e_j, e_ij = 0;
	float cos_alpha, sin_alpha, cot_alpha = 0;

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(F.maxCoeff()+1);

	Eigen::VectorXf d(F.maxCoeff() + 1);

	auto get_cot = [&](int face, int a, int b, int c)
	{
		e_i = l(face, a);
		e_j = l(face, b);
		e_ij = l(face, c);
		cos_alpha = (pow(e_i, 2) + pow(e_j, 2) - pow(e_ij, 2)) / (2 * e_i*e_j);
		sin_alpha = sqrt(1 - pow(cos_alpha, 2));
		cot_alpha = cos_alpha / sin_alpha;
		return cot_alpha /= 2;
	};

	for (int i = 0; i < F.rows(); i++)
	{
		auto cot = get_cot(i, 0, 1, 2);
		tripletList.push_back(T(F(i,0), F(i,1), cot_alpha));
		auto cot = get_cot(i, 1, 2, 0);
		tripletList.push_back(T(F(i, 1), F(i, 2), cot_alpha));
		auto cot = get_cot(i, 2, 0, 1);
		tripletList.push_back(T(F(i, 2), F(i, 0), cot_alpha));

		d(F(i, 0)) = l(i, 0);
		d(F(i, 1)) = l(i, 1);
		d(F(i, 2)) = l(i, 2);
	}

	for (int i = 0; i < F.maxCoeff() + 1; i++)
	{
		tripletList.push_back(T(i,i, -d.sum()+d(i,i)));
	}

	L.setFromTriplets(tripletList.begin(), tripletList.end());
}

