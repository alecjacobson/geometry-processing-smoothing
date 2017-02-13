#include "cotmatrix.h"

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L) //     for triangles, columns correspond to edges [1,2],[2,0],[0,1]
{
	auto cot_ij = [](auto a, auto b, auto c) 
	{ 
		double cos_ij = (a*a + b*b - c*c) / (2 * a*b);
		double sin_ij = sqrt(1 - cos_ij * cos_ij);
		return cos_ij / sin_ij;
	};

	int numVert = F.maxCoeff() + 1;

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(numVert);
	L.resize(numVert, numVert);

	for (int fIndex = 0; fIndex < F.rows(); ++fIndex)
	{
		int v0 = F(fIndex, 0);
		int v1 = F(fIndex, 1);
		int v2 = F(fIndex, 2);

		float cot_01 = 0.5*cot_ij(l(fIndex, 0), l(fIndex, 1), l(fIndex, 2));
		float cot_12 = 0.5*cot_ij(l(fIndex, 1), l(fIndex, 2), l(fIndex, 0));
		float cot_20 = 0.5*cot_ij(l(fIndex, 0), l(fIndex, 2), l(fIndex, 1));

		//off-diagonal
		tripletList.push_back(T(v0, v1, cot_01));
		tripletList.push_back(T(v1, v0, cot_01));

		tripletList.push_back(T(v1, v2, cot_12));
		tripletList.push_back(T(v2, v1, cot_12));

		tripletList.push_back(T(v2, v0, cot_20));
		tripletList.push_back(T(v0, v2, cot_20));

		//on diagonal
		tripletList.push_back(T(v0, v0, -cot_01));
		tripletList.push_back(T(v0, v0, -cot_20));

		tripletList.push_back(T(v1, v1, -cot_01));
		tripletList.push_back(T(v1, v1, -cot_12));

		tripletList.push_back(T(v2, v2, -cot_12));
		tripletList.push_back(T(v2, v2, -cot_20));
	}	

	L.setFromTriplets(tripletList.begin(), tripletList.end());
}

