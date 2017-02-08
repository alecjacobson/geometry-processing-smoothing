#include "cotmatrix.h"
#include <set>
#include <cmath>

using namespace Eigen;

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
	double a, b, c, cosine, halfCot;
	int v1, v2;

	int n = F.maxCoeff() + 1;
	int f = F.rows();
	int m = 3 * f / 2;	//for closed surface, each edge is repeated twice
	int numNan = 0;
	std::vector<Triplet<double>> lVal;
	lVal.reserve(9 * f);
	L.resize(n, n);

	for (int i = 0; i < f; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			//get squared lengths: c is the length of the edge being considered
			a = l(i, (j + 1) % 3);
			a = a*a;
			b = l(i, (j + 2) % 3);
			b = b*b;
			c = l(i, j);
			c = c*c;

			//2ab*cos(C)
			cosine = (a + b - c);

			// cot^2 t = cos^2 t / ( 1 - cos^2 t ), with premultiplication
			halfCot = 0.5*sqrt(cosine*cosine / (4*a*b - cosine*cosine));

			v1 = F(i, (j+1) % 3);
			v2 = F(i, (j+2) % 3);

			// only lower triangular matrix is needed for L
			if (v2 > v1)
				std::swap(v1, v2);

			lVal.push_back({ v1, v2, halfCot });
			lVal.push_back({ v1, v1, -halfCot });
			lVal.push_back({ v2, v2, -halfCot });
		}
	}

	L.setFromTriplets(lVal.begin(), lVal.end());

}

