#include "cotmatrix.h"
#include <cmath>

typedef Eigen::Triplet<double> T;

double cotc(double a, double b, double c, double area)
{
	// sinc = c / 2r and r = abc/4A
	double sin = c / ((a*b*c) / (2.0 * area));
	double cos = (a*a + b*b - c*c) / (2.0 * a*b);
	return 0.5 * (cos / sin);
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  int v = F.maxCoeff() + 1;
	L.resize(v, v);

	std::vector<T> tlist;
	tlist.reserve(v * 12);

	for (int i = 0; i < F.rows(); i++)
	{
		double a = l(i, 0);
		double b = l(i, 1);
		double c = l(i, 2);
		double s = (a + b + c) / 2.0;
		//by helon's formula 
		double area = sqrt(s * (s - a) * (s - b)* (s - c));

		tlist.push_back(T(F(i, 0), F(i, 1), cotc(a, b, c, area)));
		tlist.push_back(T(F(i, 1), F(i, 2), cotc(b, c, a, area)));
		tlist.push_back(T(F(i, 2), F(i, 0), cotc(c, a, b, area)));

		//reverse direction same value
		tlist.push_back(T(F(i, 1), F(i, 0), cotc(a, b, c, area)));
		tlist.push_back(T(F(i, 2), F(i, 1), cotc(b, c, a, area)));
		tlist.push_back(T(F(i, 0), F(i, 2), cotc(c, a, b, area)));

		//diagonal sum
		tlist.push_back(T(F(i, 0), F(i, 0), -cotc(a, b, c, area)));
		tlist.push_back(T(F(i, 0), F(i, 0), -cotc(c, a, b, area)));
		tlist.push_back(T(F(i, 1), F(i, 1), -cotc(a, b, c, area)));
		tlist.push_back(T(F(i, 1), F(i, 1), -cotc(b, c, a, area)));
		tlist.push_back(T(F(i, 2), F(i, 2), -cotc(b, c, a, area)));
		tlist.push_back(T(F(i, 2), F(i, 2), -cotc(c, a, b, area)));
	}
  
	L.setFromTriplets(tlist.begin(), tlist.end());
}

