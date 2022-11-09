#include "cotmatrix.h"
#include <igl/doublearea.h>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here

	// compute area for all facets
	Eigen::VectorXd doublearea(F.rows());
	igl::doublearea(l, doublearea);

	// compute 1/2*cot for angle inside each facet
	typedef Eigen::Triplet<double> T;
	typedef Eigen::SparseMatrix<double> SpMat;
	std::vector<T> triplets;
	triplets.reserve(F.rows() * 4);
	for (int i = 0; i < F.rows(); i++) {
		double a_2 = l(i, 2)*l(i, 2);
		double b_2 = l(i, 0)*l(i, 0);
		double c_2 = l(i, 1)*l(i, 1);
		double douArea_inv = 1 / doublearea(i);
		double cot_a = 0.5*(b_2 + c_2 - a_2) * douArea_inv;
		double cot_b = 0.5*(a_2 + c_2 - b_2) * douArea_inv;
		double cot_c = 0.5*(a_2 + b_2 - c_2) * douArea_inv;
		triplets.push_back(T(F(i, 0), F(i, 1), 0.5*cot_a));
		triplets.push_back(T(F(i, 1), F(i, 2), 0.5*cot_b));
		triplets.push_back(T(F(i, 2), F(i, 0), 0.5*cot_c));
	}
	int V_count = F.maxCoeff() + 1;
	SpMat cotMat(V_count, V_count);
	cotMat.setFromTriplets(triplets.begin(), triplets.end());

	// compute 1/2*cot alpha_ij + 1/2*cot beta_ij for each edge ij that exists
	SpMat cotMat_t = cotMat.transpose();
	L = cotMat_t + cotMat;

	// compute diagonal part
	for (int i = 0; i < V_count; i++) {
		L.insert(i, i) = -L.row(i).sum();
	}
}

