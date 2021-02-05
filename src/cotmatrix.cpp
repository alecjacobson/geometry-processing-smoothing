#include "cotmatrix.h"
typedef Eigen::Triplet<double> T;

void setV2F(Eigen::MatrixXi & V2F1, Eigen::MatrixXi & V2F2, int i, int j,
		int k) {

	if (V2F1(i, j) == -1) {
		V2F1(i, j) = k;
		V2F1(j, i) = k;
	} else {
		V2F2(i, j) = k;
		V2F2(j, i) = k;
	}
	return;
}

void cotmatrix(const Eigen::MatrixXd & l, const Eigen::MatrixXi & F,
		Eigen::SparseMatrix<double> & L) {
	// Add your code here

	int vnum = F.maxCoeff() + 1;
	int fnum = F.rows();

	// calculate area of faces
	Eigen::VectorXd FA(fnum);
	for (int i = 0; i < fnum; i++) {
		int v1 = F(i, 0);
		int v2 = F(i, 1);
		int v3 = F(i, 2);
		double e1len = l(v1, v2);
		double e2len = l(v1, v3);
		double e3len = l(v2, v3);
		double s = e1len + e2len + e3len;
		s = s / 2;
		double A2 = s * (s - e1len) * (s - e2len) * (s - e3len);
		double A = sqrt(A2);
		FA[i] = A;
	}

	// build point2face connection
	// 2 point i, j will link to 2 faces
	Eigen::MatrixXi V2F1, V2F2;
	V2F1 = -Eigen::MatrixXi::Ones(vnum, vnum);
	V2F2 = -Eigen::MatrixXi::Ones(vnum, vnum);
	for (int i = 0; i < fnum; i++) {
		int v1 = F(i, 0);
		int v2 = F(i, 1);
		int v3 = F(i, 2);
		setV2F(V2F1, V2F2, v1, v2, i);
	}

	std::vector<T> coef;
	for (int i = 0; i < vnum; i++) {
		// diag
		double ii = 0;
		for (int j = 0; j < vnum; j++) {

			if (i == j)
				continue;

			// is i and j connected?
			int k1 = V2F1(i, j);
			int k2 = V2F2(i, j);
			// cot c = (a^2 + b^2 - c^2) / 4A

			int k = k1;
			if (k != -1) {
				int v1 = F(k, 0);
				int v2 = F(k, 1);
				int v3 = F(k, 2);
				double e1len = l(v1, v2);
				double e2len = l(v1, v3);
				double e3len = l(v2, v3);
				double A = FA[k];
				double cot = e1len * e1len + e2len * e2len - e3len * e3len;
				cot = cot / (4 * A);
				T tmp(i, j, cot);
				coef.push_back(tmp);
				ii += cot;
			}

			k = k2;
			if (k != -1) {
				int v1 = F(k, 0);
				int v2 = F(k, 1);
				int v3 = F(k, 2);
				double e1len = l(v1, v2);
				double e2len = l(v1, v3);
				double e3len = l(v2, v3);
				double A = FA[k];
				double cot = e1len * e1len + e2len * e2len - e3len * e3len;
				cot = cot / (4 * A);
				T tmp(i, j, cot);
				coef.push_back(tmp);
				ii += cot;
			}
		}
		T tmp(i, i, ii);
		coef.push_back(tmp);
	}

	L.setFromTriplets(coef.begin(), coef.end());
}

