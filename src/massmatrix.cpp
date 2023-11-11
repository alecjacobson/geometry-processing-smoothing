#include "massmatrix.h"

void massmatrix(const Eigen::MatrixXd & l, const Eigen::MatrixXi & F,
		Eigen::DiagonalMatrix<double, Eigen::Dynamic> & M) {
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

	Eigen::MatrixXi V2F = Eigen::MatrixXi::Zero(vnum, fnum);
	for (int i = 0; i < fnum; i++) {
		int v1 = F(i, 0);
		int v2 = F(i, 1);
		int v3 = F(i, 2);
		V2F(v1, i) = 1;
		V2F(v2, i) = 1;
		V2F(v3, i) = 1;
	}
	Eigen::VectorXd diag(vnum);
	for (int i = 0; i < vnum; i++) {
		double ii = 0;
		for (int j = 0; j < fnum; j++) {
			if (V2F(i, j) == 1) {
				ii += FA[j];
			}
		}

		diag[i] = ii / 3;
	}

	M = diag.asDiagonal();
}

