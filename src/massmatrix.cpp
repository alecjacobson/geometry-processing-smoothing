#include "massmatrix.h"

void massmatrix(
		const Eigen::MatrixXd & l,
		const Eigen::MatrixXi & F,
		Eigen::DiagonalMatrix<double, Eigen::Dynamic> & M) {

	int vertexCount = F.maxCoeff() + 1;
	int triangleCount = F.rows();

	M.resize(vertexCount);
	for (int i = 0; i < vertexCount; i++) {

		double totalArea = 0;
		for (int t = 0; t < triangleCount; t++) {

			// Check if the vertex i is contained in the triangle t.
			if (F(t, 0) != i && F(t, 1) != i && F(t, 2) != i) {
				continue;
			}

			// Compute area of the triangle from the side lengths.
			// Use Heron's formula.
			double s = (l(i, 0) + l(i, 1) + l(i, 2)) / 2;
			double area = sqrt(s * (s - l(i, 0)) * (s - l(i, 1)) * 
					(s - l(i, 2)));
			totalArea += area;
		}
		M.diagonal()[i] = totalArea;
	}

	M = M * (1 / 3.0);

}

