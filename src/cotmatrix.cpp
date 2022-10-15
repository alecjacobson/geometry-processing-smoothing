#include "cotmatrix.h"

/**
* @brief Solves with the cosine law to produce the cosine of the angle
* opposite to the edge c.
*/
double cosineLaw(double a, double b, double c) {
	double top = a * a + b * b - c * c;
	double bottom = 2 * a * b;
	return top / bottom;
}

/**
* @brief Computes the cotagenent of the angle opposite to the edges with the
* given lengths in a triangle.
*
* @param lengths input array of 3 edge lengths.
* @param angles output array of 3 angles.
*/
void lengthsToCotangents(double* lengths, double* cots) {

	// First we apply cosine law.
	double cosines[3];
	cosines[0] = cosineLaw(lengths[1], lengths[2], lengths[0]);
	cosines[1] = cosineLaw(lengths[2], lengths[0], lengths[1]);
	cosines[2] = cosineLaw(lengths[0], lengths[1], lengths[2]);

	// Now we apply pythagorean and cotangent identities to produce the
	// cotangent.
	for (int i = 0; i < 3; i++) {

		// Pythaogrean identity is valid because the angle is known to be 
		// less than PI; so sine wouldn't have gone negative in there.
		double cosine = cosines[i];
		double sine = sqrt(1 - cosine * cosine);
		double cotangent = cosine / sine;
		cots[i] = cotangent;
	}
}

void cotmatrix(
		const Eigen::MatrixXd & l,
		const Eigen::MatrixXi & F,
		Eigen::SparseMatrix<double> & L) {

	int vertexCount = F.maxCoeff() + 1;
	int faceCount = l.rows();

	L.resize(vertexCount, vertexCount);
	for (int i = 0; i < faceCount; i++) {

		// Find cotangents opposite to each edge.
		double lengths[] = { l(i, 0), l(i, 1), l(i, 2) };
		double cots[3];
		lengthsToCotangents(lengths, cots);

		// For each edge, fill in the appropriate matrix stuff.
		L.coeffRef(F(i, 0), F(i, 1)) += 0.5 * cots[2];
		L.coeffRef(F(i, 1), F(i, 2)) += 0.5 * cots[0];
		L.coeffRef(F(i, 2), F(i, 0)) += 0.5 * cots[1];

		// And its symmetric after all...
		L.coeffRef(F(i, 1), F(i, 0)) += 0.5 * cots[2];
		L.coeffRef(F(i, 2), F(i, 1)) += 0.5 * cots[0];
		L.coeffRef(F(i, 0), F(i, 2)) += 0.5 * cots[1];
	}

	// Now set the diagonal
	for (int i = 0; i < vertexCount; i++) {
		L.insert(i, i) = -L.row(i).sum();
	}


}

