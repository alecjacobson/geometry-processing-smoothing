#include "cotmatrix.h"

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
	// Add your code here

	//Notes: Cool thing about triplets, when inserted into a sparse matrix 
	//the elements that have the same index will be automatically summed.
	//This means we can just deal on a triangle-by-triangle basis, calculating
	//a single cotagent and when we are done processing the final cotangent
	//weight will be available in the correct index. 

	//Also note on igl::edge_lengths() it returns the length of the edges in 
	//the following order [1,2],[2,0],[0,1]
	// -> https://github.com/libigl/libigl/blob/master/include/igl/edge_lengths.h

	//Allocate space to store the cotangent laplacian
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(3 * 4 * F.rows()); // 3 vertices, 4 permutations, m triangles

	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < 3; j++) {

			double a = l(i, j);
			double b = l(i, (j + 1) % 3);
			double c = l(i, (j + 2) % 3);
			
			//Combine law-of-cosines with sin^2 + cos^2 = 1 and cot = cos/sin
			double cosVal = ((a*a) + (b*b) - (c*c)) / (2 * a * b);
			double cotVal = 0.5 * (cosVal / std::sqrt(1 - cosVal*cosVal));

			//need this weird offset cause of the way igl::edge_lengths returns the lengths
			int vi = F(i, (j + 1) % 3);
			int vj = F(i, (j + 2) % 3);

			triplets.emplace_back(vi, vj, cotVal);
			triplets.emplace_back(vj, vi, cotVal);
			triplets.emplace_back(vj, vj, -cotVal);
			triplets.emplace_back(vi, vi, -cotVal);
		}
	}

	//Assign the sparse cotangent laplacian matrix
	L.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);
	L.setFromTriplets(triplets.begin(), triplets.end());	

}

