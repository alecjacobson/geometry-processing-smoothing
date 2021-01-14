#include "cotmatrix.h"
#include <igl/doublearea.h>
#include <iostream>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{

	std::cout << "Constructing laplacian..." << std::flush;

	// Analytical formulae are easily derived for cot(alpha) and cot(beta) in terms of
	// edge lengths only (see attached scan "laplacian_derivation_Shashwat.jpg"), and these are implemented below.

	L.resize(F.maxCoeff()+1, F.maxCoeff()+1);

	double sum_L = 0.0; // To sum off-diagonal entries to enter into the diagonal entries later on
	std::vector<double> sum_vec (F.maxCoeff()+1);

	// Compute areas*2.0
	Eigen::VectorXd dblA (F.rows());
	igl::doublearea(l, 0.0, dblA);

	// Use triplets to store values into the sparse matrix 
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplet_list;

	// Loop through each face
	for (int ii = 0; ii < F.rows(); ii++)
	{
		std::vector<int> V_ii (3); // Vertices that belong to this face
		V_ii[0] = F(ii, 0);
		V_ii[1] = F(ii, 1);
		V_ii[2] = F(ii, 2);

		// Sort the vertices to make life easier later on
		std::sort(V_ii.begin(), V_ii.end());

		// For each face, find the other faces that share edges with it.
		// To make things a bit faster, use the knowledge that there can only
		// be at most 3 faces touching this face.
		int count = 0;
		for (int jj = 0; jj < F.rows(); jj++)
		{
			if (count == 3)
				break;

			if (ii == jj) // Same triangle; handled separately
			{
				continue;
			}
			else
			{
				std::vector<int> V_jj (3); // Vertices that belong to this face
				V_jj[0] = F(jj, 0);
				V_jj[1] = F(jj, 1);
				V_jj[2] = F(jj, 2);

				// Sort the vertices to make life easier later on
				std::sort(V_jj.begin(), V_jj.end());

				// Now if any two vertices are the same, then those two triangles share an edge
				std::vector<int> common_vertices;
				std::set_intersection(V_ii.begin(), V_ii.end(), V_jj.begin(), V_jj.end(), std::back_inserter(common_vertices));

				if (common_vertices.size() == 2) // Face ii and face jj share an edge
				{
					count++;

					// Find the common edge and the uncommon edges of each face
					int common_edge1, common_edge2; // The index of the common edge
					std::vector<int> other_edges_ii, other_edges_jj;

					for (int kk = 0; kk < 3; kk++)
					{
						if (F(ii, kk) == common_vertices[0] || F(ii, kk) == common_vertices[1])
							other_edges_ii.push_back(kk);
						else
							common_edge1 = kk;

						if (F(jj, kk) == common_vertices[0] || F(jj, kk) == common_vertices[1])
							other_edges_jj.push_back(kk);
						else
							common_edge2 = kk;
					}

					double a1 = l(ii, common_edge1);
					double a2 = l(jj, common_edge2);

					// std::cout << common_edge1 << ", " << other_edges_ii[0] << ", " << other_edges_ii[1] << std::endl;
					double b = l(ii, other_edges_ii[0]);
					double c = l(ii, other_edges_ii[1]);
					double d = l(jj, other_edges_jj[0]);
					double e = l(jj, other_edges_jj[1]);

					double cos_alpha = (-a1*a1 + b*b + c*c);
					double cos_beta = (-a2*a2 + e*e + d*d);

					double sin_alpha = 2.0*dblA(ii);
					double sin_beta = 2.0*dblA(jj);

					// Note: a factor of (2*b*c) is common to both the cos and the
					// sin terms, so would cancel out. So, to save some computation,
					// it is neglected in anticipation of it cancelling out.

					double cot_alpha = cos_alpha/sin_alpha;
					double cot_beta = cos_beta/sin_beta;

					int row = common_vertices[0];
					int col = common_vertices[1];

					double val = 0.5*(cot_alpha + cot_beta);
					triplet_list.push_back(T(row, col, val));
					//triplet_list.push_back(T(col, row, val));

					//sum_L += val;
					sum_vec[row] += 1.0*val;
				}
			}
		}
	}

	L.setFromTriplets(triplet_list.begin(), triplet_list.end());

	sum_L = L.sum();

	// Now insert diagonal entries
	for (int ii = 0; ii < F.maxCoeff()+1; ii++)
		L.insert(ii, ii) = -sum_vec[ii];

	std::cout << "done." << std::endl;

}




