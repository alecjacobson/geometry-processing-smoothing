#include "cotmatrix.h"
#include <cmath>
#include <vector>
#include <iostream>

void cotmatrix(
  const Eigen::MatrixXd & l, // F by 3, lengths corresponding to each edge defined in F
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  double e_i, e_j, e_ij, s, area, cosine, sine, halfCot;
  int i, j;
  std::vector<Eigen::Triplet<double>> triplets;

  for (int faceInd = 0; faceInd < l.rows(); faceInd++) {
    for (int edgeInd = 0; edgeInd < 3; edgeInd++) { // 0 gives length of [1,2], 1 gives length of [2,0], 2 gives length of [0,1]
      
      e_ij = l(faceInd, edgeInd);
      e_i = l(faceInd, (edgeInd + 1) % 3);
      e_j = l(faceInd, (edgeInd + 2) % 3);
      
      // half-cotangent calculation
      s = (e_i + e_j + e_ij) / 2.0;
      area = std::sqrt(s * (s - e_i) * (s - e_j) * (s - e_ij));
      cosine = (std::pow(e_i, 2) + std::pow(e_j, 2) - std::pow(e_ij, 2)) / (2.0 * e_i * e_j);
      sine = (2.0 * area) / (e_i * e_j);
      halfCot = cosine / (2.0 * sine);
      
      // indices of vertices, to plug result into L matrix
      i = F(faceInd, (edgeInd + 1) % 3);
      j = F(faceInd, (edgeInd + 2) % 3);
      
      // the triplet will ADD to this entry of L! so it'll sum the cotangents for us w/ each half edge, and take care of the i = j case too
      triplets.push_back(Eigen::Triplet<double>(i, j, halfCot)); 
      triplets.push_back(Eigen::Triplet<double>(j, i, halfCot)); 
      triplets.push_back(Eigen::Triplet<double>(i, i, - halfCot));
      triplets.push_back(Eigen::Triplet<double>(j, j, - halfCot));
    }
  }
  L.setFromTriplets(triplets.begin(), triplets.end());
}

