#ifndef COTMATRIX_H
#define COTMATRIX_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cmath>
// Construct the "cotangent Laplacian" for a mesh with edge lengths `l`. Each
// entry in the output sparse, symmetric matrix `L` is given by:
//
// \\[
// L_{ij} = \begin{cases}
//          ½ \cot{α_{ij}} + ½ \cot{β_{ij}}  & \text{if edge $ij$ exists} \\\\\
//         -∑_{j≠i} L_{ij}                   & \text{if $i = j} \\\\\
//          0                                & \text{otherwise}
//          \end{cases}
// \\]
//
// where $α_{ij}$ and $β_{ij}$ are the angles opposite the edge between
// vertices i and j.
//
// Inputs:
//   l  #F by 3 matrix so that l(f,c) is the length of the "half-edge"
//     **_across_** from the cth of the fth face. 
//   F  #F by 3 list of triangle indices into a vertex list `V`. You may assume
//     that the number of vertices #V = F.maxCoeff()+1
// Outputs:
//   L  #V by #V sparse symmetric Laplacian matrix
//
void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L);
#endif
