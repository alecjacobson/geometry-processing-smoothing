#ifndef MASSMATRIX_H
#define MASSMATRIX_H
#include <Eigen/Core>
#include <Eigen/Sparse>
// Construct the diagonal(ized) mass matrix for a mesh with edge lengths `l`.
// Each enetry in the output sparse, symmetric matrix `M` is given by:
//
// \\[
// M_{ij} = \begin{cases}
//          ∑ ⅓ A_{f}  & \text{if face $f$ containing $i$ and $j$ exists} \\\\\
//          0          & \text{otherwise}
//          \end{cases}
// \\]
//
// where $A_f$ is the area of the fth face.
//
// Inputs:
//   l  #F by 3 matrix so that l(f,c) is the length of the "half-edge"
//     **_across_** from the cth of the fth face. 
//   F  #F by 3 list of triangle indices into a vertex list `V`. You may assume
//     that the number of vertices #V = F.maxCoeff()+1
// Outputs:
//   M  #V by #V diagonal matrix
//
void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M);
#endif
