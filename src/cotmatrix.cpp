#include "cotmatrix.h"
#include <iostream>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  int num_faces = F.rows();
  int V = F.maxCoeff() + 1;
  L.resize(V, V);
  L.setZero();

  Eigen::VectorXd sum_off_diagonals = Eigen::VectorXd::Zero(V);
  std::vector< Eigen::Triplet<double> > tripletList;
  for(int f_idx = 0; f_idx < num_faces; f_idx++){
    int i = F(f_idx, 0); int j = F(f_idx, 1); int k = F(f_idx, 2);
    // EDGE IJ
    double half_cot;
    half_cot = 0.5 * cotangent_triangle(i,j,f_idx, l, F);
    tripletList.push_back(Eigen::Triplet<double>(i,j, half_cot));
    tripletList.push_back(Eigen::Triplet<double>(j,i, half_cot));
    sum_off_diagonals(i) += half_cot; sum_off_diagonals(j) += half_cot;
    // EDGE JK
    half_cot = 0.5 * cotangent_triangle(j,k,f_idx, l, F);
    tripletList.push_back(Eigen::Triplet<double>(j,k, half_cot));
    tripletList.push_back(Eigen::Triplet<double>(k,j, half_cot));
    sum_off_diagonals(j) += half_cot; sum_off_diagonals(k) += half_cot;
    // EDGE KI
    half_cot = 0.5 * cotangent_triangle(k,i,f_idx, l, F);
    tripletList.push_back(Eigen::Triplet<double>(k,i, half_cot));
    tripletList.push_back(Eigen::Triplet<double>(i,k, half_cot));
    sum_off_diagonals(k) += half_cot; sum_off_diagonals(i) += half_cot;
  }
  // Diagonal entries
  for(int v = 0; v < V; v++){
    tripletList.push_back(Eigen::Triplet<double>(v,v, -sum_off_diagonals(v)));
  }

   L.setFromTriplets(tripletList.begin(), tripletList.end());
}

double cotangent_triangle(
  const int i,
  const int j,
  const int f_idx,
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F)
  {
    // find i's index in F -- this is a lil ugly
    int i_pos; int j_pos; int k_pos;
    for(int x = 0; x < 3; x++){
      if(F(f_idx, x) == i){
        i_pos = x;
        j_pos = (x + 1) % 3;
        k_pos = (x + 2) % 3;
      }
    }
    // Find edge length of edge ij
    double c = l(f_idx, k_pos);
    double a = l(f_idx, i_pos);
    double b = l(f_idx, j_pos);

    // Herons formula + solve for diameter of triangle's circumference
    double sp = (a + b + c) / 2; // semiperimeter
    double area = sqrt(sp * (sp - a) * (sp - b) * (sp - c));
    double d = (a * b * c) / (2 * area);

    // From sin and cos laws
    double sin_theta = c / d;
    double cos_theta = ((a*a) + (b*b) - (c * c)) / (2 * a * b);
    return cos_theta / sin_theta;
  }
