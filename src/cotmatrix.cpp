#include "cotmatrix.h"
#include <vector>
#include <cmath> // only for sqrt(.) and pow(.)

// computes the cot of angle_ij given adj side lengths ki, kj, and the opposite side ij to angle_ij
double compute_cot_from_lengths(double kj, double ki, double ij);

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here

  std::vector<Eigen::Triplet<double> > L_vector;

  int num_vertices = F.maxCoeff() + 1;

  // duplicates get summed up when setting from triplets
  // so, no need to interate on all possible edges

  // iterate on the face matrix F to get all the existing edges
  for (int f = 0; f < F.rows(); f++) {


    // Vertex - index: i - 0; j - 1; k - 2

    double len_ij = l(f, 2);
    double len_ki = l(f, 1);
    double len_kj = l(f, 0);

    // handle the case where the edge e_ij exists
    double cot_ij = compute_cot_from_lengths(len_kj, len_ki, len_ij);
    double cot_kj = compute_cot_from_lengths(len_ki, len_ij, len_kj);
    double cot_ki = compute_cot_from_lengths(len_ij, len_kj, len_ki);

    L_vector.push_back(Eigen::Triplet<double>(F(f, 0), F(f, 1), 0.5*cot_ij));
    L_vector.push_back(Eigen::Triplet<double>(F(f, 1), F(f, 0), 0.5*cot_ij));

    L_vector.push_back(Eigen::Triplet<double>(F(f, 1), F(f, 2), 0.5*cot_kj));
    L_vector.push_back(Eigen::Triplet<double>(F(f, 2), F(f, 1), 0.5*cot_kj));

    L_vector.push_back(Eigen::Triplet<double>(F(f, 0), F(f, 2), 0.5*cot_ki));
    L_vector.push_back(Eigen::Triplet<double>(F(f, 2), F(f, 0), 0.5*cot_ki));

    // handle the case where i=j and edge exists 
    // need negative of sum of all L_ik (for k != i) when edge exists between k and i
    // exploting the property that duplicate entry gets added finally
    L_vector.push_back(Eigen::Triplet<double>(F(f, 0), F(f, 0), -1.0*(0.5*cot_ij + 0.5*cot_ki)));
    L_vector.push_back(Eigen::Triplet<double>(F(f, 1), F(f, 1), -1.0*(0.5*cot_ij + 0.5*cot_kj)));
    L_vector.push_back(Eigen::Triplet<double>(F(f, 2), F(f, 2), -1.0*(0.5*cot_kj + 0.5*cot_ki)));

    // rest of the cases, it's zero. so nothing to do.
    
  }

  L.resize(num_vertices, num_vertices);
  L.setFromTriplets(L_vector.begin(), L_vector.end());

}

// computes the cot of angle_ij given adj side lengths ki, kj, and the opposite side ij to angle_ij
double compute_cot_from_lengths(double kj, double ki, double ij) {

  // first compute the area using Heron's formula
  double semi_perimeter = (kj + ki + ij)/2.0;
  double area_sq = semi_perimeter * (semi_perimeter - kj) * (semi_perimeter - ki) * (semi_perimeter - ij) * 1.0;
  double area = sqrt(area_sq);

  // use area = 0.5 ab Sin(C) to find Sine angles
  double sine_ij = (2.0*area)/(kj*ki*1.0);

  // use cosine rule to compute cos(angle_ij)
  double cosine_ij = (kj*kj + ki*ki - ij*ij)/(2.0*kj*ki);

  // cot_ij
  double cot_ij = cosine_ij/sine_ij;

  return cot_ij;
}

