#include "cotmatrix.h"
#include <math.h>
#include <iostream>

using namespace std;

double cotangent(double opposite, 
    double neighbor_1, 
    double neighbor_2){

    double s = (neighbor_1 + neighbor_2 + opposite) / 2.0;
    double x = s - neighbor_1;
    double y = s - neighbor_2;
    double z = s - opposite;

    double numerator = (neighbor_1 * neighbor_1) + (neighbor_2 * neighbor_2) - (opposite * opposite);
    double denominator = 4 * sqrt(s * x * y * z);

    if (denominator >= - 0.00000001 and denominator <= 0.00000001) {
     cout << "Warning! Triangle quarter area very low: " << denominator << endl;
     if (denominator == 0.0){
      denominator = 0.000000001;
     }
    }
    return numerator / denominator;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  for (int faceIndex = 0; faceIndex < F.rows(); faceIndex++){
    auto vertices = F.row(faceIndex);
    // Indices
    int v1 = vertices[0];
    int v2 = vertices[1];
    int v3 = vertices[2];

    auto lengths = l.row(faceIndex);
    // Side lengths
    double s1 = lengths[0];
    double s2 = lengths[1];
    double s3 = lengths[2];
    
    // Side 1 and 2
    double cot12 = cotangent(s3, s1, s2) / 2.0;
    // by symmetry update both indices
    L.coeffRef(v1, v2) += cot12;
    L.coeffRef(v2, v1) += cot12;

    // Side 2 and 3
    double cot23 = cotangent(s1, s2, s3) / 2.0;
    L.coeffRef(v2, v3) += cot23;
    L.coeffRef(v3, v2) += cot23;

    // Side 3 and 1
    double cot31 = cotangent(s2, s3, s1) / 2.0;
    L.coeffRef(v1, v3) += cot31;
    L.coeffRef(v3, v1) += cot31;
  }

  // Set up Laplacian equality: what leaves a node, enters it
  for (int diagIndex = 0; diagIndex < L.rows(); diagIndex++){
    // for some reason cannot actually edit the .diagonal()
    L.coeffRef(diagIndex, diagIndex) = -1 * L.row(diagIndex).sum();
  }
}

