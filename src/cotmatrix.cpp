#include "cotmatrix.h"

typedef Eigen::Triplet<double> tri;

// finds the cotangent of the angle opposite edge a
double cotangent(double a, double b, double c)
{
  double cosA = (a*a + b*b + c*c)/(2*c*b);
  double sinA = sqrt(1 - cosA*cosA);
  return cosA/sinA;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Given edges l, and the Tets F, construct the "cotangent Laplacian":
  /*
           / ½cotαᵢⱼ + ½cotβᵢⱼ  if edge ij exists
    Lᵢⱼ =  | ∑ᵢ≠ⱼ Lᵢⱼ           if i = j
           \ 0                  otherwise

    the angles are coming from the edge lengths in the l matrix

    angle A is opposite side a
    using the cosine law:
    cos(A) = (a² + b² + c²)/2bc
    then we know that sin(A) = sqrt(1 - cos²(A))
    then the cotangent is given by cos(A)/sin(A)

    L is of size #V x #V and is sparse (fill with triplets in eigen)
    l is of size  3 x #F, where l(row, i) corresponds to the edge length opposite F(row, i)
   */
  std::vector<tri> tripletList;

  for( int32_t i = 0; i < F.rows(); i++ )
  {
    int32_t v0 = F(i, 0);
    int32_t v1 = F(i, 1);
    int32_t v2 = F(i, 2);

    double cot0 = cotangent(l(i, 0), l(i, 1), l(i, 2)); // cotangent opposite edge l₀
    double cot1 = cotangent(l(i, 1), l(i, 2), l(i, 0));
    double cot2 = cotangent(l(i, 2), l(i, 0), l(i, 1));

    // angle at v0 is in cot0.
    // for off diagonals:
    // L(vi, vj) = ½cotk
    // the component from the angle opposite will be added in a later iteration of loop (½cot(β))
    tripletList.push_back(tri(v0, v1, 0.5*cot2));
    tripletList.push_back(tri(v1, v0, 0.5*cot2));
    tripletList.push_back(tri(v1, v2, 0.5*cot0));
    tripletList.push_back(tri(v2, v1, 0.5*cot0));
    tripletList.push_back(tri(v2, v0, 0.5*cot1));
    tripletList.push_back(tri(v0, v2, 0.5*cot1));
    // diagonals:
    // L(vi, vi) = - (cotj + cotk)
    tripletList.push_back(tri(v0, v0, -cot1-cot2));
    tripletList.push_back(tri(v1, v1, -cot0-cot2));
    tripletList.push_back(tri(v2, v2, -cot1-cot0));
  }
  int32_t max = F.maxCoeff() + 1;
  L.resize(max, max);
  L.setFromTriplets(tripletList.begin(), tripletList.end());
}

