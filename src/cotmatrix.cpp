#include "cotmatrix.h"

double cot3(double l1, double l2, double l3, double area){
  // sin3 = 2 * area / l1*l2  
  // cos3 = l1^2 + l2^2 - l3^2 / 2*l1*l2
  // cot3 = l1^2 + l2^2 - l3^2 / 4 * area

  double c3 = (l1*l1 + l2*l2 - l3*l3)/(4.0*area);
  return c3;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  int num_v = F.maxCoeff() + 1;
  L.resize(num_v, num_v);
  // L is V x V
  
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;

  Eigen::VectorXd areas(F.rows()), diag(num_v);
  igl::doublearea(l, 0., areas);
  areas = areas/2.0;
  // Interesting thing to figure out,
  // passing areas(i)/2.0 directly in the function call
  // did not work, but dividing by 2 here worked

  for (int i=0; i<F.rows(); i++){
    for (int j=0; j<3; j++){
      // For triangle meshes
      double cot = cot3(l(i,(j+1)%3), l(i,(j+2)%3), l(i,j), areas(i));
      cot = cot/2.0; // since we need to add 0.5*cot
 
      int idx1 = F(i, (j+1)%3);
      int idx2 = F(i, (j+2)%3);
      assert(idx1 != idx2);

      // twice for symmetric matrix
      tripletList.push_back(T(idx1, idx2, cot));
      tripletList.push_back(T(idx2, idx1, cot));
    
      // diagonals
      tripletList.push_back(T(idx1, idx1, -cot));
      tripletList.push_back(T(idx2, idx2, -cot));
    }
  }

  L.setFromTriplets(tripletList.begin(), tripletList.end());
  // setFromTriplets adds if the same indices are repeated
}
