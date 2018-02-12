#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
    int numV = F.maxCoeff() + 1;
    
    M.resize(numV,numV);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(F.rows()*3);
    
    double area, s;
    for (int i = 0; i < F.rows(); i ++) {
        s = F.row(i).sum() / 2.0;
        area = sqrt(s * (s - l(i,0)) * (s - l(i,1)) * (s - l(i,2)));
        for (int j = 0; j < 3; j ++) {
            tripletList.push_back(F(i,j), F(i,j), area / 3.0);
        }
    }
    //Might cause some issues
    M.setFromTriplets(tripletList.begin(), tripletList.end());
}

