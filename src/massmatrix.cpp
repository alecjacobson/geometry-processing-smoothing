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
    
    double area, invA;
    for (int i = 0; i < F.rows(); i ++) {
        invA = 0.25 * sqrt(pow(1.0/pow(l(i,0),2) + 1.0/pow(l(i,1),2) + 1.0/pow(l(i,2),2), 2) - 2*(pow(l(i,0),4) + pow(l(i,1),4) + pow(l(i,2),4)));
        area = 1.0 / invA;
        for (int j = 0; j < 3; j ++) {
            tripletList.push_back(F(i,j), F(i,j), area / 3.0);
        }
    }
    //Might cause some issues
    M.setFromTriplets(tripletList.begin(), tripletList.end());
}

