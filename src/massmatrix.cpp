#include "massmatrix.h"
#include <iostream>

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
    int numV = F.maxCoeff() + 1;
    
    Eigen::VectorXd curValues;
    curValues.resize(numV);
    curValues.setZero();
    
    double area, s;
    for (int i = 0; i < F.rows(); i ++) {
        //Compute area using Heron's formula
        s = l.row(i).sum() / 2.0;
        area = sqrt(s * (s - l(i,0)) * (s - l(i,1)) * (s - l(i,2)));
        for (int j = 0; j < 3; j ++) {
            curValues(F(i,j)) += area / 3.0;
        }
    }
    
    //Update elements in M
    M.diagonal() = curValues;

}

