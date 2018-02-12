#include "cotmatrix.h"
#include <iostream>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(F.rows()*4);
    
    Eigen::MatrixXd vals(3,2);
    
    //Makes it easier to reference other vertices
    vals(0,0) = 1;
    vals(0,1) = 2;
    vals(1,0) = 0;
    vals(1,1) = 2;
    vals(2,0) = 0;
    vals(2,1) = 1;
    

    double cosVal, sinVal, cotVal;
    
    int numV = F.maxCoeff() + 1;
    
    L.resize(numV,numV);
    
    for (int i = 0; i < F.rows(); i ++) {
        for (int j = 0; j < 3; j ++) {
            //Cosine Law
            cosVal = 1.0/ (2.0 * l(i,vals(j,0)) * l(i,vals(j,1)));
            cosVal = cosVal * (pow(l(i,vals(j,0)),2.0) + pow(l(i,vals(j,1)),2.0) - pow(l(i,j),2.0));
            
            //We use the fact that sine is always positive.
            sinVal = sqrt(1.0 - pow(cosVal,2.0));
            
            cotVal = cosVal / sinVal;
            tripletList.push_back(T(F(i,vals(j,0)), F(i,vals(j,1)), cotVal / 2.0));
            tripletList.push_back(T(F(i,vals(j,1)), F(i,vals(j,0)), cotVal / 2.0));
            
            tripletList.push_back(T(F(i,vals(j,0)), F(i,vals(j,0)), -cotVal / 2.0));
            tripletList.push_back(T(F(i,vals(j,1)), F(i,vals(j,1)), -cotVal / 2.0));
            
        }
        
        
    }
    L.setFromTriplets(tripletList.begin(), tripletList.end());
}

