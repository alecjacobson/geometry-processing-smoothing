#include "cotmatrix.h"

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(F.rows()*4);
    
    Eigen::MatrixXd vals(3,2);
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
            cosVal = l(i,vals(j,0)) * l(i,vals(j,1)) / 2.0;
            cosVal = cosVal * (1.0/l(i,vals(j,0)) + 1.0 / l(i,vals(j,1)) - 1.0/l(i,j));
            
            //We use the fact that sine is always positive.
            sinVal = sqrt(1.0 - pow(cosVal,2.0));
            
            cotVal = sinVal / cosVal;
            
            tripletList.push_back(F(i,vals(j,0)), F(i,vals(j,1)), cotVal / 2.0);
            tripletList.push_back(F(i,vals(j,1)), F(i,vals(j,0)), cotVal / 2.0);
            
            tripletList.push_back(F(i,vals(j,0)), F(i,vals(j,0)), -cotVal / 2.0);
            tripletList.push_back(F(i,vals(j,1)), F(i,vals(j,1)), -cotVal / 2.0);
            
        }
        
        
    }
    L.setFromTriplets(tripletList.begin(), tripletList.end());
}

