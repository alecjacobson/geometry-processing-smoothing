#include "cotmatrix.h"
#include <cmath>

typedef Eigen::Triplet<double> T;

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
    std::vector<Eigen::Triplet<double>> trips;
    L.resize(F.maxCoeff()+1, F.maxCoeff()+1);
    
    for (int i = 0; i < F.rows(); i++) {
        for (int j_ = 0; j_ < 3; j_++) {
            int j[] = {j_, (j_+1)%3, (j_+2)%3};
            double a = l(i,j[2]);
            double b = l(i,j[0]);
            double c = l(i,j[1]);
            
            //Heron's formula to calculate area
            double s = (a+b+c) / 2.0;
            double A = sqrt(s * (s-a) * (s-b) * (s-c));
            //half the cotangent = cos / sin / 2
            double cot = (((a*a-b*b-c*c)/(-2*b*c)) / ((2*A)/(b*c))) / 2;
            
            //put to matrix
            trips.push_back(T(F(i,j[0]), F(i,j[1]), cot));
            trips.push_back(T(F(i,j[1]), F(i,j[0]), cot));
            trips.push_back(T(F(i,j[0]), F(i,j[0]), -cot));
            trips.push_back(T(F(i,j[1]), F(i,j[1]), -cot));
        }
    }
    
    L.setFromTriplets(trips.begin(), trips.end());
}

