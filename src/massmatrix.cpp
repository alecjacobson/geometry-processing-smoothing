#include "massmatrix.h"

double Herons_formula(double a, double b, double c){
    double s = (a+b+c)/2.0;
    return 0.25*std::sqrt(s*(s-a)*(s-b)*(s-c));
}

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
    int n = F.maxCoeff() + 1;
    M.resize(n);
    for(int i = 0; i < n; i++){
        double area = 0;
        for(int f = 0; f < F.rows(); f++){
            if(F(f,0) == i || F(f,1) == i || F(f,2) == i){
                area += Herons_formula(l(f,0),l(f,1),l(f,2));
            }
        }
        M.diagonal()[i] = area/3.0;
    }
}

