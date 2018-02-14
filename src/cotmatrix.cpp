#include "cotmatrix.h"
#include <cmath>
#include <vector>
#include <iostream>
//calculate cos C
double cos_abc(double a, double b, double c){
    return (a * a + b * b - c * c) / (2 * a * b);
}

//calculate sin C
double sin_abc(double a, double b, double c){
    double s = (a + b + c)/2;
    double area = sqrt(s * (s - a) * (s -b) * (s - c));
    //return c / (a*b*c/(2 * area));
    return 2 * area/(a * b);
}


void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
    L.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);
    for (int row = 0; row < F.rows(); ++row){
        for (int col = 0; col < F.cols(); ++col){ //3 edges [0, 1], [1, 2] [2, 0]
            int i = F(row, col);
            int j = F(row, (col + 1) % F.cols());
            double a = l(row, col);
            double b = l(row, (col +1) % l.cols());
            double c = l(row, (col +2) % l.cols());
            L.coeffRef(i,j) += 0.5 * cos_abc(a, b, c)/sin_abc(a, b, c);
            //std::cout <<i << " " << j << " " <<0.5 * cos_abc(a, b, c)/sin_abc(a, b, c) << std::endl;
        } 
        for (int col = 0; col < F.cols(); ++col){ //3 edges [1, 0], [2, 1] [0, 2] symmetric
            int j = F(row, col);
            int i = F(row, (col + 1) % F.cols());
            L.coeffRef(i,j) = L.coeffRef(j,i);
            //std::cout <<i << " " << j << " " <<0.5 * cos_abc(a, b, c)/sin_abc(a, b, c) << std::endl;
        }
        for (int col = 0; col < F.cols(); ++col){ //diag
            int i = F(row, col);
            double a = l(row, col); //[1,2]
            double b = l(row, (col +1) % l.cols()); //[0,2]
            double c = l(row, (col +2) % l.cols()); //[0,1]
            L.coeffRef(i,i) -= 0.5 * cos_abc(a, b, c)/sin_abc(a, b, c);
            
            L.coeffRef(i,i) -= 0.5 * cos_abc(a, c, b)/sin_abc(a, c, b);

        }
    }
}

