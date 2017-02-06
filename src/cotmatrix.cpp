#include "cotmatrix.h"
#include <vector>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> trips;

    //cot(A) = cos(A) / sin(A)
    //
    //sin(A) = a/d
    //d = abc/2*area
    //area = sqrt(s(s-a)(s-b)(s-c))
    //s = (a+b+c)/2
    //cos(A) = (b^2+c^2-a^2)/(2bc)
    
    int mV = F.maxCoeff()+1;
    std::vector<double> diag(mV,0);
    for(int fi = 0; fi < F.rows(); ++fi) {
        auto&& le = l.row(fi).array();
        auto&& f = F.row(fi).array();

        double s = (le.sum()) / 2.0;
        double area = std::sqrt(s * (s-le).prod());
        double d = le.prod() / (2 * area);

        auto Sa = le/d;

        for(int i = 0; i < 3; ++i) {
            double a = le((i+0)%3);
            double b = le((i+1)%3);
            double c = le((i+2)%3);
            double ca = ( b*b + c*c - a*a ) / ( 2 * b*c );
            double v = .5 * ca / Sa(i);


            //int A = f((i+0)%3);
            int B = f((i+1)%3);
            int C = f((i+2)%3);
            trips.emplace_back(B,C,v);
            trips.emplace_back(C,B,v);
            diag[B] += -v;
            diag[C] += -v;
        }
        
    }
    for(int i = 0; i < diag.size(); ++i) {
        trips.emplace_back(i,i,diag[i]);
    }

    L = Eigen::SparseMatrix<double>(mV,mV);
    L.setFromTriplets(trips.begin(),trips.end());
}

