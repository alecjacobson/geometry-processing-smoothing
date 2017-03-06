#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{

    int mV = F.maxCoeff()+1;
    M.resize(mV);
    M.setZero();

    for(int fi = 0; fi < F.rows(); ++fi) {
        auto&& le = l.row(fi).array();
        auto&& f = F.row(fi).array();

        double s = (le.sum())/2.0;
        double area = std::sqrt(s * (s-le).prod());
        for(int i = 0; i < 3; ++i) {
            M.diagonal()(f(i)) += area;
        }

    }
    M.diagonal().noalias() = M.diagonal() / 3.0;
}

