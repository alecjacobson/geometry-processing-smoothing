#include "massmatrix.h"

void massmatrix(
    const Eigen::MatrixXd & l,
    const Eigen::MatrixXi & F,
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
    Eigen::VectorXd dblA(F.rows());
    igl::doublearea(l, 0, dblA);

    Eigen::VectorXd areas = Eigen::VectorXd::Zero(F.maxCoeff() + 1);
    
    for (int i = 0; i < F.rows(); i++) 
    {
        Eigen::Vector3i f = F.row(i);
        double area = dblA[i];

        for (int j = 0; j < f.size(); j++) 
        {
            int v = f[j];
            areas[v] += area;

        }
    }

    M.resize(F.maxCoeff() + 1);
    M.diagonal() = 1 / (double)6 * areas;
}

