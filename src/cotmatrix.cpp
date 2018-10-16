#include "cotmatrix.h"
#include <igl/doublearea.h>
#include <cmath>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
    L.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);

    Eigen::VectorXd area;
    igl::doublearea(l, area);

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(F.rows() * F.cols());

    for (int i = 0; i < F.rows(); ++i) {
        const Eigen::RowVector3d edge = l.row(i);

        const double cot0 = (edge(1) * edge(1) + edge(2) * edge(2) - edge(0) * edge(0)) / (8 * area(i));
        const double cot1 = (edge(0) * edge(0) + edge(2) * edge(2) - edge(1) * edge(1)) / (8 * area(i));
        const double cot2 = (edge(1) * edge(1) + edge(0) * edge(0) - edge(2) * edge(2)) / (8 * area(i));

        tripletList.push_back(Eigen::Triplet<double>(F(i, 0), F(i, 1), cot2));
        tripletList.push_back(Eigen::Triplet<double>(F(i, 1), F(i, 2), cot0));
        tripletList.push_back(Eigen::Triplet<double>(F(i, 2), F(i, 0), cot1));
        tripletList.push_back(Eigen::Triplet<double>(F(i, 1), F(i, 0), cot2));
        tripletList.push_back(Eigen::Triplet<double>(F(i, 2), F(i, 1), cot0));
        tripletList.push_back(Eigen::Triplet<double>(F(i, 0), F(i, 2), cot1));
    }

    L.setFromTriplets(tripletList.begin(), tripletList.end());

    for (int i = 0; i < L.rows(); ++i) {
        // See issue #15
        L.coeffRef(i, i) -= L.row(i).sum();
    }
}

