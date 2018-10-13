#include "cotmatrix.h"

#include "igl/adjacency_list.h"

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{

    // Compute area of triangle for i-th face
    auto triangle_area = [&l](int i){
        double e0 = l(i, 0);
        double e1 = l(i, 1);
        double e2 = l(i, 2);
        double s = (e0 + e1 + e2) / 2;
        return sqrt(s*(s-e0)*(s-e1)*(s-e2));
    };

    int nv = F.maxCoeff() + 1;
    L.resize(nv, nv);
    Eigen::SparseMatrix<double> cots(nv, nv);
    std::vector<Eigen::Triplet<double>> triplets(F.rows()*3);

    int a, b, c;
    double ea, eb, ec;
    double A, cot;

    for (int i = 0; i < F.rows(); ++i) {
        A = triangle_area(i);
        for (int j = 0; j < 3; ++j) {
            a = F(i, j);
            b = F(i, (j+1)%3);
            c = F(i, (j+2)%3);
            // Length of edge crossing vertex a,b,c
            ea = l(i, j);
            eb = l(i, (j+1)%3);
            ec = l(i, (j+2)%3);
            // Compute cotangents for edge(a,b)
            cot = (ea*ea + eb*eb - ec*ec) / (2*A);
            triplets.emplace_back(a, b, cot);
        }
    }

    cots.setFromTriplets(triplets.begin(), triplets.end());

    triplets.clear();
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            a = F(i, j);
            b = F(i, (j+1)%3);
            triplets.emplace_back(a, b, cots.coeff(a, b)/2 + cots.coeff(b, a)/2);
        }
    }

    double sum = 0;
    for (int i = 0; i < nv; ++i) {
        for (int j = 0; j < nv; ++j) {
            if (j != i)  {
                sum += cots.coeff(i, j)/2 + cots.coeff(j, i)/2;
            }
        }
        triplets.emplace_back(i, i, -sum);
        sum = 0;
    }

    L.setFromTriplets(triplets.begin(), triplets.end());
}
