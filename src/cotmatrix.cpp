#include "cotmatrix.h"

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  int V_num = F.maxCoeff() + 1;
  L.resize(V_num, V_num);

  auto cot = [](double a, double b, double c){
    double cosA = (b * b + c * c - a * a) * 1.0 / (2 * b * c);
    double sinA = sqrt(1 - cosA * cosA);
    return cosA / sinA;
  };

  for (int i = 0; i < F.rows(); i++) {
    double edge0 = l(i, 0);
    double edge1 = l(i, 1);
    double edge2 = l(i, 2);

    double cot0 = cot(edge0, edge1, edge2);
    double cot1 = cot(edge1, edge0, edge2);
    double cot2 = cot(edge2, edge0, edge1);

    int v0 = F(i, 0), v1 = F(i, 1), v2 = F(i, 2);
    L.coeffRef(v0, v1) += 0.5 * cot2;
    L.coeffRef(v1, v0) += 0.5 * cot2;
    L.coeffRef(v1, v2) += 0.5 * cot0;
    L.coeffRef(v2, v1) += 0.5 * cot0;
    L.coeffRef(v0, v2) += 0.5 * cot1;
    L.coeffRef(v2, v0) += 0.5 * cot1;

    L.coeffRef(v0, v0) -= 0.5 * cot1;
    L.coeffRef(v0, v0) -= 0.5 * cot2;
    L.coeffRef(v1, v1) -= 0.5 * cot0;
    L.coeffRef(v1, v1) -= 0.5 * cot2;
    L.coeffRef(v2, v2) -= 0.5 * cot0;
    L.coeffRef(v2, v2) -= 0.5 * cot1;
  }

  L.makeCompressed();
}

