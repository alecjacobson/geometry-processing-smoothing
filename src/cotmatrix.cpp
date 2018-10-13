#include "cotmatrix.h"
#include "iostream"
#include <vector>
#include <math.h>

double cosineLaw(double a, double b, double c) {
    // https://en.wikipedia.org/wiki/Law_of_cosines
    double cosA = (-a * a + b * b + c * c) / (2 * c * b);
    return cosA;
}

double cot(double a, double b, double c) {
    //computing the cosA
    double cosA = cosineLaw(a, b, c);
    //since sin^2a + cos^2a=1.... then..
    double sinA = sqrt(1 - cosA * cosA);
    return cosA / sinA;
}

double cot(const Eigen::MatrixXd &l, int idx, int offset) {
    //overloaded function for simplicity.
    //columns correspond to edges [1,2],[2,0],[0,1] offset does a shift here.
    assert(l.cols() == 3);
    double b = l(idx, (offset + 0) % 3);
    double c = l(idx, (offset + 1) % 3);
    double a = l(idx, (offset + 2) % 3);
    return cot(a, b, c);
}

void to_L_sym(Eigen::SparseMatrix<double> &L, int i, int j, double val) {
    //This insert in symmetric fashion, and we need to sum Tb if Talpha already exists.
    L.coeffRef(i, j) += val;
    if (i != j)
        L.coeffRef(j, i) += val;
}

void cotmatrix(
        const Eigen::MatrixXd &l,
        const Eigen::MatrixXi &F,
        Eigen::SparseMatrix<double> &L) {

    int N = F.maxCoeff() + 1;
    L.resize(N, N);

    //I am going to use coeffRef on the SparseMatrix for simplicity(I am going to use the same matrix as an accumulator)
    //If speed were to be a  major concern instead of using coeffRef, we could fill the SparseMatrix with triplets as explained here:
    //https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling

    //In anycase, let's reserve some memory before hand to avoid multiple reallocations.
    L.reserve(F.rows() * 3 * 2);

    //this should be already zero,since it is a newly initialized sparse matrix but ...
    L.setZero();

    for (int idx = 0; idx < F.rows(); idx++) {
        int i = F(idx, 0);
        int j = F(idx, 1);
        int k = F(idx, 2);

        double cotA = 0.5 * cot(l, idx, 0);
        double cotB = 0.5 * cot(l, idx, 1);
        double cotC = 0.5 * cot(l, idx, 2);

        //to_L_sym-> add to existing value, symmetrically [eg. (i,j) (j,i)]
        to_L_sym(L, i, j, cotA);
        to_L_sym(L, j, k, cotB);
        to_L_sym(L, k, i, cotC);

    }

    //diagonal
    for (int idx = 0; idx < N; idx++) {
        L.coeffRef(idx, idx) = -L.row(idx).sum();
    }

    L.makeCompressed(); //suppresses the remaining empty space(if any) and transforms the matrix into a compressed column storage.

}

