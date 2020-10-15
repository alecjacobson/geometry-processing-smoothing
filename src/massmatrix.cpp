#include "massmatrix.h"
#include "igl/doublearea.h"

void massmatrix(
        const Eigen::MatrixXd &l,
        const Eigen::MatrixXi &F,
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> &M) {


    int N = F.maxCoeff() + 1;
    M.resize(N);
    M.setZero();

    Eigen::MatrixXd areas;
    igl::doublearea(l, areas);

    //igl::doublearea computes the area of each triangle twice.
    areas = areas / 2.0;

    for (int idx = 0; idx < F.rows(); idx++) {
        //let's get the Vs that are part of this triangle.
        int v1 = F(idx, 0);
        int v2 = F(idx, 1);
        int v3 = F(idx, 2);

        //the area of this triangle F(idx) contributes 1/3 to each of this vs(v1,v2,v3), so..
        M.diagonal()(v1) += areas(idx) / 3.0;
        M.diagonal()(v2) += areas(idx) / 3.0;
        M.diagonal()(v3) += areas(idx) / 3.0;
    }

}

