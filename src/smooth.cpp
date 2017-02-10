#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "igl\edge_lengths.h"


void sparseify(const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & M, Eigen::SparseMatrix<double> & S) {
	int m = M.diagonal().size();
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(m);
	Eigen::VectorXd diag = M.diagonal();
	for (int i = 0; i < m; ++i) {
		tripletList.push_back(Eigen::Triplet<double>(i, i, diag(i)));
	}
	S.setFromTriplets(tripletList.begin(), tripletList.end());
}

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
	int m = F.rows();

	int nV = F.maxCoeff() + 1;
	Eigen::SparseMatrix<double> L(nV,nV);
	L.setZero();
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> M(nV);
	M.setZero();

	Eigen::SparseMatrix<double> SM(nV, nV);
	Eigen::SparseMatrix<double> A;
	Eigen::MatrixXd AL;
	Eigen::MatrixXd l(m, 3);

	igl::edge_lengths(V, F, l);
	cotmatrix(l, F, L);
	massmatrix(l, F, M);
	sparseify(M, SM);
	A = SM - lambda*L;
	//Eigen::LLT<Eigen::MatrixXd> lltOfA(A);
	//AL = lltOfA.matrixL();
		
	Eigen::MatrixXd b = M*G;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(A);
	U = chol.solve(b);

	std::cout << "DONE" << std::endl;
  

}

