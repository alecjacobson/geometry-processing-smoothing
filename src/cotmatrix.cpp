#include "cotmatrix.h"
#include <math.h> 

double myCot(double a, double b, double c){
	double cosC = (a*a + b*b - c*c) / (2*a*b);
	double sinC = sqrt(1-cosC*cosC);
	return cosC / sinC;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  int numV = F.maxCoeff()+1;
  L.resize(numV,numV);
  L.setZero();
  int aIdx, bIdx, cIdx;
  int Va, Vb, Vc;
  double cota, cotb, cotc;
	for (int ii = 0; ii < F.rows(); ii++){
		aIdx = 0;
		bIdx = 1;
		cIdx = 2;
		Va = F(ii, aIdx);
		Vb = F(ii, bIdx);
		Vc = F(ii, cIdx);
		cota = myCot(l(ii,bIdx),l(ii,cIdx),l(ii,aIdx));
		cotb = myCot(l(ii,aIdx),l(ii,cIdx),l(ii,bIdx));
		cotc = myCot(l(ii,aIdx),l(ii,bIdx),l(ii,cIdx));
		L.coeffRef(Va,Vb) += cotc / 2.0;
		L.coeffRef(Vb,Va) += cotc / 2.0;
		L.coeffRef(Vb,Vc) += cota / 2.0;
		L.coeffRef(Vc,Vb) += cota / 2.0;
		L.coeffRef(Va,Vc) += cotb / 2.0;
		L.coeffRef(Vc,Va) += cotb / 2.0;
	}

	for (int ii = 0; ii < numV; ii++){
		L.coeffRef(ii,ii) = -L.row(ii).sum();
	}
}

