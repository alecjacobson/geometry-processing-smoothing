#include "massmatrix.h"
#include "math.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  int f=F.rows();
  M.resize(F.maxCoeff()+1);
  Eigen::VectorXd area=Eigen::VectorXd::Zero(F.maxCoeff()+1);
  double s,a;
  for (int i=0; i<f; i++){
  	s=(l(i,0)+l(i,1)+l(i,2))/2;
  	a=sqrt(s*(s-l(i,0))*(s-l(i,1))*(s-l(i,2))); 
  	for (int j=0; j<3; j++){
  	  area(F(i,j))+=a;
	}
  }
  M=1.0/3*area.asDiagonal();
}

