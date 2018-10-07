#include "cotmatrix.h"
#include "math.h"

double cot(double a,double b,double c){
	double s=(a+b+c)/2;
	double area=sqrt(s*(s-a)*(s-b)*(s-c));
	double sin=area*2/b/c;
	double cos=(b*b+c*c-a*a)/2/b/c;
	return cos/sin;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  L.resize(F.maxCoeff()+1,F.maxCoeff()+1);
  int f=F.rows();
  double c;
  typedef Eigen::Triplet<double> T;
  std::vector<T> list;
  for (int i=0; i<f; i++){
  	for (int j=0; j<3; j++){
  	  c=cot(l(i,j),l(i,(j+1)%3),l(i,(j+2)%3))/2;
  	  list.push_back(T(F(i,(j+1)%3),F(i,(j+2)%3),c));
  	  list.push_back(T(F(i,(j+2)%3),F(i,(j+1)%3),c));
  	  list.push_back(T(F(i,(j+1)%3),F(i,(j+1)%3),-c));
  	  list.push_back(T(F(i,(j+2)%3),F(i,(j+2)%3),-c));
	}
  }
  L.setFromTriplets(list.begin(), list.end());
}

