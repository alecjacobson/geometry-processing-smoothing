#include "cotmatrix.h"
#include "igl/doublearea.h"
#include <iostream>
#include <cmath>

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  
  // Compute areas of each face
  Eigen::VectorXd A;
  igl::doublearea(l, A);

  // Loop through faces to determine cotangent matrix elements
  // Using Eigen's functionality that triplets at duplicate (i, j)
  // are summed up when creating a sparse matrix
  typedef Eigen::Triplet<double> T;
  std::vector<T> triples;
  triples.reserve(F.rows() * 12);
  for (int i = 0; i < F.rows(); i++) {
        // Edge e = (F(j), F(j+1)) is opposite of F(j+2) mod 3
	// cotangent of the opposite angle is (-e^2 + l_1^2 + l_2^2) / (4 A)
	// where l_1, l_2 are the other sides and A is the area of current face. 
	double len_01 = ((pow(l(i, 0), 2) + pow(l(i,1), 2) - pow(l(i,2), 2))) / (4 * A(i)); 
	double len_12 = ((pow(l(i, 1), 2) + pow(l(i,2), 2) - pow(l(i,0), 2))) / (4 * A(i)); 
	double len_02 = ((pow(l(i, 0), 2) + pow(l(i,2), 2) - pow(l(i,1), 2))) / (4 * A(i)); 
        double entry_01 = len_01 * 0.5;
        double entry_12 = len_12 * 0.5;
  	double entry_02 = len_02 * 0.5;
        // Insert entries
        triples.push_back(T(F(i,0), F(i,1), entry_01)); 
        triples.push_back(T(F(i,1), F(i,0), entry_01)); 
        triples.push_back(T(F(i,0), F(i,0), -entry_01)); 
        triples.push_back(T(F(i,1), F(i,1), -entry_01)); 

        triples.push_back(T(F(i,0), F(i,2), entry_02)); 
        triples.push_back(T(F(i,2), F(i,0), entry_02)); 
        triples.push_back(T(F(i,0), F(i,0), -entry_02)); 
        triples.push_back(T(F(i,2), F(i,2), -entry_02)); 

        triples.push_back(T(F(i,1), F(i,2), entry_12)); 
        triples.push_back(T(F(i,2), F(i,1), entry_12)); 
        triples.push_back(T(F(i,1), F(i,1), -entry_12)); 
        triples.push_back(T(F(i,2), F(i,2), -entry_12)); 
  }
  L.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);
  L.setFromTriplets(triples.begin(), triples.end()); 
}

