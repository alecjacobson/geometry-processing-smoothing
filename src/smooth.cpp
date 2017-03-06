#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"
#include <igl/edge_lengths.h>
#include <igl/massmatrix.h>
#include <Eigen/SparseCholesky>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  /*
    Using the cotmatrix and massmatrix functions, smooth the mesh in V & F

    the mesh update rule is performed:

    Mᵗ⁺¹ Vᵗ = (Mᵗ - λLᵗ) Vᵗ⁺¹

    where the t is the current smoothing "time"
    We can assume that the changes in V have a negligible effect on L and M, and can discretize explicitly through

    Mᵗ Vᵗ = (Mᵗ - λLᵗ) Vᵗ⁺¹

    We can also generalize this to have an effect via smoothing of a mesh, or via smoothing a function along a mesh.
    This is done by using G as the function to be smoothed, and can be either vertices of the mesh or function along the surface:

    Mᵗ Gᵗ = (Mᵗ - λLᵗ) Vᵗ⁺¹

    Which can be sparsely solved
   */

  // get edge lengths
  Eigen::MatrixXd l(F.rows(), 3);
  igl::edge_lengths(V, F, l);

  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
  massmatrix(l, F, M);

  Eigen::SparseMatrix<double> L;
  cotmatrix(l, F, L);
  
  Eigen::SparseMatrix<double> A(M.rows(), M.rows());
  A.setIdentity();
 
  if(Eigen::MatrixXd(L).hasNaN())
    printf("L has nan\n");
  
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  U = solver.compute(A*M - lambda*L).solve(M*G);
  
  if(U.maxCoeff()  != U.maxCoeff())
    printf("NaN in solution\n");
}
