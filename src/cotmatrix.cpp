#include "cotmatrix.h"
#include <Eigen/Sparse>
#include <math.h>
#include <ppl.h>
#include <concurrent_vector.h>
inline void vertex_index(const int i, int& il, int& ir)
{
  il = i - 1 < 0 ? 2 : i - 1;
  ir = i + 1 > 2 ? 0 : i + 1;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  typedef Eigen::Triplet<double> tuple;
  Concurrency::concurrent_vector<tuple> tuple_list;
  const int n = F.rows();
  tuple_list.reserve(12 * n);
  Concurrency::parallel_for(size_t(0),size_t(n), [&l,&F,n,&tuple_list](const int m)
  {
    int il, ir;
    int v[3];
    double e[3];
    double s = 0;
    for (int i = 0; i < 3; i++)
    {
      v[i] = F(m, i);
      e[i] = l(m, i);
      s += 0.5*e[i];
    }
    double A = sqrt(s*(s - e[0])*(s - e[1])*(s - e[2]));
    double d = 2 * A / (e[0] * e[1] * e[2]);

    double cosine[3];
    double sine[3];
    double cot[3];

    for (int i = 0; i < 3; i++)
    {
      vertex_index(i, il, ir);
      cosine[i] = (e[il] * e[il] + e[ir] * e[ir] - e[i] * e[i]) / (2 * e[il] * e[ir]);
      sine[i] = e[i] * d;
      cot[i] = cosine[i] / sine[i];
    }

    for (int i = 0; i < 3; i++)
    {
      vertex_index(i, il, ir);
      tuple_list.push_back(tuple(v[il], v[ir], 0.5*cot[i]));
      tuple_list.push_back(tuple(v[ir], v[il], 0.5*cot[i]));
      tuple_list.push_back(tuple(v[il], v[il], -0.5*cot[i]));
      tuple_list.push_back(tuple(v[ir], v[ir], -0.5*cot[i]));
    } 
  }, Concurrency::auto_partitioner());
  L.setFromTriplets(tuple_list.begin(), tuple_list.end());
}