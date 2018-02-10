#include "cotmatrix.h"

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  std::vector< Eigen::Triplet<double> > tripletList;
  // Preprocess: create a hash of edge ij to face index
  std::map<std::string, int> edge_to_f;
  int num_faces = F.rows();
  for(int f_idx= 0; f_idx < num_faces; f_idx++){
    int v0 = F(f_idx, 0); int v1 = F(f_idx, 1); int v2 = F(f_idx, 2);
    // Create string keys
    std::string key01 = std::to_string(v0) + std::to_string(v1);
    std::string key12 = std::to_string(v1) + std::to_string(v2);
    std::string key20 = std::to_string(v2) + std::to_string(v0);
    // Add mappings to edge_to_f
    edge_to_f.insert(make_pair(key01,f_idx));
    edge_to_f.insert(make_pair(key12,f_idx));
    edge_to_f.insert(make_pair(key20,f_idx));
  }
  // Create triplets of L
  int V = F.maxCoeff() + 1;
  double sum_off_diagonals = 0;
  for(int i = 0; i < V; i ++){
    for(int j = 0; j < V; j ++){
      // Find indicies into F of two faces that are incident on edge i j
      std::string keyij = std::to_string(i) + std::to_string(j);
      std::string keyji = std::to_string(j) + std::to_string(i);
      int fa_idx = edge_to_f[keyij];
      int fb_idx = edge_to_f[keyji];
      if(fa_idx != 0 && fb_idx != 0){
        // found the edge, compute cotangent a_ij and b_ij
        double cotangent_a = cotangent_triangle(i,j, fa_idx, l, F);
        double cotangent_b = cotangent_triangle(j,i, fb_idx, l, F);
        double cotangent_sum = (0.5 * cotangent_a) + (0.5 * cotangent_b);
        sum_off_diagonals += cotangent_sum;
        tripletList.push_back(Eigen::Triplet<double>(i, j, cotangent_sum));
      }
    }
  }
  // Diagonal entries
  for(int v = 0; v < V; v++){
    tripletList.push_back(Eigen::Triplet<double>(v, v, -sum_off_diagonals));
  }
  // Populate L
  L.setFromTriplets(tripletList.begin(), tripletList.end());
}

double cotangent_triangle(
  int i,
  int j,
  int f_idx,
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F)
  {
    // find i's index in F -- this is a lil ugly
    int i_pos; int j_pos; int k_pos;
    for(int x = 0; x < 3; x++){
      if(F(f_idx, x) == i){
        i_pos = x;
        j_pos = (x + 1) % 3;
        k_pos = (x + 2) % 3;
      }
    }
    // Find edge length of edge ij
    double l_ij = l(f_idx,i_pos);
    double l_jk = l(f_idx, j_pos);
    double l_ki = l(f_idx, k_pos);

    // Herons formula + solve for diameter of triangle's circumference
    double sp = (l_ij + l_jk + l_ki) / 2; // semiperimeter
    double area = sqrt(sp * (sp - l_jk) * (sp - l_ki) * (sp * l_ij));
    double d = (l_ij * l_jk * l_ki) / (2 * area);

    // From sin and cos laws
    double sin_theta = l_ij / d;
    double cos_theta = ((l_jk*l_jk) + (l_ki*l_ki) - (l_ij * l_ij)) / (2 * l_jk * l_ki);
    return 0.0;
  }
