#include "cotmatrix.h"

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
    int n = F.maxCoeff() + 1;
    L.resize(n,n);
    for(int i = 0; i < n; i++){
        for(int f = 0; f < F.rows(); f++){
            int face_index = -1;
            for(int k = 0; k < 3; k++){
                //find which vertex on face Ff is the our vertex Vi
                if(F(f,k) == i){
                    face_index = k;
                }
            }
            if(face_index >= 0){
                //We have two half edges involving vi,
                //The edge from F(f,face_index) to F(f,face_index+1)
                //and from F(f,face_index-1) to F(f,face_index)
                double c = l(f,(face_index+2)%3); //length of edge from i to F(f,(face_index+1)%3)
                double a = l(f,(face_index+1)%3);
                double b = l(f,face_index);
                double cosine = (b*b + a*a - c*c) / (2*a*b);
                double sine = std::sqrt(1-cosine*cosine);
                L.coeffRef(i,F(f,(face_index+1)%3)) += 0.5 * cosine/sine;
                
                c = l(f,(face_index+1)%3); //length of edge from i to F(f,(face_index+2)%3)
                a = l(f,(face_index+2)%3);
                b = l(f,face_index);
                cosine = (b*b + a*a - c*c) / (2*a*b);
                sine = std::sqrt(1-cosine*cosine);
                L.coeffRef(i,F(f,(face_index+2)%3)) += 0.5 * cosine/sine;
            }
        }
    }
    Eigen::VectorXd rowsum = L * Eigen::VectorXd::Ones(L.rows());
    
    for(int i = 0; i < L.rows(); i++){
        L.coeffRef(i,i) -= rowsum(i);
    }
}

