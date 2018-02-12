#include "cotmatrix.h"

// computes the cotangent of the angle opposite the side with length c, in a triangle with 
// side lengths a, b, and c.
double cotangent(double a, double b, double c) {

    // law of cosines
    double cos = (a*a + b*b - c*c)/(2*a*b);
    
    // Heron's formula for area of the triangle
    double s = (a + b + c)/2;
    double area = std::sqrt(s*(s-a)*(s-b)*(s-c));
    
    // diameter of the circumcircle
    double d = (a*b*c)/(2*area);
    
    // law of sines
    double sin = c/d;
    
    // cot = 1/tan = 1/(sin/cos) = cos/sin
    return cos/sin;
}

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L) {
    
    // compute the number of vertices and resize L accordingly.
    int n = F.maxCoeff() + 1;
    L.resize(n,n);
    
    // build a triplet list to populate L.
    typedef Eigen::Triplet<double> T;
    std::vector<T> list;
        
    for (int i = 0; i < F.rows(); i++) {
        
        // extract the vertex indices from F
        int v0 = F(i,0);
        int v1 = F(i,1);
        int v2 = F(i,2);
        
        // extract the edges lengths from l
        double e0 = l(i,2);
        double e1 = l(i,0);
        double e2 = l(i,1);
        
        // compute the half cotangent values
        double cot0 = 0.5*cotangent(e1,e2,e0);
        double cot1 = 0.5*cotangent(e0,e2,e1);
        double cot2 = 0.5*cotangent(e0,e1,e2);
        
        // add the non-diagonal entries. Because of the symmetry between edges, we add an 
        // entry for edge ij in both ij and ji. (Note that when populating the sparse matrix 
        // from the triplet list, duplicate entries are automatically summed. This means 
        // the 0.5*cot(beta) term will be added whenever we reach that face in the loop. 
        // Exploiting this feature also allows us to keep track of the rowwise sum in the 
        // diagonal entry for each row.)
        list.push_back(T(v0,v1,cot0));
        list.push_back(T(v1,v0,cot0));
        
        list.push_back(T(v1,v2,cot1));
        list.push_back(T(v2,v1,cot1));
        
        list.push_back(T(v2,v0,cot2));
        list.push_back(T(v0,v2,cot2));
                
        // add the diagonal elements
        list.push_back(T(v0,v0,-cot0));
        list.push_back(T(v1,v1,-cot0));
        
        list.push_back(T(v1,v1,-cot1));
        list.push_back(T(v2,v2,-cot1));
        
        list.push_back(T(v0,v0,-cot2));
        list.push_back(T(v2,v2,-cot2));

    } 
    
    // set L from the triplet list
    L.setFromTriplets(list.begin(), list.end());      
}

