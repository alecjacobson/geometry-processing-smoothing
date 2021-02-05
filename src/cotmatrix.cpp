#include "cotmatrix.h"

// Calcutes 1/2 * cot(theta) in a triangle with side lengths
// a, b, and c, where theta is the angle opposite c.
double halfcot(
    double a,
    double b, 
    double c) 
{
    // Heron's formula for area
    double s = (a + b + c) / 2;
    double A = std::sqrt(s * (s - a) * (s - b) * (s - c));

    // cos(theta) / sin(theta)
    return (std::pow(c, 2) - std::pow(a, 2) - std::pow(b, 2)) / (-8 * A);
}

void cotmatrix(
    const Eigen::MatrixXd & l,
    const Eigen::MatrixXi & F,
    Eigen::SparseMatrix<double> & L)
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(F.size() * 4);

    for (int i = 0; i < F.rows(); i++) 
    {
        Eigen::Vector3i f = F.row(i);
        Eigen::Vector3d lengths = l.row(i);

        for (int j = 0; j < f.size(); j++) 
        {
            // lengths[j] corresponds to the length of the edge 
            // between vertices f[(j + 1) % 3] and f[(j + 2) % 3]
            double c = halfcot(lengths[(j + 1) % 3], lengths[(j + 2) % 3], lengths[j]);
            triplets.push_back(T(f[(j + 1) % 3], f[(j + 2) % 3], c));
            triplets.push_back(T(f[(j + 2) % 3], f[(j + 1) % 3], c));

            // populate diagonal with negative sum of off-diagonal entries
            triplets.push_back(T(f[(j + 1) % 3], f[(j + 1) % 3], -c));
            triplets.push_back(T(f[(j + 2) % 3], f[(j + 2) % 3], -c));
        }
    }

    L.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);
    L.setFromTriplets(triplets.begin(), triplets.end());
}

