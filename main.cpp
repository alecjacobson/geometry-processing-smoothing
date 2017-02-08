#include "smooth.h"
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/parula.h>
#include <igl/viewer/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
  // Load input meshes
  Eigen::MatrixXd OV,V,U;
  Eigen::MatrixXi F;
  double lambda = 1e-5;
  igl::read_triangle_mesh(
    (argc>1?argv[1]:"../shared/data/sphere-noisy.obj"),OV,F);
  // Load data into MatrixXd rather than VectorXd for simpler `smooth` API
  // Just use y-coordinates as data to be smoothed
  Eigen::MatrixXd G = OV.col(1);
  if(argc>2)
  {
    if(argv[2][0] == 'n')
    {
      // Corrupt with a bit of noise
      G += 0.1*(G.maxCoeff()-G.minCoeff())*
        Eigen::MatrixXd::Random(G.rows(),G.cols());
    }else
    {
      igl::readDMAT(argv[2],G);
    }
  }

  // Create a libigl Viewer object to toggle between point cloud and mesh
  igl::viewer::Viewer viewer;
  std::cout<<R"(
  D,d  smooth data
  K    decamate(?) lambda
  k    decimate lambda
  L    Toggle lighting
  M,m  smooth mesh geometry
  R,r  reset mesh geometry and data
  L    lighting
)";
  const auto & update = [&]()
  {
    if((V.array() != V.array()).any())
    {
      std::cout<<"Too degenerate to keep smoothing. Better reset"<<std::endl;
    }
    viewer.data.set_mesh(V,F);
    viewer.data.compute_normals();
    Eigen::MatrixXd C;
    igl::parula(U,G.minCoeff(),G.maxCoeff(),C);
    viewer.data.set_colors(C);
  };
  const auto & reset = [&]()
  {
    V = OV;
    U = G;
  };
  viewer.callback_key_pressed = 
    [&](igl::viewer::Viewer &, unsigned int key, int)
  {
    switch(key)
    {
      case 'D':
      case 'd':
        //////////////////////////////////////////////////////////////////////
        // Smooth data
        //////////////////////////////////////////////////////////////////////
        // Use copy constructor to fake in-place API (may be overly
        // conservative depending on your implementation)
        smooth(V,F,Eigen::MatrixXd(U),lambda,U);
        break;
      case 'K':
      case 'k':
        lambda = (key=='K'?10.0:0.1)*lambda;
        std::cout<<"lambda: "<<lambda<<std::endl;
        break;
      case 'L':
        // Toggle lighting
        viewer.core.lighting_factor = 1.0- viewer.core.lighting_factor;
        break;
      case 'M':
      case 'm':
      {
        //////////////////////////////////////////////////////////////////////
        // Smooth mesh geometry. 
        //////////////////////////////////////////////////////////////////////
        // "Linearize" simply by conducting smooth step assuming that vertex
        // data is a signal defined over current surface: copy is needed to
        // prevent memory "aliasing"
        Eigen::MatrixXd Vcpy(V);
        smooth(Vcpy,F,Vcpy,lambda,V);
        break;
      }
      case 'R':
      case 'r':
        reset();
        break;
      default:
        return false;
    }
    update();
    return true;
  };

  reset();
  update();
  viewer.launch();

  return EXIT_SUCCESS;
}
