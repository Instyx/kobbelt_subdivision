#include <iostream>
#include "../include/kobbelt_subdivision.hpp"
#include <geometrycentral/surface/meshio.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

namespace gcs = geometrycentral::surface;

std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh, newmesh;
std::unique_ptr<gcs::VertexPositionGeometry> geometry, newgeo;

int k = 1;
double w = 1;
double min_w = 0.0;
double max_w = 2.*(std::sqrt(5)-1);
bool subdivide;
void callback() {
  ImGui::PushItemWidth(216);

  ImGui::SameLine();
  ImGui::InputInt("subdivision level", &k);
  ImGui::InputDouble("parameter w", &w);
  ImGui::Checkbox("subdivide", &subdivide);
  if(subdivide){
      auto res = k_kobbelt_subdivision(std::move(mesh), std::move(geometry), w, k);
      mesh = std::move(std::get<0>(res));
      geometry = std::move(std::get<1>(res));

      polyscope::registerSurfaceMesh("subdivided", geometry->vertexPositions, mesh->getFaceVertexList());
      subdivide=false;
  }

}
int main(int argc, char *argv[]){
    polyscope::init();
    if(argc==2){
        std::tie(mesh, geometry) = gcs::readManifoldSurfaceMesh(argv[1]);
        polyscope::registerSurfaceMesh("my mesh", geometry->vertexPositions, mesh->getFaceVertexList());
        // for(auto f : mesh->faces()){
            // std::cout << f.degree() << std::endl;
        // }
        // std::tie(newmesh, newgeo) = kobbelt_subdivision(*mesh, *geometry, 1);
        // polyscope::registerSurfaceMesh("subdivided", newgeo->vertexPositions, newmesh->getFaceVertexList());
    }
    polyscope::state::userCallback = callback;
    polyscope::show();

    return 0;
}
