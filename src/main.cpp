#include <iostream>
#include "../include/kobbelt_subdivision.hpp"
#include <geometrycentral/surface/meshio.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

namespace gcs = geometrycentral::surface;

int main(int argc, char *argv[]){
    polyscope::init();
    std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
    std::unique_ptr<gcs::VertexPositionGeometry> geometry;
    if(argc==2){
        std::tie(mesh, geometry) = gcs::readManifoldSurfaceMesh(argv[1]);
        polyscope::registerSurfaceMesh("my mesh", geometry->vertexPositions, mesh->getFaceVertexList());
    }
    polyscope::show();

    return 0;
}
