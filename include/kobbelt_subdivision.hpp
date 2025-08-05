#pragma once

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace gcs = geometrycentral::surface;

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>, std::unique_ptr<gcs::VertexPositionGeometry> > compute_edgepoints(gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &geo, double w);
