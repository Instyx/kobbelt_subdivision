#include "../include/kobbelt_subdivision.hpp"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <vector>

namespace gcs = geometrycentral::surface;

void find_link(gcs::SurfaceMesh &mesh, gcs::Vertex v, std::vector<gcs::Vertex> &adj, std::vector<gcs::Vertex> &across){
    adj.clear();
    across.clear();
    for (gcs::Halfedge he : v.outgoingHalfedges()) {
        // Get the face associated with this halfedge.
        gcs::Face f = he.face();

        // If the halfedge is on the boundary, it has no face. Skip it.
        if (!he.isInterior()){
            continue;
        }

        assert(f.degree() == 4 && "This function requires quad faces.");
        adj.push_back(he.tipVertex());
        // The "across" vertex is the tip of the *next* halfedge in the face loop.
        //   v --> neighbor1 --> across_v --> neighbor2 --> v
        //   ^        ^             ^            ^
        //  he   he.next()  he.next().next() ...
        gcs::Halfedge next_he = he.next();
        gcs::Vertex across_v = next_he.tipVertex();

        across.push_back(across_v);
    }
}

void compute_edgepoints(gcs::SurfaceMesh &mesh, gcs::VertexPositionGeometry &geo, double w){
    std::vector<gcs::Vertex> adj_vertices;
    std::vector<gcs::Vertex> across_vertices;
    for(gcs::Edge e : mesh.edges()){
       gcs::Vertex p1 = e.firstVertex();
       gcs::Vertex p2 = e.secondVertex();
       // virtual vertices
       std::vector<gcs::Vertex> adj1, adj2, across1, across2;
       find_link(mesh, p1, adj1, across1);
       find_link(mesh, p2, adj2, across2);
    }
}
