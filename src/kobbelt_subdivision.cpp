#include "../include/kobbelt_subdivision.hpp"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <geometrycentral/surface/surface_mesh_factories.h>
#include <Eigen/src/Core/Matrix.h>
#include <vector>
#include <memory>

namespace gcs = geometrycentral::surface;

geometrycentral::Vector3 compute_virtualpoint(gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &geo, gcs::Vertex v1, gcs::Vertex v2, double w){
    std::vector<gcs::Vertex> adj;
    std::vector<gcs::Vertex> across;
    for (gcs::Halfedge he : v1.outgoingHalfedges()) {
        // Get the face associated with this halfedge.
        gcs::Face f = he.face();

        // If the halfedge is on the boundary, it has no face. Skip it.
        if (!he.isInterior()){
            continue;
        }
        std::cout << f.degree() << std::endl;
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
    int idx = -1;
    for(int i = 0; i<adj.size(); ++i){
        if(adj[i]==v2){
            idx=i;
            break;
        }
    }
    geometrycentral::Vector3 v1_pos{0.,0.,0.};
    int n1 = v1.degree();
    for(auto v : adj){
        v1_pos += geo.vertexPositions[v];
    }
    v1_pos = (4. / n1) * v1_pos;
    v1_pos -= (geo.vertexPositions[adj[(idx-1+n1)%n1]] + geo.vertexPositions[adj[idx]] + geo.vertexPositions[adj[(idx+1)%n1]]);
    geometrycentral::Vector3 temp = w/(4.+w) * (geo.vertexPositions[across[(idx-2+n1)%n1]] + geo.vertexPositions[across[(idx-1+n1)%n1]] +
                                                geo.vertexPositions[across[idx]] + geo.vertexPositions[across[(idx+1)%n1]]);
    v1_pos += temp;
    temp = geometrycentral::Vector3{0.,0.,0.};
    for(auto v : across){
        temp += geo.vertexPositions[v];
    }
    v1_pos -= (4.*w / ((4.+w)*n1)) * temp;

    return v1_pos;

}


geometrycentral::Vector3 compute_extrapolatedpoint(gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &geo, gcs::Vertex b){
    std::vector<gcs::Vertex> inner_adj;
    for(gcs::Vertex n : b.adjacentVertices()){
        if(n.isBoundary()) continue;
        inner_adj.push_back(n);
    }
    geometrycentral::Vector3 temp{0.,0.,0.};
    for(auto n : inner_adj){
        temp -= geo.vertexPositions[n];
    }
    temp = temp / inner_adj.size();
    return 2*geo.vertexPositions[b] + temp;
}

geometrycentral::Vector3 compute_extrapolatedcorner(gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &geo, gcs::Vertex v1, gcs::Vertex v2){
    //v1 should be the corner to be extrapolated
    return 2*geo.vertexPositions[v1] - geo.vertexPositions[v2];
}


geometrycentral::Vector3 four_point_rule(geometrycentral::Vector3 p0, geometrycentral::Vector3 p1, geometrycentral::Vector3 p2, geometrycentral::Vector3 p3, double w){
    return (8.+w)/16 * (p1+p2) - w/16 * (p0+p3);
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>, std::unique_ptr<gcs::VertexPositionGeometry> > kobbelt_subdivision(gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &geo, double w){
    std::vector<std::vector<size_t>> newFaces;
    std::vector<geometrycentral::Vector3> newVertices(mesh.nVertices() + mesh.nFaces() + mesh.nEdges());

    gcs::VertexData<bool> isOrigVert(mesh, true);
    gcs::EdgeData<bool> isOrigEdge(mesh, true);
    geo.requireVertexPositions();

    gcs::EdgeData<geometrycentral::Vector3> edge_points(mesh);
    if(mesh.hasBoundary()){
        //mark corner vertices
        gcs::VertexData<bool> marked(mesh, false);
        for(auto he : mesh.exteriorHalfedges()){
            if(he.face()==he.next().face()){
                marked[he.tipVertex()] = true;
            }
        }

        for(auto he : mesh.exteriorHalfedges()){
            if(!marked[he.tailVertex()] && !marked[he.tipVertex()]){
                gcs::Vertex p0;
                gcs::Vertex p1 = he.tailVertex();
                gcs::Vertex p2 = he.tipVertex();
                gcs::Vertex p3 = he.next().tipVertex();
                for(auto o_he : p1.outgoingHalfedges()){
                    if(o_he.tipVertex()!=p2 && o_he.tipVertex().isBoundary()){
                        p0 = o_he.tipVertex();
                        break;
                    }
                }
                edge_points[he.edge()] = four_point_rule(geo.vertexPositions[p0], geo.vertexPositions[p1], geo.vertexPositions[p2], geo.vertexPositions[p3], w);
            }
        }

    for(gcs::Edge e : mesh.edges()){
        gcs::Vertex p1 = e.firstVertex();
        gcs::Vertex p2 = e.secondVertex();
        geometrycentral::Vector3 v1_pos = compute_virtualpoint(mesh, geo, p1, p2, w);
        geometrycentral::Vector3 v2_pos = compute_virtualpoint(mesh, geo, p2, p1, w);
        geometrycentral::Vector3 new_point = four_point_rule(v1_pos, geo.vertexPositions[p1], geo.vertexPositions[p2], v2_pos, w);
        edge_points[e] = new_point;
    }

    for(gcs::Edge e : mesh.edges()){
        if (!isOrigEdge[e]) continue;
        gcs::Vertex new_v = mesh.insertVertexAlongEdge(e).vertex();
        isOrigVert[new_v] = false;
        geo.vertexPositions[new_v] = edge_points[e];
        for (gcs::Edge ee : new_v.adjacentEdges()) {
            isOrigEdge[ee] = false;                  // mark the new edges
           // gcs::Vertex otherV = e.otherVertex(new_v);    // other side of edge
       }

    }

    for(auto v: mesh.vertices()){
        newVertices[v.getIndex()] = geo.vertexPositions[v];
    }
    int total_v = mesh.nVertices();
    int idx = 0;
    for(gcs::Face f : mesh.faces()){
        std::vector<gcs::Vertex> verts;
        for(gcs::Vertex v : f.adjacentVertices()){
            verts.push_back(v);
        }
        int deg = verts.size();
        for(int i = 0; i < verts.size(); ++i){
            std::vector<size_t> face;
            if(isOrigVert[verts[i]]){
                face.push_back(verts[(i+deg-1)%deg].getIndex());
                face.push_back(verts[(i)%deg].getIndex());
                face.push_back(verts[(i+1)%deg].getIndex());
                face.push_back(total_v+idx);
                newFaces.push_back(face);
            }
        }
        //TODO: create 2 seperate functions for open and closed nets. a lot of checks are not needed in the closed =? faster?
        for(int i = 0; i< verts.size();++i){
            if(isOrigVert[verts[i]]) continue;
            gcs::Vertex opposite = verts[(i+4)%deg];
            if(isOrigVert[opposite]){
                std::cout << "NANDATOOOOO " << std::endl;
            }
            //TODO: check for boundary cases and extrapolate the newly added mid points in the boundary
            // 3 case, either
            //                     only 1 edge of the face is adjacent to boundary (compute 4 point rule for the other direction)
            //                     2 edges (extrapoalte mid points in both directions and average)
            //                     3 edges (like 2 edges but no average since there is only one possibility?)
            bool flag = false;
            for(gcs::Edge adj_e : verts[i].adjacentEdges() ){
                if(adj_e.isBoundary()){
                    flag = true;
                    break;
                }
            }
            if(flag) continue;
            for(gcs::Edge adj_e : opposite.adjacentEdges() ){
                if(adj_e.isBoundary()){
                    flag = true;
                    break;
                }
            }
            if(flag) continue;
            gcs::Vertex v1, v2;
            for(auto adj_f : verts[i].adjacentFaces()){
                if(adj_f == f) continue;
                std::vector<gcs::Vertex> temp;
                for(auto v : adj_f.adjacentVertices()){
                    temp.push_back(v);
                }
                int curr;
                for(int j = 0; j<temp.size(); ++j){
                    if(temp[j]==verts[i]){
                        curr = j;
                        break;
                    }
                }
                v1 = temp[(curr+4)%deg];
                break;
            }
            for(auto adj_f : opposite.adjacentFaces()){
                if(adj_f == f) continue;
                std::vector<gcs::Vertex> temp;
                for(auto v : adj_f.adjacentVertices()){
                    temp.push_back(v);
                }
                int curr;
                for(int j = 0; j<temp.size(); ++j){
                    if(temp[j]==opposite){
                        curr = j;
                        break;
                    }
                }
                v2 = temp[(curr+4)%deg];
                break;
            }
            geometrycentral::Vector3 face_point = four_point_rule(geo.vertexPositions[v1], geo.vertexPositions[verts[i]], geo.vertexPositions[opposite], geo.vertexPositions[v2], w);
            newVertices[total_v+idx] = face_point;
            break;
        }
        idx++;
    }
    std::cout << "end " << std::endl;
    for(auto a : newFaces)
        std::cout << a.size() << std::endl;
    return gcs::makeManifoldSurfaceMeshAndGeometry(newFaces, newVertices);
}
