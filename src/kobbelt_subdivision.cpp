#include "../include/kobbelt_subdivision.hpp"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
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

geometrycentral::Vector3 four_point_rule(geometrycentral::Vector3 p0, geometrycentral::Vector3 p1, geometrycentral::Vector3 p2, geometrycentral::Vector3 p3, double w){
    return (8.+w)/16 * (p1+p2) - w/16 * (p0+p3);
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>, std::unique_ptr<gcs::VertexPositionGeometry> > kobbelt_subdivision(gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &geo, double w){
    std::vector<std::vector<size_t>> newFaces(mesh.nFaces()*4);
    std::vector<geometrycentral::Vector3> newVertices(mesh.nVertices() + mesh.nFaces() + mesh.nEdges());

    gcs::VertexData<gcs::Vertex> oldToNewVertexMap(mesh);
    gcs::EdgeData<gcs::Vertex> oldEdgeToNewVertexMap(mesh);
    gcs::FaceData<gcs::Vertex> oldFaceToNewVertexMap(mesh);
    gcs::VertexData<bool> isOrigVert(mesh, true);
    gcs::EdgeData<bool> isOrigEdge(mesh, true);
    geo.requireVertexPositions();

    for(gcs::Edge e : mesh.edges()){
       if (!isOrigEdge[e]) continue;
       gcs::Vertex p1 = e.firstVertex();
       gcs::Vertex p2 = e.secondVertex();
       geometrycentral::Vector3 v1_pos = compute_virtualpoint(mesh, geo, p1, p2, w);
       geometrycentral::Vector3 v2_pos = compute_virtualpoint(mesh, geo, p2, p1, w);
       geometrycentral::Vector3 new_point = four_point_rule(v1_pos, geo.vertexPositions[p1], geo.vertexPositions[p2], v2_pos, w);
       gcs::Vertex new_v = mesh.insertVertexAlongEdge(e).vertex();
       isOrigVert[new_v] = false;
       geo.vertexPositions[new_v] = new_point;
       oldEdgeToNewVertexMap[e] = new_v;
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
        for(int i = 0; i< verts.size();++i){
            if(isOrigVert[verts[i]]) continue;
            gcs::Vertex opposite = verts[(i+4)%deg];
            if(isOrigVert[opposite]){
                std::cout << "NANDATOOOOO " << std::endl;
            }
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
    return gcs::makeManifoldSurfaceMeshAndGeometry(newFaces, newVertices);
}
