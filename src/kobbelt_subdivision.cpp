#include "../include/kobbelt_subdivision.hpp"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include <geometrycentral/surface/surface_mesh_factories.h>
#include <Eigen/src/Core/Matrix.h>
#include <utility>
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

// same as finding the extrapolatedcorcer and applying 4 point rule
geometrycentral::Vector3 special_rule5(gcs::ManifoldSurfaceMesh &mesh, gcs::VertexPositionGeometry &geo, gcs::Vertex v0, gcs::Vertex v1, gcs::Vertex v2, double w){
    return (8.-w)/16. * geo.vertexPositions[v0] + (8.+2.*w)/16. * geo.vertexPositions[v1] - w/16. * geo.vertexPositions[v2];
}

geometrycentral::Vector3 four_point_rule(geometrycentral::Vector3 p0, geometrycentral::Vector3 p1, geometrycentral::Vector3 p2, geometrycentral::Vector3 p3, double w){
    return (8.+w)/16 * (p1+p2) - w/16 * (p0+p3);
}


std::pair<gcs::Vertex, gcs::Vertex> find_opposite_points(gcs::ManifoldSurfaceMesh &mesh, gcs::Face f, gcs::Vertex v, gcs::Vertex opp){
    gcs::Vertex v1, v2;
    int deg = f.degree();
    for(auto adj_f : v.adjacentFaces()){
        if(adj_f == f) continue;
        std::vector<gcs::Vertex> temp;
        for(auto v : adj_f.adjacentVertices()){
            temp.push_back(v);
        }
        int curr;
        for(int j = 0; j<temp.size(); ++j){
            if(temp[j]==v){
                curr = j;
                break;
            }
        }
        v1 = temp[(curr+4)%deg];
        break;
    }
    for(auto adj_f : opp.adjacentFaces()){
        if(adj_f == f) continue;
        std::vector<gcs::Vertex> temp;
        for(auto v : adj_f.adjacentVertices()){
            temp.push_back(v);
        }
        int curr;
        for(int j = 0; j<temp.size(); ++j){
            if(temp[j]==opp){
                curr = j;
                break;
            }
        }
        v2 = temp[(curr+4)%deg];
        break;
    }
    return {v1,v2};
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
        // corner vertices are the boundary vertices which is not connected to any inner vertex
        gcs::VertexData<bool> marked(mesh, false);
        for(auto he : mesh.exteriorHalfedges()){
            if(he.face()==he.next().face()){
                marked[he.tipVertex()] = true;
            }
        }
        // compute boundary edge points
        for(auto he : mesh.exteriorHalfedges()){
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

            if(!marked[p1] && !marked[p2]){
               edge_points[he.edge()] = four_point_rule(geo.vertexPositions[p0], geo.vertexPositions[p1], geo.vertexPositions[p2], geo.vertexPositions[p3], w);
            }
            if(marked[p1]){
                geometrycentral::Vector3 vv = compute_extrapolatedcorner(mesh, geo, p1, p2);
                edge_points[he.edge()] = four_point_rule(vv, geo.vertexPositions[p1], geo.vertexPositions[p2], geo.vertexPositions[p3], w);
            }
            if(marked[p2]){
                geometrycentral::Vector3 vv = compute_extrapolatedcorner(mesh, geo, p2, p1);
                edge_points[he.edge()] = four_point_rule(vv, geo.vertexPositions[p2], geo.vertexPositions[p1], geo.vertexPositions[p0], w);
            }
        }
    }

    // compute inner edge points
    for(gcs::Edge e : mesh.edges()){
        if(e.isBoundary()) continue;
        gcs::Vertex p1 = e.firstVertex();
        gcs::Vertex p2 = e.secondVertex();
        geometrycentral::Vector3 v1_pos = p1.isBoundary() ? compute_extrapolatedpoint(mesh, geo, p1) : compute_virtualpoint(mesh, geo, p1, p2, w);
        geometrycentral::Vector3 v2_pos = p2.isBoundary() ? compute_extrapolatedpoint(mesh, geo, p2) : compute_virtualpoint(mesh, geo, p2, p1, w);
        geometrycentral::Vector3 new_point = four_point_rule(v1_pos, geo.vertexPositions[p1], geo.vertexPositions[p2], v2_pos, w);
        edge_points[e] = new_point;
    }

    // insert the edge points to the mesh
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
        // create 4 faces for each face
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
        // compute the positiion of the face point
        //TODO: create 2 seperate functions for open and closed nets. a lot of checks are not needed in the closed =? faster?
        int count_b_edges = 0;
        for(auto e : f.adjacentEdges()){
            if(e.isBoundary()) ++count_b_edges;
        }
        // each edge counted twice, because of newly added edge points
        geometrycentral::Vector3 face_point{0.,0.,0.};
        count_b_edges/=2;
        if(count_b_edges==0){
            for(int i = 0; i< verts.size();++i){
                if(isOrigVert[verts[i]]) continue;
                gcs::Vertex opposite = verts[(i+4)%deg];
                std::pair<gcs::Vertex, gcs::Vertex> v_pair = find_opposite_points(mesh, f, verts[i], opposite);
                face_point = four_point_rule(geo.vertexPositions[v_pair.first], geo.vertexPositions[verts[i]], geo.vertexPositions[opposite], geo.vertexPositions[v_pair.second], w);
                break;
            }
        }
        else if(count_b_edges==1){
            for(int i = 0; i< verts.size();++i){
                if(isOrigVert[verts[i]]) continue;
                if(verts[i].isBoundary()) continue;
                gcs::Vertex opposite = verts[(i+4)%deg];
                if(opposite.isBoundary()) continue;
                std::pair<gcs::Vertex, gcs::Vertex> v_pair = find_opposite_points(mesh, f, verts[i], opposite);
                face_point = four_point_rule(geo.vertexPositions[v_pair.first], geo.vertexPositions[verts[i]], geo.vertexPositions[opposite], geo.vertexPositions[v_pair.second], w);
                break;
            }
        }
        else if(count_b_edges==2){
            int tt=0;
            for(int i = 0; i< verts.size();++i){
                if(tt==2) break;
                if(isOrigVert[verts[i]]) continue;
                gcs::Vertex opposite = verts[(i+4)%deg];
                std::pair<gcs::Vertex, gcs::Vertex> v_pair = find_opposite_points(mesh, f, verts[i], opposite);
                geometrycentral::Vector3 temp;
                if(verts[i].isBoundary()){
                    temp = special_rule5(mesh, geo, verts[i], opposite, v_pair.second, w);
                }
                else{
                    temp = special_rule5(mesh, geo, opposite, verts[i], v_pair.first, w);
                }
                face_point+=temp;
                ++tt;
            }
            face_point /= 2;
        }
        else if(count_b_edges==3){
            for(int i = 0; i< verts.size();++i){
                if(isOrigVert[verts[i]]) continue;
                gcs::Vertex opposite = verts[(i+4)%deg];
                if(verts[i].isBoundary() && opposite.isBoundary()) continue;
                std::pair<gcs::Vertex, gcs::Vertex> v_pair = find_opposite_points(mesh, f, verts[i], opposite);
                geometrycentral::Vector3 face_point;
                if(verts[i].isBoundary()){
                    face_point = special_rule5(mesh, geo, verts[i], opposite, v_pair.second, w);
                }
                else{
                    face_point = special_rule5(mesh, geo, opposite, verts[i], v_pair.first, w);
                }
                break;
            }
        }

        else{
            std::cout << "DISCONNECTED FACE!!" << std::endl;
        }

        newVertices[total_v+idx] = face_point;
        idx++;
    }
    return gcs::makeManifoldSurfaceMeshAndGeometry(newFaces, newVertices);
}



std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>, std::unique_ptr<gcs::VertexPositionGeometry> > k_kobbelt_subdivision(std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh, std::unique_ptr<gcs::VertexPositionGeometry> geo,
                                                                                                                           double w, int k)
{
    if (k <= 0) {
        return std::make_tuple(std::move(mesh), std::move(geo));
    }

    // Use the input pointers as the starting point for the loop.
    auto current_mesh = std::move(mesh);
    auto current_geo = std::move(geo);

    for (int i = 0; i < k; ++i) {
        std::tie(current_mesh, current_geo) = kobbelt_subdivision(*current_mesh, *current_geo, w);
    }

    return std::make_tuple(std::move(current_mesh), std::move(current_geo));
}
