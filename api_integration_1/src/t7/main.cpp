//
// Created by Ben Gavan on 15/07/2021.
//

#include <set>
#include <cmath>
#include <iostream>
#include "../../libs/gmsh.h"

using namespace std;

// Mesh sizes can be specified very accurately by providing a background mesh,
// i.e., a post-processing view that contains the target mesh sizes.

int main() {

    gmsh::initialize();

    // Merge a list-based post-processing view containing the target mesh sizes:
    try {
        gmsh::merge("t7_bgmesh.pos");
    } catch(...) {
        gmsh::logger::write("Could not load background mesh: bye!", "error");
        gmsh::finalize();
        return 0;
    }

    vector<int> tags;
    gmsh::view::getTags(tags);
    cout << "Tags: ";
    for (int i : tags)
    {
        cout << i << ", ";
    }
    cout << endl;

    gmsh::model::add("t7");

    // Create a simple rectangular geometry:
    double lc = 1e-2;
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(.1, 0, 0, lc, 2);
    gmsh::model::geo::addPoint(.1, .3, 0, lc, 3);
    gmsh::model::geo::addPoint(0, .3, 0, lc, 4);
    gmsh::model::geo::addLine(1, 2, 1);
    gmsh::model::geo::addLine(3, 2, 2);
    gmsh::model::geo::addLine(3, 4, 3);
    gmsh::model::geo::addLine(4, 1, 4);
    gmsh::model::geo::addCurveLoop({4, 1, -2, 3}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);
    gmsh::model::geo::synchronize();

    // Add the post-processing view as a new size field:
    int bg_field = gmsh::model::mesh::field::add("PostView");
    gmsh::model::mesh::field::setNumber(bg_field, "ViewIndex", 0);

    gmsh::model::mesh::field::setAsBackgroundMesh(bg_field);

    // In order to compute the mesh sizes from the background mesh only, and
    // disregard any other size constraints, one can set:
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0); // 0 => mesh generation not based on Mesh.x
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    gmsh::model::mesh::generate(2);

    gmsh::write("../../out/t7.msh");

    gmsh::finalize();

    return 0;
}



