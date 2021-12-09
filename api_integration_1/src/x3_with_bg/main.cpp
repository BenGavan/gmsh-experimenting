//
// Created by Ben Gavan on 18/07/2021.
//

#include <iostream>
#include <string>
#include "../../libs/gmsh.h"

using namespace std;

int main()
{

    gmsh::initialize();

    int t1 = gmsh::view::add("view 1");

    gmsh::open("t7.msh");

    int t2 = gmsh::model::mesh::field::add("PostView");

    gmsh::open("t7_bgmesh.pos");

    gmsh::model::mesh::field::setAsBackgroundMesh(t2);

    // In order to compute the mesh sizes from the background mesh only, and
    // disregard any other size constraints, one can set:
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0); // 0 => mesh generation not based on Mesh.x
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);


    gmsh::model::mesh::generate(2);


    gmsh::view::write(t1, "v1.msh");
    gmsh::view::write(t2, "t2.pos");

    gmsh::finalize();

    return 0;
}